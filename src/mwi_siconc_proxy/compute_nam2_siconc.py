"""compute_nam2_siconc

Compute the sea ice concentration for NOAA AMSR2 swath files

"""

import os
import sys
from pathlib import Path
import numpy as np
from netCDF4 import Dataset

#import pm_icecon.bt.compute_bt_ic as bt
from pm_icecon.bt.compute_bt_ic import (
    tb_data_mask,
    _get_wx_params,
    get_water_mask,
    bootstrap,
    coastal_fix,
    apply_invalid_ice_mask,
)
from pm_icecon.config.models.bt import BootstrapParams
from pm_icecon.config.models.bt import TbSetParams, WeatherFilterParams, WeatherFilterParamsForSeason
from pm_icecon.tests.regression.test_bt import _original_f18_example

from mwi_siconc_proxy.set_siconc_parameters import get_bootstrap_parameters
from mwi_siconc_proxy.grid_am2mbt import (
    read_am2mbt_config,
    parse_am2mbt_filename,
    GRIDID_AREADEFS,
)
from mwi_siconc_proxy.set_siconc_parameters import (
    read_btparams_file,
)

CRS_WKTS = {
    'e2n': 'PROJCRS["WGS 84 / NSIDC EASE-Grid 2.0 North",BASEGEOGCRS["WGS 84",ENSEMBLE["World Geodetic System 1984 ensemble",MEMBER["World Geodetic System 1984 (Transit)"],MEMBER["World Geodetic System 1984 (G730)"],MEMBER["World Geodetic System 1984 (G873)"],MEMBER["World Geodetic System 1984 (G1150)"],MEMBER["World Geodetic System 1984 (G1674)"],MEMBER["World Geodetic System 1984 (G1762)"],MEMBER["World Geodetic System 1984 (G2139)"],MEMBER["World Geodetic System 1984 (G2296)"],ELLIPSOID["WGS 84",6378137,298.257223563,LENGTHUNIT["metre",1]],ENSEMBLEACCURACY[2.0]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4326]],CONVERSION["US NSIDC EASE-Grid 2.0 North",METHOD["Lambert Azimuthal Equal Area",ID["EPSG",9820]],PARAMETER["Latitude of natural origin",90,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8801]],PARAMETER["Longitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["False easting",0,LENGTHUNIT["metre",1],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["easting (X)",south,MERIDIAN[90,ANGLEUNIT["degree",0.0174532925199433]],ORDER[1],LENGTHUNIT["metre",1]],AXIS["northing (Y)",south,MERIDIAN[180,ANGLEUNIT["degree",0.0174532925199433]],ORDER[2],LENGTHUNIT["metre",1]],USAGE[SCOPE["Environmental science - used as basis for EASE grid."],AREA["Northern hemisphere."],BBOX[0,-180,90,180]],ID["EPSG",6931]]',
    'e2s': 'PROJCRS["WGS 84 / NSIDC EASE-Grid 2.0 South",BASEGEOGCRS["WGS 84",ENSEMBLE["World Geodetic System 1984 ensemble",MEMBER["World Geodetic System 1984 (Transit)"],MEMBER["World Geodetic System 1984 (G730)"],MEMBER["World Geodetic System 1984 (G873)"],MEMBER["World Geodetic System 1984 (G1150)"],MEMBER["World Geodetic System 1984 (G1674)"],MEMBER["World Geodetic System 1984 (G1762)"],MEMBER["World Geodetic System 1984 (G2139)"],MEMBER["World Geodetic System 1984 (G2296)"],ELLIPSOID["WGS 84",6378137,298.257223563,LENGTHUNIT["metre",1]],ENSEMBLEACCURACY[2.0]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4326]],CONVERSION["US NSIDC EASE-Grid 2.0 South",METHOD["Lambert Azimuthal Equal Area",ID["EPSG",9820]],PARAMETER["Latitude of natural origin",-90,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8801]],PARAMETER["Longitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["False easting",0,LENGTHUNIT["metre",1],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["easting (X)",north,MERIDIAN[90,ANGLEUNIT["degree",0.0174532925199433]],ORDER[1],LENGTHUNIT["metre",1]],AXIS["northing (Y)",north,MERIDIAN[0,ANGLEUNIT["degree",0.0174532925199433]],ORDER[2],LENGTHUNIT["metre",1]],USAGE[SCOPE["Environmental science - used as basis for EASE grid."],AREA["Southern hemisphere."],BBOX[-90,-180,0,180]],ID["EPSG",6932]]',  # noqa
}


def get_run_label(path):
    """Determine the runtime label from a Path"""
    label = path.stem

    return label

def get_date_from_swath_list(path):
    """Extract date from last endtime in path"""
    try:
        with open(path) as f:
            lines = f.readlines()
            last_filename = lines[-1].rstrip()
            fn_parts = parse_am2mbt_filename(last_filename)
            date = fn_parts['timestamp_end']
    except ValueError:
        raise RuntimeError(f'Could not get date from swath list: {path}')

    return date


def get_runtime_info_siconc(
    swath_list_path,
    config,
    btparams,
):
    """Determine the information needed here"""
    info = {}
    info['gridid'] = config['gridid']
    info['area_def'] = GRIDID_AREADEFS[info['gridid']]
    info['tb_list'] = ['h37', 'v19', 'v22', 'v37']
    info['ancillary_dir'] = config['ancillary_dir']
    info['output_dir'] = config['output_dir']
    info['siconc_label'] = get_run_label(swath_list_path)

    info['tb_fn_template'] = config['tb_fn_template']
    info['sic_fnstem_template'] = config['sic_fnstem_template']
    tb_dat_fname_template = config['tb_fn_template']

    # These are probably duplicates...
    info['outputdir'] = config['output_dir']
    info['prefix'] = info['siconc_label']

    info['tb_dat_fnames'] = {}
    for tb in info['tb_list']:
        info['tb_dat_fnames'][tb] = tb_dat_fname_template.format(tb=tb, **info)
        if not Path(info['tb_dat_fnames'][tb]).is_file():
            print(f'Warning no such tb dat file: {info["tb_dat_fnames"][tb]}')

    info['date'] = get_date_from_swath_list(swath_list_path)
    info['xdim'] = info['area_def']['width']
    info['ydim'] = info['area_def']['height']

    return info


def get_tbs_from_swath_list(info, dtype=np.float32):
    """Read in the appropriate TB files"""
    # This will hold the tb data arrays
    tbs = {}

    #xdim = info['area_def']['width']
    #ydim = info['area_def']['height']
    xdim = info['xdim']
    ydim = info['ydim']
    for tb in info['tb_list']:
        print(f'tb: {tb}')
        print(f'  fn: {info["tb_dat_fnames"][tb]}')
        tb_fn = info["tb_dat_fnames"][tb]
        tbs[tb] = np.fromfile(tb_fn, dtype=dtype).reshape(ydim, xdim)

    return tbs  # dict of np-arrays with tb name as key


def interpret_oldstyle_mask(masktype, fn, xdim, ydim):
    """Get this type of mask from the old-style land/ocean mask file"""
    data = np.fromfile(fn, dtype=np.int16).reshape(ydim, xdim)

    if masktype == 'nonocean_mask':
        return np.array(data > 3)
    if masktype == 'valid_ice':
        return np.array(data > 0)
    if masktype == 'invalid_ice':
        return np.array(data == 0)
    else:
        raise RuntimeError(f'dont know how to get: {masktype}')

def get_land_mask(info):
    """Return the land mask for this gridid"""
    gridid = info['gridid']
    anc_dir = info['ancillary_dir']
    if not Path(anc_dir).is_dir():
        raise RuntimeError(f'No such anc dir: {anc_dir}')

    # Note: This is specific to the old-style NOAA sicond data format
    try:
        land_mask_fname = f'{anc_dir}/{gridid}_mask.dat'
        land_mask = interpret_oldstyle_mask(
            'nonocean_mask',
            land_mask_fname,
            xdim=info['xdim'],
            ydim=info['ydim'],
        )
    except FileNotFoundError:
        raise RuntimeError(f'Unknown land_mask_fname for gridid: {gridid}\n{land_mask_fname}')

    return land_mask

    
def get_invalid_ice_mask(info):
    """Return the invalid mask for this gridid"""
    gridid = info['gridid']
    anc_dir = info['ancillary_dir']
    if not Path(anc_dir).is_dir():
        raise RuntimeError(f'No such anc dir: {anc_dir}')
    # Note: mask uses date of file, not of each individual swath
    date = info['date']
    month = date.month
    xdim = info['xdim']
    ydim = info['ydim']

    ## Note: This is specific to the old-style NOAA sicond data format
    #if gridid == 'e2n10':
    #    sst_mask_fn = f'{anc_dir}/e2_sst_n_{month:02d}.dat'
    #    sst_mask = interpret_oldstyle_mask('invalid_ice', sst_mask_fn, xdim=xdim, ydim=ydim)
    #elif gridid == 'e2s10':
    #    sst_mask_fn = f'{anc_dir}/e2_sst_s_{month:02d}.dat'
    #    sst_mask = interpret_oldstyle_mask('invalid_ice', sst_mask_fn, xdim=xdim, ydim=ydim)
    #else:
    #    raise RuntimeError(f'Unknown sst for gridid: {gridid}, month: {month}')

    # Note: This is specific to the old-style NOAA sicond data format
    try:
        sst_mask_fname = f'{anc_dir}/{gridid}_sst_{month:02d}.dat'
        sst_mask = interpret_oldstyle_mask(
            'invalid_ice',
            sst_mask_fname,
            xdim=info['xdim'],
            ydim=info['ydim'],
        )
    except FileNotFoundError:
        raise RuntimeError(f'Unknown land_mask_fname for gridid: {gridid}')

    return sst_mask


def fill_pole_hole(gridid, conc):
    """Simple fill the pole hole routine"""

    # if the middle of the land mask is land, don't need to fill pole hole
    pole_radius = None
    if gridid == 'e2n10':
        pole_radius = 15
        polehole_ileft = 513
        polehole_iright = 538
        polehole_jtop = 513
        polehole_jbottom = 538
    elif gridid == 'e2n20':
        pole_radius = 8
        polehole_ileft = 254
        polehole_iright = 271
        polehole_jtop = 254
        polehole_jbottom = 271

    # Return if the gridid indicates no pole hole filling needed
    if pole_radius is None:
        print('No pole hole calculated')
        return conc

    # Note: this is a view, so the values will update when this
    #       subset is updated
    near_pole_conc = conc[
        polehole_jtop:polehole_jbottom,
        polehole_ileft:polehole_iright
    ]

    is_pole_hole = (near_pole_conc < 0.01) | (near_pole_conc > 100)
    if np.any(np.isfinite(near_pole_conc[is_pole_hole])):
        near_pole_mean = np.nanmean(near_pole_conc[~is_pole_hole])
        near_pole_conc[is_pole_hole] = near_pole_mean
    else:
        print(f'No suitable values found near pole hole; not filled')


    return conc


def compute_nam2_siconc(
    tbs,
    btparams,
    land_mask,
    invalid_ice_mask,
    date,
    gridid,
    land_value=254,
    missing_value=255,
):
    """Compute the NOAA AMSR2 sea ice conc using Bootstrap algorithm
    This is modeled after the 'goddard_bootstrap' routine in pm_icecon"""
    return_fields = {}

    tb_mask = tb_data_mask(
        tbs=(
            tbs['v37'],
            tbs['h37'],
            tbs['v19'],
            tbs['v22'],
        ),
        min_tb=btparams.mintb,
        max_tb=btparams.maxtb,
    )
    return_fields['tb_mask'] = tb_mask

    season_params = _get_wx_params(
        date=date,
        weather_filter_seasons=btparams.weather_filter_seasons,
    )

    # Note: weather_mask is called water_mask in BT code
    weather_mask = get_water_mask(
        v37=tbs['v37'],
        h37=tbs['h37'],
        v22=tbs['v22'],
        v19=tbs['v19'],
        land_mask=land_mask,
        tb_mask=tb_mask,
        ln1=btparams.vh37_params.lnline,
        date=date,
        wintrc=season_params.wintrc,
        wslope=season_params.wslope,
        wxlimt=season_params.wxlimt,
    )
    return_fields['weather_mask'] = weather_mask

    siconc = bootstrap(
        tb_v37=tbs['v37'],
        tb_h37=tbs['h37'],
        tb_v19=tbs['v19'],
        params=btparams,
        land_mask=land_mask,
        tb_mask=tb_mask,
        water_mask=weather_mask,
        missing_flag_value=missing_value,
    )
    return_fields['siconc_init'] = siconc.copy()

    siconc[weather_mask] = 0
    siconc[tb_mask] = missing_value
    siconc[land_mask] = land_value

    siconc_prior = siconc.copy()
    siconc = apply_invalid_ice_mask(
        conc=siconc,
        missing_flag_value=missing_value,
        land_flag_value=land_value,
        invalid_ice_mask=invalid_ice_mask,
    )
    applied_invalid = siconc != siconc_prior
    return_fields['applied_invalid'] = applied_invalid
    return_fields['invalid_ice'] = invalid_ice_mask

    siconc_prior = siconc.copy()
    siconc = coastal_fix(
        conc=siconc,
        missing_flag_value=missing_value,
        land_mask=land_mask,
        minic=btparams.minic,
    )
    applied_spillover = siconc_prior != siconc
    return_fields['applied_spillover'] = applied_spillover

    # Here, pole hole is filled with missing_value (255)
    siconc_prior = siconc.copy()
    siconc = fill_pole_hole(gridid, siconc)
    #applied_polehole = (siconc_prior != siconc) & np.isfinite(siconc)
    applied_polehole = siconc != siconc_prior
    return_fields['applied_polehole'] = applied_polehole

    siconc[siconc < btparams.minic] = 0

    # Add a few extra fields to keep things neat
    return_fields['land_mask'] = land_mask

    return siconc, return_fields

def get_siconc_path(info, fn_extension='.dat'):
    """Determine the output file neme from the info dict"""
    fn_wo_ext = info['sic_fnstem_template'].format(**info)
    filename = fn_wo_ext + fn_extension

    return filename

def calculate_flag_field(siconc, return_fields):
    """Determine the flag values
    Note: the order of operations matters"""
    flag_field = np.zeros(siconc.shape, dtype=np.uint8)

    flag_field[np.isnan(siconc)] = 255
    flag_field[return_fields['tb_mask']] = 128

    flag_field[return_fields['weather_mask']] = 8
    flag_field[return_fields['applied_spillover']] = 16

    flag_field[return_fields['invalid_ice']] = 4
    flag_field[return_fields['applied_polehole']] = 64
    flag_field[return_fields['land_mask']] = 120

    flag_meanings = '0: Good, 4: InvalidSeaiceClimatology, 8: Weather, 16: Spillover, 64: FilledPoleHole, 120: Land, 128: MissingData'

    return flag_field, flag_meanings


def save_to_netcdf(
    info,
    siconc,
    flag_field,
    flag_field_meanings,
    return_fields,
):
    """Save the file to output netCDF file"""
    # Get the file name
    nc_fn = info['sic_fnstem_template'].format(**info) + '.nc'
    ds = Dataset(nc_fn, 'w')

    xdim = info['xdim']
    ydim = info['ydim']

    gridid = info['gridid']
    xleft = info['area_def']['area_extent'][0]
    xright = info['area_def']['area_extent'][2]
    ytop = info['area_def']['area_extent'][3]
    ybottom = info['area_def']['area_extent'][1]

    # Note: Assuming x and y resolution are the same
    resolution = (xright - xleft) / xdim
    halfres = resolution / 2

    # Set the x variable (dimension)
    x = ds.createDimension('x', xdim)  # noqa
    xs = ds.createVariable('x', np.float32, ('x',))
    xs_values = np.linspace(
        xleft + halfres,
        xright - halfres,
        num=xdim, 
        dtype=np.float32)
    xs[:] = xs_values[:]
    xs.standard_name = 'projection_x_coordinate'
    xs.long_name = 'x coordinate of projection'  
    xs.units = 'meters'
    xs.axis = 'X'

    # Set the y variable (dimension)
    y = ds.createDimension('y', ydim)  # noqa
    ys = ds.createVariable('y', np.float32, ('y',))
    # Note that y values start from "top" (positive y)
    ys_values = np.linspace(
        ytop - halfres,
        ybottom + halfres,
        num=ydim,
        dtype=np.float32)
    ys[:] = ys_values[:]
    ys.standard_name = 'projection_y_coordinate'
    ys.long_name = 'y coordinate of projection'
    ys.units = 'meters'
    ys.axis = 'Y'

    # Create the coordinate reference system variable
    if gridid == 'e2n10':
        long_name = "NSIDC_EASE2_N10km_subset"
        lat_proj_orig = 90.
        crs_wkt = CRS_WKTS['e2n']
    elif gridid == 'e2s10':
        long_name = "NSIDC_EASE2_S10km_subset"
        lat_proj_orig = -90.
        crs_wkt = CRS_WKTS['e2s']
    elif gridid == 'e2n20':
        long_name = "NSIDC_EASE2_N20km_subset"
        lat_proj_orig = 90.
        crs_wkt = CRS_WKTS['e2n']
    elif gridid == 'e2s20':
        long_name = "NSIDC_EASE2_S20km_subset"
        lat_proj_orig = -90.
        crs_wkt = CRS_WKTS['e2s']

    crs = ds.createVariable('crs', 'i4')
    crs.grid_mapping_name = 'lambert_azimuthal_equal_area'
    crs.longitude_of_projection = 0.
    crs.false_easting = 0.
    crs.false_northing = 0.
    crs.semi_major_axis = 6378137.
    crs.inverse_flattening = 298.257223563
    crs.GeoTransform = "-5250000 10000 0 5250000 0 -10000"
    crs.long_name = long_name
    crs.latitude_of_projection_origin = lat_proj_orig
    crs.crs_wkt = crs_wkt

    # Note: this is a 2d variable (no time dimension)
    ubyte_siconc = siconc.copy()
    ubyte_siconc[ubyte_siconc < 0] = 0
    ubyte_siconc[ubyte_siconc > 255] = 255
    ubyte_siconc = ubyte_siconc.astype(np.uint8)

    conc = ds.createVariable('seaice_conc', np.uint8, ('y', 'x'), zlib=True)  # noqa
    conc[:, :] = ubyte_siconc[:, :]
    conc.coverage_content_type = 'image'

    conc.grid_mapping = 'crs'
    conc.units = '1'
    conc.long_name = 'Bootstrap Sea Ice Concentration'
    conc.standard_name = 'sea_ice_area_fraction'
    conc.valid_range = np.array((0, 100))
    conc.flag_values = np.array((254, 255))
    conc.flag_meanings = 'land missing'
    conc.packing_convention = 'netCDF'
    conc.packing_convention_description = 'unpacked = scale_factor * packed + add_offset'
    conc.scale_factor = 0.01
    conc.add_offset = 0.

    flags = ds.createVariable('quality_flags', np.uint8, ('y', 'x'), zlib=True)  # noqa
    flags[:, :] = flag_field[:, :]
    flags.long_name = 'Quality Flags'
    flags.grid_mapping = 'crs'
    flags.units = 1
    flags.comment = flag_field_meanings

    # Note: this is a 2d variable (no time dimension)
    # Note: No multiyear sea ice algorithm is implemented yet
    multiyear_ice = np.zeros(siconc.shape, dtype=np.float32)

    myi_conc = ds.createVariable('multiyear_seaice_conc', np.uint8, ('y', 'x'), zlib=True)  # noqa
    myi_conc[:, :] = multiyear_ice[:, :]
    myi_conc.coverage_content_type = 'image'

    myi_conc.grid_mapping = 'crs'
    myi_conc.units = '1'
    myi_conc.long_name = 'Multiyear Sea Ice Concentration'
    myi_conc.standard_name = 'sea_ice_area_fraction'
    myi_conc.valid_range = np.array((0, 100))
    myi_conc.flag_values = np.array((254, 255))
    myi_conc.flag_meanings = 'land missing'
    myi_conc.packing_convention = 'netCDF'
    myi_conc.packing_convention_description = 'unpacked = scale_factor * packed + add_offset'
    myi_conc.scale_factor = 0.01
    myi_conc.add_offset = 0.
    myi_conc.comment = 'The multiyear sea ice algorithm is not yet implemented'
    
    ds.close()
    print(f'Wrote: {nc_fn}')


def main_compute_nam2_siconc():
    """Check that we're were"""
    try:
        swath_list_filename = sys.argv[1]
        swath_list_path = Path(swath_list_filename)
        assert swath_list_path.is_file()
        run_label = get_run_label(swath_list_path)
    except IndexError:
        raise RuntimeError('First arg of compute_nam2_siconc.py should be a swath list filename, whose base-name is used as the label for this run.')
    except AssertionError:
        raise RuntimeError(f'No such swath file list: {swath_list_filename}')

    try:
        config_filename = sys.argv[2]
        config = read_am2mbt_config(config_filename)
    except IndexError:
        raise RuntimeError('Second arg of compute_nam2_siconc.py should be config file name')

    try:
        btparams_filename = sys.argv[3]
        btparams = read_btparams_file(btparams_filename)
    except IndexError:
        raise RuntimeError('Third arg of compute_nam2_siconc.py should be bootstrap parameters file')

    # Initialize info structure
    info = get_runtime_info_siconc(
        swath_list_path,
        config,
        btparams,
    )

    tbs = get_tbs_from_swath_list(info)
    land_mask = get_land_mask(info)
    invalid_ice_mask = get_invalid_ice_mask(info)
    siconc, return_fields = compute_nam2_siconc(
        tbs=tbs,
        btparams=btparams,
        land_mask=land_mask,
        invalid_ice_mask=invalid_ice_mask,
        date=info['date'],
        gridid=info['gridid'],
    )

    output_filename = get_siconc_path(info)
    siconc.tofile(output_filename)

    flag_field, flag_field_meanings = calculate_flag_field(siconc, return_fields)

    save_to_netcdf(
        info,
        siconc,
        flag_field,
        flag_field_meanings,
        return_fields,
    )


if __name__ == '__main__':
    main_compute_nam2_siconc()
