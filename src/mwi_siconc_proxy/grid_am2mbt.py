""" Grid NOAA AMSR2 MBT to NSIDC NH EASE2 10km Arctic subset

See get_usage_string() for example usage
"""

from netCDF4 import Dataset
import numpy as np
import os
import datetime as dt
import yaml
from pathlib import Path
from pyresample.geometry import AreaDefinition, SwathDefinition
from pyresample.bilinear import NumpyBilinearResampler
import warnings


'''Options for NumpyBilinearResampler:
def resample_bilinear(data, source_geo_def, target_area_def, radius=50e3,
                      neighbours=32, nprocs=1, fill_value=0,
                      reduce_data=True, segments=None, epsilon=0):
    """Resample using bilinear interpolation.

    data : numpy array
        Array of single channel data points or
        (source_geo_def.shape, k) array of k channels of datapoints
    source_geo_def : object
        Geometry definition of source data
    target_area_def : object
        Geometry definition of target area
    radius : float, optional
        Cut-off distance in meters
    neighbours : int, optional
        Number of neighbours to consider for each grid point when
        searching the closest corner points
    nprocs : int, optional
        Number of processor cores to be used for getting neighbour info
    fill_value : {int, None}, optional
        Set undetermined pixels to this value.
        If fill_value is None a masked array is returned with undetermined
        pixels masked
    reduce_data : bool, optional
        Perform initial coarse reduction of source dataset in order
        to reduce execution time
    segments : int or None
        Number of segments to use when resampling.
        If set to None an estimate will be calculated
    epsilon : float, optional
        Allowed uncertainty in meters. Increasing uncertainty
        reduces execution time
'''
EXAMPLE_AM2MBT_FNAME = 'AMSR2-MBT_v2r2_GW1_s202501060005067_e202501060144046_c202501060156059.nc'
EXAMPLE_CONFIG_FNAME = 'am2mbt_e2n10_config.yaml'
DEFINED_AM2MBT_ERRORS = [
    'nofile',
]
DEFAULT_CONFIG_AM2MBT = {
    'var_fn_template': '{dirname}/{gridid}_{varname}_{tstamp}.dat',
    'gridid': 'e2n10',
    'output_dir': './output',
    'intermed_dir': './intermed',
    'overwrite': False,
    # nc_vars maps the to-be-gridded fields to their AMSR2-MBT file var names
    'nc_vars': {
        # tbs
        'v18': 'Brightness_Temperature_18_GHzV',
        'v23': 'Brightness_Temperature_23_GHzV',
        'h36': 'Brightness_Temperature_36_GHzH',
        'v36': 'Brightness_Temperature_36_GHzV',

        # other fields
        'landflag': 'Land_Ocean_Flag_6_to_36',

        # lats
        'latlow': 'Latitude_for_Low_Resolution',
        'lat18': 'Latitude_for_18',
        'lat23': 'Latitude_for_23',
        'lat36': 'Latitude_for_36',

        # lons
        'lonlow': 'Longitude_for_Low_Resolution',
        'lon18': 'Longitude_for_18',
        'lon23': 'Longitude_for_23',
        'lon36': 'Longitude_for_36',
    },
    # grid_vars specifies the data and lat/lon fields used to geolocate the data
    'grid_vars': {
        'v19': ['v18', 'lat18', 'lon18'],
        'v22': ['v23', 'lat23', 'lon23'],
        'h37': ['h36', 'lat36', 'lon36'],
        'v37': ['v36', 'lat36', 'lon36'],
    },
}

GRIDID_AREADEFS = {
    'latlon_5deg': {
        'area_id': 'latlon_5deg_grid',
        'description': 'Simple_latlon_grid_at_10_deg',
        'proj_id': 'sample_latlon',
        'projection': 'EPSG:4326',
        'width': 72,
        'height': 18,
        'area_extent': (-180, -90, 180, 90),
    },
    'e2n10': {
        'area_id': 'ease2_nh_arctic_10km',
        'description': 'Arctic Sea Ice subset of EASE2 Northern Hemisphere',
        'proj_id': 'ease2_nh_arctic_10km',
        'projection': 'EPSG:6931',
        'width': 1050,
        'height': 1050,
        'area_extent': (-5250000, -5250000, 5250000, 5250000),
    },
    'e2n20': {
        'area_id': 'ease2_nh_arctic_20km',
        'description': 'Arctic Sea Ice subset of EASE2 Northern Hemisphere',
        'proj_id': 'ease2_nh_arctic_20km',
        'projection': 'EPSG:6931',
        'width': 525,
        'height': 525,
        'area_extent': (-5250000, -5250000, 5250000, 5250000),
    },
    'e2s10': {
        'area_id': 'ease2_sh_antarctic_10km',
        'description': 'Antarctic Sea Ice subset of EASE2 Southern Hemisphere',
        'proj_id': 'ease2_sh_antarctic_10km',
        'projection': 'EPSG:6932',
        'width': 840,
        'height': 840,
        'area_extent': (-4200000, -4200000, 4200000, 4200000),
    },
    'e2s20': {
        'area_id': 'ease2_sh_antarctic_20km',
        'description': 'Antarctic Sea Ice subset of EASE2 Southern Hemisphere',
        'proj_id': 'ease2_sh_antarctic_20km',
        'projection': 'EPSG:6932',
        'width': 420,
        'height': 420,
        'area_extent': (-4200000, -4200000, 4200000, 4200000),
    },
    'e2n10_NP_4x4': {
        'area_id': 'E2N10_NP_test',
        'description': 'Test grid near North Pole',
        'proj_id': 'e2n20_np_test',
        'projection': 'EPSG:6931',
        'width': 4,
        'height': 4,
        'area_extent': (-20000, -20000, 20000, 20000),
    },
    'e2n10_NP_10x10': {
        'area_id': 'E2N10_NP_test',
        'description': 'Test grid near North Pole',
        'proj_id': 'e2n20_np_test',
        'projection': 'EPSG:6931',
        'width': 10,
        'height': 10,
        'area_extent': (-50000, -50000, 50000, 50000),
    },
}


def xwm(m='exiting in xwm()'):
    # Note: this is only for development and should be used in production
    raise SystemExit(m)  # pragma: no cover


def get_pyresample_area_def(gridid):
    """Return the pyresample module's AreaDefinition for a specified grid"""
    try:
        assert gridid in GRIDID_AREADEFS.keys()
    except AssertionError:
        raise ValueError(f'No AreaDefinition information for: {gridid}')

    gridid_area_def = GRIDID_AREADEFS[gridid]
    area_def = AreaDefinition(
        gridid_area_def['area_id'],
        gridid_area_def['description'],
        gridid_area_def['proj_id'],
        gridid_area_def['projection'],
        gridid_area_def['width'],
        gridid_area_def['height'],
        gridid_area_def['area_extent'],
    )

    return area_def


def get_pyresample_swathdef(lons, lats):
    """Return the pyresample module's SwathDefinition for a set of lons, lats"""
    swathdef = SwathDefinition(lons=lons, lats=lats)

    return swathdef


def get_usage_string():
    bfn = os.path.basename(__file__)

    usage_str = f'''
    Usage:
        {bfn} <AMSR2_MBT_ncfile> [<config_yaml_fname>]
    eg:
        {bfn} {EXAMPLE_AM2MBT_FNAME}
        {bfn} {EXAMPLE_AM2MBT_FNAME} {EXAMPLE_CONFIG_FNAME}
    '''

    return usage_str


def get_error_message(grid_mbt_error):
    error_str = ''
    if grid_mbt_error == 'nofile':
        error_str += '\n'
        error_str += 'No input file provided'
        error_str += '\n'
        error_str += get_usage_string()

    return error_str


def parse_swath_timestamp(ts, ts_type='start'):
    """Extract the timestamp from a std AMSR2-MBT file name part
    This part will be 1 letter followed by 15 digits"""
    try:
        assert ts_type in ('start', 'end', 'creation')
    except AssertionError:
        ts_type = 'unknown'

    try:
        if ts_type != 'unknown':
            assert ts_type[0] == ts[0]
    except AssertionError:
        raise ValueError(f'timestamp type {ts_type} mismatch first char {ts[0]}')

    try:
        year=int(ts[1:5])
        month=int(ts[5:7])
        day=int(ts[7:9])
        hour=int(ts[9:11])
        minute=int(ts[11:13])
        seconds=int(ts[13:15])
        microseconds=int(ts[15:16]) * 100000
        timestamp_datetime = \
            dt.datetime(year, month, day, hour, minute, seconds, microseconds)
    except ValueError:
        raise ValueError(f'Could not parse {ts_type} timestamp: {ts}')

    return timestamp_datetime


def parse_am2mbt_filename(fn):
    """Return dict with standard parts of AMSR2 MBT file name
    Note: A typical file is EXAMPLE_AM2MBT_FNAME:
    AMSR2-MBT_v2r2_GW1_s202501060005067_e202501060144046_c202501060156059.nc
    which has parts:
      AMSR2-MBT
      v2r2
      GW1
      s202501060005067
      e202501060144046
      c202501060156059
    """
    fn_path = Path(fn)
    fn_stem = fn_path.stem

    am2mbt_parts = {}
    am2mbt_parts['is_ncfile'] = fn_path.suffix == '.nc'

    fn_parts = fn_stem.split('_')

    am2mbt_parts['is_standard'] = True
    if len(fn_parts) != 6:
        am2mbt_parts['is_standard'] = False
    else:
        # Parse the filename as expected
        # prod_name, prod_ver, prod_plat, ts_start, ts_end, ts_gen = \
        am2mbt_parts['product_name'] = fn_parts[0]
        am2mbt_parts['product_version'] = fn_parts[1]
        am2mbt_parts['product_platform'] = fn_parts[2]
        try:
            am2mbt_parts['timestamp_start'] = parse_swath_timestamp(fn_parts[3], 'start')
            am2mbt_parts['timestamp_end'] = parse_swath_timestamp(fn_parts[4], 'end')
            am2mbt_parts['timestamp_creation'] = parse_swath_timestamp(fn_parts[5], 'creation')
        except ValueError:
            print(f'Some timestamps not parsed')
            am2mbt_parts['is_standard'] = False

    if not am2mbt_parts['is_ncfile']:
        print(f'WARNING: this is not an ncfile:\n  {fn}')
        am2mbt_parts['is_standard'] = False

    return am2mbt_parts


def get_filename_timestamp(fn):
    fn_parts = parse_am2mbt_filename(fn)

    try:
        return fn_parts['timestamp_start']
    except KeyError:
        raise ValueError(f'No timestamp_start found in: {fn}')


def verify_am2mbt_dataset(am2ds, config):
    '''Confirm that the netCDF4 Dataset has suitable variables'''
    data_vars = [varname for varname in am2ds.variables.keys()]

    nc_vars = [varname for varname in config['nc_vars'].keys()]
    grid_vars = [varname for varname in config['grid_vars'].keys()]

    # Verify that each var needed for grid_var is in nc_vars
    for grid_var in grid_vars:
        contrib_vars = config['grid_vars'][grid_var]
        for contrib_var in contrib_vars:
            try:
                assert contrib_var in nc_vars
            except AssertionError:
                raise ValueError(f'Could not find {contrib_var} of {contrib_vars} for {grid_var} in {nc_vars}')

    # Verify that each contributing var is in the netCDF dataset
    for nc_var in nc_vars:
        expected_data_var = config['nc_vars'][nc_var]
        try:
            assert expected_data_var in data_vars
        except AssertionError:
            raise ValueError(f'Could not find {expected_data_var} in {data_vars}')

    # End of verify_am2mbt_dataset().
    # If control reaches here, dataset is verified.

    ### DEBUG and PREPARATION ###
    # Methodology for saving the config file as it exists here
    #    out_fn = 'am2mbt_e2n10_config.yaml'
    #    with open(out_fn, 'w') as f:
    #        yaml.dump(config, f)
    #    print(f'Wrote config to: {out_fn}')


def verify_am2mbt_file(
    fn,
    config=DEFAULT_CONFIG_AM2MBT,
    default_timestamp=dt.datetime(1900, 1, 1, 1, 1, 1, 100000)
):
    '''Confirm that the netCDF4 file contains suitable variables'''

    # Set timestamp to default if can't be calculated
    try:
        filename_tstamp = get_filename_timestamp(fn)
    except ValueError:
        filename_tstamp = default_timestamp

    try:
        am2ds = Dataset(fn)
    except OSError:
        raise OSError(f'Could not open input file as netCDF: {fn}')

    verify_am2mbt_dataset(am2ds, config)

    filename_info = parse_am2mbt_filename(fn)
    return filename_info


def read_am2mbt_config(fn):
    from yaml.scanner import ScannerError

    try:
        with open(fn, 'r') as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f'Could not file config file: {fn}')
    except ScannerError:
        raise ValueError(f'Bad yaml file: {fn}')

    if config is None:
        raise ValueError(f'config is None read from: {fn}')

    return config


# TODO: handle radius
def grid_swath_data(vals, lons, lats, areadef, radius=20e3, nanval=0, label=None):
    # Grid the data to an areadef

    swathdef = get_pyresample_swathdef(lons=lons, lats=lats)

    resampler = NumpyBilinearResampler(swathdef, areadef, radius)

    if label is not None:
        print(f'resampling swath data ({label})...', end='', flush=True)
    with warnings.catch_warnings(record=True) as warning_list:
        gridded = resampler.resample(vals).astype(np.float32)

        for warning in warning_list:
            message = warning.message

            # Ignore proj messages about using CRS codes
            if isinstance(message, UserWarning) and 'converting to a PROJ string' in str(message):
                pass
            else:
                print(f'warning: {message}')

    if label is not None:
        print('done', flush=True)

    # Set NaN values
    if nanval is not None:
        gridded[gridded == nanval] = np.nan

    return gridded


def get_runtime_info(fn, config, overwrite=None):
    """Determine the information needed to create the data files"""
    fn_info = parse_am2mbt_filename(fn)
    runtime_info = {}
    if overwrite is not None:
        runtime_info['overwrite'] = overwrite
    else:
        runtime_info['overwrite'] = config['overwrite']

    runtime_info['timestamp_str'] = fn_info['timestamp_start'].strftime('%Y%m%d%H%M')
    runtime_info['tb_source_str'] = f"{fn_info['product_name']}_{fn_info['product_version']}"
    runtime_info['gridid'] = config['gridid']
    runtime_info['var_fn_template'] = config['var_fn_template']
    runtime_info['output_dir'] = config['output_dir']
    runtime_info['intermed_dir'] = config['intermed_dir']

    return runtime_info


def set_up_directories(runtime_info):
    # Ensure that directories specified in runtime_info exist
    os.makedirs(runtime_info['output_dir'], exist_ok=True)
    os.makedirs(runtime_info['intermed_dir'], exist_ok=True)


def get_intermed_var_filename(varname, runtime_info):
    # Determine the intermediate filename for gridded swath variable
    fn_template = runtime_info['var_fn_template']
    intermed_dir = runtime_info['intermed_dir']
    timestamp_str = runtime_info['timestamp_str']
    gridid = runtime_info['gridid']

    #'var_fn_template': 'DIRNAME/VAR_TSTAMP.dat',
    #fn_template = '{dirname}/{varname}_{tstamp}.dat'
    intermed_fn = fn_template.format(
        dirname=intermed_dir,
        varname=varname,
        tstamp=timestamp_str,
        gridid=gridid,
    )

    return intermed_fn


def process_am2mbt_file(fn, config, process_overwrite=None):  # pragma: no cover
    """Extract information from an AMSR2 MBT netCDF file"""
    #print(f'Processing:')
    #print(f'  file: {fn}')
    #print(f'  config: {config}')
    print(f'Processing file: {fn}')

    runtime_info = get_runtime_info(fn, config, overwrite=process_overwrite)
    overwrite = runtime_info['overwrite']

    set_up_directories(runtime_info)

    # 'grid_vars': {
    #     'v19': ['v18', 'lat18', 'lon18'],
    # 'nc_vars': {
    #     'v18': 'Brightness_Temperature_18_GHzV',
    ds = Dataset(fn, 'r')

    # TODO: Also grid time field
    for grid_var in config['grid_vars']:
        grid_var_filename = get_intermed_var_filename(grid_var, runtime_info)
        if not os.path.isfile(grid_var_filename):
            print(f'Creating: {grid_var_filename}')
        else:
            if overwrite:
                print(f'Overwriting: {grid_var_filename}')
            else:
                print(f'File exists: {grid_var_filename}')
                continue

        data_var, lat_var, lon_var = config['grid_vars'][grid_var]
        data_nc_varname = config['nc_vars'][data_var]
        lat_nc_varname = config['nc_vars'][lat_var]
        lon_nc_varname = config['nc_vars'][lon_var]

        data_arr = np.array(ds.variables[data_nc_varname])
        lat_arr = np.array(ds.variables[lat_nc_varname])
        lon_arr = np.array(ds.variables[lon_nc_varname])

        areadef = get_pyresample_area_def(config['gridid'])
        gridded = grid_swath_data(data_arr, lon_arr, lat_arr, areadef, label=grid_var)

        gridded.tofile(grid_var_filename)
        print(f'  Wrote: {grid_var_filename}', flush=True)


def main():
    """
    checks for existence of AMSR2-MBT swath file
    if config file is specified,
      verifies existence of config file
      verifies valididy of config file
    Then calls the processing routine
    """

    import sys

    # Check for a second argument with name of YAML with config info
    try:
        cfn = sys.argv[2]
        assert os.path.isfile(cfn)
        config = read_am2mbt_config(cfn)
    except IndexError:
        # No config file specified
        #config = DEFAULT_CONFIG_AM2MBT
        raise RuntimeError(f'No config file specified')
    except AssertionError:
        raise ValueError(f'No such specified config file: {cfn}')
    except ValueError:
        raise ValueError(f'Could not parse specified config file: {cfn}')

    # Get name of file to process from first cmdline arg
    #   and do some minimal validation of the file
    try:
        ifn = sys.argv[1]
        assert os.path.isfile(ifn)
    except IndexError:
        raise ValueError(f'{get_error_message("nofile")}')
    except AssertionError:
        raise FileNotFoundError(f'Input file does not exist: {ifn}')

    # NB: These are not tested by unit tests because tested separately
    verify_am2mbt_file(ifn, config=config)  # pragma: no cover
    process_am2mbt_file(ifn, config)  # pragma: no cover


if __name__ == '__main__':
    main()  # pragma: no cover
