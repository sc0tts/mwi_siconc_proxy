# Tests for test_grid_am2mbt.py

import datetime as dt
import pytest
from netCDF4 import Dataset
import yaml
import numpy as np
from pyresample.geometry import AreaDefinition, SwathDefinition
import warnings

from mwi_siconc_proxy import grid_am2mbt
from mwi_siconc_proxy.grid_am2mbt import (
    EXAMPLE_AM2MBT_FNAME,
    DEFINED_AM2MBT_ERRORS,
    DEFAULT_CONFIG_AM2MBT,
    EXAMPLE_CONFIG_FNAME,
    get_usage_string,
    get_error_message,
    parse_swath_timestamp,
    parse_am2mbt_filename,
    get_filename_timestamp,
    verify_am2mbt_dataset,
    verify_am2mbt_file,
    read_am2mbt_config,
    main,
    get_pyresample_area_def,
    get_pyresample_swathdef,
    GRIDID_AREADEFS,
    grid_swath_data,
    get_runtime_info,
)


EXPECTED_DATETIME = dt.datetime(2025, 1, 2, 20, 30, 40, 500000)
YmdHMSf_STR = EXPECTED_DATETIME.strftime('%Y%m%d%H%M%S%f')[:15]

LATLON_AREADEF = {
    'area_id': 'sample_latlon_grid',
    'description': 'Simple_latlon_grid_at_5_deg',
    'proj_id': 'sample_latlon',
    'projection': 'EPSG:4326',
    'width': 72,
    'height': 18,
    'area_extent': (-180, 90, 180, -90),
},



"""Note: an AMSR2 MBT swath file with only the variables needed for this
gridding would have headers like:

netcdf reduced {
dimensions:
        Number_of_Scans = UNLIMITED ; // (3960 currently)
        Number_of_low_rez_FOVs = 243 ;
        Number_of_low_rez_channels = 6 ;
        Number_of_hi_rez_FOVs = 486 ;
variables:
        float Brightness_Temperature_18_GHzH(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Brightness_Temperature_18_GHzV(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Brightness_Temperature_23_GHzV(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Brightness_Temperature_36_GHzH(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Brightness_Temperature_36_GHzV(Number_of_Scans, Number_of_low_rez_FOVs) ;
        short Land_Ocean_Flag_6_to_36(Number_of_low_rez_channels, Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Latitude_for_23(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Latitude_for_36(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Latitude_for_Low_Resolution(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Longitude_for_18(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Longitude_for_23(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Longitude_for_36(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Longitude_for_Low_Resolution(Number_of_Scans, Number_of_low_rez_FOVs) ;
        short Pixel_Data_Quality_6_to_36(Number_of_Scans, Number_of_hi_rez_FOVs) ;
        float RFI_NOAA_OCEAN_18h(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float RFI_NOAA_OCEAN_18v(Number_of_Scans, Number_of_low_rez_FOVs) ;
        char quality_information ;

        float Brightness_Temperature_18_GHzH(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Brightness_Temperature_18_GHzV(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Brightness_Temperature_23_GHzV(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Brightness_Temperature_36_GHzH(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Brightness_Temperature_36_GHzV(Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Latitude_for_23               (Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Latitude_for_36               (Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Latitude_for_Low_Resolution   (Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Longitude_for_18              (Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Longitude_for_23              (Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Longitude_for_36              (Number_of_Scans, Number_of_low_rez_FOVs) ;
        float Longitude_for_Low_Resolution  (Number_of_Scans, Number_of_low_rez_FOVs) ;
        float RFI_NOAA_OCEAN_18h            (Number_of_Scans, Number_of_low_rez_FOVs) ;
        float RFI_NOAA_OCEAN_18v 
"""

def gen_temp_file_ref(dpath, fname):
    # Creates a Path (file) from a PathDir and a filename
    dpath.mkdir()
    fpath = dpath / fname

    return fpath


def gen_dummy_am2mbt_dataset(path):
    """ Create a dummy am2mbt_dataset/file for testing
    The path needs to be writeable.  Can be a tmp_path"""
    with Dataset(path, 'w') as ds:
        # Create the dimensions...
        ds.createDimension("Number_of_Scans", None)  # Unlimited
        ds.createDimension("Number_of_low_rez_FOVs", 2)
        ds.createDimension("Number_of_low_rez_channels", 6)
        ds.createDimension("Number_of_hi_rez_FOVs", 4)  # only needed for Pix...

        # Create the variables with dims nscans, nlowrez_fovs
        swath_vars = {
            'Brightness_Temperature_18_GHzH': ('Number_of_Scans', 'Number_of_low_rez_FOVs'),
            'Brightness_Temperature_18_GHzV': ('Number_of_Scans', 'Number_of_low_rez_FOVs'),
            'Brightness_Temperature_23_GHzV': ('Number_of_Scans', 'Number_of_low_rez_FOVs'),
            'Brightness_Temperature_36_GHzH': ('Number_of_Scans', 'Number_of_low_rez_FOVs'),
            'Brightness_Temperature_36_GHzV': ('Number_of_Scans', 'Number_of_low_rez_FOVs'),
            'Latitude_for_18': ('Number_of_Scans', 'Number_of_low_rez_FOVs'),               
            'Latitude_for_23': ('Number_of_Scans', 'Number_of_low_rez_FOVs'),               
            'Latitude_for_36': ('Number_of_Scans', 'Number_of_low_rez_FOVs'),               
            'Latitude_for_Low_Resolution': ('Number_of_Scans', 'Number_of_low_rez_FOVs'),   
            'Longitude_for_18': ('Number_of_Scans', 'Number_of_low_rez_FOVs'),              
            'Longitude_for_23': ('Number_of_Scans', 'Number_of_low_rez_FOVs'),              
            'Longitude_for_36': ('Number_of_Scans', 'Number_of_low_rez_FOVs'),              
            'Longitude_for_Low_Resolution': ('Number_of_Scans', 'Number_of_low_rez_FOVs'),  
            'RFI_NOAA_OCEAN_18h': ('Number_of_Scans', 'Number_of_low_rez_FOVs'),            
            'RFI_NOAA_OCEAN_18v': ('Number_of_Scans', 'Number_of_low_rez_FOVs'), 
        }
        for swath_var_name in swath_vars:
            swath_var_dtype = "f4"
            nscans, nfovs = swath_vars[swath_var_name]
            var_ref = ds.createVariable(swath_var_name, swath_var_dtype, (nscans, nfovs))

        # Create the other variables individually
        var_name = 'Pixel_Data_Quality_6_to_36'
        var_dtype= 'i2'
        var_ref = ds.createVariable(
            var_name,
            var_dtype,
            ('Number_of_Scans', 'Number_of_hi_rez_FOVs'),
        )

        var_name = 'Land_Ocean_Flag_6_to_36'
        var_dtype= 'i2'
        var_ref = ds.createVariable(
            var_name,
            var_dtype,
            ('Number_of_low_rez_channels', 'Number_of_Scans', 'Number_of_low_rez_FOVs'),
        )

        var_name = 'quality_information'
        var_dtype= 'S1'
        var_ref = ds.createVariable(
            var_name,
            var_dtype,
            (),
        )

    # End of context for writing this netCDF file


def test_default_config_has_keys():
    config = DEFAULT_CONFIG_AM2MBT
    expected_keys = ("nc_vars", "grid_vars", "var_fn_template")
    for key in expected_keys:
        assert key in config.keys()


def test_usage_string():
    # Verify that usage string is returned
    usage_string = get_usage_string()
    assert 'Usage' in usage_string
    

def test_error_messages():
    # Verify each error message
    for error in DEFINED_AM2MBT_ERRORS:
        error_message = get_error_message(error)
        assert error_message is not None
        assert error_message != ''


def test_parse_swath_timestamp_types():
    # Ensure that start, end, and creation values can be computed
    for leading_char, ts_type in zip(('s', 'e', 'c'), ('start', 'end', 'creation')):
        test_ts_string = f'{leading_char}{YmdHMSf_STR}'
        ts_calculated = parse_swath_timestamp(test_ts_string, ts_type)
        assert ts_calculated == EXPECTED_DATETIME

def test_parse_swath_timestamp_bad_type():
    # Verify that a bad timestamp type raises ValueError
    leading_char = 'Z'
    test_ts_string = f'{leading_char}{YmdHMSf_STR}'
    with pytest.raises(ValueError):
        parse_swath_timestamp(test_ts_string, 'start')


def test_parse_swath_timestamp_unknown():
    # Verify that a bad timestamp type raises ValueError
    leading_char = 'Z'
    test_ts_string = f'{leading_char}{YmdHMSf_STR}'

    # This should be work without error
    calculated_timestamp = parse_swath_timestamp(test_ts_string, 'unknown')
    assert calculated_timestamp == EXPECTED_DATETIME

    # Partial string should raise ValueError
    with pytest.raises(ValueError):
        calculated_timestamp = parse_swath_timestamp(test_ts_string[:10], 'unknown')


def test_get_filename_timestamp():
    # Verify extraction of timestamp from typical AMSR2 swath file name
    # Note: EXAMPLE_AM2MBT_FNAME is defined in grid_am2mbt.py as:
    #  AMSR2-MBT_v2r2_GW1_s202501060005067_e202501060144046_c202501060156059.nc
    calc_timestamp = get_filename_timestamp(EXAMPLE_AM2MBT_FNAME)
    expected_timestamp = dt.datetime(2025, 1, 6, 0, 5, 6, 700000)

    assert expected_timestamp == calc_timestamp


def test_parse_swath_timestamp_bad_timestamp():
    # Verify bad timestamp fails
    # Note that this bad timestamp begins with 's' ('e' or 'c' would also work)
    bad_timestamp = 's8adk43k'
    with pytest.raises(ValueError):
        parse_swath_timestamp(bad_timestamp)


def test_verify_am2mbt_file_fn_must_exist():
    # Ensure that filename must exist if given
    nonexistent_filename = 'thisisnotafile.nc'
    with pytest.raises(OSError):
        verify_am2mbt_file(nonexistent_filename)


def test_verify_am2mbt_file_sample_fn(tmp_path):
    # Ensure that sample netCDF file passes verify() routine
    dpath = tmp_path / "ncdir"
    dummy_am2mbt_ncname = EXAMPLE_AM2MBT_FNAME
    dummy_am2mbt_path = gen_temp_file_ref(dpath, dummy_am2mbt_ncname)
    gen_dummy_am2mbt_dataset(dummy_am2mbt_path)

    verify_am2mbt_file(dummy_am2mbt_path)


def test_verified_fn_detects_nc_suffix(tmp_path):
    # If a file ends in '.nc' then is netCDF file
    dpath = tmp_path / "ncdir"
    dummy_am2mbt_ncname = EXAMPLE_AM2MBT_FNAME
    filename_parts = parse_am2mbt_filename(dummy_am2mbt_ncname)
    assert dummy_am2mbt_ncname[-3:] == '.nc'
    assert filename_parts['is_ncfile']


def test_verified_fn_detects_not_nc_suffix(tmp_path):
    # If a file doesn't ends in '.nc' then should still get some parts
    dpath = tmp_path / "ncdir"
    dummy_am2mbt_ncname = EXAMPLE_AM2MBT_FNAME
    dummy_am2mbt_ncname = dummy_am2mbt_ncname.replace('.nc', '.notnc')
    filename_parts = parse_am2mbt_filename(dummy_am2mbt_ncname)
    #breakpoint()
    assert dummy_am2mbt_ncname[-3:] != '.nc'
    assert not filename_parts['is_ncfile']
    assert not filename_parts['is_standard']


def test_verify_am2mbt_file_fn_must_yield_timestamp(tmp_path):
    # Ensure that filename must exist if given
    dpath = tmp_path / "ncdir"
    bad_timestamp_filename = EXAMPLE_AM2MBT_FNAME.replace('_s', '_x')
    bad_timestamp_path = gen_temp_file_ref(dpath, bad_timestamp_filename)
    # This causes bad_timestamp_path to be a "valid" netcdf file
    gen_dummy_am2mbt_dataset(bad_timestamp_path)

    default_timestamp = dt.datetime(2000, 2, 2, 2, 2, 2, 200000)
    verify_am2mbt_file(bad_timestamp_path, default_timestamp=default_timestamp)


def test_verify_am2mbt_file_must_be_ncfile(tmp_path):
    dpath = tmp_path / "badnc_dir"
    dpath.mkdir()
    badnc_path = dpath / EXAMPLE_AM2MBT_FNAME
    badnc_path.write_text('this is not a valid netCDF file!', encoding='utf-8')

    with pytest.raises(OSError):
        verify_am2mbt_file(badnc_path)


def test_verify_am2mbt_dataset_min_config(tmp_path):
    dpath = tmp_path / "ncdir"
    dummy_nc_path = gen_temp_file_ref(dpath, EXAMPLE_AM2MBT_FNAME)
    gen_dummy_am2mbt_dataset(dummy_nc_path)
    ds = Dataset(dummy_nc_path, 'r')

    min_config = {
        'nc_vars': {
            'v18': 'Brightness_Temperature_18_GHzH',
            'lat18': 'Latitude_for_18',
            'lon18': 'Longitude_for_18',
        },
        'grid_vars': {
            'v19': ['v18', 'lat18', 'lon18'],
        },
    }
    verify_am2mbt_dataset(ds, min_config)


def test_verify_am2mbt_dataset_must_have_consistent_vars(tmp_path):
    dpath = tmp_path / "ncdir"
    dummy_nc_path = gen_temp_file_ref(dpath, EXAMPLE_AM2MBT_FNAME)
    gen_dummy_am2mbt_dataset(dummy_nc_path)
    ds = Dataset(dummy_nc_path, 'r')

    config_without_contrib_var = {
        'nc_vars': {
            'v18': 'Brightness_Temperature_18_GHzH',
            'lat18': 'Latitude_for_18',
            #'lon18': 'Longitude_for_18',
        },
        'grid_vars': {
            'v19': ['v18', 'lat18', 'lon18'],
        },
    }
    with pytest.raises(ValueError):
        verify_am2mbt_dataset(ds, config_without_contrib_var)


def test_verify_am2mbt_dataset_unexpected_nc_var(tmp_path):
    dpath = tmp_path / "ncdir"
    dummy_nc_path = gen_temp_file_ref(dpath, EXAMPLE_AM2MBT_FNAME)
    gen_dummy_am2mbt_dataset(dummy_nc_path)
    ds = Dataset(dummy_nc_path, 'r')

    config_with_unexpected_nc_var = {
        'nc_vars': {
            'v18': 'Brightness_Temperature_18_GHzH',
            'lat18': 'Latitude_for_18',
            'lon18': 'Longitude_for_18',
            'missing18': 'data_field_not_in_netCDF_file',
        },
        'grid_vars': {
            'v19': ['v18', 'lat18', 'lon18'],
        },
    }
    with pytest.raises(ValueError):
        verify_am2mbt_dataset(ds, config_with_unexpected_nc_var)


def test_sample_config_file(tmp_path):
    # Test the use of a sconfiguration file by writing the default config to tmp file
    dpath = tmp_path / 'config'
    dpath.mkdir()
    config_path = dpath / EXAMPLE_CONFIG_FNAME

    with open(config_path, 'w') as ymlfile:
        yaml.safe_dump(DEFAULT_CONFIG_AM2MBT, ymlfile)

    with open(config_path, 'r') as configfile:
        config = yaml.safe_load(configfile)

    # Assure that the tmpfile writing worked
    assert config == DEFAULT_CONFIG_AM2MBT

    read_config = read_am2mbt_config(config_path)
    assert read_config == DEFAULT_CONFIG_AM2MBT


def test_read_am2mbt_config_nofile(tmp_path):
    # Test the use of a sconfiguration file by writing the default config to tmp file
    dpath = tmp_path / 'config'
    dpath.mkdir()
    config_path = dpath / 'nosuchfile'

    with pytest.raises(FileNotFoundError):
        read_config = read_am2mbt_config(config_path)


def test_read_am2mbt_config_not_yaml_file(tmp_path):
    # Test the use of a sconfiguration file by writing the default config to tmp file
    dpath = tmp_path / 'config'
    dpath.mkdir()
    config_path = dpath / 'notayamlfile'
    config_path.write_text('*&#*(#-not a YAML file', encoding='utf-8')

    with pytest.raises(ValueError):
        read_config = read_am2mbt_config(config_path)


def test_read_am2mbt_config_empty_yaml_file(tmp_path):
    # Test the use of a sconfiguration file by writing the default config to tmp file
    dpath = tmp_path / 'config'
    dpath.mkdir()
    config_path = dpath / 'nosuchfile'
    config_path.write_text('', encoding='utf-8')

    with pytest.raises(ValueError):
        read_config = read_am2mbt_config(config_path)


def test_main_args_exist(monkeypatch, tmp_path):
    """Test erroneous cmdline args calling as __main__()"""
    import sys

    with monkeypatch.context() as m:
        # No arg1 means no file was specified
        m.setattr(sys, 'argv', ['grid_am2mbt',])
        with pytest.raises(ValueError) as error_no_arg1:
            main_output = main()
        assert 'No input file provided' in str(error_no_arg1.value)

        # No file (arg1) found is an error
        m.setattr(sys, 'argv', ['grid_am2mbt', 'isnotafile'])
        with pytest.raises(FileNotFoundError) as error_no_such_arg1:
            main_output = main()
        assert 'Input file does not exist' in str(error_no_such_arg1.value)

    # Create temporary file
    dpath = tmp_path / 'run'
    dpath.mkdir()
    fn_path = dpath / 'sample.nc'
    fn_path.write_text('', encoding='utf-8')

    bad_config_path = dpath / 'config.yaml'
    bad_config_path.write_text('', encoding='utf-8')

    with monkeypatch.context() as m:
        # No arg1 means no file was specified
        m.setattr(sys, 'argv', ['grid_am2mbt', fn_path, 'not_a_config_file'])
        with pytest.raises(ValueError) as error_no_such_arg2:
            main_output = main()
        assert 'No such specified config file' in str(error_no_such_arg2.value)

        # arg2 that is not valid config file is error
        m.setattr(sys, 'argv', ['grid_am2mbt', fn_path, bad_config_path])
        with pytest.raises(ValueError) as error_bad_config_path:
            main_output = main()
        assert 'Could not parse' in str(error_bad_config_path.value)


def test_get_get_pyresample_area_def_requires_valid_gridid():
    with pytest.raises(ValueError):
        get_pyresample_area_def('bad_gridid')


def test_get_get_pyresample_area_def_has_gridid_information():
    for gridid in GRIDID_AREADEFS.keys():
        area_def = get_pyresample_area_def(gridid)
        assert type(area_def) == AreaDefinition


def test_get_get_pyresample_area_def_notices_missing_gridid():
    gridid = 'not_a_valid_gridid'
    with pytest.raises(ValueError) as error_bad_config_path:
        area_def = get_pyresample_area_def(gridid)


def test_get_pyresample_swathdef():
    # Ensure that we can define a pyresample SwathDefinition
    lats = np.ones((3, 2))
    lons = np.ones((3, 2))
    swath_def = get_pyresample_swathdef(lons=lons, lats=lats)
    assert type(swath_def) == SwathDefinition


def test_get_pyresample_swathdef():
    # Ensure that we get a valid SwathDefinition from lons and lats
    lat_vals = np.linspace(80, 20, num=7, endpoint=True)
    lon_vals = np.linspace(-10, 10, num=3, endpoint=True)
    lons, lats = np.meshgrid(lon_vals, lat_vals)

    swathdef = get_pyresample_swathdef(lons, lats)
    assert isinstance(swathdef, SwathDefinition)

def test_grid_swath_data():
    # Ensure we can grid a simple swath

    # Define 10-deg separated swath
    lat_vals = np.linspace(80, 20, num=7, endpoint=True)
    lon_vals = np.linspace(-10, 10, num=3, endpoint=True)
    lons, lats = np.meshgrid(lon_vals, lat_vals)
    
    # simple values are lat+lon
    #vals = lats
    vals = lons
    #vals = lats + lons

    #areadef = AreaDefinition('latlon_5deg')
    areadef = get_pyresample_area_def('latlon_5deg')

    # There are pyproj UserWarning and numpy resampler warnings to ignore
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # Default radius of 20e3 yields all zeros
        #gridded = grid_swath_data(vals, lons, lats, areadef)
        gridded = grid_swath_data(vals, lons, lats, areadef, radius=20e7)

    assert isinstance(gridded, np.ndarray)


def test_get_runtime_info():
    # Verify we get expected info from typical input file and config
    test_config = DEFAULT_CONFIG_AM2MBT
    test_am2mbt_ncname = EXAMPLE_AM2MBT_FNAME

    runtime_info = get_runtime_info(test_am2mbt_ncname, test_config)

    assert 'timestamp_str' in runtime_info.keys()
    assert 'tb_source_str' in runtime_info.keys()
