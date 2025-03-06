""" set_siconc_parameters.py

Provide parameters needed for sea ice concentration algorithms

Sample usage:

     python mwi_siconc_proxy/src/mwi_siconc_proxy/set_siconc_parameters.py nam2sic_nh btparams_nam2sic_nh.yaml
"""

import sys
from pathlib import Path
import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from pm_icecon.config.models.bt import (
    BootstrapParams,
    WeatherFilterParams,
    WeatherFilterParamsForSeason,
    TbSetParams,
)


NOAA_AMSR2_BT_PARAMS = {
    'nam2sic_nh': BootstrapParams(
        add1=0.0,
        add2=-2.0,
        minic=10.0,
        maxic=1.0,
        mintb=10.0,
        maxtb=320.0,
        weather_filter_seasons=[
            WeatherFilterParamsForSeason(
                start_month=6,
                start_day=None,
                end_month=10,
                end_day=15,
                weather_filter_params=WeatherFilterParams(
                    wintrc=89.2,
                    wslope=0.50375,
                    wxlimt=21.0
                )
            ),
            WeatherFilterParamsForSeason(
                start_month=10,
                start_day=16,
                end_month=5,
                end_day=None,
                weather_filter_params=WeatherFilterParams(
                    wintrc=87.6467,
                    wslope=0.517333,
                    wxlimt=14.0
                )
            ),
        ],
        vh37_params=TbSetParams(
            water_tie_point_set=(201.916, 132.815),
            ice_tie_point_set=(255.67, 241.713),
            lnline={'offset': -73.5471, 'slope': 1.21104}
        ),
                v1937_params=TbSetParams(
            water_tie_point_set=(201.916, 178.771),
            ice_tie_point_set=(255.67, 258.341),
            lnline={'offset': 47.0061, 'slope': 0.809335}
        ),
    ),
}

def get_bootstrap_parameters(bt_parameter_label, param_set=NOAA_AMSR2_BT_PARAMS):
    """Return the BootstrapParams for this bt parameter set"""
    try:
        params = param_set[bt_parameter_label]
    except KeyError:
        raise ValueError(f'Missing Bootstrap parameters:\n  {bt_parameter_label} is not in {param_set}:\n    {param_set.keys()}')

    return params


def write_btparams_file(params, path):
    # Unfortunately, can't use yaml.safe_dump()
    #   because BootstrapParams uses non-standard structures
    with open(path, 'w') as f:
        yaml.dump(params, f)


def read_btparams_file(path):
    """Read a BootstrapParams structure from the path"""
    with open(path, 'r') as f:
        params = yaml.load(f, Loader=Loader)

    return params


def main_set_siconc_parameters():
    # Save the parameters listed here to a yaml file
    idx = sys.argv[1]
    fn = sys.argv[2]
    path = Path(fn)
    if path.is_file():
        raise RuntimeError(f'output file exists, quitting: {fn}')

    assert fn[-4:] == '.yml' or fn[-5:] == '.yaml'

    params = NOAA_AMSR2_BT_PARAMS[idx]

    write_btparams_file(params, path)

    # Check that it's the same
    testparams = read_btparams_file(path)
    assert params == testparams


if __name__ == '__main__':
    main_set_siconc_parameters()
