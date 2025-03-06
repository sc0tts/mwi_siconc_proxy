"""overlay_gridded.py

Overlay a series of gridded files
Later files in the list overwrite values placed on the grid from earlier files
"""

import os
import sys
import numpy as np
from pathlib import Path
from mwi_siconc_proxy.grid_am2mbt import (
    DEFAULT_CONFIG_AM2MBT,
    parse_swath_timestamp,
    parse_am2mbt_filename,
    get_intermed_var_filename,
    GRIDID_AREADEFS,
    read_am2mbt_config,
)
from mwi_siconc_proxy.grid_am2mbt import get_runtime_info as get_runinfo_swath


def get_runtime_info(run_label, config):
    # Return the information needed for the runtime
    runtime_info = {}
    runtime_info['gridid'] = config['gridid']
    runtime_info['prefix'] = run_label
    runtime_info['tb_list'] = ['h37', 'v19', 'v22', 'v37']
    runtime_info['outputdir'] = config['output_dir']
    runtime_info['tb_fn_template'] = config['tb_fn_template']

    return runtime_info


def overlay_gridded(listfile, run_label, swath_config):
    """Overlay the files from the files in the listfile"""

    print(f'Overlaying grids for: {run_label}')
    #swath_config = DEFAULT_CONFIG_AM2MBT

    # TODO: need to add gridid to intermediate file name
    runtime_info = get_runtime_info(run_label, config)
    gridid = runtime_info['gridid']
    grid_areadef = GRIDID_AREADEFS[gridid]

    xdim = grid_areadef['width']
    ydim = grid_areadef['height']

    tb_grids = {}
    for tb in runtime_info['tb_list']:
        tb_grids[tb] = np.zeros((ydim, xdim), dtype=np.float32)
        tb_grids[tb][:] = np.nan

    with open(listfile, 'r') as f:
        for line in f.readlines():
            fn = line.rstrip()
            print(f'file: {fn}')
            fn_parts = parse_am2mbt_filename(fn)
            swath_info = get_runinfo_swath(fn, swath_config)

            for tb in runtime_info['tb_list']:
                print(f'tb: {tb}')
                intermed_tb_fn = get_intermed_var_filename(tb, swath_info)
                gridded_swath = np.fromfile(intermed_tb_fn, dtype=np.float32).reshape(ydim, xdim)
                is_valid = np.isfinite(gridded_swath)
                tb_grids[tb][is_valid] = gridded_swath[is_valid]

    # Here, the tb_grids structure has a filled tb field for each tb
    for tb in runtime_info['tb_list']:

        fn_template = runtime_info['tb_fn_template']
        fn = fn_template.format(
            outputdir=runtime_info['outputdir'],
            prefix=runtime_info['prefix'],
            gridid=gridid,
            tb=tb,
        )
        tb_grids[tb].tofile(fn)
        print(f'Wrote: {fn}')

    print(f'Finished with overlay_gridded')


if __name__ == '__main__':
    try:
        swath_file_list = sys.argv[1]
        swath_list_path = Path(swath_file_list)
        #assert os.path.isfile(swath_file_list)
        assert swath_list_path.is_file()
    except IndexError:
        raise ValueError('No swath file list given')
    except AssertionError:
        raise ValueError(f'No such swath file list: {swath_file_list}')

    # Check for a second argument with name of YAML with config info
    try:
        config_fn = sys.argv[2]
        assert os.path.isfile(config_fn)
        config = read_am2mbt_config(config_fn)
    except IndexError:
        # No config file specified
        #config = DEFAULT_CONFIG_AM2MBT
        raise RuntimeError(f'No config file specified')
    except AssertionError:
        raise ValueError(f'No such specified config file: {cfn}')
    except ValueError:
        raise ValueError(f'Could not parse specified config file: {cfn}')

    run_label = swath_list_path.stem
    overlay_gridded(swath_file_list, run_label, config)
