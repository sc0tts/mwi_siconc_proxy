# MWI Bootstrap Sea Ice Concentration (proxy)
Compute Bootstrap-algorithm sea ice concentration for MWI data.  This currently uses proxy data (AMSR2 swaths data from CLASS) until actual or sample MWI data is available to refine the algorithm.

# Installation

Note: Users can adapt these instructions as appropriate on their systems.  For instance, any of the following directories can be a symbolic link.  And with changes to the runtime scripts or configuration files, any of the specific directory names can be changed.


## Clone this repository:

  `git clone git@github.com:sc0tts/mwi_siconc_proxy.git`

## Create a directory with the ancillary data:

  Ancillary data is available from:

    `https://www.dropbox.com/scl/fo/bpxya67aw2fw74hcc2n6d/AAog12InU1DpTE7f5uNkN2k?rlkey=zdt7mglbpzlowvddl8aaw9owg&st=7eda60d7&dl=0`

    and is the directory:

       `mwi_ancillary/`

  Download this directory at the same level as `mwi_siconc_proxy/`

## Ensure that data are available:

  A day's worth of proxy data (NOAA swath data) are available at:

    `https://www.dropbox.com/scl/fo/bpxya67aw2fw74hcc2n6d/AAog12InU1DpTE7f5uNkN2k?rlkey=zdt7mglbpzlowvddl8aaw9owg&st=7eda60d7&dl=0`

    in the directory:

      `mwi_proxy_data/`

  Download this directory at the same level as `mwi_siconc_proxy/`

## Copy (or create a link) to the runtime scripts in this directory:

  Script to generate Northern Hemisphere 20km EASE2 data:

    `mwi_siconc_proxy/scripts/run_mwi_proxy_nh.sh`

  Script to generate Southern Hemisphere 20km EASE2 data:

    `mwi_siconc_proxy/scripts/run_mwi_proxy_sh.sh`

## Verify that runtime directory is appropriate

Using the `tree` comment, the directory structure of this runtime environment will initially look like:

```
$ tree -L 1
.
├── mwi_ancillary
├── mwi_proxy_data
├── mwi_siconc_proxy
├── run_mwi_proxy_nh.sh
└── run_mwi_proxy_sh.sh
```

Note: ensure that the execute-permission bit is set for the .sh scripts.

## Create a Python environment needed to run the code

There are different methods for installing a python package into a working environment.  The following is one way; it involves creating a conda environment and then installing this code as a local package in editable mode.  Other methods will also work.

An environment.yml file provides the recipe for setting up a conda environment that will run this code.  Other methods of incorporating this package into a Python environment will also work.

To create a conda environment named `mwi_siconc_proxy`, the environment.yml file can be used:

  `conda env create --file=./mwi_siconc_proxy/environment.yml`

Note: Some users prefer the command `mamba` intead of `conda.  Either will work.

This will create a conda environment named: `mwi_siconc_proxy`.  (See the top of the `environment.yml` file.)

To install this repository as an editable module, the method we have found to be least problematic is:
  - Activate the conda environment:
    - `conda activate mwi_siconc_proxy`
  - Change directory to the root of this code repository, i.e.:
    - `cd ./mwi_siconc_proxy/`
  - From the code repository directory, install this module using pip:
    - `python -m pip install --no-build-isolation --no-deps -e .`
    - Note that `.` at the end of that command; it indicates the directory in which the command is run

## Run the sample scripts

If the code, ancillary, and proxy data are installed as described above, two sample scripts can be run which will generate netCDF files of sea ice concentration; there is one script demonstrating generation for the Northern Hemisphere, and one script for the Southern Hemisphere.

Note that in general, three files are needed for a run:
  - config: A configuration file specifying runtime values such as intermediate and output directories as well as filename templates.  This file also indicates how data variables in the swath input files are mapped to the brightness temperature variables needed for the sea ice concentration algorithm
  - parameter file:  A file giving the parameters for the Bootstrap sea ice concentration algorithm
  - swathfile list:  A file containing a list of input swath files.  Currently, the name of this swathfile is used to name the final sea ice concentration output file.

Note: Intermediate brightness temperature files are named based on the grid, variable name (ie TB channel), and timestamp.  This allows the same file to be used without re-generation if multiple runs would use the same input file(s).

Note: intermediate files are not currently removed by these scripts.  This allows files to be re-used by later script runs.

Generate the Northern Hemisphere sea ice concentration:

  `./run_mwi_proxy_nh.sh`

Generate the Southern Hemisphere sea ice concentration:

  `./run_mwi_proxy_sh.sh`
