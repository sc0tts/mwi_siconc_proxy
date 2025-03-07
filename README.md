# MWI Bootstrap Sea Ice Concentration (proxy)
Compute Bootstrap-algorithm sea ice concentration for MWI data.  This currently uses proxy data (AMSR2 swaths data from CLASS) until actual or sample MWI data is available to refine the algorithm.

## Installation

Note: Users can adapt these instructions as appropriate on their systems.  For instance, any of the following directories can be a symbolic link.  And with changes to the runtime scripts or configuration files, any of the specific directory names can be changed.


#Clone this repository:

  git clone git@github.com:sc0tts/mwi_siconc_proxy.git

#Create a directory with the ancillary data:

  Ancillary data is available from:
    https://www.dropbox.com/scl/fo/bpxya67aw2fw74hcc2n6d/AAog12InU1DpTE7f5uNkN2k?rlkey=zdt7mglbpzlowvddl8aaw9owg&st=7eda60d7&dl=0
    and is the directory:
       mwi_ancillary/
  Download this directory at the same level as mwi_siconc_proxy/

#Ensure that data are available:

  A day's worth of proxy data (NOAA swath data) are available at
    https://www.dropbox.com/scl/fo/bpxya67aw2fw74hcc2n6d/AAog12InU1DpTE7f5uNkN2k?rlkey=zdt7mglbpzlowvddl8aaw9owg&st=7eda60d7&dl=0
    in the directory:
      mwi_proxy_data/
  Download this directory at the same level as mwi_siconc_proxy/

#Copy (or create a link) to the runtime scripts in this directory:
  Script to generate Northern Hemisphere 20km EASE2 data:
    mwi_siconc_proxy/scripts/run_mwi_proxy_nh.sh
  Script to generate Southern Hemisphere 20km EASE2 data:
    mwi_siconc_proxy/scripts/run_mwi_proxy_sh.sh

#Verify that runtime directory is appropriate

Using the `tree` comment, the directory structure of this runtime environment will initially look like:

$ tree -L 1
.
├── mwi_ancillary
├── mwi_proxy_data
├── mwi_siconc_proxy
├── run_mwi_proxy_nh.sh
└── run_mwi_proxy_sh.sh


# Create the conda environment needed to run the code


