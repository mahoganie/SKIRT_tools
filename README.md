# SKIRT_tools
Tools to help streamline the process of running the SKIRT radiative transfer code on FIRE galaxies. Includes functions to convert data to a SKIRT-ready format, easily configure and run SKIRT simulations, and measure fluxes and luminosities from SKIRT outputs.

These tools can be used with any simulation as long as the data is loaded into the correct particle format. However, as the three-band images are layered using FIRE_studio, the mock images tend to turn out best with galaxies from the FIRE simulations.

These tools were inspired by Courtney Klein's mockobservation_tools: 

And they make use of FIRE studio, a code created by Alex Gurvich https://ascl.net/2202.006 \


## To install SKIRT and PTS:
SKIRT and PTS are required to use these tools. Please navigate to the SKIRT documentation and follow the steps to install them: https://skirt.ugent.be/root/_home.html

## To install FIRE_studio and SKIRT_tools:
Navigate to the directory where you want the packages to be installed and git clone the repositories:

git clone https://github.com/courtk32/FIRE_studio.git \
git clone https://github.com/mahoganie/SKIRT_tools.git

Then nativigate into the cloned directories and run:

pip install .


## How to cite:
Citation for FIRE_studio (Gurvich 2022): https://ui.adsabs.harvard.edu/abs/2022ascl.soft02006G/abstract \
Other citations: \
Starburst99: (Leitherer et al. 1999) https://iopscience.iop.org/article/10.1086/313233 \
(Weingartner & Draine 2001) https://iopscience.iop.org/article/10.1086/318651


## Available bands
The available bands and their center wavelengths in microns are:

Band / Central wavelength
                    
GALEX_FUV / 0.1535079263
GALEX_NUV / 0.2300784642

SDSS_u / 0.3556457276
SDSS_g / 0.4702444002
SDSS_r / 0.6175540451
SDSS_i / 0.748994523
SDSS_z / 0.8946683031

2MASS_J  / 1.239286141
2MASS_H  / 1.649435452
2MASS_Ks / 2.163824108

WISE_W1 / 3.389693501
WISE_W2 / 4.640627423
WISE_W3 / 12.56759007
WISE_W4 / 22.31422389

HERSCHEL_PACS70   / 70.76977781
HERSCHEL_PACS100  / 100.80266
HERSCHEL_PACS160  / 161.892196
HERSCHEL_SPIRE250 / 252.5495116
HERSCHEL_SPIRE350 / 354.2725745
HERSCHEL_SPIRE500 / 515.3633266
