###########################################################################################################################################################################################################
TITIPy (Topside Ionosphere Turbulence Indices with Python)
###########################################################################################################################################################################################################
This Python tool calculates RODI, ROTI, ROTEI and other physical parameters by downloading and analyzing electron density, electron temperature, and TEC data measured by Langmuir Probes and POD antennas on-board ESA Swarm satellites.

Developed by Alessio Pignalberi for ESA INTENS project. For questions please ask to alessio.pignalberi@ingv.it
###########################################################################################################################################################################################################
Requirements:
OS=Linux
Language=Python 3.7.6+

Python libraries needed (not built-in): 
numpy 1.18.1+
cdflib 0.3.18+
basemap 1.2.0+
basemap-data-hires 1.2.0+
apexpy 1.0.3+
patool 1.12+
spacepy 0.2.1+

It is suggested to install Python through Anaconda (https://www.anaconda.com/products/individual). 
After that, Python needed libraries can be installed with the following commands from terminal:
conda install -c conda-forge basemap
conda install basemap-data-hires
pip install apexpy
pip install patool
pip install cdflib
pip install spacepy
###########################################################################################################################################################################################################
TITIPy folder contains the folowing files:
README.md (this file), containing useful information for installation and usage of this tool.
LICENSE, the license file.

Main.py, this is the script to launch a TITIPy run. It imports and runs the other scripts of the tool.
Terminal_interface.py, it contains the terminal interface to allow the user to input the year, month, day, and Swarm satellite.
Downloading_Swarm_data.py, this script downloads and decompresses Langmuir Probes and TEC files (cdf format) from ESA Swarm FTP. Downloaded files are put in the folder ../Downloaded_data/.
Reading_Swarm_data_cdf.py, this script reads the downloaded cdf files, calculates several quantities, and saves data in simple txt files. Data are put in the folder ../Organized_data/.
Parameters_calculation_and_mapping.py, this script calculates ionospheric indices (RODI, ROTEI, ROTI) and save the results in txt files. If activated, corresponding scatter plot maps and polar magnetic plots are created. Data and figures are put in the folder ../Output/.
Functions.py, it contains several functions used for indices calculation and binning.

TITIPY_input_parameters.txt, through this file some parameters, other than those defined through the terminal interface, can be defined.
logo.png, TITIPy logo figure image.
User_credentials.txt, this file contains the username and password to access the ESA Swarm dissemination center. Please replace USERNAME and PASSWORD with your own credentials. To obtain credentials see https://earth.esa.int/eogateway/faq/how-do-i-access-swarm-data
###########################################################################################################################################################################################################
TITPy tool can be run from terminal simply by navigating to the corresponding folder tool and by typing the following command:
python Main.py

Alternatively, any Python IDLE (Integrated Development and Learning Environment), like Spyder or Jupyter, can be used to run the script Main.py

Results of the TITIPy run, data and figures (if wanted), are put in the folder (created during the run) named YYYYMMDDS, where YYYY=year, MM=month, DD=day of the month, and S=Swarm satellite, chosen through the terminal interface.

For further info on TITIPy, please refers to the publication: Pignalberi, A. (2021) TITIPy: A Python tool for the calculation and mapping of topside ionosphere turbulence indices. Computers & Geosciences,148, 104675. doi: 10.1016/j.cageo.2020.104675
