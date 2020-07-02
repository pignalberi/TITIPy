# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 15:38:32 2020

@author: pigna
"""

import os
import logging

logging.basicConfig(filename='warning.log', level=logging.WARNING)


print('\n\n*************************   TITIPy   *************************')
print('******Topside Ionosphere Turbulence Indices with Python******\n')
print("""This Python module calculates RODI, ROTI, ROTEI and other 
physical parameters by downloading and analyzing electron 
density, electron temperature, and TEC data measured by 
Langmuir Probes and POD antennas on-board ESA Swarm satellites.

Developed by Alessio Pignalberi for ESA INTENS project
For questions please ask to alessio.pignalberi@ingv.it\n""")

###############################################################################
#Importing scripts
###############################################################################
from Terminal_interface import *

#creating tree folders
main_folder=os.getcwd()

path_initial = os.path.join(main_folder,str(YEAR).zfill(4)+str(MONTH).zfill(2)+str(DOM).zfill(2)+str(SAT))
try:
    os.mkdir(path_initial)
    os.makedirs(os.path.join(path_initial,'Downloaded_data','LP')) 
    os.makedirs(os.path.join(path_initial,'Downloaded_data','TEC')) 
    os.makedirs(os.path.join(path_initial,'Organized_data','LP')) 
    os.makedirs(os.path.join(path_initial,'Organized_data','TEC'))
    os.makedirs(os.path.join(path_initial,'Output','data','LP'))
    os.makedirs(os.path.join(path_initial,'Output','data','TEC'))
    os.makedirs(os.path.join(path_initial,'Output','figures','LP'))
    os.makedirs(os.path.join(path_initial,'Output','figures','TEC'))
except FileExistsError:
    print("Directory " , path_initial ,  " already exists")

path_downloaded_data=os.path.join(path_initial,'Downloaded_data')
path_organized_data=os.path.join(path_initial,'Organized_data')
path_output=os.path.join(path_initial,'Output')


import Downloading_Swarm_data

import Reading_Swarm_data_cdf

import Parameters_calculation_and_mapping

#writing run information file
filename_parameters=os.path.join(main_folder,'TITIPy_input_parameters.txt')
file=open(filename_parameters,'r')
lines=list(file) 
file.close()

filename_run='TITIPy_run_info_'+str(YEAR).zfill(4)+str(MONTH).zfill(2)+str(DOM).zfill(2)+str(SAT)+'.txt'
f=open(os.path.join(path_initial,filename_run),'w')
 
f.write("""***************************************************************************************************
**************************************TITIPy run info**********************************************
***************************************************************************************************\n""")

f.write('%4i\t\tYear\n%2i\t\tMonth\n%2i\t\tDay\n%1s\t\tSwarm satellite\n' % \
         (YEAR,MONTH,DOM,SAT))

for line in lines[3:6]:    
    f.write(line)

f.close()
    



print('**********************************************')
print("""Done\n  
Please find output files and figures (if wanted) in the
following folder: .../"""+str(YEAR).zfill(4)+str(MONTH).zfill(2)+str(DOM).zfill(2)+str(SAT)+"""/Output/ """)
print('\nThank you for using. Goodbye')
print('**********************************************')


