# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 15:36:15 2020

@author: alessio
"""

from ftplib import FTP
import os
import patoolib
import warnings

warnings.filterwarnings("ignore") 

from Terminal_interface import *

main_folder=os.getcwd()

path_initial = os.path.join(main_folder,str(YEAR).zfill(4)+str(MONTH).zfill(2)+str(DOM).zfill(2)+str(SAT))
path_downloaded_data_LP=os.path.join(path_initial,'Downloaded_data','LP')
path_downloaded_data_TEC=os.path.join(path_initial,'Downloaded_data','TEC')

print('\n**********Swarm data download module**********')

##################################################################################
#Downloading Langmuir Probes data files from Swarm FTP
##################################################################################
print('Downloading data from Swarm FTP...')

try:
    for sat in SAT:
        ftp = FTP('swarm-diss.eo.esa.int')   # connect to host, default port
        ftp.login()               # user anonymous, passwd anonymous
        
        ftp.cwd('/Level1b/Latest_baselines/EFIx_LP/Sat_'+str(sat))
        
        listing=[]
        ftp.retrlines("LIST", listing.append)
        
        filenames=[]
        for index in listing:
            words=index.split(None, 8)
            if(words[-1].lstrip()[-4:]=='.ZIP'):
                filenames.append(words[-1].lstrip())
        
        flag=False
        for filename in filenames:
            if(filename[19:23]==str(YEAR).zfill(4) and filename[23:25]==str(MONTH).zfill(2) and filename[25:27]==str(DOM).zfill(2)):
                file_founded=filename
                flag=True
                break
        
    
        os.chdir(path_downloaded_data_LP)

        print('    Downloading LP file for Swarm '+str(sat)+' for '+str(YEAR).zfill(4)+'/'+str(MONTH).zfill(2)+'/'+str(DOM).zfill(2)+' ...')
        
        with open(file_founded, 'wb' ) as file :
                ftp.retrbinary('RETR %s' % file_founded, file.write)
                file.close()                

        if(not flag):
            raise Exception('    File not found for Swarm '+str(sat)+' for '+str(YEAR).zfill(4)+'/'+str(MONTH).zfill(2)+'/'+str(DOM).zfill(2))
            
        os.chdir(main_folder) 

    ftp.close()

except RuntimeError:
    raise Exception('Problems in connecting to or downloading from Swarm FTP')


##################################################################################
#Downloading TEC data files from Swarm FTP
##################################################################################

try:
    for sat in SAT:
        ftp = FTP('swarm-diss.eo.esa.int')   # connect to host, default port
        ftp.login()               # user anonymous, passwd anonymous
        
        ftp.cwd('/Level2daily/Latest_baselines/TEC/TMS/Sat_'+str(sat))
        
        listing=[]
        ftp.retrlines("LIST", listing.append)
        
        filenames=[]
        for index in listing:
            words=index.split(None, 8)
            if(words[-1].lstrip()[-4:]=='.ZIP'):
                filenames.append(words[-1].lstrip())
        
        flag=False
        for filename in filenames:
            if(filename[19:23]==str(YEAR).zfill(4) and filename[23:25]==str(MONTH).zfill(2) and filename[25:27]==str(DOM).zfill(2)):
                file_founded=filename
                flag=True
                break
        
    
        os.chdir(path_downloaded_data_TEC)

        print('    Downloading TEC file for Swarm '+str(sat)+' for '+str(YEAR).zfill(4)+'/'+str(MONTH).zfill(2)+'/'+str(DOM).zfill(2)+' ...')
        
        with open(file_founded, 'wb' ) as file :
                ftp.retrbinary('RETR %s' % file_founded, file.write)
                file.close()                

        if(not flag):
            raise Exception('    File not found for Swarm '+str(sat)+' for '+str(YEAR).zfill(4)+'/'+str(MONTH).zfill(2)+'/'+str(DOM).zfill(2))
            
        os.chdir(main_folder) 

    ftp.close()

except RuntimeError:
    raise Exception('Problems in connecting to or downloading from Swarm FTP')
        
##################################################################################
#Decompressing downloaded files 
##################################################################################
print('Decompression of downloaded files...')

#LP
os.chdir(path_downloaded_data_LP)
swarm_files=os.listdir(path_downloaded_data_LP)
swarm_files.sort()
for file in swarm_files:
    patoolib.extract_archive(file, outdir=path_downloaded_data_LP,verbosity=-1)
    os.remove(file)

files=os.listdir()
for file in files:
    if(file[-4:]=='.cdf'):
        pass
    else:
        os.remove(file)

os.chdir(main_folder)

#TEC
os.chdir(path_downloaded_data_TEC)
swarm_files=os.listdir(path_downloaded_data_TEC)
swarm_files.sort()
for file in swarm_files:
    patoolib.extract_archive(file, outdir=path_downloaded_data_TEC,verbosity=-1)
    os.remove(file)

files=os.listdir()
for file in files:
    if(file[-4:]=='.cdf'):
        pass
    else:
        os.remove(file)

os.chdir(main_folder)    
