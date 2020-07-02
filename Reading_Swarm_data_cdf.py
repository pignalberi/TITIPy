# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 11:03:02 2018

@author: alessio
"""

import numpy as np
import os
import datetime
import cdflib
from apexpy import Apex
import warnings

warnings.filterwarnings("ignore") 

from Terminal_interface import *

main_folder=os.getcwd()

path_initial = os.path.join(main_folder,str(YEAR).zfill(4)+str(MONTH).zfill(2)+str(DOM).zfill(2)+str(SAT))
path_downloaded_data_LP=os.path.join(path_initial,'Downloaded_data','LP')
path_downloaded_data_TEC=os.path.join(path_initial,'Downloaded_data','TEC')

path_organized_data_LP=os.path.join(path_initial,'Organized_data','LP')
path_organized_data_TEC=os.path.join(path_initial,'Organized_data','TEC')


MEAN_EARTH_RADIUS=6371.007 #mean Earth radius IUGG

print('Reading and converting downloaded data...')

###############################################################################
#Working on Langmuir Probes data
###############################################################################

filename_swarm=os.listdir(path_downloaded_data_LP)
filename_swarm.sort()

file=filename_swarm[0]
for file in filename_swarm:
    
    SAT=file[11:12]
    YEAR=file[19:23]
    MONTH=file[23:25]
    DOM=file[25:27]
    
    filename_swarm_out='Swarm_LP_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DOM+'_data.txt'
    f=open(os.path.join(path_organized_data_LP,filename_swarm_out),'w')
    f.write('year	      month	       day	    	hour		min	    	sec	    	msec		doy	         hourUT		  hourLT		     MLT		    LatGeo		    LonGeo		Radius [m]	    	Height [km]	      	  LatMag_QD	      	  LonMag_QD		Ne [cm^-3]	  	Te [K]		Flag_LP      Flag_Ne	     Flag_Te\n')
       
    cdf = cdflib.CDF(os.path.join(path_downloaded_data_LP,file))
        
    #reading variables
    Timestamp=cdf.varget('Timestamp')   
    Timestamp=cdflib.cdfepoch.breakdown(Timestamp)
    Latitude=cdf.varget('Latitude')
    Longitude=cdf.varget('Longitude')
    Radius=cdf.varget('Radius')
    Ne=cdf.varget('Ne')
    Te=cdf.varget('Te')
    Flags_LP=cdf.varget('Flags_LP')
    Flags_Ne=cdf.varget('Flags_Ne')
    Flags_Te=cdf.varget('Flags_Te')
            
    cdf.close()
    
    height=Radius/1000.-MEAN_EARTH_RADIUS
    
    date=datetime(Timestamp[0][0],Timestamp[0][1],Timestamp[0][2])  
    refh=np.nanmean(height)
    A=Apex(date,refh)

    year=[]
    month=[]
    day=[]
    hour=[]
    minute=[]
    second=[]
    millisecond=[]
    hourUT=[]
    doy=[]
    mlat=[]
    mlon=[]
    MLT=[]
    for i in range(len(Timestamp)):
        year.append(Timestamp[i][0])
        month.append(Timestamp[i][1])
        day.append(Timestamp[i][2])
        hour.append(Timestamp[i][3])
        minute.append(Timestamp[i][4])
        second.append(Timestamp[i][5])
        millisecond.append(Timestamp[i][6])            
        doy.append(int(datetime(year[i],month[i],day[i]).strftime('%j')))
        hourUT.append((3600.*hour[i]+60*minute[i]+second[i]+0.001*millisecond[i])/3600.)
        
        #Magnetic coordinates        
        mag_lat,mag_lon=A.convert(Latitude[i],Longitude[i],'geo','qd',height[i])  
        mlt=A.mlon2mlt(mag_lon,datetime(year[i],month[i],day[i],hour[i],minute[i],second[i])) 
        mlat.append(mag_lat)
        mlon.append(mag_lon)
        MLT.append(mlt)

        
        
    year=np.array(year)
    month=np.array(month)
    day=np.array(day)
    hour=np.array(hour)
    minute=np.array(minute)
    second=np.array(second)
    millisecond=np.array(millisecond)     
    doy=np.array(doy)
    hourUT=np.array(hourUT)
    hourLT=(hourUT+Longitude/15.)%24      
    mlat=np.array(mlat)    
    mlon=np.array(mlon)    
    MLT=np.array(MLT)       
    
    
    for i in range(len(year)):
        f.write('%4i\t\t%2i\t\t%2i\t\t%2i\t\t%2i\t\t%2i\t\t%3i\t\t%3i\t\t%8.5f\t\t%8.5f\t\t%8.5f\t\t%10.5f\t\t%10.5f\t\t%10.2f\t\t%10.2f\t\t%10.5f\t\t%10.5f\t\t%8i\t\t%6i\t\t%3i\t\t%3i\t\t%3i\n' % \
                 (year[i],month[i],day[i],hour[i],minute[i],second[i],millisecond[i],doy[i],hourUT[i],hourLT[i],MLT[i],Latitude[i],Longitude[i],Radius[i],height[i],mlat[i],mlon[i],Ne[i],Te[i],Flags_LP[i],Flags_Ne[i],Flags_Te[i]))                
    
    f.close()

###############################################################################
#Working on TEC data
###############################################################################
    
filename_swarm=os.listdir(path_downloaded_data_TEC)
filename_swarm.sort()

PRN_index=np.linspace(1,32,32,dtype=int)

for file in filename_swarm:
    
    SAT=file[11:12]
    YEAR=file[19:23]
    MONTH=file[23:25]
    DOM=file[25:27]

    filename_swarm_out=[]    
    for prn in PRN_index:
        filename_swarm_out.append('Swarm_TEC_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DOM+'_PRN'+str(prn).zfill(2)+'_data.txt')
    
    f=[]
    for i in range(len(PRN_index)):  
        f.append(open(os.path.join(path_organized_data_TEC,filename_swarm_out[i]),'w'))
        f[i].write('year          month	   	day	  	hour	   	min	   	sec	   	msec	    	doy	      	  hourUT		  hourLT	   	   MLT		    	LatGeo_LEO		LonGeo_LEO		Radius [m]	       Height [km]	         LatMag_QD		 LonMag_QD	      Abs_STEC [TECU]       Abs_VTEC [TECU] 	     Rel_STEC [TECU] 	  Rel_STEC_RMS [TECU]   		Elev_Angle     	PRN         	   GPS_pos_X [m]                  GPS_pos_Y [m]                  GPS_pos_Z [m]                  LEO_pos_X [m]                  LEO_pos_Y [m]                  LEO_pos_Z [m]\n')
    
    cdf = cdflib.CDF(os.path.join(path_downloaded_data_TEC,file))
        
    #reading variables
    Timestamp=cdf.varget('Timestamp')   
    Timestamp=cdflib.cdfepoch.breakdown(Timestamp)
    Latitude=cdf.varget('Latitude')
    Longitude=cdf.varget('Longitude')
    Radius=cdf.varget('Radius')
    PRN=cdf.varget('PRN')
    Absolute_STEC=cdf.varget('Absolute_STEC')
    Absolute_VTEC=cdf.varget('Absolute_VTEC')
    Relative_STEC=cdf.varget('Relative_STEC')
    Relative_STEC_RMS=cdf.varget('Relative_STEC_RMS')
    Elevation_Angle=cdf.varget('Elevation_Angle')
    GPS_Position=cdf.varget('GPS_Position')
    LEO_Position=cdf.varget('LEO_Position')
    
        
    cdf.close()
    
    height=Radius/1000.-MEAN_EARTH_RADIUS
    
    date=datetime(Timestamp[0][0],Timestamp[0][1],Timestamp[0][2])  
    refh=np.nanmean(height)
    
    A=Apex(date,refh)

    GPS_Position_X=np.array(GPS_Position[:,0])
    GPS_Position_Y=np.array(GPS_Position[:,1])
    GPS_Position_Z=np.array(GPS_Position[:,2])
    LEO_Position_X=np.array(LEO_Position[:,0])
    LEO_Position_Y=np.array(LEO_Position[:,1])
    LEO_Position_Z=np.array(LEO_Position[:,2])    
    
        
    year=[]
    month=[]
    day=[]
    hour=[]
    minute=[]
    second=[]
    millisecond=[]
    hourUT=[]
    doy=[]
    mlat=[]
    mlon=[]
    MLT=[]
    for i in range(len(Timestamp)):
        year.append(Timestamp[i][0])
        month.append(Timestamp[i][1])
        day.append(Timestamp[i][2])
        hour.append(Timestamp[i][3])
        minute.append(Timestamp[i][4])
        second.append(Timestamp[i][5])
        millisecond.append(Timestamp[i][6])            
        doy.append(int(datetime(year[i],month[i],day[i]).strftime('%j')))
        hourUT.append((3600.*hour[i]+60*minute[i]+second[i]+0.001*millisecond[i])/3600.)
        
        #Magnetic coordinates        
        mag_lat,mag_lon=A.convert(Latitude[i],Longitude[i],'geo','qd',height[i])  
        mlt=A.mlon2mlt(mag_lon,datetime(year[i],month[i],day[i],hour[i],minute[i],second[i])) 
        mlat.append(mag_lat)
        mlon.append(mag_lon)
        MLT.append(mlt)
        

    year=np.array(year)
    month=np.array(month)
    day=np.array(day)
    hour=np.array(hour)
    minute=np.array(minute)
    second=np.array(second)
    millisecond=np.array(millisecond)     
    doy=np.array(doy)
    hourUT=np.array(hourUT)
    hourLT=(hourUT+Longitude/15.)%24      
    mlat=np.array(mlat)    
    mlon=np.array(mlon)    
    MLT=np.array(MLT)       
    
        
    for i in range(len(year)):
        f[(PRN[i]-1)].write('%4i\t\t%2i\t\t%2i\t\t%2i\t\t%2i\t\t%2i\t\t%3i\t\t%3i\t\t%8.5f\t\t%8.5f\t\t%8.5f\t\t%10.5f\t\t%10.5f\t\t%10.2f\t\t%10.2f\t\t%10.5f\t\t%10.5f\t\t%10.5f\t\t%10.5f\t\t%10.5f\t\t%10.5f\t\t\t%7.3f\t\t%2i\t\t%16.5f\t\t%16.5f\t\t%16.5f\t\t%16.5f\t\t%16.5f\t\t%16.5f\n' % \
         (year[i],month[i],day[i],hour[i],minute[i],second[i],millisecond[i],doy[i],hourUT[i],hourLT[i],MLT[i],\
          Latitude[i],Longitude[i],Radius[i],height[i],mlat[i],mlon[i],\
          Absolute_STEC[i],Absolute_VTEC[i],Relative_STEC[i],Relative_STEC_RMS[i],Elevation_Angle[i],PRN[i],GPS_Position_X[i],GPS_Position_Y[i],\
          GPS_Position_Z[i],LEO_Position_X[i],LEO_Position_Y[i],LEO_Position_Z[i]))
    
    for i in range(len(PRN_index)):
        f[i].close()