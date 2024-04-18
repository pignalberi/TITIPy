# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 15:22:41 2019

@author: alessio
"""

import os
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.basemap import Basemap
import datetime
from spacepy import coordinates
from spacepy.time import Ticktock
from apexpy import Apex
    
from Functions import *
from Terminal_interface import *

import warnings

warnings.filterwarnings("ignore") 

print('\n**********Swarm data analysis module**********')

plt.ioff()

MEAN_EARTH_RADIUS=6371.007 #mean Earth radius IUGG
SHELL_HEIGHT=400          #height of the thin shell (above Swarm satellite) in km

###############################################################################
#Defining paths
###############################################################################

main_folder=os.getcwd()

path_initial = os.path.join(main_folder,str(YEAR).zfill(4)+str(MONTH).zfill(2)+str(DOM).zfill(2)+str(SAT))
path_organized_data_LP=os.path.join(path_initial,'Organized_data','LP')
path_organized_data_TEC=os.path.join(path_initial,'Organized_data','TEC')

path_output=os.path.join(path_initial,'Output')

###############################################################################
#Working on Langmuir Probes data
###############################################################################

filename_swarm=os.listdir(path_organized_data_LP)
filename_swarm.sort()


#Parameters for indices calculation
FLAG=np.nan 
FREQUENCY=2                              #frequency of the input time series
WINDOW=int(WINDOW_WIDTH*FREQUENCY+1)     #number of points in the windows (it must be odd)
GAP=int(WINDOW/2)                        #number of missing data allowed in the window


for filename in filename_swarm:
    
    SAT=filename[9:10]
    YEAR=filename[11:15]
    MONTH=filename[16:18]
    DAY=filename[19:21]
    
    print('Working on Swarm '+SAT+' '+YEAR+'/'+MONTH+'/'+DAY+' LP data...')
    
    ###############################################################################
    #Reading variables in files
    ###############################################################################
    #Time
    year,month,day,hour,minute,second,msec,doy=\
    np.genfromtxt(os.path.join(path_organized_data_LP,filename),unpack=True,usecols=(0,1,2,3,4,5,6,7),dtype=int,skip_header=1)
    
    hourUT,hourLT,MLT=\
    np.genfromtxt(os.path.join(path_organized_data_LP,filename),unpack=True,usecols=(8,9,10),dtype=float,skip_header=1)
    
    #Position
    lat_geo,lon_geo,radius,height,lat_mag,lon_mag=\
    np.genfromtxt(os.path.join(path_organized_data_LP,filename),unpack=True,usecols=(11,12,13,14,15,16),dtype=float,skip_header=1)
    
    #Values
    Ne,Te,flag_LP,flag_Ne,flag_Te=\
    np.genfromtxt(os.path.join(path_organized_data_LP,filename),unpack=True,usecols=(17,18,19,20,21),dtype=int,skip_header=1)
    

    Ne=list(Ne)
    #Ne filter based on flags
    for i in range(len(Ne)):
        if(flag_LP[i]!=1 or flag_Ne[i]>29):
            Ne[i]=np.nan
    
    Te=list(Te)
    #Te filter based on flags
    for i in range(len(Te)):
        if(flag_LP[i]!=1 or (flag_Te[i]!=10 and flag_Te[i]!=20)):
            Te[i]=np.nan
    
   
    ###############################################################################
    #RODI, ROTEI, and other LP parameteres calculation
    ###############################################################################
    date=[]
    for i in range(len(Ne)):
        if(msec[i]>=0 and msec[i]<500):
            msec_2=0
        if(msec[i]>=500 and msec[i]<=999):
            msec_2=500  
        date.append(datetime(year[i],month[i],day[i],hour[i],minute[i],second[i],msec_2*1000))
    
    date_list=[]
    for i in range(len(date)):
        #time in seconds from the beginning of the considered month
        time=(date[i].day-1)*86400+date[i].hour*3600+date[i].minute*60+date[i].second+date[i].microsecond/1e6
        date_list.append(time)
    
    print('    Calculating RODI, ROTEI, and other LP parameters...')
    
    #Ne
    rod_parameters=rod_rodi_calculation(Ne,date_list,WINDOW,FREQUENCY,GAP,FLAG)
    gradient=gradient_calculation(Ne,date_list,radius,lat_geo,lon_geo,FREQUENCY,FLAG)
    
    #Te
    roTe_parameters=rod_rodi_calculation(Te,date_list,WINDOW,FREQUENCY,GAP,FLAG)
    gradient_Te=gradient_calculation(Te,date_list,radius,lat_geo,lon_geo,FREQUENCY,FLAG)

    
    print('    Saving output file...')
    
    filename_output='Swarm_LP_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_output.txt'
    f=open(os.path.join(path_output,'data','LP',filename_output),'w')
    f.write('year          month		day	    	hour		min	    	sec	    	msec		doy	         hourUT		  hourLT		     MLT		    LatGeo		    LonGeo		Radius [m]	    	Height [km]	      	  LatMag_QD	      	  LonMag_QD		Ne [cm^-3]	  	Flag_LP      Flag_Ne	   	    ROD [cm^-3/s]       	   RODI [cm^-3/s]    		mean_ROD [cm^-3/s]		Ne_grad [cm^-3/km]		   mean_Ne [cm^-3]			  Te [K]	      Flag_Te	  	 	ROTE [K/s]    	   		ROTEI [K/s]   	 		mean_ROTE [K/s]		Te_grad [K/km]			mean_Te [K]\n')
    for i in range(len(Ne)):
        f.write('%4i\t\t%2i\t\t%2i\t\t%2i\t\t%2i\t\t%2i\t\t%3i\t\t%3i\t\t%8.5f\t\t%8.5f\t\t%8.5f\t\t%10.5f\t\t%10.5f\t\t%10.2f\t\t%10.2f\t\t%10.5f\t\t%10.5f\t\t%9.0f\t\t%3i\t\t%3i\t\t\t%9.1f\t\t\t%9.1f\t\t\t%9.1f\t\t\t%9.1f\t\t\t%9.1f\t\t\t%8.0f\t\t%3i\t\t\t%9.1f\t\t\t%9.1f\t\t\t%9.1f\t\t\t%9.1f\t\t\t%9.1f\n' % \
                (year[i],month[i],day[i],hour[i],minute[i],second[i],msec[i],doy[i],hourUT[i],hourLT[i],MLT[i],lat_geo[i],lon_geo[i],radius[i],height[i],lat_mag[i],lon_mag[i],Ne[i],flag_LP[i],flag_Ne[i],rod_parameters[1][i],rod_parameters[2][i],rod_parameters[3][i],gradient[1][i],rod_parameters[5][i],Te[i],flag_Te[i],roTe_parameters[1][i],roTe_parameters[2][i],roTe_parameters[3][i],gradient_Te[1][i],roTe_parameters[5][i]))
    
    f.close()


    
    if(FIGURE):
        
        ###########################################################################
        #Binning in Magnetic coordinates
        ###########################################################################
        print('    Binning data in magnetic coordinates...')
        
        #North Pole
        lat_QD_NP=np.ma.masked_less(lat_mag,40)
        MLT_map_NP=np.ma.masked_array(MLT,lat_QD_NP.mask)#/24.*(2*np.pi)
        lat_QD_NP=-1.*(lat_QD_NP-90.)
        
        for i in range(len(lat_QD_NP)):
            if(lat_QD_NP[i]>=50):
                lat_QD_NP[i]=np.nan
                MLT_map_NP[i]=np.nan
               
        for i in range(len(MLT_map_NP)):
            if(MLT_map_NP[i]>=24):
                MLT_map_NP[i]=np.nan
                lat_QD_NP[i]=np.nan
            
        binned_rodi_QD_MLT_NP=binning_2D(MLT_map_NP,lat_QD_NP,rod_parameters[2],0,24,0.25,0,50,2)  
        binned_roTei_QD_MLT_NP=binning_2D(MLT_map_NP,lat_QD_NP,roTe_parameters[2],0,24,0.25,0,50,2)  
        binned_Ne_QD_MLT_NP=binning_2D(MLT_map_NP,lat_QD_NP,Ne,0,24,0.25,0,50,2)  
        binned_Te_QD_MLT_NP=binning_2D(MLT_map_NP,lat_QD_NP,Te,0,24,0.25,0,50,2)  
    
        
        #South Pole
        lat_QD_SP=np.ma.masked_greater(lat_mag,-40)
        MLT_map_SP=np.ma.masked_array(MLT,lat_QD_SP.mask)#/24.*(2*np.pi)
        lat_QD_SP=lat_QD_SP+90.
        
        for i in range(len(lat_QD_SP)):
            if(lat_QD_SP[i]>=50):
                lat_QD_SP[i]=np.nan
                MLT_map_SP[i]=np.nan
        
        for i in range(len(MLT_map_SP)):
            if(MLT_map_SP[i]>=24):
                MLT_map_SP[i]=np.nan
                lat_QD_SP[i]=np.nan
                
        binned_rodi_QD_MLT_SP=binning_2D(MLT_map_SP,lat_QD_SP,rod_parameters[2],0,24,0.25,0,50,2)  
        binned_roTei_QD_MLT_SP=binning_2D(MLT_map_SP,lat_QD_SP,roTe_parameters[2],0,24,0.25,0,50,2)  
        binned_Ne_QD_MLT_SP=binning_2D(MLT_map_SP,lat_QD_SP,Ne,0,24,0.25,0,50,2)  
        binned_Te_QD_MLT_SP=binning_2D(MLT_map_SP,lat_QD_SP,Te,0,24,0.25,0,50,2)  
        
    
        binned_rodi_QD_MLT=binning_2D(MLT,lat_mag,rod_parameters[2],0,24,0.25,-90,90,2) 
        binned_roTei_QD_MLT=binning_2D(MLT,lat_mag,roTe_parameters[2],0,24,0.25,-90,90,2)  
        binned_Ne_QD_MLT=binning_2D(MLT,lat_mag,Ne,0,24,0.25,-90,90,2) 
        binned_Te_QD_MLT=binning_2D(MLT,lat_mag,Te,0,24,0.25,-90,90,2) 
    
    
    
        ###############################################################################
        #RODI figures
        ###############################################################################
        logo = plt.imread(os.path.join(main_folder,'logo.png'))
        
        print('    Making output figures...')
        
        #Map scatter lat-lon Global
        fig=plt.figure(figsize=(18,10))
        
        m=Basemap(projection='robin',lon_0=0,llcrnrlat=-90,urcrnrlat=90,\
         llcrnrlon=-180,urcrnrlon=180,resolution='l')     
        
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,30),labels=[1,1,0,0],fontsize=15)
        m.drawmeridians(np.arange(-180,180,60),labels=[0,0,1,1],fontsize=15)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)
               
        x, y = m(lon_geo,lat_geo)
        
        fig=m.scatter(x,y,s=50,c=np.log10(rod_parameters[2]),edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=3,vmax=5)
            
        cbar=m.colorbar(fig,location='right',pad="5%",size="2.3%")#,format='%3.1e')
        cbar.set_label('$log_{10}(RODI)\ (log_{10}(cm^{-3}/s))$', fontsize=15)
        cbar.ax.tick_params(labelsize=15)
    
        plt.figimage(logo,0,0,zorder=-1,alpha=1)
    
        plt.title('Scatter RODI, Global, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','RODI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_Global_scatter.png'), format='png',dpi=150)
        
        plt.clf()
        plt.close()
        
        #Map scatter lat-lon North Pole
        fig=plt.figure(figsize=(12,10))
        
        m = Basemap(lon_0=0, boundinglat=40, projection='npstere',round=True,resolution='l')
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,10),labels=[1,1,1,1],fontsize=15)
        m.drawmeridians(np.arange(0,360,30),labels=[1,1,1,1],fontsize=15)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)
        
        x, y = m(lon_geo,lat_geo)
    
        fig=m.scatter(x,y,s=50,c=np.log10(rod_parameters[2]),edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=3,vmax=5)
    
        cbar=m.colorbar(fig,location='right',pad="10%",size="5%")#,format='%3.1e')
        cbar.set_label('$log_{10}(RODI)\ (log_{10}(cm^{-3}/s))$', fontsize=15)
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Scatter RODI, North Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','RODI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_North_Pole_scatter.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        #Map scatter lat-lon South Pole
        fig=plt.figure(figsize=(12,10))

        m = Basemap(lon_0=0, boundinglat=-40, projection='spstere',round=True,resolution='l')
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,10),labels=[1,1,1,1],fontsize=15)
        m.drawmeridians(np.arange(0,360,30),labels=[1,1,1,1],fontsize=15)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)        
               
        x, y = m(lon_geo,lat_geo)
    
        fig=m.scatter(x,y,s=50,c=np.log10(rod_parameters[2]),edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=3,vmax=5)
        
        cbar=m.colorbar(fig,location='right',pad="10%",size="5%")#,format='%3.1e')
        cbar.set_label('$log_{10}(RODI)\ (log_{10}(cm^{-3}/s)])$', fontsize=15)       
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Scatter RODI, South Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','RODI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_South_Pole_scatter.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        ###########################################################################
        #QD-MLT polar plots
        ###########################################################################
        
        #North Pole
        fig=plt.figure(figsize=(12,12))
        
        rad=np.linspace(0,50,26, endpoint='True')
        azm=np.linspace(0,24,97, endpoint='True')/24.*(2*np.pi)
        r, th = np.meshgrid(rad, azm)
        z=np.log10(binned_rodi_QD_MLT_NP[0].T)
    
        ax=plt.subplot(projection="polar")
        plt.pcolormesh(th-np.pi/2., r, z,cmap=plt.cm.jet,vmin=3,vmax=5)
        plt.plot(azm, r, color='k', ls='none') 
        plt.grid()
        plt.xticks(np.linspace(0,2*np.pi,24,endpoint=False),(np.concatenate((np.linspace(6,23,18,dtype=int),np.linspace(0,5,6,dtype=int)))),fontsize=15)
    
        ax.set_rorigin(0)
        ax.set_theta_zero_location('E', offset=0)
        plt.yticks(np.linspace(0,40,5),['90','80','70','60','50'],fontsize=15,alpha=1,color='black')
    
        plt.xlabel('MLT-Magnetic Latitude ', fontsize=20)
        cbar=plt.colorbar(orientation='horizontal',fraction=0.05)
        cbar.set_label('$log_{10}$(RODI) ($log_{10}$(el/cm$^{3}$/s))', fontsize=15)        
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Magnetic polar plot RODI, North Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','RODI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_North_Pole_magnetic_polar_plot.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        
        #South Pole
        fig=plt.figure(figsize=(12,12))
        
        rad=np.linspace(0,50,26, endpoint='True')
        azm=np.linspace(0,24,97, endpoint='True')/24.*(2*np.pi)
        r, th = np.meshgrid(rad, azm)
        z=np.log10(binned_rodi_QD_MLT_SP[0].T)
        
        ax=plt.subplot(projection="polar")
        plt.pcolormesh(th-np.pi/2., r, z,cmap=plt.cm.jet,vmin=3,vmax=5)
        plt.plot(azm, r, color='k', ls='none') 
        plt.grid()
        plt.xticks(np.linspace(0,2*np.pi,24,endpoint=False),(np.concatenate((np.linspace(6,23,18,dtype=int),np.linspace(0,5,6,dtype=int)))),fontsize=15)
        
        ax.set_rorigin(0)
        ax.set_theta_zero_location('E', offset=0)
        plt.yticks(np.linspace(0,40,5),['-90','-80','-70','-60','-50'],fontsize=15,alpha=1,color='black')
        
        plt.xlabel('MLT-Magnetic Latitude ', fontsize=20)
        cbar=plt.colorbar(orientation='horizontal',fraction=0.05)
        cbar.set_label('$log_{10}$(RODI) ($log_{10}$(el/cm$^{3}$/s))', fontsize=15)
        
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Magnetic polar plot RODI, South Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','RODI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_South_Pole_magnetic_polar_plot.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        
        #Global        
        fig=plt.figure(figsize=(18,10))
        
        y=np.linspace(-90,90,91, endpoint='True')
        x=np.linspace(0,24,97,endpoint='True')
        X,Y = np.meshgrid(x,y)
        z=np.log10(binned_rodi_QD_MLT[0])
    
        plt.pcolormesh(X,Y,z,cmap=plt.cm.jet,vmin=3,vmax=5)
        plt.grid()
        plt.xticks(np.linspace(0,24,25,endpoint=True),fontsize=15)
        plt.yticks(np.linspace(-90,90,7),fontsize=15)
        plt.xlabel('MLT', fontsize=20)
        plt.ylabel('Magnetic Latitude', fontsize=20)
        cbar=plt.colorbar(orientation='vertical',fraction=0.05)
        cbar.set_label('$log_{10}$(RODI) ($log_{10}$(el/cm$^{3}$/s))', fontsize=20)
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Magnetic coordinates plot RODI, Global, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','RODI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_Global_magnetic_coordinates.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
    
    
        ###############################################################################
        #ROTEI figures
        ###############################################################################    
        #Map scatter lat-lon Global
        fig=plt.figure(figsize=(18,10))
        
        m=Basemap(projection='robin',lon_0=0,llcrnrlat=-90,urcrnrlat=90,\
         llcrnrlon=-180,urcrnrlon=180,resolution='l')     
        
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,30),labels=[1,1,0,0],fontsize=15)
        m.drawmeridians(np.arange(-180,180,60),labels=[0,0,1,1],fontsize=15)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)
               
        x, y = m(lon_geo,lat_geo)
        
        fig=m.scatter(x,y,s=50,c=np.log10(roTe_parameters[2]),edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=0,vmax=3.5)
                
        cbar=m.colorbar(fig,location='right',pad="5%",size="2.3%")#,format='%3.1e')
        cbar.set_label('$log_{10}(ROTEI)\ (log_{10}(K/s))$', fontsize=15)
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Scatter ROTEI, Global, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','ROTEI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_Global_scatter.png'), format='png',dpi=150)
        
        plt.clf()
        plt.close()
        
        #Map scatter lat-lon North Pole
        fig=plt.figure(figsize=(12,10))
        
        m = Basemap(lon_0=0, boundinglat=40, projection='npstere',round=True,resolution='l')
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,10),labels=[1,1,1,1],fontsize=15)
        m.drawmeridians(np.arange(0,360,30),labels=[1,1,1,1],fontsize=15)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)
               
        x, y = m(lon_geo,lat_geo)
    
        fig=m.scatter(x,y,s=50,c=np.log10(roTe_parameters[2]),edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=0,vmax=3.5)
        
        cbar=m.colorbar(fig,location='right',pad="10%",size="5%")#,format='%3.1e')
        cbar.set_label('$log_{10}(ROTEI)\ (log_{10}(K/s))$', fontsize=15)
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Scatter ROTEI, North Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','ROTEI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_North_Pole_scatter.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        #Map scatter lat-lon South Pole
        fig=plt.figure(figsize=(12,10))
        
        m = Basemap(lon_0=0, boundinglat=-40, projection='spstere',round=True,resolution='l')
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,10),labels=[1,1,1,1],fontsize=15)
        m.drawmeridians(np.arange(0,360,30),labels=[1,1,1,1],fontsize=15)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)
               
        x, y = m(lon_geo,lat_geo)
    
        fig=m.scatter(x,y,s=50,c=np.log10(roTe_parameters[2]),edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=0,vmax=3.5)
        
        cbar=m.colorbar(fig,location='right',pad="10%",size="5%")#,format='%3.1e')
        cbar.set_label('$log_{10}(ROTEI)\ (log_{10}(K/s)])$', fontsize=15)        
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Scatter ROTEI, South Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','ROTEI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_South_Pole_scatter.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        ###########################################################################
        #QD-MLT polar maps
        ###########################################################################
        
        #North Pole
        fig=plt.figure(figsize=(12,12))
        
        rad=np.linspace(0,50,26, endpoint='True')
        azm=np.linspace(0,24,97, endpoint='True')/24.*(2*np.pi)
        r, th = np.meshgrid(rad, azm)
        z=np.log10(binned_roTei_QD_MLT_NP[0].T)
    
        ax=plt.subplot(projection="polar")
        plt.pcolormesh(th-np.pi/2., r, z,cmap=plt.cm.jet,vmin=0,vmax=3.5)
        plt.plot(azm, r, color='k', ls='none') 
        plt.grid()
        plt.xticks(np.linspace(0,2*np.pi,24,endpoint=False),(np.concatenate((np.linspace(6,23,18,dtype=int),np.linspace(0,5,6,dtype=int)))),fontsize=15)
    
        ax.set_rorigin(0)
        ax.set_theta_zero_location('E', offset=0)
        plt.yticks(np.linspace(0,40,5),['90','80','70','60','50'],fontsize=15,alpha=1,color='black')
    
        plt.xlabel('MLT-Magnetic Latitude ', fontsize=20)
        cbar=plt.colorbar(orientation='horizontal',fraction=0.05)
        cbar.set_label('$log_{10}$(ROTEI) ($log_{10}$(K/s))', fontsize=15)        
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Magnetic polar plot ROTEI, North Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','ROTEI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_North_Pole_magnetic_polar_plot.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        
        #South Pole
        fig=plt.figure(figsize=(12,12))
        
        rad=np.linspace(0,50,26, endpoint='True')
        azm=np.linspace(0,24,97, endpoint='True')/24.*(2*np.pi)
        r, th = np.meshgrid(rad, azm)
        z=np.log10(binned_roTei_QD_MLT_SP[0].T)
        
        ax=plt.subplot(projection="polar")
        plt.pcolormesh(th-np.pi/2., r, z,cmap=plt.cm.jet,vmin=0,vmax=3.5)
        plt.plot(azm, r, color='k', ls='none') 
        plt.grid()
        plt.xticks(np.linspace(0,2*np.pi,24,endpoint=False),(np.concatenate((np.linspace(6,23,18,dtype=int),np.linspace(0,5,6,dtype=int)))),fontsize=15)
        
        ax.set_rorigin(0)
        ax.set_theta_zero_location('E', offset=0)
        plt.yticks(np.linspace(0,40,5),['-90','-80','-70','-60','-50'],fontsize=15,alpha=1,color='black')
        
        plt.xlabel('MLT-Magnetic Latitude ', fontsize=20)
        cbar=plt.colorbar(orientation='horizontal',fraction=0.05)
        cbar.set_label('$log_{10}$(ROTEI) ($log_{10}$(K/s))', fontsize=15)        
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Magnetic polar plot ROTEI, South Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','ROTEI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_South_Pole_magnetic_polar_plot.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        
        #Global        
        fig=plt.figure(figsize=(18,10))
        
        y=np.linspace(-90,90,91, endpoint='True')
        x=np.linspace(0,24,97,endpoint='True')
        X,Y = np.meshgrid(x,y)
        z=np.log10(binned_roTei_QD_MLT[0])
    
        plt.pcolormesh(X,Y,z,cmap=plt.cm.jet,vmin=0,vmax=3.5)
        plt.grid()
        plt.xticks(np.linspace(0,24,25,endpoint=True),fontsize=15)
        plt.yticks(np.linspace(-90,90,7),fontsize=15)
        plt.xlabel('MLT', fontsize=20)
        plt.ylabel('Magnetic Latitude', fontsize=20)
        
        cbar=plt.colorbar(orientation='vertical',fraction=0.05)
        cbar.set_label('$log_{10}$(ROTEI) ($log_{10}$(K/s))', fontsize=20)
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)   
        
        plt.title('Magnetic coordinates plot ROTEI, Global, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','ROTEI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_Global_magnetic_coordinates.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
    
        
        ###############################################################################
        #Ne figures
        ###############################################################################    
        #Map scatter lat-lon Global
        fig=plt.figure(figsize=(18,10))
        
        m=Basemap(projection='robin',lon_0=0,llcrnrlat=-90,urcrnrlat=90,\
         llcrnrlon=-180,urcrnrlon=180,resolution='l')     
        
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,30),labels=[1,1,0,0],fontsize=15)
        m.drawmeridians(np.arange(-180,180,60),labels=[0,0,1,1],fontsize=15)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)
               
        x, y = m(lon_geo,lat_geo)
        
        fig=m.scatter(x,y,s=50,c=np.log10(Ne),edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=3,vmax=6)
            
        cbar=m.colorbar(fig,location='right',pad="10%",size="2.3%")#,format='%3.1e')
        cbar.set_label('$log_{10}$(Ne) ($log_{10}$(el/cm$^{3}$))', fontsize=15)
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Scatter Ne, Global, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','Ne_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_Global_scatter.png'), format='png',dpi=150)
        
        plt.clf()
        plt.close()
        
        #Map scatter lat-lon North Pole
        fig=plt.figure(figsize=(12,10))
        
        m = Basemap(lon_0=0, boundinglat=40, projection='npstere',round=True,resolution='l')
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,10),labels=[1,1,1,1],fontsize=15)
        m.drawmeridians(np.arange(0,360,30),labels=[1,1,1,1],fontsize=15)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)
               
        x, y = m(lon_geo,lat_geo)
    
        fig=m.scatter(x,y,s=50,c=np.log10(Ne),edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=3,vmax=6)
    
        cbar=m.colorbar(fig,location='right',pad="10%",size="5%")#,format='%3.1e')
        cbar.set_label('$log_{10}$(Ne) ($log_{10}$(el/cm$^{3}$))', fontsize=15)
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Scatter Ne, North Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','Ne_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_North_Pole_scatter.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        #Map scatter lat-lon South Pole
        fig=plt.figure(figsize=(12,10))
        
        m = Basemap(lon_0=0, boundinglat=-40, projection='spstere',round=True,resolution='l')
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,10),labels=[1,1,1,1],fontsize=15)
        m.drawmeridians(np.arange(0,360,30),labels=[1,1,1,1],fontsize=15)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)
               
        x, y = m(lon_geo,lat_geo)
    
        fig=m.scatter(x,y,s=50,c=np.log10(Ne),edgecolors='None', linewidths=None,cmap=plt.cm.jet)#,vmin=3,vmax=5)
        
        cbar=m.colorbar(fig,location='right',pad="10%",size="5%")#,format='%3.1e')
        cbar.set_label('$log_{10}$(Ne) ($log_{10}$(el/cm$^{3}$))', fontsize=15)        
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Scatter Ne, South Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','Ne_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_South_Pole_scatter.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        ###########################################################################
        #QD-MLT polar plots
        ###########################################################################
        
        #North Pole
        fig=plt.figure(figsize=(12,12))
        
        rad=np.linspace(0,50,26, endpoint='True')
        azm=np.linspace(0,24,97, endpoint='True')/24.*(2*np.pi)
        r, th = np.meshgrid(rad, azm)
        z=np.log10(binned_Ne_QD_MLT_NP[0].T)
    
        ax=plt.subplot(projection="polar")
        plt.pcolormesh(th-np.pi/2., r, z,cmap=plt.cm.jet,vmin=3,vmax=6)
        plt.plot(azm, r, color='k', ls='none') 
        plt.grid()
        plt.xticks(np.linspace(0,2*np.pi,24,endpoint=False),(np.concatenate((np.linspace(6,23,18,dtype=int),np.linspace(0,5,6,dtype=int)))),fontsize=15)
    
        ax.set_rorigin(0)
        ax.set_theta_zero_location('E', offset=0)
        plt.yticks(np.linspace(0,40,5),['90','80','70','60','50'],fontsize=15,alpha=1,color='black')
    
        plt.xlabel('MLT-Magnetic Latitude ', fontsize=20)
        cbar=plt.colorbar(orientation='horizontal',fraction=0.05)
        cbar.set_label('$log_{10}$(Ne) ($log_{10}$(el/cm$^{3}$))', fontsize=15)        
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Magnetic polar plot Ne, North Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','Ne_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_North_Pole_magnetic_polar_plot.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        
        #South Pole
        fig=plt.figure(figsize=(12,12))
        
        rad=np.linspace(0,50,26, endpoint='True')
        azm=np.linspace(0,24,97, endpoint='True')/24.*(2*np.pi)
        r, th = np.meshgrid(rad, azm)
        z=np.log10(binned_Ne_QD_MLT_SP[0].T)
        
        ax=plt.subplot(projection="polar")
        plt.pcolormesh(th-np.pi/2., r, z,cmap=plt.cm.jet,vmin=3,vmax=6)
        plt.plot(azm, r, color='k', ls='none') 
        plt.grid()
        plt.xticks(np.linspace(0,2*np.pi,24,endpoint=False),(np.concatenate((np.linspace(6,23,18,dtype=int),np.linspace(0,5,6,dtype=int)))),fontsize=15)
        
        ax.set_rorigin(0)
        ax.set_theta_zero_location('E', offset=0)
        plt.yticks(np.linspace(0,40,5),['-90','-80','-70','-60','-50'],fontsize=15,alpha=1,color='black')
        
        plt.xlabel('MLT-Magnetic Latitude ', fontsize=20)
        cbar=plt.colorbar(orientation='horizontal',fraction=0.05)
        cbar.set_label('$log_{10}$(Ne) ($log_{10}$(el/cm$^{3}$))', fontsize=15)        
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Magnetic polar plot Ne, South Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','Ne_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_South_Pole_magnetic_polar_plot.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        
        #Global        
        fig=plt.figure(figsize=(18,10))
        
        y=np.linspace(-90,90,91, endpoint='True')
        x=np.linspace(0,24,97,endpoint='True')
        X,Y = np.meshgrid(x,y)
        z=np.log10(binned_Ne_QD_MLT[0])
    
        plt.pcolormesh(X,Y,z,cmap=plt.cm.jet,vmin=3,vmax=6)
        plt.grid()
        plt.xticks(np.linspace(0,24,25,endpoint=True),fontsize=15)
        plt.yticks(np.linspace(-90,90,7),fontsize=15)
        plt.xlabel('MLT', fontsize=20)
        plt.ylabel('Magnetic Latitude', fontsize=20)
        
        cbar=plt.colorbar(orientation='vertical',fraction=0.05)
        cbar.set_label('$log_{10}$(Ne) ($log_{10}$(el/cm$^{3}$))', fontsize=20)
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Magnetic coordinates plot Ne, Global, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','Ne_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_Global_magnetic_coordinates.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
    
    
        ###############################################################################
        #Te figures
        ###############################################################################    
        #Map scatter lat-lon Global
        fig=plt.figure(figsize=(18,10))
        
        m=Basemap(projection='robin',lon_0=0,llcrnrlat=-90,urcrnrlat=90,\
         llcrnrlon=-180,urcrnrlon=180,resolution='l')     
        
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,30),labels=[1,1,0,0],fontsize=15)
        m.drawmeridians(np.arange(-180,180,60),labels=[0,0,1,1],fontsize=15)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)
               
        x, y = m(lon_geo,lat_geo)
        
        fig=m.scatter(x,y,s=50,c=Te,edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=0,vmax=6000)
        
        cbar=m.colorbar(fig,location='right',pad="10%",size="2.3%")#,format='%3.1e')
        cbar.set_label('Te (K)', fontsize=15)
        cbar.ax.tick_params(labelsize=15)
 
        plt.figimage(logo,0,0,zorder=-1,alpha=1)
            
        plt.title('Scatter Te, Global, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','Te_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_Global_scatter.png'), format='png',dpi=150)
        
        plt.clf()
        plt.close()
        
        #Map scatter lat-lon North Pole
        fig=plt.figure(figsize=(12,10))
        
        m = Basemap(lon_0=0, boundinglat=40, projection='npstere',round=True,resolution='l')
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,10),labels=[1,1,1,1],fontsize=15)
        m.drawmeridians(np.arange(0,360,30),labels=[1,1,1,1],fontsize=15)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)
               
        x, y = m(lon_geo,lat_geo)
    
        fig=m.scatter(x,y,s=50,c=Te,edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=0,vmax=6000)
        
        cbar=m.colorbar(fig,location='right',pad="10%",size="5%")#,format='%3.1e')
        cbar.set_label('Te (K)', fontsize=15)
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Scatter Te, North Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','Te_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_North_Pole_scatter.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        #Map scatter lat-lon South Pole
        fig=plt.figure(figsize=(12,10))
        
        m = Basemap(lon_0=0, boundinglat=-40, projection='spstere',round=True,resolution='l')
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,10),labels=[1,1,1,1],fontsize=15)
        m.drawmeridians(np.arange(0,360,30),labels=[1,1,1,1],fontsize=15)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)
               
        x, y = m(lon_geo,lat_geo)
    
        fig=m.scatter(x,y,s=50,c=Te,edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=0,vmax=6000)
        
        cbar=m.colorbar(fig,location='right',pad="10%",size="5%")#,format='%3.1e')
        cbar.set_label('Te (K)', fontsize=15)        
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Scatter Te, South Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','Te_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_South_Pole_scatter.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        ###########################################################################
        #QD-MLT polar plots
        ###########################################################################
        
        #North Pole
        fig=plt.figure(figsize=(12,12))
        
        rad=np.linspace(0,50,26, endpoint='True')
        azm=np.linspace(0,24,97, endpoint='True')/24.*(2*np.pi)
        r, th = np.meshgrid(rad, azm)
        z=binned_Te_QD_MLT_NP[0].T
    
        ax=plt.subplot(projection="polar")
        plt.pcolormesh(th-np.pi/2., r, z,cmap=plt.cm.jet,vmin=0,vmax=6000)
        plt.plot(azm, r, color='k', ls='none') 
        plt.grid()
        plt.xticks(np.linspace(0,2*np.pi,24,endpoint=False),(np.concatenate((np.linspace(6,23,18,dtype=int),np.linspace(0,5,6,dtype=int)))),fontsize=15)
    
        ax.set_rorigin(0)
        ax.set_theta_zero_location('E', offset=0)
        plt.yticks(np.linspace(0,40,5),['90','80','70','60','50'],fontsize=15,alpha=1,color='black')
    
        plt.xlabel('MLT-Magnetic Latitude ', fontsize=20)
        cbar=plt.colorbar(orientation='horizontal',fraction=0.05)
        cbar.set_label('Te (K)', fontsize=15)        
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Magnetic polar plot Te, North Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','Te_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_North_Pole_magnetic_polar_plot.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        
        #South Pole
        fig=plt.figure(figsize=(12,12))
        
        rad=np.linspace(0,50,26, endpoint='True')
        azm=np.linspace(0,24,97, endpoint='True')/24.*(2*np.pi)
        r, th = np.meshgrid(rad, azm)
        z=binned_Te_QD_MLT_SP[0].T
        
        ax=plt.subplot(projection="polar")
        plt.pcolormesh(th-np.pi/2., r, z,cmap=plt.cm.jet,vmin=0,vmax=6000)
        plt.plot(azm, r, color='k', ls='none') 
        plt.grid()
        plt.xticks(np.linspace(0,2*np.pi,24,endpoint=False),(np.concatenate((np.linspace(6,23,18,dtype=int),np.linspace(0,5,6,dtype=int)))),fontsize=15)
        
        ax.set_rorigin(0)
        ax.set_theta_zero_location('E', offset=0)
        plt.yticks(np.linspace(0,40,5),['-90','-80','-70','-60','-50'],fontsize=15,alpha=1,color='black')
        
        plt.xlabel('MLT-Magnetic Latitude ', fontsize=20)
        cbar=plt.colorbar(orientation='horizontal',fraction=0.05)
        cbar.set_label('Te (K)', fontsize=15)       
        cbar.ax.tick_params(labelsize=15)
 
        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Magnetic polar plot Te, South Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','Te_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_South_Pole_magnetic_polar_plot.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        
        #Global        
        fig=plt.figure(figsize=(18,10))
        
        y=np.linspace(-90,90,91, endpoint='True')
        x=np.linspace(0,24,97,endpoint='True')
        X,Y = np.meshgrid(x,y)
        z=binned_Te_QD_MLT[0]
    
        plt.pcolormesh(X,Y,z,cmap=plt.cm.jet,vmin=0,vmax=6000)
        plt.grid()
        plt.xticks(np.linspace(0,24,25,endpoint=True),fontsize=15)
        plt.yticks(np.linspace(-90,90,7),fontsize=15)
        plt.xlabel('MLT', fontsize=20)
        plt.ylabel('Magnetic Latitude', fontsize=20)
        
        cbar=plt.colorbar(orientation='vertical',fraction=0.05)
        cbar.set_label('Te (K)', fontsize=15)
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Magnetic coordinates plot Te, Global, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','LP','Te_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_Global_magnetic_coordinates.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()    
        
        

###############################################################################
#Working on TEC data
###############################################################################

PRN_index=np.linspace(1,32,32,dtype=int)

filename_swarm=os.listdir(path_organized_data_TEC)
filename_swarm.sort()


#Parameters for indices calculation
FLAG=np.nan 
if(FREQUENCY_CHANGE_TEC):
    FREQUENCY=0.1               #frequency of the input time series
    WINDOW=int((WINDOW_WIDTH*10)*FREQUENCY+1)     #number of points in the windows (it must be odd)

else:
    FREQUENCY=1                 #frequency of the input time series
    WINDOW=int(WINDOW_WIDTH*FREQUENCY+1)     #number of points in the windows (it must be odd)
GAP=int(WINDOW/2)             #number of missing data allowed in the window


lat_mag_complete=[]
MLT_complete=[]
TEC_for_ROTI_complete=[]
ROTI_complete=[]
for filename in filename_swarm:
    
    SAT=filename[10:11]
    YEAR=filename[12:16]
    MONTH=filename[17:19]
    DAY=filename[20:22]
    PRN=filename[26:28]
    
    print('Working on Swarm '+SAT+' '+YEAR+'/'+MONTH+'/'+DAY+' PRN='+PRN+' TEC data...')

    if(os.stat(os.path.join(path_organized_data_TEC,filename)).st_size<1000):
        print('    GPS satellite outside FOV')
        continue
        
    ###############################################################################
    #Reading variables in files
    ###############################################################################
    #Time
    year,month,day,hour,minute,second,msec,doy=\
    np.genfromtxt(os.path.join(path_organized_data_TEC,filename),unpack=True,usecols=(0,1,2,3,4,5,6,7),dtype=int,skip_header=1)
    
    hourUT,hourLT,MLT=\
    np.genfromtxt(os.path.join(path_organized_data_TEC,filename),unpack=True,usecols=(8,9,10),dtype=float,skip_header=1)
    
    #Position
    lat_geo,lon_geo,radius,height,lat_mag,lon_mag=\
    np.genfromtxt(os.path.join(path_organized_data_TEC,filename),unpack=True,usecols=(11,12,13,14,15,16),dtype=float,skip_header=1)
    
    #Values
    abs_sTEC,abs_vTEC,rel_sTEC,rel_sTEC_RMS,elev_angle=\
    np.genfromtxt(os.path.join(path_organized_data_TEC,filename),unpack=True,usecols=(17,18,19,20,21),dtype=float,skip_header=1)
    
    #PRN
    prn=np.genfromtxt(os.path.join(path_organized_data_TEC,filename),unpack=True,usecols=(22),dtype=int,skip_header=1)    

    #GPS and LEO position
    GPS_Position_X,GPS_Position_Y,GPS_Position_Z,LEO_Position_X,LEO_Position_Y,LEO_Position_Z=\
    np.genfromtxt(os.path.join(path_organized_data_TEC,filename),unpack=True,usecols=(23,24,25,26,27,28),dtype=float,skip_header=1)
    
    
    #Control on prn
    for i in range(len(prn)):
        if(prn[i]!=int(PRN)):
           raise Exception('PRN not coincidents')
    
    #lat_mag_complete.append(lat_mag)
    #MLT_complete.append(MLT)
    
    ###############################################################################
    #ROTI and other TEC paramaters calculation
    ###############################################################################
    date=[]
    date_str=[]
    for i in range(len(prn)):
        date.append(datetime(year[i],month[i],day[i],hour[i],minute[i],second[i]))
        date_str.append(str(datetime(year[i],month[i],day[i],hour[i],minute[i],second[i])))
    
    date_list=[]
    for i in range(len(date)):
        #time in seconds from the beginning of the considered month
        time=(date[i].day-1)*86400+date[i].hour*3600+date[i].minute*60+date[i].second
        date_list.append(time)
    
    print('    Calculating ROTI...')
    
    if(TEC_INPUT==1):
        TEC_for_ROTI=abs_sTEC
        label_TEC='Absolute sTEC (TECU)'
    if(TEC_INPUT==2):
        TEC_for_ROTI=abs_vTEC
        label_TEC='Absolute vTEC (TECU)'
    if(TEC_INPUT==3):
        TEC_for_ROTI=rel_sTEC     
        label_TEC='Relative sTEC (TECU)'
        
    TEC_for_ROTI_complete.append(TEC_for_ROTI)
    
    
    #cartesian geodetic WGS84 coordinates
    GPS_position_wgs84 = coordinates.Coords(np.matrix((GPS_Position_X,GPS_Position_Y,GPS_Position_Z)).T, 'GEO', 'car')
    GPS_position_wgs84.ticks = Ticktock(date_str, 'ISO')
    LEO_position_wgs84 = coordinates.Coords(np.matrix((LEO_Position_X,LEO_Position_Y,LEO_Position_Z)).T, 'GEO', 'car') 
    LEO_position_wgs84.ticks = Ticktock(date_str, 'ISO')
    
    #converting in geographic spherical coordinates
    GPS_position_geo = GPS_position_wgs84.convert('GEO', 'sph')
    LEO_position_geo = LEO_position_wgs84.convert('GEO', 'sph')
    

    #changing coordinate system from global to local
    X_diff=GPS_Position_X-LEO_Position_X
    Y_diff=GPS_Position_Y-LEO_Position_Y
    Z_diff=GPS_Position_Z-LEO_Position_Z
    dist=np.sqrt(X_diff**2+Y_diff**2+Z_diff**2)
    lat_LEO=np.deg2rad(LEO_position_geo.lati)
    lon_LEO=np.deg2rad(LEO_position_geo.long)
    
    X_prime=-X_diff*np.sin(lat_LEO)*np.cos(lon_LEO)-Y_diff*np.sin(lat_LEO)*np.sin(lon_LEO)+Z_diff*np.cos(lat_LEO)
    Y_prime=-X_diff*np.sin(lon_LEO)+Y_diff*np.cos(lon_LEO)
    Z_prime=X_diff*np.cos(lat_LEO)*np.cos(lon_LEO)+Y_diff*np.cos(lat_LEO)*np.sin(lon_LEO)+Z_diff*np.sin(lat_LEO)
    
    #GPS look angles
    #Az=np.arctan(Y_prime/X_prime)
    Az=np.arctan2(Y_prime,X_prime)
    Ze=np.arccos(Z_prime/dist)
    El=np.pi/2-Ze
        
    #calculating GPS coordinate on the thin sspherical hell height
    psi=np.pi/2-El-np.arcsin((MEAN_EARTH_RADIUS+height)/(MEAN_EARTH_RADIUS+height+SHELL_HEIGHT)*np.cos(El))
    
    lat_IPP=np.arcsin(np.sin(lat_LEO)*np.cos(psi)+np.cos(lat_LEO)*np.sin(psi)*np.cos(Az))      
    lon_IPP=lon_LEO+np.arcsin(np.sin(psi)*np.sin(Az)/np.cos(lat_IPP))
    
    lat_geo_IPP=np.rad2deg(lat_IPP)
    lon_geo_IPP=np.rad2deg(lon_IPP)   

    
    #magnetic coordinates
    A=Apex(datetime(year[0],month[0],day[0]),np.nanmean(height)+SHELL_HEIGHT)
    
    mlat_GPS=[]
    mlon_GPS=[]
    MLT_GPS=[]
    for i in range(len(GPS_position_geo)):
        mag_lat,mag_lon=A.convert(lat_geo_IPP[i],lon_geo_IPP[i],'geo','qd',height[i]+SHELL_HEIGHT)  
        mlt=A.mlon2mlt(mag_lon,datetime(year[i],month[i],day[i],hour[i],minute[i],second[i])) 
        mlat_GPS.append(mag_lat)
        mlon_GPS.append(mag_lon)
        MLT_GPS.append(mlt)
  
    lat_mag_complete.append(mlat_GPS)
    MLT_complete.append(MLT_GPS)   
    
    #TEC
    rot_parameters=rod_rodi_calculation(TEC_for_ROTI,date_list,WINDOW,FREQUENCY,GAP,FLAG)
    gradient=gradient_calculation(TEC_for_ROTI,date_list,radius+SHELL_HEIGHT*1000,lat_geo_IPP,lon_geo_IPP,FREQUENCY,FLAG)
    
    ROTI_complete.append(rot_parameters[2])
    
    print('    Saving output file...')
    
    filename_output='Swarm_TEC_2Hz_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_PRN'+str(PRN).zfill(2)+'_output.txt'
    f=open(os.path.join(path_output,'data','TEC',filename_output),'w')
    f.write('year          month	       day	  	hour	   	min	   	sec	   	msec	    	doy	      	  hourUT		  hourLT	   	   MLT		    	LatGeo_LEO		LonGeo_LEO		Radius [m]	       Height [km]	         LatMag_QD		 LonMag_QD	      Abs_STEC [TECU]       Abs_VTEC [TECU] 	     Rel_STEC [TECU] 	  Rel_STEC_RMS [TECU]   		Elev_Angle     	PRN         	LatGeo_IPP             LonGeo_IPP          Height_IPP [km]                  ROT [TECU/s]    		ROTI [TECU/s]  	   ROT_mean [TECU/s] 		   TEC_grad [TECU/km]   	   mean_TEC [TECU]\n')
    for i in range(len(prn)):
        f.write('%4i\t\t%2i\t\t%2i\t\t%2i\t\t%2i\t\t%2i\t\t%3i\t\t%3i\t\t%8.5f\t\t%8.5f\t\t%8.5f\t\t%10.5f\t\t%10.5f\t\t%10.2f\t\t%10.2f\t\t%10.5f\t\t%10.5f\t\t%10.5f\t\t%10.5f\t\t%10.5f\t\t%10.5f\t\t\t%7.3f\t\t%2i\t\t%10.5f\t\t%10.5f\t\t%9.2f\t\t\t%9.5f\t\t\t%9.5f\t\t\t%9.5f\t\t\t%9.5f\t\t\t%8.3f\n' % \
         (year[i],month[i],day[i],hour[i],minute[i],second[i],msec[i],doy[i],hourUT[i],hourLT[i],MLT[i],\
          lat_geo[i],lon_geo[i],radius[i],height[i],lat_mag[i],lon_mag[i],\
          abs_sTEC[i],abs_vTEC[i],rel_sTEC[i],rel_sTEC_RMS[i],El[i],prn[i],lat_geo_IPP[i],lon_geo_IPP[i],\
          height[i]+SHELL_HEIGHT,rot_parameters[1][i],rot_parameters[2][i],rot_parameters[3][i],gradient[1][i],rot_parameters[5][i]))
    
    f.close()
    
    
    
    if(FIGURE):
          
        ###############################################################################
        #ROTI figures
        ###############################################################################
        print('    Making output figures...')
        
        #Map scatter lat-lon Global
        fig=plt.figure(figsize=(18,18))
            
        plt.subplot(2,1,1) 
        m=Basemap(projection='robin',lon_0=0,llcrnrlat=-90,urcrnrlat=90,\
         llcrnrlon=-180,urcrnrlon=180,resolution='l')     
        
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,30),labels=[1,1,0,0],fontsize=20)
        m.drawmeridians(np.arange(-180,180,60),labels=[0,0,1,1],fontsize=20)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)
        
        x, y = m(lon_geo_IPP,lat_geo_IPP)#m(GPS_position_geo.long,GPS_position_geo.lati)#m(lon_geo,lat_geo)
        fig=m.scatter(x,y,s=50,c=np.log10(rot_parameters[2]),edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=-3,vmax=-1)
    
        
        plt.subplot(2,1,2)
        plt.scatter(lat_geo_IPP,np.rad2deg(El),s=50,c=np.log10(rot_parameters[2]),edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=-3,vmax=-1)
        plt.ylabel('Elevation Angle (°)', fontsize=20)
        plt.xlabel('Geographic Latitude', fontsize=20)
        plt.grid()
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlim(-90,90)
        plt.ylim(15,95)
        
        cbar=m.colorbar(fig,location='right',pad="5%",size="2.3%")
        cbar.set_label('$log_{10}(ROTI)\ (log_{10}(TECU/s))$', fontsize=20)
        cbar.ax.tick_params(labelsize=20)
 
        plt.figimage(logo,0,0,zorder=-1,alpha=1)
            
        plt.suptitle('Scatter ROTI, Global, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+', PRN='+str(PRN).zfill(2)+'\n',fontsize=30,fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','TEC','ROTI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_PRN'+str(PRN).zfill(2)+'_Global_scatter.png'), format='png',dpi=150)
    
        plt.close()
        plt.clf()
    
    
        #Map scatter lat-lon North Pole
        fig=plt.figure(figsize=(25,10))
        
        plt.subplot(1,2,1) 
        m = Basemap(lon_0=0, boundinglat=40, projection='npstere',round=True,resolution='l')
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,10),labels=[1,1,1,1],fontsize=20)
        m.drawmeridians(np.arange(0,360,30),labels=[1,1,1,1],fontsize=20)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)
        
        x, y = m(lon_geo_IPP,lat_geo_IPP)#m(GPS_position_geo.long,GPS_position_geo.lati)#m(lon_geo,lat_geo)
        fig=m.scatter(x,y,s=50,c=np.log10(rot_parameters[2]),edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=-3,vmax=-1)#,latlon=True)
      
        
        plt.subplot(1,2,2)
        plt.scatter(lat_geo_IPP,np.rad2deg(El),s=50,c=np.log10(rot_parameters[2]),edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=-3,vmax=-1)#,latlon=True)
        plt.ylabel('Elevation Angle (°)', fontsize=20)
        plt.xlabel('Geographic Latitude', fontsize=20)
        plt.grid()
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlim(40,90)
        plt.ylim(15,95)

        cbar=m.colorbar(fig,location='right',pad="5%",size="2.3%")#,format='%3.1e')
        cbar.set_label('$log_{10}(ROTI)\ (log_{10}(TECU/s))$', fontsize=20)
        cbar.ax.tick_params(labelsize=20)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.suptitle('Scatter ROTI, North Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+', PRN='+str(PRN).zfill(2)+'\n',fontsize=30,fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','TEC','ROTI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_PRN'+str(PRN).zfill(2)+'_North_Pole_scatter.png'), format='png',dpi=150)
    
        plt.close()
        plt.clf()
    
    
        
        #Map scatter lat-lon South Pole
        fig=plt.figure(figsize=(25,10))
            
        plt.subplot(1,2,1) 
        m = Basemap(lon_0=0, boundinglat=-40, projection='spstere',round=True,resolution='l')
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,10),labels=[1,1,1,1],fontsize=20)
        m.drawmeridians(np.arange(0,360,30),labels=[1,1,1,1],fontsize=20)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)    
        
        x, y = m(lon_geo_IPP,lat_geo_IPP)#m(GPS_position_geo.long,GPS_position_geo.lati)#m(lon_geo,lat_geo)
        fig=m.scatter(x,y,s=50,c=np.log10(rot_parameters[2]),edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=-3,vmax=-1)#,latlon=True)

        
        plt.subplot(1,2,2)
        plt.scatter(lat_geo_IPP,np.rad2deg(El),s=50,c=np.log10(rot_parameters[2]),edgecolors='None', linewidths=None,cmap=plt.cm.jet,vmin=-3,vmax=-1)#,latlon=True)
        plt.ylabel('Elevation Angle (°)', fontsize=20)
        plt.xlabel('Geographic Latitude', fontsize=20)
        plt.grid()
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlim(-90,-40)
        plt.ylim(15,95)
        
        cbar=m.colorbar(fig,location='right',pad="5%",size="2.3%")#,format='%3.1e')
        cbar.set_label('$log_{10}(ROTI)\ (log_{10}(TECU/s))$', fontsize=20)
        cbar.ax.tick_params(labelsize=20)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.suptitle('Scatter ROTI, South Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+', PRN='+str(PRN).zfill(2)+'\n',fontsize=30,fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','TEC','ROTI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_PRN'+str(PRN).zfill(2)+'_South_Pole_scatter.png'), format='png',dpi=150)
        
        plt.close()
        plt.clf()
    
                  
        ###############################################################################
        #TEC figures
        ###############################################################################    
        #Map scatter lat-lon Global
        fig=plt.figure(figsize=(18,10))
        
        m=Basemap(projection='robin',lon_0=0,llcrnrlat=-90,urcrnrlat=90,\
         llcrnrlon=-180,urcrnrlon=180,resolution='l')     
        
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,30),labels=[1,1,0,0],fontsize=15)
        m.drawmeridians(np.arange(-180,180,60),labels=[0,0,1,1],fontsize=15)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)
               
        x, y = m(lon_geo_IPP,lat_geo_IPP)#m(GPS_position_geo.long,GPS_position_geo.lati)#m(lon_geo,lat_geo)
        
        fig=m.scatter(x,y,s=50,c=TEC_for_ROTI,edgecolors='None', linewidths=None,cmap=plt.cm.jet)#,vmin=0,vmax=150)
        
        cbar=m.colorbar(fig,location='right',pad="5%",size="2.3%")#,format='%3.1e')
        cbar.set_label(label_TEC, fontsize=15)
        cbar.ax.tick_params(labelsize=15)
    
        plt.title('Scatter '+label_TEC[:-7]+', Global, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+', PRN='+str(PRN).zfill(2)+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','TEC',label_TEC[:-7].split()[0]+'_'+label_TEC[:-7].split()[1]+'_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_PRN'+str(PRN).zfill(2)+'_Global_scatter.png'), format='png',dpi=150)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.clf()
        plt.close()
        
        #Map scatter lat-lon North Pole
        fig=plt.figure(figsize=(12,10))
        
        m = Basemap(lon_0=0, boundinglat=40, projection='npstere',round=True,resolution='l')
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,10),labels=[1,1,1,1],fontsize=15)
        m.drawmeridians(np.arange(0,360,30),labels=[1,1,1,1],fontsize=15)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)
        
        x, y = m(lon_geo_IPP,lat_geo_IPP)#m(GPS_position_geo.long,GPS_position_geo.lati)#m(lon_geo,lat_geo)
    
        fig=m.scatter(x,y,s=50,c=TEC_for_ROTI,edgecolors='None', linewidths=None,cmap=plt.cm.jet)#,vmin=0,vmax=50)#,latlon=True)
    
        cbar=m.colorbar(fig,location='right',pad="10%",size="5%")#,format='%3.1e')
        cbar.set_label(label_TEC, fontsize=15)
        cbar.ax.tick_params(labelsize=15)

        plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
        plt.title('Scatter '+label_TEC[:-7]+', North Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+', PRN='+str(PRN).zfill(2)+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','TEC',label_TEC[:-7].split()[0]+'_'+label_TEC[:-7].split()[1]+'_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_PRN'+str(PRN).zfill(2)+'_North_Pole_scatter.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()
        
        #Map scatter lat-lon South Pole
        fig=plt.figure(figsize=(12,10))
        
        m = Basemap(lon_0=0, boundinglat=-40, projection='spstere',round=True,resolution='l')
        m.drawcoastlines(linewidth=1)
        m.drawparallels(np.arange(-90,90,10),labels=[1,1,1,1],fontsize=15)
        m.drawmeridians(np.arange(0,360,30),labels=[1,1,1,1],fontsize=15)
        m.drawcountries()
        m.shadedrelief(scale=0.25,alpha=0.5)  
               
        x, y = m(lon_geo_IPP,lat_geo_IPP)#m(GPS_position_geo.long,GPS_position_geo.lati)#m(lon_geo,lat_geo)
    
        fig=m.scatter(x,y,s=50,c=TEC_for_ROTI,edgecolors='None', linewidths=None,cmap=plt.cm.jet)#,vmin=0,vmax=50)#,latlon=True)
        
     
        cbar=m.colorbar(fig,location='right',pad="10%",size="5%")#,format='%3.1e')
        cbar.set_label(label_TEC, fontsize=15)    
        cbar.ax.tick_params(labelsize=15)
 
        plt.figimage(logo,0,0,zorder=-1,alpha=1)
            
        plt.title('Scatter '+label_TEC[:-7]+', South Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+', PRN='+str(PRN).zfill(2)+'\n\n',fontsize=20, fontweight='bold')    
        plt.savefig(os.path.join(path_output,'figures','TEC',label_TEC[:-7].split()[0]+'_'+label_TEC[:-7].split()[1]+'_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_PRN'+str(PRN).zfill(2)+'_South_Pole_scatter.png'), format='png',dpi=150)
    
        plt.clf()
        plt.close()



###########################################################################
#All PRNs merged together
###########################################################################
if(FIGURE):
    
    print('Working on Swarm '+SAT+' '+YEAR+'/'+MONTH+'/'+DAY+' all PRNs TEC data...')
    
    lat_mag=np.concatenate(lat_mag_complete)
    MLT=np.concatenate(MLT_complete)
    ROTI=np.concatenate(ROTI_complete)
    ###########################################################################
    #Binning in Magnetic coordinates
    ###########################################################################
    print('    Binning data in magnetic coordinates...')
    
    lat_QD_NP=np.ma.masked_less(lat_mag,40)
    MLT_map_NP=np.ma.masked_array(MLT,lat_QD_NP.mask)#/24.*(2*np.pi)
    lat_QD_NP=-1.*(lat_QD_NP-90.)
    
    for i in range(len(lat_QD_NP)):
        if(lat_QD_NP[i]>50):
            lat_QD_NP[i]=np.nan
            MLT_map_NP[i]=np.nan
           
    binned_roti_QD_MLT_NP=binning_2D(MLT_map_NP,lat_QD_NP,ROTI,0,24,0.25,0,50,2)  
        
    
    lat_QD_SP=np.ma.masked_greater(lat_mag,-40)
    MLT_map_SP=np.ma.masked_array(MLT,lat_QD_SP.mask)#/24.*(2*np.pi)
    lat_QD_SP=lat_QD_SP+90.
    
    for i in range(len(lat_QD_SP)):
        if(lat_QD_SP[i]>50):
            lat_QD_SP[i]=np.nan
            MLT_map_SP[i]=np.nan
            
    binned_roti_QD_MLT_SP=binning_2D(MLT_map_SP,lat_QD_SP,ROTI,0,24,0.25,0,50,2)  
        
    binned_roti_QD_MLT=binning_2D(MLT,lat_mag,ROTI,0,24,0.25,-90,90,2)    
    
     
    ###########################################################################
    #QD-MLT polar plots
    ###########################################################################
    
    #North Pole
    fig=plt.figure(figsize=(12,12))
    
    rad=np.linspace(0,50,26, endpoint='True')
    azm=np.linspace(0,24,97, endpoint='True')/24.*(2*np.pi)
    r, th = np.meshgrid(rad, azm)
    z=np.log10(binned_roti_QD_MLT_NP[0].T)
    
    ax=plt.subplot(projection="polar")
    plt.pcolormesh(th-np.pi/2., r, z,cmap=plt.cm.jet,vmin=-3,vmax=-1)
    plt.plot(azm, r, color='k', ls='none') 
    plt.grid()
    plt.xticks(np.linspace(0,2*np.pi,24,endpoint=False),(np.concatenate((np.linspace(6,23,18,dtype=int),np.linspace(0,5,6,dtype=int)))),fontsize=15)
    
    ax.set_rorigin(0)
    ax.set_theta_zero_location('E', offset=0)
    plt.yticks(np.linspace(0,40,5),['90','80','70','60','50'],fontsize=15,alpha=1,color='black')
    
    plt.xlabel('MLT-Magnetic Latitude ', fontsize=20)
    cbar=plt.colorbar(orientation='horizontal',fraction=0.05)
    cbar.set_label('$log_{10}$(ROTI) ($log_{10}$(TECU/s))', fontsize=15)
    
    cbar.ax.tick_params(labelsize=15)

    plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
    plt.title('Magnetic polar plot ROTI, North Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n',fontsize=20, fontweight='bold')    
    plt.savefig(os.path.join(path_output,'figures','TEC','ROTI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_all_PRNs_North_Pole_magnetic_polar_plot.png'), format='png',dpi=150)
    
    plt.clf()
    plt.close()
    
    
    #South Pole
    fig=plt.figure(figsize=(12,12))
    
    rad=np.linspace(0,50,26, endpoint='True')
    azm=np.linspace(0,24,97, endpoint='True')/24.*(2*np.pi)
    r, th = np.meshgrid(rad, azm)
    z=np.log10(binned_roti_QD_MLT_SP[0].T)
    
    ax=plt.subplot(projection="polar")
    plt.pcolormesh(th-np.pi/2., r, z,cmap=plt.cm.jet,vmin=-3,vmax=-1)
    plt.plot(azm, r, color='k', ls='none') 
    plt.grid()
    plt.xticks(np.linspace(0,2*np.pi,24,endpoint=False),(np.concatenate((np.linspace(6,23,18,dtype=int),np.linspace(0,5,6,dtype=int)))),fontsize=15)
    
    ax.set_rorigin(0)
    ax.set_theta_zero_location('E', offset=0)
    plt.yticks(np.linspace(0,40,5),['-90','-80','-70','-60','-50'],fontsize=15,alpha=1,color='black')
    
    plt.xlabel('MLT-Magnetic Latitude ', fontsize=20)
    cbar=plt.colorbar(orientation='horizontal',fraction=0.05)
    cbar.set_label('$log_{10}$(ROTI) ($log_{10}$(TECU/s))', fontsize=15)
    
    cbar.ax.tick_params(labelsize=15)

    plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
    plt.title('Magnetic polar plot ROTI, South Pole, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n',fontsize=20, fontweight='bold')    
    plt.savefig(os.path.join(path_output,'figures','TEC','ROTI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_all_PRNs_South_Pole_magnetic_polar_plot.png'), format='png',dpi=150)
    
    plt.clf()
    plt.close()
    
    
    #Global        
    fig=plt.figure(figsize=(18,10))
    
    y=np.linspace(-90,90,91, endpoint='True')
    x=np.linspace(0,24,97,endpoint='True')
    X,Y = np.meshgrid(x,y)
    z=np.log10(binned_roti_QD_MLT[0])
    
    plt.pcolormesh(X,Y,z,cmap=plt.cm.jet,vmin=-3,vmax=-1)
    plt.grid()
    plt.xticks(np.linspace(0,24,25,endpoint=True),fontsize=15)
    plt.yticks(np.linspace(-90,90,7),fontsize=15)
    plt.xlabel('MLT', fontsize=20)
    plt.ylabel('Magnetic Latitude', fontsize=20)
    cbar=plt.colorbar(orientation='vertical',fraction=0.05)
    cbar.set_label('$log_{10}$(ROTI) ($log_{10}$(TECU/s))', fontsize='xx-large')
    cbar.ax.tick_params(labelsize=15)

    plt.figimage(logo,0,0,zorder=-1,alpha=1)
        
    plt.title('Magnetic coordinates plot ROTI, Global, Swarm '+SAT+', '+YEAR+'/'+MONTH+'/'+DAY+'\n',fontsize=20, fontweight='bold')    
    plt.savefig(os.path.join(path_output,'figures','TEC','ROTI_Swarm_'+SAT+'_'+YEAR+'_'+MONTH+'_'+DAY+'_all_PRNs_Global_magnetic_coordinates.png'), format='png',dpi=150)
    
    plt.clf()
    plt.close()
    
   
    ############################################################################### 
    


