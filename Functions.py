# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 11:34:33 2020

@author: pigna
"""


import numpy as np

###############################################################################
#Definition of functions used
###############################################################################       


def rod_rodi_calculation(Ne_list,date_list,window_lenght,frequency,gap,flag):
    rod_list=[]
    rodi_list=[]
    index_list=[]
    rod_mean_list=[]
    rod_abs_mean_list=[]
    Ne_mean_list=[]

    for ijk in range(len(date_list)):
        #RODI calculation
        index_window_center=ijk
        window_center=date_list[ijk]

        if(ijk<int(window_lenght/2) or ijk>len(date_list)-(int(window_lenght/2)+1)):
            indices_values_windowed=[flag]
        else:            
            indices_values_windowed=[]
            for i in range(-(int(window_lenght/2)),(int(window_lenght/2))+1):
                if(np.abs(date_list[ijk+i]-window_center)<=(int(window_lenght/2))/frequency):
                    indices_values_windowed.append(ijk+i)
                    
        if((window_lenght-len(indices_values_windowed))<=gap):            
            rod=[]
            Ne=[]
            for j in range(len(indices_values_windowed)-1):
                delta_time=date_list[indices_values_windowed[j+1]]-date_list[indices_values_windowed[j]]
                delta_Ne=Ne_list[indices_values_windowed[j+1]]-Ne_list[indices_values_windowed[j]]
                Ne.append(Ne_list[indices_values_windowed[j]])
                if(delta_time==(1./frequency)):
                    rod.append(float(delta_Ne)/float(delta_time))
            
            rod_mean=np.nanmean(rod)
            rodi=np.nanstd(rod,ddof=1)
            rod_abs_mean=np.nanmean(np.abs(rod))
            
            rodi_list.append(rodi)
            index_list.append(index_window_center)
            rod_mean_list.append(rod_mean)
            rod_abs_mean_list.append(rod_abs_mean)
            Ne_mean_list.append(np.nanmean(Ne))
        else:
            rodi_list.append(flag)
            rod_mean_list.append(flag)
            index_list.append(index_window_center)
            rod_abs_mean_list.append(flag)
            Ne_mean_list.append(flag)
        
        #ROD calculation
        if(ijk<len(date_list)-1):
            delta_time=date_list[ijk+1]-date_list[ijk]
            delta_Ne=Ne_list[ijk+1]-Ne_list[ijk]
            if(delta_time==(1./frequency)):
                rod_list.append(float(delta_Ne)/float(delta_time))
            else:
                rod_list.append(flag)
            
        else:   
            rod_list.append(flag)
        

    
    return (index_list,rod_list,rodi_list,rod_mean_list,rod_abs_mean_list,Ne_mean_list)



def spherical_distance(radius,lat1,lon1,lat2,lon2):
    dLat=np.radians(lat2-lat1);
    dLon=np.radians(lon2-lon1);
    lat1=np.radians(lat1);
    lat2=np.radians(lat2);
    a=np.sin(dLat/2)**2+np.sin(dLon/2)**2*np.cos(lat1)*np.cos(lat2)
    c = 2*np.arctan2(np.sqrt(a),np.sqrt(1-a))
    
    return radius*c



def gradient_calculation(Ne_list,date_list,radius,lat_geo,lon_geo,frequency,flag):
    gradient=[]
    index_gradient=[]
    for i in range(len(Ne_list)):
        if(i<len(Ne_list)-1):            
            if(date_list[i+1]-date_list[i]==1./frequency):
                delta_Ne=Ne_list[i+1]-Ne_list[i]
                delta_distance=spherical_distance(radius[i]/1000.,lat_geo[i],lon_geo[i],lat_geo[i+1],lon_geo[i+1])
                gradient.append(float(delta_Ne)/float(delta_distance))
                index_gradient.append(i)
            else:
                gradient.append(flag)
                index_gradient.append(i)
        else:
            gradient.append(flag)
            index_gradient.append(i)
    
    return (index_gradient,gradient)



def init_list_of_objects(size):
    list_of_objects = list()
    for i in range(0,size):
        list_of_objects.append( list() ) 
    return list_of_objects


def binning_2D(binning_fun_x,binning_fun_y,fun_to_bin,bin_min_x,bin_max_x,step_bin_x,bin_min_y,bin_max_y,step_bin_y):
    num_bin_x=int((bin_max_x-bin_min_x)/step_bin_x)
    num_bin_y=int((bin_max_y-bin_min_y)/step_bin_y)
    binned=init_list_of_objects(num_bin_x*num_bin_y)
    counter_bin=num_bin_x*num_bin_y*[np.nan]
    for j in range(len(fun_to_bin)):
        if (not (np.isnan(binning_fun_x[j]) or np.isnan(binning_fun_y[j]) or np.ma.is_masked(binning_fun_x[j]) or np.ma.is_masked(binning_fun_y[j]))):
            index_x=int((binning_fun_x[j]-bin_min_x)/step_bin_x)
            index_y=int((binning_fun_y[j]-bin_min_y)/step_bin_y)
            if (not np.isnan(fun_to_bin[j])):
                binned[index_x+index_y*num_bin_x].append(fun_to_bin[j])

    for k in range(num_bin_y):
        for i in range(num_bin_x):
            counter_bin[i+k*num_bin_x]=len(binned[i+k*num_bin_x])
    
    binned=np.array(binned)
    binned=np.reshape(binned,(num_bin_y,num_bin_x))
    counter_bin=np.array(counter_bin)
    counter_bin=np.reshape(counter_bin,(num_bin_y,num_bin_x))
    binned_mean=np.empty((num_bin_y,num_bin_x))*np.nan
    binned_dev=np.empty((num_bin_y,num_bin_x))*np.nan
    binned_median=np.empty((num_bin_y,num_bin_x))*np.nan
    binned_upper_perc=np.empty((num_bin_y,num_bin_x))*np.nan
    binned_lower_perc=np.empty((num_bin_y,num_bin_x))*np.nan
    for i in range(num_bin_y):
        for j in range(num_bin_x):
            binned_mean[i][j]=np.nanmean(MAD(binned[i][j]))
            binned_dev[i][j]=np.nanstd(MAD(binned[i][j]))
            binned_median[i][j]=np.nanmedian(MAD(binned[i][j]))
            binned_upper_perc[i][j]=np.nanpercentile(MAD(binned[i][j]),75)
            binned_lower_perc[i][j]=np.nanpercentile(MAD(binned[i][j]),25)
    return (binned_mean,binned_dev,counter_bin,binned_median,binned_lower_perc,binned_upper_perc)





def MAD(array):
    SIGMA_MAD=1
    median=np.nanmedian(array)
    mad_values=[]
    for i in range(len(array)):
        mad_values.append(np.abs(array[i]-median))
    MAD=1.4826*np.nanmedian(mad_values)
           
    for i in range(len(array)):
        if(array[i]>= median-SIGMA_MAD*MAD and array[i]<= median+SIGMA_MAD*MAD):
            continue
        else:
            array[i]=np.nan
                 
    return array