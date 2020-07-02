# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 16:21:27 2020

@author: pigna
"""

import os
from datetime import datetime, timedelta
import warnings

warnings.filterwarnings("ignore") 

###############################################################################
#Terminal interface for input parameters
###############################################################################

today = datetime.today()
last_day=today-timedelta(days=6)
last_day_string = last_day.strftime("%d/%m/%Y")

print('**********Input time and satellite**********')
print('Valid dates are from 05/12/2013 to '+last_day_string)


YEAR=int(input('Insert year: '))

if(YEAR<2013 or YEAR>today.year+1):
    raise RuntimeError('No valid year!')

MONTH=int(input('Insert month: '))
if(MONTH<1 or MONTH>12):
    raise RuntimeError('No valid month!')

DOM=int(input('Insert day of the month: '))
try:
    datetime(YEAR,MONTH,DOM)
except:
    raise RuntimeError('No valid day!')

SAT=str(input('Insert Swarm satellite (A,B,C): ')).upper()
if(not (SAT in ('A','B','C'))):
    raise RuntimeError('No valid satellite!')


###############################################################################
#Reading of input parameters file
###############################################################################

filename_parameters=os.path.join(os.getcwd(),'TITIPy_input_parameters.txt')
file=open(filename_parameters,'r')

lines=list(file)    

FIGURE=str.rsplit(lines[3])[0].upper()
if(not (FIGURE in ('Y','N'))):
    raise RuntimeError('No valid input for plotting!')
if(FIGURE=='Y'):
    FIGURE=True
else:
    FIGURE=False

WINDOW_WIDTH=int(str.rsplit(lines[4])[0])
if(datetime(YEAR,MONTH,DOM)<datetime(2014,7,15)):
    FREQUENCY_CHANGE_TEC=True
else:
    FREQUENCY_CHANGE_TEC=False

TEC_INPUT=int(str.rsplit(lines[5])[0])
if(not (TEC_INPUT in (1,2,3))):
    raise RuntimeError('No valid TEC input for ROTI calculation!')

file.close()


   