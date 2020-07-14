##### HERE WERE DOING IMPORTS
import os
import sys
import glob
import datetime as dt
from datetime import datetime
from datetime import timedelta
from dateutil.relativedelta import relativedelta

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
##from matplotlib import style
import numpy as np
import numpy.ma as ma
import pandas as pd
import scipy
import scipy.signal 
import spacepy
from spacepy.pycdf import CDF
from spacepy.plot import applySmartTimeTicks
from spacepy import omni, time, pybats
from spacepy import pybats
import spacepy.plot
from spacepy.plot import add_arrows, style
import spacepy.irbempy as ib
import spacepy.time as spt
import spacepy.coordinates as spc
import tailwag as tw
#######################################################################################

#### Prepare to plot:
# Create an output directory:
##outdir = 'fuze_plots/'
##if not os.path.exists(outdir): os.mkdir(outdir)



###################################################################################   

def main():
    print("THIS IS THE BEGINNING OF THE PROGRAM:  \n")   
    print(spacepy.plot.available(returnvals = False))
    
    
    
    
    ####  swiping the event points from the excel file that were gonna use
    Event_Data = pd.read_excel("Event_Points.xlsx", header = 0)


    ####  getting three arrays: Point_list as strings, Date_List as datetimes, and the used cluster satellite
    Point_List = Event_Data['Narrowed Point']
    Date_List  = [dt.datetime.strptime(date,'%Y-%m-%d %H:%M:%S.%f') for date in Point_List]
    Sat_List   = Event_Data['USED SAT']
    
       
    for (a, b) in zip(Date_List, Sat_List):
        fig, ax1, ax2, ax3, ax4, ax5, n = tw.fusion(a, b, 12, outdir = 'fusion_plots/')
        
        
        
#########################################################################################################
# If this file is used via IPython's "run" magic command,
# Evaluate on all figures.
if __name__ == '__main__':
    main()
       
    
    
    

    
 
