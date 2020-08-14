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
import scatter
import math



#######################################################################################

#### Prepare to plot:
# Create an output directory:
##outdir = 'fuze_plots/'
##if not os.path.exists(outdir): os.mkdir(outdir)



###################################################################################   

def main():
    print("THIS IS THE BEGINNING OF THE PROGRAM:  \n")   
    print(spacepy.plot.available(returnvals = False))
    
    print()
    
    
    ####  swiping the event points from the excel file that were gonna use
    Event_Data = pd.read_excel("Event_Points.xlsx", header = 0)
    

    
    
    
   


    ####  getting three arrays: Point_list as strings, Date_List as datetimes, and the used cluster satellite
    Point_List = Event_Data['Narrowed Point']
    Date_List  = [dt.datetime.strptime(date,'%Y-%m-%d %H:%M:%S.%f') for date in Point_List]
    Sat_List   = Event_Data['USED SAT']
    Op_List    = Event_Data['Operation']
   
   
    ##for (a, b, c) in zip(Date_List, Sat_List, Op_List):
       ## if c == "SKIP":
        ##    continue
      ##  else:======
          ##  fig, ax1, ax2, ax3, ax4, ax5, n = tw.fusion(a, b, 12, outdir = 'fusion_plots/')
    
    
    
    
    ## SET UP THE ARRAY THAT WERE GOING TO NEED FOR PLOTTING
    Vz_89Values       = []
    Vz_96Values       = []
    Vz_01Values       = []
    T89_diff_Values   = []
    T96_diff_Values   = []
    T01_diff_Values   = []
    
    
        
    
    for (a, b, c) in zip(Date_List, Sat_List, Op_List):
         if c == "SKIP":
            continue
         else:
             print(type(a))
             print(a)   
             
             
             ##### STEP 1 GET ALL THE DATA THAT WERE GONNA NEED ############
             ### FETCH TSYG DATA HERE!  
              ######## STEP 2 START APPENDING THE VALUES TO THE ARRAYS #############
             my_clus, my_ttimes, bork_after, bork_before = scatter.get_crossing_info(a, debug=False)
             Vz_val = tw.Vz_finder_range(str(my_clus[0]), my_clus[0], 2)
             print("bork")
             
             
             print("BRO THIS IS THE cluster date:  " + str(my_clus[0]))
             print("YO BRO THIS IS THE DATE FROM T89:  " + str(my_ttimes['T89']))
             print("YO BRO THIS IS THE DATE FROM T96:  " + str(my_ttimes['T96']))
             print("YO BRO THIS IS THE DATE FROM T01storm:  " + str(my_ttimes['T01STORM']))
             
             if isinstance(my_ttimes['T89'], dt.datetime):
                 print("AWWWEEE SHEEEET DAWG 89 is A DATETIME")
                 timediff89 = my_clus[0] - my_ttimes['T89']
                 minutediff89 = int(timediff89.total_seconds() / 60)                
                 print("The time difference (in minutes) for T89 is:   " + str(minutediff89))
                 T89_diff_Values.append(minutediff89)
                 Vz_89Values.append(Vz_val)
                 
             else:
                 print("THIS IS ACTUALLY A NAN")
                 
             if isinstance(my_ttimes['T96'], dt.datetime):
                 print("AWWWEEE SHEEEET DAWG 89 is A DATETIME")
                 timediff96 = my_clus[0] - my_ttimes['T96']
                 minutediff96 = int(timediff96.total_seconds() / 60)
                 print("The time difference (in minutes) for T89 is:   " + str(minutediff96))
                 T96_diff_Values.append(minutediff96)
                 Vz_96Values.append(Vz_val)
                 
             else:
                 print("THIS IS ACTUALLY A NAN")
    
             if isinstance(my_ttimes['T01STORM'], dt.datetime):
                 print("AWWWEEE SHEEEET DAWG 89 is A DATETIME")
                 timediff01 = my_clus[0] - my_ttimes['T01STORM']
                 minutediff01 = int(timediff01.total_seconds() / 60)          
                 print("The time difference (in minutes) for T89 is:   " + str(minutediff01))
                 T01_diff_Values.append(minutediff01)
                 Vz_01Values.append(Vz_val)
                 
                 
             else:
                 print("THIS IS ACTUALLY A NAN")
    
    
    
    for (a, b) in zip(Vz_89Values, T89_diff_Values):
        print(str(a) + "          " + str(b))
        
    ######## ###########STEP 3 LETS MAKE THAT GRAPH!!!!!! ###############################
             
   
    plt.scatter(Vz_89Values, T89_diff_Values, color='green',  label='T89', alpha=0.5)    
    plt.scatter(Vz_96Values, T96_diff_Values, color='orange', label='T96', alpha=0.5)
    plt.scatter(Vz_01Values, T01_diff_Values, color='crimson',label='T01', alpha=0.5)
    
    trend89 = np.polyfit(Vz_89Values, T89_diff_Values, 1)
    trendpoly89 = np.poly1d(trend89) 
    plt.plot(Vz_89Values,trendpoly89(Vz_89Values), color = 'green', lw=.5)
    
    trend96 = np.polyfit(Vz_96Values, T96_diff_Values, 1)
    trendpoly96 = np.poly1d(trend96) 
    plt.plot(Vz_96Values,trendpoly96(Vz_96Values), color = 'orange', lw=.5)
    
    trend01 = np.polyfit(Vz_01Values, T01_diff_Values, 1)
    trendpoly01 = np.poly1d(trend01) 
    plt.plot(Vz_01Values,trendpoly01(Vz_01Values), color = 'crimson', lw=.5)
    
    
    
    
    plt.title('Vz and time Differential graph')
    plt.xlabel('Vz Values')
    plt.ylabel('Time differentials')
    plt.legend(loc='best')
    plt.show()
             
             
           


             
             
        
        
        
        
        
  
 
    
        

#########################################################################################################
# If this file is used via IPython's "run" magic command,
# Evaluate on all figures.
if __name__ == '__main__':
    main()
       
    
    
     
    
  
    
    

    
 
