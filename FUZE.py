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
import numpy as np
import numpy.ma as ma
import pandas as pd
import scipy
import scipy.signal
import spacepy
from spacepy.pycdf import CDF
from spacepy.plot import applySmartTimeTicks, style
from spacepy import omni, time, pybats
from spacepy import pybats
import spacepy.plot
from spacepy.plot import add_arrows
import spacepy.irbempy as ib
import spacepy.time as spt
import spacepy.coordinates as spc
import tailwag as tw
#######################################################################################

#### Prepare to plot:
# Create an output directory:
outdir = 'fuze_plots/'
if not os.path.exists(outdir): os.mkdir(outdir)



###################################################################################   

def main():
    print("THIS IS THE BEGINNING OF THE PROGRAM:  \n")   
    
    
    
    
    ####  swiping the event points from the excel file that were gonna use
    Event_Data = pd.read_excel("Event_Points.xlsx", header = 0)


    ####  getting three arrays: Point_list as strings, Date_List as datetimes, and the used cluster satellite
    Point_List = Event_Data['Narrowed Point']
    Date_List  = [dt.datetime.strptime(date,'%Y-%m-%dT%H:%M:%S.%f') for date in Point_List]
    Sat_List   = Event_Data['USED SAT']
    
    ##for (a, b) in zip(Point_List, Date_List):
        ##print(type(a), type(b))
    
    
    

    
    
    for (a, b, c) in zip(Date_List, Sat_List, Point_List):
    ###  Print to see which date we are currently using
        print("Entering the for loop with the following event date: " + c + "\n\n")
    
    
    
        #####   Here were grabbing our cluster data that we will need for this event
        #####   If we need several months the function will swipe that shit so dont worry boomer
        #####   change that last number to be the amount of hours you want on each side of your interval
        Data_Dict = tw.fetch_cluster_data(a, b, 12)
        
        
        # FETCH TSYG DATA HERE!
        
        TSYG_T89_time, TSYG_T89_B  = tw.gen_sat_tsyg(Data_Dict, extMag='T89', dbase='qd1min');  
        TSYG_T96_time, TSYG_T96_B = tw.gen_sat_tsyg(Data_Dict, extMag='T96', dbase='qd1min');  
        TSYG_T01_time, TSYG_T01_B = tw.gen_sat_tsyg(Data_Dict, extMag='T01STORM', dbase='qd1min');  
      
        
        
        print("#############################################################")
        print("The following stuff is what is inside the T89 function......\n")
        ##for thing in TSYG_T89:
            ##print(thing)
    
    
        ####LOOP THAT WILL ITERATE 3 TIMES TO GIVE US OUR 3 GRAPHS PER EPOCH........#####
        for n in range(2, 13, 2):
            print("WERE NOW INSIDE THE LOOPING FUNCTION TO SNAG VALUES......\n")
            ####  Set up start and stop times so that we can narrow our data down early
            start_Time = a - timedelta(hours = n)
            end_Time   = a + timedelta(hours = n)
            #### Set up a filter to narrow down the stuff from the massive dictionary 
            #### That we just made
            filter1 = (Data_Dict['cis_time']>=start_Time) & (Data_Dict['cis_time']<=end_Time)    
            CIS_time = Data_Dict['cis_time'][filter1]
            maskpro = Data_Dict['dens_h'][filter1]
            maskoxy = Data_Dict['dens_o'][filter1]
            
            filter2 = (Data_Dict['fgm_time']>=start_Time) & (Data_Dict['fgm_time']<=end_Time)
            FGM_time = Data_Dict['fgm_time'][filter2]
            Bx = Data_Dict['b'][filter2,0]
            x = Data_Dict['xyz'][filter2,0]
            y = Data_Dict['xyz'][filter2,1]
            z = Data_Dict['xyz'][filter2,2]

            
            
            
            T89_time = TSYG_T89_time[filter2]
            T89_Bx = TSYG_T89_B[filter2, 0]
            
            T96_time = TSYG_T96_time[filter2]
            T96_Bx = TSYG_T96_B[filter2, 0]
            
            T01_time = TSYG_T01_time[filter2]
            T01_Bx = TSYG_T01_B[filter2, 0]
            
            
            print("WERE NOW HITTING DAT FOR LOOP.....\n")
            # FILTER TSYG DATA HERE.
            ##for (a, b) in zip(T89_time,T89_Bx):
                ##print (a, b)
            
            
            

            

            ind_x = x[round(len(FGM_time)/2)]
            ind_y = y[round(len(FGM_time)/2)]
            ind_z = z[round(len(FGM_time)/2)]
        
            print("The ind_x is:   " + str(ind_x) + "\n")
            print("The ind_y is:   " + str(ind_y) + "\n")
            
            print("The length of the FGM array is:   " + str(len(FGM_time)) + "\n")
            print("The length of half FGM array is:   " + str(round(len(FGM_time)/2)) + "\n")
            
             

###################################################################################
#                          PLOTTING SECTION
###################################################################################
        
            fig = plt.figure(figsize=(8.5,11))
            fig.suptitle(c + " with interval hours: " +  str(n), fontsize=16)
        
        
            ####Creating the 4 subplots were going to need
            ax1, ax2 = fig.add_subplot(321), fig.add_subplot(322)
            ax3, ax4 = fig.add_subplot(312), fig.add_subplot(313)

            # Axes 1: Orbit in x-y plane:
            line1 = ax1.plot(x, y)
            add_arrows(line1, n = 5, size = 18, style = '->')
        
            ax1.set(xlabel='X in R$_E$', ylabel='Y in R$_E$',
                    title=' X vs Y position')
            ax1.annotate('Cross', xy=(ind_x,ind_y), xytext = (ind_x + 0.2, ind_y),
                         arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),)  
            
            # Axes 2: Orbit in x-z plane:
            line2 = ax2.plot(x, z)
            ax2.set(xlabel='X in R$_E$', ylabel='Z in R$_E$',
                    title=' X vs Z position')
            add_arrows(line2, n = 5, size = 18, style = '->')
            ax2.annotate('Cross', xy=(ind_x,ind_z), xytext = (ind_x + 0.2, ind_z),
                         arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),) 
    
            # Axes 3: Magnetic field x-component:
            ax3.plot(FGM_time, Bx, lw = .5)        
            ax3.set(xlabel='time', ylabel='$B_X$ in nT',
                    title='$B_X$ Magnitude')
            ##ax3.set_xlim([start_Time, end_Time])
            ax3.hlines(0, FGM_time[0], FGM_time[-1], lw=1, color='black', linestyle = '--')
            ax3.axvline(x = a, ymin = 0, ymax = 1, lw=1, color = 'black', linestyle = '--')
            applySmartTimeTicks(ax3, [start_Time, end_Time])
            
            ### Extra Plots### Magnetic field x-compnent from T89, T96, T01STORM
            ax3.plot(T89_time, T89_Bx, color = 'r', lw=.5)
            ax3.plot(T96_time, T96_Bx, color = 'g', lw=.5)
            ax3.plot(T01_time, T01_Bx, color = 'y', lw=.5)
            
            
            
            
           
          
            # Axes 4: Composition - Oxygen
            ax4.set_xlabel('time')
            ax4.set_title('Proton and Oxygen Densities')
        
            ax4.plot(CIS_time, maskoxy, color = 'g', lw = .5, alpha=0.7)
            ax4.set_ylabel('Oxygen per $cm^{3}$', color = 'g')        
            ## ax4.set_xlim([start_Time, end_Time])
            ax4.grid(b = None, which='major', axis='both', color = 'g', alpha= 0.2, linestyle = '--')
            ax4.axvline(x = a, ymin = 0, ymax = 1, lw=1, color = 'black', linestyle = '--')

            # Axes 4: Composition - Protons
            ax5 = ax4.twinx()         
            ax5.plot(CIS_time, maskpro, color = 'r', lw=.5)
            ax5.set_ylabel('Protons per $cm^{3}$', color = 'r')          
            ##ax5.set_xlim([start_Time, end_Time])               
            ax5.grid(b = None, which='major', axis='both', color = 'r', linestyle = '--', alpha= 0.2 )
            applySmartTimeTicks(ax4, [start_Time, end_Time])  
    
        
            ##### With these 2 blocks were adding planet earth to our postion graphs
            ax1 = fig.add_subplot(321)
            ax2 = fig.add_subplot(322)
       
            spacepy.pybats.add_planet(ax1)
            spacepy.pybats.add_planet(ax2)
        
        
            #### Here were making sure everything stays nice and cool looking
            plt.tight_layout(rect=[0, 0, 1, .95])
            plt.show()
        
        
            ####  HERE WERE SAVING THE GRAPHS TO A FOLDER AND GIVING THE FILE A TITLE THAT HAVE THE DATE
            #####  AND THE INTERVAL HOURS USED
            fig.savefig(outdir + f'fuze_T{a:%Y%m%d_%H%M%S}_int{n:02d}.png')
            
        
#########################################################################################################
# If this file is used via IPython's "run" magic command,
# Evaluate on all figures.
if __name__ == '__main__':
    main()  
