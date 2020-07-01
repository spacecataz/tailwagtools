


        



    

##### HERE WERE DOING IMPORTS
import datetime as dt
from datetime import datetime
from datetime import timedelta
from dateutil.relativedelta import relativedelta
import glob
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import numpy.ma as ma
import os
import pandas as pd
import scipy
import scipy.signal
import spacepy
from spacepy.pycdf import CDF
import tailwag as tw

from spacepy.plot import applySmartTimeTicks, style

from spacepy import omni, time
from spacepy import pybats

    
###################################################################################   
print("THIS IS THE BEGINNING OF THE PROGRAM:  \n")   


##### READ THE DATA FROM THE EXCEL FILE
Event_Data = pd.read_excel("Event_Points.xlsx", header = 0)

for line in Event_Data:
    print(line)

    
       
##### HERE WE GRAB THE DATA THAT WERE GOING TO NEED FROM THE EXCEL COLUMNS 
##### The first two lines grab the dates and convert them to datetime objects
##### The last two lines grab the used cluster satellites
Point_List = Event_Data['Narrowed Point']
Date_List = [dt.datetime.strptime(date,'%Y-%m-%dT%H:%M:%S.%f') for date in Point_List]
Sat_List = Event_Data['USED SAT']



## Print Statements to check the var types in the arrays:
## Date list is datetime.datetime, point list is string
## CIS and FGM are both int 
'''
for (a,b, c) in zip(Date_List, Sat_List):

    print(type(a), type(b)    
''' 
        


    
        
print("We are now starting to for loop to iterate through the Date_List array: \n")

### This loop will go through datimes in the excel file and grab the satellite
### and satellite that we will pull data from
for (a, b, c) in zip(Date_List, Sat_List, Point_List):
    ###  Print to see which date we are currently using
    print(a)
    
    #### get all of the stuff from these cluster satellites and turn it into a dicrtionary
    #### a feeds the date, b feeds the satellite integer
    #### 12 is a time interval on each side of the used date
    Data_Dict = tw.fetch_cluster_data(a, b, 12)
    
    

    
    
    ####LOOP THAT WILL ITERATE 3 TIMES TO GIVE US OUR 3 GRAPHS PER EPOCH........#####
    for n in range(2, 13, 2):
    
    
    
    
    
  
        ####  Set up start and stop times so that we can narrow our data down early
        start_Time = a - timedelta(hours = n)
        end_Time   = a + timedelta(hours = n)
    
        
        
        
    
    
     

        #### Set up a filter to narrow down the stuff from the massive dictionary 
        #### That we just made
        filter1 = (Data_Dict['cis_time']>=start_Time) & (Data_Dict['cis_time']<=end_Time)    
        CIS_time = Data_Dict['cis_time'][filter1]
        maskpro = ma.masked_values(scipy.signal.medfilt(Data_Dict['dens_h'][filter1], kernel_size = 11), -1.00000E+31)
        maskoxy = ma.masked_values(scipy.signal.medfilt(Data_Dict['dens_o'][filter1], kernel_size = 11), -1.00000E+31)
    
    
    
        #### 
        filter2 = (Data_Dict['fgm_time']>=start_Time) & (Data_Dict['fgm_time']<=end_Time)
        FGM_time = Data_Dict['fgm_time'][filter2]
        Bx = Data_Dict['b'][filter2,0]
        x = Data_Dict['xyz'][filter2,0]/6378.16
        y = Data_Dict['xyz'][filter2,1]/6378.16
        z = Data_Dict['xyz'][filter2,2]/6378.16
    
    
   
        # Convert to ticktocks:
        ticks = time.Ticktock(FGM_time, 'ISO')
       # Fetch minute-level omni:
        omni_data = omni.get_omni(ticks, dbase='qd1min') 
        
        for thing in omni_data:
            print(thing)
        
      
        
###################################################################################
#                          PLOTTING SECTION
###################################################################################
        
        fig = plt.figure(figsize=(8.5,11))
        fig.suptitle(c + " with interval hours: " +  str(n), fontsize=16)
        
        
        ####Creating the 4 subplots were going to need
        ax1, ax2 = fig.add_subplot(321), fig.add_subplot(322)
        ax3, ax4 = fig.add_subplot(312), fig.add_subplot(313)
        
        ax1.plot(x, y)
        ax1.set(xlabel='X in R$_E$', ylabel='Y in R$_E$',
                title=' X vs Y position')
            
        ax2.plot(x, z)
        ax2.set(xlabel='X in R$_E$', ylabel='Z in R$_E$',
                title=' X vs Z position')
    
        ax3.plot(FGM_time, Bx)        
        ax3.set(xlabel='time', ylabel='$B_X$ in nT',
                title='$B_X$ Magnitude')
        ##ax3.set_xlim([start_Time, end_Time])
        ax3.hlines(0, FGM_time[0], FGM_time[-1], lw=2, color='black', linestyle = '--')
        applySmartTimeTicks(ax3, [start_Time, end_Time])
        
        ##ax3b = ax3.twinx()
        ##ax3b.plot(FGM_time, omni_data, color = 'r', lw=1)
        ##ax3b.set_ylabel('Protons per $cm^{3}$')          
        ##ax3b.set_xlim([start_Time, end_Time])          
       ## ax3b.grid(b = None, which='major', axis='both')
        
        
        ax4.set_xlabel('time')
        ax4.set_title('Proton and Oxygen Densities')
        
        ax4.plot(CIS_time, maskoxy, 'g-', lw=1, alpha=0.7)
        ax4.set_ylabel('Oxygen per $cm^{3}$')        
        ax4.set_xlim([start_Time, end_Time])
        ax4.grid(b = None, which='major', axis='both')
        applySmartTimeTicks(ax4, [start_Time, end_Time])  
        
        
        
        #### Here we have a 5th plot but its going to be on the 4th subplot so were good fam
        ax5 = ax4.twinx()
          
        ax5.plot(CIS_time, maskpro, color = 'r', lw=1)
        ax5.set_ylabel('Protons per $cm^{3}$')          
        ax5.set_xlim([start_Time, end_Time])
                
        ax5.grid(b = None, which='major', axis='both')
    
        
        
        
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
        #### Probably want to give the files a better name than me        
        fig.savefig("/home/doge/Research/TAILWAG/tailwagtools/Output/*" + c + " with interval hours: " +  str(n)+ ".jpg") 
   

   
    
     
    
   

    
    
   
   
print("GAME OVER, THE USER WINS")
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
  
    
    
    
   