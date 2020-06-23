# -*- coding: utf-8 -*-



##### HERE WERE DOING IMPORTS

import datetime as dt
import glob
import matplotlib
import numpy as np
import numpy.ma as ma
import os
import scipy
import scipy.signal
import spacepy
import pandas as pd
from datetime import datetime
from datetime import timedelta
from dateutil.relativedelta import relativedelta




os.environ["CDF_LIB"] = "/home/doge/Research/TAILWAG/tailwagtools/CDF/lib"
from spacepy.pycdf import CDF

path = '/home/doge/Research/TAILWAG/tailwagtools/DATA/Monthly_Data'


#### You need this line for the strptime function bc the geniuses that wrote 
#### the datetime module named their class datetime 
#### because they're FUCKING MORONS!!!!
##from datetime import datetime


##  SET UP SOME DUMMY VARIABLES FOR CHANGING AND CHECKING LATER!!!!!!
dummy_date = dt.datetime(2001, 1, 1, 1, 1, 1, 0)
dummy_CIS = "CIS9"
dummy_FGM = "FGM9"


#####  Here are some print statements to check the dummy outputs
##print (dummy_date)
##print (dummy_CIS)
##print (dummy_FGM)
##print ("End of dummy test statements \n")



##### READ THE DATA FROM THE EXCEL FILE
Event_Data = pd.read_excel("Event_Points.xlsx", header = 0)


##### HERE WE GRAB THE DATA THAT WERE GOING TO NEED FROM THE EXCEL COLUMNS 
Point_List = Event_Data['Narrowed Point']
CIS_CAPS_List = Event_Data['USED CIS']
FGM_CAPS_List = Event_Data['USED FGM']

#### MAKING A LIST THAT CONVERTS THE USED CLUSTER 
CIS_STRING_List = [clus.lower() + "_pps_cis_" for clus in CIS_CAPS_List]
FGM_STRING_List = [clus.lower() + "_cps_fgm_spin_" for clus in FGM_CAPS_List]

for thing in CIS_STRING_List:
    print(thing)

##[i for i in os.listdir(path) if os.path.isfile(os.path.join(path,i)) and \
         ##'001_MN_DX' in i]

##### Section where we convert out strings to datetime objects
Date_List = [dt.datetime.strptime(date,'%Y-%m-%dT%H:%M:%S.%f') for date in Point_List]
Start_Date = [date - dt.timedelta(hours = 12) for date in Date_List]
End_Date = [date + dt.timedelta(hours = 12) for date in Date_List]


#####  print statesments for checking if things are correct
'''
print("Now were gonna check the start, date, and end date lists:" + "\n")
for (a, b, c) in zip (Start_Date, Date_List, End_Date):
    print (a, b, c)
    
for (a, b) in zip(CIS_List, FGM_List):
    print(a, b)
'''   


'''
for (a, b) in zip (Start_Date, End_Date):
    if ((a.month == b.month) and (a.year == b.year)):
        print("badger")
    else:
        print("mushroom")
'''


#####  BIG FOR LOOP THATS GONNA BE FOR EVERY EVENT SORRY I DONT REMEMBER DICTIONARIES ATM AND IM REALLY FUCKING LAZY
for (point,date, start, end, cis, fgm, CISPICK, FGMPICK) in zip(Point_List, Date_List, Start_Date, End_Date, CIS_STRING_List, FGM_STRING_List, CIS_CAPS_List, FGM_CAPS_List):
    
   
    
    
    #####  SET UP ALL OF THE ARRAYS HERE THAT WE WILL USE TO THROW OUR VALUES FROM THE FILES INTO
    CIS_time = []
    FGM_time = []
    x = []
    y = []
    z = []
    Vx = []
    Bx = []
    protons = []
    oxygens = []   
    
    month_diff = end.month - start.month
    
    
    print("################" + "\n")
    print("We're now doing the following point: " + point + "\n")
    print("The start point is:    " + str(start) + "\n")
    print("The end point is:    " + str(end)+ "\n")
    print("We are now entering the for loop with the month diff as: " + str(month_diff) + "\n")
    
    
    for i in range(month_diff + 1):
        
        print ("The iterative month is : " + str(i) + "\n")
        
        
        ##WERE GONNA USE this variable to iterate through the various months taht were going to have....
        ###  What we need to do here is add the i months to the start date
        diff_date = start + relativedelta(months=i)
        print("Our differential date is:   " + str(diff_date) + "\n")
        
        
        
        
        
        
        diff_str = diff_date.strftime("%Y%m")
        diff_cis = cis + diff_str
        diff_fgm = fgm + diff_str
       
        
        print("The diff string is:   " + str(diff_str))
        print("The cis string is:   " + str(diff_cis))
        print("The fgm string is:   " + str(diff_fgm) + "\n")
        
        
        
        CISFILE = [i for i in os.listdir(path) if os.path.isfile(os.path.join(path,i)) and diff_cis in i]
        FGMFILE = [i for i in os.listdir(path) if os.path.isfile(os.path.join(path,i)) and diff_fgm in i]
        
        print(CISFILE)
        print(FGMFILE)
        
        
        cispath = "/" + str(CISFILE).strip('[]').strip("'")
        fgmpath = "/" + str(FGMFILE).strip('[]').strip("'")
        
        print(cispath)
        print(fgmpath)
        
        print(path + cispath)
        print(path + fgmpath)
        
        
        CISDATA = CDF(path + cispath)
        FGMDATA = CDF(path + fgmpath)
        
        ####  USE THIS TO CHECK WHATS ACTUALLY INSIDE THE CIS FILE
        ##for line in CISDATA:
            ##print (line)
        
        ####  USE THIS TO CHECK WHATS ACTUALLY INSIDE THE FGM FILE
        ##for line in FGMDATA:
            ##print (line)
            
            
            
         #####  ARRAYS FOR THE TIME TAGS   
        ##CIS_time_temp = CISDATA['Epoch__' + CISPICK + '_PP_CIS']
        ##FGM_time_temp = FGMDATA['time_tags__' + FGMPICK + '_CP_FGM_SPIN']
            
        
        
        
        
        ##### ARRAYS TO WITH X, Y AND Z COORDS FROM FGM CONVERTED TO Re
        x_temp = FGMDATA['sc_pos_xyz_gse__' + FGMPICK + '_CP_FGM_SPIN'][0:,0]/6378.16
        y_temp = FGMDATA['sc_pos_xyz_gse__' + FGMPICK + '_CP_FGM_SPIN'][0:,1]/6378.16
        z_temp = FGMDATA['sc_pos_xyz_gse__' + FGMPICK + '_CP_FGM_SPIN'][0:,2]/6378.16
       
       
        ####HERE WERE OBVIOUSLY GRABBING Bx, Vx AND Vz .........########
        Bx_temp = FGMDATA['B_vec_xyz_gse__' + FGMPICK + '_CP_FGM_SPIN'][0:,0]   
        Vx_temp = CISDATA['V_p_xyz_gse__' + CISPICK + '_PP_CIS'][0:,0]
        
        ##protons_temp = CISDATA['N_p__' + CISPICK + '_PP_CIS']
        ##oxygens_temp = CISDATA['N_O1__' + CISPICK + '_PP_CIS']
        
        

        
       
     
        
'''      
        for i in range(len(CIS_time_temp)):
            if (CIS_time_temp[i]>= start) & (CIS_time_temp[i] <= end):
                CIS_time.append(CIS_time_temp[i])
                Vx.append(Vx_temp)
                protons.append(protons_temp[i])
                oxygens.append(oxygens_temp[i])
                
        for i in range(len(FGM_time_temp)):
            if (FGM_time_temp[i]>= start) & (FGM_time_temp[i] <= end):
                FGM_time.append(FGM_time_temp[i])
                x.append(x_temp)
                y.append(y_temp)
                z.append(z_temp)
                Bx.append(Bx_temp)
'''                
           
        
        
        
       
        
        ##protons = scipy.signal.medfilt(CISDATA['N_p__' + CISPICK + '_PP_CIS'], kernel_size = 11)
        ##oxygens = scipy.signal.medfilt(CISDATA['N_O1__' + CISPICK + '_PP_CIS'], kernel_size = 11)
        
       
        
        ##################################################################################################
        
  
    #################################################################################################################   
   
    
    
    #### Here were masking the three arrays to filter out the crazy values that ruin the graphs
    ##maskVx = ma.masked_values(Vx, -1.00000E+31)
    ##maskpro = ma.masked_values(protons, -1.00000E+31)
    ##maskoxy = ma.masked_values(oxygens, -1.00000E+31)
      
    for line in CIS_time:
        print(line)
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ###print("\n" + "We're now calculating and creating the graphs for the following date:    " + point)


    ####LOOP THAT WILL ITERATE TO GIVE US OUR EPOCH GRAPHS#####
    ###or n in range(2, 12, 2):



        '''
        
        ####Creating the 4 subplots were going to need
        ax1, ax2 = fig.add_subplot(321), fig.add_subplot(322)
        ax3, ax4 = fig.add_subplot(312), fig.add_subplot(313)
        '''
        
        
        
        
        '''
        ax1.plot(orb_x, orb_y)
        ax1.set(xlabel='X in R$_E$', ylabel='Y in R$_E$',
                title=' X vs Y position')
        
        ##ax1.arrow(orb_x[half_len], orb_y[half_len], orb_x[half_len + 1] - orb_x[half_len], orb_y[half_len + 1] - orb_y[half_len], shape='full', lw=.1, length_includes_head=True, head_width=.1)
        ##axa = fig.add_subplot(111)
        ##add_planet(axa)
        
        
        
        
        
        ax2.plot(orb_x, orb_z)
        ax2.set(xlabel='X in R$_E$', ylabel='Z in R$_E$',
                title=' X vs Z position')
        ##ax2.arrow(orb_x[half_len], orb_z[half_len], orb_x[half_len + 1] - orb_x[half_len], orb_z[half_len + 1] - orb_z[half_len], shape='full', lw=.1, length_includes_head=True, head_width=.1)
        ##ax2 = fig.add_subplot(322)
        ##add_planet(ax2)
        
        
        
        
        
        ax3.plot(FGM_time, Bx)        
        ax3.set(xlabel='time', ylabel='$B_X$ in nT',
                title='$B_X$ Magnitude')
        ax3.set_xlim([start_Time, end_Time])
        ax3.hlines(0, FGM_time[0], FGM_time[-1], lw=2, color='black', linestyle = '--')
        applySmartTimeTicks(ax3, [start_Time, end_Time])
        
        
        ax4.set_xlabel('time')
        ax4.set_title('Proton and Oxygen Densities')
        
        ax4.plot(CIS_time, maskoxy, 'g-', lw=1, alpha=0.7)
        ax4.set_ylabel('Oxygen per $cm^{3}$')        
        ax4.set_xlim([start_Time, end_Time])
        ax4.grid(b = None, which='major', axis='both')
        applySmartTimeTicks(ax4, [start_Time, end_Time])  ##spacepy.add_planet(ax1)
        ##spacepy.add_planet(ax2)https://www.howtoforge.com/how-to-install-microsoft-teams-linux-on-ubuntu-and-centos/
        
        
        
        #### Here we have a 5th plot but its going to be on the 4th subplot so were good fam
        ax5 = ax4.twinx()
          
        ax5.plot(CIS_time, maskpro, color = 'r', lw=1)
        ax5.set_ylabel('Protons per $cm^{3}$')          
        ax5.set_xlim([start_Time, end_Time])
                
        ax5.grid(b = None, which='major', axis='both')
        applySmartTimeTicks(ax4, [start_Time, end_Time])
        
        
        
        ##### With these 2 blocks were adding planet earth to our postion graphs
        ax1 = fig.add_subplot(321)
        ##ax2 = fig.add_subplot(322)
       
        ##spacepy.pybats.add_planet(ax1)
        ##spacepy.pybats.add_planet(ax2)
        
        
        
        #### Here were making sure everything stays nice and cool looking
        plt.tight_layout(rect=[0, 0, 1, .95])
        plt.show()
        
        
        ####  HERE WERE SAVING THE GRAPHS TO A FOLDER AND GIVING THE FILE A TITLE THAT HAVE THE DATE
        #####  AND THE INTERVAL HOURS USED
        #### Probably want to give the files a better name than me        
        fig.savefig("/home/doge/Research/TAILWAG/tailwagtools/Output_Folder/*""The graph for " +  Title + " with interval hours_ " +  str(n) + ".jpg")  
        '''
        
        
        
print("GAME OVER, THE USER WINS")
        
   
        
        
        
        
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
  
    
    
    
   