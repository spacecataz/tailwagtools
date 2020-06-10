#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Hello researcher, I want to play a game...This program uses things like GLOB and CDF reader from NASA and spacepy
in order to read CDF files and create graphs of around a specified time (that you can pull from the title of the file) a
nd create 3 graphs with various time intervals on each side of the specified time (up to 6 hours on each side).  y
ou will analyze for analyzing the earth's magnetotail.  This will require y=ou to strip CIS, FGM and datetime data
from these CDF files and create a single page of multiple graphs  with 2 position graphs (X vs Y and X vs Z) and graphs
for densities of oxygen and hydrogen and magnetic field mangitudes, some of these might change depending on what supreme
leader dan says.  dont worry it should be easy to change things around.  Once you get that figured out, save
these graphs to a folder and make sure you properly label them so we know which is which..........

"""
'''
THINGS TO DO STILL:
    
    -- Make things look snazzy with Latex
    -- Make the program create a folder per each event and put the 2, 4, and 6 hour interval graphs in that folder
    -- Make sure there are legends on the graphs
    -- Get the Tsyganenko Bx values on the Bx magnitude graph (the one in the middle)
    -- FIX THOSE FUCKING ARROWS JAMIE!!!!!



'''


####HERES WHERE WE DO ALL OF OUR IMPORTS, SOME OF THEM DONT WORK ON THIS COMPUTER FOR SOME REASON???
############fIGURE THIS STUFF OUT FOR YOUR OWN MACHINE
import glob
import os
 
 

##import spacepy.irbempy 

from spacepy import pybats
import spacepy
from spacepy.pycdf import CDF
import spacepy.omni as om
import spacepy.coordinates as spc 
##import spacepy.irbempy as ib 
from spacepy.pycdf import CDF
import matplotlib.pyplot as pltapplySmartTimeTicks
from spacepy.plot import applySmartTimeTicks, style
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import numpy.ma as ma
from spacepy.plot import style
import scipy
import scipy.signal



import datetime


##style()


####    d8888b. d8888b.  .d88b.   d888b  d8888b.  .d8b.  .88b  d88.      .d8888. d888888b  .d8b.  d8888b. d888888b .d8888.      db   db d88888b d8888b. d88888b
####    88  `8D 88  `8D .8P  Y8. 88' Y8b 88  `8D d8' `8b 88'YbdP`88      88'  YP `~~88~~' d8' `8b 88  `8D `~~88~~' 88'  YP      88   88 88'     88  `8D 88'
####    88oodD' 88oobY' 88    88 88      88oobY' 88ooo88 88  88  88      `8bo.      88    88ooo88 88oobY'    88    `8bo.        88ooo88 88ooooo 88oobY' 88ooooo
####    88~~~   88`8b   88    88 88  ooo 88`8b   88~~~88 88  88  88        `Y8b.    88    88~~~88 88`8b      88      `Y8b.      88~~~88 88~~~~~ 88`8b   88~~~~~
####    88      88 `88. `8b  d8' 88. ~8~ 88 `88. 88   88 88  88  88      db   8D    88    88   88 88 `88.    88    db   8D      88   88 88.     88 `88. 88.
####    88      88   YD  `Y88P'   Y888P  88   YD YP   YP YP  YP  YP      `8888Y'    YP    YP   YP 88   YD    YP    `8888Y'      YP   YP Y88888P 88   YD Y88888P




####HERE, WE USE GLOB IN ORDER TO PARSE THROUGH A FOLDER OF CDF FILES, WILL NEED TO CHANGE THE PATH
############### FOR YOUr OWN MACHINE
######### MAKE SURE THAT THE TIME INTERVAL HAS 6 HOURS ON EACH SIDE OF IT 
for path in glob.glob("/home/doge/python/tailwagtools/Monkey_Test_Folder/*"):
    
  
    ####THESE TWO STATEMENTS GRAB THE TWO FILES YOU WILL NEED AS OBJECTS
    CISdata = CDF(path + '/CIS_FILE.cdf')
    FGMdata = CDF(path + '/FGM_FILE.cdf')
    
  

    ####WHEREAS THESE TWO GET THE TIME TAG ARRAYS YOU WILL NEED A LITTLE LATER...#######
    CIS_time = CISdata['Epoch__C1_PP_CIS']
    FGM_time = FGMdata['time_tags__C1_CP_FGM_SPIN']
  
  
    ####HERE WERE GETTING OUR X,Y AND Z COORDINATE POSITIONS############
    #########WE'RE DIVIDING IN ORDER TO CONVERT TO Re COORDINATES
    x = FGMdata['sc_pos_xyz_gse__C1_CP_FGM_SPIN'][0:,0]/6378.16
    y = FGMdata['sc_pos_xyz_gse__C1_CP_FGM_SPIN'][0:,1]/6378.16
    z = FGMdata['sc_pos_xyz_gse__C1_CP_FGM_SPIN'][0:,2]/6378.16

    

    ####HERE WERE OBVIOUSLY GRABBING Bx, Vx AND Vz .........########
    Bx = FGMdata['B_vec_xyz_gse__C1_CP_FGM_SPIN'][0:,0]   
    Vx = CISdata['V_p_xyz_gse__C1_PP_CIS'][0:,0]
    Vz = CISdata['V_p_xyz_gse__C1_PP_CIS'][0:,2]
   
    #### AND HERE WERE GETTING THE ARRAYS SETUP FOR OUR DENSITIES########
    ####  YOU DONT NECESSARILY NEED THIS (COULD BE MORE EFFICIENT) BUT ITS NICE FOR DEBUGGING, TRUST ME.....
    protons = scipy.signal.medfilt(CISdata['N_p__C1_PP_CIS'], kernel_size = 11)
    oxygens = scipy.signal.medfilt(CISdata['N_O1__C1_PP_CIS'], kernel_size = 11)
    
    ####HERE WE APPLY A MASK TO THE ARRAYS TO FILTER OUT THE CRAZY VALUES THAT 
    ####WILL FUCK THE GRAPHS UP!!!!!!!!!!!
    maskVx = ma.masked_values(Vx, -1.00000E+31)
    maskpro = ma.masked_values(protons, -1.00000E+31)
    maskoxy = ma.masked_values(oxygens, -1.00000E+31)
    
    
    ####GETTING THE FOLDER NAME AS A STRING AND CONVERTING IT INTO A 
    ####DATETIME OBJECT,/MUST INCLUDE THE STRING FOR COMPARISONS
    ####TITLE WILL BE USED LATER IN THE CODE BUT LETS KEEP IT HERE FOR NOW
    Title = os.path.basename(path)  
    Filename = os.path.basename(path + ".00.0")      
    Filedate = datetime.datetime.strptime(Filename, '%Y-%m-%d-%H-%M.%S.%f')
    
    #### Print statement to help figure out wtf is coming out..............
    print("We're now calculating and creating the graphs for the following date:    " + Title)
    
    
    ######PRINT STATEMENTS FOR TESTING OBVIOUSLY........########### 
    '''
    print ("The path is: " + path)
    print("The Filename is: " + Filename)
    print("AND THE FILEDATE IS: " )
    print(Filedate)
    
    for thing in FGM_time:
        if thing <= Filedate:
            print ("LLLLLLLAAAAVVVAAAAAA")
        else:
            print("NYOOOOOOOOOOOOOOOOOOOOOOO")
    '''
    #### ABOVE, THE STUFF IN THE FGM TIME CAN BE PROPERLY COMPARED TO THE 
    #### FILEDATE, AND WERE GETTING LAVA AND NYOO
    
    
    
   
    
    
    
    
    ####LOOP THAT WILL ITERATE 3 TIMES TO GIVE US OUR 3 GRAPHS PER EPOCH........#####
    for n in range(2, 7, 2):
    
        ####HERE WE CREATE THE START AND END TIMES AS DATETIME OBJECTS, THEY WILL DECREASE AND INCREASE WITH THE LOOP BY N HOURS#######
        start_Time = Filedate - datetime.timedelta(hours = n)
        end_Time   = Filedate + datetime.timedelta(hours = n)
        
        
        
        #### PRINT STATEMENTS TO TEST WHAT THE START AND STOP TIMES ARE, 
        ####  AND TO ALSO TEST WHAT THEIR TYPES ARE.......  ##spacepy.add_planet(ax1)
        ##spacepy.add_planet(ax2)
        #### ALL THREE HAVE THE SAME TYPE AT THE MOMENT......
        '''+++
        print (start_Time)
        print (end_Time)
        
        print type(start_Time)
        print type(end_Time)
        print type(Filedate)
        '''
        
        
        #### ANOTHER TEST TO MAKE SURE THE START, END, AND FGM ARE ALL THE 
        #### SAME TYPES.....AND YES THEY ARE!!!!!!
        '''
        for line in FGM_time:
            if (type(line) == type(start_Time)) & (type(line) == type(end_Time)):
                print ("They are the same types dawg")
            else:
                print("These are not the same types,python median filter or the droids you're looking for......")
        '''
        '''
        print("Here the the first few thigns in x position...")
        for line in x[0:10]:
            print line
        '''
        
        ####  HERE ARE SOME ARRAYS THAT YOU NEED TO CREATE FOR GRAPHS LATER.....
        orb_x = []
        orb_y = []
        orb_z = []
        orb_time = []
        
        '''
        for stamp in FGM_time:
            if (stamp >= start_Time) & (stamp <= end_Time):
                print FGM_time[0]
                print ("NYANCAT")
        '''      
                
         #####  HERE WERE ADDING STUFF TO THOSE ARRAYS FROM ABOVE AND USING THE NEW INTERVALS WHILE MAKING SURE THAT E
        # ####### VERYTHING PROPERLY LINES UP "TIMEWISE"
        for i in range(len(FGM_time)):
            
            if (FGM_time[i] >= start_Time) & (FGM_time[i] <= end_Time):
                orb_x.append(x[i])
                orb_y.append(y[i])
                orb_z.append(z[i])
                orb_time.append(FGM_time[i])

            
       
   
        
        
        
        
        
        '''    
        print("The lengths of x, y, z and FGM time are:   ")
        print len(x)
        print len(y)
        print len(z)
        print len(FGM_time)
        
        print("The lengths of orbx, orby, orbz and orbtime are:   ")
        print len(orb_x)
        print len(orb_y)
        print len(orb_z)
        print len(orb_time)
        '''        
        
       ###### 
        half_len = len(orb_time)/2
        print (half_len)
        
        fig = plt.figure(figsize=(8.5,11))
        fig.suptitle( Title + " with interval hours: " +  str(n), fontsize=16)
        
        
        ####Creating the 4 subplots were going to need
        ax1, ax2 = fig.add_subplot(321), fig.add_subplot(322)
        ax3, ax4 = fig.add_subplot(312), fig.add_subplot(313)
        
        
        
        
        
    
       
        ax1.plot(orb_x, orb_y)
        ax1.set(xlabel='X in R$_E$', ylabel='Y in R$_E$',
                title=' X vs Y position')
        
        ax1.arrow(orb_x[half_len], orb_y[half_len], orb_x[half_len + 1] - orb_x[half_len], orb_y[half_len + 1] - orb_y[half_len], shape='full', lw=.1, length_includes_head=True, head_width=.1)
        ##axa = fig.add_subplot(111)
        ##add_planet(axa)
        
        
        
        
        
        ax2.plot(orb_x, orb_z)
        ax2.set(xlabel='X in R$_E$', ylabel='Z in R$_E$',
                title=' X vs Z position')
        ax2.arrow(orb_x[half_len], orb_z[half_len], orb_x[half_len + 1] - orb_x[half_len], orb_z[half_len + 1] - orb_z[half_len], shape='full', lw=.1, length_includes_head=True, head_width=.1)
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
        ##spacepy.add_planet(ax2)
        
        
        
        #### Here we have a 5th plot but its going to be on the 4th subplot so were good fam
        ax5 = ax4.twinx()
          
        ax5.plot(CIS_time, maskpro, color = 'r', lw=1)
        ax5.set_ylabel('Protons per $cm^{3}$')          
        ax5.set_xlim([start_Time, end_Time])
                
        ax5.grid(b = None, which='major', axis='both')
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
        #### Probably want to give the files a better name than me        
        fig.savefig("/home/doge/python/tailwagtools/Final_Graph_Folder/*""The graph for " +  Title + " with interval hours: " +  str(n) + ".jpg")
        

    
####  STATEMENT TO KNOW THAT WEVE REACHED THE END OF THE PROGRAM
print ("GAME OVER, THE USER WINS!!")



    
      
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
        
        