'''
    GREETINGS USER!!!! this is the script that you can use to generate fancy 
graphs to determine satellite positions, Bx values for cluster AND
the Tsyganenko models, AND GET proton and Oxygen outflow....

    In order to use this script, you must use the magic excel file called 
    Event_Points...it will read 3 columns and kick out some graphs using 
    functions from the tailwag.py.....its pretty simply, check a look!!!
    
    Problems that need to be fixed still:
        
        - There's a lot of bs related to reading average outflow levels 
            from 2005 to later...this is necessarily this programs fault
            but keep an eye out for that K?


'''
##### Here we're doing our imports
import datetime as dt
from datetime import datetime
import pandas as pd
import tailwag as tw


    
##START HERE!!!    
def main():
    
    ##  Just some fancy statements to start this thing off.... 
    ##print(spacepy.plot.available(returnvals = False))
    print()
    print("Hello there, I wanna play a game")
    print()

    ####  Reading the excel file and putting it in memeory, need the columns...
    Event_Data = pd.read_excel("Event_Points.xlsx", header = 0)
    
    ####  getting three arrays: Point_list as strings, Sat_List will give the 
    ####    cluster satellite used as an int, Op_List has the string "SKIP"
    ####     to determine if we will skip the date, Date_List will be an array
    ####       of datetime objects that are the same/made from Point_List.....
    Point_List = Event_Data['Narrowed Point']   
    Sat_List   = Event_Data['USED SAT']
    Op_List    = Event_Data['Operation']
    Date_List  = [dt.datetime.strptime(date,'%Y-%m-%d %H:%M:%S.%f') for date in Point_List]
    
    #### Loop to determine whether we skip the date, or make a graph with it....
    ####  Go read the function to see what all the cool inputs and 
    ####    conditions do...HINT: THIS IS CAN BE USED TO CHECK OUTFLOW LEVELS!!
    for (a, b, c) in zip(Date_List, Sat_List, Op_List):
        if c == "SKIP":
            continue
        else:
            tw.fusion(a, b, 12, add_tsyg=True, add_scatter = True, outdir = 'check_plots/')
            
    #### We reached the end, so heres a print statement to let us know....
    print("game over, the user wins")
    
#########################################################################################################
# If this file is used via IPython's "run" magic command,
# Evaluate on all figures.
if __name__ == '__main__':
    main()
       
    
    
     
    
  
    
    

    
 
