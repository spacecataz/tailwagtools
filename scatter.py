#import things
import matplotlib.pyplot as plt; plt.ion() #interactive matplotlib
import pandas as pd
import numpy as np
import tailwag
import datetime as dt
import scipy
from scipy import stats
import math


    
def to_pickle(data_dic, picklename):
   
   import pickle

   pickle_out = open(picklename, "wb")
   pickle.dump(data_dic, pickle_out)
   pickle_out.close()




def read_Event_Points(filename):
    '''
    For the Event_Points data file, read the Excel file via Pandas and convert the list of points into an epoch.
    
    Parameters
    ----------
    filename : string
       Name of the Excel file storing the event points.
       
    Returns
    -------
    events['epoch'] : datetime.datetime
        The time of the cluster crossing of the plasmasheet.
        
    Example
    _______
    >>> # From inside the tailwagtools directory...
    >>> import scatter
    >>> epochs = scatter.read_Event_Points("Event_Points.xlsx")
    >>> print(epochs)
    [datetime.datetime(2001, 8, 15, 9, 9, 38), ... ]
    
    '''

    #Read filename as an Excel file using Pandas
    event_data = pd.read_excel(filename, header = 0)
    #Example: event_data = pd.read_excel("tailwagtools/Event_Points.xlsx", header = 0)
    
    #Read the points as strings as a string to decide what to do with the point
    point_list = event_data['Narrowed Point']
    op_List    = event_data['Operation']
    
    
    ###create the epoch array that were gonna need for the datetimes
    epoch = []
    
    ##loop through the dates in point_list and decide whether to add them
    ##or not based on the op_list variable...needs to be "SKIP".....
    for (date, op) in zip(point_list, op_List):
        if op == "SKIP":
            continue
        else:
            epoch.append(dt.datetime.strptime(date,'%Y-%m-%d %H:%M:%S.%f'))
    
    #Convert point to datetime epoch
    ##epoch  = [dt.datetime.strptime(date,'%Y-%m-%d %H:%M:%S.%f') for date in point_list]
    
    return epoch


def get_crossing_info(epoch, debug=False):
    '''
    For a given crossing epoch, calculate and return:
    - the crossing time for all Tsgyanenko models under consideration
    - the average density before and after the cluster crossing
      for both H+ and O+

    Parameters
    ----------
    epoch : datetime.datetime
       The time of the cluster crossing of the plasmasheet.

    Other Parameters
    ----------------
    debug : boolean
       Turn on debug information.  Defaults to False

    Returns
    -------
    t_cluster : datetime.datetime
        Time of crossing observed in Cluster
    t_times : array of datetime.datetime
        Time of crossing observed in T89, T96, AND T01STORM
        T89 =  t_times[0], T96 = t_times[1], AND T01STORM = t_times[2]      
    cluster_before : array of float values
        An array of values that yield the average value before 
        the time crossing.. ...
        cluster_before[0] = H value, cluster_before[1] = O value
    cluster_after : array of float values
        An array of values that yield the average value after 
        the time crossing.. ...
        cluster_after[0] = H value, cluster_after[1] = O value


    Example
    _______
    >>> import datetime as dt
    >>> import scatter
    >>> epoch = dt.datetime(2001,8,19,20,0,0)
    >>> t_times = scatter.get_crossing_info(epoch)
    >>> print(t_times)
    {'T89': dmarray([datetime.datetime(...)], dtype=object),
        'T96': dmarray([datetime.datetime(...)], dtype=object),
        'T01STORM': dmarray([datetime.datetime(...)], dtype=object)
        }
    '''
    
    #Fetch cluster data given an epoch
    data = tailwag.fetch_cluster_data(epoch) #get Cluster data
    loc = np.abs(data['b'][:,0])==np.abs(data['b'][:,0]).min() #location of Cluster crossing
    t_cluster = data['fgm_time'][loc] #Cluster crossing time
    
    #Fetch cluster plasma data
    t_cis = data['cis_time'] #get Cluster CIS times
    h_dens = data['dens_h']
    o_dens = data['dens_o']
    
    #Location t_cis closest to t_cluster
    #MXB notes: should we change this into tricky indexing instead of np.where? this is probably too chonky
    index_cluster = np.where(t_cis<=t_cluster)[0][-1] #locate CIS time <= cluster crossing time
    tdelt_before = t_cluster-t_cis[index_cluster] #timedelta using t_cis less than t_cluster
    tdelt_after = t_cis[index_cluster+1]-t_cluster #timedelta using t_cis more than t_cluster
    if tdelt_before > tdelt_after: #compare the two, if the timedelta of t_cis less > t_cis more, then change index_cluster to the location of t_cis more
        index_cluster = index_cluster+1

    # Print some debug information:
    if debug:
        print('----------DEBUG----------')
        print(f'\tCrossing epoch from file: {epoch}')
        print(f'\tCrossing time from CIS epochs: {t_cis[index_cluster]}')
        
     
    
        #Calculate Cluster average density after
    loc = t_cis>t_cis[index_cluster] #location after Cluster crossing
    cluster_after = h_dens[loc].mean() + o_dens[loc].mean() #average H+ density + average O+ density
    
        #Calculate Cluster average density before
    loc = t_cis<t_cis[index_cluster] #location before Cluster crossing
    cluster_before = h_dens[loc].mean() + o_dens[loc].mean() #average H+ density + average O+ density
    

        #Get crossing time and densities for each Tsyg model
    t_times = {} #container for times
    index_Tsyg = {}

    for vers in ['T89','T96', 'T01STORM']:
            #Obtain Tsyg crossing times
        t, b_Tsyg = tailwag.gen_sat_tsyg(data, extMag = vers) #get Tsyg data
        loc = np.abs(b_Tsyg[:,0])==np.abs(b_Tsyg[:,0]).min() #location of crossing
            # It's possible to have no results from a given Tsyg model.
            # If that's the case, return NaNs.
        if t[loc].size > 0:
            t_times[vers] = t[loc][0] #Tsyg crosing time
        else:
                t_times[vers] = np.nan
    
    return(t_cluster, t_times, cluster_after, cluster_before)


def split_crossing_info(epoch, sat, interval_hour, debug=False):
    '''
    For a given crossing epoch, calculate and return:
    - the crossing time for all Tsgyanenko models under consideration
    - the average density before and after the cluster crossing
      for both H+ and O+

    Parameters
    ----------
    epoch : datetime.datetime
       The time of the cluster crossing of the plasmasheet.
    sat : an int
        Cluster satellite 1
    interval_hour : int
        A int that will give us the number of hours on each side of our POI
    
    Other Parameters
    ----------------
    debug : boolean
       Turn on debug information.  Defaults to False

    Returns
    -------
    t_cluster : datetime.datetime
        Time of crossing observed in Cluster
    t_times : array of datetime.datetime
        Time of crossing observed in T89, T96, AND T01STORM
        T89 =  t_times[0], T96 = t_times[1], AND T01STORM = t_times[2]      
    cluster_before : array of float values
        An array of values that yield the average value before 
        the time crossing.. ...
        cluster_before[0] = H value, cluster_before[1] = O value
    cluster_after : array of float values
        An array of values that yield the average value after 
        the time crossing.. ...
        cluster_after[0] = H value, cluster_after[1] = O value


    Example
    _______
    >>> import datetime as dt
    >>> import scatter
    >>> epoch = dt.datetime(2001,8,19,20,0,0)
    >>> t_times = scatter.get_crossing_info(epoch)
    >>> print(t_times)
    {'T89': dmarray([datetime.datetime(...)], dtype=object),
        'T96': dmarray([datetime.datetime(...)], dtype=object),
        'T01STORM': dmarray([datetime.datetime(...)], dtype=object)
        }
    '''
    
    #Fetch cluster data given an epoch
    data = tailwag.fetch_cluster_data(epoch, sat=sat, tspan=interval_hour) #get Cluster data
    
    #Fetch cluster plasma data
    t_cis = data['cis_time'] #get Cluster CIS times
    h_dens = data['dens_h']
    o_dens = data['dens_o']
    
    #Location t_cis closest to t_cluster
    # MXB notes: should we change this into tricky indexing instead of np.where?
    # this is probably too chonky
    # DTW NOTES: I think this should be epoch.
    index_cluster = np.where(t_cis<=epoch)[0][-1] #locate CIS time <= cluster crossing time
    
    ##print("The tdeltbefore is the following:     " + str(tdelt_before))
    ##print("The tcis index cluster + 1 is the following:     " + str(t_cis[index_cluster+1]))
    ##print("The tcluster is the following:     " + str(t_cluster))
   
    # Print some debug information:
    if debug:
        print('----------DEBUG----------')
        print(f'\tCrossing epoch from file:      {epoch}')
        print(f'\tCrossing time from CIS epochs: {t_cis[index_cluster]}')
     
    cluster_before = []
    cluster_after  = []
        #Calculate Cluster average density after
    loc = t_cis>t_cis[index_cluster] #location after Cluster crossing
    cluster_after.append(h_dens[loc].mean()) #average H+ density + average O+ density
    cluster_after.append(o_dens[loc].mean())
        #Calculate Cluster average density before
    loc = t_cis<t_cis[index_cluster] #location before Cluster crossing
    cluster_before.append(h_dens[loc].mean()) #average H+ density + average O+ density
    cluster_before.append(o_dens[loc].mean())
    

        #Get crossing time and densities for each Tsyg model
    t_times = {} #container for times
    index_Tsyg = {}

    for vers in ['T89','T96', 'T01STORM']:
            #Obtain Tsyg crossing times
        t, b_Tsyg = tailwag.gen_sat_tsyg(data, extMag = vers) #get Tsyg data
        loc = np.abs(b_Tsyg[:,0])==np.abs(b_Tsyg[:,0]).min() #location of crossing
            # It's possible to have no results from a given Tsyg model.
            # If that's the case, return NaNs.
        if t[loc].size > 0:
            t_times[vers] = t[loc][0] #Tsyg crosing time
        else:
                t_times[vers] = np.nan
    

    return(t_cluster, t_times, cluster_after, cluster_before)


if __name__ == '__main__':
    '''
    If run as script, the following commands will execute.
    '''

    # Let's set up a quick demo for now:
    # Goal: have entire plotting/statistical analysis run here.

    # Get crossing times from Excel file.
    all_epochs = read_Event_Points("Event_Points.xlsx")

    # Create a container for the data.  GOAL: have vectors to easily plot!
    # We can keep adding entries here for things like Vz, F10.7, etc. etc.
    data = {'epoch':[], 'dtT89':[], 'dDensH':[], 'dDens':[], 'dDensO':[],
            'dtT96':[], 'dtT01STORM':[]}

    # Now, loop through all events and gather information:
    for i, t in enumerate(all_epochs[:60]):

        print(f'****WORKING ON CROSSING #{i} AT {t}*****')
        # Get crossing info from cool function:
        t_clus, t_times, c_after, c_before = get_crossing_info(t)
        
        # Store info into our big data container:
        data['epoch'].append(t)
        for vers in ['T89','T96', 'T01STORM']:
            if type( t_times[vers] ) != type(t):
                data['dt'+vers].append(np.nan)
            else:
                data['dt'+vers].append( (t - t_times[vers]).total_seconds()/60. ) #in minutes!
        data['dDens'].append( c_before - c_after )
    
    to_pickle(data,"tastypickle")


    # If we're smart, we're saving this data to an external file.
    # We then make the plot and analysis stuff SEPARATELY because
    # processing the data takes waaayyy too long.
    # Save the data as a pickle.
    # see here: https://nbviewer.jupyter.org/url/www-personal.umich.edu/~dwelling/python/notebooks/primer03_fileIO.ipynb

    
    # Add verification plotting here!
