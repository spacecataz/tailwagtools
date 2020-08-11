#import things
import matplotlib.pyplot as plt; plt.ion() #interactive matplotlib
import pandas as pd
import numpy as np
from tailwagtools import tailwag
import datetime as dt
import scipy
from scipy import stats
import math

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
    >>> from tailwagtools import scatter
    >>> events['epoch'] = scatter.read_Event_Points("tailwagtools/Event_Points.xlsx")
    >>> print(events['epoch'])
    {'epoch' : [datetime.datetime(...),...]}
    
    Returns array of datetime objects stored in events['epoch']
    
    '''

    #Read filename as an Excel file using Pandas
    event_data = pd.read_excel(filename, header = 0)
    #Example: event_data = pd.read_excel("tailwagtools/Event_Points.xlsx", header = 0)
    
    #Separate points list
    point_list = event_data['Narrowed Point']
    
    #Convert point to datetime epoch
    epoch  = [dt.datetime.strptime(date,'%Y-%m-%d %H:%M:%S.%f') for date in point_list]
    
    return epoch


def get_crossing_info(epoch):
    '''
    For a given crossing epoch, calculate and return:
    - the crossing time for all Tsgyanenko models under consideration
    - the average density before and after the cluster crossing
      for both H+ and O+

    Parameters
    ----------
    epoch : datetime.datetime
       The time of the cluster crossing of the plasmasheet.

    Returns
    -------
    t_cluster : datetime.datetime
        Time of crossing observed in Cluster
    t_times : datetime.datetime
        Time of crossing observed in T89, T96, AND T01STORM


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
    
    
    #Calculate Cluster average density after
    loc = t_cis>t_cis[index_cluster] #location after Cluster crossing
    cluster_after = h_dens[loc].mean() + o_dens[loc].mean() #average H+ density + average O+ density
    
    #Calculate Cluster average density before
    loc = t_cis<t_cis[index_cluster] #location before Cluster crossing
    cluster_before = h_dens[loc].mean() + o_dens[loc].mean() #average H+ density + average O+ density
    
    
    #Get crossing time and densities for each Tsyg model
    t_times = {} #container for times
    index_Tsyg = {}
    Tsyg_after = {} #container for densities after crossing
    Tsyg_before = {} #container for densities before crossing

    for vers in ['T89','T96', 'T01STORM']:
        #Obtain Tsyg crossing times
        t, b_Tsyg = tailwag.gen_sat_tsyg(data, extMag = vers) #get Tsyg data
        loc = np.abs(b_Tsyg[:,0])==np.abs(b_Tsyg[:,0]).min() #location of crossing
        t_times[vers] = t[loc] #Tsyg crosing time
        
        #Obtain average plasma densities before/after for each Tsyg crossing time
        index_Tsyg[vers] = np.where(t_cis<=t_times[vers])[0][-1] #locate CIS time <= Tsyg crossing time
        tdelt_before = t_times[vers]-t_cis[index_Tsyg[vers]] #timedelta using t_cis less than t_cluster
        tdelt_after = t_cis[index_Tsyg[vers]+1]-t_times[vers] #timedelta using t_cis more than t_cluster
        if tdelt_before > tdelt_after: #compare the two, if the timedelta of t_cis less > t_cis more, then change index_cluster to the location of t_cis more
            index_Tsyg[vers] = index_Tsyg[vers]+1
        
        #THIS PART DOES NOT WORK ):
        #Calculate Tsyg average density after
        loc = t_cis>t_cis[index_Tsyg[vers]] #location after Cluster crossing
        Tsyg_after[vers] = h_dens[loc].mean() + o_dens[loc].mean() #average H+ density + average O+ density
        
        #Calculate Tsyg average density before
        loc = t_cis<t_cis[index_Tsyg] #location before Cluster crossing
        Tsyg_before[vers] = h_dens[loc].mean() + o_dens[loc].mean() #average H+ density + average O+ density

    return(t_cluster, t_times, cluster_after, cluster_before, Tsyg_before, Tsyg_after)


