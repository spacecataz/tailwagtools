.#!/usr/bin/env python
'''
A module for performing analysis of outflow wagging the tail.
'''

# Set install directory:
install_dir = '/'.join(__loader__.path.split('/')[:-1])+'/'

# A map of the header for a CSV info file:
table_map = {'start':4,'end':5,'epoch':3,'kp':7,'f107':9,'cis':10,'fgm':11}

def _parse_table_line(parts):
    '''
    Take one line from table, parse into useful objects, return a dictionary
    of values.  Uses the "table_map" to find correct columns and map to
    dictionary keys.
    '''

    import re
    from dateutil.parser import parse

    raise DepreciationWarning('This function is no longer useful.')
    
    out = {}
    
    # Get times:
    for x in ('start', 'end', 'epoch'):
        out[x] = parse(parts[ table_map[x] ])

    # Get Kp:
    out['kp'] = float(re.search('\d\.?\d{0,3}', parts[table_map['kp']]).group())

    # Get F107
    if parts[ table_map['f107'] ]:
        out['f107'] = float(parts[ table_map['f107'] ])
    else:
        out['f107'] = 0.0
    
    # Save file names.
    for x in ('cis', 'fgm'):
         out[x] = table_map[x]
    
    return out

def parse_event_table(filename='default'):
    '''
    Load the table of Cluster tail passes for use in other scripts.
    Default behavior is to read the CSV file provided with this package.
    This can be overridden via the **filename** kwarg.
    
    The returned object is a list, one entry per event.  Within each entry
    is a dictionary of useful information with the following entries:

    *start* - datetime of the interval start date/time.
    *end*   - datetime of the interval end date/time.
    *epoch* - datetime of the central epoch of the event.
    *kp*    - Kp value at epoch.
    *f107*  - f107 value at epoch.
    *cis*   - string containing name of Cluster CIS file.
    *fgm*   - string containing name of Cluster FGM file.
    '''

    raise DepreciationWarning('This function is no longer useful.')

    # Set file name path:
    if filename == 'default':
        filename=install_dir+'event_info.csv'

    # Open file; slurp contents into list.
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    # Get header:
    head = lines.pop(0).split(',')

    # Create data structure (list of events)
    data = []

    # Loop through lines, skipping year-only rows:
    for l in lines:
        parts = l.split(',')
        if parts[0]: continue # If only year, skip.
        data.append(_parse_table_line(parts))

    return data
    
def gen_sat_tsyg(cluster_data, extMag='T89', dbase='qd1min'):
    '''
    Given a cluster data dict that contains both time and position of a 
    satellite, create magnetic field along that orbit via one of the 
    Tsyganenko empirical models.

    The time range, coordinate system, and other details are pulled 
    automatically from the input data dictionary.

    Parameters
    ==========
    cluster_file : dictionary
        A dictionary of cluster values from fetch_cluster_data.

    Other Parameters
    ================
    extMag : string
        Name of external magnetic field model.  Choices can be found in 
        documentation for 
        [Spacepy's ONERA interface](https://pythonhosted.org/SpacePy/autosummary/spacepy.irbempy.get_Bfield.html)
        Default is 'T89', or the Tsyganenko '89 model.
    dbase : string
        Set the OMNI/QD database source for obtaining solar wind
        and other empirical input values.  Defaults to 'QDHourly'.

    Returns
    =======
    time : array
        An array of datetime objects.
    mag_gsm : array
        A nTime-by-3 array of the magnetic field data.

    OPEN QUESTIONS/TO-DO

    -verification of results- compare to CCMC test?
    
    '''

    import numpy as np
    
    from spacepy.pycdf import CDF
    import spacepy.omni as om
    import spacepy.time as spt 
    import spacepy.coordinates as spc 
    import spacepy.irbempy as ib 

    # Convert CDF values to spacepy-specific classes.
    # These enable coordinate transforms and are required for
    # the ONERA interface.
    time = spt.Ticktock(cluster_data['fgm_time'], 'ISO')
    pos  = spc.Coords( cluster_data['xyz'], cluster_data['coords'], 'car')

    # Get input solar wind and QinDenton values:
    QD = om.get_omni(time, dbase=dbase)

    # get_Bfield expects certain labels for all variables...
    QD['dens'] = QD['Den_P']
    QD['velo'] = QD['Vsw']

    # NOTE: Irbem has a hard limit of 100,000 points.  We frequently
    # exceed that.  We have to build new data containers and loop over
    # points to prevent problems.
    npts, nTimeMax = len(time), 100000
    b_tsyg = np.zeros((npts,3))
    
    # Get b-field along orbit, repeating as required to overcome
    # limits of Irbemlib.
    for i in range(0,npts,nTimeMax):
        iStop = min(npts, i+nTimeMax)
        b_now = ib.get_Bfield(time[i:iStop], pos[i:iStop],
                              extMag=extMag, omnivals=QD)
        b_tsyg[i:iStop,:] = b_now['Bvec']
        
    # Convert magnetic field from GEO to GSM.
    if cluster_data['coords']!='GEO':
        b_geo = spc.Coords( b_tsyg, 'GEO', 'car', ticks=time)
        b_out = b_geo.convert(cluster_data['coords'], 'car')

    # Return time and b-field data.  Strip away extraneous information
    # added by spacepy (via .UTC and .data, which give values only.)
    return(time.UTC, b_out.data)


def get_cluster_filename(epoch, sat=1):
    '''
    Given an epoch, return the CIS and FGM Cluster files that contain the
    data associated with that event.

    Defaults to Cluster 1 spacecraft, use *sat* kwarg to change.

    
    Parameters
    ----------
    epoch : datetime.datetime
         A datetime object corresponding to the time of interest.

    Other Parameters
    ----------------
    sat : int
        The number, 1-4, indicating which satellite to use.

    Returns
    -------
    A two element tuple containing the FGM and CIS filenames in that order.

    '''

    from glob import glob

    # Build path to folder with data:
    path = install_dir + '/DATA/Monthly_Data/'

    # Search for FGM file; stop kindly if not found:
    try:
        fgm=glob(path+f'c{sat}_cps_fgm_spin_{epoch:%Y%m}*_{epoch:%Y%m}*.cdf')[0]
    except IndexError:
        raise FileNotFoundError('No matching FGM file for this time.')

    # Search for CIS file; stop kindly if not found:
    try:
        cis=glob(path+f'c{sat}_pps_cis_{epoch:%Y%m}*_{epoch:%Y%m}*.cdf')[0]
    except IndexError:
        raise FileNotFoundError('No matching FGM file for this time.')

    return (fgm, cis)

    
def fetch_cluster_data(epoch, sat=1, tspan=12, kernel_size=7,
                       coords='GSM', debug=False, prune=True):
    '''
    For a given *epoch* and time span, *tspan*, about *epoch*, open and
    re-organize data from Cluster CIS and FGM files into a ready-to-use 
    dictionary.

    The key-value pairs inside the returned dictionary are described below.
    Note that there are two time vectors as the CIS and FGM data report
    at different cadences.

    | Key      | Value                                             |
    |----------|---------------------------------------------------|
    |fgm_time  | Array of datetimes for the FGM data.              |
    |cis_time  | Array of datetimes for the CIS data.              |
    |dens_h    | Array of density values, protons only.            |
    |dens_o    | Array of density values, oxygen ions only.        |
    |v_h       | 3D velocity of cold protons DO NOT USE YET.       |
    |b         | An Nx3 array of Bx, By, Bz values                 |
    |xyz       | An Nx3 array of the spacecraft's location         |
    |coords    | A string giving the coordinate system of the data |

    Parameters
    ----------
    epoch : datetime.datetime
         A datetime object corresponding to the central time of interest.

    Other Parameters
    ----------------
    tspan : int
        Time span to cover, in hours.  Data returned will always include
        epoch +/- tspan.
    sat : int
        The number, 1-4, indicating which satellite to use.
    kernel_size : int
        Size of median filter window size for smoothing density data.
        MUST be a ODD NUMBER.
    coords : str
        Set the output coordinate system, defaults to GSM.
    debug : bool
        Turn on verbose debug info.
    prune : bool
        Prune data to only epoch +/- tspan.

    Returns
    -------
    data : dict
        A dictionary of values pulled from the Cluster CDFs.  Values 
        are described above.

    '''
    
    import datetime as dt
    import numpy as np
    from spacepy.pycdf import CDF
    from scipy.signal import medfilt
    import spacepy.time as spt 
    import spacepy.coordinates as spc
    
    # Convert tspan to a timedelta:
    tspan = dt.timedelta(hours=tspan)
    
    # Start by grabbing file names.  If epoch+/-tspan crosses the month
    # boundary, add more files.
    fgm_files, cis_files = [], []
    for t in [epoch-tspan, epoch, epoch+tspan]:
        fgm, cis = get_cluster_filename(t, sat=sat)
        # Only add files if the name is unique in the list:
        if fgm not in fgm_files:
            fgm_files.append(fgm)
            cis_files.append(cis)

    if debug:
        print(f'Found {len(fgm_files)} FGM and {len(cis_files)} CIS files')
        print('Data files located: ')
        for f, c in zip(fgm_files,cis_files):
           print(f'\t{f}\n\t{c}\n')

    # Based on satellite number, build variable maps:
    cis_map = {'cis_time':f'Epoch__C{sat}_PP_CIS',
               'dens_h'  :f'N_p__C{sat}_PP_CIS',
               'dens_o'  :f'N_O1__C{sat}_PP_CIS',
               'v_h'     :f'V_p_xyz_gse__C{sat}_PP_CIS'}
    fgm_map = {'fgm_time':f'time_tags__C{sat}_CP_FGM_SPIN',
               'xyz':f'sc_pos_xyz_gse__C{sat}_CP_FGM_SPIN',
               'b'  :f'B_vec_xyz_gse__C{sat}_CP_FGM_SPIN'}
        
    # Build empty data container:
    data = {} # Empty dict
    allkeys = list(cis_map.keys())+list(fgm_map.keys())
    for k in allkeys:
        # Create empty numpy arrays for each variable:
        dtype = object if 'time' in k else None
        data[k] = np.zeros(0,dtype=dtype)

    # Load data into containers:
    for f, c in zip(fgm_files, cis_files):
        # Open CDF files:
        fgm=CDF(f)
        cis=CDF(c)
        # Append data to data container:
        for k in cis_map:
            data[k] = np.append(data[k],cis[cis_map[k]][...])
        for k in fgm_map:
            data[k] = np.append(data[k],fgm[fgm_map[k]][...])

    # Reshape 3D data that got flattened on append:
    for k in ['xyz', 'b', 'v_h']:
        data[k] = data[k].reshape( ( int(data[k].size/3), 3) )

    # Prune data down to only what we need!
    if prune:
        start_Time = epoch - tspan
        end_Time   = epoch + tspan

        # Filter down data to only within time range:
        mask = (data['cis_time']>=start_Time) & (data['cis_time']<=end_Time)    
        for x in ['cis_time', 'dens_h', 'dens_o']:
            data[x] = data[x][mask]
        data['v_h'] = data['v_h'][mask,:]
        
        # Do the same thing as the above and set up a second filter for FGM data   
        mask = (data['fgm_time']>=start_Time) & (data['fgm_time']<=end_Time)
        data['fgm_time'] = data['fgm_time'][mask]
        data['xyz']      = data['xyz'][mask,:]
        data['b']        = data['b'][mask,:]
        
    # Unit conversions: km -> RE
    data['xyz'] /= 6378.16

    # Mask bad data points, apply median filter to density:
    data['dens_h'] = np.ma.masked_values(
        medfilt(data['dens_h'], kernel_size=kernel_size),
        cis[cis_map['dens_h']].attrs['FILLVAL'])
    data['dens_o'] = np.ma.masked_values(
        medfilt(data['dens_o'], kernel_size=kernel_size),
        cis[cis_map['dens_o']].attrs['FILLVAL'])
    data['b' ] = np.ma.masked_values(data['b'],
                                     fgm[fgm_map['b']].attrs['FILLVAL'])

    # Coordinate transformations!
    # Start with time as ticktock object, values as coords objects:
    time = spt.Ticktock(data['fgm_time'], 'ISO')
    xyz  = spc.Coords( data['xyz'], 'GSE', 'car', ticks=time)
    b    = spc.Coords( data['b'],   'GSE', 'car', ticks=time)
    # Now, rotate:
    data['xyz'] = xyz.convert(coords, 'car').data
    data['b'  ] = b.convert(  coords, 'car').data

    # Store the coordinate system we are using:
    data['coords'] = coords
    
    return data

def get_crossing_info(epoch):
    '''
    For a given crossing epoch, calculate and return the following:

    - The crossing time for all Tsgyanenko models under consideration
    - The average density before and after the cluster crossing
      for both H+ and O+.

    Parameters
    ----------
    epoch : datetime.datetime
       The time of the cluster crossing of the plasmasheet.

    Other Parameters
    ----------------
    bsens : float
       The maximum value of Bx in Tsyg results that is considered
       to be a crossing.  Defaults to 1 nT.

    Returns
    -------

    Examples
    --------
    '''

    # Get cluster data associated with epoch:
    data = fetc_cluster_data(epoch, bsens=1)

    # Get crossing time for each Tsyg model.
    # This is a dictionary of datetimes.
    t_times = {} # Container for results.

    for vers in ['T89','T96', 'T01STORM']:
        # DTW notes: this may fail if a model is not available.
        t, b01 = tailwag.gen_sat_tsyg(data, extMag=vers) # get Tsyg data
        loc = np.abs(b01[:,0])==np.abs(b01[:,0]).min() # location of crossing
        # Test to see if crossing data is legit:
        # If there are no points found, then there was no Tsyg model data
        # If the minimum b_field is not reasonably close to zero, then
        # there is a data gap over the crossing.
        if loc.size == 0 or np.abs(b01[loc,0])>bsens:
            t_times[vers] = np.nan
            

    return t_times


def kp_finder_event(date):
    '''
    This function will take in a single datetime 
    object and return the kp index value for it...
    
    Parameters
    ==========
    date : datetime
        Your single datetime event you will use to find the kp index



    Returns
    =======
    Kpindex : float64
        returned value from the array...
    '''
    import spacepy.omni as om
    import spacepy.time as spt 
   
  

    timeticks = spt.Ticktock(date, 'ISO')
    d = om.get_omni(timeticks, dbase='qd1min')

    Kpindex = d['Kp'][0]
    

    return Kpindex


def kp_finder_range(t1, t2):
    '''
    This function takes in two datetime objects
    and finds the range of kp indices as an array
    
    Parameters
    ==========
    t1 : datetime
        The initial date as a datetime object
        
    t2 : datetime
        The final date as a datetime object

    Returns
    =======
    Kprange : array
       the values of the kp in an array...
       
    UTCrnge : array
       the values of the time  in an array...
    '''
    import datetime as dt
    import spacepy.omni as om
    import spacepy.time as spt 
   

    ##################This is the list built with list comprehension that will build the list of the proper 
    ###############dates to be feed into timeticks below....to get all of our data.......
    isotime = [(t1 + dt.timedelta(minutes = x)).isoformat() 
                                                        for x in range(int((t2-t1).total_seconds()//60))]

    timeticks = spt.Ticktock(isotime, 'ISO')

    d = om.get_omni(timeticks, dbase='qd1min')
    
    Kprange = d['Kp'][...]
    UTCrnge = d['UTC'][...]
    
    return Kprange, UTCrnge


def Vz_finder_range(date_str, date_obj, int_hours):
    '''
    This function takes in two datetime objects
    and finds the average value of an array of Vz 
    values over a time interval of specified hours 
    an array
    
    Parameters
    ==========
    date_str : string
        your initial date point as a string object
        
    date_obj : datetime
        your initial date point as a datetime object
        
    int_hours : int
        your interval hours as an int

    Returns
    =======
    avg : float
       the average value of the VZ array
    
    '''

    import datetime as dt
    from datetime import datetime  
    import numpy as np
    import omni as omni
 
   
    ##  Create and print the filename that is going to be passed to get
    ##  OMNI data...
    filename = "OMNI_1min_" + date_str[:7].replace("-", "_")
    print(filename)
    
    
    ####HERE WE CREATE THE START AND END TIMES AS DATETIME OBJECTS, THEY WILL DECREASE AND INCREASE WITH THE LOOP BY N HOURS#######
    start_Time = date_obj - dt.timedelta(hours = int_hours)
    end_Time   = date_obj + dt.timedelta(hours = int_hours)
        
        
    ##  Creating the OMNI dictionary.....
    OMNI_dic = omni.read_ascii(filename)
    
    
    ## Create an empty dictionary that we are going to use to fill
    ## with our usefull masked values......
    data_dic = {}
        
        
    # Create the mask that we are about to use to filter values...
    mask = (OMNI_dic['time']>=start_Time) & (OMNI_dic['time']<=end_Time) 
        
       
    ##getting our values with the mask using list comprehension....
    for x in ['time', 'Vz Velocity']:
        data_dic[x] = OMNI_dic[x][mask]
        data_dic['Vz Velocity'] = OMNI_dic['Vz Velocity'][mask]
    ##compute the average of the Vz Velocity array      
    avg = np.mean(data_dic['Vz Velocity'])   
    
    return avg


def read_ascii(filename):
    '''
    Load data into dictionary.
    '''
    
    import re
    import datetime as dt
    import numpy as np
    from spacepy.datamodel import dmarray

    if filename[-4:] == '.lst':
        datafile   = filename
        formatfile = filename[:-4]+'.fmt'
    elif filename[-4:] == '.fmt':
        formatfile = filename
        datafile   = filename[:-4]+'.lst'
    else:
        formatfile = filename+'.fmt'
        datafile   = filename+'.lst'

    try:
        fmt = open(formatfile, 'r')
    except:
        fmt = False

    info = {}
    var  = []
    flag = []
        
    if fmt:
        raw = fmt.readlines()
        fmt.close()
        # Determine time resolution: hour or minute:
        tres = 'Hour'
        for line in raw:
            if 'Minute' in line: tres='Minute'
        
        # Skip header.
        while True:
            raw.pop(0)
            if tres in raw[0]: break
        last = raw.pop(0)

        # Parse rest of header:
        for l in raw:
            # This regular expression skips the initial digit & space,
            # Lazily grabs the variable name (? turns off greed in this context)
            #
            x = re.search('\d+\s+(.+?)\s+([FI]\d+\.?\d*)', l)
            #x=re.search('\d+\s+(.+?),\s*(\S+)', l) # Old, breaks w/o units.
            if x:
                grps = x.groups()
                # If units available, split and keep:
                var_units = grps[0].split(',')
                if len(var_units)==1: var_units.append('None')
                info[var_units[0]] = ','.join(var_units[1:])
                var.append(var_units[0])
                fmt_code = grps[-1]
                if 'F' in fmt_code:
                    i1, i2 = int(fmt_code[1]), -1*int(fmt_code[-1])
                    flag.append(10**(i1+i2-2) - 10**(i2))
                elif 'I' in fmt_code:
                    i1 = int(fmt_code[-1])-1
                    flag.append(10**i1 -1)
            else:
                raise ValueError(f'Cannot parse header line: {l}')
        
    # Read in data.
    raw = open(datafile, 'r').readlines()
    nLines = len(raw)

    # Create container.
    data = {}
    data['time'] = np.zeros(nLines, dtype=object)
    for k in info:
        data[k] = np.zeros(nLines)
        #data[k] = dmarray(np.zeros(nLines), attrs={'units':info[k]})

    # Now save data.
    t_index = 3 + (tres == 'Minute')*1
    for i, l in enumerate(raw):
        parts = l.split()
        doy = dt.timedelta(days=int(parts[1])-1)
        minute = int(float(parts[3]))*(tres == 'Minute')
        data['time'][i]=dt.datetime(int(parts[0]), 1, 1, int(parts[2]),
                                    minute, 0) + doy
        for j, k in enumerate(var):
            data[k][i] = parts[j+t_index]

    # Mask out bad data.
    for i, k in enumerate(var):
        data[k] = np.ma.masked_greater_equal(data[k], flag[i])

            
    return data

def to_pickle(data_dic, picklename):
    '''
    Parameters
    ----------
    data_dic : dictionary
        DESCRIPTION.
    picklename : TYPE
        DESCRIPTION.

    Returns
    -------
    None....It just makes a pickle bro, -high five-......

    '''
   
    import pickle

    pickle_out = open(picklename, "wb")
    pickle.dump(data_dic, pickle_out)
    pickle_out.close()



 
def fusion(Date_date, Sat_int, interval_hours, add_tsyg=True, add_scatter = True, outdir = 'fusion_plots/'):
    '''
    Parameters
    ----------
    Date_date : datetime object
        Our POI as a datetime object to be fed into various functions
    Sat_int : an int
        The Cluster satellite we will use from 1-4
    interval_hours : int
        A int that will give us the number of hours on each side of our POI

    Returns
    -------
    fig : plot
       
    ax1 : a subplot
       A subplot containg X vs Y positions
    ax2 : a subplot
       A subplot containg X vs Z positions
    ax3 : a subplot
       A subplot containing various Bx magnitudes from cluster and TSG models
    ax4 : a subplot
       A subplot for Oxygen densities 
    ax5 : a subplot
       A subplot for Hydrogen densities 

    '''
    from datetime import timedelta 
    from matplotlib import pyplot as plt
    import os
    import spacepy
    import scatter
    
    from spacepy.plot import applySmartTimeTicks, style, add_arrows
    from spacepy import pybats

    style('spacepy')
    if not os.path.exists(outdir): os.mkdir(outdir)

    # Get string with date & time:
    Date_str = f'{Date_date:%Y-%m-%d %H:%M:%S}'
    
    ### Grabbing the cluster data from the given satellite with our 
    ###   selected interval hours.....
    Data_Dict = fetch_cluster_data(Date_date, Sat_int, interval_hours)
    
    ### FETCH TSYG DATA HERE!       
    if add_tsyg:
        TSYG_T89_time, TSYG_T89_B = gen_sat_tsyg(Data_Dict, extMag='T89', dbase='qd1min');  
        TSYG_T96_time, TSYG_T96_B = gen_sat_tsyg(Data_Dict, extMag='T96', dbase='qd1min');  
        TSYG_T01_time, TSYG_T01_B = gen_sat_tsyg(Data_Dict, extMag='T01STORM', dbase='qd1min'); 
        # Create "helper variable" names:
        T89_time = TSYG_T89_time
        T89_Bx = TSYG_T89_B[:, 0]
        
        T96_time = TSYG_T96_time
        T96_Bx = TSYG_T96_B[:, 0]
        
        T01_time = TSYG_T01_time
        T01_Bx = TSYG_T01_B[:, 0]

    ### Trim down data to largest time range.
    ####  Set up start and stop times so that we can narrow our data down early
    start_Time = Date_date - timedelta(hours = interval_hours)
    end_Time   = Date_date + timedelta(hours = interval_hours)
    
    #### Set up a filter for CIS data to narrow down the stuff from the 
    #### massive dictionary that we just made
    CIS_time = Data_Dict['cis_time']
    maskpro = Data_Dict['dens_h']
    maskoxy = Data_Dict['dens_o']
        
    #### Do the same thing as the above and set up a second   filter for FGM data   
    FGM_time = Data_Dict['fgm_time']
    Bx = Data_Dict['b'][:,0]
    x = Data_Dict['xyz'][:,0]
    y = Data_Dict['xyz'][:,1]
    z = Data_Dict['xyz'][:,2]
    
    #### These values are what we will use to plot the crossing point
    #### On our graphs that we are about to make......
    ind_x = x[round(len(FGM_time)/2)]
    ind_y = y[round(len(FGM_time)/2)]
    ind_z = z[round(len(FGM_time)/2)]
    
    
    ### Adding a section here to get the densitiies for before and after the crossing....
    ### H is first, O is second...
    if add_scatter:
        t_cluster, cross_times, cluster_after, cluster_before = scatter.get_crossing_info(Date_date)
        initial_h = cluster_before[0]
        initial_o = cluster_before[1]
        final_h   = cluster_after[0]
        final_o   = cluster_after[1]
        cross_T89 = cross_times['T89']
        cross_T96 = cross_times['T96']
        cross_T01 = cross_times['T01STORM']
        
       
            
    
    

    ###################################################################################
    #        PLOTTING SECTION
    ###################################################################################
    
    ####LOOP THAT WILL ITERATE 3 TIMES TO GIVE US OUR 3 GRAPHS PER EPOCH........#####
    for n in range(interval_hours, 1, -2):
        print(f"WERE NOW INSIDE THE LOOPING FUNCTION TO SNAG VALUES...... n={n}\n")

        ####  Set up start and stop times so that we can narrow our data down early
        start_Time = Date_date - timedelta(hours = n)
        end_Time   = Date_date + timedelta(hours = n)
        
        #### Create our figure to plot on, give the title and make sure its the 
        #### size of a piece of standard  paper homie
        fig = plt.figure(figsize=(8.5,11))
        fig.suptitle(Date_str + f" $\\pm${n} Hours", fontsize=16)
        
        
        ####Creating the 4 subplots were going to need
        ax1, ax2 = fig.add_subplot(321), fig.add_subplot(322)
        ax3, ax4 = fig.add_subplot(312), fig.add_subplot(313)

        # Cut down orbits to ONLY period of interest +/- N hours
        filter2 = (Data_Dict['fgm_time']>=start_Time) & (Data_Dict['fgm_time']<=end_Time)
        x = Data_Dict['xyz'][filter2,0]
        y = Data_Dict['xyz'][filter2,1]
        z = Data_Dict['xyz'][filter2,2]

        # Axes 1: Orbit in x-y plane:
        line1 = ax1.plot(x, y)        
        ax1.set(xlabel='X in R$_E$', ylabel='Y in R$_E$',
                    title=' X vs Y position')
        add_arrows(line1, n = 5, size = 18, style = '->')     
        ax1.annotate('Tail Crossing', xy=(ind_x,ind_y), xytext = (ind_x + 0.2, ind_y),
                         arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),)  
        ax1.plot(ind_x, ind_y, 'o', ms=3)   
        
        # Axes 2: Orbit in x-z plane:
        line2 = ax2.plot(x, z)
        ax2.set(xlabel='X in R$_E$', ylabel='Z in R$_E$',
                title=' X vs Z position')
        add_arrows(line2, n = 5, size = 18, style = '->')
        ax2.annotate('Tail Crossing', xy=(ind_x,ind_z), xytext = (ind_x + 0.2, ind_z),
                     arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),) 
        ax2.plot(ind_x, ind_z, 'o', ms=3)
    
        # Axes 3: Magnetic field x-component:
        ax3.plot(FGM_time, Bx, 'b-', lw = .5, label=f'Cluster {Sat_int} B$_X$')
        ax3.set(xlabel='time', ylabel='$B_X$ in nT',
                    title='$B_X$ Magnitude')
        ax3.hlines(0, start_Time, end_Time, lw=1, color='black', linestyle = '--')
        ax3.axvline(x = Date_date, ymin = 0, ymax = 1, lw=1, color = 'black', linestyle = '--')
        applySmartTimeTicks(ax3, [start_Time, end_Time])
            
        
        ### Extra Plots### Magnetic field x-compnent from T89, T96, T01STORM
        if add_tsyg:
            ax3.plot(T89_time, T89_Bx, lw=.5, color='green',   label='T89')
            ax3.plot(T96_time, T96_Bx, lw=.5, color='orange',  label='T96')
            ax3.plot(T01_time, T01_Bx, lw=.5, color='crimson', label='T01S')
            ax3.legend(loc='best')
            
        # Axes 4: Composition - Oxygen
        ax4.set_xlabel('time')
        ax4.set_title('Proton and Oxygen Densities')       
        ax4.plot(CIS_time, maskoxy, color = 'g', lw = .5, alpha=0.7)
        ax4.set_ylabel('Oxygen per $cm^{3}$', color = 'g')
        ax4.tick_params(axis='y', labelcolor='g')
        ax4.grid(axis='y')
        #ax4.grid(which='major', axis='both', color = 'g', alpha= 0.2, linestyle = '--')
        ax4.axvline(x = Date_date, ymin = 0, ymax = 1, lw=1, color = 'black', linestyle = '--')

        # Axes 5: Composition - Protons
        ax5 = ax4.twinx()         
        ax5.plot(CIS_time, maskpro, color = 'r', lw=.5)
        ax5.set_ylabel('Protons per $cm^{3}$', color = 'r')
        ax5.tick_params(axis='y', labelcolor='r')
        ax5.grid(axis='y')
        #ax5.grid(b = None, which='major', axis='both', color = 'r', linestyle = '--', alpha= 0.2 )
        applySmartTimeTicks(ax4, [start_Time, end_Time])  
        
        
        if add_scatter:
            ax4.hlines(initial_o, start_Time, Date_date, lw=1, color='g', linestyle = '--')
            ax4.hlines(final_o, Date_date, end_Time, lw=1, color='g', linestyle = '--')
            ax5.hlines(initial_h, start_Time, Date_date, lw=1, color='r', linestyle = '--')
            ax5.hlines(final_h, Date_date, end_Time, lw=1, color='r', linestyle = '--')
            ax3.axvline(x = cross_T89, ymin = 0, ymax = 1, lw=1, color = 'green', linestyle = '--')
            ax3.axvline(x = cross_T96, ymin = 0, ymax = 1, lw=1, color = 'orange', linestyle = '--')
            ax3.axvline(x = cross_T01, ymin = 0, ymax = 1, lw=1, color = 'crimson', linestyle = '--')
                  
        
        
        
        
    
        
        ##### With these 2 blocks were adding planet earth to our postion graphs
        ##### and adding a cool style  to our plots
        ax1 = fig.add_subplot(321)
        ax2 = fig.add_subplot(322)
       
        spacepy.pybats.add_planet(ax1)
        spacepy.pybats.add_planet(ax2)
        
        ####  Doing a tight layout because it will do weird shit otherwise   
        fig.tight_layout(rect=[0, 0, 1, .95])
        plt.show()
        
        ####  HERE WERE SAVING THE GRAPHS TO A FOLDER AND GIVING THE FILE A TITLE THAT HAVE THE DATE
        #####  AND THE INTERVAL HOURS USED
        fig.savefig(outdir + f'fuze_T{Date_date:%Y%m%d_%H%M%S}_int{n:02d}.png')
        
    #### Return dat shit
    return fig, ax1, ax2, ax3, ax4, ax5, n
    
# This next block only executes if the code is run as a
# script and not a module (e.g., "run tailway" vs. "import tailwag").
if __name__ == '__main__':

    import matplotlib.pyplot as plt
    from spacepy.plot import applySmartTimeTicks, style
    from spacepy.pycdf import CDF

    # Turn on good plot style via spacepy:
    style()
    
    # Let's run a quick visual test to make sure things are working:
    # Each call uses a different Tsyganenko model.
    time, mag_t89 = gen_sat_tsyg('sample_data/example_cluster.cdf',
                                 coords='GSE', dbase='qd1min')
    time, mag_t01 = gen_sat_tsyg('sample_data/example_cluster.cdf',
                                 extMag='T01STORM', coords='GSE', dbase='qd1min')

    # Open our real data:
    data = CDF('sample_data/example_cluster.cdf')
    mag_cls = data['B_vec_xyz_gse__C1_CP_FGM_SPIN'][...]
    
    # Create our plot and axes:
    fig  = plt.figure( figsize=[8,8] )
    axes = fig.subplots(3,1)

    # Plot each field component:
    for i,x in enumerate('xyz'):
        ax = axes[i]
        ax.plot(time, mag_cls[:,i], label='Cluster 1', lw=2.0)
        ax.plot(time, mag_t89[:,i], label='T89')
        ax.plot(time, mag_t01[:,i], label='T01s')

        # Apply better time ticks; only label bottom-most axes:
        applySmartTimeTicks(ax, time, dolabel=ax==axes[-1])
        # Add y-labels:
        ax.set_ylabel('B$_{}$ ($nT$)'.format(x.upper()))

    # Add plot title and legend:
    axes[0].legend(loc='best', ncol=3)
    axes[0].set_title('Cluster vs. T89 vs. T01s', size=18)
    
    # Sharpen things up.
    fig.tight_layout()

    # Show the plot:
    plt.show()
