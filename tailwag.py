#!/usr/bin/env python
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
                       coords='GSM', debug=False):
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
    
    
def Fusion(Date_date, Sat_int, interval_hours, add_tsyg=True):
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
    import spacepy
    from spacepy.plot import applySmartTimeTicks, style, add_arrows
    from spacepy import pybats

    style('spacepy')

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
    
    
    ####LOOP THAT WILL ITERATE 3 TIMES TO GIVE US OUR 3 GRAPHS PER EPOCH........#####
    for n in range(2, interval_hours + 1, 2):
        print("WERE NOW INSIDE THE LOOPING FUNCTION TO SNAG VALUES......\n")
        ####  Set up start and stop times so that we can narrow our data down early
        start_Time = Date_date - timedelta(hours = n)
        end_Time   = Date_date + timedelta(hours = n)
        
        
        #### Set up a filter for CIS data to narrow down the stuff from the 
        #### massive dictionary that we just made
        filter1 = (Data_Dict['cis_time']>=start_Time) & (Data_Dict['cis_time']<=end_Time)    
        CIS_time = Data_Dict['cis_time'][filter1]
        maskpro = Data_Dict['dens_h'][filter1]
        maskoxy = Data_Dict['dens_o'][filter1]
        #### Do the same thing as the above and set up a second   filter for FGM data   
        filter2 = (Data_Dict['fgm_time']>=start_Time) & (Data_Dict['fgm_time']<=end_Time)
        FGM_time = Data_Dict['fgm_time'][filter2]
        Bx = Data_Dict['b'][filter2,0]
        x = Data_Dict['xyz'][filter2,0]
        y = Data_Dict['xyz'][filter2,1]
        z = Data_Dict['xyz'][filter2,2]
        
        #### Here we are creating the arrays that we are going to be using in 
        #### our graphs, they are require the FGM filter to match things correctly       
        if add_tsyg:
            T89_time = TSYG_T89_time[filter2]
            T89_Bx = TSYG_T89_B[filter2, 0]
            
            T96_time = TSYG_T96_time[filter2]
            T96_Bx = TSYG_T96_B[filter2, 0]
            
            T01_time = TSYG_T01_time[filter2]
            T01_Bx = TSYG_T01_B[filter2, 0]

        #### These values are what we will use to plot the crossing point
        #### On our graphs that we are about to make......
        ind_x = x[round(len(FGM_time)/2)]
        ind_y = y[round(len(FGM_time)/2)]
        ind_z = z[round(len(FGM_time)/2)]
            
        ###################################################################################
        #        PLOTTING SECTION
        ###################################################################################
            
        #### Create our figure to plot on, give the title and make sure its the 
        #### size of a piece of standard  paper homie
        fig = plt.figure(figsize=(8.5,11))
        fig.suptitle(Date_str + " with interval hours: " +  str(n), fontsize=16)
        
        
        ####Creating the 4 subplots were going to need
        ax1, ax2 = fig.add_subplot(321), fig.add_subplot(322)
        ax3, ax4 = fig.add_subplot(312), fig.add_subplot(313)


        # Axes 1: Orbit in x-y plane:
        line1 = ax1.plot(x, y)        
        ax1.set(xlabel='X in R$_E$', ylabel='Y in R$_E$',
                    title=' X vs Y position')
        add_arrows(line1, n = 5, size = 18, style = '->')     
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
        ax3.hlines(0, FGM_time[0], FGM_time[-1], lw=1, color='black', linestyle = '--')
        ax3.axvline(x = Date_date, ymin = 0, ymax = 1, lw=1, color = 'black', linestyle = '--')
        applySmartTimeTicks(ax3, [start_Time, end_Time])
            
        
        ### Extra Plots### Magnetic field x-compnent from T89, T96, T01STORM
        if add_tsyg:
            ax3.plot(T89_time, T89_Bx, lw=.5)
            ax3.plot(T96_time, T96_Bx, lw=.5)
            ax3.plot(T01_time, T01_Bx, lw=.5)
        
        
         # Axes 4: Composition - Oxygen
        ax4.set_xlabel('time')
        ax4.set_title('Proton and Oxygen Densities')       
        ax4.plot(CIS_time, maskoxy, color = 'g', lw = .5, alpha=0.7)
        ax4.set_ylabel('Oxygen per $cm^{3}$', color = 'g')        
        ax4.grid(b = None, which='major', axis='both', color = 'g', alpha= 0.2, linestyle = '--')
        ax4.axvline(x = Date_date, ymin = 0, ymax = 1, lw=1, color = 'black', linestyle = '--')

        # Axes 5: Composition - Protons
        ax5 = ax4.twinx()         
        ax5.plot(CIS_time, maskpro, color = 'r', lw=.5)
        ax5.set_ylabel('Protons per $cm^{3}$', color = 'r')          
        ax5.grid(b = None, which='major', axis='both', color = 'r', linestyle = '--', alpha= 0.2 )
        applySmartTimeTicks(ax4, [start_Time, end_Time])  
    
        
        ##### With these 2 blocks were adding planet earth to our postion graphs
        ##### and adding a cool style  to our plots
        ax1 = fig.add_subplot(321)
        ax2 = fig.add_subplot(322)
       
        spacepy.pybats.add_planet(ax1)
        spacepy.pybats.add_planet(ax2)
        
        ####  Doing a tight layout because it will do weird shit otherwise   
        fig.tight_layout(rect=[0, 0, 1, .95])
        
        #### Return dat shit
        return fig, ax1, ax2, ax3, ax4, ax5
    
    
    
    
    
    
    
    
    
    
    
    

