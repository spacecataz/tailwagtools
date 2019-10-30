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
    
def gen_sat_tsyg(cluster_file, extMag='T89', coords='GSM', dbase='QDhourly'):
    '''
    Given a CDF input file that contains both time and position of a satellite, 
    create magnetic field along that orbit via one of the Tsyganenko empirical
    models.

    Parameters
    ==========
    cluster_file : string
        Path of CDF-formatted cluster input file with orbit location.

    Other Parameters
    ================
    extMag : string
        Name of external magnetic field model.  Choices can be found in 
        documentation for 
        [Spacepy's ONERA interface](https://pythonhosted.org/SpacePy/autosummary/spacepy.irbempy.get_Bfield.html)
        Default is 'T89', or the Tsyganenko '89 model.
    coords : string
        Set the output coordinate system.  Defaults to 'GSM'.
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
    -high-resolution OMNI.
    -verification of results- compare to CCMC test?
    
    '''

    from spacepy.pycdf import CDF
    import spacepy.omni as om
    import spacepy.time as spt 
    import spacepy.coordinates as spc 
    import spacepy.irbempy as ib 

    # Open CDF file, get required variables:
    clus = CDF(cluster_file)
    time = clus['time_tags__C1_CP_FGM_SPIN'][...]      # Extract time
    pos  = clus['sc_pos_xyz_gse__C1_CP_FGM_SPIN'][...] # Extract orbit as nTime x 3D array.

    # Convert CDF values to spacepy-specific classes.
    # These enable coordinate transforms and are required for
    # the ONERA interface.
    time = spt.Ticktock(time, 'ISO')
    pos  = spc.Coords( pos/6371, 'GSE', 'car') # Km->Re.

    # Get input solar wind and QinDenton values:
    QD = om.get_omni(time, dbase=dbase)

    # get_Bfield expects certain labels for all variables...
    QD['dens'] = QD['Den_P']
    QD['velo'] = QD['Vsw']
    
    # Get b-field along orbit:
    b_tsyg = ib.get_Bfield(time, pos, extMag=extMag, omnivals=QD)

    # Convert magnetic field from GEO to GSM.
    if coords!='GEO':
        b_geo = spc.Coords( b_tsyg['Bvec'], 'GEO', 'car', ticks=time)
        b_out = b_geo.convert(coords, 'car')

    # Return time and b-field data.  Strip away extraneous information
    # added by spacepy (via .UTC and .data, which give values only.)
    return(time.UTC, b_out.data)


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
