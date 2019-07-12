#!/usr/bin/env python
'''
A module for performing analysis of outflow wagging the tail.
'''


def gen_sat_tsyg(cluster_file, extMag='T89', coords='GSM'):
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

    # Get b-field along orbit:
    b_tsyg = ib.get_Bfield(time, pos, extMag=extMag)

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
    time, mag_t89 = gen_sat_tsyg('sample_data/example_cluster.cdf', coords='GSE')
    time, mag_t01 = gen_sat_tsyg('sample_data/example_cluster.cdf',
                                 extMag='T01STORM', coords='GSE')

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
