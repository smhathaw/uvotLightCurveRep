#!/usr/bin/python

# import all from pylab
from matplotlib import pylab
from pylab import *

# Script specifically to manipulate 3C454 data after it has been reduced
# for B and J magnitudes

# set object variable
object = '3C454'

# set initial and final observation ID numbers
observation = '31018'
initial = '017'
final = '026'

# define directory to get files from
initial_final = '%s_%s' % (initial, final)
directory = "/Users/Stephanie/Swift/uvotLightCurveOutput/%s/%s" % (object, initial_final)

# define a function to get data from reduced Swift files that will be plot.
# Swift data returns magnitudes, time or error depening on the input variable

def SwiftData(filter, which_array):

    # import fits from astropy module to read data from fits files
    from astropy.io import fits
    
    # set the filter name
    if filter == 'B':
        ff = 'bb'

    # open the fits file that contains all of the magnitudes output
    # using uvotsource on each observation for filter
    hdulist = fits.open("%s/%s_photOut.fits" % (directory, ff))
    # set variable to table of data in hdulist
    tbdata = hdulist[1].data

    # get magnitude values for y values in plot
    mag = []
    for magnitude in tbdata.field('AB_MAG'):
        mag.append(magnitude)

    # get the time that the exposure was taken from the file
    time = []
    for time_seconds_MET in tbdata.field('MET'):
        # convert seconds to days
        time_days_MET = time_seconds_MET/86400
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        # convert to MJD
        time_MJD = 51910 + time_days_MET
        time.append(time_MJD)
        

    # get the y error
    error = []
    for mag_err in tbdata.field('AB_MAG_ERR'):
        error.append(mag_err)
        
    # close the fits file
    hdulist.close()

    if which_array == 'time':
        return time
    elif which_array == 'magnitude':
        return mag
    elif which_array == 'error':
        return error
    else:
        return 'error, no array specified'

# function that inputs J or B and outputs the array of Smarts data specified 
# (magnitude, time or error)

def SmartsData(filter, which_array):
    
    # import urllib2 module to fetch URLs
    import urllib2
    # set the correct url
    url = 'http://www.astro.yale.edu/smarts/glast/tables/%s.tab' % (object)
    # create a data file to write the data to
    open('data.txt','wb').write(urllib2.urlopen(url).read())
    file = open('data.txt','r')
    
    mag = []
    smarts_time = []
    error = []
    row = 0

    # set column numbers for J or B magnitudes
    if filter == 'B':
        # time column
        t = 1
        # magnitude column
        m = 2
        # error column
        e = 3
        # magnitude correction
        c = 0.163
        
    if filter == 'J':
        # time column
        t = 10
        # magnitude column
        m = 11
        # error column
        e = 12
        # magnitude correction
        c = 0
        
    # call SwiftData function to get time array so as to only use 
    # corresponding SMARTS data
    time_range = SwiftData('B','time')
    for line in file:
        # Divide data table (one long string into mangageable lines/lists
        line = line.strip()
        columns = line.split()
        if row > 2:
            # get filter data from table and convert to float
            start = float(columns[t])
            # convert to time from jd to mjd
            mjd = start - 2400000.5
            # conver magnitudes from Jonsons-cousins to AB
            JCmag = float(columns[m])
            ABmag = JCmag - c
            mag_error =  float(columns[e])
    
            # append filter data in date range to lists
            if mjd >= time_range[0] and mjd <= time_range[len(time_range)-1]:
                smarts_time.append(mjd)
                mag.append(ABmag)
                error.append(mag_error)
        row = row + 1 
    file.close

    if which_array == 'time':
        return smarts_time
    elif which_array == 'magnitude':
        return mag
    elif which_array == 'error':
        return error
    else:
        return 'SMARTS data error'

# define a function that will plot a light curve for Swift and Smarts data 
# with Fermi data plotted below.

def LightCurve():
    
    from astropy.io import fits
    import matplotlib.pyplot as plt
    from matplotlib import pylab
    from matplotlib.ticker import ScalarFormatter

    # open and set the fermi data
    # set file name
    if object == '3C454':
        fermiFile = '3C454.3_86400.lc.txt'

    hdulist = fits.open("/Users/Stephanie/Fermi/%s/%s" % (object, fermiFile))
    tbdata = hdulist[1].data
    time = []
    mag = []
    for time_seconds_MET in tbdata.field('START'):
        # convert MET to MJD
        time_days_MET = time_seconds_MET/86400
        time_MJD = 51910 + time_days_MET 
        time.append(time_MJD)
    for flux in tbdata.field('FLUX_100_300000'):
        mag.append(flux)
    errorLAT = []
    for error in tbdata.field('ERROR_100_300000'):
        errorLAT.append(error)

    # call the SwiftData function to get the correct times
    # define filter
    swift_filter1 = 'B'
    plotSWIFTtime1 = SwiftData(swift_filter1,'time')
    
    lat = 0
    plotLATtime = []
    plotLATmag = []
    plotLATerror = []
    for times in time:
        if times >= plotSWIFTtime1[0] and times <= plotSWIFTtime1[len(plotSWIFTtime1)-1]:
            plotLATtime.append(times)
            plotLATmag.append(mag[lat])
            plotLATerror.append(errorLAT[lat])
        lat = lat + 1
    hdulist.close

    # set filters
    smarts_filter1 = 'B'
    smarts_filter2 = 'J'

    # call SwiftData and SmartsData functions to get magnitudes, times and errors
    plotSWIFTmag1 = SwiftData(swift_filter1, 'magnitude')
    plotSWIFTerror1 = SwiftData(swift_filter1, 'error')
    plotSMARTStime1 = SmartsData(smarts_filter1, 'time')
    plotSMARTSmag1 = SmartsData(smarts_filter1, 'magnitude')
    plotSMARTSerror1 = SmartsData(smarts_filter1, 'error')
    plotSMARTStime2 = SmartsData(smarts_filter2, 'time')
    plotSMARTSmag2 = SmartsData(smarts_filter2,'magnitude')
    plotSMARTSerror2 = SmartsData(smarts_filter2,'error')

    # plot the light curves

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)
    ax0.errorbar(plotSWIFTtime1, plotSWIFTmag1, yerr=plotSWIFTerror1, color='red', fmt='o', label='Swift %s' % (swift_filter1))
    ax0.errorbar(plotSMARTStime1, plotSMARTSmag1, yerr=plotSMARTSerror1, color='blue',fmt='o', label='SMARTS %s' % (smarts_filter1))
    ax0.errorbar(plotSMARTStime2, plotSMARTSmag2, yerr=plotSMARTSerror2, color='purple', fmt='o', label='SMARTS %s' % (smarts_filter2))
    
    ax0.set_title(object)
    ax0.legend(bbox_to_anchor=(0.95, 1), loc=2)
    ax0.set_ylabel('Magnitude')
    ax0.set_ylim([16,11])
    
    ax1.errorbar(plotLATtime, plotLATmag, yerr=plotLATerror, color='gray', fmt='o')
    ax1.set_ylabel('.1 - 300 GeV  ' + r'$photons/cm^2/s$')
    plt.subplots_adjust(bottom = 0.1, right = 0.8, top = 0.9)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    xticklabels = ax1.get_xticklabels()
    plt.setp(xticklabels, visible=True)
    fmt=ScalarFormatter(useOffset=False)
    fmt.set_scientific(False)
    gca().xaxis.set_major_formatter(fmt)
    plt.show()


# define a function to plot J mag vs B mag for 3C454
def BJmagPlot():

    # call madules to plot
    import matplotlib.pyplot as plt
    from matplotlib import pylab
    from matplotlib.ticker import ScalarFormatter
    import numpy as np

    # define filters
    filter1 = 'B'
    filter2 = 'J'
    
    # call the functions to get the Smarts data
    smartsMag1 = SmartsData(filter1, 'magnitude')
    smartsMag2 = SmartsData(filter2, 'magnitude')
    smartsTime1 = SmartsData(filter1, 'time')
    smartsTime2 = SmartsData(filter2, 'time')
    smartsError1 = SmartsData(filter1, 'error')
    smartsError2 = SmartsData(filter2, 'error')

    # set arrays so that B,J are plotted for the
    # same times
    plot1smarts = []
    plot2smarts = []
    plotSmartsError1 = []
    plotSmartsError2 = []
    for x1 in range(len(smartsTime1)):
        for x2 in range(len(smartsTime2)):
            if int(smartsTime1[x1]) == int(smartsTime2[x2]):
                plot1smarts.append(smartsMag1[x1])
                plot2smarts.append(smartsMag2[x2])
                plotSmartsError1.append(smartsError1[x1])
                plotSmartsError2.append(smartsError2[x2])
    
    # get the swift data
    swiftMag1 = SwiftData(filter1, 'magnitude')
    swiftTime1 = SwiftData(filter1, 'time')
    swiftError1 = SwiftData(filter1, 'error')
    
    # Set the arrays so that swift B and Smarts J are plotted together
    plot1swift = []
    plot2smarts_swift = []
    plotSwiftError1 = []
    plotSmartsError2_swift = []
    for x1 in range(len(swiftTime1)):
        for x2 in range(len(smartsTime2)):
            if int(swiftTime1[x1]) == int(smartsTime2[x2]):
                plot1swift.append(swiftMag1[x1])
                plot2smarts_swift.append(smartsMag2[x2])
                plotSwiftError1.append(swiftError1[x1])
                plotSmartsError2_swift.append(smartsError2[x2])


    # line with slope of 1
    s = []
    h = []
    for num in range(3):
        s.append(num + 13.6)
        h.append(s[num] - 2.7)

    # make plots
    fig, ax = plt.subplots(nrows = 1)
    ax.errorbar(plot1smarts, plot2smarts, xerr=plotSmartsError1, yerr=plotSmartsError2, color='blue', fmt='o', label='SMARTS %s' % (filter1))
    ax.errorbar(plot1swift, plot2smarts_swift, xerr=plotSwiftError1, yerr=plotSmartsError2_swift, color='red', fmt='o', label='Swift %s' % (filter1))
    ax.set_ylabel('%s Magnitude'  % (filter2))
    ax.set_xlabel('%s Magnitude'  % (filter1))
    ax.set_title(object)
    fit_np = np.polyfit(plot1smarts, plot2smarts, 1)
    ax.plot(plot1smarts, np.polyval(fit_np, plot1smarts),'b-', lw = 2, label='best fit')
    ax.plot( s, h, "g-", lw = 2, label='slope = 1')
    plt.axis([16, 13.5, 13.5, 11])
    axis = plt.gca()
    axis.set_autoscale_on(False)
    legend(loc=2)
    plt.show()

#LightCurve()
BJmagPlot()
