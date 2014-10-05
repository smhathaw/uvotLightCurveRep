#!/usr/bin/python

# Reduce Swift data to finished light curves. This script follows my guide for
# uvot light curve reduction which follows the information on Swift's wedbsite

# User inputs the object and observations that will be reduced
object = raw_input("Input object name: ") 
# set observation name variables based on the object input
if object == '3C279':
    observation = '35019'
    initial = '136'
    initial_int = int(initial)
    final = '170'
    final_int = int(final)
elif object == '3C454':
    observation = '31018'
    initial = '017'
    final = '026'
    initial_int = int(initial)
    final_int = int(final)
elif object == '2155-304':
    observation = '30795'
    initial = '126'
    final = '150'
    initial_int = int(initial)
    final_int = int(final)

# Set the number of filters to 6
filter_num_int = 6

# define a function (called inside) to name the filters
def filterName(f):
    if f == 0:
        return 'bb'
    elif f == 1:
        return 'm2'
    elif f == 2:
        return 'uu'
    elif f == 3:
        return 'vv'
    elif f == 4:
        return 'w1'
    elif f == 5:
        return 'w2'

# define a function that returns the obsID number as a string
def observationID(obs):
    if len(str(obs)) == 1:
        return "000" + observation + "00" + str(obs)
    elif len(str(obs)) == 2:
        return "000" + observation + "0" + str(obs)
    else:
        return "000" + observation + str(obs)
# define a function (called inside) to name files
def fileName(obsID, file_type, filter):
    if file_type == 'sky_file':
        return "sw" + obsID + "u" + filter + "_sk.img"
    elif file_type == 'exp_file':
        return "sw" + obsID + "u" + filter + "_ex.img"
    elif file_type == 'sky_file_gz':
        return "sw" + obsID + "u" + filter + "_sk.img.gz"
    elif file_type == 'exp_file_gz':
        return "sw" + obsID + "u" + filter + "_ex.img.gz"

# define a function (called outside) to create an output directory
def directoryOut():
    # create a directory for output
    import os # import os module to use system for command line input
    initial_final = "%s_%s" % (initial, final) # named for final and initial obs
    os.system("mkdir /Users/Stephanie/Swift/uvotLightCurveOutput/%s" % (object))
    directory_out = "/Users/Stephanie/Swift/uvotLightCurveOutput/%s/%s" % (object, initial_final)
    os.system("mkdir %s" % (directory_out))

# define directory out outside of function
initial_final = "%s_%s" % (initial, final)
directory_out = "/Users/Stephanie/Swift/uvotLightCurveOutput/%s/%s" % (object, initial_final)

# define a function to unzip files 
def unzipFun():
    # For loop to execute once for each filter
    for x in range(filter_num_int): 
        # call filterName function to name filter
        filter = filterName(x)

        # Loop over all of the observation ID numbers
        final_int_1 = final_int + 1 # Add 1 to the final observation number
        for y in range(initial_int,final_int_1):
            # call the function to get the observation ID number
            obsID = observationID(y)
            # call fileName function to name files            
            sky_file = fileName(obsID, 'sky_file', filter)
            exp_file = fileName(obsID, 'exp_file', filter)
            sky_file_gz = fileName(obsID, 'sky_file_gz', filter)
            exp_file_gz = fileName(obsID, 'exp_file_gz', filter)
        
            # Run gunzip on files
            import os # call the os module
            # unzip sky image file
            os.system("gunzip /Users/Stephanie/Swift/%s/%s/uvot/image/%s" % (object, obsID, sky_file_gz)) # call os.system to execute the gunzip command
            # unzip the exposure map file
            os.system("gunzip /Users/Stephanie/Swift/%s/%s/uvot/image/%s" % (object, obsID, exp_file_gz)) # call os.system to execute the gunzip command
        

# function to name the output file for uvotsource
def uvotsourceOutFile(filter):
    return "%s_photOut.fits" % (filter)


# define a function to reduce data (called outside)
def reduction():
    # For loop to execute once for each filter
    for x in range(filter_num_int):
        # call filterName function to name filter
        filter = filterName(x)
    
        # Loop over all of the observation ID numbers
        final_int_1 = final_int + 1 # Add 1 to the final observation number
        for y in range(initial_int,final_int_1):
            # call function to get observation ID
            obsID = observationID(y)
            # call fileName function to name files
            sky_file = fileName(obsID, 'sky_file', filter)
            exp_file = fileName(obsID, 'exp_file', filter)

            # Check if the file that will have data reduction done exists
            import os.path # call path from os module
            file_exists = os.path.isfile("/Users/Stephanie/Swift/%s/%s/uvot/image/%s" % (object, obsID, sky_file))
        
            # if the file exists begin data reduction
            if file_exists == True:

                # First check the aspect correction (if aspect correction
                # is not equal to direct then those extensions should be excluded)

                # import the fits file reader from astropy
                from astropy.io import fits
                # open the file to check ASPCORR header
                hdulist1 = fits.open('/Users/Stephanie/Swift/%s/%s/uvot/image/%s' % (object, obsID, sky_file))
        
                extensions = [] # create an empty list to fill with extensions
                for extension in hdulist1: # loop over extensions
                    extensions.append(extension) # add extension to list
                num_extensions = len(extensions) # get the number of extensions

                exclude_ext = [] # create empty list for exluded extensions
                for ext in range(1, num_extensions): # loop through extensions
                    # Set aspcorr to value for ASPCORR from the file
                    aspcorr = hdulist1[ext].header['ASPCORR']
                    # put extensions with ASPCORR not set to DIRECT on the
                    # exclude list
                    if aspcorr != "DIRECT":
                        ext_str = str(ext) # converst extension to string
                        exclude_ext.append(ext_str)
                # make exlude list a string
                exclude_str = ", ".join(exclude_ext) # make exclude list a string
                hdulist1.close() # close the file

                # Run uvotimsum to coadd extensions in a single file
                import os # import the os module
                # set the name of summed sky image file
                usum_sky_file = "usumu" + filter + "_" + str(y) + "_sky.fits"
                # set the name of the summed exposure map file
                usum_ex_file = "usumu" + filter + "_" + str(y) + "_ex.fits"
            
                # use system fom os module to run the ftool command uvotimsum
                # first coadd the sky image file
                os.system("uvotimsum infile = /Users/Stephanie/Swift/%s/%s/uvot/image/%s outfile = %s/%s exclude = %s clobber = yes" % (object, obsID, sky_file, directory_out, usum_sky_file, exclude_str))
                # then coadd the exposure map image file
                os.system("uvotimsum infile = /Users/Stephanie/Swift/%s/%s/uvot/image/%s outfile = %s/%s exclude = %s clobber = yes" % (object, obsID, exp_file, directory_out, usum_ex_file, exclude_str))

                # Run uvotsource to do photometry using data from the sky image
                # file the exposure map file and calibration data from the
                # CALBD datbase to output magnitudes from the source at specific
                # times (in MET)

                # call uvotsourceOutFile to name the output file for uvotsource
                photometry_out = uvotsourceOutFile(filter)

                # use system to run the ftool command uvotsource to do photometry
                # Note that source region and background region files should have
                # been created and saved in the Swift/"object" directory. To find
                # the correct coordinates uvotdetect should be run a summed file.
                # I compared the source coordinates to those named in the Swift
                # database.
            
                # uvotsource runs on one extension of a sky image file, however
                # here I run it on my summed file for each observation of a single
                # filter (there is only one extension). Every output table is put
                # into the same file for every observation of a single filter
                os.system("uvotsource image = %s/%s srcreg = /Users/Stephanie/Swift/%s/source5.reg bkgreg = /Users/Stephanie/Swift/%s/backgr.reg sigma = 5 expfile = %s/%s outfile = %s/%s" % (directory_out, usum_sky_file, object, object,directory_out, usum_ex_file, directory_out, photometry_out))


from matplotlib import pylab
from pylab import *

# Function to plot single light curves for a filter (called outside)
def plotSingle(filter):
    
    # Determine the actual filter name
    if filter == 'bb':
        filter_correct = 'B'
    elif filter == 'uu':
        filter_correct = 'U'
    elif filter == 'vv':
        filter_correct = 'V'
    elif filter == 'm2':
        filter_correct = 'UVM2'
    elif filter == 'w1':
        filter_correct = 'UVW1'
    elif filter == 'w2':
        filter_correct = 'UVW2'

    # Plot the magnitudes vs swift time (mission elapsed time) using the
    # magnitude and MET output by uvotsource 
        
    # import fits to read fits files from astropy.io
    from astropy.io import fits 
    # import pyplot from matplotlib to plot data
    import matplotlib.pyplot as plt
    # import ScalarFormatter from matplotlib.ticker to correct the 
    # x-axis ticker values
    from matplotlib.ticker import ScalarFormatter
    
    # call the function to name the output file by uvotsource
    photometry_out = uvotsourceOutFile(filter)

    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for a single filter
    hdulist2 = fits.open("%s/%s" % (directory_out, photometry_out))
    tbdata =  hdulist2[1].data # set variable to table of data in hdulist
    y = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdata.field('AB_MAG'):
        y.append(magnitude) # append magnitude to y list
    x = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdata.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        ##### NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        x.append(time_MJD) # append time to x list
        
    y_error = [] # empty list for y error
    for mag_err in tbdata.field('AB_MAG_ERR'):
        y_error.append(mag_err)
        
    hdulist2.close() # close the fits file
        
    # Plot magnitude vs time (x vs y) in a scatter plot
    fig, ax = plt.subplots()
    ax.errorbar(x,y, yerr=y_error, color='purple', fmt='o')
    plt.xlabel('Time in MJD')
    plt.ylabel('%s Magnitude' % (filter_correct))
    plt.title(object)
    fmt=matplotlib.ticker.ScalarFormatter(useOffset=False)
    fmt.set_scientific(False)
    gca().xaxis.set_major_formatter(fmt)
    plt.show()
        
# define a function to plot all filters
def plotAllFilters():
    
    # Plot the magnitudes vs swift time (mission elapsed time) using the    
    # magnitude and MET output by uvotsource                                
    
    # import fits to read fits files from astropy.io                        
    from astropy.io import fits
    # import pyplot from matplotlib to plot data                            
    import matplotlib.pyplot as plt
    # import ScalarFormatter from matplotlib.ticker to correct the          
    # x-axis ticker values                                                  
    from matplotlib.ticker import ScalarFormatter
    from matplotlib import pylab
    
    # open the fits file that contains all of the magnitudes ouput          
    # using uvotsource on each observation for b filter              
    hdulistB = fits.open("%s/bb_photOut.fits" % (directory_out))
    tbdataB =  hdulistB[1].data # set variable to table of data in hdulist   
    yB = [] # empty list for y values in plot                                
        # get magnitude values from table in fits file                          
    for magnitude in tbdataB.field('AB_MAG'):
        yB.append(magnitude - 1) # append magnitude to y list
    xB = [] # empty list for x values
    # get the time that the exposure was taken from the file                
    for time_seconds_MET in tbdataB.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days    
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD                    
        xB.append(time_MJD) # append time to x list
    y_errorB = [] # empty list for y error
    for mag_err in tbdataB.field('AB_MAG_ERR'):
        y_errorB.append(mag_err) # append y_error to y_error list
    hdulistB.close() # close the fits file

    print 'b mag = ' , yB
    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for m2 filter                    
    hdulistM2 = fits.open("%s/m2_photOut.fits" % (directory_out))
    tbdataM2 =  hdulistM2[1].data # set variable to table of data in hdulist  
    yM2 = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataM2.field('AB_MAG'):
        yM2.append(magnitude + 1) # append magnitude to y list
    xM2 = [] # empty list for x values
    # get the time that the exposure was taken from the file     
    for time_seconds_MET in tbdataM2.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days   
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)   
        time_MJD = 51910 + time_days_MET # conver to MJD                   
        xM2.append(time_MJD) # append time to x list                         
    y_errorM2 = [] # empty list for y error                  
    for mag_err in tbdataM2.field('AB_MAG_ERR'):
        y_errorM2.append(mag_err) # append y_error to y_error list             
    hdulistM2.close() # close the fits file                               
    
    # open the fits file that contains all of the magnitudes ouput         
    # using uvotsource on each observation for u filter                      
    hdulistU = fits.open("%s/uu_photOut.fits" % (directory_out))
    tbdataU =  hdulistU[1].data # set variable to table of data in hdulist 
    yU = [] # empty list for y values in plot                                
    # get magnitude values from table in fits file                           
    for magnitude in tbdataU.field('AB_MAG'):
        yU.append(magnitude) # append magnitude to y list                   
    xU = [] # empty list for x values                                    
    # get the time that the exposure was taken from the file             
    for time_seconds_MET in tbdataU.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days  
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)   
        time_MJD = 51910 + time_days_MET # conver to MJD                
        xU.append(time_MJD) # append time to x list                         
    y_errorU = [] # empty list for y error                            
    for mag_err in tbdataU.field('AB_MAG_ERR'):
        y_errorU.append(mag_err) # append y_error to y_error list
    hdulistU.close() # close the fits file                               
    
    # open the fits file that contains all of the magnitudes ouput            
    # using uvotsource on each observation for V filter                     
    hdulistV = fits.open("%s/vv_photOut.fits" % (directory_out,))
    tbdataV =  hdulistV[1].data # set variable to table of data in hdulist
    yV = [] # empty list for y values in plot                                
    # get magnitude values from table in fits file        
    for magnitude in tbdataV.field('AB_MAG'):
        yV.append(magnitude - 2) # append magnitude to y list
    xV = [] # empty list for x values    
    # get the time that the exposure was taken from the file           
    for time_seconds_MET in tbdataV.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days    
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD  
        xV.append(time_MJD) # append time to x list 
    y_errorV = [] # empty list for y error          
    for mag_err in tbdataV.field('AB_MAG_ERR'):
        y_errorV.append(mag_err) # append y_error to y_error list  
    hdulistV.close() # close the fits file                               
    
    # open the fits file that contains all of the magnitudes ouput            
    # using uvotsource on each observation for w1 filter                 
    hdulistW1 = fits.open("%s/w1_photOut.fits" % (directory_out))
    tbdataW1 =  hdulistW1[1].data # set variable to table of data in hdulist  
    yW1 = [] # empty list for y values in plot                             
    # get magnitude values from table in fits file                   
    for magnitude in tbdataW1.field('AB_MAG'):
        yW1.append(magnitude) # append magnitude to y list                
    xW1 = [] # empty list for x values           
    # get the time that the exposure was taken from the file           
    for time_seconds_MET in tbdataW1.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days    
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD) 
        time_MJD = 51910 + time_days_MET # conver to MJD 
        xW1.append(time_MJD) # append time to x list   
    y_errorW1 = [] # empty list for y error                                   
    for mag_err in tbdataW1.field('AB_MAG_ERR'):
        y_errorW1.append(mag_err) # append y_error to y_error list        
    hdulistW1.close() # close the fits file
                               
    # open the fits file that contains all of the magnitudes ouput         
    # using uvotsource on each observation for w2 filter                
    hdulistW2 = fits.open("%s/w2_photOut.fits" % (directory_out,))
    tbdataW2 =  hdulistW2[1].data # set variable to table of data in hdulist
    yW2 = [] # empty list for y values in plot                             
    # get magnitude values from table in fits file            
    for magnitude in tbdataW2.field('AB_MAG'):
        yW2.append(magnitude + 2) # append magnitude to y list
    xW2 = [] # empty list for x values                                   
    # get the time that the exposure was taken from the file                   
    for time_seconds_MET in tbdataW2.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days   
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD     
        xW2.append(time_MJD) # append time to x list  
    y_errorW2 = [] # empty list for y error                                    
    for mag_err in tbdataW2.field('AB_MAG_ERR'):
        y_errorW2.append(mag_err) # append y_error to y_error list       
    hdulistW2.close() # close the fits file                               
    
    # Get B and V Data from SMARTS
    import urllib2 # import module for fetching URLs  
    # set the correct url
    url = 'http://www.astro.yale.edu/smarts/glast/tables/%s.tab' % (object)    
    # create a data file to right the tables from url to 
    open('data.txt', 'wb').write(urllib2.urlopen(url).read())
    f = open('data.txt', 'r') 
    # create empty lists
    smartsB_y = []
    smartsB_x = []
    smartsB_error = []
    b = 0 # set b to 0 
    smartsV_y = []
    smartsV_x = []
    smartsV_error = []
    v = 0 # set v to 0  
    for line in f:
        # divide data table (one long string) into managable lines/lists
        line = line.strip()
        columns = line.split()
        if b > 2:
            # get B filter data from the table and convert to float
            smarts_B_Start = float(columns[1])
            # convert time from jd to mjd
            smarts_B_mjd = smarts_B_Start - 2400000.5
            # convert magnitudes from Johnsons to AB
            smarts_B_mag = float(columns[2])
            smarts_B_ABmag = smarts_B_mag - 0.163
            smarts_B_err = float(columns[3])
        # append B filter data in date range to above list
            if smarts_B_mjd >= xB[0] and smarts_B_mjd <= xB[len(xB) - 1]:
                smartsB_x.append(smarts_B_mjd)
                smartsB_y.append(smarts_B_ABmag - 1)
                smartsB_error.append(smarts_B_err)
        b = b + 1
        if v > 2:
            # get V filter data from the table and convert to float    
            smarts_V_Start = float(columns[4])
            # convert time from jd to mjd 
            smarts_V_mjd = smarts_V_Start - 2400000.5
            smarts_V_mag = float(columns[5])
            # convert from Johnson to AB magnitude
            smarts_V_ABmag = smarts_V_mag - 0.044
            smarts_V_err = float(columns[6])
        # append V filter data in date range to above list                     
            if smarts_V_mjd >= xV[0] and smarts_V_mjd <= xV[len(xV) - 1]:
                smartsV_x.append(smarts_V_mjd)
                smartsV_y.append(smarts_V_ABmag - 2)
                smartsV_error.append(smarts_V_err)
        v = v + 1
    
    f.close() # close file
    
    # extract Fermi data from previously downloaded files
    # open the fits file previously downloaded from the Fermi website
    
    # set the file name depending on the object
    if object == '3C279':
        filename = 'J123939+044409_604800.lc.txt'
    elif object == '3C454':
        filename = '3C454.3_86400.lc.txt'
    elif object == '2155-304':
        filename = 'PKS2155-304_86400.lc.txt'

    # open first file and add energy and time data to list
    hdulistLAT = fits.open("/Users/Stephanie/Fermi/%s/%s" % (object, filename))
    tbdataLAT = hdulistLAT[1].data # set variable to table of data in hdulist
    xLAT = [] # empty list for the x values (time)
    yLAT = [] # empty list for y values (magnitudes)
    for time_seconds_MET in tbdataLAT.field('START'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xLAT.append(time_MJD)
    for flux in tbdataLAT.field('FLUX_100_300000'):
        yLAT.append(flux)
    # get the errors
    errorLAT = []
    for error  in tbdataLAT.field('ERROR_100_300000'):
        errorLAT.append(error)

    # Plot only the same portion as Swift and SMARTS data
    lat = 0 # iteration variable
    plotLATx = []
    plotLATy = []
    plotLATerror = []
    for times in xLAT:
        if times >= xM2[0] and times <= xM2[len(xM2) - 1]:
            plotLATx.append(times)
            plotLATy.append(yLAT[lat])
            plotLATerror.append(errorLAT[lat])
        lat = lat + 1
    hdulistLAT.close()

    #print 'xLAT = ', plotLATx, 'yLAT = ', plotLATy, 'errorLAT = ' , plotLATerror
    print 'timeB = ', xB, ' B = ', xB, ' timeV = ' , xV, ' Y = ', yV

    # Plot magnitude vs time (x vs y) in a single scatter plot if prompted
    #fig, (ax0, ax1) = plt.subplots(nrows=2)
    f = figure()

    ax0 = subplot(211)
    ax0.errorbar(xB,yB, yerr=y_errorB, color='blue', fmt='o', label='B')
    ax0.errorbar(xM2,yM2, yerr=y_errorM2, color='purple', fmt='o', label='UVM2')
    ax0.errorbar(xU,yU, yerr=y_errorU, color='green', fmt='o', label='U')
    ax0.errorbar(xV,yV, yerr=y_errorV, color='red', fmt='o', label='V')
    ax0.errorbar(xW1,yW1, yerr=y_errorW1, color='cyan', fmt='o', label='UVW1')
    ax0.errorbar(xW2,yW2, yerr=y_errorW2, color='orange', fmt='o', label='UVW2')
    ax0.errorbar(smartsB_x,smartsB_y, yerr=smartsB_error, color='black', fmt='o', label='Smarts B')
    ax0.errorbar(smartsV_x,smartsV_y, yerr=smartsV_error, color='magenta', fmt='o', label='Smarts V')
    ax0.set_title(object)
    ax0.legend(bbox_to_anchor=(0.95, 1.1), loc=2)
    ax0.set_ylabel('Magnitude')
    #xticklabels = ax1.get_xticklabels()

    ax1 = subplot(212)
    ax1.errorbar(plotLATx, plotLATy, yerr = plotLATerror, color = 'gray', fmt='o')
    ax1.set_ylabel('.1 - 300 GeV ' + r'$photons/cm^2/s$')
    ax1.set_xlabel('MJD')

    plt.subplots_adjust(bottom = 0.1, right = 0.8, top = 0.9)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    show()

# function to plot the b and v filters
def plotBvV():
    # Plot the magnitudes vs swift time (mission elapsed time) using the
    # magnitude and MET output by uvotsource
    
    # import fits to read fits files from astropy.io
    from astropy.io import fits
    # import pyplot from matplotlib to plot data
    import matplotlib.pyplot as plt
    # import ScalarFormatter from matplotlib.ticker to correct the
    # x-axis ticker values
    from matplotlib.ticker import ScalarFormatter
    
    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for b filter
    hdulistB = fits.open("%s/bb_photOut.fits" % (directory_out))
    tbdataB =  hdulistB[1].data # set variable to table of data in hdulist
    yB = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataB.field('AB_MAG'):
        yB.append(magnitude - 1) # append magnitude to y list
    xB = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataB.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xB.append(time_MJD) # append time to x list
    y_errorB = [] # empty list for y error
    for mag_err in tbdataB.field('AB_MAG_ERR'):
        y_errorB.append(mag_err) # append y_error to y_error list
    hdulistB.close() # close the fits file

    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for V filter
    hdulistV = fits.open("%s/vv_photOut.fits" % (directory_out,))
    tbdataV =  hdulistV[1].data # set variable to table of data in hdulist
    yV = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataV.field('AB_MAG'):
        yV.append(magnitude) # append magnitude to y list
    xV = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataV.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xV.append(time_MJD) # append time to x list
    y_errorV = [] # empty list for y error
    for mag_err in tbdataV.field('AB_MAG_ERR'):
        y_errorV.append(mag_err) # append y_error to y_error list
    hdulistV.close() # close the fits file

    #   bvt = 0
    # for item in xB:
    #    print 'Btime = ', xB[bvt],  ' Vtime = ', xV[bvt], ' B= ', yB[bvt], ' V = ', yV[bvt]
    #    bvt = bvt + 1

    # make lists of what will be plot
    plotB=[]
    plotV=[]
    plotBerror=[]
    plotVerror=[]
    for x in range(len(xB)):
        for y in range(len(xV)):
            if int(xB[x]) == int(xV[y]):
                plotB.append(yB[x])
                plotV.append(yV[y])
                plotBerror.append(y_errorB[x])
                plotVerror.append(y_errorV[y])


    # print 'B = ', plotB, ' V = ', plotV
    
    s = []
    h = []
    for num in range(3):
        s.append(num + 13)
        h.append(s[num] + .65)
    
    # plot B magnitude vs V magnitude
    fig, ax = plt.subplots(nrows = 1)
    ax.errorbar(plotB, plotV, xerr=plotBerror, yerr=plotVerror, color='blue', fmt='o', label='B v V')
    ax.set_ylabel('V Magnitude')
    ax.set_xlabel('B Magnitude')
    ax.set_title(object)
    xticklabels = ax.get_xticklabels()
    plt.setp(xticklabels, visible=True)
    fmt=matplotlib.ticker.ScalarFormatter(useOffset=True)
    fmt.set_scientific(False)
    gca().xaxis.set_major_formatter(fmt)
    fit_np = np.polyfit(plotB, plotV, 1)
    ax.plot(plotB, np.polyval(fit_np, plotB), "r-", lw = 2, label='best fit')
    ax.plot( s, h, "g-", lw = 2, label='one to one line')
    plt.legend(loc=2)
    show()

# function to plot the b and v filters
def plotBvM2():
    # Plot the magnitudes vs swift time (mission elapsed time) using the
    # magnitude and MET output by uvotsource
    
    # import fits to read fits files from astropy.io
    from astropy.io import fits
    # import pyplot from matplotlib to plot data
    import matplotlib.pyplot as plt
    # import ScalarFormatter from matplotlib.ticker to correct the
    # x-axis ticker values
    from matplotlib.ticker import ScalarFormatter
    
    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for b filter
    hdulistB = fits.open("%s/bb_photOut.fits" % (directory_out))
    tbdataB =  hdulistB[1].data # set variable to table of data in hdulist
    yB = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataB.field('AB_MAG'):
        yB.append(magnitude - 1) # append magnitude to y list
    xB = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataB.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xB.append(time_MJD) # append time to x list
    y_errorB = [] # empty list for y error
    for mag_err in tbdataB.field('AB_MAG_ERR'):
        y_errorB.append(mag_err) # append y_error to y_error list
    hdulistB.close() # close the fits file
    
    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for m2 filter
    hdulistM2 = fits.open("%s/m2_photOut.fits" % (directory_out))
    tbdataM2 =  hdulistM2[1].data # set variable to table of data in hdulist
    yM2 = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataM2.field('AB_MAG'):
        yM2.append(magnitude) # append magnitude to y list
    xM2 = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataM2.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xM2.append(time_MJD) # append time to x list
    y_errorM2 = [] # empty list for y error
    for mag_err in tbdataM2.field('AB_MAG_ERR'):
        y_errorM2.append(mag_err) # append y_error to y_error list
    hdulistM2.close() # close the fits file

    # make lists of what will be plot
    plotB=[]
    plotM2=[]
    plotBerror=[]
    plotM2error=[]
    for x in range(len(xB)):
        for y in range(len(xM2)):
            if int(xB[x]) == int(xM2[y]):
                plotB.append(yB[x])
                plotM2.append(yM2[y])
                plotBerror.append(y_errorB[x])
                plotM2error.append(y_errorM2[y])

    s = []
    h = []
    for num in range(2):
        s.append(num + 12)
        h.append(s[num] + 1.7)


    # plot B magnitude vs M2 magnitude
    fig, ax = plt.subplots(nrows = 1)
    ax.errorbar(plotB, plotM2, xerr=plotBerror, yerr=plotM2error, color='blue', fmt='o', label='B v M2')
    ax.set_ylabel('M2 Magnitude')
    ax.set_xlabel('B Magnitude')
    ax.set_title(object)
    xticklabels = ax.get_xticklabels()
    plt.setp(xticklabels, visible=True)
    fmt=matplotlib.ticker.ScalarFormatter(useOffset=True)
    fmt.set_scientific(False)
    gca().xaxis.set_major_formatter(fmt)
    fit_np = np.polyfit(plotB, plotM2, 1)
    ax.plot(plotB, np.polyval(fit_np, plotB), "r-", lw = 2, label='best fit')
    ax.plot( s, h, "g-", lw = 2, label='line of slope 1')
    legend(loc = 2)
    show()

#### Function to plot all data, except the Fermi Data
def plotAllOptUVFilters():
    
    # Plot the magnitudes vs swift time (mission elapsed time) using the
    # magnitude and MET output by uvotsource
    
    # import fits to read fits files from astropy.io
    from astropy.io import fits
    # import pyplot from matplotlib to plot data
    import matplotlib.pyplot as plt
    # import ScalarFormatter from matplotlib.ticker to correct the
    # x-axis ticker values
    from matplotlib.ticker import ScalarFormatter
    
    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for b filter
    hdulistB = fits.open("%s/bb_photOut.fits" % (directory_out))
    tbdataB =  hdulistB[1].data # set variable to table of data in hdulist
    yB = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataB.field('AB_MAG'):
        yB.append(magnitude - 1) # append magnitude to y list
    xB = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataB.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xB.append(time_MJD) # append time to x list
    y_errorB = [] # empty list for y error
    for mag_err in tbdataB.field('AB_MAG_ERR'):
        y_errorB.append(mag_err) # append y_error to y_error list
    hdulistB.close() # close the fits file
    
    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for m2 filter
    hdulistM2 = fits.open("%s/m2_photOut.fits" % (directory_out))
    tbdataM2 =  hdulistM2[1].data # set variable to table of data in hdulist
    yM2 = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataM2.field('AB_MAG'):
        yM2.append(magnitude + 1) # append magnitude to y list
    xM2 = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataM2.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xM2.append(time_MJD) # append time to x list
    y_errorM2 = [] # empty list for y error
    for mag_err in tbdataM2.field('AB_MAG_ERR'):
        y_errorM2.append(mag_err) # append y_error to y_error list
    hdulistM2.close() # close the fits file
    
    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for u filter
    hdulistU = fits.open("%s/uu_photOut.fits" % (directory_out))
    tbdataU =  hdulistU[1].data # set variable to table of data in hdulist
    yU = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataU.field('AB_MAG'):
        yU.append(magnitude) # append magnitude to y list
    xU = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataU.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xU.append(time_MJD) # append time to x list
    y_errorU = [] # empty list for y error
    for mag_err in tbdataU.field('AB_MAG_ERR'):
        y_errorU.append(mag_err) # append y_error to y_error list
    hdulistU.close() # close the fits file
    
    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for V filter
    hdulistV = fits.open("%s/vv_photOut.fits" % (directory_out,))
    tbdataV =  hdulistV[1].data # set variable to table of data in hdulist
    yV = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataV.field('AB_MAG'):
        yV.append(magnitude - 2) # append magnitude to y list
    xV = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataV.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xV.append(time_MJD) # append time to x list
    y_errorV = [] # empty list for y error
    for mag_err in tbdataV.field('AB_MAG_ERR'):
        y_errorV.append(mag_err) # append y_error to y_error list
    hdulistV.close() # close the fits file
    
    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for w1 filter
    hdulistW1 = fits.open("%s/w1_photOut.fits" % (directory_out))
    tbdataW1 =  hdulistW1[1].data # set variable to table of data in hdulist
    yW1 = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataW1.field('AB_MAG'):
        yW1.append(magnitude) # append magnitude to y list
    xW1 = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataW1.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xW1.append(time_MJD) # append time to x list
    y_errorW1 = [] # empty list for y error
    for mag_err in tbdataW1.field('AB_MAG_ERR'):
        y_errorW1.append(mag_err) # append y_error to y_error list
    hdulistW1.close() # close the fits file
    
    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for w2 filter
    hdulistW2 = fits.open("%s/w2_photOut.fits" % (directory_out,))
    tbdataW2 =  hdulistW2[1].data # set variable to table of data in hdulist
    yW2 = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataW2.field('AB_MAG'):
        yW2.append(magnitude + 2) # append magnitude to y list
    xW2 = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataW2.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xW2.append(time_MJD) # append time to x list
    y_errorW2 = [] # empty list for y error
    for mag_err in tbdataW2.field('AB_MAG_ERR'):
        y_errorW2.append(mag_err) # append y_error to y_error list
    hdulistW2.close() # close the fits file
    
    # Get B and V Data from SMARTS
    import urllib2 # import module for fetching URLs
    # set the correct url
    url = 'http://www.astro.yale.edu/smarts/glast/tables/%s.tab' % (object)
    # create a data file to right the tables from url to
    open('data.txt', 'wb').write(urllib2.urlopen(url).read())
    f = open('data.txt', 'r')
    # create empty lists
    smartsB_y = []
    smartsB_x = []
    smartsB_error = []
    b = 0 # set b to 0
    smartsV_y = []
    smartsV_x = []
    smartsV_error = []
    v = 0 # set v to 0
    for line in f:
        # divide data table (one long string) into managable lines/lists
        line = line.strip()
        columns = line.split()
        if b > 2:
            # get B filter data from the table and convert to float
            smarts_B_Start = float(columns[1])
            # convert time from jd to mjd
            smarts_B_mjd = smarts_B_Start - 2400000.5
            # convert magnitudes from Johnsons to AB
            smarts_B_mag = float(columns[2])
            smarts_B_ABmag = smarts_B_mag - 0.163
            smarts_B_err = float(columns[3])
            # append B filter data in date range to above list
            if smarts_B_mjd >= xB[0] and smarts_B_mjd <= xB[len(xB) - 1]:
                smartsB_x.append(smarts_B_mjd)
                smartsB_y.append(smarts_B_ABmag - 1)
                smartsB_error.append(smarts_B_err)
        b = b + 1
        if v > 2:
            # get V filter data from the table and convert to float
            smarts_V_Start = float(columns[4])
            # convert time from jd to mjd
            smarts_V_mjd = smarts_V_Start - 2400000.5
            smarts_V_mag = float(columns[5])
            # convert from Johnson to AB magnitude
            smarts_V_ABmag = smarts_V_mag - 0.044
            smarts_V_err = float(columns[6])
            # append V filter data in date range to above list
            if smarts_V_mjd >= xV[0] and smarts_V_mjd <= xV[len(xV) - 1]:
                smartsV_x.append(smarts_V_mjd)
                smartsV_y.append(smarts_V_ABmag - 2)
                smartsV_error.append(smarts_V_err)
        v = v + 1
    
    f.close() # close file

    fig, ax0 = plt.subplots(nrows=1)
    ax0.errorbar(xB,yB, yerr=y_errorB, color='blue', fmt='o', label='B')
    ax0.errorbar(xM2,yM2, yerr=y_errorM2, color='purple', fmt='o', label='UVM2')
    ax0.errorbar(xU,yU, yerr=y_errorU, color='green', fmt='o', label='U')
    ax0.errorbar(xV,yV, yerr=y_errorV, color='red', fmt='o', label='V')
    ax0.errorbar(xW1,yW1, yerr=y_errorW1, color='cyan', fmt='o', label='UVW1')
    ax0.errorbar(xW2,yW2, yerr=y_errorW2, color='orange', fmt='o', label='UVW2')
    ax0.errorbar(smartsB_x,smartsB_y, yerr=smartsB_error, color='black', fmt='o', label='Smarts B')
    ax0.errorbar(smartsV_x,smartsV_y, yerr=smartsV_error, color='magenta', fmt='o', label='Smarts V')
    ax0.set_title(object)
    ax0.legend(bbox_to_anchor=(0.98, 1.02), loc=2)
    ax0.set_ylabel('Magnitude')
    xticklabels = ax0.get_xticklabels()
    plt.setp(xticklabels, visible=True)
    fmt=matplotlib.ticker.ScalarFormatter(useOffset=True)
    fmt.set_scientific(False)
    gca().xaxis.set_major_formatter(fmt)
    plt.subplots_adjust(bottom = 0.1, right = 0.8, top = 0.9)
    plt.show()

def plotBvJ():
    
    # import fits to read fits files from astropy.io
    from astropy.io import fits
    # import pyplot from matplotlib to plot data
    import matplotlib.pyplot as plt
    # import ScalarFormatter from matplotlib.ticker to correct the
    # x-axis ticker values
    from matplotlib.ticker import ScalarFormatter
    
    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for b filter
    hdulistB = fits.open("%s/bb_photOut.fits" % (directory_out))
    tbdataB =  hdulistB[1].data # set variable to table of data in hdulist
    yB = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataB.field('AB_MAG'):
        yB.append(magnitude) # append magnitude to y list
    xB = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataB.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xB.append(time_MJD) # append time to x list
    y_errorB = [] # empty list for y error
    for mag_err in tbdataB.field('AB_MAG_ERR'):
        y_errorB.append(mag_err) # append y_error to y_error list
    hdulistB.close() # close the fits file
    
    # get J data from SMARTS
    import urllib2 # import module for fetching URLs
    # set the correct url
    url = 'http://www.astro.yale.edu/smarts/glast/tables/%s.tab' % (object)
    # create a data file to right the tables from url to
    open('data.txt', 'wb').write(urllib2.urlopen(url).read())
    f = open('data.txt', 'r')
    # create empty lists
    smartsJ_y = []
    smartsJ_x = []
    smartsJ_error = []
    j = 0 # set j to 0
    for line in f:
        # divide data table (one long string) into managable lines/lists
        line = line.strip()
        columns = line.split()
        if j > 2:
            # get J filter data from the table and convert to float
            smarts_J_Start = float(columns[10])
            # convert time from jd to mjd
            smarts_J_mjd = smarts_J_Start - 2400000.5
            # convert magnitudes from Johnsons to AB
            smarts_J_mag = float(columns[11])
            smarts_J_ABmag = smarts_J_mag
            smarts_J_err = float(columns[12])
            # append B filter data in date range to above list
            if smarts_J_mjd >= xB[0] and smarts_J_mjd <= xB[len(xB) - 1]:
                smartsJ_x.append(smarts_J_mjd)
                smartsJ_y.append(smarts_J_ABmag)
                smartsJ_error.append(smarts_J_err)
        j = j + 1
    f.close() # close file
    
    print object
    b1 = 0
    for b in yB:
        print 'B = ', yB[b1]
        print 'time = ', xB[b1]
        b1 = b1 + 1
    j1 = 0
    for w in smartsJ_y:
        print 'J = ', smartsJ_y[j1]
        print 'time = ', smartsJ_x[j1]
        j1 = j1 + 1

    s = []
    h = []
    for num in range(2):
        s.append(num + 15.4)
        h.append(s[num] - 2.7)

    
    # make lists of what will be plot
    plotB=[]
    plotJ=[]
    plotBerror=[]
    plotJerror=[]
    for x in range(len(xB)):
        for y in range(len(smartsJ_x)):
            if int(xB[x]) == int(smartsJ_x[y]):
                plotB.append(yB[x] + 0.163)
                plotJ.append(smartsJ_y[y])
                plotBerror.append(y_errorB[x])
                plotJerror.append(smartsJ_error[y])

    b = 0
    for x in plotB:
        print 'B = ', plotB[b]
        b = b + 1
    j = 0
    for y in plotJ:
        print 'J = ', plotJ[j]
        j = j + 1


    # plot B magnitude vs J magnitude
    fig, ax = plt.subplots(nrows = 1)
    ax.errorbar(plotB, plotJ, xerr=plotBerror, yerr=plotJerror, color='blue', fmt='o', label='J v B')
    ax.set_ylabel('J Magnitude')
    ax.set_xlabel('B Magnitude')
    ax.set_title(object)
    xticklabels = ax.get_xticklabels()
    plt.setp(xticklabels, visible=True)
    fmt=matplotlib.ticker.ScalarFormatter(useOffset=True)
    fmt.set_scientific(False)
    gca().xaxis.set_major_formatter(fmt)
    yticklabels = ax.get_yticklabels()
    plt.setp(yticklabels, visible=True)
    fmt=matplotlib.ticker.ScalarFormatter(useOffset=True)
    fmt.set_scientific(False)
    gca().yaxis.set_major_formatter(fmt)
    fit_np = np.polyfit(plotB, plotJ, 1)
    ax.plot(plotB, np.polyval(fit_np, plotB), "r-", lw = 2, label='best fit')
    ax.plot( s, h, "g-", lw = 2, label='one to one line')
    legend(loc=2)
    show()


def plotM2vsJ():
    # import fits to read fits files from astropy.io
    from astropy.io import fits
    # import pyplot from matplotlib to plot data
    import matplotlib.pyplot as plt
    # import ScalarFormatter from matplotlib.ticker to correct the
    # x-axis ticker values
    from matplotlib.ticker import ScalarFormatter
    
    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for m2 filter
    hdulistM2 = fits.open("%s/m2_photOut.fits" % (directory_out))
    tbdataM2 =  hdulistM2[1].data # set variable to table of data in hdulist
    yM2 = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataM2.field('AB_MAG'):
        yM2.append(magnitude) # append magnitude to y list
    xM2 = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataM2.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xM2.append(time_MJD) # append time to x list
    y_errorM2 = [] # empty list for y error
    for mag_err in tbdataM2.field('AB_MAG_ERR'):
        y_errorM2.append(mag_err) # append y_error to y_error list
    hdulistM2.close() # close the fits file

    # get J data from SMARTS
    import urllib2 # import module for fetching URLs
    # set the correct url
    url = 'http://www.astro.yale.edu/smarts/glast/tables/%s.tab' % (object)
    # create a data file to right the tables from url to
    open('data.txt', 'wb').write(urllib2.urlopen(url).read())
    f = open('data.txt', 'r')
    # create empty lists
    smartsJ_y = []
    smartsJ_x = []
    smartsJ_error = []
    j = 0 # set j to 0
    for line in f:
        # divide data table (one long string) into managable lines/lists
        line = line.strip()
        columns = line.split()
        if j > 2:
            # get J filter data from the table and convert to float
            smarts_J_Start = float(columns[10])
            # convert time from jd to mjd
            smarts_J_mjd = smarts_J_Start - 2400000.5
            # convert magnitudes from Johnsons to AB
            smarts_J_mag = float(columns[11])
            smarts_J_ABmag = smarts_J_mag - 0.163
            smarts_J_err = float(columns[12])
            # append B filter data in date range to above list
            if smarts_J_mjd >= xM2[0] and smarts_J_mjd <= xM2[len(xM2) - 1]:
                smartsJ_x.append(smarts_J_mjd)
                smartsJ_y.append(smarts_J_ABmag)
                smartsJ_error.append(smarts_J_err)
        j = j + 1
    f.close() # close file
    
    # make lists of what will be plot
    plotM2=[]
    plotJ=[]
    plotM2error=[]
    plotJerror=[]
    for x in range(len(xM2)):
        for y in range(len(smartsJ_x)):
            if int(xM2[x]) == int(smartsJ_x[y]):
                plotM2.append(yM2[x])
                plotJ.append(smartsJ_y[y])
                plotM2error.append(y_errorM2[x])
                plotJerror.append(smartsJ_error[y])
# m2j = 0
        #   for x in xM2:
        #print 'time_m2 = ', xM2[m2j], ' time_J = ', smartsJ_x[m2j]
        #m2j = 1 + m2j

    s = []
    h = []
    for num in range(2):
        s.append(num + 16.4)
        h.append(s[num] - 3.9)

    # plot J magnitude vs M2 magnitude
    fig, ax = plt.subplots(nrows = 1)
    ax.errorbar(plotM2, plotJ, xerr=plotM2error, yerr=plotJerror, color='blue', fmt='o', label='J v M2')
    ax.set_ylabel('J Magnitude')
    ax.set_xlabel('M2 Magnitude')
    ax.set_title(object)
    xticklabels = ax.get_xticklabels()
    plt.setp(xticklabels, visible=True)
    fmt=matplotlib.ticker.ScalarFormatter(useOffset=True)
    fmt.set_scientific(False)
    gca().xaxis.set_major_formatter(fmt)
    yticklabels = ax.get_yticklabels()
    plt.setp(yticklabels, visible=True)
    fmt=matplotlib.ticker.ScalarFormatter(useOffset=True)
    fmt.set_scientific(False)
    gca().yaxis.set_major_formatter(fmt)
    fit_np = np.polyfit(plotM2, plotJ, 1)
    ax.plot(plotM2, np.polyval(fit_np, plotM2), "r-", lw = 2, label='best fit')
    ax.plot( s, h, "g-", lw = 2, label='one to one line')
    legend(loc=2)
    show()

def plotBvFermi():

    # import fits to read fits files from astropy.io
    from astropy.io import fits
    # import pyplot from matplotlib to plot data
    import matplotlib.pyplot as plt
    # import ScalarFormatter from matplotlib.ticker to correct the
    # x-axis ticker values
    from matplotlib.ticker import ScalarFormatter

    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for b filter
    hdulistB = fits.open("%s/bb_photOut.fits" % (directory_out))
    tbdataB =  hdulistB[1].data # set variable to table of data in hdulist
    yB = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataB.field('AB_MAG'):
        yB.append(magnitude) # append magnitude to y list
    xB = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataB.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xB.append(time_MJD) # append time to x list
    y_errorB = [] # empty list for y error
    for mag_err in tbdataB.field('AB_MAG_ERR'):
        y_errorB.append(mag_err) # append y_error to y_error list
    hdulistB.close() # close the fits file

    # extract Fermi data from previously downloaded files
    # open the fits file previously downloaded from the Fermi website

    # set the file name depending on the object
    if object == '3C279':
        filename = 'J123939+044409_604800.lc.txt'
    elif object == '3C454':
        filename = '3C454.3_86400.lc.txt'
    elif object == '2155-304':
        filename = 'PKS2155-304_86400.lc.txt'

    # open first file and add energy and time data to list
    hdulistLAT = fits.open("/Users/Stephanie/Fermi/%s/%s" % (object, filename))
    tbdataLAT = hdulistLAT[1].data # set variable to table of data in hdulist
    xLAT = [] # empty list for the x values (time)
    yLAT = [] # empty list for y values (magnitudes)
    for time_seconds_MET in tbdataLAT.field('START'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xLAT.append(time_MJD)
    for flux in tbdataLAT.field('FLUX_100_300000'):
        yLAT.append(flux)
    # get the errors
    errorLAT = []
    for error  in tbdataLAT.field('ERROR_100_300000'):
        errorLAT.append(error)

    # Plot only the same portion as Swift and SMARTS data
    lat = 0 # iteration variable
    plotLATx = []
    plotLATy = []
    plotLATerror = []
    for times in xLAT:
        if times >= xB[0] and times <= xB[len(xB) - 1]:
            plotLATx.append(times)
            plotLATy.append(yLAT[lat])
            plotLATerror.append(errorLAT[lat])
        lat = lat + 1
    hdulistLAT.close()

    # make lists of what will be plot
    plotB=[]
    plotLAT=[]
    plotBerror=[]
    plotLATerror1=[]
    for x in range(len(xB)):
        for y in range(len(plotLATx)):
            if int(xB[x]) == int(plotLATx[y]):
                plotB.append(math.pow(10,yB[x]))
                plotLAT.append(plotLATy[y])
                plotBerror.append(math.pow(10,y_errorB[x]))
                plotLATerror1.append(plotLATerror[y])

#blat = 0
#   for x in xB:
#       print 'time_LAT = ', plotLATx[blat], ' time_B = ', xB[blat]
#       blat = 1 + blat
    s = []
    h = []
    for num in range(3):
        s.append(num + 13)
        h.append(num + 13)

    # plot J magnitude vs M2 magnitude
    fig, ax = plt.subplots(nrows = 1)
    ax.errorbar(plotB, plotLAT, xerr=plotBerror, yerr=plotLATerror1, color='blue', fmt='o', label='.1-300GeV vs B')
    ax.set_ylabel('.1-300GeV ' + r'$photons/cm^.s$' )
    ax.set_xlabel('B Magnitude')
    ax.set_title(object)
    xticklabels = ax.get_xticklabels()
    plt.setp(xticklabels, visible=True)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#fit_np = np.polyfit(plotB, plotLAT, 1)
#   ax.plot(plotB, np.polyval(fit_np, plotB), "r-", lw = 2, label='best fit')
#   ax.plot( s, h, "g-", lw = 2, label='one to one line')
#   legend(loc=2)
    show()

def plotM2vFermi():
    
    # import fits to read fits files from astropy.io
    from astropy.io import fits
    # import pyplot from matplotlib to plot data
    import matplotlib.pyplot as plt
    # import ScalarFormatter from matplotlib.ticker to correct the
    # x-axis ticker values
    from matplotlib.ticker import ScalarFormatter
    import numpy as np
    
    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for b filter
    hdulistM2 = fits.open("%s/m2_photOut.fits" % (directory_out))
    tbdataM2 =  hdulistM2[1].data # set variable to table of data in hdulist
    yM2 = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataM2.field('AB_MAG'):
        yM2.append(magnitude) # append magnitude to y list
    xM2 = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataM2.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xM2.append(time_MJD) # append time to x list
    y_errorM2 = [] # empty list for y error
    for mag_err in tbdataM2.field('AB_MAG_ERR'):
        y_errorM2.append(mag_err) # append y_error to y_error list
    hdulistM2.close() # close the fits file
    
    # extract Fermi data from previously downloaded files
    # open the fits file previously downloaded from the Fermi website
    
    # set the file name depending on the object
    if object == '3C279':
        filename = 'J123939+044409_604800.lc.txt'
    elif object == '3C454':
        filename = '3C454.3_86400.lc.txt'
    elif object == '2155-304':
        filename = 'PKS2155-304_86400.lc.txt'
    
    # open first file and add energy and time data to list
    hdulistLAT = fits.open("/Users/Stephanie/Fermi/%s/%s" % (object, filename))
    tbdataLAT = hdulistLAT[1].data # set variable to table of data in hdulist
    xLAT = [] # empty list for the x values (time)
    yLAT = [] # empty list for y values (magnitudes)
    for time_seconds_MET in tbdataLAT.field('START'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xLAT.append(time_MJD)
    for flux in tbdataLAT.field('FLUX_100_300000'):
        yLAT.append(flux)
    # get the errors
    errorLAT = []
    for error  in tbdataLAT.field('ERROR_100_300000'):
        errorLAT.append(error)
    
    # Plot only the same portion as Swift and SMARTS data
    lat = 0 # iteration variable
    plotLATx = []
    plotLATy = []
    plotLATerror = []
    for times in xLAT:
        if times >= xM2[0] and times <= xM2[len(xM2) - 1]:
            plotLATx.append(times)
            plotLATy.append(yLAT[lat])
            plotLATerror.append(errorLAT[lat])
        lat = lat + 1
    hdulistLAT.close()
    
    # make lists of what will be plot
    plotM2=[]
    plotLAT=[]
    plotM2error=[]
    plotLATerror1=[]
    for x in range(len(xM2)):
        for y in range(len(plotLATx)):
            if int(xM2[x]) == int(plotLATx[y]):
                plotM2.append(math.pow(10,yM2[x]))
                plotLAT.append(plotLATy[y])
                plotM2error.append(math.pow(10,y_errorM2[x]))
                plotLATerror1.append(plotLATerror[y])
    
#m2lat = 0
#   print len(xM2)
#   for x in xM2:
#       print 'time_LAT = ', plotLATx[m2lat], ' time_m2 = ', xM2[m2lat]
#       m2lat = 1 + m2lat
#       print m2lat
    s = []
    h = []
    for num in range(3):
        s.append(num + 13)
        h.append(num + 13)

    # plot M2 magnitude vs Fermi magnitude
    fig, ax = plt.subplots(nrows = 1)
    ax.errorbar(plotM2, plotLAT, xerr=plotM2error, yerr=plotLATerror1, color='blue', fmt='o', label='.1-300 GeV vs M2')
    ax.set_ylabel('.1-300GeV' + ' ' + r'$photons/cm^2/s$')
    ax.set_xlabel('UVM2 Magnitude')
    ax.set_title(object)
    xticklabels = ax.get_xticklabels()
    plt.setp(xticklabels, visible=True)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#   fit_np = np.polyfit(plotM2, plotLAT, 1)
#   ax.plot(plotM2, np.polyval(fit_np, plotM2), "r-", lw = 2, label='best fit')
#   ax.plot( s, h, "g-", lw = 2, label='one to one line')
#   legend(loc=2)
    plt.show()


def plotSMARTSBJ():
    # import fits to read fits files from astropy.io
    from astropy.io import fits
    # import pyplot from matplotlib to plot data
    import matplotlib.pyplot as plt
    # import ScalarFormatter from matplotlib.ticker to correct the
    # x-axis ticker values
    from matplotlib.ticker import ScalarFormatter

    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for m2 filter
    hdulistM2 = fits.open("%s/m2_photOut.fits" % (directory_out))
    tbdataM2 =  hdulistM2[1].data # set variable to table of data in hdulist
    yM2 = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataM2.field('AB_MAG'):
        yM2.append(magnitude) # append magnitude to y list
    xM2 = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataM2.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
    # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xM2.append(time_MJD) # append time to x list
    y_errorM2 = [] # empty list for y error
    for mag_err in tbdataM2.field('AB_MAG_ERR'):
        y_errorM2.append(mag_err) # append y_error to y_error list
    hdulistM2.close() # close the fits file

    # Get B and V Data from SMARTS
    import urllib2 # import module for fetching URLs
    # set the correct url
    url = 'http://www.astro.yale.edu/smarts/glast/tables/%s.tab' % (object)
    # create a data file to right the tables from url to
    open('data.txt', 'wb').write(urllib2.urlopen(url).read())
    f = open('data.txt', 'r')
    # create empty lists
    smartsB_y = []
    smartsB_x = []
    smartsB_error = []
    b = 0 # set b to 0
    smartsJ_y = []
    smartsJ_x = []
    smartsJ_error = []
    j = 0 # set v to 0
    for line in f:
        # divide data table (one long string) into managable lines/lists
        line = line.strip()
        columns = line.split()
        if b > 2:
            # get B filter data from the table and convert to float
            smarts_B_Start = float(columns[1])
            # convert time from jd to mjd
            smarts_B_mjd = smarts_B_Start - 2400000.5
            # convert magnitudes from Johnsons to AB
            smarts_B_mag = float(columns[2])
            smarts_B_ABmag = smarts_B_mag - 0.163
            smarts_B_err = float(columns[3])
            # append B filter data in date range to above list
            if smarts_B_mjd >= xM2[0] and smarts_B_mjd <= xM2[len(xM2) - 1]:
                smartsB_x.append(smarts_B_mjd)
                smartsB_y.append(smarts_B_ABmag)
                smartsB_error.append(smarts_B_err)
        b = b + 1
        if j > 2:
            # get V filter data from the table and convert to float
            smarts_J_Start = float(columns[10])
            # convert time from jd to mjd
            smarts_J_mjd = smarts_J_Start - 2400000.5
            smarts_J_mag = float(columns[11])
            # convert from Johnson to AB magnitude
            smarts_J_ABmag = smarts_J_mag - 0.044
            smarts_J_err = float(columns[12])
            # append V filter data in date range to above list
            if smarts_J_mjd >= xM2[0] and smarts_J_mjd <= xM2[len(xM2) - 1]:
                smartsJ_x.append(smarts_J_mjd)
                smartsJ_y.append(smarts_J_ABmag)
                smartsJ_error.append(smarts_J_err)
        j = j + 1
    
    f.close() # close file


    # make lists of what will be plot
    plotB=[]
    plotJ=[]
    plotBerror=[]
    plotJerror=[]
    for x in range(len(smartsB_x)):
        for y in range(len(smartsJ_x)):
            if int(smartsB_x[x]) == int(smartsJ_x[y]):
                plotB.append(smartsB_y[x])
                plotJ.append(smartsJ_y[y])
                plotBerror.append(smartsB_error[x])
                plotJerror.append(smartsJ_error[y])
#  m2j = 0
#   for x in xM2:
#print 'time_m2 = ', xM2[m2j], ' time_J = ', smartsJ_x[m2j]
#m2j = 1 + m2j

    # open the fits file that contains all of the magnitudes ouput
    # using uvotsource on each observation for b filter
    hdulistB = fits.open("%s/bb_photOut.fits" % (directory_out))
    tbdataB =  hdulistB[1].data # set variable to table of data in hdulist
    yB = [] # empty list for y values in plot
    # get magnitude values from table in fits file
    for magnitude in tbdataB.field('AB_MAG'):
        yB.append(magnitude) # append magnitude to y list
        xB = [] # empty list for x values
    # get the time that the exposure was taken from the file
    for time_seconds_MET in tbdataB.field('MET'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # NOTE: Swift MET starts from Jan. 1st, 2001 (MMDDYY) = 51910 (MJD)
        time_MJD = 51910 + time_days_MET # conver to MJD
        xB.append(time_MJD) # append time to x list
    y_errorB = [] # empty list for y error
    for mag_err in tbdataB.field('AB_MAG_ERR'):
        y_errorB.append(mag_err) # append y_error to y_error list
    hdulistB.close() # close the fits file


    # make lists of what will be plot
    plotB2=[]
    plotJ2=[]
    plotB2error=[]
    plotJ2error=[]
    plotBtime = []
    plotJtime = []
    for x in range(len(xB)):
        for y in range(len(smartsJ_x)):
            if int(xB[x]) == int(smartsJ_x[y]):
                plotB2.append(yB[x] + 0.163)
                plotJ2.append(smartsJ_y[y])
                plotB2error.append(y_errorB[x])
                plotJ2error.append(smartsJ_error[y])
                plotBtime.append(xB[x])
                plotJtime.append(smartsJ_x[y])

    m2j = 0
    for x in plotBtime:
        print 'time_B = ', plotBtime[m2j], ' time_J = ', plotJtime[m2j]
        print 'B = ', plotB2[m2j] , ' J = ' , plotJ2[m2j]
        m2j = m2j + 1
    l = 0
    for y in plotJtime:
        print 'time_Bsmarts = ', smartsB_x[l], ' time_J = ', smartsJ_x[l]
        l = 1 + l

    s = []
    h = []
    for num in range(3):
        s.append(num + 13.6)
        h.append(s[num] - 2.5 )


    # plot J magnitude vs M2 magnitude
    fig, ax = plt.subplots(nrows = 1)
    ax.errorbar(plotB, plotJ, xerr=plotBerror, yerr=plotJerror, color='blue', fmt='o', label='SMARTs JvB')
    ax.errorbar(plotB2, plotJ2, xerr=plotB2error, yerr=plotJ2error, color='red', fmt='o', label=' J v Swift B')
    ax.set_ylabel('J Magnitude')
    ax.set_xlabel('B Magnitude')
    ax.set_title(object)
    xticklabels = ax.get_xticklabels()
    plt.setp(xticklabels, visible=True)
    fmt=matplotlib.ticker.ScalarFormatter(useOffset=True)
    fmt.set_scientific(False)
    gca().xaxis.set_major_formatter(fmt)
    yticklabels = ax.get_yticklabels()
    plt.setp(yticklabels, visible=True)
    fmt=matplotlib.ticker.ScalarFormatter(useOffset=True)
    fmt.set_scientific(False)
    gca().yaxis.set_major_formatter(fmt)
    fit_np = np.polyfit(plotB, plotJ, 1)
    ax.plot(plotB, np.polyval(fit_np, plotB), "r-", lw = 2, label='best fit')
    ax.plot( s, h, "g-", lw = 2, label='one to one line')
    legend(loc=2)
    show()




###################
## Call Functions
#directoryOut()
#unzipFun()
#reduction()

###################
### to run plotSingle(filter) ask user which filter to plot
#filter = raw_input("Which filter should be plot? ")
#plotSingle(filter)
#plotAllOptUVFilters()
#plotAllFilters()

###################
#plotBvV()
#plotBvM2()
#plotM2vsJ()
#plotBvJ()
#plotBvFermi()
plotM2vFermi()
#plotSMARTSBJ()

