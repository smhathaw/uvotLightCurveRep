#!/usr/bin/python

# Reduce Swift data to finished light curves. This script follows my guide for
# uvot light curve reduction which follows the information on Swift's wedbsite

# User inputs the object and observations that will be reduced
object = raw_input("Input object name: ") 
# set observation name variables based on the object input
if object == '3C279':
    observation = '35019'
    initial = '130'
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
    final = '151'
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
    
    # extract Fermi data from previously downloaded files
    # open the fits file previously downloaded from the Fermi website
    
    # set the file name depending on the object
    if object == '3C279':
        filename1 = 'L140725162056B40DE2BC75_PH00.fits'
        filename2 = 'L140725162056B40DE2BC75_PH01.fits'
    elif object == '3C454':
        filename1 = 'L140725164157B40DE2BC86_PH00.fits'
        filename2 = 'L140725164157B40DE2BC86_PH01.fits'
    elif object == '2155-304':
        filename1 = 'L140725160209B40DE2BC27_PH00.fits'
        filename2 = 'L140725160209B40DE2BC27_PH01.fits'
        
    # open first file and add energy and time data to list
    hdulistLAT1 = fits.open("/Users/Stephanie/Fermi/%s/%s" % (object, filename1))
    tbdataLAT1 = hdulistLAT1[1].data # set variable to table of data in hdulist
    yLAT = [] # empty list for y values (magnitudes)
    for energy in tbdataLAT1.field('ENERGY'):
        yLAT.append(energy)
    xLAT = [] # empty list for the x values (time)
    for time_seconds_MET in tbdataLAT1.field('TIME'):
        time_days_MET = time_seconds_MET/86400 # convert seconds to days
        # MET begins at midnight Jan 1, 2014
        time_MJD = 51910 + time_days_MET # convert to MJD
        xLAT.append(time_MJD)
    hdulistLAT1.close()
    
    # open second file and energya and time data to list
    hdulistLAT2 = fits.open("/Users/Stephanie/Fermi/%s/%s" % (object, filename2))
    tbdataLAT2 = hdulistLAT2[1].data # set variable to table of data in hdulist
    for energy2 in tbdataLAT2.field('ENERGY'):
        yLAT.append(energy2)
    for time_seconds_MET2 in tbdataLAT2.field('TIME'):
        time_days_MET2 = time_seconds_MET2/86400 # convert seconds to days
        # MET begins at midnisht Jan 1, 2014
        time_MJD2 = 51910 + time_days_MET2 # convert to MJD
        xLAT.append(time_MJD2)
    hdulistLAT2.close()
    
    
    # Plot magnitude vs time (x vs y) in a single scatter plot if prompted
    fig, (ax0, ax1) = plt.subplots(nrows=2)
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
    xticklabels = ax0.get_xticklabels()
    plt.setp(xticklabels, visible=True)
    fmt=matplotlib.ticker.ScalarFormatter(useOffset=True)
    fmt.set_scientific(False)
    gca().xaxis.set_major_formatter(fmt)
    
    ax1.scatter(xLAT, yLAT, color = 'gray')
    plt.ylabel('MeV')
    plt.xlabel('MJD')
    plt.subplots_adjust(bottom = 0.1, right = 0.8, top = 0.9)
#   xticklabels = ax0.get_xticklabels()
#   plt.setp(xticklabels, visible=True)
#   fmt=matplotlib.ticker.ScalarFormatter(useOffset=False)
#   fmt.set_scientific(False)
#   gca().xaxis.set_major_formatter(fmt)
#   plt.show()

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


    # make losts of what will be plot
    plotB=[]
    plotV=[]
    plotBerror=[]
    plotVerror=[]
    minVB = min(len(xB), len(xV))
    for x in range(minVB):
        if int(xB[x]) == int(xV[x]):
            plotB.append(yB[x])
            plotV.append(yV[x])
            plotBerror.append(y_errorB[x])
            plotVerror.append(y_errorV[x])

    # plot B magnitude vs V magnitude
    fig, ax = plt.subplots(nrows = 1)
    ax.errorbar(plotB, plotV, xerr=plotBerror, yerr=plotVerror, color='blue', fmt='o')
    ax.set_ylabel('B Magnitude')
    ax.set_xlabel('V Magnitude')
    ax.set_title(object)
    xticklabels = ax.get_xticklabels()
    plt.setp(xticklabels, visible=True)
    fmt=matplotlib.ticker.ScalarFormatter(useOffset=True)
    fmt.set_scientific(False)
    gca().xaxis.set_major_formatter(fmt)
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
    
    # make losts of what will be plot
    plotB=[]
    plotM2=[]
    minM2B = min(len(xB), len(xM2))
    for x in range(minM2B):
        if int(xB[x]) == int(xM2[x]):
            plotB.append(yB[x])
            plotM2.append(yM2[x])
    # plot B magnitude vs V magnitude
    fig, ax = plt.subplots(nrows = 1)
    ax.plot(plotM2,plotB, 'bo' )
    ax.set_ylabel('B Magnitude')
    ax.set_xlabel('M2 Magnitude')
    ax.set_title(object)
    xticklabels = ax.get_xticklabels()
    plt.setp(xticklabels, visible=True)
    fmt=matplotlib.ticker.ScalarFormatter(useOffset=True)
    fmt.set_scientific(False)
    gca().xaxis.set_major_formatter(fmt)
    show()


###################
## Call Functions
#directoryOut()
#unzipFun()
#reduction()
### to run plotSingle(filter) ask user which filter to plot
#filter = raw_input("Which filter should be plot? ")
#plotSingle(filter)
#plotAllFilters()
plotBvV()
#plotBvM2()
