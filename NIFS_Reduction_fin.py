import pyfits as fits
import numpy as np
from scipy.optimize import curve_fit, leastsq
import sys
import getopt
import os
import time
import shutil


'''
Last changed: 
July 29. 2014 - Morten Stostad - Improved readability

Gemini NIFS data reduction script - PYRAF version
Reduction for:  GENERAL BASELINE CALIBRATIONS

DESCRIPTION   	         	                                  
                                                                         
This script is a pyraf-compatible version of the original NIFS scripts  
written by Tracy Beck. Original pyraf-compatible code by R. McDermid, with
further modifications by T. Do and M. Stostad. 

This is a module written to reduce all NIFS-data, put together by combining
three different modules by R. McDermid and T. Beck (Step 2, 3 and 4) with
improved telluric correction (Step 3.2). The following is the minimum needed
to be able to reduce data:

Science files taken on the night of observation
Calibration files taken on the night of observation
Text files of arc lines (or telluric lines if lamp calib files unavailable)
Template of a typical A-star (e.g. Vega)
Lists of all the calibration files (e.g. skylist, arclist, etc)
Pyraf, numpy, pyfits and scipy installed

There are four major steps in the reduction procedure.

1. Creating separate files that contain lists of the files; there 
   should be one skylist that has a list of all the sky files, one
   arclist that has all the arc files, etc. THIS STEP IS NOT IN THE PIPELINE.
   Use mkNifsFilesList.py to create most of them - some have to be created
   manually, however..
2. Preparing the base calibration files (ronchi flats, bad pixel mask,  
   wavelength solutions, etc). Yields five files; a shift reference
   file, a flat field, a flat BPM, an arc frame and a ronchi flat.
3. Preparing for telluric correction. Yields a telluric correction fits file.
4. Reducing the science data by using the files created in 1) and 2). 

The following are the only parameters the user should have to change.
Of course, it doesn't hurt to understand the rest of the code, but 
ideally the program should work fine by just changing these parameters
and running run_reduction() in the python shell.

DISCLAIMER: This is a very basic Python pipeline, and it is not 
particularly well made. We publish it mainly because the previous
pipeline from the Gemini website had obvious flaws, and could not
be used in Python without a lot of extra work. Anyone that uses
this pipeline should be aware that although it should work for
general reduction of NIFS data, its main purpose was reduction of
the data from Stostad+ 2014, and has only been tested on that data.
'''

#Change a value to False if the step should be skipped. The steps are
#ordered [Prepare base calibration files, prepare telluric file, 
#reduce science data]. For instance, if the telluric correction file 
#is already acquired, change to [True, False, True].

steps = [True, True, True]

# The folders for the calibration and science files, and the folder for the reduced
# files. (raw_data and red_data are just the same variables in step 3&4).
# dat_dir = name of top folder
# raw_dir = location of calibration files
# reduce_dir = backup location
# reduced_all = folder for reduced files

datDir = "/Users/username/nifs/date/"
raw_dir = datDir+"calib/"
reduce_dir= datDir+"backup/"
reduced_all = datDir + 'reduced/'
raw_data = raw_dir
red_data   = reduce_dir


# No parameters for step 1; the user has to create the filelists themselves.
# If the script mkNifsFilesList.py is available, the user could use that.

#################### PARAMETERS FOR STEP 2 ####################


# Folders of the raw calibration data, the output folder for reduced files and
# the folder where the lamp reference files are. These are the only parameters
# every user will have to change for step 1. The backup directory for the 
# finished calibration files can also be changed.

arcDir = "/Users/username/nifs/arclamps/data/"  # lamp reference files folder
backupDir = reduce_dir+'calib_backup/' # another backup folder


useSkyLines = False      # set to True if no lamp calibrations were taken and
                         # telluric lines should be used for wavelength calibration.
                         # In almost all cases this should be False.

# It's assumed that the list of files (e.g. flatlist) are already in 
# datDir/calib/ from Step 1.

# Set the file names for the lists of the calibration file names. 
# These are the names of the list created in Step 1.

flatlist = "flatlist"
flatdarklist = "flatdarklist"
ronchilist = "ronchilist"

if useSkyLines:
    arcStr = 'skylist'
    arcStrDrk = 'skydarklist'
else:
    arcStr = 'arclist'
    arcStrDrk = 'arcdarklist'

# Change steps_basec if not all the substeps of the base calibration should be run
# [find first shift file, make flats, find wavelength solution, make ronchi flat]
steps_basec = [True,True,True,True]


#################### PARAMETERS FOR STEP 3 ####################

# The spectral resolution of the science data spectra.

R = 5000

# These are the file names of processed calibration files, which were created by
# step 2. If they're not manually modified, simply put ManualNames = False,
# and you won't have to worry about them. If you have changed the names since
# step 2 created them, put ManualNames = True and change the file names here.

ManualNames = False
# calflat    = "flat"
# bpm = "flat_bpm.pl"
# arc        = "wsoln"
# ronchiflat = "rflat"
# shiftimage = "shiftFile"

# Specify the name of the text file containing the name of the A star
# fits file and the corresponding sky file.
# 
# Example: 
# tellistfile = 'astar_hd155379'
# astarskylistfile = 'astar_hd155379_sky'

tellistfile = ""
astarskylistfile = ""

# The high-resolution stellar template:
fits_template = '/Users/username/nifs/templates/vega_all.fits'

# The output name of the final telluric file:
telluric_norm = 'telluric_norm_weighted.fits'



#################### PARAMETERS FOR STEP 4 ####################


rmReducedData = False

# Same as above. Don't worry about the names (they will be
# overwritten) if ManualNames = False.
ManualNames = False

calflat = "flat"
arc = "wrgnN20120505S0384"
ronchiflat = "rgnN20120505S0542"
shiftimage = "shiftFile"
bpm = "flat_bpm.pl"


# Name of the text files with the list of science/sky filenames.
scilistfile = 'sciencename'
skylistfile = 'skylist'

# optionally set to use a specific sky file (if not set to None).
# otherwise, will combine all the files in the skyfile list into 1 and
# use that sky. This is in contrast to the procedure in the original
# script that requires one sky for every science frame.
useSkyFile = None

telluric  = telluric_norm


# END EDITING HERE
#--------------------------------------------------------------------------

#STEP 1:
#The user must create the lists of calibration files outside this script.


######################################################################
######################################################################
######################################################################
# STEP 2:
# PREPARING BASE CALIBRATION FILES
#
# To do the main part of the calibration in step 4 we must create certain
# calibration files. There are five of these:
#
# flat.fits
# flat_bpm.pl
# wsoln.fits
# rflat.fits
# shiftFile.fits
#
# The following code was written by R. McDermid. It should create all
# these files in the current directory.
######################################################################
# Gemini NIFS data reduction script - PYRAF version
# Reduction for:  GENERAL BASELINE CALIBRATIONS
#
###########################################################################
# DESCRIPTION   	         	                                  #
#                                                                         #
# This script is a pyraf-compatible version of the original NIFS scripts  #
# written by Tracy Beck.                                                  #
#                                                                         #
#                                                                         #
###########################################################################
# Current limitations:  The NIFS Baseline calibration reductions have     #
# not been well tested on data that was obtained with non-standard        #
# wavelength configurations.  (i.e., a different central wavelength       #
# setting than the default Z, J, H and K-band values of 1.05, 1.25, 1.65  #
# and 2.20 microns).                                                      #
#                                                                         #
###########################################################################

###########################################################################
# STEP 2.1:  Prepare IRAF  		                                  #
###########################################################################

# Import the pyraf module and relevant packages
from pyraf import iraf


def prep_nifs(log):
    iraf.gemini()
    iraf.nifs()
    iraf.gnirs()
    iraf.gemtools()

    # Unlearn the used tasks
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs,iraf.nifs)

    iraf.set(stdimage='imt2048')

    # Prepare the package for NIFS
    iraf.nsheaders("nifs",logfile=log)

    # Set clobber to 'yes' for the script. This still does not make the gemini tasks
    # overwrite files, so you will likely have to remove files if you re-run the script.
    user_clobber=iraf.envget("clobber")
    iraf.reset(clobber='yes')


###########################################################################
# STEP 2.2:  Determine the shift to the MDF file                          #
###########################################################################

def mk_mdf_shift(infile,raw_dir,log,outpref="s",overwrite=True):
    outfile = outpref+infile
    if os.path.exists(outfile) & overwrite:
        print '% mk_mdf_shift: file exists, will now overwrite: '+outfile
        os.remove(outfile)
    # find the MDF shift, by default prefixes the file with an "s"
    print '% NIFS_Basecalib: Determine the shift to the MDF file:'
    iraf.nfprepare(infile,rawpath=raw_dir,outpref=outpref, shiftx='INDEF',  
                   shifty='INDEF',fl_vardq='no',fl_corr='no',fl_nonl='no', logfile=log)
    
    # return the name of the shiftfile
    return outfile

###########################################################################
# STEP 2.3:  Make the Flat field and BPM                                  #
###########################################################################

def mk_flat(flatlist,flatdarklist,raw_dir,shift_image,log):
    print '% NIFS_Basecalib: Make the flat field and BPM:'

    # nfprepare preps the data by validating the data, read the array
    # data, look at saturation and linear limit levels, optionally
    # subtract reference pixels, correct for non-linearity (not on by
    # default), calculate the variance, calculate the data quality, and
    # detect cosmic rays (not on by default).
    #
    # Importantly, nfprepare also adds in the shifts from flexture (MDF
    # shifts) to the data. In this case, it was calculated using one file
    # above and used for the rest here as 'shiftim'.

    # note: the keyword 'fl_cut' in nsreduce was only recently changed
    # in Dec. 2012 from 'fl_nscut', so there may be bugs with this
    # keyword depending on which version of the Gemini IRAF routines
    # you are using
    calflat=str(open(flatlist, "r").readlines()[0]).strip()
    flatdark=str(open(flatdarklist, "r").readlines()[0]).strip()

    iraf.nfprepare("@"+flatlist,rawpath=raw_dir,shiftim=shift_image,
                   fl_vardq='yes',fl_int='yes',fl_corr='no',fl_nonl='no', logfile=log)

    iraf.nfprepare("@"+flatdarklist,rawpath=raw_dir,shiftim=shift_image,
                   fl_vardq='yes',fl_int='yes',fl_corr='no',fl_nonl='no', logfile=log)


    # gemcombine combines multiple frames into one
    iraf.gemcombine("n//@"+flatlist,output="gn"+calflat,fl_dqpr='yes',
                    fl_vardq='yes',masktype="none",logfile=log)
    iraf.gemcombine("n//@"+flatdarklist,output="gn"+flatdark,fl_dqpr='yes',
                    fl_vardq='yes',masktype="none",logfile=log,combine="median")

    # nsreduce does sky subtraction and flattens the spectroscopic data
    iraf.nsreduce("gn"+calflat,fl_cut='yes',fl_nsappw='yes',fl_vardq='yes',
                   fl_sky='no',fl_dark='no',fl_flat='no',logfile=log,outprefix='r')
    iraf.nsreduce("gn"+flatdark,fl_cut='yes',fl_nsappw='yes',fl_vardq='yes',
                  fl_sky='no',fl_dark='no',fl_flat='no',logfile=log,outprefix='r')


    # creating flat image, final name = rnN....._sflat.fits
    print '% NIFS_Basecalib: creating flat image, final name = rgnN....._sflat.fits'

    iraf.nsflat("rgn"+calflat,darks="gn"+os.path.splitext(flatdark)[0],
                flatfile="rgn"+os.path.splitext(calflat)[0]+"_sflat",
                darkfile="rgn"+os.path.splitext(flatdark)[0]+"_dark",
                fl_save_dark='yes',process="fit",
                thr_flo=0.15,thr_fup=1.55,fl_vardq='yes',logfile=log) 


    #rectify the flat for slit function differences - make the final flat.
    print '% NIFS_Basecalib: rectify the flat for slit function differences - make the final flat.'
    iraf.nsslitfunction("rgn"+calflat,"rgn"+os.path.splitext(calflat)[0]+"_flat",
                        flat="rgn"+os.path.splitext(calflat)[0]+"_sflat",
                        dark="rgn"+os.path.splitext(flatdark)[0]+"_dark",
                        combine="median",
                        order=3,fl_vary='no',logfile=log)

    # return [flat file, dark for the flat, bad pixel mask]
    return ("rgn"+os.path.splitext(calflat)[0]+"_flat","rgn"+os.path.splitext(flatdark)[0]+"_dark",
            "rgn"+os.path.splitext(calflat)[0]+"_sflat_bpm.pl")

###########################################################################
# STEP 2.4:  Reduce the Arc and determine the wavelength solution         #
###########################################################################

def mk_wavelength_solution(arcStr,arcStrDrk,raw_dir,shift_image,
                           flat_image,bpm,log,arcDir,useSkyLines=False):
    arc=str(open(arcStr, "r").readlines()[0]).strip()
    arcdark=str(open(arcStrDrk, "r").readlines()[0]).strip()

    iraf.nfprepare('@'+arcStr, rawpath=raw_dir, shiftimage=shift_image,
                   bpm=bpm,fl_vardq='yes',
                   fl_corr='no',fl_nonl='no', logfile=log)

    iraf.nfprepare('@'+arcStrDrk, rawpath=raw_dir, shiftimage=shift_image,
                   bpm=bpm,fl_vardq='yes',
                   fl_corr='no',fl_nonl='no', logfile=log)

    # Determine the number of input arcs and arc darks so that the
    # routine runs automatically for single or multiple files.
    nfiles = len(open(arcStr).readlines())
    if nfiles > 1: 
       iraf.gemcombine("n//@"+arcStr,output="gn"+os.path.splitext(arc)[0],
          fl_dqpr='yes',fl_vardq='yes',masktype="none",logfile=log)
    else:
       iraf.copy("n"+os.path.splitext(arc)[0]+".fits",
                 "gn"+os.path.splitext(arc)[0]+".fits")

    nfiles = len(open(arcStrDrk).readlines())
    if nfiles > 1:
       iraf.gemcombine("n//@"+arcStrDrk,output="gn"+os.path.splitext(arcdark)[0],
          fl_dqpr='yes',fl_vardq='yes',masktype="none",logfile=log)
    else:
       iraf.copy("n"+os.path.splitext(arcdark)[0]+".fits",
                 "gn"+os.path.splitext(arcdark)[0]+".fits")

    iraf.nsreduce("gn"+arc,outpr="r",darki="gn"+os.path.splitext(arcdark)[0],
                  flati=flat_image,
                  fl_vardq='no', fl_cut='yes', fl_nsappw='yes', fl_sky='no',
                  fl_dark='yes',fl_flat='yes', 
                  logfile=log)

    ###########################################################################
    #  DATA REDUCTION HINT -                                                  # 
    # For the nswavelength call, the different wavelength settings            #
    # use different vaues for some of the parameters. For optimal auto        #
    # results, use:                                                           #
    #                                                                         #
    # K-band: thresho=50.0, cradius=8.0   -->  (gives rms of 0.1 to 0.3)      # 
    # H-band: thresho=100.0, cradius=8.0  -->  (gives rms of 0.05 to 0.15)    #
    # J-band: thresho=100.0               -->  (gives rms of 0.03 to 0.09)    # 
    # Z-band: Currently not working very well for non-interactive mode        #
    #                                                                         # 
    # Note that better RMS fits can be obtained by running the wavelength     #
    # calibration interactively and identifying all of the lines              #
    # manually.  Tedious, but will give more accurate results than the        #
    # automatic mode (i.e., fl_inter-).  Use fl_iner+ for manual mode.        #
    #                                                                         # 
    ###########################################################################

    # Determine the wavelength of the observation and set the arc coordinate
    # file.  If the user wishes to change the coordinate file to a different
    # one, they need only to change the "clist" variable to their line list
    # in the coordli= parameter in the nswavelength call.
    hdulist = fits.open("rgn"+os.path.splitext(arc)[0]+".fits")
    band = hdulist[0].header['GRATING'][0:1]


    if useSkyLines:
        # optionally use the sky lines instead of lamps
        clist=arcDir+"ohlines.dat"
        my_thresh=50
        nfound = 5
        nlost = 1
    else:
        nfound = 10
        nlost = 10
        if band == "Z":
            clist=arcDir+"ArXe_Z.dat"
            my_thresh=100.0
        elif band == "K":
            clist=arcDir+"ArXe_K.dat"
            my_thresh=50.0
        else:
            clist=arcDir+"argon.dat"
            my_thresh=100.0


    iraf.nswavelength("rgn"+arc, coordli=clist, nsum=10, thresho=my_thresh, 
       trace='yes',fwidth=2.0,match=-6,cradius=8.0,fl_inter='no',nfound=nfound,
       nlost=nlost, logfile=log)

    # return the name of the file with the wavelength solution
    return 'wrgn'+arc

###########################################################################
# STEP 2.5:                                                               #
#  Trace the spatial curvature and spectral distortion in the Ronchi flat #
###########################################################################

def mk_ronchi_flat(ronchilist,raw_dir,shift_image,flat_image,flat_dark,bpm,log):
    ronchiflat=str(open(ronchilist, "r").readlines()[0]).strip()
    
    iraf.nfprepare("@"+ronchilist,rawpath=raw_dir, shiftimage=shift_image,
                   bpm=bpm,
                   fl_vardq='yes',fl_corr='no',fl_nonl='no',logfile=log)

    # Determine the number of input Ronchi calibration mask files so that
    # the routine runs automatically for single or multiple files.
    nfiles = len(open(ronchilist).readlines())
    if nfiles > 1:
        iraf.gemcombine("n//@"+ronchilist,output="gn"+ronchiflat,fl_dqpr='yes',
                        masktype="none",fl_vardq='yes',logfile=log)
    else:
        iraf.copy("n"+ronchiflat+".fits","gn"+ronchiflat+".fits")

    iraf.nsreduce("gn"+ronchiflat, outpref="r",  dark=flat_dark,
                  flatimage=flat_image,
                  fl_cut='yes', fl_nsappw='yes', 
                  fl_flat='yes', fl_sky='no', fl_dark='yes', fl_vardq='no',
                  logfile=log)

    iraf.nfsdist("rgn"+ronchiflat, 
                 fwidth=6.0, cradius=8.0, glshift=2.8, 
                 minsep=6.5, thresh=2000.0, nlost=3, 
                 fl_inter='no',logfile=log)
    
    return "rgn"+ronchiflat

## ###########################################################################
## # Reset to user defaults                                                  #
## ###########################################################################

def nifs_finish():
    if user_clobber == "no":
        iraf.set(clobber='no')


def run_basecalib(datDir, raw_dir, reduce_dir, backupDir, arcDir, useSkyLines, flatlist,\
               flatdarklist, ronchilist, arcStr, arcStrDrk, steps_basec):

    
    flatfile1=str(open(flatlist, "r").readlines()[0]).strip()
    flatdark1=str(open(flatdarklist, "r").readlines()[0]).strip()
    ronchifile1=str(open(ronchilist, "r").readlines()[0]).strip()
    arcfile1=str(open(arcStr, "r").readlines()[0]).strip()


    # Create a log file and back up the previous one if it already exists
    log = 'Basecalib.log'
    if os.path.exists(log):
        t = time.localtime()
        app = "_"+str(t[0])+str(t[1]).zfill(2)+str(t[2]).zfill(2)+'_'+ \
            str(t[3]).zfill(2)+':'+str(t[4]).zfill(2)+':'+str(t[5]).zfill(2)
        shutil.move(log,log+app)

    # Reduce data in sequence

    prep_nifs(log)  # reset things
    
    # find the shift from the first file
    if steps_basec[0]:
        print 'Using first file to find the initial shift: '+flatfile1
        shift_image = mk_mdf_shift(flatfile1,raw_dir,log)
    else:
        shift_image = "s"+flatfile1

    # make the flat files
    if steps_basec[1]:
        flat_image, flat_dark, bpm = mk_flat(flatlist,flatdarklist,raw_dir,shift_image,log)
    else:
        flat_image = "rgn"+os.path.splitext(flatfile1)[0]+"_flat.fits"
        bpm = "rgn"+os.path.splitext(flatfile1)[0]+"_sflat_bpm.pl"
        flat_dark = "rgn"+os.path.splitext(flatdark1)[0]+"_dark.fits"

    # find the wavelength solution
    if steps_basec[2]:
        wave_file = mk_wavelength_solution(arcStr,arcStrDrk,raw_dir,shift_image,
                                           flat_image,bpm,log,
                                           useSkyLines=useSkyLines,
                                           arcDir=arcDir)
    else:
        wave_file = 'wrgn'+arcfile1

    # make the ronchi flats
    if steps_basec[3]:
        ronchi_file = mk_ronchi_flat(ronchilist,raw_dir,shift_image,flat_image,
                                     flat_dark,bpm,log)
    else:
        ronchi_file = "rgn"+ronchifile1

    # make the backup directory for the calibration files, in case
    # they get erased.
    backupDir = reduce_dir+'calib_backup/'

    if not(os.path.isdir(reduce_dir)):
        print "Directory not found, making directory:"+reduce_dir
        os.mkdir(reduce_dir)


    if not(os.path.isdir(backupDir)):
        print "Directory not found, making directory:"+backupDir
        os.mkdir(backupDir)

    # copy the reduced files


    print 'backing up ronchi flat: '+ronchi_file
    shutil.copyfile(ronchi_file,reduce_dir+ronchi_file) # ronchi flat
    shutil.copyfile(ronchi_file,backupDir+ronchi_file)
    
    # shutil.copyfile("wrgn"+arc,reduce_dir+'wavelength_solution.fits')
    print 'backing up wavelength solution: '+wave_file
    shutil.copyfile(wave_file,reduce_dir+wave_file)
    shutil.copyfile(wave_file,backupDir+wave_file) 

    # delete the old database if it's there, then insert the new one
    try:  
        shutil.copytree('database',reduce_dir+'database')  
    except OSError:
        shutil.rmtree(reduce_dir+'database')
        shutil.copytree('database',reduce_dir+'database')
    try:
        shutil.copytree('database',backupDir+'database')
    except OSError:
        shutil.rmtree(backupDir+'database')
        shutil.copytree('database',backupDir+'database')
    
    # bad pixel map
    print 'backing up bad pixel mask: '+'flat_bpm.pl'
    shutil.copyfile(bpm,reduce_dir+'flat_bpm.pl')
    shutil.copyfile(bpm,backupDir+'flat_bpm.pl')
    
    # flat
    print 'backing up flat: '+flat_image+'.fits'
    shutil.copyfile(flat_image+'.fits',reduce_dir+'flat.fits')
    shutil.copyfile(flat_image+'.fits',backupDir+'flat.fits')
    
    # shift file
    print "backing up shift file: shiftFile.fits"
    shutil.copyfile(shift_image,reduce_dir+'shiftFile.fits')
    shutil.copyfile(shift_image,backupDir+'shiftFile.fits')


    # copy the names of the files out into a file
    fileTypes = ['ronchiflat','badPixelMask','flat','shiftImage','waveSolution']
    fileTypeNames = [ronchi_file,bpm,flat_image,shift_image,wave_file]
    newfileTypeNames = [ronchi_file,'flat_bpm.pl','flat.fits','shiftFile.fits',wave_file]
    output = open(reduce_dir+'calibfiles','w')
    for ii in np.arange(len(fileTypes)):
        output.write('%s \t %s \t %s\n' % (fileTypes[ii],fileTypeNames[ii],\
                                           newfileTypeNames[ii]))
    output.close()
    
## ###########################################################################
## #		End of the Baseline Calibration reduction                    #
## ###########################################################################
## #	                                                                     #
## #  The final output files created from this script for later science      #
## #  reduction have prefixes and file names of:                             #
## #     1. Shift reference file:  "s"+calflat                               #
## #     2. Flat field:  "rn"+calflat+"_flat"                                #
## #     3. Flat BPM (for DQ plane generation):  "rn"+calflat+"_flat_bpm.pl" #
## #     4. Wavelength referenced Arc:  "wrgn"+arc                           #
## #     5. Spatially referenced Ronchi Flat:  "rgn"+ronchiflat              #
## #     6. A database for some info on the ronchi flat and arc: "database"  #
## #     For this reduction,                                                 #
## #        Shift ref. file =  sN20060210S0195.fits                          #
## #        Flat field      =  rgnN20060210S0195_flat.fits                   #
## #        Flat BPM        =  rgnN20060210S0195_sflat_bpm.pl                #
## #        Arc frame       =  wrgnN20060210S0191.fits                       #
## #        Ronchi flat     =  rgnN20060210S0389.fits                        #
## #        database        =  database     (a folder)                       #
## #	                                                                     #
## #  Because of the shutil calls at the end of the script, these required   #
## #  files are also copied to reduce_dir with new (more convenient) names.  #
## #  These names are:                                                       #
## #                                                                         #
## #        Shift ref. file =  shiftFile                                     #
## #        Flat field      =  flat                                          #
## #        Flat BPM        =  flat_bpm.pl                                   #
## #        Arc frame       =  wrgnN20060210S0191.fits                       #
## #        Ronchi flat     =  rgnN20060210S0389.fits                        #
## #        database        =  database                                      #
## #                                                                         #
## #  The arc frame and ronchi flat have to have the same names as the       #
## #  database references them with their original names.                    #
## #                                                                         #
## #  A file with the name "calibfiles" is also written to reduce_dir, so    #
## #  the program can find the names of the different calibration files.     #
## #                                                                         #
## #  These files are all the new files you need after having finished the   #
## #  base calibration.                                                      #
## #                                                                         #
## ###########################################################################












######################################################################
######################################################################
######################################################################
#STEP 3:
#PREPARING FOR TELLURIC CORRECTION
#
#To do the telluric correction in step 4 we must have a normalized telluric
#absorption spectrum. To get this we need to do the following:
#
#1. Combining all the different exposures of the A-star taken on the
#   night of observation to find a spectrum for the A-star observed.
#   (code originally written by R.McDermid, slightly modified)
#2. Changing the stellar template (in our case Vega, but other templates
#   can just as easily be used) into a usable form. There are three
#   things needed to be done. 
#   2a) Convolving the much higher-resolution template data to be
#       smooth enough for comparison to the A-star observation data.
#       This step also normalizes the template data.
#   2b) Fixing for any velocity the A-star in observation might have.
#   2c) Do a linear interpolation to force the data points from the
#       stellar template to be at the same wavelengths as the A-star
#       measurements.
#3. Dividing the observation A-star spectrum by the fixed template
#   spectrum, then normalizing again.
#
#####################################################################


###### STEP 3.1 (Creating A-star spectrum) ######



# Some gaussian fitters have to be defined for later use:

def gauss_function(x, offset, a, x0, sigma):
    return offset+a*np.exp(-(x-x0)**2/(2*sigma**2))

def test_fit(k, g):
    '''A simple 1d gaussian fitter.
    '''
    m = np.argmax(g)
    center = k[m]
    maxPt = np.max(g)
    dx = k[1]-k[0]
    integral = (g*dx).sum()
    if maxPt != 0.0:
        sigma = integral/maxPt/2.0
    else:
        sigma = 1
    popt, pcov = curve_fit(gauss_function, k, g, p0 = \
                           [np.min(g), maxPt, center, sigma])
    return popt

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a 2d gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = leastsq(errorfunction, params)
    return p


# Finished defining the gaussian fitters.





def astar_spec(raw_data, red_data, ManualNames, tellistfile, astarskylistfile):

    ###########################################################################
    # Gemini NIFS data reduction script
    # Reduction for:  TELLURIC STANDARD CALIBRATIONS
    #
    # DESCRIPTION:
    #
    # This script is a pyraf-compatible version of the original NIFS scripts
    # written by Tracy Beck. 
    #
    # Notes:
    # - A sky frame is constructed by median combining sky frames.  Editing
    # the way the sky subtraction is done should be easy if telluric data
    # were obtained by offsetting to the sky (In this case, just look at the
    # way the NIFS Science reduction is constructed).
    #
    #
    #
    # VERSION:
    #  1.0: Adapted from original cl scripts by Tracy Beck. R.McDermid, 05Sep2012
    #  2.0: Modified to be automatic with use of a template A-star, step 3.2. 
    #  M.Stostad, July2013
    #  3.0: Added to rest of pipeline. M.Stostad, July2013
    # AUTHOR: R.McDermid (rmcdermid@gemini.edu)
    # Modified by: M.Stostad (morten.stostad@mail.utoronto.ca)
    ###########################################################################

    ###########################################################################
    # STEP 3.1.1:  Prepare IRAF                                               #
    ###########################################################################

    # Import the pyraf module and relevant packages
    from pyraf import iraf
    iraf.gemini()
    iraf.nifs()
    iraf.gnirs()
    iraf.gemtools()
    import pyfits as fits

    # Unlearn the used tasks
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs,iraf.nifs)

    # Create a log file and back up the previous one if it already exists
    log = 'Telluric.log'
    if os.path.exists(log):
        t = time.localtime()
        app = "_"+str(t[0])+str(t[1]).zfill(2)+str(t[2]).zfill(2)+'_'+ \
            str(t[3]).zfill(2)+':'+str(t[4]).zfill(2)+':'+str(t[5]).zfill(2)
        shutil.move(log,log+app)

    iraf.set(stdimage='imt2048')

    # Prepare the package for NIFS
    iraf.nsheaders("nifs",logfile=log)

    # Set clobber to 'yes' for the script. This still does not make the gemini tasks
    # overwrite files, so you will likely have to remove files if you re-run the script.
    user_clobber=iraf.envget("clobber")
    iraf.reset(clobber='yes')


    # Use the first telluric frame as the base name for the combined telluric spectrum
    telluric=str(open(tellistfile, "r").readlines()[0]).strip()


    ###########################################################################
    # STEP 3.1.2:  Reduce the Telluric Standard                               #
    ###########################################################################

    # Get the names of the base calibration files created in step 2 (assuming the
    # file calibfiles was created as it should have been):

    if ManualNames == False:
        calib_names = open(red_data+'calibfiles','r')
        ronchiflat  = calib_names.readline().split()[2].strip()  #'rgn'+"N20..."+'.fits'
        bpm         = calib_names.readline().split()[2].strip()  #'flat_bpm.pl'
        calflat     = calib_names.readline().split()[2].strip()  #'flat'
        shiftimage  = calib_names.readline().split()[2].strip()  #'shiftFile'
        arc         = calib_names.readline().split()[2].strip()  #'wrgn'+"N20..."+'.fits'

    # Prepare the data
    iraf.nfprepare("@"+tellistfile,rawpath=raw_data,shiftim=red_data+shiftimage,
      bpm=red_data+bpm,fl_vardq='yes',
      fl_int='yes',fl_corr='no',fl_nonl='no')

    astarskyfile = open(astarskylistfile)
    astar_sky  = astarskyfile.readline().strip()

    iraf.nfprepare("@"+astarskylistfile,rawpath=raw_data,shiftim=red_data+shiftimage,
      bpm=red_data+bpm,fl_vardq='yes',
      fl_int='yes',fl_corr='no',fl_nonl='no')

    # Make a median combined sky from the offset frames, which uses the telluric
    # frame as the name, with _sky appended.
    iraf.gemcombine("n//@"+astarskylistfile,output="g"+astar_sky,
      fl_dqpr='yes',fl_vardq='yes',masktype="none",logfile=log)

    # Do the sky subtraction on all the individual frames. Read the list (then
    # get rid of stupid '\n' character returns) first.
    telluriclist=open(tellistfile, "r").readlines()
    telluriclist=[word.strip() for word in telluriclist]
    for image in telluriclist:
        iraf.gemarith ("n"+image, "-", "g"+astar_sky, "sn"+image, fl_vardq="yes", 
                         logfile=log)

    #reduce and flat field the data
    iraf.nsreduce("sn@"+tellistfile,outpref="r", 
                  flatim=red_data+calflat,
                  fl_cut='yes',fl_nsappw='yes',fl_vardq='yes',fl_sky='no',
                  fl_dark='no',fl_flat='yes',logfile=log)

    #fix bad pixels from the DQ plane
    iraf.nffixbad("rsn@"+tellistfile,outpref="b",logfile=log)
    print "The arc file is", arc
    print "The ronchi flat is", ronchiflat
    #derive the 2D to 3D spatial/spectral transformation
    iraf.nsfitcoords("brsn@"+tellistfile,outpref="f", fl_int='no',
                     lamptr=arc, sdisttr=ronchiflat,
                     logfile=log,lxorder=4,syorder=4)

    #apply the transformation determined in the nffitcoords step
    iraf.nstransform("fbrsn@"+tellistfile,outpref="t",logfile=log)

    # make cubes
    # Reformat the data into a 3-D datacube
    iraf.nifcube ("tfbrsn@"+tellistfile, logfile=log)


    # Extract 1D spectra from the 3D data:

    outpref = 'gx'
    summed_specs = None
    telcubes = open(tellistfile, 'r').readlines()

    # A spectrum has to be extracted for each data cube of the A star. This
    # spectrum is extracted by fitting a 2d gaussian to the median flux close
    # to the center of the spatial directions, and then taking the spectra of
    # all the pixels within the FWHM of this 2d gaussian. Then a median 
    # background ring is subtracted to create the summed spectrum. 
    for file in telcubes:
        file = "ctfbrsn" + file.strip()
        tel_data = fits.getdata(file)

        # Create the median brightness fits file:
        brightness = np.zeros((tel_data.shape[1],tel_data.shape[2])) #stand-in
        for i in range(tel_data.shape[2]):
            for j in range(tel_data.shape[1]):
                brightness[j,i] = np.median(tel_data[:,j,i])  #double loop over 
                                                          #every pixel

        hdulist = fits.open(file)
        hdu = hdulist[1]
        hdu.data = brightness

        try:
            hdulist.writeto('m'+file)
        except IOError:
            os.remove('m'+file)
            hdulist.writeto('m'+file)

        med_data = brightness

        ylength = len(med_data[:,0])
        xlength = len(med_data[0,:])
        # Doing a 2d gaussian fit around the center of the image:
        y = int(ylength/2)
        x = int(xlength/2)

        print "Finding position of the A star in file:", file
        g2d_range = y*4.0/5.0 # How many pixels should the radius of the 2d fitter be?

        print "Fit between spatial coordinates:", \
            (max(0, y - g2d_range), max(0, x - g2d_range)), \
            (min(ylength, y + g2d_range), min(xlength, x + g2d_range))
        fit = fitgaussian(med_data[max(0, y - g2d_range) : \
                                           min(ylength, y + g2d_range), \
                                           max(0, x - g2d_range): \
                                           min(xlength, x + g2d_range)])
        y = y - g2d_range + fit[1] 
        x = x - g2d_range + fit[2] # The fitter assumes the minimum point is (0,0).

        r = 2.0*(fit[3] + fit[4])/2 # Using a rough avg width to find the "radius"
                    # of the star.

        print "The star is at y =", y, "x =", x, ", and has a radius of", r, "pixels."

        hdr = fits.getheader(file, 1)
        spec_len = hdr['NAXIS3']


        # Find the list of locations within the radius specified earlier, then
        # extract the spectra from these locations and put them in spec_list.
        star_locs = []
        for j in range(int(y - r), int(y + r + 2)):    
            for i in range(int(x - r), int(x + r + 2)): 
                          #range only uses ints, better with too many
                          #than with too few (some extra calculations, though)

                distance = np.sqrt((j - y)**2 + (i - x)**2)
                if (distance <= r) and j >= 0 and i >= 0:
                    star_locs.append([j, i])
        spec_list = np.zeros((len(star_locs), spec_len))
        indx = 0
        for loc in star_locs:
             if loc[0] < tel_data.shape[1] and loc[1] < tel_data.shape[2]:
                 spec_list[indx] = np.array(tel_data[:, loc[0], loc[1]])
                 indx += 1


        rinner = r*1.5    # To find a median background spectrum from a ring
        rout = r*2.0    # around the star we need an inner and outer radius


        # Now do the exact same as above, but with the bacground ring instead
        # of the inner circle. This yields a list of background spectra which 
        # we can take the median background spectrum from.
        bg_locs = []
        for j in range(int(y - rout), int(y + rout + 2)):    
            for i in range(int(x - rout), int(x + rout + 2)): 
                distance = np.sqrt((j - y)**2 + (i - x)**2)
                if (distance > rinner) and (distance <= rout) and j >= 0 and i >= 0:
                    bg_locs.append([j, i])
        bg_locs = np.array(bg_locs)

        bgspec_list = np.zeros((len(bg_locs), spec_len))
        indx = 0
        for loc in bg_locs:
                bgspec_list[indx] = np.array(tel_data[:, loc[0], loc[1]])
                # Sometimes the pipeline fails here if the fit of the
                # A-star is bad, for whatever reason. If so you won't
                # have to run the first step of the pipeline again,
                # but you should check the A-star files to see why
                # the 2d Gaussian fitter didn't work.
                indx += 1
        median = np.median(bgspec_list, axis = 0)


        # Finally subtract the median background spectrum (multiplied by the number
        # of spectra summed over when finding the summed star spectrum) from the 
        # summed star spectrum.
        sum1 = np.sum(spec_list, 0) - median*len(spec_list)
        sum1 = sum1 / np.median(sum1)
        # summed_specs is the list of finished spectra found from this loop, one
        # for each .fits star file.
        if summed_specs == None:  
            summed_specs = sum1
        else:
            summed_specs = np.vstack((summed_specs, sum1))

    med_spec = np.zeros(len(summed_specs[0,:]))
    for column in range(len(summed_specs[0,:])):
        med_spec[column]= np.median(summed_specs[:,column])

    # Creating a new fits file with a new header (to reflect the change to 1d)
    # and our new data:
    cubename = "ctfbrsn" + open(tellistfile, 'r').readline().strip()
    oldhdr = fits.getheader(cubename, 1)

    hdu = fits.PrimaryHDU()
    hdr = hdu.header
    hdr['BITPIX'] = oldhdr['BITPIX']
    hdr['NAXIS'] = 1
    hdr['NAXIS1'] = len(med_spec)
    hdr['CTYPE'] = 'LINEAR'
    hdr['CRVAL1'] = oldhdr['CRVAL3']
    hdr['CRPIX1'] = oldhdr['CRPIX3']
    hdr['CDELT1'] = oldhdr['CD3_3']
    hdr['CD1_1'] = oldhdr['CD3_3']

    try:
        fits.writeto(outpref+cubename, med_spec, hdr)
        print 'Writing new fits file', outpref+cubename
    except IOError:
        print 'Removing old telluric 1d-file:', outpref+cubename
        os.remove(outpref+cubename)
        print 'Writing new fits file', outpref+cubename
        fits.writeto(outpref+cubename, med_spec, hdr)


    ###########################################################################
    # Reset to user defaults                                                  #
    ###########################################################################
    if user_clobber == "no":
        iraf.set(clobber='no')

    ###########################################################################
    #          End of the Telluric Calibration Data Reduction                 #
    #                                                                         #
    #  The output of this reduction script is a 1-D spectrum used for         #
    # telluric calibration of NIFS science data.  For this particular         #
    # reduction the output file name is "gxctfbrsn"+telluric, or:             #
    # gxctfbrsnN20100401S0138. The file prefixes are described below.         #
    #                                                                         #
    # g = gemcombined/gemarithed   n=nfprepared  s=skysubtracted              #
    # r=nsreduced  b = bad pixel corrected  f= run through nffitcoords        # 
    # t = nftransformed   x = extracted to a 1D spectrum  m = flat 3d cube    #
    #                                                                         #
    # This script is meant to be a guideline as a method of a typical data    #
    # reduction for NIFS frames.  Of course, NIFS PIs can add or skip steps   #
    # in this reduction as they deem fit in order to reduce their particular  #
    # datasets.                                                               #
    #                                                                         #
    # version. 3.0: The output name of the 1d telluric file is now set in the #
    # beginning of the pipeline. By default it is                             #
    # "telluric_norm_weighed.fits"                                            #
    #                                                                         #
    ###########################################################################

    return outpref+cubename


###### STEP 3.2 (Functions for changing the stellar template) #######



def velocity_fix(fits_template, fits_data, vel_template = 'v_astartemplate.fits'):
    '''The A star will have some velocity. Because the telluric absorption
    will remain at the correct place regardless, we will have to shift the 
    spectrum from the template to simulate the template having the same velocity.
    We find the velocity of the star by making a gaussian fit to both the data
    and the template at the Brackett-gamma line, then measuring the difference.
    '''

    # We find the wavelength for the Brackett-gamma line of the template first:
    dat = fits.getdata(fits_template)
    header = fits.getheader(fits_template)
    start = header['CRVAL1']
    spec_delt = header['CD1_1']
    no_bins = header['NAXIS1']

    # Making a linear fit from the areas surrounding the Brackett-gamma line,
    # as we have to correct for the downwards slope in the spectrum before fitting
    # a gaussian. We need to find the x-values corresponding to these areas - I have
    # chosen the wavelengths to be 21450-21870 angstrom, excluding 21550-21770 as the
    # domain of the Brackett-gamma. These wavelengths are semi-arbitrary and can
    # be changed at will.
    
    linfitmin = int((21050 - start) / spec_delt)
    xmin = int((21450 - start) / spec_delt)
    xmax = int((21870 - start) / spec_delt)
    linfitmax = int((22270 - start) / spec_delt)  
    xval = np.append(np.arange(linfitmin, xmin), np.arange(xmax, linfitmax))
    lindata = np.append(dat[linfitmin:xmin], dat[xmax:linfitmax])
    m, b = np.polyfit(xval, lindata, 1)
    fitted_dat = dat/(m*np.arange(no_bins) + b)
    # Now that we've fitted the data to a linear curve we can run the gaussian fitter:
    gaussfit = test_fit(np.arange(xmin,xmax), fitted_dat[xmin:xmax])
    x = gaussfit[2]
    wlength = start + x*spec_delt  # Done! Found the emission wavelength for the template.

    # Now the exact same procedure for the data:
    dat1 = fits.getdata(fits_data)
    dat1 = dat1/np.median(dat1)
    header1 = fits.getheader(fits_data)
    start1 = header1['CRVAL1']
    spec_delt1 = header1['CD1_1']
    no_bins1 = header1['NAXIS1']


    # The wavelengths have been adjusted slightly to account for some telluric absorption
    # lines that would have clouded the linear/gaussian fit.
    linfitmin1 = int((21350 - start1) / spec_delt1)
    xmin1 = int((21550 - start1) / spec_delt1)
    xmax1 = int((21770 - start1) / spec_delt1)
    linfitmax1 = int((21940 - start1) / spec_delt1)  
    xval1 = np.append(np.arange(linfitmin1, xmin1), np.arange(xmax1, linfitmax1))
    lindata1 = np.append(dat1[linfitmin1:xmin1], dat1[xmax1:linfitmax1])
    m1, b1 = np.polyfit(xval1, lindata1, 1)
    fitted_dat1 = dat1/(m1*np.arange(no_bins1) + b1)
    gaussfit1 = test_fit(np.arange(xmin1,xmax1), fitted_dat1[xmin1:xmax1])
    x1 = gaussfit1[2]
    wlength1 = start1 + x1*spec_delt1  # Found the wavelength for the template - now our data

    print "Template emission wavelength is", wlength,\
        ", data emission wavelength is", wlength1

    diff = wlength1 - wlength
    velocity = diff * 2.99792458E8 / wlength
    
    print "The velocity of the A star is", velocity, "m s^-1"
    return velocity


def prepare_data(fits_template, fits_data, R = 5000):
    '''Prepares the sigma for the convolution. Assumes all header info
    is in units of angstrom.
    '''
    #Find delta lambda (i.e. sigma) from the data info and R
    header_data = fits.getheader(fits_data)
    len_spec = header_data['NAXIS1']
    start = header_data['CRVAL1']
    bin_size = header_data['CD1_1']
    
    sigma_wv = (start + (len_spec*bin_size)/2)/R 

    #Now convert this into units of pixels of the star template spectrum
    header_template = fits.getheader(fits_template)
    bin_size_template = header_template['CD1_1']
    
    sigma = sigma_wv / bin_size_template

    return sigma

def convolve(fits_template, sigma, smooth_template = 'c_astartemplate.fits'):
    '''Convolute the 1d data from a fits file with the sigma specified. Sigma must be
    given in number of bins. This is done on the high resolution template spectrum so
    that it's smoothed out for the interpolation done later. 
    '''
    dat = fits.getdata(fits_template)
    header = fits.getheader(fits_template)
    len_spec = header['NAXIS1']
    start = header['CRVAL1']
    spec_delt = header['CD1_1']
    wx = np.arange(start, start + len_spec*spec_delt, spec_delt)

    # Use a gaussian of arbitrary amplitude - as long as the spectrum
    # is normalized after the convolution the amplitude is irrelevant
    con = np.convolve(gauss_function(np.arange(-100,101), \
                            0, 1/(sigma*np.sqrt(2*np.pi)), 0, sigma), dat)
    con1 = con[100:-100]

    # Creates a .fits file with the smoothed data.
    hdulist = fits.open(fits_template)
    hdu = hdulist[0]
    hdu.data = con1/np.median(con1)
    try:
        hdulist.writeto(smooth_template)
    except IOError:
        os.remove(smooth_template)
        hdulist.writeto(smooth_template)


def lin_interpol(fits_template, fits_data, velocity, fixed_template):
    '''Does a linear interpolation of the fits data, fixing it for velocity and
    making the data points conform to the format of the data. Because the template
    has a very high resolution a simple linear interpolation doesn't lose significant
    data.
    '''

    # To do the interpolation we need the x-axis values for the template and our
    # A star data. I call these x-axis values hres_spec (for the template) and
    # lres_spec (for our data).
    indx = 0
    dat = fits.getdata(fits_template)
    header = fits.getheader(fits_template)
    len_spec = header['NAXIS1']
    start = header['CRVAL1']
    spec_delt = header['CD1_1']
    hres_spec = np.arange(start, start + len_spec*spec_delt, spec_delt)
    indx = 0

    # Fixing the template spectrum for velocity:
    for wlength in hres_spec:
        wlength = wlength + wlength*velocity/(2.99792458E8)
        hres_spec[indx] = wlength
        indx += 1

    # Add _d at the end to signify that the variable is for the data and not 
    # for the template:
    indx_d = 0
    dat_d = fits.getdata(fits_data)
    header_d = fits.getheader(fits_data)
    len_spec_d = header_d['NAXIS1']
    start_d = header_d['CRVAL1']
    spec_delt_d = header_d['CD1_1']
    lres_spec = np.linspace(start_d, start_d + (len_spec_d-1)*spec_delt_d, \
                            len_spec_d)

    # Doing the interpolation. The new data will be the template data, but
    # interpolated so that every data point corresponds to a wavelength value
    # from our original lres_spec. An arbitrary data point would be
    # (lres_spec[i], new_data[i]) for any i.
    new_data = np.zeros(len_spec_d)
    indx = 0
    print 'Fitting stellar template to A star data format..'
    for x in lres_spec:
        for indxx1 in range(len(hres_spec)):
            if hres_spec[indxx1] > x:
                break
        indxx0 = indxx1 - 1
        x0 = hres_spec[indxx0]
        x1 = hres_spec[indxx1]
        y0 = dat[indxx0] 
        y1 = dat[indxx1]
        y = y0 + ((y1 - y0)*(x - x0))/(x1 - x0)
        new_data[indx] = y
        indx += 1
    print 'Stellar template data successfully modified.'

    # Creating a new fits file and changing the headers to their new values.
    hdulist = fits.open(fits_template)
    hdu = hdulist[0]
    hdu.header['NAXIS1'] = len(new_data)
    hdu.header['CRVAL1'] = lres_spec[0]
    hdu.header['CDELT1'] = lres_spec[1] - lres_spec[0]
    hdu.header['CD1_1'] = lres_spec[1] - lres_spec[0]
    hdu.data = new_data
    try:
        hdulist.writeto(fixed_template)
    except IOError:
        os.remove(fixed_template)
        hdulist.writeto(fixed_template)

    return new_data


def telluric_fits(fin_template, fits_data, telluric_norm = 'telluric_norm1.fits'):
    '''Creates the fits file of the final normalized telluric spectrum.
    '''
    # Because the finished template data should have the same wavelength values
    # and indices as the data, creating the actual telluric file is very simple.

    hdulist = fits.open(fits_data)
    hdu = hdulist[0]
    data = hdu.data/(fin_template)
    hdu.data = data/(np.median(data))
    try:
        hdulist.writeto(telluric_norm)
    except IOError:
        os.remove(telluric_norm)
        hdulist.writeto(telluric_norm)

    print "backing up telluric norm:", telluric_norm
    shutil.copyfile(telluric_norm,reduce_dir+telluric_norm)
    shutil.copyfile(telluric_norm,backupDir+telluric_norm)



def run_telluric(raw_data, red_data, ManualNames, tellistfile, astarskylistfile,\
                    fits_template, R, telluric_norm):
    '''Modifies the template, then creates the telluric spectrum.
    '''

    fits_data = astar_spec(raw_data, red_data, ManualNames, tellistfile, astarskylistfile)

    smooth_template = 'c_astartemplate.fits' # name after template is smoothed
    fixed_template = 'v' + smooth_template # template is velocity fixed and in right format

    # The A star will have some velocity - this function finds it and stores it for 
    # interpolation later.
    velocity = velocity_fix(fits_template, fits_data)

    # Finding the sigma (in pixel size) for the convolution and performs the convolution,
    # creating a fits file with name specified in smooth_template above.
    sigma = prepare_data(fits_template, fits_data, R)
    print "The sigma for the convolution (in pixels of high-res spectrum) is", sigma
    convolve(fits_template, sigma, smooth_template)

    # Does a linear interpolation such that the indexed data points in the template are
    # equal to the indexed data points from the A star, also fixing for velocity.
    fin_template = lin_interpol(smooth_template, fits_data, velocity, fixed_template)

    # Finally creates the .fits file using the spectrum from the A star
    telluric_fits(fin_template, fits_data, telluric_norm = telluric_norm)




######################################################################
######################################################################
######################################################################
# STEP 4:
# REDUCTION OF SCIENCE DATA
# 
# This is the final step that the others have all been preparing for.
# Having the final calibration files (calflat, arc, ronchiflat, shiftimage, 
# bpm and telluric_norm), this portion of the script creates the final 3D data
# cubes.
#
######################################################################
######################################################################
######################################################################

def run_science(raw_data, red_data, rmReducedData, calflat, arc, ronchiflat, \
                shiftimage, bpm, scilistfile, skylistfile, useSkyFile, telluric, \
                ManualNames, reduced_all):
    '''
    Gemini NIFS data reduction script
    Reduction for:  SCIENCE DATA
     
    
    - Merging of data cubes is not implemented here.
    
    '''
    
    ###############################################################
    # STEP 4.1:  PREPARE IRAF                                     #
    ###############################################################

    # Import some useful python utilities
    import sys
    import getopt
    import os
    import time
    import shutil
    # Import the pyraf module and relevant packages
    from pyraf import iraf
    iraf.gemini()
    iraf.nifs()
    iraf.gnirs()
    iraf.gemtools()
    import pyfits
    import numpy as np
    import glob

    # Unlearn the used tasks
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs,iraf.nifs)

    # Create a log file and back up the previous one if it already exists
    log = 'Science.log'
    if os.path.exists(log):
        t = time.localtime()
        app = "_"+str(t[0])+str(t[1]).zfill(2)+str(t[2]).zfill(2)+'_'+ \
            str(t[3]).zfill(2)+':'+str(t[4]).zfill(2)+':'+str(t[5]).zfill(2)
        shutil.move(log,log+app)

    iraf.set(stdimage='imt2048')

    # Prepare the package for NIFS
    iraf.nsheaders("nifs",logfile=log)

    ###############################################################
    # STEP 4.2: SET REDUCTION FILE NAMES AND PATHS                #
    ###############################################################

    # Set clobber to 'yes' for the script. This still does not make the gemini tasks
    # overwrite files, so you will likely have to remove files if you re-run the script.
    user_clobber=iraf.envget("clobber")
    iraf.reset(clobber='yes')


    if ManualNames == False:
        calib_names = open(red_data+'calibfiles','r')
        ronchiflat  = calib_names.readline().split()[2].strip()  #'rgn'+"N20..."+'.fits'
        bpm         = calib_names.readline().split()[2].strip()  #'flat_bpm.pl'
        calflat     = calib_names.readline().split()[2].strip()  #'flat'
        shiftimage  = calib_names.readline().split()[2].strip()  #'shiftFile'
        arc         = calib_names.readline().split()[2].strip()  #'wrgn'+"N20..."+'.fits'


    if rmReducedData:
        # remove the data
        print 'Removing previously reduced files'
        sciFiles = np.loadtxt(scilistfile,dtype='str')
        for k in np.arange(len(sciFiles)):
            existingFiles = glob.glob('*?'+sciFiles[k])
            for inFile in existingFiles:
                print 'REMOVING: '+inFile
                os.remove(inFile)

    ###########################################################################
    # STEP 4.3:  Reduce the Science Data                                      #
    ###########################################################################

    iraf.nfprepare("@"+scilistfile, rawpath=raw_data,
       shiftimage=red_data+shiftimage,fl_vardq='yes',
       bpm=red_data+bpm,
       logfile=log)

    if useSkyFile == None:
        iraf.nfprepare("@"+skylistfile, rawpath=raw_data,
                       shiftimage=red_data+shiftimage,fl_vardq='yes',
                       bpm=red_data+bpm,
                       logfile=log)
        nfiles = len(open(skylistfile).readlines())
        skyoutfile = str(open(skylistfile, "r").readlines()[0]).strip()
        if nfiles > 1:
            iraf.gemcombine("n//@"+skylistfile,output="gn"+os.path.splitext(skyoutfile)[0],
                            fl_dqpr='yes',fl_vardq='yes',masktype="none",logfile=log)
            useSkyFile = "gn"+skyoutfile
        else:
            useSkyFile = "n"+skyoutfile


    #############################################################
    #  DATA REDUCTION HINT -                                    #
    #                                                           #
    # At the present time, we found that there are problems     #
    # with the WCS coordinates.  The automatic sky              #
    # frame ID and subtraction does not work very well in       #
    # "nsreduce" for NIFS data reductions if the number of sky  #
    # frames does not equal the number if science frames.  As   #
    # a result, the sky subtraction in this script was set up   #
    # to work outside of the "nsreduce" call.  This should work #
    # for most modes of science acquisition.  However you do    #
    # have to ensure that each frame in 'scilist' has a         #
    # corresponding frame in the "skylist" file. If you share   #
    # sky frames between differerent science exposures, you will#
    # have to duplicate those in the skylist.                   #
    #############################################################

    # Read in the frame lists (removing '\n' line breaks from the strings)
    scilist=open(scilistfile, "r").readlines()
    scilist=[word.strip() for word in scilist]


    for i in range(len(scilist)):
        iraf.gemarith ("n"+scilist[i], "-", useSkyFile, "gn"+scilist[i], fl_vardq="yes", 
                         logfile=log)

    # Flat field and cut the data
    iraf.nsreduce("gn@"+scilistfile, fl_cut='yes', fl_nsappw='yes', fl_dark='no', fl_sky='no', 
       fl_flat='yes', flatimage=red_data+calflat,
       fl_vardq='yes',logfile=log)

    # Interpolate over bad pixels flagged in the DQ plane
    iraf.nffixbad("rgn@"+scilistfile,logfile=log)

    # Derive the 2D to 3D spatial/spectral transformation
    iraf.nsfitcoords("brgn@"+scilistfile,lamptransf=arc, 
                     sdisttransf=ronchiflat,logfile=log, fl_int='no', lxorder=4, syorder=4)

    # Apply the transformation determined in the nffitcoords step
    iraf.nstransform("fbrgn@"+scilistfile, logfile=log)



    #iraf.nftelluric("tfbrgn@"+scilistfile, telluric, logfile=log)

    # Reformat the data into a 3-D datacube
    iraf.nifcube ("tfbrgn@"+scilistfile, logfile=log)



    # correct the data for telluric absorption features

    telluric_data = fits.getdata(telluric)
    cubes = open(scilistfile, 'r').readlines()


    if not(os.path.isdir(reduced_all)):
        print "Directory not found, making directory:"+reduced_all
        os.mkdir(reduced_all)

    # Copy the reduced files
    for inx in range(len(cubes)):
        print 'backing up science files: '+"ctfbrgn"+cubes[inx].strip()
        shutil.copyfile("ctfbrgn"+cubes[inx].strip(), reduced_all + \
                        "ctfbrgn" + cubes[inx].strip())

    print "The telluric file is", telluric
    
    available = telluric_data > 0.1

    for file1 in cubes:

        file1 = "ctfbrgn" + file1.strip()
        hdulist = fits.open(file1)
        hdu = hdulist[1]
        cube_data = hdu.data

        print "Removing telluric absorption from", file1
        for j in range(len(cube_data[0,:,0])):
            for i in range(len(cube_data[0,0,:])):
                cube_data[available,j,i] = cube_data[available,j,i]/telluric_data[available]

        hdu.data = cube_data
        try:
            hdulist.writeto(reduced_all+'a'+file1)
        except IOError:
            os.remove(reduced_all+'a'+file1)
            hdulist.writeto(reduced_all+'a'+file1)





    ###########################################################################
    # Reset to user defaults                                                  #
    ###########################################################################
    if user_clobber == "no":
        iraf.set(clobber='no')

    ###########################################################################
    #          End of the Science Data Reduction                              #
    #                                                                         #
    # The output of this reduction is a set of 3-D data cubes that have been  #
    # sky subtracted, flat fielded, cleaned for bad pixels, telluric          #
    # corrected and rectified into a cohesive datacube format.  In the case   #
    # of this reduction, the final output files are called: actfbrgn+science, #
    # or: actfbrgnN20100401S0182.fits                                         #
    #     actfbrgnN20100401S0184.fits                                         #
    #     actfbrgnN20100401S0186.fits                                         #
    #     actfbrgnN20100401S0188.fits                                         #
    #                                                                         #
    # The meaning of the output prefixes are described below:                 #
    #                                                                         #
    # g = gemcombined   n=nfprepared  s=skysubtracted   r=nsreduced           #
    # b = bad pixel corrected  f= run through nffitcoords                     # 
    # t = nftransformed   a = corrected for telluric absorption features      #
    # c = rectified to a 3D datacube                                          #
    #                                                                         #
    # This script is meant to be a guideline as a method of a typical data    #
    # reduction for NIFS frames.  Of course, NIFS PIs can add or skip steps   #
    # in this reduction as they deem fit in order to reduce their particular  #
    # datasets.                                                               #
    #                                                                         #
    ###########################################################################







def run_reduction():

    if steps[0] == True:
        run_basecalib(datDir, raw_dir, reduce_dir, backupDir, arcDir, useSkyLines, \
                flatlist, flatdarklist, ronchilist, arcStr, arcStrDrk, steps_basec)

    if steps[1] == True:
        run_telluric(raw_data, red_data, ManualNames, tellistfile, \
                         astarskylistfile, fits_template, R, telluric_norm)

    if steps[2] == True:
        run_science(raw_data, red_data, rmReducedData, calflat, arc, ronchiflat, \
                    shiftimage, bpm, scilistfile, skylistfile, useSkyFile, telluric, \
                    ManualNames, reduced_all)

