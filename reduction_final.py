#!/usr/bin/env python

from pyraf import iraf
import os
import sys
import glob
try:
    from astropy.io import fits as pf
except:
    import pyfits as pf
import logging

file_list=glob.glob('IRCA*.fits')
file_list.sort()

#logger
fname = pf.open(file_list[0])[0].header['FRAMEID']+'_reduction.log'
logging.basicConfig(filename=fname, level=logging.INFO)

ch1='[1:256,1:1024]' 	#o
ch2='[257:512,1:1024]' 	#e

distcorr_file = 'distcor.cl'
iraf.task(distcor = distcorr_file)
geomap='ircs+ao188_20mas_distmap_20131118.dbs'
flat = 'calflat.fits'

#dark_folder = 'data/imaging/ircs_UH30B/CALIB/DARK'
dark_dir = '../CALIB/DARK/master_darks'
if not os.path.exists(dark_dir):
    print('Dark folder does is empty: {}'.format(dark_folder))
    logging.warning('dark file does not exist')
    sys.exit()



for filename in file_list:
  #get header
  hdr=pf.open(filename)[0].header
  #get specific dark
  #exptime='master_dark_'+str(hdr['EXP1TIME'])+'s.fits'
  #dark=os.path.join(dark_dir,exptime)
  #dark_out = filename[:-5]+'d.fits'
  #try:
  #  iraf.imar(filename,'-',dark,dark_out) #should save header
  #except Exception as e:
  #  print(e)
  #  logging.warning(e)
  #  sys.exit()

  #flatfielding
  flat_out = filename[:-5]+'f.fits'
  try:
    iraf.imar(filename,'/',flat, flat_out)
  except Exception as e:
    print(e)
    logging.warning(e)
    sys.exit()
  logging.info('(1) Flat-fielding: {}'.format(flat))
  #geometric distortion correction
  geom_out = filename[:-5]+'fg.fits'
  try:
    iraf.distcor(flat_out,geom_out,geomap) #should save header
  except Exception as e:
    print(e)
    logging.warning(e)
    sys.exit()
  logging.info('(2) Distortion correction:\n{0} & {1}'.format(distcorr_file,geomap))
  #cropping
  a1=geom_out[:-5]+ch1
  b1=filename[:-5]+'fg_ch1'
  a2=geom_out[:-5]+ch2
  b2=filename[:-5]+'fg_ch2'
  try:
    iraf.imcopy(a1,b1) #should save header
    iraf.imcopy(a2,b2) #should save header
  except Exception as e:
    print(e)
    logging.warning(e)
    sys.exit()
  logging.info('(3) Cropped into o ({0}) and e channgels ({1})'.format(ch1,ch2))
  #remove intermediate files
  #iraf.imdel(dark_out)
  iraf.imdel(flat_out)
  iraf.imdel(geom_out)
