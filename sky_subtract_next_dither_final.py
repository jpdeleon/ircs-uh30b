#!/usr/bin/env python
'''
Last run: 7/1 after using revised reduction_(7.1).py
nframes = 3; ndither = 5
'''
import glob
try:
  from astropy.io import fits as pf
except:
  import pyfits as pf

ch1=glob.glob('IRCA*fg_ch1.fits')
ch1.sort()
ch2=glob.glob('IRCA*fg_ch2.fits')
ch2.sort()

nframes = 3
wpangle = 4
ndither = 5
one_set=nframes*wpangle

#ch1
for m,i in enumerate(range(0,one_set)):
    im_1=pf.getdata(ch1[i]) #dith1
    im_2=pf.getdata(ch1[i+one_set]) #dith2
    head1 = pf.getheader(ch1[i])
    head2 = pf.getheader(ch1[i+one_set])
    # sky subtraction: dith1 - dith2
    im1 = im_1 - im_2 
    print('{0} - {1}'.format(head1['I_DTHPOS'], head2['I_DTHPOS']))
    pf.writeto(ch1[i][:-9]+'s_ch1.fits',im1,head1)
    # sky subtraction: dith2 - dith1
    im2 = im_2 - im_1
    print('{1} - {0}'.format(head1['I_DTHPOS'], head2['I_DTHPOS']))
    pf.writeto(ch1[i+one_set][:-9]+'s_ch1.fits',im2,head2)

for n,j in enumerate(range(2*one_set,3*one_set)):
    im_1=pf.getdata(ch1[j]) #dither 3 
    im_2=pf.getdata(ch1[j+one_set]) #dither 4
    head1 = pf.getheader(ch1[j])
    head2 = pf.getheader(ch1[j+one_set])
    # sky subtraction: dith3- dith4
    im1 = im_1 - im_2 
    print('{0} - {1}'.format(head1['I_DTHPOS'],head2['I_DTHPOS']))
    pf.writeto(ch1[j][:-9]+'s_ch1.fits',im1,head1)
    # sky subtraction: dith4- dith3
    im2 = im_2 - im_1 
    print('{1} - {0}'.format(head1['I_DTHPOS'],head2['I_DTHPOS']))
    pf.writeto(ch1[j+one_set][:-9]+'s_ch1.fits',im2,head2)

for o,j in enumerate(range(2*one_set,3*one_set)): #dith3
    im_1=pf.getdata(ch1[j+2*one_set]) #dither 5
    im_2=pf.getdata(ch1[j]) #dither 3
    head1 = pf.getheader(ch1[j+2*one_set])
    head2 = pf.getheader(ch1[j])
    # sky subtraction: dith5- dith3
    im = im_1 - im_2 
    print('{0} - {1}'.format(head1['I_DTHPOS'],head2['I_DTHPOS']))
    pf.writeto(ch1[j+2*one_set][:-9]+'s_ch1.fits',im,head1)

##################

#ch2
for m,i in enumerate(range(0,one_set)):
    im_1=pf.getdata(ch2[i]) #dith1
    im_2=pf.getdata(ch2[i+one_set]) #dith2
    head1 = pf.getheader(ch2[i])
    head2 = pf.getheader(ch2[i+one_set])
    # sky subtraction: dith1 - dith2
    im1 = im_1 - im_2 
    print('{0} - {1}'.format(head1['I_DTHPOS'], head2['I_DTHPOS']))
    pf.writeto(ch2[i][:-9]+'s_ch2.fits',im1,head1)
    # sky subtraction: dith2 - dith1
    im2 = im_2 - im_1
    print('{1} - {0}'.format(head1['I_DTHPOS'], head2['I_DTHPOS']))
    pf.writeto(ch2[i+one_set][:-9]+'s_ch2.fits',im2,head2)

for n,j in enumerate(range(2*one_set,3*one_set)):
    im_1=pf.getdata(ch2[j]) #dither 3 
    im_2=pf.getdata(ch2[j+one_set]) #dither 4
    head1 = pf.getheader(ch2[j])
    head2 = pf.getheader(ch2[j+one_set])
    # sky subtraction: dith3- dith4
    im1 = im_1 - im_2 
    print('{0} - {1}'.format(head1['I_DTHPOS'],head2['I_DTHPOS']))
    pf.writeto(ch2[j][:-9]+'s_ch2.fits',im1,head1)
    # sky subtraction: dith4- dith3
    im2 = im_2 - im_1 
    print('{1} - {0}'.format(head1['I_DTHPOS'],head2['I_DTHPOS']))
    pf.writeto(ch2[j+one_set][:-9]+'s_ch2.fits',im2,head2)

for o,j in enumerate(range(2*one_set,3*one_set)): #dith3
    im_1=pf.getdata(ch2[j+2*one_set]) #dither 5
    im_2=pf.getdata(ch2[j]) #dither 3
    head1 = pf.getheader(ch2[j+2*one_set])
    head2 = pf.getheader(ch2[j])
    # sky subtraction: dith5- dith3
    im = im_1 - im_2 
    print('{0} - {1}'.format(head1['I_DTHPOS'],head2['I_DTHPOS']))
    pf.writeto(ch2[j+2*one_set][:-9]+'s_ch2.fits',im,head1)



