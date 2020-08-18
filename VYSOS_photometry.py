
# Pipeline to do photometry on VYSOS targets
# Vanshree Bhalotia, 2020

################################################################################

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import photutils
from photutils import aperture_photometry
from photutils import CircularAperture
from photutils import CircularAnnulus

from astropy import units as u
from astropy.coordinates import SkyCoord

from astropy.io import fits
from astropy.wcs import WCS

################################################################################

# load images & astrometry outputs 
#######################################
path = '/usr/local/Cellar/astrometry-net/0.82/data/aug08-00' # <---
filters = ['V','R','i'] # <---

imagelist = np.array([])
sourcelist = np.array([])
for i in range(1,10):# <---
    imagename = path+'0'+str(i)+str(filters[0])+'.fit' # <---
    xyname = path+'0'+str(i)+str(filters[0])+'-indx.xyls'
    #print(imagename)
    #print(xyname)
    imagelist = np.append(imagelist,imagename)
    sourcelist = np.append(sourcelist,xyname)


for i in range(10,16):# <---
    imagename = path+str(i)+str(filters[0])+'.fit' # <---
    xyname = path+str(i)+str(filters[0])+'-indx.xyls'
    #print(imagename)
    #print(xyname)
    imagelist = np.append(imagelist,imagename)
    sourcelist = np.append(sourcelist,xyname)


print(imagelist)
print(sourcelist)


#   
#######################################

for imagename in imagelist:
    image_data = fits.getdata(imagename)

for xyname in sourcelist:
    hdulxy = fits.open(xyname)
    xycoords = hdulxy[1].data
    xcoord = [x for x,y in xycoords]
    ycoord = [y for x,y in xycoords]


#  load bias frames & take median
#######################################


#  perform  bias subtraction
#######################################

        ### show image

        ### show bias subtracted image

#  perform  aperture photometry
#######################################


#  convert to magnitudes 
#######################################


#  save outputs in a .txt file 
#######################################

