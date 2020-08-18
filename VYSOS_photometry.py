
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

rootlist = np.array([])
for i in range(1,10):# <---
    rootname = path+'0'+str(i)+str(filters[0])# <---
    rootlist = np.append(rootlist,rootname)

for i in range(10,16):# <---
    rootname = path+str(i)+str(filters[0])
    rootlist = np.append(rootlist,rootname)

print(rootlist)

#   
#######################################

for rootname in rootlist:
    imagename = rootname+'.fit'
    xysrcname = rootname+'-indx.xyls'
    
    image_data = fits.getdata(imagename)
    hdulxy = fits.open(xysrcname)
    xycoords = hdulxy[1].data
    xcoord = [x for x,y in xycoords]
    ycoord = [y for x,y in xycoords]

    #print(xysrcname)


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

