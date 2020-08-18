
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

#  load bias frames & take median
#######################################
biasfilenames = []
for i in range(1,10):
    biasfilename_i = '/Users/vanshree/cosmo/bias/bias-000'+str(i)+'.fit'
    biasfilenames = np.append(biasfilenames,biasfilename_i)
biasfilenames = np.append(biasfilenames,'/Users/vanshree/cosmo/bias/bias-0010.fit')


# fig = plt.figure(figsize=(20,10))
# plt.set_cmap('viridis')
# plt.show()

biasarr = []

for n in range(len(biasfilenames)):
    file = str(biasfilenames[n])
    hdubias = fits.open(file)
    imbias = hdubias[0].data
    hdbias = hdubias[0].header
        
    # f1 = fig.add_subplot(4,3,n+1)
    # vmedian = np.nanmedian(imbias)
    
    # f1.imshow(np.arcsinh(imbias))
    # plt.title(biasfilenames[n][10::])

    biasarr.append(imbias)

# median bias frame
plt.set_cmap('viridis')

print("Shape of stacked bias-frames = ",np.shape(biasarr)) #should be 5 x 80 x 512
biasmed = np.median(biasarr, axis=0)
print("Shape of the median image generated = ",np.shape(biasmed))
vmedian = np.nanmedian(biasmed)

# plt.title("Combined Median bias frame")
# plt.imshow(biasmed,vmin=0.5*np.abs(vmedian), vmax=1.5*np.abs(vmedian), aspect=2)
# plt.colorbar()
# plt.show()
print("Mean of "+ "combined bias frame" + " is " + str(np.mean(biasmed)) + " counts ")


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

    # ########## plotting 
    # theta = np.arange(21)/10*np.pi 
    # ct = np.cos(theta) 
    # st = np.sin(theta) 
    # r = 25

    # fig = plt.figure(figsize=(12,12))
    
    # ax = plt.subplot()
    # ax.imshow(image_data, cmap='viridis', origin='lower')
    # ax.set_xlabel('Right Ascension')
    # ax.set_ylabel('Declination')

    # for i in range(len(xcoord)):
    #     ax.plot(xcoord[i]+r*ct, ycoord[i]+r*st, color='red', lw=2)
    #     ax.text(xcoord[i]+30, ycoord[i]+30, str(i+1), color='red', fontsize=14)

    # plt.show()


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


input(':')
