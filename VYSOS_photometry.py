
# Pipeline to do photometry on VYSOS targets
# Vanshree Bhalotia, 2020

#############################################################################################################

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

#############################################################################################################

#  load bias frames & take median
#######################################
biasfilenames = []
for i in range(1,10):
    biasfilename_i = '/Users/vanshree/cosmo/bias/bias-000'+str(i)+'.fit'
    biasfilenames = np.append(biasfilenames,biasfilename_i)
biasfilenames = np.append(biasfilenames,'/Users/vanshree/cosmo/bias/bias-0010.fit')

biasarr = []

for n in range(len(biasfilenames)):
    file = str(biasfilenames[n])
    hdubias = fits.open(file)
    imbias = hdubias[0].data
    hdbias = hdubias[0].header

    biasarr.append(imbias)

# print("Shape of stacked bias-frames = ",np.shape(biasarr)) #should be 5 x 80 x 512
biasmed = np.median(biasarr, axis=0)
# print("Shape of the median image generated = ",np.shape(biasmed))

print("Mean of "+ "combined bias frame" + " is " + str(np.mean(biasmed)) + " counts ")

#############################################################################################################

# load images & astrometry outputs

#############################################################################################################

path = '/usr/local/Cellar/astrometry-net/0.82/data/aug08-00' # <---
filters = ['V','R','i'] # <---

rootlist = np.array([])
for i in range(1,10):# <---
    rootname = path+'0'+str(i)+str(filters[0])# <---
    rootlist = np.append(rootlist,rootname)

for i in range(10,16):# <---
    rootname = path+str(i)+str(filters[0])
    rootlist = np.append(rootlist,rootname)

#######################################
# calculate src values for various images
#######################################

# sourcelist = np.array([])

# for rootname in rootlist:
#     xysrcname = rootname+'.rdls'
#     sourcelist = np.append(sourcelist,xysrcname)
    
#     hdulxy = fits.open(xysrcname)
#     xycoords = hdulxy[1].data
#     xcoord = [x for x,y in xycoords]
#     ycoord = [y for x,y in xycoords]

#     for ra in xcoord: 
    

#######################################
# image/sources one by one 
#######################################

starlist = np.array([])

for rootname in rootlist:
    
    imagename = rootname+'.fit'
    xysrcname = rootname+'.rdls' # loading RA and Dec
    print(xysrcname)

    hdul = fits.open(imagename)
    image_data = fits.getdata(imagename)
    im_header = hdul[0].header

    hdulxy = fits.open(xysrcname)
    xycoords = hdulxy[1].data
    xcoord = np.array([x for x,y in xycoords])
    ycoord = np.array([y for x,y in xycoords])

    for i in range(len(rootlist)):
        ra_i = xcoord[i]
        dec_i = ycoord[i]
        
        radiff = abs(ra_i-xcoord)
        decdiff = abs(dec_i-ycoord)

        ragap = 0.001666666667 #<-- 6 arseconds in degrees
        decgap = 0.001666666667 #<--
        ra_good = abs(ra_i - xcoord)<ragap
        dec_good = abs(dec_i - ycoord)<decgap

        print(xycoords[ra_good & dec_good])
        simstar = xycoords[ra_good & dec_good]
        starlist = np.append(starlist,simstar)
        # if abs(ra_i-xcoord)>0.5: #and abs(dec_i-ycoord)>0.5:
        #     continue


#  perform  bias subtraction
#######################################

    #print('image mean',np.mean(image_data))

    biassubtractedimage = image_data - biasmed

    #print('bias subtracted image mean',np.mean(biassubtractedimage))

#   
#######################################
    

#  perform  aperture photometry
#######################################
    xycoords = np.array(list(xycoords))

    aperture = CircularAperture(xycoords,r=5)
    annulus_aperture = CircularAnnulus(xycoords, r_in=8., r_out=10.)
    apers = [aperture, annulus_aperture]

    phot_table = aperture_photometry(biassubtractedimage,apers)
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g' 
    # print(phot_table)

    ### subtract background 
    bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
    bkg_sum = bkg_mean * aperture.area
    final_sum = phot_table['aperture_sum_0'] - bkg_sum
    phot_table['residual_aperture_sum'] = final_sum
    phot_table['residual_aperture_sum'].info.format = '%.8g'

    #print('\n')
    #print(phot_table['residual_aperture_sum']) 


#  convert to magnitudes  
#######################################
    exptime = im_header['EXPOSURE']
    # print('exptime',exptime)

    counts_array = np.array(phot_table['residual_aperture_sum'])
    time_array = np.ones(len(counts_array))*exptime

    inst_mags_array = -2.5*np.log10(counts_array/time_array)
    # print('inst_mags_array = ',inst_mags_array)


#  save outputs in a .txt file 
#######################################
    # print('xcenter=',phot_table['xcenter'])
    # print('ycenter=',phot_table['ycenter'])

    # print('aperture_sum_0 =',np.array(phot_table['aperture_sum_0']))
    # print('aperture_sum_1 =',np.array(phot_table['aperture_sum_1']))

    # print('residual_aperture_sum = ',np.array(phot_table['residual_aperture_sum']))
    # print('inst_mags_array =',inst_mags_array)


input(':')
