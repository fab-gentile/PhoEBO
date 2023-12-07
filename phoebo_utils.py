from astropy.io import fits
import numpy as np
from reproject import reproject_interp


def resize(image, start_hdr, final_hdr):
    hdu=fits.PrimaryHDU()
    hdu.header=start_hdr
    hdu.data=image
    array, footprint = reproject_interp(hdu, final_hdr)
    array[np.isnan(array)]=0
    return array

def crop_center(img,cropx,cropy):
    y,x = img.shape
    startx = x//2-(cropx//2)
    starty = y//2-(cropy//2)    
    return img[starty:starty+cropy,startx:startx+cropx]

def adu_to_ab(adu,zpt):
    adu=adu.astype(float)
    return -2.5*np.log10(adu) + zpt

def ab_to_ujy(ab):
    return 10**((23.9-ab)/2.5)

def adu_to_ab_err(adu_err,adu):
    adu=adu.astype(float)
    adu_err=adu_err.astype(float)
    return abs(2.5*(adu_err)/adu/np.log(10))

def ab_to_ujy_err(ab_err,ujy):
    return abs(ujy*ab_err/2.5*np.log(10))
