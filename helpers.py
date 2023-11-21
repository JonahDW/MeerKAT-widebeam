import numpy as np

import casacore.images as pim
from astropy.io import fits
from astropy.coordinates import angular_separation

import katbeam

def open_fits_casa(file):
    '''
    Open an image in fits or CASA format and return the image
    '''
    if '.fits' in file.lower():
        imagedata = fits.open(file)
    else:
        image = pim.image(file)
        image.putmask(False)
        image.tofits('temp.fits', velocity=False)
        imagedata = fits.open('temp.fits')
        # Clean up
        os.system('rm temp.fits')

    return imagedata

def makeXYDegreesDistanceMaps(image, wcs, RADeg, decDeg):
    """
    Returns an array of distance along x, y axes in degrees from given position.
    Inspired by original function from: https://github.com/simonsobs/nemo/blob/main/nemo/nemoCython.pyx
    """
    x0, y0 = wcs.all_world2pix(RADeg, decDeg)
    ra0, dec0 = RADeg, decDeg
    ra1, dec1 = wcs.all_pix2world(x0+1, y0+1)
    xPixScale = angular_separation(ra0, dec0, ra1, dec0)
    yPixScale = angular_separation(ra0, dec0, ra0, dec1)

    Y = image.shape[0]
    X = image.shape[1]

    # Real space map of angular distance in degrees
    xDegreesMap = np.ones([Y, X], dtype=np.float64)
    yDegreesMap = np.ones([Y, X], dtype=np.float64)

    xx, yy = np.meshgrid(np.arange(X), np.arange(Y))
    xDegreesMap = (xx-x0)*xPixScale
    yDegreesMap = (yy-y0)*yPixScale

    return xDegreesMap, yDegreesMap

def trim_image(trimg, beam_model, thresh, trim):

    # Assumes N is at the top, E at the left
    y, x = np.where(beam_model >= thresh)
    yMin, yMax = y.min(), y.max()
    xMin, xMax = x.min(), x.max()

    clipped_image = trimg[0].data[0, 0, yMin:yMax, xMin:xMax]
    trimg[0].data = np.zeros([1, 1, clipped_image.shape[0], clipped_image.shape[1]])
    trimg[0].data[0, 0] = clipped_image

    oldCRPIX1 = trimg[0].header['CRPIX1']
    oldCRPIX2 = trimg[0].header['CRPIX2']
    trimg[0].header['NAXIS1'] = clipped_image.shape[1]
    trimg[0].header['NAXIS2'] = clipped_image.shape[0]
    trimg[0].header['CRPIX1'] = oldCRPIX1 - xMin
    trimg[0].header['CRPIX2'] = oldCRPIX2 - yMin

    return trimg

def write_image_fits(image, filename, model):
    '''
    Write image to file
    '''
    # Add beam model to header
    if model == 'katbeam':
        image[0].header['PBCOR'] = 'katbeam-'+katbeam.__version__
    else:
        image[0].header['PBCOR'] = model

    if not filename.endswith('.fits'):
        filename += '.fits'
    image.writeto(filename, overwrite=True)