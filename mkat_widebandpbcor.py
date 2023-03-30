import os
import sys
import copy
import glob
import numpy as np

import astropy.io.fits as pyfits
import astropy.stats as apyStats

import matplotlib.pyplot as plt

import argparse
from pathlib import Path

import katbeam

from katbeam import JimBeam
from astLib import *

import casacore.images as pim

WCSKeys=['NAXIS', 'CTYPE', 'CRPIX', 'CRVAL', 'CDELT', 'CUNIT']

def open_fits_casa(file):
    '''
    Open an image in fits or CASA format and return the image
    '''
    if '.fits' in file.lower():
        imagedata = pyfits.open(file)
    else:
        image = pim.image(file)
        image.putmask(False)
        image.tofits('temp.fits', velocity=False)
        imagedata = pyfits.open('temp.fits')
        # Clean up
        os.system('rm temp.fits')

    return imagedata

def makeXYDegreesDistanceMaps(image, wcs, RADeg, decDeg):
    """
    Returns an array of distance along x, y axes in degrees from given position.
    Inspired by original function from: https://github.com/simonsobs/nemo/blob/main/nemo/nemoCython.pyx
    """
    x0, y0 = wcs.wcs2pix(RADeg, decDeg)
    ra0, dec0 = RADeg, decDeg
    ra1, dec1 = wcs.pix2wcs(x0+1, y0+1)
    xPixScale = astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
    yPixScale = astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)

    Y = image.shape[0]
    X = image.shape[1]

    # Real space map of angular distance in degrees
    xDegreesMap = np.ones([Y, X], dtype=np.float64)
    yDegreesMap = np.ones([Y, X], dtype=np.float64)

    xx, yy = np.meshgrid(np.arange(X), np.arange(Y))
    xDegreesMap = (xx-x0)*xPixScale
    yDegreesMap = (yy-y0)*yPixScale

    return xDegreesMap, yDegreesMap

def trim_image(trimg, wcs, freqMHz, band, thresh, trim):

    try:
        beam = JimBeam('MKAT-AA-%s-JIM-2020' % band.upper())
    except ValueError:
        print('Band selected not recognized, check input carefully')
        sys.exit()

    shape = trimg[0].data[0, 0].shape
    RADeg, decDeg = wcs.getCentreWCSCoords()
    xDegMap, yDegMap = makeXYDegreesDistanceMaps(np.ones(shape, dtype = np.float64), wcs, RADeg, decDeg)

    central_I = beam.I(xDegMap, yDegMap, freqMHz)
    trimg[0].data[0,0][central_I < thresh] = np.nan

    if trim:
        # Assumes N is at the top, E at the left
        y, x = np.where(central_I >= thresh)
        yMin, yMax = y.min(), y.max()
        xMin, xMax = x.min(), x.max()
        _, decMin = wcs.pix2wcs(xMax, yMin)
        _, decMax = wcs.pix2wcs(xMin, yMax)
        decMinMax = np.array([decMin, decMax])
        yDecAbsMax = np.array([yMin, yMax])[np.argmax(abs(decMinMax))]
        RAMin, _ = wcs.pix2wcs(xMax, yDecAbsMax)
        RAMax, _ = wcs.pix2wcs(xMin, yDecAbsMax)
        clipDict = astImages.clipUsingRADecCoords(trimg[0].data[0, 0], wcs, RAMin, RAMax, decMin, decMax)
        trimg[0].data = np.zeros([1, 1, clipDict['data'].shape[0], clipDict['data'].shape[1]])
        trimg[0].data[0, 0] = clipDict['data']

        for i in range(1, 3):
            for k in WCSKeys:
                trimg[0].header['%s%d' % (k, i)]=clipDict['wcs'].header['%s%d' % (k, i)]

    return trimg

def calculate_widebandpb(in_image, model, model_images, nterms, band, freqAxis, freqs, write_beams):
    '''
    Calculate the wideband PB correction
    '''
    tt0_img = open_fits_casa(in_image+'.image.tt0')
    cfreq = tt0_img[0].header['CRVAL'+str(freqAxis)]/1e6 #MHz
    dfreq = tt0_img[0].header['CDELT'+str(freqAxis)]

    temp_img = tt0_img
    im_shape = temp_img[0].data[0, 0].shape

    if model == 'katbeam':
        try:
            beam = JimBeam('MKAT-AA-%s-JIM-2020' % band.upper())
        except ValueError:
            print('Band selected not recognized, check input carefully')
            sys.exit()

        if len(freqs) < nterms:
            print('Number of frequencies is lower than number of Taylor terms, fit cannot be performed')
            sys.exit()

        # Get coordinates from image header and create coordinate map
        wcs = astWCS.WCS(tt0_img[0].header, mode='pyfits').copy()
        wcs.updateFromHeader()
        RADeg, decDeg = wcs.getCentreWCSCoords()
        xDegMap, yDegMap = makeXYDegreesDistanceMaps(np.ones(im_shape, dtype = np.float64), wcs, RADeg, decDeg)

        freqs = np.array(freqs, dtype=float)
        pbs = np.zeros(shape=(len(freqs),len(xDegMap)*len(yDegMap)))
        for i, freqMHz in enumerate(freqs):
            pbs[i,:] = beam.I(xDegMap, yDegMap, freqMHz).flatten()

    if model == 'plumber':
        pb_images = model_images

        pbs = np.zeros(shape=(len(pb_images),im_shape[0]*im_shape[1]))
        freqs = np.zeros(len(pb_images))
        for i, pb_im in enumerate(pb_images):
            imagedata = open_fits_casa(pb_im)
            pbs[i,:] = imagedata[0].data[0,0].flatten()

            freqs[i] = float(pb_im.split('_')[-1].split('M')[0])

    if model == 'holo':
        pb_images = model_images

        pbs = np.zeros(shape=(len(pb_images),im_shape[0]*im_shape[1]))
        freqs = np.zeros(len(pb_images))
        for i, pb_im in enumerate(pb_images):
            imagedata = open_fits_casa(pb_im)

            pbs[i,:] = imagedata[0].data[0,0].flatten()
            freqs[i] = (imagedata[0].header['CRVAL3']-imagedata[0].header['CRPIX3']*imagedata[0].header['CDELT3'])/1e6 #MHz

    # Make sure frequencies are sorted
    sorted_freq = freqs.argsort()
    freqs = freqs[sorted_freq]
    pbs = pbs[sorted_freq,:]

    y = pbs
    x = (freqs-cfreq)/cfreq
    pb_polyfit = np.polynomial.polynomial.polyfit(x, y, deg=nterms-1)

    # Renormalize tt0 to correct any small errors
    pb_polyfit[0,:] /= np.nanmax(pb_polyfit[0,:])

    if write_beams:
        # Write PB images
        temp_img[0].data[0,0,:,:] = pb_polyfit[0,:].reshape(im_shape[0], im_shape[1])
        temp_img.writeto('pb_tt0.fits', overwrite=True)
        if nterms > 1:
            temp_img[0].data[0,0,:,:] = pb_polyfit[1,:].reshape(im_shape[0], im_shape[1])
            temp_img.writeto('pb_tt1.fits', overwrite=True)

            temp_img[0].data[0,0,:,:] = (pb_polyfit[1,:].reshape(im_shape[0], im_shape[1])
                                        / pb_polyfit[0,:].reshape(im_shape[0], im_shape[1]))
            temp_img.writeto('pb_alpha.fits', overwrite=True)

    return pb_polyfit[:nterms,:]

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    in_image = args.image_name
    model = args.model
    model_images = args.model_images
    nterms = args.nterms
    write_beams = args.write_beams
    freqs = args.freqs
    band = args.band
    thresh = args.thresh
    trim = args.trim
    alpha_thresh = args.alpha_thresh

    # Open tt0 image for header information
    tt0_img = open_fits_casa(in_image+'.image.tt0')

    # Determine frequency axis
    # Order depends on if e.g. CASA or WSClean image
    polAxis = None
    freqAxis = None
    for i in range(1, 5):
        if tt0_img[0].header['CTYPE%d' % (i)] == 'STOKES':
            polAxis=i
        elif tt0_img[0].header['CTYPE%d' % (i)] == 'FREQ':
            freqAxis=i
    assert(polAxis is not None and freqAxis is not None)
    freqMHz = tt0_img[0].header['CRVAL%d' % (freqAxis)]/1e6

    wcs = astWCS.WCS(tt0_img[0].header, mode='pyfits').copy()
    wcs.updateFromHeader()

    print("Image central Frequency = %.3f MHz" % (freqMHz))
    print("Selected %s-band" % band.upper())

    prefix = "pbcorr"
    outDir=os.path.split(os.path.abspath(in_image))[0]
    fileBaseName = os.path.basename(in_image).split('.')[0]

    # Calculate wideband Taylor terms
    pb_tt = calculate_widebandpb(in_image, model, model_images, nterms, band, freqAxis, freqs, write_beams)

    tt0_img = open_fits_casa(in_image+'.image.tt0')
    pb_tt0 = pb_tt[0,:].reshape(tt0_img[0].data.shape[2],tt0_img[0].data.shape[3])

    # Apply PB correction
    tt0_img[0].data[0,0,:,:] = tt0_img[0].data[0,0,:,:]/pb_tt0
    tt0_untrimmed = copy.deepcopy(tt0_img)

    # Add beam model to header and write to file
    if model == 'katbeam':
        tt0_img[0].header['PBCOR'] = 'katbeam-'+katbeam.__version__
    else:
        tt0_img[0].header['PBCOR'] = model

    tt0_img = trim_image(tt0_img, wcs, freqMHz, band, thresh, trim)
    outFileName = os.path.join(outDir,fileBaseName+"_"+prefix+"_tt0.fits")
    tt0_img.writeto(outFileName, overwrite = True)

    if nterms > 1:
        tt1_img = open_fits_casa(in_image+'.image.tt1')
        pb_tt1 = pb_tt[1,:].reshape(tt1_img[0].data.shape[2],tt1_img[0].data.shape[3])

        #Apply PB correction
        tt1_img[0].data[0,0,:,:] = (tt1_img[0].data[0,0,:,:] - pb_tt1*tt0_untrimmed[0].data[0,0,:,:])/pb_tt0
        tt1_untrimmed = copy.deepcopy(tt1_img)

        # Create alpha image and apply PB correction
        alpha = tt1_untrimmed[0].data[0,0,:,:]/tt0_untrimmed[0].data[0,0,:,:]
        mask = tt0_untrimmed[0].data[0,0,:,:] < alpha_thresh
        alpha[mask] = np.nan

        # Open image so fits data structure is already there
        tt0_untrimmed[0].data[0,0,:,:] = alpha
        if model == 'katbeam':
            tt0_untrimmed[0].header['PBCOR'] = 'katbeam-'+katbeam.__version__
        else:
            tt0_untrimmed[0].header['PBCOR'] = model

        tt0_untrimmed = trim_image(tt0_untrimmed, wcs, freqMHz, band, thresh, trim)
        outFileName = os.path.join(outDir,fileBaseName+"_"+prefix+"_alpha.fits")
        tt0_untrimmed.writeto(outFileName, overwrite = True)

def new_argument_parser():

    parser = argparse.ArgumentParser()

    parser.add_argument("image_name", type=str, help="""An image to correct. Input assumes 
                        CASA file structure, i.e. image_name.image.tt0, image_name.image.tt1, etc.""")
    parser.add_argument("--model", type=str, default="katbeam",
                        help="""Which primary beam model to use, options are katbeam, plumber, 
                                and holo(graphic) (default=katbeam).""")
    parser.add_argument("--model_images", nargs='+', default=None,
                        help="""If using plumber or holo model, specify files with PB images to 
                                to fit wideband primary beam.""")
    parser.add_argument("--nterms", default=2, type=int, help="Number of Taylor coefficients")
    parser.add_argument("--write_beams", action='store_true',
                        help="""Write derived beams to fits files (default=do not write files).""")
    parser.add_argument("--freqs", nargs='+', default=1285,
                        help="""If using katbeam, which frequencies (in MHz) to use to generate the primary beam.
                                The number of frequencies should be equal or greater than the number of Taylor terms,
                                and for more accurate results should resemble the frequency structure of the data
                                used for your imaging (default=1285)""")
    parser.add_argument("-b", "--band", dest="band", default="L", help="""If using katbeam, specify band 
                            for which primary beam model will be used, can be L, UHF, or S (default = L).""")
    parser.add_argument("-t", "--threshold", dest="thresh", default=0.3, type=float,
                        help="""Threshold (at central frequency) below which image pixels will be 
                                set to blank values (nan). Use to remove areas where the primary beam
                                correction is large (default=0.3).""")
    parser.add_argument("-T", "--trim", dest="trim", help="""Trim image outside valid region (set by
                        --threshold) to reduce size.""", default=False, action='store_true')
    parser.add_argument("--alpha_thresh", default=0, type=float,
                        help="""Mask all pixels below this flux level in the spectral index image (default=0).""")

    return parser

if __name__ == '__main__':
    main()