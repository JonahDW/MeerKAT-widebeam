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
import nemoCython

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

def trim_image(trimg, wcs, freqMHz, band, thresh):

    try:
        beam = JimBeam('MKAT-AA-%s-JIM-2020' % band.upper())
    except ValueError:
        print('Band selected not recognized, check input carefully')
        sys.exit()

    shape = trimg[0].data[0, 0].shape
    RADeg, decDeg = wcs.getCentreWCSCoords()
    maxRDeg = 5.0
    xDegMap, yDegMap = nemoCython.makeXYDegreesDistanceMaps(np.ones(shape, dtype = np.float64), wcs, RADeg, decDeg, maxRDeg)

    central_I = beam.I(xDegMap, yDegMap, freqMHz)
    trimg[0].data[0,0][central_I < thresh] = np.nan

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

def calculate_widebandpb(in_image, model, nterms, band, freqAxis, freqs):
    '''
    Calculate the wideband PB correction
    '''
    tt0_img = open_fits_casa(in_image+'.image.tt0')
    cfreq = tt0_img[0].header['CRVAL'+str(freqAxis)]
    dfreq = tt0_img[0].header['CDELT'+str(freqAxis)]
    im_shape = tt0_img[0].data[0, 0].shape

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
        maxRDeg = 5.0
        xDegMap, yDegMap = nemoCython.makeXYDegreesDistanceMaps(np.ones(im_shape, dtype = np.float64), wcs, RADeg, decDeg, maxRDeg)

        pbs = np.zeros(shape=(3,len(xDegMap),len(yDegMap)))
        for i, freqMHz in enumerate(freqs):
            pbs[i,:,:] = beam.I(xDegMap, yDegMap, freqMHz)

    else:
        pb_images = model

        pbs = np.zeros(shape=(len(pb_images),im_shape[0]*im_shape[1]))
        freqs = np.zeros(len(pb_images))
        for i, pb_im in enumerate(pb_images):
            imagedata = open_fits_casa(pb_im)

            pbs[i,:] = imagedata[0].data[0,0].flatten()
            freqs[i] = imagedata[0].header['CRVAL3']-imagedata[0].header['CRPIX3']*imagedata[0].header['CDELT3']

    x = (freqs-cfreq)/cfreq
    pb_polyfit  = np.polynomial.polynomial.polyfit(x, pbs, deg=nterms-1)

    return pb_polyfit[:nterms,:]

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    in_image = args.image_name
    model = args.model
    nterms = args.nterms
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
    pb_tt = calculate_widebandpb(in_image, model, nterms, band, freqAxis)

    tt0_img = open_fits_casa(in_image+'.image.tt0')
    pb_tt0 = pb_tt[0,:].reshape(tt0_img[0].data.shape[2],tt0_img[0].data.shape[3])

    # Apply PB correction
    tt0_img[0].data[0,0,:,:] /= pb_tt0
    tt0_untrimmed = copy.deepcopy(tt0_img)

    # Add beam model to header and write to file
    if model == 'katbeam':
        tt0_img[0].header['PBCOR'] = 'katbeam-'+katbeam.__version__
    else:
        tt0_img[0].header['PBCOR'] = 'HOLO'
    tt0_img = trim_image(tt0_img, wcs, freqMHz, band, thresh)
    outFileName = os.path.join(outDir,fileBaseName+"_"+prefix+"_tt0.fits")
    tt0_img.writeto(outFileName, overwrite = True)

    if nterms > 1:
        tt1_img = open_fits_casa(in_image+'.image.tt1')
        pb_tt1 = pb_tt[1,:].reshape(tt1_img[0].data.shape[2],tt1_img[0].data.shape[3])

        #Apply PB correction
        tt1_img[0].data[0,0,:,:] = (tt1_img[0].data[0,0,:,:] - pb_tt1*tt0_untrimmed[0].data[0,0,:,:])/pb_tt0
        tt1_untrimmed = copy.deepcopy(tt1_img)

        # Add beam model to header and write to file
        if model == 'katbeam':
            tt1_img[0].header['PBCOR'] = 'katbeam-'+katbeam.__version__
        else:
            tt1_img[0].header['PBCOR'] = 'HOLO'
        tt1_img = trim_image(tt1_img, wcs, freqMHz, band, thresh)
        outFileName = os.path.join(outDir,fileBaseName+"_"+prefix+"_tt1.fits")
        tt1_img.writeto(outFileName, overwrite = True)

        # Create alpha image
        alpha = tt1_untrimmed[0].data[0,0,:,:]/tt0_untrimmed[0].data[0,0,:,:]
        mask = tt0_untrimmed[0].data[0,0,:,:] < alpha_thresh
        alpha[mask] = np.nan

        # Open image so fits data structure is already there
        tt0_untrimmed[0].data[0,0,:,:] = alpha
        tt0_untrimmed = trim_image(tt0_untrimmed, wcs, freqMHz, band, thresh)
        outFileName = os.path.join(outDir,fileBaseName+"_"+prefix+"_alpha.fits")
        tt0_untrimmed.writeto(outFileName, overwrite = True)

def new_argument_parser():

    parser = argparse.ArgumentParser()

    parser.add_argument("image_name", type=str, help="""An image to correct. Input assumes 
                        CASA file structure, i.e. image_name.image.tt0, image_name.image.tt1, etc.""")
    parser.add_argument("--model", nargs='+', default="katbeam")
    parser.add_argument("--nterms", default=2, type=int, help="Number of Taylor coefficients")
    parser.add_argument("--freqs", nargs='+', default=1285,
                        help="""If using katbeam, which frequencies (in MHz) to use to generate the primary beam.
                                The number of frequencies should be equal or greater than the number of Taylor terms,
                                and for more accurate results should resemble the frequency structure of the data
                                used for your imaging""")
    parser.add_argument("-b", "--band", dest="band", default="L", help="""If using katbeam, 
                        specify band for which primary beam model will be used, can be L, UHF, or S (default = L).""")
    parser.add_argument("-t", "--threshold", dest="thresh", help="""Threshold (at central frequency) below which 
                        image pixels will be set to blank values (nan). Use to remove areas where the primary beam
                        correction is large.""", default=0.3, type=float)
    parser.add_argument("-T", "--trim", dest="trim", help="""Trim image outside valid region (set by
                        --threshold) to reduce size.""", default=False, action='store_true')
    parser.add_argument("--alpha_thresh", default=0, type=float,
                        help="""Mask all pixels below this flux level in the spectral index image.""")

    return parser

if __name__ == '__main__':
    main()