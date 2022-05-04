#! /usr/bin/env python3

"""

Command line tool for applying primary beam correction to images

Gets the frequency to use from the image header

"""

import os
import sys
import numpy as np

import astropy.io.fits as pyfits
import astropy.stats as apyStats

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

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    inFileName = str(args.image_file)
    band = args.band
    thresh = args.thresh
    trim = args.trim

    # If the filename does not in .fits we shall attempt to make it conform
    # and assume it is a CASA image.
    img = open_fits_casa(inFileName)

    prefix = "pbcorr"
    if trim == True:
        prefix=prefix+"_trim"
    outDir=os.path.split(os.path.abspath(inFileName))[0]
    fileBaseName = os.path.basename(inFileName).split('.')[0]
    outFileName = os.path.join(outDir,fileBaseName+"_"+prefix+".fits")

    # Order depends on if e.g. CASA or WSClean image
    polAxis = None
    freqAxis = None
    for i in range(1, 5):
        if img[0].header['CTYPE%d' % (i)] == 'STOKES':
            polAxis=i
        elif img[0].header['CTYPE%d' % (i)] == 'FREQ':
            freqAxis=i
    assert(polAxis is not None and freqAxis is not None)
    shape = img[0].data[0, 0].shape
    wcs = astWCS.WCS(img[0].header, mode = 'pyfits').copy()
    freqMHz = img[0].header['CRVAL%d' % (freqAxis)]/1e6
    wcs.updateFromHeader()

    print("Frequency = %.3f MHz" % (freqMHz))
    print("Selected %s-band" % band.upper())

    try:
        beam = JimBeam('MKAT-AA-%s-JIM-2020' % band.upper())
    except ValueError:
        print('Band selected not recognized, check input carefully')
        sys.exit()

    RADeg, decDeg = wcs.getCentreWCSCoords()
    maxRDeg = 5.0
    xDegMap, yDegMap = nemoCython.makeXYDegreesDistanceMaps(np.ones(shape, dtype = np.float64), wcs, RADeg, decDeg, maxRDeg)
    I = beam.I(xDegMap, yDegMap, freqMHz)
    img[0].data[0, 0] = img[0].data[0, 0]/I
    img[0].data[0, 0][I < thresh] = np.nan

    if trim == True:
        # Assumes N is at the top, E at the left
        y, x = np.where(I >= thresh)
        yMin, yMax = y.min(), y.max()
        xMin, xMax = x.min(), x.max()
        _, decMin = wcs.pix2wcs(xMax, yMin)
        _, decMax = wcs.pix2wcs(xMin, yMax)
        decMinMax = np.array([decMin, decMax])
        yDecAbsMax = np.array([yMin, yMax])[np.argmax(abs(decMinMax))]
        RAMin, _ = wcs.pix2wcs(xMax, yDecAbsMax)
        RAMax, _ = wcs.pix2wcs(xMin, yDecAbsMax)
        clipDict = astImages.clipUsingRADecCoords(img[0].data[0, 0], wcs, RAMin, RAMax, decMin, decMax)
        img[0].data = np.zeros([1, 1, clipDict['data'].shape[0], clipDict['data'].shape[1]])
        img[0].data[0, 0] = clipDict['data']

        for i in range(1, 3):
            for k in WCSKeys:
                img[0].header['%s%d' % (k, i)]=clipDict['wcs'].header['%s%d' % (k, i)]

    # We may as well report RMS while we're at it
    d = img[0].data[np.isnan(img[0].data) == False]
    sigma = 1e6
    for i in range(10):
        mask = np.logical_and(np.greater(d, d.mean()-3*sigma), np.less(d, d.mean()+3*sigma))
        sigma = np.std(d[mask])
    print("3-sigma clipped stdev image RMS estimate = %.1f uJy" % (sigma*1e6))
    sbi = apyStats.biweight_scale(d, c = 9.0, modify_sample_size = True)
    print("Biweight scale image RMS estimate = %.1f uJy" % (sbi*1e6))

    # Add katbeam version to header
    img[0].header['PBCOR'] = 'katbeam-'+katbeam.__version__
    img.writeto(outFileName, overwrite = True)

def new_argument_parser():

    parser = argparse.ArgumentParser("mkat_primary_beam_correct")

    parser.add_argument("image_file", help="""An image to correct.""")
    parser.add_argument("-b", "--band", dest="band", default="L", help="""Band for which primary beam model 
                        will be used, can be L, UHF, or S (default = L).""")
    parser.add_argument("-t", "--threshold", dest="thresh", help="""Threshold below which image pixels
                        will be set to blank values (nan). Use to remove areas where the primary beam
                        correction is large.""", default = 0.3, type = float)
    parser.add_argument("-T", "--trim", dest="trim", help="""Trim image outside valid region (set by
                        --threshold) to reduce size.""", default = False, action = 'store_true')

    return parser

if __name__ == '__main__':
    main()