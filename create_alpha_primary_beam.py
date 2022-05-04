"""

Calculate spectral index image using taylor term images

"""
import sys
import os
import argparse
import numpy as np

from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clip

import casacore.images as pim
from scipy import ndimage

import matplotlib.pyplot as plt

def alpha_beam(alpha, header):
    # Determine frequency axis:
    for i in range(1,5):
        if 'FREQ' in header['CTYPE'+str(i)]:
            freq_idx = i
            break
    freq = float(header['CRVAL'+str(freq_idx)])/1e9 #GHz

    fwhm = 57.5*(1.5/freq)/60 #degrees

    center = np.array([header['CRPIX1'],header['CRPIX2']])
    dist = np.linalg.norm(np.indices(alpha.shape) - center[:,None,None] + 0.5, axis=0)
    dist *= max(header['CDELT1'],header['CDELT2'])

    pb_alpha = -8*np.log(2)*(dist/fwhm)**2

    return pb_alpha

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    tt0_imfile = args.tt0
    tt1_imfile = args.tt1
    threshold = float(args.threshold)
    pbcor = args.pbcor

    tt0 = pim.image(tt0_imfile)
    tt0.putmask(False)
    tt0.tofits('tt0.fits', velocity=False)
    tt0 = fits.open('tt0.fits')

    tt1 = pim.image(tt1_imfile)
    tt1.putmask(False)
    tt1.tofits('tt1.fits', velocity=False)
    tt1 = fits.open('tt1.fits')

    alpha = tt1[0].data[0,0,:,:]/tt0[0].data[0,0,:,:]
    mask = tt0[0].data[0,0,:,:] < threshold
    alpha[mask] = np.nan

    if pbcor:
        inband_alpha = alpha_beam(alpha, tt0[0].header)
        alpha -= inband_alpha

    tt0[0].data[0,0,:,:] = alpha

    outfile = os.path.join(tt0_imfile.split('.')[0]+'_alpha.fits')
    tt0.writeto(outfile, overwrite=True)

    # Clean up
    os.system('rm tt0.fits')
    os.system('rm tt1.fits')

def new_argument_parser():

    parser = argparse.ArgumentParser()

    parser.add_argument("tt0", help="""0th order taylor term image.""")
    parser.add_argument("tt1", help="""1st order taylor term image.""")
    parser.add_argument("-t", "--threshold", default=0,
                        help="Mask all pixels below this flux level.")
    parser.add_argument("--pbcor", action='store_true',
                        help="""Correct for in band spectral index from
                                primary beam. Input should contain FWHM
                                of the primary beam in degrees.""")

    return parser

if __name__ == '__main__':
    main()