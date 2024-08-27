import os
import sys
import copy
import glob
import numpy as np

import argparse
from pathlib import Path

from katbeam import JimBeam
from astropy.wcs import WCS

import helpers

def calculate_beams(image, model, model_images, band, freqs, outdir):
    '''
    Get beams at different frequencies
    '''
    im_shape = image[0].data[0, 0].shape

    if model == 'katbeam':
        try:
            beam = JimBeam('MKAT-AA-%s-JIM-2020' % band.upper())
            print("Selected katbeam %s-band" % band.upper())
        except ValueError:
            print('Band selected not recognized, check input carefully')
            sys.exit()

        # Get coordinates from image header and create coordinate map
        wcs = WCS(image[0].header, naxis=2).copy()
        RADeg = image[0].header['CRVAL1']
        decDeg = image[0].header['CRVAL2']
        xDegMap, yDegMap = helpers.makeXYDegreesDistanceMaps(np.ones(im_shape, dtype = np.float64), wcs, RADeg, decDeg)

        freqs = np.array(freqs, dtype=float)
        pbs = np.memmap(os.path.join(outdir,'temp_beams.dat'), dtype=np.float32, mode='w+',
                        shape=(len(freqs),len(xDegMap)*len(yDegMap)))
        for i, freqMHz in enumerate(freqs):
            pbs[i,:] = beam.I(xDegMap, yDegMap, freqMHz).flatten()

    if model == 'images':
        assert len(freqs) == len(model_images), "Different number of frequencies and model images!"

        pb_images = model_images
        pbs = np.memmap(os.path.join(outdir,'temp_beams.dat'), dtype=np.float32, mode='w+',
                        shape=(len(pb_images),im_shape[0]*im_shape[1]))
        freqs = np.array(freqs, dtype=float)
        for i, pb_im in enumerate(pb_images):
            imagedata = helpers.open_fits_casa(pb_im)
            pbs[i,:] = np.squeeze(imagedata[0].data).flatten()

    # Make sure frequencies are sorted
    sorted_freq = freqs.argsort()
    freqs = freqs[sorted_freq]
    pbs = pbs[sorted_freq,:]

    return pbs, freqs

def casa_widebandpbcor(in_image, pbs, nterms, freqAxis, freqs, thresh, alpha_thresh, trim, write_beams, model, outdir):
    '''
    Calculate the wideband PB correction for CASA images
    '''
    imagename = os.path.basename(in_image)
    out_image = os.path.join(outdir, imagename)

    tt0_img = helpers.open_fits_casa(in_image+'.image.tt0')
    wcs = WCS(tt0_img[0].header).copy()
    temp_img = copy.deepcopy(tt0_img)

    # Fit nth order polynomial to primary beams
    cfreq = tt0_img[0].header['CRVAL'+str(freqAxis)]/1e6 #MHz
    x = (freqs-cfreq)/cfreq

    pb_polyfit = np.polynomial.polynomial.polyfit(x, pbs, deg=nterms-1)
    # Renormalize tt0 to correct any small errors
    pb_polyfit[0,:] /= np.nanmax(pb_polyfit[0,:])

    # Apply PB for Taylor term 0 image
    pb_tt0 = pb_polyfit[0,:].reshape(tt0_img[0].data.shape[2],tt0_img[0].data.shape[3])
    tt0_img[0].data[0,0,:,:] = tt0_img[0].data[0,0,:,:]/pb_tt0
    tt0_untrimmed = copy.deepcopy(tt0_img)

    # Trim and write image
    tt0_img[0].data[0, 0][pb_tt0 < thresh] = np.nan
    if trim:
        tt0_img = helpers.trim_image(tt0_img, pb_tt0, thresh, trim)
    print(f"--> Saving primary beam corrected image '{out_image+'-pbcor_tt0.fits'}'")
    helpers.write_image_fits(tt0_img, out_image+'-pbcor_tt0', model)

    if write_beams:
        # Write PB images
        temp_img[0].data[0,0,:,:] = pb_tt0
        print(f"--> Saving tt0 primary beam 'pb_tt0.fits'")
        helpers.write_image_fits(temp_img, 'pb_tt0', model)

    if nterms > 1:
        tt1_img = helpers.open_fits_casa(in_image+'.image.tt1')
        pb_tt1 = pb_polyfit[1,:].reshape(tt1_img[0].data.shape[2],tt1_img[0].data.shape[3])

        #Apply PB correction
        tt1_img[0].data[0,0,:,:] = (tt1_img[0].data[0,0,:,:] - pb_tt1*tt0_untrimmed[0].data[0,0,:,:])/pb_tt0
        tt1_untrimmed = copy.deepcopy(tt1_img)

        # Create alpha image and apply PB correction
        alpha = tt1_untrimmed[0].data[0,0,:,:]/tt0_untrimmed[0].data[0,0,:,:]
        mask = tt0_untrimmed[0].data[0,0,:,:] < alpha_thresh
        alpha[mask] = np.nan

        # Open image so fits data structure is already there
        tt0_untrimmed[0].data[0,0,:,:] = alpha
        if trim:
            tt0_untrimmed = helpers.trim_image(tt0_untrimmed, pb_tt0, thresh, trim)
        print(f"--> Saving primary beam corrected spectral index image '{out_image+'-pbcor_alpha.fits'}'")
        helpers.write_image_fits(tt0_untrimmed, out_image+'-pbcor_alpha', model)

        if write_beams:
            temp_img[0].data[0,0,:,:] = pb_tt1
            print(f"--> Saving tt1 primary beam 'pb_tt1.fits'")
            helpers.write_image_fits(temp_img, 'pb_tt1', model)

            temp_img[0].data[0,0,:,:] = pb_tt1 / pb_tt0
            print(f"--> Saving alpha primary beam 'pb_alpha.fits'")
            helpers.write_image_fits(temp_img, 'pb_alpha', model)

def polynomial_widebandpbcor(in_image, pbs, nterms, freqAxis, freqs, thresh, trim, write_beams, model, outdir):
    '''
    Perform primary beam correction on channel images and combine them in a weighted average
    '''
    imagename = os.path.basename(in_image.rsplit('.',1)[0])
    out_image = os.path.join(outdir, imagename)

    mfs_img = helpers.open_fits_casa(in_image)
    wcs = WCS(mfs_img[0].header).copy()

    img_shape = np.squeeze(mfs_img[0].data).shape
    temp_img = copy.deepcopy(mfs_img)

    # Fit nth order polynomial to primary beams
    cfreq = mfs_img[0].header['CRVAL'+str(freqAxis)]/1e6 #MHz
    x = (freqs-cfreq)/cfreq

    pb_polyfit = np.polynomial.polynomial.polyfit(x, pbs, deg=nterms-1)
    # Renormalize tt0 to correct any small errors
    pb_polyfit[0,:] /= np.nanmax(pb_polyfit[0,:])

    # Apply PB for Taylor term 0 image
    mfs_pb = pb_polyfit[0,:].reshape(img_shape)
    mfs_img[0].data[0,0,:,:] = mfs_img[0].data[0,0,:,:]/mfs_pb

    if trim:
        mfs_img = helpers.trim_image(mfs_img, mfs_pb, thresh, trim)
    print(f"--> Saving primary beam corrected image '{out_image+'-pbcor.fits'}'")
    helpers.write_image_fits(mfs_img, out_image+'-pbcor', model)

    if write_beams:
        # Write PB images
        temp_img[0].data[0,0,:,:] = mfs_pb
        print(f"--> Saving mfs primary beam 'mfs_pb.fits'")
        helpers.write_image_fits(temp_img, 'mfs_pb', model)

def weighted_widebandpbcor(in_image, pbs, freqs, weights, thresh, trim, write_beams, model, outdir):
    '''
    Calculate the wideband PB correction with weighted images
    '''
    imagename = os.path.basename(in_image.rsplit('.',1)[0])
    out_image = os.path.join(outdir, imagename)

    mfs_img = helpers.open_fits_casa(in_image)
    wcs = WCS(mfs_img[0].header).copy()

    img_shape = np.squeeze(mfs_img[0].data).shape
    temp_img = copy.deepcopy(mfs_img)

    weights = np.array(weights).astype(float)
    mfs_pb = np.average(pbs, axis=0, weights=weights)
    mfs_pb = mfs_pb.reshape(img_shape)

    mfs_img[0].data[0,0,:,:] = mfs_img[0].data[0,0,:,:]/mfs_pb
    mfs_img[0].data[0, 0][mfs_pb < thresh] = np.nan
    if trim:
        mfs_img = helpers.trim_image(mfs_img, mfs_pb, thresh, trim)
    print(f"--> Saving primary beam corrected image '{out_image+'-pbcor.fits'}'")
    helpers.write_image_fits(mfs_img, out_image+'-pbcor', model)

    if write_beams:
        # Write PB images
        temp_img[0].data[0,0,:,:] = mfs_pb
        print(f"--> Saving mfs primary beam 'mfs_pb.fits'")
        helpers.write_image_fits(temp_img, 'mfs_pb', model)

def channel_pbcor(in_image, pb, freq, thresh, trim, model, outdir):
    '''
    Do primary beam correction on channel imagem, using single frequency
    '''
    imagename = os.path.basename(in_image.rsplit('.',1)[0])
    out_image = os.path.join(outdir, imagename)

    channel_img = helpers.open_fits_casa(in_image)
    wcs = WCS(channel_img[0].header).copy()

    img_shape = np.squeeze(channel_img[0].data).shape
    pb = pb.reshape(img_shape)

    channel_img[0].data[0,0,:,:] = channel_img[0].data[0,0,:,:]/pb
    return_image = channel_img[0].data[0,0,:,:]

    # Threshold and trim
    channel_img[0].data[0, 0][pb < thresh] = np.nan
    if trim:
        img = helpers.trim_image(channel_img, pb, thresh, trim)
    print(f"--> Saving primary beam corrected channel image '{out_image+'-pbcor.fits'}'")
    helpers.write_image_fits(channel_img, out_image+'-pbcor', model)

    return return_image

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    in_image = args.image_name
    mfs_mode = args.mfs_mode
    model = args.model
    band = args.band
    thresh = args.thresh
    trim = args.trim
    model_images = args.model_images
    freqs = args.freqs
    weights = args.weights
    invert_weights = args.invert_weights
    nterms = args.nterms
    alpha_thresh = args.alpha_thresh
    write_beams = args.write_beams
    outdir = args.outdir

    # Open mfs image for header information
    if 'wsclean' in mfs_mode.lower():
        mfs_image_file = in_image+'-MFS-image.fits'
    elif mfs_mode.lower() == 'casa':
        mfs_image_file = in_image+'.image.tt0'
    else:
        mfs_image_file = in_image

    if os.path.exists(mfs_image_file):
        mfs_img = helpers.open_fits_casa(mfs_image_file)
    else:
        print(f'Input image {mfs_image_file} not found')
        sys.exit(0)

    if outdir is None:
        outdir = os.path.dirname(in_image)

    # Determine frequency axis
    # Order depends on if e.g. CASA or WSClean image
    polAxis = None
    freqAxis = None
    for i in range(1, 5):
        if mfs_img[0].header['CTYPE%d' % (i)] == 'STOKES':
            polAxis=i
        elif mfs_img[0].header['CTYPE%d' % (i)] == 'FREQ':
            freqAxis=i
    assert(polAxis is not None and freqAxis is not None)
    freqMHz = mfs_img[0].header['CRVAL%d' % (freqAxis)]/1e6

    print("Image central Frequency = %.3f MHz" % (freqMHz))

    # In wsclean, determine frequencies and weights of channels
    if 'wsclean' in mfs_mode.lower():
        print("Reading frequencies and weights assuming wsclean outputs")

        freqs = []
        weights = []
        channel_images = sorted(glob.glob(in_image+'-0*-image.fits'))
        for channel_im_file in channel_images:
            chan_num = int(os.path.basename(channel_im_file).rsplit('-',2)[1])
            channel_im = helpers.open_fits_casa(channel_im_file)

            freq = channel_im[0].header['CRVAL'+str(freqAxis)]/1e6 #MHz
            weight = 1/channel_im[0].header['WSCIMGWG']*1e9
            if invert_weights:
                weight = 1/weight
            print(f"Channel image {chan_num} has central frequency {freq:.2f} MHz and weight {weight:.2g}")

            freqs.append(freq)
            weights.append(weight)

    # If using single frequency, determine frequency from image
    if mfs_mode.lower() == 'none':
        freqs = [freqMHz]

    # Get primary beam models
    pb_images, freqs = calculate_beams(mfs_img, model, model_images, band, freqs, outdir)

    # CASA primary beam correction
    if mfs_mode.lower() == 'casa':
        if len(freqs) < nterms:
            print("Number of frequencies is lower than number of Taylor terms, fit cannot be performed")
            sys.exit()
        casa_widebandpbcor(in_image, pb_images, nterms, freqAxis, freqs,
                           thresh, alpha_thresh, trim, write_beams, model, outdir)

    # WSClean polynomial primary beam correction
    elif mfs_mode.lower() == 'wsclean-poly':
        if len(freqs) < nterms:
            print("Number of frequencies is lower than number of Taylor terms, fit cannot be performed")
            sys.exit()
        polynomial_widebandpbcor(mfs_image_file, pb_images, nterms, freqAxis, freqs,
                                 thresh, trim, write_beams, model, outdir)

    # Single frequency primary beam correction
    elif mfs_mode.lower() == 'none':
        channel_image = channel_pbcor(in_image, pb_images[0], freqMHz, thresh, trim, model, outdir)
    # WSClean or weighted primary beam correction
    else:
        assert len(freqs) == len(weights), "Different number of frequencies and weights!"

        weighted_widebandpbcor(mfs_image_file, pb_images, freqs, weights, 
                               thresh, trim, write_beams, model, outdir)

    # Clean up
    os.remove(os.path.join(outdir,'temp_beams.dat'))


def new_argument_parser():

    parser = argparse.ArgumentParser()

    parser.add_argument("image_name", type=str, help="""An image to correct. In 'casa' mode 
                        CASA file structure is assumed, i.e. image_name.image.tt0, image_name.image.tt1.
                        In 'wsclean' mode, wsclean output will be assumed, i.e. imagename-0000-image.fits,
                        imagename-0001-image.fits for channel images and imagename-MFS-image.fits for
                        the mfs image.""")
    parser.add_argument("mfs_mode", type=str, help="""Which multi-frequency synthesis mode has been used
                        for the image to determine the primary beam model. 'casa' assumes Taylor term
                        imaging with casa file structure, while 'wsclean' assumes weighted imaging with
                        a wsclean structure. 'weighted' will assume weighted imaging but without file
                        structure, so frequencies and weights are input manually. To do primary beam 
                        correction at a single frequency without wideband mfs, choose 'none'""")
    parser.add_argument("-m", "--model", type=str, default="katbeam",
                        help="""Which primary beam model to use, options are 'katbeam' or 'images', in the latter
                                case input images must be specified by model_images parameter (default=katbeam).""")
    parser.add_argument("-b", "--band", dest="band", default="L", help="""If using katbeam, specify band 
                            for which primary beam model will be used, can be L, UHF, or S (default = L).""")
    parser.add_argument("-t", "--threshold", dest="thresh", default=0.3, type=float,
                        help="""Threshold (at central frequency) below which image pixels will be 
                                set to blank values (nan). Use to remove areas where the primary beam
                                correction is large (default=0.3).""")
    parser.add_argument("-T", "--trim", dest="trim", help="""Trim image outside valid region (set by
                        --threshold) to reduce size.""", default=False, action='store_true')
    parser.add_argument("--model_images", nargs='+', default=None,
                        help="""Specify files with PB images to to fit wideband primary beam. Corresponding
                                frequencies must be given with the freqs parameter.""")
    parser.add_argument("--freqs", nargs='+', default=None,
                        help="""Frequencies (in MHz) to use to generate the primary beam (katbeam) or frequencies
                        corresponding to input model images (images). If using Taylor terms, the number of 
                        frequencies should be equal or greater than the number of Taylor terms. 
                        For the most accurate results should resemble the frequency structure of 
                        the data used for imaging""")
    parser.add_argument("--weights", nargs='+', default=None,
                        help="""Weights associated with the input frequencies for a weighted primary beam""")
    parser.add_argument("--invert_weights", action='store_true', help="Invert the weights")
    parser.add_argument("--nterms", default=2, type=int, help="Number of Taylor coefficients")
    parser.add_argument("--alpha_thresh", default=0, type=float,
                        help="""Mask all pixels below this flux level in the spectral index image (default=0).""")
    parser.add_argument("--write_beams", action='store_true',
                        help="""Write derived beams to fits files (default=do not write files).""")
    parser.add_argument("--outdir", default=None, type=str, help="Output directory of images.")

    return parser

if __name__ == '__main__':
    main()