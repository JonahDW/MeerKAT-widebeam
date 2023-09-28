import os
import argparse
import numpy as np

import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

from astropy.io import fits

import helpers

def get_beam_pol(beams, pols, stokes='I'):
    '''
    Get specified polarisation of beam
    '''
    HH = pols == 'HH'
    HV = pols == 'HV'
    VH = pols == 'VH'
    VV = pols == 'VV'

    if stokes == 'I':
        beams_out = 0.5 * (np.sum(np.abs(beams)**2, axis=0))
    if stokes == 'Q':
        beams_out = 0.5 * (np.abs(beams[HH,:,:,:])**2
                         + np.abs(beams[HV,:,:,:])**2
                         - np.abs(beams[VH,:,:,:])**2
                         - np.abs(beams[VV,:,:,:])**2)
    if stokes == 'U':
        beams_out = np.real(beams[HH,:,:,:] 
                    * np.conjugate(beams[VH,:,:,:])
                    + beams[HV,:,:,:]
                    * np.conjugate(beams[VV,:,:,:]))
    if stokes == 'V':
        beams_out = np.imag(beams[HH,:,:,:] 
                    * np.conjugate(beams[VH,:,:,:])
                    + beams[HV,:,:,:]
                    * np.conjugate(beams[VV,:,:,:]))

    return beams_out

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    holodir = args.holodir
    freqs = args.freqs
    outdir = args.outdir
    stokes = args.stokes
    celldeg = args.cell/3600.
    imsize = args.imsize

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Read in files from directory
    beam_file = os.path.join(holodir,'beam.npy')
    freq_file = os.path.join(holodir,'freq_MHz.npy')
    pols_file = os.path.join(holodir,'pols.npy')
    margin_file = os.path.join(holodir,'margin_deg.npy')

    beams = np.load(beam_file, 'r')
    margin = np.load(margin_file)
    pols = np.load(pols_file).astype(str)
    in_freqs = np.load(freq_file)

    out_freqs = np.array(freqs).astype(float)
    for freq in out_freqs:
        nearest_freq = np.argmin(np.abs(in_freqs - float(freq)))
        print(f'Creating beam for {freq} MHz, nearest sampled frequency {in_freqs[nearest_freq]:.2f} MHz')

        for stoke in stokes:
            # Get beam and average over all antennas
            out_beams = get_beam_pol(beams[:,:,nearest_freq,:,:], pols, stokes=stoke)
            beam = np.mean(out_beams, axis=0)
            # Make sure beam is normalised
            beam /= beam.max()

            # Rescale beam to specified cell and image size
            scaling_func = RegularGridInterpolator((margin, margin), beam)
            new_margin = np.linspace(-imsize*celldeg/2, imsize*celldeg/2, imsize)
            XX, YY = np.meshgrid(new_margin, new_margin)
            scaled_beam = scaling_func((XX, YY))

            # Create and write to fits
            hdu = fits.PrimaryHDU(scaled_beam)
            hdul = fits.HDUList([hdu])
            outfile = os.path.join(outdir, f'{stoke}_{in_freqs[nearest_freq]:.2f}MHz.fits')
            hdul.writeto(outfile, overwrite=True)

def new_argument_parser():

    parser = argparse.ArgumentParser()

    parser.add_argument("holodir", type=str, help="""Directory containing holographic images, assuming
                        npz format from data by de Villiers (2023).""")
    parser.add_argument("--freqs", nargs='+', default=None, help="Frequencies of output images, in MHz.")
    parser.add_argument("--outdir", default='.', type=str, help="Output directory of images.")
    parser.add_argument("--stokes", nargs='+', default='I', help="Beam polarization(s).")
    parser.add_argument("--cell", default=2.0, type=float, help="Cell size, in arcsec.")
    parser.add_argument("--imsize", default=2000, type=int, help="Image size, in pixels")

    return parser

if __name__ == '__main__':
    main()