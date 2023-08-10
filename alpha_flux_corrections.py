import os
import sys
import json

import numpy as np
from astropy.table import Table

import matplotlib.pyplot as plt
from argparse import ArgumentParser
from pathlib import Path

class Catalog:

    def __init__(self, catalog_file):
        self.dirname = os.path.dirname(catalog_file)
        self.cat_name = os.path.basename(catalog_file).rsplit('.',1)[0]
        print(f'Reading in catalog: {self.cat_name}')

        self.table = Table.read(catalog_file)

        # Parse meta
        header = self.table.meta

        # Determine frequency axis:
        for i in range(1,5):
            if 'FREQ' in header['CTYPE'+str(i)]:
                freq_idx = i
                break
        self.freq = float(header['CRVAL'+str(freq_idx)])/1e9 #GHz
        self.dfreq = float(header['CDELT'+str(freq_idx)])/1e9 #GHz

        self.frequency_range = np.arange(self.freq-self.dfreq/2,self.freq+self.dfreq/2, 1e-3)

    def correct_alpha(self, alpha_cols, dist_col):
        '''
        Correct spectral indices as a function of distance with a previously obtained fit
        '''
        dist = self.table[dist_col]

        correction = alpha_correction(dist, self.frequency_range, self.freq)

        correction_col = 'Spectral_index_correction'
        # Go through alpha columns and apply correction
        for col in alpha_cols:
            self.table[col] -= correction

        if correction_col in self.table.colnames:
            sys.exit('Correction column is already present, indicating that a correction has already been applied. Exiting..')
        else:
            alpha_col_index = self.table.colnames.index(alpha_cols[0])
            self.table.add_column(correction,
                                  name=correction_col,
                                  index=alpha_col_index+1)

    def correct_flux(self, flux_cols, alpha, dist_col):
        '''
        Correct fluxes with primary beam, using spectral index
        '''
        dist = self.table[dist_col]
        if alpha:
            try:
                alpha = float(alpha)
            except ValueError:
                alpha = self.table[alpha]
        else:
            alpha = -0.75

        correction = flux_correction(dist, self.frequency_range, self.freq, alpha)

        correction_col = 'Flux_correction'
        # Go through flux columns and apply correction
        for col in flux_cols:
            self.table[col] /= correction
            # Also apply to error column
            if 'E_'+col in self.table.colnames:
                self.table['E_'+col] /= correction

        if correction_col in self.table.colnames:
            sys.exit('Correction column is already present, indicating that a correction has already been applied. Exiting..')
        else:
            flux_col_index = self.table.colnames.index(flux_cols[0])
            self.table.add_column(correction,
                                  name=correction_col,
                                  index=flux_col_index+1)

        return correction

def primary_beam_1D(beam_fwhm, b, freq, offset):
    '''
    MeerKAT L-band primary beam from Mauch et al 2019

    Keyword arguments:
    a (float) -- Constant scaling for the angle (in degrees)
    b (float) -- Constant scaling for the ratio offset/angle
    freq (float) -- Central frequency in GHz
    offset (array) -- Offsets for which to calculate the beam amplitude
    '''
    theta_beam = beam_fwhm*(1.5/freq)
    x = offset/theta_beam
    a_beam = (np.cos(b*np.pi*x)/(1-4*(b*x)**2))**2

    return a_beam

def flux_correction(sep_pc, frequency_range, central_freq, alpha):
    # Correct flux modified beam pattern induced by spectral shape
    meerkat_beam_fwhm = 57.5/60 #degrees

    # Make arrays and set them up for broadcasting
    freqs = frequency_range.reshape(-1,1)
    sep_pc = sep_pc.reshape(1,-1)

    if hasattr(alpha, '__len__'):
        alpha = alpha.filled(-0.8)
        alpha.reshape(1,-1)
    flux = freqs**(alpha)

    # Calculate reference flux at pointing center
    ref_beam = primary_beam_1D(meerkat_beam_fwhm, 1.189, freqs, 0)
    ref_flux = np.trapz(flux*ref_beam, freqs, axis=0)

    beam = primary_beam_1D(meerkat_beam_fwhm, 1.189, freqs, sep_pc)
    total_flux = np.trapz(flux*beam, freqs, axis=0)/ref_flux

    attenuation = primary_beam_1D(meerkat_beam_fwhm, 1.189, central_freq, sep_pc)

    correction = total_flux/attenuation

    return np.squeeze(correction)

def alpha_correction(sep_pc, frequency_range, central_freq):
    # Correct flux modified beam pattern induced by spectral shape
    meerkat_beam_fwhm = 57.5/60*(1.5/central_freq) #degrees

    # Make arrays and set them up for broadcasting
    freqs = frequency_range.reshape(1,-1)
    sep_pc = sep_pc.reshape(-1,1)

    alpha_int = np.trapz(-8*np.log(2)*(sep_pc/meerkat_beam_fwhm)**2*(freqs/1.27)**2,freqs,axis=1)
    alpha_basic = -8*np.log(2)*(sep_pc/meerkat_beam_fwhm)**2

    correction = alpha_int - alpha_basic[:,0]

    return correction

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    catalog_file = args.catalog
    dist_col = args.dist_col
    flux_cols = args.flux_cols
    alpha_cols = args.alpha_cols
    alpha = args.alpha

    # Initialize catalog
    catalog = Catalog(catalog_file)
    if alpha_cols:
        catalog.correct_alpha(alpha_cols, dist_col)

    if flux_cols:
        catalog.correct_flux(flux_cols, alpha, dist_col)

    catalog.table.write(catalog_file, overwrite=True)


def new_argument_parser():

    parser = ArgumentParser()

    parser.add_argument("catalog",
                        help="Input catalog to correct spectral index and flux")
    parser.add_argument("-d", "--dist_col", default='Sep_PC',
                        help="""Catalog column containing the distance of the source
                                to the pointing center (default = 'Sep_PC')""")
    parser.add_argument("-f", "--flux_cols", nargs='*',
                        help="""Flux column for to apply flux correction on
                                (default = none, don't apply flux corrections)""")
    parser.add_argument("-a", "--alpha_cols", nargs='*',
                        help="""Spectral index column(s) to apply corrections on.
                                (default = none, don't apply spectral index corrections)""")
    parser.add_argument("--alpha", default=None,
                        help="""Spectral index to calculate flux corrections for.
                                Can be either float, or column in the catalog (default = -0.8)""")

    return parser

if __name__ == '__main__':
    main()