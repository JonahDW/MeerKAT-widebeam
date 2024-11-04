# MeerKAT-widebeam

This module contains several files that are all in one way or another related to correcting the effects of the primary beam on MeerKAT radio images. 

## meerkat_widebandpbcor.py

Calculate and apply primary beam correction for wideband data using primary beams and multiple frequencies. These can be generated using katbeam, or a series of primary beam images can be given at different frequencies. Several modes are supported, chosen through the parameter `mfs_mode`, related to how the image was generated and the file structure. Several usage examples are given in the `examples` directory.

The `mfs_mode` options are: 
- `none` : A simple primary beam correction at a single (central) frequency. 
- `casa` : Primary beam models are fit with the a Taylor polynomial and applied to the corresponding Taylor term images. If the number of Taylor terms is two or more, a spectral index image is automatically generated. The procedure to do the primary beam correction is meant to emulate the procedure as it is done in the CASA [`widebandpbcor`](https://casa.nrao.edu/docs/taskref/widebandpbcor-task.html) task. 
- `wsclean`: Assumes the images were created using the [WSClean](https://wsclean.readthedocs.io/en/latest/) MFS imaging mode, and a weighted primary beam is created using the weights from the individual channel images. 
- `wsclean-poly`: Just as `wsclean`, assumes the images were created using WSClean MFS imaging mode, but with `-fit-spectral-pol` option engaged, fitting a polynomial to the data. The primary beam is generated is thus generated in a similar manner to the `casa` option.
- `weighted`: Create a weighted primary beam model (like `wsclean` mode), manually specifying the frequencies and corresponding weights. 

Two different types of primary beam models are supported in the `model` keyword, either `katbeam` (default) which uses [katbeam](https://github.com/ska-sa/katbeam) models, or `images` where the primary beam model images on file can be specified. This latter option is useful, for instance, when you have generated holographic primary beam models with the `make_holo_beams.py` script (see below). Beware that when you supply lists of model images, weights, and frequencies, the order must be correct! The script can not and will not check beyond whether the lists are the same length.

```
usage: meerkat_widebandpbcor.py [-h] [-m MODEL] [-b BAND] [-t THRESHOLD] [-T]
                                [--freqs FREQS [FREQS ...]]
                                [--weights WEIGHTS [WEIGHTS ...]]
                                [--model_images MODEL_IMAGES [MODEL_IMAGES ...]]
                                [--nterms NTERMS] [--correct_channels]
                                [--alpha_thresh ALPHA_THRESH] [--write_beams]
                                [--outdir OUTDIR]
                                image_name mfs_mode

positional arguments:
  image_name            An image to correct. In 'casa' mode CASA file
                        structure is assumed, i.e. image_name.image.tt0,
                        image_name.image.tt1. In 'wsclean' mode, wsclean
                        output will be assumed, i.e.
                        imagename-0000-image.fits, imagename-0001-image.fits
                        for channel images and imagename-MFS-image.fits for
                        the mfs image.
  mfs_mode              Which multi-frequency synthesis mode has been used for
                        the image to determine the primary beam model. 'casa'
                        assumes Taylor term imaging with casa file structure,
                        while 'wsclean' assumes weighted imaging with a
                        wsclean structure. 'weighted' will assume weighted
                        imaging but without file structure, so frequencies and
                        weights are input manually. To do primary beam
                        correction at a single frequency without wideband mfs,
                        choose 'none'

optional arguments:
  -h, --help            show this help message and exit
  -m MODEL, --model MODEL
                        Which primary beam model to use, options are 'katbeam'
                        or 'images', in the latter case input images must be
                        specified by model_images parameter (default=katbeam).
  -b BAND, --band BAND  If using katbeam, specify band for which primary beam
                        model will be used, can be L, UHF, or S (default = L).
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold (at central frequency) below which image
                        pixels will be set to blank values (nan). Use to
                        remove areas where the primary beam correction is
                        large (default=0.3).
  -T, --trim            Trim image outside valid region (set by --threshold)
                        to reduce size.
  --freqs FREQS [FREQS ...]
                        Frequencies (in MHz) to use to generate the primary
                        beam (katbeam) or frequencies corresponding to input
                        model images (images). If using Taylor terms, the
                        number of frequencies should be equal or greater than
                        the number of Taylor terms. For the most accurate
                        results should resemble the frequency structure of the
                        data used for imaging
  --weights WEIGHTS [WEIGHTS ...]
                        Weights associated with the input frequencies for a
                        weighted primary beam
  --model_images MODEL_IMAGES [MODEL_IMAGES ...]
                        Specify files with PB images to to fit wideband
                        primary beam. Corresponding frequencies must be given
                        with the freqs parameter.
  --nterms NTERMS       Number of Taylor coefficients for polynomial fit
  --correct_channels    In wsclean mode, also correct all channel images at
                        the same time.
  --alpha_thresh ALPHA_THRESH
                        Mask all pixels below this flux level in the spectral
                        index image (default=0).
  --write_beams         Write derived beams to fits files (default=do not
                        write files).
  --outdir OUTDIR       Output directory of images.
```

## make_holo_beams.py

Script to create MeerKAT holographic primary beam models for different frequencies and polarisations based on data from [de Villiers (2023)](https://archive-gw-1.kat.ac.za/public/repository/10.48479/wdb0-h061/index.html). These can be used as input model images in `meerkat_widebandpbcor.py`.

```
usage: make_holo_beams.py [-h] [--freqs FREQS [FREQS ...]] [--outdir OUTDIR]
                          [--stokes STOKES [STOKES ...]] [--cell CELL]
                          [--imsize IMSIZE]
                          holodir

positional arguments:
  holodir               Directory containing holographic images, assuming npz
                        format from data by de Villiers (2023).

optional arguments:
  -h, --help            show this help message and exit
  --freqs FREQS [FREQS ...]
                        Frequencies of output images, in MHz.
  --outdir OUTDIR       Output directory of images.
  --stokes STOKES [STOKES ...]
                        Beam polarization(s).
  --cell CELL           Cell size, in arcsec.
  --imsize IMSIZE       Image size, in pixels
```

## wideband_corrections.py

This script assumes a catalog created from an image that was corrected with a primary beam which did not take into account wideband effects. It takes a catalog of sources and calculates and applies residual corrections to the fluxes and spectral indices of sources. If imaging products are not available or wideband primary beam correction in the image plane is too much of a hassle, these residual corrections can do a decent job of correcting for the effects of the wideband primary beam. The derivation and motivation behind these corrections is given in [Wagenveld et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023A%26A...673A.113W/abstract).

```
usage: alpha_flux_corrections.py [-h] [-d DIST_COL]
                                 [-f [FLUX_COLS [FLUX_COLS ...]]]
                                 [-a [ALPHA_COLS [ALPHA_COLS ...]]]
                                 [--alpha ALPHA]
                                 catalog

positional arguments:
  catalog               Input catalog to correct spectral index and flux

optional arguments:
  -h, --help            show this help message and exit
  -d DIST_COL, --dist_col DIST_COL
                        Catalog column containing the distance of the source
                        to the pointing center (default = 'Sep_PC')
  -f [FLUX_COLS [FLUX_COLS ...]], --flux_cols [FLUX_COLS [FLUX_COLS ...]]
                        Flux column for to apply flux correction on (default =
                        none, don't apply flux corrections)
  -a [ALPHA_COLS [ALPHA_COLS ...]], --alpha_cols [ALPHA_COLS [ALPHA_COLS ...]]
                        Spectral index column(s) to apply corrections on.
                        (default = none, don't apply spectral index
                        corrections)
  --alpha ALPHA         Spectral index to calculate flux corrections for. Can
                        be either float, or column in the catalog (default =
                        -0.8)
```
