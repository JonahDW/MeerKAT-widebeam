# MeerKAT-widebeam

This module contains several files that are all in one way or another related to correcting the effects of the primary beam on MeerKAT radio images. 

## meerkat_widebandpbcor.py

Calculate and apply primary beam correction for wideband data using primary beams and multiple frequencies. These can be generated using katbeam, or a series of primary beam images can be given at different frequencies. Several modes are supported, related to how the image was generated and the file structure. In `casa` mode, primary beams are fit with the a Taylor polynomial and applied to the corresponding Taylor term images. If the number of Taylor terms is two or more, a spectral index image is automatically generated. The procedure to do the primary beam correction is meant to emulate the procedure as it is done in the CASA [`widebandpbcor`](https://casa.nrao.edu/docs/taskref/widebandpbcor-task.html) task. If using the `wsclean` mode, it's assumed the images were created using the [WSClean](https://wsclean.readthedocs.io/en/latest/) MFS imaging mode, and a weighted primary beam is created using the weights from the individual spectral windows. The `weighted` mode assumes a weighted primary beam as well, but the frequencies and corresponding weights can be specified manually.

```
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
                        weights are input manually.

optional arguments:
  -h, --help            show this help message and exit
  -m MODEL, --model MODEL
                        Which primary beam model to use, options are katbeam,
                        plumber, and holo(graphic) (default=katbeam).
  -b BAND, --band BAND  If using katbeam, specify band for which primary beam
                        model will be used, can be L, UHF, or S (default = L).
  -t THRESH, --threshold THRESH
                        Threshold (at central frequency) below which image
                        pixels will be set to blank values (nan). Use to
                        remove areas where the primary beam correction is
                        large (default=0.3).
  -T, --trim            Trim image outside valid region (set by --threshold)
                        to reduce size.
  --model_images MODEL_IMAGES [MODEL_IMAGES ...]
                        If using plumber or holo model, specify files with PB
                        images to to fit wideband primary beam.
  --freqs FREQS [FREQS ...]
                        If using katbeam, which frequencies (in MHz) to use to
                        generate the primary beam. If using Taylor terms, the
                        number of frequencies should be equal or greater than
                        the number of Taylor terms. For the most accurate
                        results should resemble the frequency structure of the
                        data used for imaging
  --weights WEIGHTS [WEIGHTS ...]
                        Weights associated with the input frequencies for a
                        weighted primary beam
  --nterms NTERMS       Number of Taylor coefficients
  --alpha_thresh ALPHA_THRESH
                        Mask all pixels below this flux level in the spectral
                        index image (default=0).
  --write_beams         Write derived beams to fits files (default=do not
                        write files).
```

## wideband_corrections.py

This script assumes a catalog created from an image that was corrected with a primary beam which did not take into account wideband effects. It takes a catalog of sources and calculates and applies residual corrections to the fluxes and spectral indices of sources, assuming a predefined spectral index. If imaging products are not available or primary beam correction in the image plane is too much of a hassle, these residual corrections can do a decent job of correcting for the effects of the wideband primary beam.

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
