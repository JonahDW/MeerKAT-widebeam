# PB-correction

This module contains several files that are all in one way or another related to correcting the effects of the primary beam on MeerKAT radio images. 

## mkat_widebandpbcor.py

Calculate and apply primary beam correction for wideband data using primary beams and multiple frequencies. These can be generated using katbeam, or a series of primary beam images can be given at different frequencies. In standard mode, these are fit with the desired number of Taylor terms and applied to the corresponding Taylor term images. If the number of Taylor terms is two or more, a spectral index image is automatically generated. The procedure to do the primary beam correction is meant to emulate the procedure as it is done in the CASA [`widebandpbcor`](https://casa.nrao.edu/docs/taskref/widebandpbcor-task.html) task. If using the `wsclean_mfs` option, it's assumed the images were created using the [WSClean](https://wsclean.readthedocs.io/en/latest/) MFS imaging mode, and a weighted primary beam is created using the weights from the individual spectral windows.

```
usage: mkat_widebandpbcor.py [-h] [--wsclean_mfs] [--model MODEL]
                             [--model_images MODEL_IMAGES [MODEL_IMAGES ...]]
                             [--nterms NTERMS] [--write_beams]
                             [--freqs FREQS [FREQS ...]] [-b BAND] [-t THRESH]
                             [-T] [--alpha_thresh ALPHA_THRESH]
                             image_name

positional arguments:
  image_name            An image to correct. Standard input assumes CASA file
                        structure, i.e. image_name.image.tt0,
                        image_name.image.tt1. If wsclean-mfs is set to true,
                        wsclean output will be assumed, i.e.
                        imagename-0000-image.fits, imagename-0001-image.fits
                        for channel images and imagename-MFS-image.fits for
                        the mfs image. In this case frequencies will be
                        obtained from the channel images and the freqs option
                        ignored.

optional arguments:
  -h, --help            show this help message and exit
  --wsclean_mfs         Assume wsclean mfs weighting scheme instead of Taylor
                        term images. Currently only works with katbeam model.
  --model MODEL         Which primary beam model to use, options are katbeam,
                        plumber, and holo(graphic) (default=katbeam).
  --model_images MODEL_IMAGES [MODEL_IMAGES ...]
                        If using plumber or holo model, specify files with PB
                        images to to fit wideband primary beam.
  --nterms NTERMS       Number of Taylor coefficients
  --write_beams         Write derived beams to fits files (default=do not
                        write files).
  --freqs FREQS [FREQS ...]
                        If using katbeam, which frequencies (in MHz) to use to
                        generate the primary beam. The number of frequencies
                        should be equal or greater than the number of Taylor
                        terms, and for more accurate results should resemble
                        the frequency structure of the data used for your
                        imaging (default=1285)
  -b BAND, --band BAND  If using katbeam, specify band for which primary beam
                        model will be used, can be L, UHF, or S (default = L).
  -t THRESH, --threshold THRESH
                        Threshold (at central frequency) below which image
                        pixels will be set to blank values (nan). Use to
                        remove areas where the primary beam correction is
                        large (default=0.3).
  -T, --trim            Trim image outside valid region (set by --threshold)
                        to reduce size.
  --alpha_thresh ALPHA_THRESH
                        Mask all pixels below this flux level in the spectral
                        index image (default=0).
```

## mkat_primary_beam_correct.py

This script, courtesy of Kenda Knowles, uses functionality from the [`katbeam`](https://github.com/ska-sa/katbeam) module to correct images with primary beam model from MeerKAT. Currently, L-, S-, and UHF-band are supported. The input image can be in FITS or CASA format.

```
usage: mkat_primary_beam_correct [-h] [-b BAND] [-t THRESH] [-T] image_file

positional arguments:
  image_file            An image to correct.

optional arguments:
  -h, --help            show this help message and exit
  -b BAND, --band BAND  Band for which primary beam model will be used, can be
                        L, UHF, or S (default = L).
  -t THRESH, --threshold THRESH
                        Threshold below which image pixels will be set to
                        blank values (nan). Use to remove areas where the
                        primary beam correction is large.
  -T, --trim            Trim image outside valid region (set by --threshold)
                        to reduce size.
```

## create_alpha_primary_beam.py

Creates a spectral index image from the CASA `tt0` and `tt1` images, and optionally applies primary beam correction to the images

```
usage: create_alpha_primary_beam.py [-h] [-t THRESHOLD] [--pbcor] tt0 tt1

positional arguments:
  tt0                   0th order taylor term image.
  tt1                   1st order taylor term image.

optional arguments:
  -h, --help            show this help message and exit
  -t THRESHOLD, --threshold THRESHOLD
                        Mask all pixels below this flux level.
  --pbcor               Correct for in band spectral index from primary beam.
                        Input should contain FWHM of the primary beam in
                        degrees.
```


## alpha_flux_corrections.py

Take a catalog of sources and calculates and applies residual corrections to the fluxes and spectral indices of sources. These residual corrections are necessary if the images from which they were derived had a large bandwidth.

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
