#!/bin/bash
set -e

# Directory containing data in .npz format from de Villiers (2023)
holo_dir='Holo-dir'
# Output directory for primary beam images
pb_images='Holo-PB-images'
# Image(s) to be primary beam corrected
image='Images/my-image'

# Create holographic primary beams
freqs=(2652.24 2706.92 2761.61 2816.30 2870.99 2925.67 2980.36 3035.05 3089.74 3144.42 3199.11 3253.80 3308.49 3363.17 3417.86)
python $pb_path'make_holo_beams.py' $holo_dir --freqs ${freqs[@]} --outdir $pb_images --cell 0.7 --imsize 8192

# Perform primary beam correction
python meerkat_widebandpbcor.py $image 'wsclean' --model 'images' --model_images ${pb_images[@]} --freqs ${freqs[@]} --write_beams -t 0.05 -T
