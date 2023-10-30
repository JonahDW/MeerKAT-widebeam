#!/bin/bash
set -e

# Image to be primary beam corrected
image='Images/my-image'

# Perform primary beam correction
python meerkat_widebandpbcor.py $image 'wsclean' --model katbeam -t 0.05 -T
