#!/bin/bash
set -e

# Image to be primary beam corrected
image='Images/my-image'

# Perform primary beam correction
freqs=(908.04 952.34 996.65 1043.46 1092.78 1144.61 1198.94 1255.79 1317.23 1381.18 1448.05 1519.94 1593.92 1656.20)
weights=(0.0378 0.0248 0.0580 0.0669 0.0683 0.0733 0.0259 0.0271 0.1170 0.1389 0.1463 0.0938 0.0368 0.0850)
python meerkat_widebandpbcor.py $image 'weighted' --model katbeam -t 0.05 -T --freqs ${freqs[@]} --weights ${weights[@]}
