#!/bin/bash
set -e

images='*.image.tt0'

# Define frequencies in MHz
freqs=(882.8 936.3 989.8 1043.3 1096.8 1150.3 1203.8 1257.3 1310.8 1364.3 1417.8 1471.3 1524.8 1578.3 1631.8 1685.3)

for i in $images
do
	source=$(echo $i | cut -d'.' -f 1)
	python 'mkat_widebandpbcor.py' $source --model katbeam --freqs ${freqs[@]} -t 0.05 --trim --alpha_thresh 5e-5

done

