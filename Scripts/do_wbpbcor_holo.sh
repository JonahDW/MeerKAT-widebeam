#!/bin/bash
set -e

casa_python=/soft/astro/casa/casa-6.1.0-118/bin/python3
images='*.image.tt0'

# Regrid PB images
pb_images='StokesI_holo/*.chan32.im'

for i in $images
do
	for pbi in $pb_images
	do
		$casa_python StokesI_holo/imregrid.py $i $pbi -f
	done
	scaled_pb=(StokesI_holo/*.chan32_scaled.im)

	source=$(echo $i | cut -d'.' -f 1)
	python 'mkat_widebandpbcor.py' $source --model holo --model_images ${scaled_pb[@]} -t 0.05 --trim --alpha_thresh 5e-5
	rm -rf ${scaled_pb[@]}
done

rm *.log
