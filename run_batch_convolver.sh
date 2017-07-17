#!/bin/bash

for i in {0..41}
do
	python jediconvolver.py "/Volumes/rliu/production_1_grids/trial_$i/" "trial$[i]_HST.fits" "/Volumes/rliu/production_1_grids/trial_$i/" "physics_settings/psf_decam.fits"

done

for i in {0..41}
do
    python jediconvolver.py "/Volumes/rliu/production_1_grids/90_trial_$i/" "90_trial$[i]_HST.fits" "/Volumes/rliu/production_1_grids/90_trial_$i/" "physics_settings/psf_decam.fits"

done