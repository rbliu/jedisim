#!/bin/bash
for i in {191..205}
do
	python jedimaster_grid.py "physics_settings/production_1_grids_more/config_$i.txt"
	for f in "trial_$i/distorted_*/"
	do
		tar -zcvf "trial_$i/distorted.tar.gz" $f
		rm -r $f
	done

	for f in "trial_$i/stamp_*/"
	do
		tar -zcvf "trial_$i/stamp.tar.gz" $f
		rm -r $f
	done

	for f in "90_trial_$i/distorted_*/"
	do
		tar -zcvf "90_trial_$i/distorted.tar.gz" $f
		rm -r $f
	done

	for f in "90_trial_$i/stamp_*/"
	do
		tar -zcvf "90_trial_$i/stamp.tar.gz" $f
		rm -r $f
	done

	cp -R "trial_$i/" "/Volumes/rliu/production_1_grids/trial_$i"
	rm -R "trial_$i/"
	cp -R "90_trial_$i/" "/Volumes/rliu/production_1_grids/90_trial_$i"
	rm -R "90_trial_$i/"
done
