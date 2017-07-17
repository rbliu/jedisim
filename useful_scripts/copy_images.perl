#!/usr/dist/bin/perl -w


open(OUT,">copy_images");

for ($i=0; $i<84; $i++){

$name1 = "/Volumes/rliu/jedisim_grid/trial_" . $i . "/trial" . $i . "_LSST_convolved_noise.fits";
$name2 = "/Volumes/rliu/jedisim_grid/90_trial_" . $i . "/90_trial" . $i . "_LSST_convolved_noise.fits";
$dir = "/Volumes/rliu/jedisim_grid/LSST_convolved_noise/";
print OUT "cp $name1 $dir\n";
print OUT "cp $name2 $dir\n";

}

close(OUT);