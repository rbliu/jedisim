#!/usr/dist/bin/perl -w


open(OUT,">clean_decam");

for ($i=0; $i<80; $i++){

$name1 = "/Volumes/rliu/jedisim_clusters/trial_" . $i . "/trial" . $i . "_decam_convolved*.fits";
$name2 = "/Volumes/rliu/jedisim_clusters/90_trial_" . $i . "/90_trial" . $i . "_decam_convolved*.fits";
print OUT "rm $name1 \n";
print OUT "rm $name2 \n";

}

close(OUT);