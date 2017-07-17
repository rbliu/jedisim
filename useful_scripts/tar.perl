#!/usr/dist/bin/perl -w

#to create a shell file to compress all the stamp_*/ and distorted_*/ folders and remove the original folders

open(OUT,">tar_stamps");

for ($i=0; $i<80; $i++){

$name1 = "/Volumes/rliu/jedisim_clusters/trial_" . $i . "/stamp.tar.gz";
$name2 = "/Volumes/rliu/jedisim_clusters/trial_" . $i . "/stamp_*";
$name3 = "/Volumes/rliu/jedisim_clusters/trial_" . $i . "/distorted.tar.gz";
$name4 = "/Volumes/rliu/jedisim_clusters/trial_" . $i . "/distorted_*";

$name5 = "/Volumes/rliu/jedisim_clusters/90_trial_" . $i . "/stamp.tar.gz";
$name6 = "/Volumes/rliu/jedisim_clusters/90_trial_" . $i . "/stamp_*";
$name7 = "/Volumes/rliu/jedisim_clusters/90_trial_" . $i . "/distorted.tar.gz";
$name8 = "/Volumes/rliu/jedisim_clusters/90_trial_" . $i . "/distorted_*";

print OUT "tar -zcvf $name1 $name2 \n";
print OUT "tar -zcvf $name3 $name4 \n";
print OUT "rm -rf $name2 \n";
print OUT "rm -rf $name4 \n";

print OUT "tar -zcvf $name5 $name6 \n";
print OUT "tar -zcvf $name7 $name8 \n";
print OUT "rm -rf $name6 \n";
print OUT "rm -rf $name8 \n";

}

close(OUT);