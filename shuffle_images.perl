#!/usr/dist/bin/perl -w
use List::Util qw(shuffle);

for ($i=0; $i<168; $i++) {

$inlist[$i]=$i;

}
$N=168;

@outlist = shuffle(@inlist);
open(OUT,">commands");
open(OUT2,">mappings");

for ($j=0; $j<84;$j++) {

$inname = "trial" . $j . "_LSST_convolved_noise.fits";
$outname = "grid" . $outlist[$j] . ".fits";

print OUT "cp $inname $outname \n";
print OUT2 "$inname $outname \n";

}

for ($j=84; $j<168; $j++) {
	$k= $j-84;
$inname = "90_trial" . $k . "_LSST_convolved_noise.fits";
$outname = "grid" . $outlist[$j] . ".fits";

print OUT "cp $inname $outname \n";
print OUT2 "$inname $outname \n";
}

close(OUT);
close(OUT2);
