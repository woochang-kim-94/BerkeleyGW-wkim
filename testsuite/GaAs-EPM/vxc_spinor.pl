#!/usr/bin/env perl

# double size of vxc.dat file

open(vin, "vxc_in");
open(OUT, ">vxc_out");

$nbnd = 8;
$count=0;

foreach $line (<vin>) {
  $count++;
  chomp($line);
  if($count%($nbnd+1) == 1) {
    ($null, $kx, $ky, $kz, $ndg, $nodg) = split(/\s+/, $line);
    $ndg=2*$ndg;
    $n0dg=2*$n0dg;
    printf OUT "   %.10f   %.10f   %.10f           %2i           %2i\n",$kx,$ky,$kz,$ndg,$nodg;
  }
  else {
    ($null, $ns, $nb, $rev, $imv) = split(/\s+/, $line);
    $nb1 = 2*$nb - 1;
    $nb2 = 2*$nb;
    printf OUT "  %i         %2i        %.9f           %.9f\n",$ns,$nb1,$rev,$imv;
    printf OUT "  %i         %2i        %.9f           %.9f\n",$ns,$nb2,$rev,$imv;
  }
}

close(vin);
close(OUT);
