#!/usr/bin/env perl -w
# kptlist.pl
#
# Usage: kptlist.pl < SCF_KPOINTS > kpoints
#
# Extracts a formatted list of k-points
# for use in the GW Sigma code
#
# David Prendergast, UCB, April 2007

$nkpt=<>;
chomp($nkpt);
print "number_kpoints $nkpt\n";
print "begin kpoints\n";
for($i=0; $i<$nkpt; $i++) {
  $line=<>;
  chomp($line);
  @data=split ' ', $line;
  $kz = pop @data;
  $ky = pop @data;
  $kx = pop @data;
  print "  $kx  $ky  $kz  1.0\n";
}
print "end\n";
exit;
