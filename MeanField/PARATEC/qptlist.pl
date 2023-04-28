#!/usr/bin/env perl -w
# qptlist.pl
#
# Usage: qptlist.pl shift_x shift_y shift_z < SCF_KPOINTS > qpoints
#
# Extracts a formatted list of q-points
# for use in the GW Epsilon code
#
# David Prendergast, UCB, April 2007

if( @ARGV != 3 ) {
  print "usage qptlist.pl shift_x shift_y shift_z < SCF_KPOINTS > qpoints\n";
  exit;
}

$shift_x = shift @ARGV;
$shift_y = shift @ARGV;
$shift_z = shift @ARGV;

$nkpt=<>;
chomp($nkpt);
print "number_qpoints $nkpt\n";
print "begin qpoints\n";
for($i=0; $i<$nkpt; $i++) {
  $line=<>;
  chomp($line);
  @data=split ' ', $line;
  $kz = pop @data;
  $ky = pop @data;
  $kx = pop @data;
  if( $kx==0.0 && $ky==0.0 && $kz==0.0 ) {
    $lenx=length($kx); @data=split /\./, $kx; $decx=length($data[-1]);
    $leny=length($ky); @data=split /\./, $ky; $decy=length($data[-1]);
    $lenz=length($kz); @data=split /\./, $kz; $decz=length($data[-1]);
    $kx+=$shift_x;
    $ky+=$shift_y;
    $kz+=$shift_z;
    $pstr="";
    $pstr.="  %" . $lenx . "." . $decx . "f";
    $pstr.="  %" . $leny . "." . $decy . "f";
    $pstr.="  %" . $lenz . "." . $decz . "f";
    $pstr.="  1.0  1\n";
    printf "$pstr", $kx, $ky, $kz;
  } else {
    print "  $kx  $ky  $kz  1.0  0\n";
  }
}
print "end\n";
exit;
