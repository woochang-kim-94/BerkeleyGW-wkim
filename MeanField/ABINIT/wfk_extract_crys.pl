#!/usr/bin/perl
#===============================================================================
#
# Perl routine
#
# wfk_extract_crys.pl            Written by  Tonatiuh Rangel         Last Modified 07/30/2014    
#
# Perl routine to parse crystal informations 
# from standard output of ABINIT (out file).
#
# It creates a file named "crystal.dat" which is used later by abi2bgw
#
# Developers:
#   Tonatiuh Rangel, Berkeley CA, trangel@lbl.gov
#   Felipe Jornada, Berkeley CA
#
#===============================================================================




#  type crystal
#    real(DP) :: celvol !< cell volume in real space (a.u.)
#    real(DP) :: recvol !< cell volume in reciprocal space (a.u.)
#    real(DP) :: alat !< lattice constant in real space (a.u.)
#    real(DP) :: blat !< lattice constant in reciprocal space (a.u.)
#    real(DP) :: avec(3,3) !< lattice vectors in real space (alat)
#    real(DP) :: bvec(3,3) !< lattice vectors in reciprocal space (blat)
#    real(DP) :: adot(3,3) !< metric tensor in real space (a.u.)
#    real(DP) :: bdot(3,3) !< metric tensor in reciprocal space (a.u.)
#    integer :: nat !< number of atoms
#    integer, pointer :: atyp(:) !< atomic species, atyp(1:nat)
#    real(DP), pointer :: apos(:,:) !< atomic positions, apos(1:3,1:nat) (alat)
#  end type crystal


if( @ARGV <  1 ) {
    ARGVerror();
}

$file=$ARGV[0];
$crystal_file="crystal.dat";

read_file($file);

#Extract BGW variable atyp 
extract_atyp();

#Convert rprimd+acell to alat+avec
get_alat_avec();
#Convert apos from [bohrs] to [alat]
convert_apos_to_alat();

print_crystal_file($crystal_file);

################
sub get_alat_avec(){
  $alat=$acell[1];
  for($ii=1;$ii<=3;$ii++){
    $avec[$ii][1]=$rprim[$ii][1]*$acell[$ii];
    $avec[$ii][2]=$rprim[$ii][2]*$acell[$ii];
    $avec[$ii][3]=$rprim[$ii][3]*$acell[$ii];
  }
# Convert to alat:
  for($ii=1;$ii<=3;$ii++){
    for($jj=1;$jj<=3;$jj++){
      $avec[$ii][$jj]=$avec[$ii][$jj]/$alat;
    }
  }
}
################


################
sub convert_apos_to_alat(){
  for($iat=1;$iat<=$natom;$iat++){
    for($ii=1;$ii<=3;$ii++){
      $apos[$iat][$ii]=$apos[$iat][$ii]/$alat;
    }
  }
}
################

################
sub extract_atyp(){
################
  for($iat=1;$iat<=$natom;$iat++){
    $ityp=$typat[$iat];
    #print "ityp $ityp, znucl $znucl[$ityp]\n";
    $zz=$znucl[$ityp];
    $atyp[$iat]=int($zz);
  }
}

################
sub print_crystal_file(){
################
  $file=$_[0];
  
  printf("**************************************\n");
  open OUT, ">$file" or die "Opening '$file': $!\n";
  $iline=0;
#  type crystal
#    real(DP) :: celvol !< cell volume in real space (a.u.)
  printf( OUT "celvol %20.12f\n",$celvol);
#    real(DP) :: recvol !< cell volume in reciprocal space (a.u.)
#    real(DP) :: alat !< lattice constant in real space (a.u.)
  printf( OUT "alat %20.12f\n",$alat);
#    real(DP) :: blat !< lattice constant in reciprocal space (a.u.)
#    real(DP) :: avec(3,3) !< lattice vectors in real space (alat)
  printf( OUT "avec\n");
  for( $ii=1;$ii<=3;$ii++){
    printf( OUT " %20.12f %20.12f %20.12f\n",$avec[$ii][1],$avec[$ii][2],$avec[$ii][3]);
  }
#    real(DP) :: bvec(3,3) !< lattice vectors in reciprocal space (blat)
#    real(DP) :: adot(3,3) !< metric tensor in real space (a.u.)
#    real(DP) :: bdot(3,3) !< metric tensor in reciprocal space (a.u.)
# All of the above can be obtained from the lattice parameters:
#  printf( OUT "acell %20.12f %20.12f %20.12f\n",$acell[1],$acell[2],$acell[3]);
#    integer :: nat !< number of atoms
  printf( OUT "nat %6i\n",$natom);
#    integer, pointer :: atyp(:) !< atomic species, atyp(1:nat)
  printf( OUT "atyp\n");
  for ($iat=1;$iat<=$natom;$iat++){
    printf( OUT " %6i \n",$atyp[$iat]);
  }
#    real(DP), pointer :: apos(:,:) !< atomic positions, apos(1:3,1:nat) (alat)
  printf(OUT "apos\n");
  for ($iat=1;$iat<=$natom;$iat++){
    printf( OUT " %20.12f %20.12f %20.12f\n",$apos[$iat][1],$apos[$iat][2],$apos[$iat][3]);
  }
#  end type crystal



  close OUT;
}


################
sub read_file(){
################
  $file=$_[0];
  
  printf("**************************************\n");
  printf("Lines read from OUT file:\n");
  open TEMP, "<$file" or die "Opening '$file': $!\n";
  $iline=0;

  $lntypat=0; $lnatom=0; $ipseudo=0; $iapos=0;
  $state=0;  
  while ($line = <TEMP>){
   $iline=$iline+1;

   if( $state==0){
     if( $line=~/^\s+Unit\scell\svolume\sucvol/){
       print $line;
       get_tokens();
       $celvol=$tokens[4]; 
     }if( $line=~/^\s+acell/){
       print $line;
       get_tokens();
       $acell[1]=$tokens[1]; 
       $acell[2]=$tokens[2]; 
       $acell[3]=$tokens[3]; 
     }if( $line=~/^\s+natom/){
       print $line;
       get_tokens();
       $natom=$tokens[1];
       $lnatom=1;
     }elsif( $line=~/rprim/){
       print $line;
       get_tokens();
       $rprim[1][1]=$tokens[1];
       $rprim[1][2]=$tokens[2];
       $rprim[1][3]=$tokens[3];
       $state=10;
     }elsif( $line=~/^\s+ntypat/){
       print $line;
       get_tokens();
       $ntypat=$tokens[1];
       $lntypat=1;
     }elsif( $line=~/^\s+typat/){
       print $line;
       if( $lnatom == 0 ){ 
         printf("Error at line $iline: typat found before natom\n");
         exit(1);
       }
       get_tokens();
       $itypat=0;
       for($itoken=1;$itoken<$ntokens;$itoken++){
         $itypat++;
         $typat[$itypat]=$tokens[$itoken];
       }
       if($itypat<$natom){
         $state=20;
       }
     }elsif( $line=~/pspini:\satom\stype/){
       print $line;
       get_tokens();
       $ipseudo++;
       $ipseudo_=$tokens[4];
       if($ipseudo ne $ipseudo_){
         printf("Error at line $iline: pseudofile $ipseudo expected, found: $ipseudo_\n");
         exit(1);
       }
       $state=30;
     }elsif( $line=~/^\s+xcart/){
       print $line;
       get_tokens();
       $iapos++;
       $apos[$iapos][1]=$tokens[1];
       $apos[$iapos][2]=$tokens[2];
       $apos[$iapos][3]=$tokens[3];
       if( $iapos < $natom){
          $state=40;
       }
     }
   }elsif( $state==10 ){
     print $line;
     get_tokens();
     $rprim[2][1]=$tokens[0];
     $rprim[2][2]=$tokens[1];
     $rprim[2][3]=$tokens[2];
     $state=11;
   }elsif( $state==11 ){
     print $line;
     get_tokens();
     $rprim[3][1]=$tokens[0];
     $rprim[3][2]=$tokens[1];
     $rprim[3][3]=$tokens[2];
     $state=0;
   }elsif( $state==20 ){
     print $line;
     get_tokens();
     for($itoken=0;$itoken<$ntokens;$itoken++){
       $itypat++;
       $typat[$itypat]=$tokens[$itoken];
     }
     if($itypat==$natom){
       $state=0;
     }elsif($itypat>$natom){
       printf("Error at line $iline: itypat=$itypat > natom=$natom\n");
       exit(1);
     }
   }elsif( $state==30 ){
     if( $line=~/znucl/ ){
       print $line;
       get_tokens();
       $znucl[$ipseudo]=$tokens[1];
       $state=0;
     }
   }elsif( $state==40 ){
     print $line;
     get_tokens();
     $iapos++;
     $apos[$iapos][1]=$tokens[0];
     $apos[$iapos][2]=$tokens[1];
     $apos[$iapos][3]=$tokens[2];
     if( $iapos == $natom){
        $state=0;
     }
   }    
  }
  
  close TEMP;
}
################
sub get_tokens(){
 chomp($line);
 @tokens=split(" ",$line);
 $ntokens=@tokens;
}
################


###############
sub ARGVerror {
###############
    chop ( $case = qx /  echo \$PWD \| awk -F \/ '{print\$NF} '  /  );   #directory at which the user is working
    system("clear");
    print "\e[0;31m*******************************************************\e[00m \n";
    print "The scripts needs to be run as follows:\n";
    print "\n";
    print "wfk_extract_crys.pl \e[0;34m  outfile \e[0m\n";
    print "\n";
    print "outfile is an ABINIT standard file: out[A-Z]* \n";
    print "\e[0;31m******************************************************\e[00m \n";
    exit(1);
}

# Declare the subroutines


# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
# Left trim function to remove leading whitespace
sub ltrim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	return $string;
}
# Right trim function to remove trailing whitespace
sub rtrim($)
{
	my $string = shift;
	$string =~ s/\s+$//;
	return $string;
}
