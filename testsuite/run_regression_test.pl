#!/usr/bin/env perl
#
# Copyright (C) 2005-2014 H. Appel, M. Marques, X. Andrade, D. Strubbe
#
# Originally GPL, dual-licensed under BerkeleyGW license by permission of these authors.
#
# Based on testsuite from Octopus project:
# $Id: oct-run_regression_test.pl 6053 2009-11-10 19:48:49Z dstrubbe $
# with many updates from Octopus since then.

#use warnings;
use Getopt::Std;
use File::Basename;
use File::Spec;
use Fcntl ':mode';
use Time::HiRes qw(gettimeofday tv_interval);
use Scalar::Util qw(looks_like_number);

sub usage {

  print <<EndOfUsage;

 Copyright (C) 2005-2014 H. Appel, M. Marques, X. Andrade, D. Strubbe

Usage: run_regression_test.pl [options]

    -n        dry-run
    -v        verbose
    -h        this usage
    -D        name of the directory where to look for the executables   
    -s        run everything serial
    -f        filename of testsuite [required]
    -p        preserve working directories
    -l        copy output log to current directory
    -m        run matches only (assumes there are work directories)

Exit codes:
    0         all tests passed
    1..253    number of test failures
    254       test skipped
    255       internal error

Report bugs to <BerkeleyGW\@civet.berkeley.edu>
EndOfUsage

  exit 0;
}


sub set_precision{
  my $p = $_[0];
  if($p ne "default"){
    $precnum = 1.0*$p;
  } else {
    $precnum = 0.0001
  }
}

# Check whether STDOUT is a terminal. If not, no ANSI sequences are written.
if(-t STDOUT) {
    $color_start{blue}="\033[34m";
    $color_end{blue}="\033[0m";
    $color_start{red}="\033[31m";
    $color_end{red}="\033[0m";
    $color_start{green}="\033[32m";
    $color_end{green}="\033[0m";
} else {
    $color_start{blue}="";
    $color_end{blue}="";
    $color_start{red}="";
    $color_end{red}="";
    $color_start{green}="";
    $color_end{green}="";
}

if (not @ARGV) { usage; }

getopts("nlvshD:f:pm");

# avoid warnings 'used only once: possible typo'
$useless = $opt_h;
$useless = $opt_l;

$test_succeeded = 1;
# set to 0 if anything fails

# Default values
use File::Temp qw/tempdir/;

# Handle options
$opt_h && usage;

my $exec_directory;
if($opt_D) {
 $exec_directory = $opt_D;
 if($exec_directory !~ /^\//){
  $exec_directory = get_env("PWD")."/$exec_directory";
 }
} else {
 $exec_directory = get_env("PWD")."/../../bin";
} 

if(length($opt_f) == 0) { 
    die255("ERROR: You must supply the name of a test file with the -f option.\n"); 
}

$aexec = get_env("EXEC");
if(get_env("BGW_TEST_MPI_NPROCS") ne "") {
    $global_np = get_env("BGW_TEST_MPI_NPROCS");
} elsif (get_env("BGW_TEST_NPROCS") ne "") {
    print STDERR "WARNING: BGW_TEST_NPROCS is deprecated, use BGW_TEST_MPI_NPROCS.\n";
    $global_np = get_env("BGW_TEST_NPROCS");
} else {
    $global_np = "";
}

$np = "serial";
$is_parallel = 0;

if(!$opt_s) {
# MPI stuff
    $mpiexec = get_env("MPIEXEC");
    $machinelist = get_env("MACHINELIST");
    if ("$mpiexec" eq "") { $mpiexec = `which mpiexec 2> /dev/null`; }
    chomp($mpiexec);

    if( "$mpiexec" eq "" ) { 
	print "No mpiexec found: running in serial.\n\n";
    } else {
# mpiexec without arguments (to check if it is available)
	$mpiexec_raw = $mpiexec;
	$mpiexec_raw =~ s/\ (.*)//;
	if ( ! -e "$mpiexec_raw" ) {
	    print "mpiexec command ($mpiexec_raw) does not exist: running in serial.\n\n";
	    $mpiexec = "";
	} elsif( ! -x "$mpiexec_raw" ) {
	    print "mpiexec command ($mpiexec_raw) is not executable: running in serial.\n\n";
	    $mpiexec = "";
	} else {
# default number of processors is 1
	    $np = 1;
	    $is_parallel = 1;
	}
    }
} else {
    $mpiexec = "";
}

$ENV{SHEXEC} = "";
if ( "$mpiexec" =~ /ibrun/ ) {
    $ENV{SHEXEC} = "MY_NSLOTS=1 $mpiexec $machinelist $aexec";
} elsif ( "$mpiexec" =~ /runjob/ ) {
    $ENV{SHEXEC} = "$mpiexec --np 1 $machinelist : $aexec";
} elsif ( "$mpiexec" ne "" ) {
    $ENV{SHEXEC} = "$mpiexec -n 1 $machinelist $aexec";
}

$output = "out";
$pipe_input = 0;
$cmd_line = 0;

# This variable counts the number of failed testcases.
$failures = 0;

$tempdirpath = get_env("TEMPDIRPATH");
if ("$tempdirpath" eq "") { $tempdirpath = '/tmp'; }
if (! -d $tempdirpath) { mkdir $tempdirpath; }

set_precision("default");

# testsuite
open(TESTSUITE, "<".$opt_f ) or die255("ERROR: cannot open testsuite file '$opt_f'.\n");

$return_value = 0;

while ($_ = <TESTSUITE>) {
  # remove trailing newline
  chomp;
  # remove trailing whitespace
  $_ =~ s/^\s+//;

  if($return_value != 0) { 
      print "\nSkipping subsequent steps due to nonzero exit code.\n\n";
      exit $failures;
  }

  # skip comments
  next if /^#/;

  if ( $_ =~ /^Test\s*:\s*(.*)\s*$/) {
    $test{"name"} = $1;
    print "$color_start{blue} ***** $test{\"name\"} ***** $color_end{blue} \n\n";
    print "Using test file  : $opt_f \n";

  } elsif ( $_ =~ /^Enabled\s*:\s*(.*)\s*$/) {
    %test = ();
    $enabled = $1;
    $enabled =~ s/^\s*//;
    $enabled =~ s/\s*$//;
    $test{"enabled"} = $enabled;

    if( $enabled eq "Yes" ) {
	$pwd = get_env("PWD");
	if (!$opt_m) {
	    $workdir = tempdir("$tempdirpath/BGW.XXXXXX");
	    chomp($workdir);

	    system ("rm -rf $workdir");
	    mkdir $workdir;

	    $scriptname = "$workdir/matches.sh";
	    open(SCRIPT, ">$scriptname") or die255("ERROR: could not create '$scriptname'.\n");
	    print SCRIPT "#\!/usr/bin/env bash\n\n";
	    print SCRIPT "perl $pwd/$0 -m -D $exec_directory -f $pwd/$opt_f\n";
	    close(SCRIPT);
	    chmod 0755, $scriptname;
  	} else {
	    $workdir = $pwd;
	}
        $ENV{TESTDIR} = File::Spec->rel2abs(dirname($opt_f));
        $ENV{WORKDIR} = File::Spec->rel2abs($workdir);

	print "Using workdir    : $workdir \n";
	if($opt_p) {
	    print "Workdir will be saved.\n";
	}
	
    } else {
	if ( $enabled eq "No") {
	    print "Test disabled: skipping test\n\n";
	    exit 254;
	} else {
	    die255("ERROR: Unknown option 'Enabled = $enabled' in testsuite file.\n\n");
	}
    }

  } else {

    if ( $enabled eq "") {
      die255("ERROR: Testsuite file must set Enabled tag before any others (except Test).\n\n");
    }
    if ( $_ =~ /^Executable\s*:\s*(.*)\s*$/) {
	$exe = "$exec_directory/$1";

	# beware: trailing space in executable name may cause this to fail.
	if( ! -x $exe) {
	  print "\nSkipping test: executable $exe not available.\n\n";
	  if (!$opt_p && !$opt_m && $test_succeeded) { system ("rm -rf $workdir"); }
	  if($failures == 0) {
	      exit 254;
	  } else {
	      exit $failures;
	      # if a previous step has failed, mark as failed not skipped
	  }
        }
    }

    elsif ( $_ =~ /^Processors\s*:\s*(.*)\s*$/) {
	$np = $1;
    }

    # for debugging purposes, to halt test after the part you are studying has run
    elsif ( $_ =~ "STOP TEST") {
	die255("\nSTOP TEST\n");
    }

    elsif ( $_ =~ /^TestGroups\s*:\s*(.*)\s*$/) {
	$groups = $1;
	print "Test belong to test groups: $groups\n";
    }

    elsif ( $_ =~ /^Banner\s*:\s*(.*)\s*$/) {
	$banner = $1;
	$len = length($banner);
	print "\n";
	print '+'; print '-'x($len+2); print "+\n";
	print "| $banner |\n";
	print '+'; print '-'x($len+2); print "+\n";
	print "\n";
    }

    elsif ( $_ =~ /^Copy\s*:\s*(\S*)\s*(\S*)\s*$/) {
	if($opt_m) { next; }

	$input_base = $1;
	$input_file = dirname($opt_f) . "/" . $input_base;
	$input_newname = $2;
	if($input_newname eq "") {
	  $input_newname = $input_base;
	}

	if( -f $input_file ) {
	  print "\n\nCopying file $input_file to $input_newname\n";
	  if(!$opt_n) {
	      system("cp $input_file $workdir/$input_newname");
	      # Ensure that the input file is writable so that it can
	      # be overwritten by the next test.
	      $mode = (stat "$workdir/$input_newname")[2];
	      chmod $mode|S_IWUSR, "$workdir/$input_newname";
	  }
	} else {
	  die255("ERROR: could not find input file: $input_file\n");
	}

    }

    elsif ( $_ =~ /^Command\s*:\s*(.*)\s*$/) {
	if($opt_m || $opt_n) { next; }

	$args = "$1";
	$command = "cd $workdir; $args";
	#BGQ needs runjob even for BerkeleyGW serial jobs
        if ( ("$mpiexec" =~ /runjob/) and ("$args" =~ /\.x/) ) {
	  $specify_np = "--np 1 : ";
	  $my_nslots = "";
	  $command = "cd $workdir; $my_nslots $mpiexec $specify_np $machinelist $aexec $args";
        }

	print "\nRunning command : $command\n";
	$return_value = system($command);

	if($return_value != 0) {
	    print "Command failed with exit code $return_value.\n";
	    $failures++;
	}
    }

    elsif ( $_ =~ /^Arguments\s*:\s*(.*)\s*$/) {
	if($opt_m || $opt_n) { next; }

	$args = "$1";
	$command = "cd $workdir; $aexec $exe $args";
	#BGQ needs runjob even for BerkeleyGW serial jobs
        if ( ("$mpiexec" =~ /runjob/) and ("$exe" =~ /\.x/) ) {
	  $specify_np = "--np 1 : ";
	  $my_nslots = "";
	  $command = "cd $workdir; $my_nslots $mpiexec $specify_np $machinelist $aexec $exe $args";
        }

	print "\nRunning command : $command\n";
	$return_value = system($command);

	if($return_value != 0) {
	    print "Command failed with exit code $return_value.\n";
	    $failures++;
	}
    }

    elsif ( $_ =~ /^Unpack\s*:\s*(.*)\s*$/) {
	if($opt_m || $opt_n) { next; }

	$fname = $1;
	if (("$fname" =~ /\.tar$/)) {
		# Do not unpack
		$tar_opts = "xf";
	}
	elsif (("$fname" =~ /\.gz$/) or ("$fname" =~ /\.tgz$/)) {
		# Unpack with gzip
		$tar_opts = "xzf";
	}
	elsif (("$fname" =~ /\.bz2$/) or ("$fname" =~ /\.bzip2$/)) {
		# Unpack with bzip2
		$tar_opts = "xjf";
	}
	elsif (("$fname" =~ /\.xz$/)) {
		# Unpack with xz
		$tar_opts = "xJf";
	}
	else {
		print STDERR "\nWARNING: File extension not recognized for unpacking: $fname\n";
		$tar_opts = "xzf"
	}
	$command = "tar $tar_opts " . dirname($opt_f) . "/$fname -C $workdir";

	print "\nUnpacking archive : $command\n";
	$return_value = system($command);

	if($return_value != 0) {
	    print "Unpack failed with exit code $return_value.\n";
	    $failures++;
	}
    }

    elsif ( $_ =~ /^Output\s*:\s*(.*)\s*$/) {
	$output = $1;
    }

    elsif ( $_ =~ /^Input\s*:\s*(\S*)\s*(\S*)\s*(\S*)\s*$/) {
	$input_base = $1;
	$input_none = ($input_base eq "NONE");
	$nocopy = (($2 eq "NOCOPY") or ($3 eq "NOCOPY"));
	if ( !$input_none ){
            if ( $nocopy ){
	      $input_file = $workdir . "/" . $input_base;
            } else {
	      $input_file = dirname($opt_f) . "/" . $input_base;
            }
	} else {
	    $input_file = "";
	}

	if($opt_m) {
	    print "\n\nFor input file : $input_base\n\n";
	    $return_value = 0;
	} else {
	  $input_newname = $2;
	  $pipe_input = ($input_newname eq "PIPE");
	  $cmd_line = ($input_newname eq "CMDLINE");
	  if($input_newname eq "") {
	      $input_newname = $input_base;
	  }

	  if( !$input_none ) {
	      if( -f $input_file ) {
                  if ( !$nocopy ) {
	              if( !$pipe_input && !$cmd_line) {
	                  print "\n\nUsing input file : $input_file \n";
	                  system("cp $input_file $workdir/$input_newname");
	                  # Ensure that the input file is writable so that it can
	                  # be overwritten by the next test.
	                  $mode = (stat "$workdir/$input_newname")[2];
	                  chmod $mode|S_IWUSR, "$workdir/$input_newname";
	              } else {
	                  system("cp $input_file $workdir/");
	              }
                  }
	      } else {
	          if(!$opt_n) {
	              die255("ERROR: could not find input file: $input_file\n");
	          }
	      }
	  }

	  print "\nStarting test run ...\n";

	  if ($pipe_input) {
	    $input_text = "< $input_base";
	  } elsif($cmd_line) {
	    $input_text = "$input_base";
	  } else {
	    $input_text = ""; 
	  }

	  # serial or MPI run?
	  if ( $is_parallel && $np ne "serial" && !$opt_s) {
	    if("$global_np" ne "") {
	      $np = $global_np;
	    }
	    if ( "$mpiexec" =~ /ibrun/ ) { # used by SGE parallel environment
	      $specify_np = "";
	      $my_nslots = "MY_NSLOTS=$np";
	    } elsif ("$mpiexec" =~ /runjob/) { # used by BlueGene
	      $specify_np = "--np $np : ";
	      $my_nslots = "";
	    } elsif ("$mpiexec" =~ /poe/) { # used by IBM PE
	      $specify_np = "";
	      $my_nslots = "MP_PROCS=$np";
	    } else { # for mpirun and Cray's aprun
	      $specify_np = "-n $np";
	      $my_nslots = "";
	    }
	    $command_line = "cd $workdir; $my_nslots $mpiexec $specify_np $machinelist $aexec $exe $input_text > $output";
	  } else {
	    #BGQ needs runjob even for BerkeleyGW serial jobs
	    if ( ("$mpiexec" =~ /runjob/) and ("$exe" =~ /\.x/) ) {
	      $specify_np = "--np 1 ";
	      $my_nslots = "";
	      $command_line = "cd $workdir; $my_nslots $mpiexec $specify_np $machinelist : $aexec $exe $input_text > $output";
            } else {
	      $command_line = "cd $workdir; $aexec $exe $input_text > $output";
            }
	  }

	  print "Executing: " . $command_line . "\n";

	  if(!$opt_n) {
	    $test_start = [gettimeofday];
	    $return_value = system("$command_line");
	    $test_end   = [gettimeofday];
	    
	    $elapsed = tv_interval($test_start, $test_end);
	    printf("\tElapsed time: %8.1f s\n\n", $elapsed);

            if($return_value == 0) { 
		print "Finished test run.\n\n"; 
		printf "%-40s%s", " Execution", ": \t [ $color_start{green}  OK  $color_end{green} ] \n"; 

	    } else { 
		print "\n\nTest run failed with exit code $return_value.\n\n"; 
		printf "%-40s%s", " Execution", ": \t [ $color_start{red} FAIL $color_end{red} ] \n\n"; 
		$failures++; 
		$test_succeeded = 0; 
	    } 
	    $test{"run"} = 1;
	  }
	  if ($opt_l && !$opt_m && !$opt_n) {
	    system ("cat $workdir/$output >> out.log");
	  }
	}

	$pipe_input = 0;
	$cmd_line = 0;
    }

    elsif ( $_ =~ /^Precision\s*:\s*(.*)\s*$/) {
	set_precision($1);
    }

    elsif ( $_ =~ /^match/ ) {
	if (!$opt_n && $return_value == 0) {
	    if(run_match_new($_)){
		printf "%-40s%s", "$name", ": \t [ $color_start{green}  OK  $color_end{green} ] \t (Calculated value = $value)\n";
		if ($opt_v) { print_hline(); }
	    } else {
		printf "%-40s%s", "$name", ": \t [ $color_start{red} FAIL $color_end{red} ] \n";
		print_hline();
		$test_succeeded = 0;
		$failures++;
	    }
	}
    }

    elsif (length($_) != 0) {
	die255("ERROR: Unknown command '$_'\n");
    }

  }

}

if (!$opt_p && !$opt_m && $test_succeeded) { system ("rm -rf $workdir"); }

print "\n";
close(TESTSUITE);

exit $failures;


sub run_match_new {
  die255("ERROR: Have to run before matching.\n") if !$test{"run"} && !opt_m;

  # parse match line
  my ($line, $match, $pre_command, $ref_value, $off);
  $line = $_[0];
  $line =~ s/\\;/_COLUMN_/g;
  ($match, $name, $pre_command, $ref_value) = split(/;/, $line);
  $pre_command =~ s/_COLUMN_/;/g;
  $ref_value =~ s/^\s*//;
  $ref_value =~ s/\s*$//;

  # parse command
  $pre_command =~ /\s*(\w+)\s*\((.*)\)/;

  my $func = $1;
  my $params = $2;

  # parse parameters
  $params =~ s/\\,/_COMMA_/g;
  my @par = split(/,/, $params);
  for($params=0; $params <= $#par; $params++){
    $par[$params] =~ s/_COMMA_/,/g;
    $par[$params] =~ s/^\s*//;
    $par[$params] =~ s/\s*$//;
  }

  if($func eq "SHELL"){ # function SHELL(shell code)
    check_num_args(1, 1, $#par, $func);
    $pre_command = $par[0];

  }elsif($func eq "LINE") { # function LINE(filename, line, field)
    check_num_args(3, 3, $#par, $func);
    if($par[1] < 0) { # negative number means from end of file
      $line_num = "`wc -l $par[0] | awk '{print \$1}'`";
      $pre_command = "awk -v n=$line_num '(NR==n+$par[1]+1) {printf \$$par[2]}' $par[0]";
    } else {
      $pre_command = "awk '(NR==$par[1]) {printf \$$par[2]}' $par[0]";
    }

  }elsif($func eq "GREP") { # function GREP(filename, 're', field <, offset>)
    check_num_args(3, 4, $#par, $func);
    if($#par == 3) {
        $off = $par[3];
    } else {
	$off = 0;
    }
    # -a means even if the file is considered binary due to a stray funny character, it will work
    $pre_command = "grep -a -A$off $par[1] $par[0]";
    $pre_command .= " | awk '(NR==$off+1) {printf \$$par[2]}'";
    # if there are multiple occurrences found by grep, we will only be taking the first one via awk

  }elsif($func eq "SIZE") { # function SIZE(filename)
    check_num_args(1, 1, $#par, $func);
    $pre_command = "ls -lt $par[0] | awk '{printf \$5}'";

  }else{ # error
    printf STDERR "Unknown command '$func'\n";
    return 0;
  }

  # 'set -e; set -o pipefail' (bash3 only) would make the whole pipe series gives an error if any step does; 
  # otherwise the error comes only if the last step failed.
  $value = qx(cd $workdir && $pre_command);
  # Perl gives error code shifted, for some reason.
  $exit_code = $? >> 8;
  if($exit_code) {
      print STDERR "Match command failed: $pre_command\n";
      return 0;
  }

  # extract non-numeric string (including possibility of NaN)
  if($value =~ /([0-9\-+.eEdDnNaA]+)/) {
      $value = $1;
      chomp $value;
  } else {
      $value = "";
  }

  if(length($value) == 0) {
      print STDERR "Match command returned nothing: $pre_command\n";
      return 0;
  }

  if(!looks_like_number($value)) {
      print STDERR "Match command returned non-numeric value '$value': $pre_command\n";
      return 0;
  }

  # at this point, we know that the command was successful, and returned a number. 
  $success = (abs(($value)-($ref_value)) <= $precnum);

  if(!$success || $opt_v) {
    print_hline();
    print "Match".$name.":\n\n";
    print "   Calculated value : ".$value."\n";
    print "   Reference value  : ".$ref_value."\n";
    print "   Difference       : ".abs($ref_value - $value)."\n";
    print "   Tolerance        : ".$precnum."\n\n";
  }

  return $success;
}

sub print_hline {
  print "\n-----------------------------------------\n\n";
}

# return value of environment variable (specified by string argument), or "" if not set
sub get_env {
    if(exists($ENV{$_[0]})) {
        return $ENV{$_[0]};
    } else {
        return "";
    }
}

# args: min num args, max num args, args given, function name
sub check_num_args {
    my $min_num_args   = $_[0];
    my $max_num_args   = $_[1];
    my $given_num_args = $_[2]+1;
    my $func_name      = $_[3];

    if($given_num_args < $min_num_args) {
        die255("$func_name given $given_num_args argument(s) but needs at least $min_num_args.\n");
    }
    if($given_num_args > $max_num_args) {
        die255("$func_name given $given_num_args argument(s) but can take no more than $max_num_args.\n");
    }
}

sub die255 {
    print STDERR $_[0];
    exit 255;
}
