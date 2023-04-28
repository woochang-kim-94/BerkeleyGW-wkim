#!/usr/bin/env perl

# Query a buildbot using the Octopus/APE/BerkeleyGW testsuite infrastructure for
# values obtained by each buildslave for a particular test and match.
# Parses the HTML status pages.
# Tested with BuildBot 0.7.12 and 0.8.5.

# modify these to choose the input file and match you want to search for
$inputfile = "GaAs-EPM/sig_GPP.inp";
$match = "n2, k1 CH + Static Remainder";

# options specifying setup for BerkeleyGW; will only work from civet.berkeley.edu due to https
$bbpath = "http://localhost:8080";
$shell_num = 5;

print "URL: $bbpath\n";
print "Input file: $inputfile\n";
print "Match: $match\n\n";

# get list of latest builds
if(-e "one_box_per_builder") { system ("rm one_box_per_builder"); }
system ("wget -nv $bbpath/one_box_per_builder");

if(@ARGV > 0) {
    open(ONEBOX, ">>one_box_per_builder") or die "cannot open one_box_per_builder\n";
    print ONEBOX "\n";
    foreach(@ARGV) {
	print ONEBOX "LOCAL FILENAME: $_\n";
    }
    close(ONEBOX);
}

open(ONEBOX, "<one_box_per_builder") or die "cannot open one_box_per_builder\n";

$total = 0.0;
$counts = 0;
$max = -inf;
$maxname = "";
$min = inf;
$minname = "";

while ($_ = <ONEBOX>) {
# BB 0.7.12
# <td align="center" class="LastBuild box success"><a href="builders/lascar_x86_64_gfortran_cl_intel/builds/139">10187</a><br />build<br />successful</td>
# BB 0.8.5
#<a href="builders/mauchly_x86_64_intel_openmp/builds/80">10898</a>
    if ( $_ =~ /<a href="builders\/(.*)\/builds\/(.*)">(.*)<\/a>/) {
	$builder = $1;
	$build_num = $2;
	$svn_rev = $3;
	# rebuild the URL
	$url = "builders/$builder/builds/$build_num";
	print "\nBuilder: $builder, at svn revision $svn_rev\n";

	# remove old file, or new ones will be named 'stdio.2' etc.
	if(-e "stdio") { system ("rm stdio"); }
	system ("wget -nv $bbpath/$url/steps/shell_$shell_num/logs/stdio");

	$name = $builder;
	open(TESTLOG, "<stdio") or print "cannot open test log\n";
    } elsif ( $_ =~ /LOCAL FILENAME: (.*)/) {
	print "\n$_\n";
	$name = $1;
	open(TESTLOG, "<$1");
    } else {
	next;
    }
    $match_found = 0;
    while ($_ = <TESTLOG>) {
	# do not use ~= / .. / here or $filename needs to have special characters escaped
	if(index($_, $inputfile) != -1) {
	    while ($_ = <TESTLOG>) {
		if(index($_, $match) != -1) {
		    if($_ =~ /\(Calculated value = (.*)\)/) {  # match OK
			print $_;
			$value = $1;
		    } else {  # match FAIL
			while ($_ = <TESTLOG>) {
			    if(index($_, $match) == -1) {
				if($_ !~ /^$/) { # print if not blank
				    print $_;
				    if($_ =~ /Calculated value : (.*)/) {
					$value = $1;
				    }
				}
			    } else {
				last;
			    }
			}
		    }
		    $total += $value;
		    $counts += 1;
		    if($value < $min) {
			$minname = $name;
			$min = $value;
		    }
		    if($value > $max) {
			$maxname = $name;
			$max = $value;
		    }
		    $match_found = 1;
		}
		if($_ =~ /Using input file/) { last; }
	    }
	}
    }
    close(TESTLOG);
    if($match_found == 0) {
	print "Match not found.\n";
    }
    # why not? builder down, svn or compilation failed, not in the right category of builders, etc.
}

if($counts == 0) {
    print "No matches found.\n";
} else {
    print "\n\n=== SUMMARY ===\n";
    print "Based on $counts matches found.\n";
    print "Minimum   = $min\n";
    print "    ($minname)\n";
    print "Maximum   = $max\n";
    print "    ($maxname)\n";
    print "Average   = " . ($total / $counts) . "\n\n";
    print "Center    = " . ($max + $min)/2 . "\n";
    printf "Precision = %e\n", ($max - $min)/2;
}

close(ONEBOX);
