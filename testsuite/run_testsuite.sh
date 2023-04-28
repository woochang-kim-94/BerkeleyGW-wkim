#!/bin/bash -l
#
# Copyright (C) 2005-2009 Heiko Appel, David Strubbe
#
# Originally GPL, dual-licensed under BerkeleyGW license by permission of these authors.
#
# Based on Octopus file:
# $Id: oct-run_testsuite 2730 2007-02-28 20:57:45Z xavier $


# Paths.
#prefix=/global/home/users/dstrubbe/octopus
#pkgdatadir=${prefix}/share/octopus
#testsuite="$pkgdatadir/testsuite"
testsuite="."

# Failure reporting.
failed_tests=0
skipped_tests=0
passed_tests=0
total_tests=0
failure_report=""


# Usage.
function usage() {
    cat <<EOF

 Copyright (C) 2005-2009 by Heiko Appel, David Strubbe

Usage: run_testsuite [options]
    
     -h            this message
     -n            dry-run mode (show what would be executed)
     -g LIST       comma-separated list of test groups to run
     -q            query testfiles for the given groups (no tests are run)
     -b DIR        directory where to look for the binary files [default: ../bin]
     -d DIR        directory where to look for the testsuite
     -l            local run
     -p PREFIX     installation prefix [default: /usr]
     -w            without MPI tests (do not run parallel tests)
     -c            delete all .log files and work directories after the run

Report bugs to <BerkeleyGW@civet.berkeley.edu>.
EOF
 exit 0;
}


# Find tests to run. Takes as argument the string passed to the
# -g command line option.
function find_tests() {
    groups="$1"
    if [ -n "$groups" ]; then
	groups_pat=`echo "$groups" | sed 's/ *, */|/'`
	groups_pat="^TestGroups *:.*($groups_pat)"
    else # No groups given defaults to all tests.
	groups_pat="."
    fi
    testfiles=`find $testsuite -name "*.test" | \
	xargs grep -l -E "$groups_pat" | sort -u`
    echo "$testfiles"
}


# Run all tests. Takes as first argument a list of testfile names.
function run_testsuite() {

tests="$1"

# Check for 'preserve working directories' flag.
if [ "$cleanup" == "yes" ]; then
    preserve_opt=""
else
    preserve_opt="-p"
fi

for y in $tests; do
    total_tests=$((total_tests + 1))
    ybase=`basename $y`
    xybase=${ybase%.test}
#    if [ "${local_run}" == "yes" ]; then
#	bin_directory=`pwd`/../src/main
#	runtest_directory=$testsuite
#    else
#	bin_directory=$prefix/bin
#	runtest_directory=$bin_directory
#    fi
#    if [ -n "${exec_suffix}" ]; then
#	suffix_opt="-s ${exec_suffix}"
#    else
#	suffix_opt=""
#    fi
    if [ -d "${bin_directory}" ]; then
	./run_regression_test.pl -l $opt_n \
	    $preserve_opt -D $bin_directory -f $y $serial_opt | \
	    tee ${xybase}.log
	# Look for failed/passed/skipped testcases and add
	# to failure summary.
	failures=${PIPESTATUS[0]}
	if [ "${failures}" -eq 0 ]; then
	    passed_tests=$((passed_tests + 1))
	elif [ ${failures} -eq 254 ]; then
	    skipped_tests=$((skipped_tests + 1))
	else
	    failed_tests=$((failed_tests + 1))
	    test_name=`echo $y | sed "s|$testsuite/||"`
	    name_length=`echo -n $test_name | wc -m`
	    space_length=$((50 - name_length))
	    spaces=`for((i = 1; i <= space_length; i++)); do echo -n ' '; done`
	    failure_report="$failure_report    $test_name$spaces$failures\n"
	fi
    fi
    [ -e out.log ] && mv out.log ${xybase}.out.log
done
}


# Parse command line.

# Some default settings.
query="no"
cleanup="yes"
test_groups=""
dry_run="no"
bin_directory="../bin"

# default is yes
# environment variable can set to no
if [ "x$SAVETESTDIRS" == "xyes" ]; then
    cleanup="no";
fi

while getopts "hnlg:qwd:b:" opt ; do
    case "$opt" in
        h) usage;;
        n) dry_run="yes";; 
#        p) prefix="$OPTARG"; testsuite=$prefix/share/octopus/testsuite;;
#        e) exec_suffix="$OPTARG";;
        l) local_run="yes";;
	g) test_groups="$OPTARG";;
	q) query="yes";;
        w) serial_opt="-s";;
	d) directory="$OPTARG";;
	b) bin_directory="$OPTARG";;
#	c) cleanup="yes";;
        ?) echo "Error parsing arguments"; exit 1;;
    esac
done
shift $[ OPTIND - 1 ]


# Find testfiles.
if [ -n "$directory" ]; then
    testsuite="$directory"
else
    [ "$local_run" == "yes" ] && testsuite=$(pwd)
fi

testfiles=`find_tests "$test_groups"`


# Query mode? If so, list files and exit.
if [ "$query" == "yes" ]; then
    echo "Testfiles for groups $test_groups:"
    echo ""
    for f in $testfiles; do
	echo ${f##$testsuite/}
    done
    exit 0
fi


# No testfiles found, abort.
if [ -z "$testfiles" ]; then
    echo "No testfiles for group(s) $test_groups found."

# Otherwise, start the whole machinery.
else
    # Get epoch seconds at testsuite start.
    testsuite_start=$(date +%s) 

    if [ "$dry_run" == "yes" ]; then
	opt_n="-n"
    else
	opt_n=""
    fi

    # Run testsuite.
    run_testsuite "$testfiles"

    # Failure reporting to STDOUT.
    echo -e "    Passed:  $passed_tests / $total_tests"
    echo -e "    Skipped: $skipped_tests / $total_tests"
    if [ $failed_tests -gt 0 ]; then
	echo "    Failed:  $failed_tests / $total_tests"
	echo
	echo "    testfile                                          # failed testcases"
	echo "    --------------------------------------------------------------------"
	echo -e "$failure_report"
    else
	if [ $passed_tests -gt 0 ]; then
	    echo -e "\nEverything seems to be OK"
	else
	    echo -e "\nAll tests were skipped."
	    # make sure a failure will be reported by the exit status
	    failed_tests=100
	fi
    fi
    echo

    # Clean up.
    [ "$cleanup" == "yes" ] && rm -f *.log

    # Get epoch seconds after we are done and compute time difference.
    testsuite_end=$(date +%s)
    timediff_sec=$[ testsuite_end - testsuite_start ]

    RUNTIME="Total run-time of the testsuite: \
    $(printf '%02d:%02d:%02d' $[timediff_sec / 3600] \
    $[(timediff_sec % 3600) / 60 ] $[timediff_sec % 60])"

    echo $RUNTIME
    echo ""
fi

exit $failed_tests


# Local Variables:
# mode: shell-script
# coding: utf-8
# End:
