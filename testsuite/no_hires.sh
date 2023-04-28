#!/bin/bash

# Some clusters will not have the Perl package Time::HiRes installed which is needed for
# timing in the testsuite. You can just deactivate the usage with this script.

sed -i.bak '/gettimeofday/d' run_regression_test.pl
sed -i.bak '/elapsed/d' run_regression_test.pl
