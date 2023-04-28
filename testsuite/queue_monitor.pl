#!/usr/bin/env perl
# Originally by David Strubbe, September 2010
# perl script for running testsuite with PBS or SLURM scheduler

use warnings;
use Time::HiRes qw(gettimeofday tv_interval);

$write_interval = 5 * 60;
$output_file = "test.out";
$qstat_monitor_interval = 30;
$tail_monitor_interval = 0.5;

$starttime = [gettimeofday];

if($#ARGV < 4 || $#ARGV > 5) {
    $numargs = $#ARGV + 1;
    print "Number of arguments = $numargs\n";
    die "Usage: perl queue_monitor.pl {slurm|pbs} job.scr timeout_seconds jobtitle end_delay [queue]\n";
}

$scheduler = $ARGV[0];

if(lc($scheduler) eq "pbs") {
    $field_queue = 5;
    $qstat_line_offset = 2;
    $name_option = "-N";
    $queue_option = "-q";
    $queue_job_command = "qstat";
    $queue_queue_command = "qstat";
    $kill_command = "qdel";
    $submit_command = "qsub";
} elsif(lc($scheduler) eq "slurm") {
    $field_queue = 1;
    $qstat_line_offset = 1;
    $name_option = "-J";
    $queue_option = "-p";
    $queue_job_command = "squeue -j";
    $queue_queue_command = "squeue -p";
    $kill_command = "scancel";
    $submit_command = "sbatch";
} else {
    die "Unknown scheduler choice '$scheduler'; must be 'slurm' or 'pbs'.\n";
}

$script = $ARGV[1];
$timeout_seconds = $ARGV[2];
$jobtitle = $ARGV[3];
$end_delay = $ARGV[4];
if($#ARGV == 5) {
    $queue = $queue_option . " " .$ARGV[5];
} else {
    $queue = "";
}

$stderr = $jobtitle . ".err";
$stdout = $jobtitle . ".out";
system "rm -f $stderr $stdout $output_file";
# on some systems, if these files are already present there may be an error
# also, clearing things out makes sure the timestamps will mean what we think they do

$submit_command_line = "$submit_command $name_option $jobtitle -e $stderr -o $stdout $queue $script";
print $submit_command_line . "\n";
$job_out = qx($submit_command_line);
# Perl gives error code shifted, for some reason.
$exit_code = $? >> 8;
# SLURM example stdout: Submitted batch job 3010046
print $job_out;
if($exit_code) {
    die("Job submission failed.\n");
}
chomp $job_out;
@job_out_bits = split(' ', $job_out);
$jobname = $job_out_bits[$#job_out_bits];

$last_write = 0;
$running = 0;
$writing = 0;
$test_finished = 0;
$job_finished = 0;
$fileage = 0;

do {
    sleep($qstat_monitor_interval);
    
    $qstat = qx($queue_job_command $jobname 2> qstat_err_);
    $exit_code = $? >> 8;
    $qstat_err = `cat qstat_err_`;

    # this means qstat has failed, and all bets are off
    if($exit_code && $qstat_err !~ /Unknown Job Id/ && $qstat_err !~ /Invalid job id/) {
	print "error checking job status: " . $qstat_err . "\n";
    }

    $job_status = ""; # default in case not set below
    chomp $qstat; # there may be an end of line here otherwise
    @lines = split(/\n/, $qstat);
    # search to be robust against extra junk before the info we need
    for $iline (0..($#lines-$qstat_line_offset)) {
	if(lc($lines[$iline]) =~ "job id" || $lines[$iline] =~ "JOBID") {
	    @fields = split(' ', $lines[$iline+$qstat_line_offset]);
	    $job_status = $fields[4];
	}
    }

# SLURM:
#             JOBID   PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
#           3010272 development test_BGW dstrubbe  R       0:54      1 c557-803

# PBS: (also may be 'Job id')
#
# Job ID                    Name             User            Time Use S Queue
# ------------------------- ---------------- --------------- -------- - -----
# 1195769.edique02           slurm            dstrubbe        00:00:16 R debug

    $currenttime = [gettimeofday];
    $elapsed = tv_interval($starttime, $currenttime);

    if($running eq 0 && $job_status eq "R") {
	print "Job is running.\n";
	$running = 1;
    }

    if($running eq 1 && $writing eq 0) {
	# don't tail before MPI job has started writing
	# otherwise some of file from previous run may be written out
	# sometimes this age will never go negative, only less than a second but positive
	if(-e $output_file && -M $output_file <= 1e-4) {
	    $writing = 1;
	    print "Parallel job output:
===========================================\n";
	    system "tail -n +1 --pid=$$ --sleep-interval=$tail_monitor_interval -f $output_file &";
	    # giving PID of this Perl script will make tail stop when this script does
	    # -n +1 is needed so we start at the beginning of the file
	}
    }

    if($running eq 0 && $elapsed - $last_write > $write_interval) {
	$last_write = $elapsed;
	print $qstat . "\n";
    }

    # we cannot stop until the job standard out file is present
    if(-e $stdout) {
	# if the standard error does not say the equivalent of
	#   qstat: Unknown Job Id 583497.sched-00.scs.lbl.gov
	# then qstat has failed and we don't know if the job is finished!
	# slurm_load_jobs error: Invalid job id specified
	# we will not consider the job done in this case, and continue waiting.
	if(($qstat eq "" && ($qstat_err =~ /Invalid job id/ || $qstat_err =~ /Unknown Job Id/))
	    || ($job_status ne "R" && $job_status ne "PD" && $job_status ne "Q")) {
	    $job_finished = 1;
	}

	$status = `cat $stdout`;
	if($status =~ /Exit status/) {
	    $test_finished = 1;
	}
    }

    # kill if timeout occurred; do not kill if already stopped (i.e. we are in end_delay period)
    if($elapsed > $timeout_seconds && $job_finished eq 0) {
	system("$kill_command $jobname");
	sleep($end_delay);
	print "===========================================\n";
	print $qstat . "\n";
	$queue = $fields[$field_queue];
	print "Here is the rest of the $queue queue:\n";
	system "$queue_queue_command $queue";

	system "rm qstat_err_";
	die("Job has not completed before timeout. Killing job.\n");
    }

    if(-e $output_file) {
	$fileage = time - (stat($output_file))[9];
    }
} while (($job_finished eq 0 && $test_finished eq 0) || (-e $output_file && $fileage < $end_delay));
# We stop when the job is no longer running (i.e. not in queue, or listed as completed or exiting)
# or "Exit status" has appeared in the standard out, and its standard out file exists, 
# and the job has not written to the output file in the last few seconds.

system "rm qstat_err_";

print "===========================================
Parallel job has finished.\n";

if(-e $output_file) {
    if($running eq 0) {
	# in this case, the job finished before one $qstat_monitor_interval, and so it will not have been output yet
	print "Parallel job output:
===========================================\n";
	system "cat $output_file";
    }

    if($fileage >= $end_delay + $qstat_monitor_interval) {
	print "No write to output for $fileage seconds (> $end_delay).\n";
    }
}
else {
    print "No output file was written.\n";
}

# temporary debugging info below
if($job_finished == 0) {
    print "No, job did not finish!\n";
}
if($test_finished == 0) {
    print "No, test has not finished!\n";
}
# end temp debug info

$errors = `head $stderr`;
if($errors ne "") {
    print "\nStandard error from job script:\n";
    system "cat $stderr";
    print "===========================================\n";
}

$status = `cat $stdout`;
if($status !~ /Exit status = 0/) {
    print $status . "\n";
    die("Test failed.\n");
}
