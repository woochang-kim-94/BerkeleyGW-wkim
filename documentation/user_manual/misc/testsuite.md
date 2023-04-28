# BerkeleyGW testsuite

BerkeleyGW ships with a comprehensive testsuite, which is run automatically
every time a developer changes the code on a number of different architectures
and with a number of different compilers and math libraries. Every time
BerkeleyGW is ported to a new architecture, we strongly recommend the user to
run the testsuite to ensure correctness of the result produced by the code.


## Running tests

Below are the typical steps to run the BerkeleyGW testsuite (assuming that you
are in the `testsuite` directory):

  - Build the BerkeleyGW binaries with "make all-flavors" or "make -j all-flavors".
    Testsuites will check both the real and complex versions of the binaries.

  - In serial: type `make check` in the `testsuite/` directory or in the root
    BerkeleyGW directory. `make check-save` will do the same, but retain the
    working directories for debugging.

  - In parallel without a scheduler: type `make check-parallel`. Set
    environment variable `SAVETESTDIRS=yes` to save working directories. The
    tests use 4 cores, so be sure that at least that many are available, or use
    `BGW_TEST_MPI_NPROCS` (see below). See, for instance, the script
    `singularity.scr` in the `testsuite/` directory, which is executed directly
    and without a scheduler.

  - In parallel with a scheduler (e.g., SLURM, PBS, etc.): type `make
    check-jobscript` or `make check-jobscript-save`, which will execute the job
    script for the machine you are using, as specified in `arch.mk` by the line
    `TESTSCRIPT`. See example job scripts for various machines in this
    directory, which typically have a `*.scr` extension.  For instance, to run
    the testsuite on the KNL nodes of the Cori supercomputer at NERSC, go to the
    `testsuite` directory and type `sbatch cori2.scr`. Run `sbatch cori2_haswell`
    if you want to run the testsuite using the Haswell nodes.


The following environment variables control how the testsuite gets executed,
and are typically put into a submission script:

  - `BGW_TEST_MPI_NPROCS`: set to overrule the number of processors listed in
    the test files. Useful for testing correctness of parallelization, if you
    don't have 4 cores, or to run the testsuite faster with more cores.

  - `OMP_NUM_THREADS`: if you build BGW with the `-DOMP` flag, this environment
    variable will set the number of OpenMP threads used in running the testsuite. 

  - `TEMPDIRPATH`: sets the scratch directory where temporary working
    directories will be created. Default is `/tmp`, but on some machines you
    may not have permission to write there.

  - `MPIEXEC`: sets the command for parallel runs. Default is `` `which
    mpiexec` ``. Note that `mpiexec` and `mpirun` are generally equivalent
    (except for Intel MPI, in which case `mpirun` is recommended). Set this if
    you are using a different command such as `ibrun` (SGE parallel
    environment), `runjob` (BlueGene), or `aprun` (Cray), if you need to use a
    different `mpiexec` than the one in your path, or if some options need to
    be appended to the command. Depending on the value of this variable, three
    styles of command line are used: `ibrun`, `runjob`, and `mpirun`/`aprun`.

  - `BGW_TEST_NPROCS`: deprecated, same as `BGW_TEST_MPI_NPROCS`, provided for
    compatibility with versions 1.0.x.

  - `MACHINELIST`: if set, will be added to command line for MPI runs. This is
    needed in some MPI implementations.

  - `EXEC`: if set, will be added before name of executable. Used to run code
    through `valgrind`.


If you want to run a subset of the tests, you can enable/disable individual
tests with the `toggle_tests.py` script, described below.

## How to update the testsuite manually

If you have updated one of the input files in a testsuite (e.g., Si-Wire-EPM/plotxct.inp),
or you have fixed a bug in the code somewhere that would affect the result of a specific
output file (e.g., plotxct.out), you can modify the BGW_CHECK_FLAGS environmental variable
(e.g., export BGW_CHECK_FLAGS=" -d Si-Wire-EPM/ ") in the corresponding testsuite script and run it
(e.g., sbatch cori2_haswell.scr). Check the file "test.out" inthe testsuite directory and locate
the failed test with the keyword "FAIL". The "Calculated value" and "Reference value" will be
printed out before the failed test. Then go to the corresponding directory and open the .test file
(e.g., Si-Wire-EPM/Si_wire.test). Find the reference value and replace it with the new value.
Rerun the testsuite script to double-check that all the tests in the corresponding directories are OK.

## Scripts

The following are the scripts shipped with the testsuite:

1. `run_testsuite.sh`: runs all tests from filenames ending in `*.test`
   that can be found under the current folder and which are enabled.

2. `toggle_tests.py`: enable or disable a particular test. This script is
   supposed to be executed interactively.

3. `run_regression_test.pl`: called from `run_testsuite.sh` for each individual
   test. The syntax is: `../run_regression_test.pl [-s] [-p] -D ../../bin -f
   [testname].test`.  Include the `-s` flag to run in serial; otherwise it will
   be in parallel.  Include the `-p` flag to preserve the working directory.

4. `*.scr`: job scripts for running the testsuite on various parallel machines.

5. `interactive*.sh`: scripts to run in parallel, in interactive mode. You have
   to follow the instructions in the script, not just run it.

6. `queue_monitor.pl`: a Perl script for parallel tests with a SLURM or PBS
   scheduler, which submits a job and monitors the queue to see what happens.
   Formally used by the BuildBot (deprecated).

7. `fix_testsuite.py`: (for developers) can be used to create or update match
   reference values in test files.

8. `no_hires.sh`: Some clusters will not have the Perl package `Time::HiRes`
   installed which is needed for timing in `run_regression_test.pl`. You can
   just deactivate the usage with this script.

9. `buildbot_query.pl`: (for developers) can be used to find the range of
   values obtained on the buildbot workers for a particular match, to set
   values in the test files (deprecated).


## Writing tests

If you are interested in modifying the testsuite and/or writing new tests,
please refer to the developer's manual.
