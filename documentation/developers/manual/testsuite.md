# BerkeleyGW testsuite

Before you start, please take a look at the testsuite section in the BerkeleyGW
users' manual (BerkeleyGW/documentation/user_manual/misc/testsuite.md).


## Writing tests

Below is a reference guide on how to write tests. This is intended mostly for
developers.  The test files consist of lines beginning with one of a set of
tags, parsed by the Perl script. Comment lines beginning with `#` will be
ignored.

- `Test : title`:
  Write a title to output to identify the test. Should be the first tag in the file.

- `Banner : text`:
  Write a heading in a box to identify a specific part of a test. Helpful when
  input files are generated by `sed` and do not have a distinct name.

- `Enabled : Yes/No`:
  If `Yes`, will be run; if `No`, will not be run. Use to turn off a test
  without deleting it from the repository. *Should be the second tag in the
  file*.

- `TestGroups : group-name [group-name2 [...]]`:
  The `run_testsuite.sh` script can be run with argument `-g` and a list of
  groups. Then tests will only be run if there is a match between the argument
  list and the list in `TestGroups`. Examples of groups that could be used:
  `parallel`, `serial`, `long`. This tag is actually read by `run_testsuite.sh`
  rather than `run_regression_test.pl`.

- `Executable : program.x`:
  Name of executable to be run (path is `bin/` directory). Persists for
  multiple runs, until next `Executable` tag.

- `Processors : integer`:
  Number of processors to use. Ignored if `mpiexec` is not available. May not
  exceed the number in the reservation in the job scripts. Set to special
  value `serial` to run without `mpiexec`. This flag is ignored if
  `BGW_TEST_MPI_NPROCS` is set and `Processors` is not `serial`.

- `Command : shell-command`:
  Specify a shell command, which will be executed in the working directory.
  Useful for renaming log files that would be overwritten by a subsequent run,
  or preprocessing input files by `sed` or other utilities (after putting in
  working directory with `Copy`).

- `Copy : file.in [file_newname.in]`:
  Copy a file from test directory to run directory. If new name is not
  supplied, original name will be used. Useful if you want to `sed` a file
  before it is run.

- `Unpack : data.tar.gz`:
  The file `data.tar.gz` in the test directory will have `tar xzf` run on it,
  with resulting files going to the working directory.

- `Output : file.out`:
  Output will be piped into the file named here, in the working directory.

- `Arguments : arg1 arg2`:
  The current `Executable` will be executed (in serial) with the text given
  here as command-line argument(s). No files will be copied to working
  directory. Use redirection with `>` to capture the output if necessary.

- `Input : file.in [file_newname.in/PIPE/CMDLINE/NONE/NOCOPY]`:
  The file named here will be copied from the test directory to the working
  directory, for use as input (unless `NOCOPY` is specified). If it is followed
  by `PIPE`, the file will be piped into the executable. If it is followed by
  `CMDLINE`, the file will be given as a command-line argument to the
  executable. If it is followed by `NONE`, no file will be copied (useful if an
  input file was already generated by `sed` in the working directory). If it is
  followed by another name, the file's name in the working directory will be
  the new name. This tag causes actual execution of the run that has been set
  up by previous tags.

- `Precision : 1e-5`:
  A floating point number, the tolerance for testing whether a match has passed
  or failed. Persists until next `Precision` tag. Default is `1e-4`.

- `match ; name ; COMMAND(..); reference-value`:
  Extracts a calculated number from a run and tests it against the reference
  value. The name is an identifier printed to output. The number is extracted
  as the standard output from the listed `COMMAND`, which is one of this set:

  - `GREP(filename, 'search-text', field, offset)`:
    Finds the first line in file containing `search-text` (which MAY NOT
    contain a comma), and returns the specified field of that line. `Field` is
    meant in the sense used by `awk`, i.e. the line is parsed into white-space
    separated groups, indexed starting with 1. The optional `offset` directs
    the use of the Mth line after the line containing `search-text`. No offset
    is equivalent to `offset`=0. This is the most robust of the commands and
    should be used when possible in preference to the others.

  - `SIZE(filename)`:
    Returns the size of the specified file via `ls -lt`. Useful for binary
    files whose contents cannot easily be checked. Not very robust.
  - `LINE(filename, line, field)`:
      Returns the specified field of the specified line number from the file.
      Use `GREP` instead if possible.
  - `SHELL(shell-command)`:
      The result is the standard output of the listed command. Deprecated; use
      `GREP` or `LINE` unless absolutely necessary.

- `STOP TEST`:
  For debugging purposes, to halt test after the part you are studying has run.