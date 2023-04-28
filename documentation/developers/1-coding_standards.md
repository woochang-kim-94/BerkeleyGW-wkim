Coding Standards
================

In all contributions to the code, please follow the spirit of these guidelines as much as possible, and try to rectify violations you come across.

Guidelines for Committing
-------------------------

 1. Always check over your changes with `git status` and `git diff` to be sure that you really mean to change what is there. Avoid accidentally committing changes. Remember that you can commit (or diff) files or directories individually, in case you have changes in your working directory which you do not wish to commit.
 1. Be descriptive and verbose in your commit messages. Ideally, for all but the most trivial changes, there should be a paragraph explaining the goal of the changes and how they accomplish it. Follow the style of other commits. This allows others to follow along with what you are up to and help find any problems, and indeed will be helpful for you looking back at your changes in a few months. Make comments covering all of the changes in your commit.
 1. Correct spelling is required!
 1. Make many small independent commits rather than a smaller number of large ones, which makes it easier for others to see what you are doing and what the origin of any problem that arises might be (especially with the buildbot and occasionally manually searching for the revision where a certain problem arose). If you find yourself making extensive changes about different issues that will result in a large and complicated commit, stop and make changes about just one of the issues in a separate checked version, and commit that, before continuing with the other issues.
 1. When you add new input variables, as in the `inread.f90` files, you should at the same time add their description to the sample input files (`.inp`), so that we have documentation on the possible options.
 1. Copy or repeat code as little as possible as it makes the code harder to understand and maintain. You can generally reuse the same code with if conditions or pulling out the part you want to use again as a subroutine.
 1. Don't change indentation at the same time that you make meaningful changes, or else it difficult to see what meaningful changes were actually made. Use two separate commits for clarity.
 1. If you find your local version to be in conflict, be very careful about resolving it to be sure that you are not undoing other recent work in the repository. If you have trouble or are unsure how to proceed, ask for help.
 1. If you introduce verbatim code from some external source, document carefully in the code itself and in your commit message where it came from and what its licensing provisions are (and how they are compatible with our BSD license!).
 1. Run `dev-scripts/make_deps.py` and `dev-script/mako_preprocess.py` whenever you change dependencies or Mako template files (`.f90p`) respectivelly. See further information below.

''Unless otherwise specified, the rules below apply to Fortran source. Some can be applied to C or C++ source as well.''

Required Elements
-----------------

 1. Every subroutine, function, module, and program must begin with `use global_m` and `implicit none` and `private` (for default scope). Exceptions: modules used by `global_m` cannot use `global_m`, so they should use the constituent modules of `global_m` explicitly. Module routines should not use any modules or contain `implicit none`, which should be at module scope instead. Modules containing only public parameters, types, and/or interfaces may have default `public` scope. Modules should list what ought to be publicly accessible via `public :: `, but should not list anything explicitly as `private` since this should be the default. All interfaces should contain `implicit none`. The constituent modules of `global_m` should not be used explicitly.
 1. All source files must have the line `#include "f_defs.h"`.
 1. All source files must have a header identifying the type, name, original creator, and date, and have some comments about the purpose and usage.
 1. Every subroutine and function must have `PUSH_SUB(routine_name)` as the first executable line (''i.e.'' the earliest position that will compile, after all the variable declarations) and have `POP_SUB(routine_name)` as the last line before any `return` statements. A subroutine contained within another should use this style: `PUSH_SUB(routine_name.contained_name)`. If there are multiple possible exit points from a routine, each must have a `POP_SUB`. Exceptions: 
   a. `push_pop_m` itself and any modules it uses
   b. any routine called before the debug level has been set
   c. certain functions that are called so many times as to swamp the debugging output, ''e.g.'' those in `Common/random.o`. In this case, where the `PUSH_SUB` would be, this line must be written: `! no push/pop since called too frequently`.

Program Structure
-----------------

 1. We use exactly two spaces for indentation, for the body of blocks such as `if` and `do`, and for continuation lines. Do not use tabs anywhere in source code. The largest block in the file (typically `module`, `program`, or `subroutine`) should begin in column zero. Preprocessor directives (''e.g.''  `#ifdef`, `#include`) should be in column zero. The keyword 'contains' should be at the same indentation as the beginning and end of the program unit. (See "Editor Configuration" below for help with this.)
 1. All named program blocks (''e.g.'' `program`, `subroutine`, `function`, `module`, `interface`) should repeat the name in the `end` statement. ''e.g.'' `program sigma ... end program sigma`.
 1. Intrinsics and keywords should be written in lower case. Preprocessor macros should be written in all upper case.
 1. Fortran source files that do not have the canonical extension `.f90` (''e.g.'' `.F90` or `.h`) should end in a syntax highlighting block for emacs:
{{{
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:  
}}}

Names
-----

 1. No variable name may be a single letter, since this makes it very hard to search for where it is used or modified. See https://github.com/jornada/DOENDO if you need a quick tool to rename single-letter variables.
 1. Names of variables, arguments, functions, routines, modules, and programs should be descriptive and comprehensible, and if appropriate, comments should explain what they are.
 1. All binary executables shall have the suffix `.x.`, and if they respond to the presence or absence of `-DCPLX` in `flavor.mk`, they should be named `.real.x` or `.cplx.` (in the Makefile: `$(FLAVOR).x`).
 1. Use of the term `tape` in comments, variable names, or output is banned.

Declarations
------------

 1. All subroutines and functions must list the argument types in the same order as in the call, and declare an `intent` for each. A blank line should separate the arguments from any subsequent declarations of local variables. 
 1. Variable and argument declarations should use the separator `::`, as in `integer :: ii`.
 1. Subroutines that follow a CONTAINS statement in a module or program should not specify the size of array arguments, ''e.g.'' write `gvec(:,:) !< (3, ng)` in preference to `gvec(3, ng)`. This avoids possible creation of array temporaries (wasting time and memory copying), and allows the array passed to the routine to have a larger leading dimension.
 1. Do not declare `optional` arguments to routines without an explicit interface. Make sure they are in a module if you want to do that. Calls that pass optional arguments should call them by keyword, ''i.e.'' `dont_die = .true.` which helps identify it as optional.
 1. Use Doxygen tags so that variable names and subroutine descriptions are automatically documented. See CodingStandardsDoxygen.
 1. Follow the same structure for whitespace (''e.g.'' space between commas and around an equal sign) as Python's [https://www.python.org/dev/peps/pep-0008/  PEP8].
 1. Avoid unnecessary blank lines within ''quanta of computation''. For example, if you have a 10-line piece of code to sort G-vectors, you can add blank lines before and after this piece of code, not but within it.

Modules
-------

 1. All Fortran files should be inside a module! This is allows the compiler to generate automatic interfaces, and it is also required by the automatic dependency system (`dev-scripts/make_deps.py`).
 1. Every module must be named with the suffix `_m`, ''e.g.'' `global_m`, and be in a file without that suffix, ''e.g.'' `global.f90`. This naming scheme is used by the build system in `Common/common-rules.mk` for dependencies, and also helps to keep the namespaces for modules, types, routines, etc. clearly separated.
 1. Use of `only` list for module usage is highly discouraged. Use of the functionality to change the names in an `only` list is prohibited, except for when handling a module from an external library in which case it is merely highly discouraged: change the name in our code instead.
 1. Module `use` statements should be in alphabetical order.
 1. Routines from external libraries should have interfaces defined in a module (ideally one provided by the library itself). This module must be used in any routines making the library calls. All FFTW calls should be made via the wrappers in `fftw_m`.
 1. MPI calls from Fortran should use the global integer `mpierr` from module `peinfo_m` for the (basically useless) last argument, rather than declaring its own local integer.

Data types
----------

 1. All real and complex types must be declared via the types such as `real(DP)` and `complex(DPC)` defined in `Common/nrtype_m` (in `Common/global_m`), or the template types `SCALAR` if the type varies between real and complex flavors of the code. Use `ZERO` and `ONE` for 0 and 1 values of the `SCALAR` type, and `MYCONJG` for the complex conjugate of the `SCALAR` type. For MPI calls, use data types `MPI_REAL_DP`, `MPI_COMPLEX_DPC`, and `MPI_SCALAR`. Generally, use macros from `f_defs.h` instead of `#ifdef CPLX` whenever possible. Don't use single-precision variables, and make constants double-precision too to avoid accidental loss of precision: ''e.g'' `0d0`, `dble(4)`, `CMPLX(0,1)`, `17d-4`, `0.5_dp`.
 1. To perform an explicit type conversion, use the macros in `Common/f_defs.h`. Use the constants and parameters from `Common/nrtype_m`.

Output
------

 1. All warnings and errors must be written to unit 0 (standard error) rather than unit 6 (standard output). Warnings should begin with `WARNING:`. Errors should be handled with the `die` function, never with `stop` or `exit`.
 1. Use `write` rather than `print`, and supply a format rather than using `*`. In `write` statements, do not break up string literals with continuation lines, so that the text can be `grep`'d for in the code. Start your output at the second column, and avoid extra spaces (eg: do `write(6,'(1x,a)') 'Blah'`, and not `write(6,'(2x,a)') ' blah'`)
 1. Write full sentences to the output including the period. For example, write "Interpolating matrix with super algorithm X." instead of "Super Algorithm X." Only capitalize the first letter of each sentence.
 1. To write an empty line, use `write(6,'()')`, and not `write(6,*)` nor `write(6,'(a)') ''`. The first option is the simplest way to write an empty line without a trailing space.
 1. Avoid unnecessary white spaces when writing integer and double precision numbers. Use the `i0` edit descriptor  for integer formats of automatic and minimum sizes, and `f0.2` for general floating points with 2 decimal places and the minimum possible width. For scientific notation, the automatic width doesn't work, so you'll need to manually minimize the amount of whitespace. For example, if you know that you will print a positive number, use `es9.3e2`. So, the output should look like (`i0` and `f0.1`):
{{{
 Memory necessary to store v(q) for 42 q-points: 103.6 MB.
}}}
    instead of (eg: `i6` and `f8.1`):
{{{
 Memory necessary to store v(q) for     42 q-points:    103.6 MB.
}}}
    and absolutely never like (eg: `i1` and `f4.1`):
{{{
 Memory necessary to store v(q) for * q-points: **** MB.
}}}
 1. When printing k-points, use a real edit descriptor with a fixed width: either `f9.6` or `f6.3`. The only exception is when writing k-points that are not in either the 1st BZ or the conventional [0,1) range.
 1. Don't write ugly/scary/unnecessary information to the standard output. If you need to write debugging information, enclose your write statements with `#ifdef VERBOSE / #endif`.
 1. Avoid unnecessary white spaces between sentences, and group output statements and parameters belonging to the same category with a simple itemization. Use a single dash followed by a white space for the items. Example:
{{{
 Calculation parameters:
 - Cutoff for Monte-Carlo averages (Ry): 1.000E-12
 - Cutoff for dielectric matrix (Ry): 60.00
 - Broadening (eV): .0001
 - Number of k-points: 42
 - We will not use truncation.

 Reading wavefunctions from file WFN.
 Successfully read file WFN.
}}}

Comments
--------

 1. With the exceptions of section headers, all major comments should look like the following `! JRD: Comment` -- replace "JRD" with your initials.
 1. Be sure to write explicit comments explaining the purpose and function of any unusual constructs introduced to the code, especially if they could look wrong or useless to others.
 1. To avoid preprocessor warnings (and XLF compilation failures), an odd number of single (''i.e.'' apostrophes) or double quotes may not be used in comments (use ` instead).
 1. Keep original lines as comments if you modify third-party routines such as `Common/sort`, `Common/random`, `MeanField/EPM`, or `MeanField/spglib-1.0.9`, especially if your change is an extension of the original code.

Required Macros and Routines
----------------------------

 1. Do not use `call flush` (an intrinsic function call) but instead `flush`, an intrinsic introduced in Fortran 2003.
 1. The intrinsics `allocate` and `deallocate` should never be used. Instead use the macros `SAFE_ALLOCATE` and `SAFE_DEALLOCATE` (for arrays) and `SAFE_DEALLOCATE_P` (for pointers). Danger: DO NOT use these macros in a one-line `if` statement; they must be used with `if(...) then` or else some of the generated lines are outside the condition.
 1. Use routine `open_file` instead of the intrinsic `open` when opening file units, and `close_file` instead of `close`, which will check that the operation occurs properly and write an error if not.
 1. For Fortran string concatenation, instead of the built-in `//` operator, use the overloaded `+` operator, defined in `message_m.mod`. Otherwise, problems can occur with cpp from gcc version 5 and later.

Forbidden
---------

 1. Use of the `goto` and `continue` keywords, as well as the related `err` and `end` parameters of `close`, `open`, and `read` intrinsic functions, is categorically forbidden. Use the `cycle` and `exit` keywords and the `iostat` parameter instead. Keywords `sequence` and `common` are also forbidden.
 1. Do not declare anything `external`. Use modules for functions and subroutines defined in this package, and interfaces for intrinsics and libraries.
 1. Do not assume any value for a loop counter outside its loop. This works in C, but not necessarily in Fortran. You must save the value in another variable.
 1. Continuation lines shall not begin with `&`.

Portability
-----------

 1. All code must work correctly at `-O3` with gcc, PGI, and Intel compilers. If it does not, it likely contains an error even if it seems to work without optimization. Making sure it works at `-O3` with sun, XLF, open64, etc. is strongly encouraged.
 1. Do not make explicit reference in source files to the compiler being used. Instead use the macros in `Common/compiler.h`, and add new ones if (strictly!) necessary. Even better, rely as much as possible on the `iso_c_binding` module and call C functions if you require a particular feature that is not present in the Fortran 2003 standard.
 1. All Fortran source code should be valid under the Fortran 2003 standard.  Do not rely on features from the Fortran 2008 or later standard. You can check the various standards [https://gcc.gnu.org/wiki/GFortranStandards#Fortran_2003 here].

Scripts
-------

 1. All scripts shall contain a "shebang" in the first line, ''e.g.'' `#!/usr/bin/env perl`, `#!/usr/bin/env python`, etc. Exception: always use the following "shebang" for bash scripts: `#!/bin/bash -l`. The reason is that in many supercomputers, one has to submit jobs invoking bash as a login shell, and it's not possible to have two arguments in the shebang (''i.e.'', `#!/usr/bin/env bash -l` is invalid).
 1. Python scripts shall be compatible with version 3. The utility `2to3` can be used to convert from v2 to v3.

Makefiles
---------

 1. All Makefiles must begin (using the appropriate number of `../` for the depth of the directory) with:
{{{
PREFIX=../..
include $(PREFIX)/Common/common-rules.mk
}}}
 1. All Makefiles must end with the line `include make_deps.mk`. Run `dev-scripts/make_deps.py` every time you add any file or change any dependency to automaticall create the `make_deps.mk` file.
 1. All Makefiles shall define targets `default` and `all`.
 1. Use the pre-existing macros, variables, etc. defined via `Common/common-rules.mk` as much as possible.
 1. Explicit reference to the compiler being used should only be made in `Common/common-rules.mk`.
 1. Object files (`.o.`) and modules (`_m.mod`) should have dependencies only on modules they used (`_m.mod`) or header files they include (`_inc.F90`), but should not have dependencies on object files. Dependencies for a file that creates a module should list both object file and module at the beginning: ''e.g.'' `mtxel.o mtxel_m.mod : $(COMMON)/fftw_m.mod`. In most cases, these rules are automatically written in `make_deps.mk` by `dev-scripts/make_deps.py`.
 1. All parts of the code must build correctly with parallel make (`make -j`) with any number of threads. If that fails, there is a mistake, usually missing dependencies.
 1. Every target that is not the name of a file should be in the `.PHONY` list in `Common/common-rules.mk`.

Editor Configuration
--------------------

To configure your editor to use our indentation convention:
 * `emacs` : in menus, select Options->Advanced(Customize)->Emacs->Programming->Language. Click on F90, then on Indent. Now set all the numeric values here to 2, and Save and Done. Equivalently, paste the following lines to your `~/.emacs` file:
{{{
(custom-set-variables
 '(f90-continuation-indent 2)
 '(f90-do-indent 2)
 '(f90-if-indent 2)
 '(f90-type-indent 2)
 '(vc-follow-symlinks nil))
(setq c-default-style "linux"
      c-basic-offset 3)
}}}
The last line also suppresses the warning emacs will give you when you open a symbolically linked file, which is a good way to handle `arch.mk` and `pw2bgw.f90`, among other things.
 * `vi` (recommended editor): paste the following lines to your `~/.vimrc` file:
{{{
:let fortran_free_source=1
au BufNewFile,BufRead *.f90p,*.fp set filetype=fortran
au filetype fortran set expandtab softtabstop=2 shiftwidth=2 autoindent
set backspace=2
}}}
 * `nano` or `pico`: don't bother.

