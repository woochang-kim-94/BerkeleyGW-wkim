# BerkeleyGW developers guideline

1. Follow coding guidelines listed in the file `1-coding_standards.md`

2. If you are associated with the BerkeleyGW development team, read the files
   `2.1-git_setup.md` and `2.2-git-development.md` to configure your git
   environment.

3. These codes are to be used on various machine architectures
   employing various compilers and libraries. Therefore, try not
   to write machine-specific code or quick fixes that will have to
   be reversed later.

4. Limit compilation changes as much as possible to `arch.mk` and ensure
   that this is suitable for all codes in `Epsilon/`, `Sigma/`, and `BSE/`.
   Hopefully individual Makefiles will only require the addition of
   new source-code files.

5. Test the codes using the `examples/` and please augment the `examples/` with
   your own test cases (provided they are not too computationally expensive).

6. Other documentation intended for developers are in the `manual` subfolder.
