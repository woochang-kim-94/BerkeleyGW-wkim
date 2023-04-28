#!/usr/bin/env python

# This script will automatically create the dependency graph for all directories
# and write the Makefile dependency into make_deps.mk
#
# Rules:
# - Only one "use" statement per line
# - `only` and `intrinsic` keywords are accepted
# - All modules that *don't end with _m* will be ignored, as they are assumed
#   to be system modules
# - You can specify the path of the module by appending the keyword
#   `-I "path"` after the `use module_m` keyword. Example:
#   use foo !-I /etc/foo
# - You can also manually ignore the automatic dependency by appending the keyword
#   `-nodep`. Example:
#   use foo !-nodep
# - If the path is not set, it defaults to `./` if a corresponding .f90 or .F90
#   file exists, otherwise it defaults to $(COMMON)/
#
# NOTE: we don't need to put the dependencies of _inc files here, because they
# are taken care by in common-rules.mk
#
# Felipe H. da Jornada (2015)


from __future__ import print_function, unicode_literals

from glob import glob
import re
import os


def str_sorted(x):
    return sorted(x, key=lambda x: x.encode('utf-8'))


class FortranDependency(object):
    # These regex use python bytes instead of str and are used to match the
    # module names in the Fortran file
    re_is_prog = re.compile(br'^ *program *([a-zA-Z0-9_]*)', re.MULTILINE)
    re_dep_list = re.compile(br'^ *use *(.*)$', re.MULTILINE)
    re_name = re.compile(br'[^ ,!]+')
    re_path = re.compile(br'!-I *(.*?) *$')
    re_nodep = re.compile(br'!-nodep *$')

    # Regex to determine if the file is a cpp include
    re_is_inc = re.compile(r'_inc.[Ff]90[p]$')
    # Regex used to rename the file extension
    re_f90_ext = re.compile(r'\.[fF]90')

    # Dependency tree, shared across all class instances.
    dep_tree = {}

    def __init__(self, fname, dname, dir_fnames=None, verbose=False):
        '''Initializes and parses a Fortran file for its dependencies

        Parameters
        ----------

        fname: str
            Basename of the Fortran file name, including extension.
        dname: str
            Basename of the directory containing the file.
        dir_fnames: str or None
            List containing names of all files in the current directory.
            This is required to figure out if a dependency is present in the
            current directory or under Common/. If `dir_fname` is None, then
            we automatically build the list.
        verbose: bool
            Produces more detailed output.
        '''

        self.fname = fname
        self.dname = dname
        if self.dname == 'Common':
            self.prefix_dir = '$(COMMON)'
        else:
            self.prefix_dir = None
        self.dir_fnames = dir_fnames if dir_fnames is not None else []
        self.verbose = verbose

        self.is_prog = False
        self.prog_name = None

        # Direct module dependencies requried to compute the program/module.
        self.compile_dependencies = set()
        # Full list of object dependencies required to link the program.
        self.link_dependencies = set()
        self.resolving_link_dependency = False
        self.dependency_tree_solved = False

        self.parse()

    def parse(self):
        '''Parse the Fortran file and get the dependency list.
        
        This routine is automatically called when the class is instantiated.
        '''
        fname = self.fname

        # Ignore cpp include files
        if self.re_is_inc.search(fname) is not None:
            return

        # Perform the file regex using python's bytes instead of str as it
        # avoids headaches from source files having different encodings.
        with open(fname, "rb") as f:
            txt = f.read()
            match = self.re_is_prog.search(txt)
            if match is not None:
                self.is_prog = True
                self.prog_name = match.group(1).decode('utf-8')
            raw_dep_list = self.re_dep_list.findall(txt)

        if not len(raw_dep_list):
            return
        if self.verbose:
            print('  - ' + fname)

        # Loop over all module dependencies that we found for this file
        for raw_dep_item in raw_dep_list:
            # Trim the module dependence string and get the real module name
            match = self.re_name.search(raw_dep_item)
            if match is None:
                continue
            # Note that we convert the real module name to a string here
            dep_mod_name = match.group(0).decode('utf-8')

            # All BGW modules finish with _m. So, it's safe to ignore system
            # modules that don't end with _m.
            if not dep_mod_name.endswith('_m'):
                continue

            # If there's a "!-nodep" flag, ignore this module dependency
            if self.re_nodep.search(raw_dep_item):
                continue

            # Determine the path of the module dependency
            match = self.re_path.search(raw_dep_item)
            if match is not None:
                # If there's a path flag "-I dep", use it!
                dep_path = match.group(1)
            else:
                # If not, check if the module dependency is present in the
                # current directory. If not, assume it is under $(COMMON)/
                dep_path = ''
                if not ((dep_mod_name[:-2]+'.f90' in self.dir_fnames) or
                        (dep_mod_name[:-2]+'.F90' in self.dir_fnames)):
                    # Original module file is not here. It must be in Common then.
                    dep_path = '$(COMMON)/'

            # Full dependency path
            full_dep = os.path.join(dep_path, dep_mod_name + '.mod')
            self.compile_dependencies.add(full_dep)

            if self.verbose:
                print('    - {}'.format(full_dep))
        if self.verbose:
            print()

    def get_compile_dependency_str(self):
        '''Get the make_deps.mk compilation dependency string.

        This will generate a string with all prerequisite modules that are
        required to compile the Fortran file and build the object/module
        associated with the instance. Return None if there are no dependencies.
        '''
        if len(self.compile_dependencies) == 0:
            return None

        # Target object file associated with this instance
        fname_o = self.re_f90_ext.sub('.o', self.fname)
        # Target module file associated with this instance (if not a program)
        if self.is_prog:
            fname_mod = ''
        else:
            fname_mod = self.re_f90_ext.sub('_m.mod', self.fname)

        deps = self.compile_dependencies
        deps = str_sorted(deps)
        if self.prefix_dir is not None:
            # If we need to add a prefix_dir, such as $(COMMON), update the
            # variables fname_o and fname_mod, if they are not empty.
            if len(fname_o):
                fname_o = os.path.join(self.prefix_dir, fname_o)
            if len(fname_mod):
                fname_mod = os.path.join(self.prefix_dir, fname_mod)
            # Prerequisite list
            deps = [os.path.join(self.prefix_dir, dep) for dep in deps]
        deps = ['\t' + dep for dep in deps]
        return '{} {}: \\\n{}\n'.format(fname_o, fname_mod, ' \\\n'.join(deps))

    def get_link_dependency_str(self):
        '''Get the make_deps.mk linking dependency string.

        This will generate a string containing, recursively, all prerequisite
        objects that are required to link the them and create the executable
        associated with the instance. Return None if there are no dependencies.
        '''
        if not self.dependency_tree_solved:
            self.solve_link_dependencies()

        if not self.is_prog:
            return None

        sout = ''
        link_deps = sorted(self.link_dependencies,
                key=lambda x: b'' if x=='$(GLOBALOBJS)' else x.encode('utf-8'))
        link_deps = ['\t' + link_dep for link_dep in link_deps]

        targets = []
        for suffix in '', '$(FLAVOR)':
            fname_o = self.prog_name + suffix + '.x'
            if self.prefix_dir is not None:
                fname_o = os.path.join(self.prefix_dir, fname_o)
            targets.append(fname_o)

        return '{}: \\\n{}\n'.format(' '.join(targets), ' \\\n'.join(link_deps))

    def solve_link_dependencies(self):
        '''Visits the dependency tree and collect all object dependencies.
        
        This step is required to produce a list of object files for the linking
        step. It is called automatically by `get_link_dependency_str` if there
        are unsolved dependencies.
        '''
        if self.dependency_tree_solved:
            return

        if self.resolving_link_dependency:
            raise RuntimeError('Circular dependency!')
        self.resolving_link_dependency = True

        # Add dependency from main object file
        fname_o = self.re_f90_ext.sub('.o', self.fname)
        if self.dname == 'Common':
            fname_o = os.path.join('$(COMMON)', fname_o)
        self.link_dependencies.add(fname_o)

        # Add dependency from modules
        for dep in self.compile_dependencies:
            dep_obj = dep.replace('_m.mod', '.o')
            dep_key = dep.replace('_m.mod', '')

            # Remove $(COMMON) in order to get the module key and directory
            if dep.startswith('$(COMMON)'):
                dep_dname = 'Common'
                dep_key = dep_key.replace('$(COMMON)/', '')
            else:
                dep_dname = self.dname

            # Add $(COMMON) to the object name if necessary
            if self.dname == 'Common':
                dep_obj = os.path.join('$(COMMON)', dep_obj)

            if dep_obj == '$(COMMON)/global.o':
                # Special case if the dependency is global.o
                # Note stricly necessary, but makes output files nicer
                self.link_dependencies.add('$(GLOBALOBJS)')
            else:
                # General case: recursively solve module dependencies
                self.link_dependencies.add(dep_obj)
                fort_dep = self.dep_tree[dep_dname]['files'][dep_key]
                fort_dep.solve_link_dependencies()
                self.link_dependencies.update(fort_dep.link_dependencies)

        self.dependency_tree_solved = True
        self.resolving_link_dependency = False


class DependencyParser(object):
    re_is_dir_valid = re.compile(
        r'^\s*include.*make_deps.mk\s*$', flags=re.MULTILINE)

    def parse(self, dirs=None, root=None, verbose=False):
        '''Automatically parses and generate Makefile dependencies.

        Parses dependencies recursively under the current directory, or
        optinally by specifying the `dirs` argument. Generate a `make_deps.mk`
        file under each subdirectory.

        Note that, by default, only subdirectories that contain a Makefile or a
        common-rulkes.mk file and which match against the string "include
        make_deps.mk" are included.
        '''

        self.root = root
        if self.root is None:
            self.find_root_directory()

        print('BerkeleyGW root directory:', self.root)

        os.chdir(self.root)
        if not self.is_bgw_root_dir(self.root):
            raise RuntimeError('{} is not the root directory of BerkeleyGW!'.
                    format(self.root))

        def print_step(title):
            if verbose:
                print()
                print(title)
                print('='*len(title))
            else:
                print(title)

        # If no set of directories is given, walk over all directories recursively.
        if dirs is None or len(dirs) == 0:
            dirs = self.get_dirs()
        old_dir = os.getcwd()

        dep_tree = FortranDependency.dep_tree
        print_step('Parsing directories and files for dependencies')
        for dirname in dirs:
            if verbose:
                sout = 'Directory: {}'.format(dirname)
                print()
                print(sout)
                print('-'*len(sout))

            os.chdir(os.path.join(old_dir, dirname))
            dname = os.path.basename(os.getcwd())

            fnames = str_sorted(glob('*.[fF]90'))
            if verbose:
                print('  Found {} files'.format(len(fnames)))
            if len(fnames) < 1:
                return

            dep_tree[dname] = {}
            dep_subtree = dep_tree[dname]
            dep_subtree['dirname'] = dirname
            dep_subtree['files'] = {}
            for fname in fnames:
                fort_dep = FortranDependency(fname, dname, fnames, verbose=verbose)
                dep_key = fname.replace('.F90', '').replace('.f90', '')
                dep_subtree['files'][dep_key] = fort_dep


        print_step('Solving dependency tree')
        for dname in dep_tree.keys():
            dep_subtree = dep_tree[dname]
            dirname = dep_subtree['dirname']
            os.chdir(os.path.join(old_dir, dirname))

            if verbose:
                sout = 'Directory: {}'.format(dirname)
                print()
                print('  ' + sout)
                print('  ' + '-'*len(sout))
            for dep_key in dep_subtree['files'].keys():
                if verbose:
                    print('  - ' + dep_key)
                fort_dep = dep_subtree['files'][dep_key]
                fort_dep.solve_link_dependencies()


        print_step('Writing "make_deps.mk" files')
        for dname in dep_tree.keys():
            dep_subtree = dep_tree[dname]
            dirname = dep_subtree['dirname']
            os.chdir(os.path.join(old_dir, dirname))

            if verbose:
                sout = 'Directory: {}'.format(dirname)
                print()
                print('  ' + sout)
                print('  ' + '-'*len(sout))

            fname_deps = 'make_deps.mk'
            fout = open(fname_deps, 'w')
            fout.write(('# This file was generated by make_deps.py. '
                        'Do not edit this directly.\n\n'))

            fout.write('\n#Compile dependencies:\n')
            if verbose:
                print('  Writing compile dependencies to ' + fname_deps)
            for dep_key in str_sorted(dep_subtree['files'].keys()):
                fort_dep = dep_subtree['files'][dep_key]
                sout = fort_dep.get_compile_dependency_str()
                if sout is None:
                    continue
                if verbose:
                    print('  - ' + dep_key)
                    #print(sout)
                fout.write(sout)

            fout.write('\n#Link dependencies:\n')
            if verbose:
                print('  Writing link dependencies to ' + fname_deps)
            for dep_key in str_sorted(dep_subtree['files'].keys()):
                fort_dep = dep_subtree['files'][dep_key]
                sout = fort_dep.get_link_dependency_str()
                if sout is None:
                    continue
                if verbose:
                    print('  - ' + dep_key)
                    #print(sout)
                fout.write(sout)

            if verbose:
                print('  File {} writen'.format(fname_deps))
            fout.close()

        print('Ok!')

    def is_dir_valid(self, dirname):
        '''Determines if `dirname` should be considered for automatic dependencies.

        Only directoryes that contain a file common-rules.mk or Makefile, which
        further match against the string "include make_deps.mk", are considered.
        '''
        fnames = glob(os.path.join(dirname, 'common-rules.mk'))
        fnames += glob(os.path.join(dirname, 'Makefile'))
        for fname in fnames:
            with open(fname, 'r') as f:
                data = f.read()
                if self.re_is_dir_valid.search(data) is not None:
                    return True
        return False

    def get_dirs(self, root_dir='.', ignore_patterns=None):
        '''Returns all directories recursively under `root_dir` for the dependency script.

        Recurvisely walk over subdirectories under `root` and returns directories
        matching the criterium in `is_dir_valid`. Also, ignore directories that
        look like version control (*.git, *.svn).
        '''
        if ignore_patterns is None:
            ignore_patterns = ['.git', '.svn']

        matches = []
        for root, _, _ in os.walk(root_dir):
            for pattern in ignore_patterns:
                if pattern in root:
                    break
            else:
                if self.is_dir_valid(root):
                    matches.append(root)
        return matches

    def is_bgw_root_dir(self, dirname):
        '''Determines if `dirname` is the root directory of BerkeleyGW.'''
        for subdir in 'bin', 'Common', 'config':
            cdir = os.path.join(dirname, subdir)
            if not os.path.exists(cdir):
                return False
        return True

    def find_root_directory(self):
        '''Finds the root directory of BerkeleyGW.'''
        cur_dir = os.path.abspath(os.getcwd())
        while True:
            if self.is_bgw_root_dir(cur_dir):
                self.root = cur_dir
                return

            new_dir = os.path.dirname(cur_dir)
            if new_dir==cur_dir:
                raise Exception('Cannot find root directory of BerkeleyGW!')
            cur_dir = new_dir


def main(dirs=None, root=None, verbose=False):
    dep_parser = DependencyParser()
    dep_parser.parse(dirs, root, verbose)


if __name__ == '__main__':
    import argparse
    desc = ('Automatically find module-dependencies for the make files, and '
            'store in "make_deps.mk" files. Default is to find dependencies '
            'all subdirectories, as long as the corresponding file "Makefile" '
            r'or "*.mk" contains a line matching "^include .*make_deps\.mk$".')
    parser = argparse.ArgumentParser(description=desc)
    # FHJ: not working with --root. Also, not necessary.
    #parser.add_argument('dirs', nargs='*', help=(
    #                    'Only work on directories DIRS instead of all '
    #                    'directories under current directory/subdirectories.'))
    parser.add_argument('--root', help=(
                        'Root directory of BerkeleyGW. Default is to '
                        'recursively transverse parent directories until the '
                        'Common subdirectory is found.'))
    parser.add_argument('--verbose', default=False, action='store_true')

    args = parser.parse_args()
    main(**vars(args))
