#!/usr/bin/env python
"""
Assemble the manual for mkdocs.
"""

import os
import shutil
import sys
import glob
sys.path.append('lib')

from bgw_h5_spec_to_markdown import h5_spec_to_markdown
from bgw_input_keywords_to_markdown import MarkdownParserRenderer
from bgw_fixed_input_to_markdown import fixed_keyword_to_markdown

# Directories where the manual will be built
website_dir = 'website'
mkdocs_docs_dir = os.path.join(website_dir, 'docs')
mkdocs_stylesheets_dir = os.path.join(mkdocs_docs_dir, 'stylesheets')


# Create the main directories
if not os.path.exists(website_dir):
    os.makedirs(website_dir)

if not os.path.exists(mkdocs_docs_dir):
    os.makedirs(mkdocs_docs_dir)
else:
    # Cleanup docs directory
    shutil.rmtree(mkdocs_docs_dir)

if not os.path.exists(mkdocs_stylesheets_dir):
    os.makedirs(mkdocs_stylesheets_dir)


def relative_link(src, dst):
    '''Make relative symbolic link pointing to `src` at `dst`.

    If src is 'A/B/C.txt' and dst is 'A/D/F.txt', this is equivalent to
        ln -s ../../B/C.txt A/D/F.txt
    '''
    if os.path.isdir(dst):
        # Make link relative to directory `dst`
        src_rel = os.path.relpath(src, dst)
    else:
        # Make link relative to directory containing `dst`
        src_rel = os.path.relpath(src, os.path.dirname(dst))
    if os.path.lexists(dst):
        if not os.path.samefile(src, dst):
            os.unlink(dst)
            os.symlink(src_rel, dst)
    else:
        os.symlink(src_rel, dst)


def filtered_copy(src, dst):
    '''Copy file src to dst, removing INTERNAL_ONLY flags'''

    lines = []
    internal = False
    with open(src, 'r') as f_in:
        for line in f_in:
            if line.startswith('#' + 'BEGIN_INTERNAL_ONLY'):
                internal = True
            elif line.startswith('#' + 'END_INTERNAL_ONLY'):
                internal = False
            elif not internal:
                lines += [line]
    with open(dst, 'w') as f_out:
        f_out.write(''.join(lines))


# Copy stylesheets
dirname = 'stylesheets'
for basename in os.listdir(dirname):
    src = os.path.join(dirname, basename)
    dst = os.path.join(mkdocs_stylesheets_dir, basename)
    relative_link(src, dst)


# Copy markdown files
copy_directory_content = ['img']
for dirname in copy_directory_content:
    for basename in os.listdir(dirname):
        src = os.path.join(dirname, basename)
        dst = os.path.join(mkdocs_docs_dir, basename)
        relative_link(src, dst)

copy_directory_content = ['misc', 'Overviews']
for dirname in copy_directory_content:
    for basename in os.listdir(dirname):
        src = os.path.join(dirname, basename)
        dst = os.path.join(mkdocs_docs_dir, basename)
        filtered_copy(src, dst)


# Generate keyword markdown files
inputs_markdown_title = [
    ('epsilon.inp', 'epsilon-keywords.md', '`Epsilon` code input keywords (`epsilon.inp`)'),
    ('sigma.inp', 'sigma-keywords.md', '`Sigma` code input keywords (`sigma.inp`)'),
    ('kernel.inp', 'kernel-keywords.md', '`Kernel` code input keywords (`kernel.inp`)'),
    ('absorption.inp', 'absorption-keywords.md', '`Absorption` code input keywords (`absorption.inp`)'),
    ('inteqp.inp', 'inteqp-keywords.md', '`inteqp` code input keywords (`inteqp.inp`)'),
    ('nonlinearoptics.inp', 'nonlinearoptics-keywords.md', '`nonlinearoptics` code input keywords (`nonlinearoptics.inp`)'),
    ('plotxct.inp', 'plotxct-keywords.md', '`plotxct` code input keywords (`plotxct.inp`)'),
    ('parabands.inp', 'parabands-keywords.md', '`ParaBands` code input keywords (`parabands.inp`)'),
    ('bgw2sgw.inp', 'bgw2sgw-keywords.md', '`bgw2sgw` input keywords (`bgw2sgw.inp`)')
    ]

for inp, md, title in inputs_markdown_title:

    inp_fname = os.path.join('../input_files', inp)
    out_fname = os.path.join(mkdocs_docs_dir, md)

    renderer = MarkdownParserRenderer.from_filtered_file(inp_fname, title)
    renderer.render(out_fname)


# Generate fixed input description files
fixed_inputs_markdown_desc = [
    ('abi2bgw.in', 'abi2bgw-input.md', '`abi2bgw` code input format (`abi2bgw.in`)'),
    ('epm.inp', 'epm-input.md', '`EPM` code input format (`epm.inp`)'),
    ('pw2bgw.inp', 'pw2bgw-input.md', '`pw2bgw` code input format (`pw2bgw.inp`)'),
    ('kgrid.inp', 'kgrid-input.md', '`kgrid` code input format (`kgrid.inp`)'),
    ('sig2wan.inp', 'sig2wan-input.md', '`sig2wan` code input format (`sig2wan.inp`)'),
    ('gsphere.inp', 'gsphere-input.md', '`gsphere` code input format (`gsphere.inp`)'),
    ('summarize_eigenvectors.inp', 'summarize_eigenvectors-input.md', '`summarize_eigenvectors` code input format (`summarize_eigenvectors.inp`)'),
    ('epsconv.inp', 'epsconv-input.md', '`epsascbin` or `epsbinasc` code input format (`epsconv.inp`)'),
    ('epsinvomega.inp', 'epsinvomega-input.md', '`epsinvomega` code input format (`epsinvomega.inp`)'),
    ('epsmat_intp.inp', 'epsmat_intp-input.md', '`epsmat_intp` code input format (`epsmat_intp.inp`)'),
    ('epsmat_merge.inp', 'epsmat_merge-input.md', '`epsmat_merge` code input format (`epsmat_merge.inp`)'),
    ('epsmat_old2hdf5.inp', 'epsmat_old2hdf5-input.md', '`epsmat_old2hdf5` code input format (`epsmat_old2hdf5.inp`)'),
    ('epsomega.inp', 'epsomega-input.md', '`epsomega` code input format (`epsomega.inp`)'),
    ]

for inp, md, desc in fixed_inputs_markdown_desc:

    inp_fname = os.path.join('../input_files', inp)
    out_fname = os.path.join(mkdocs_docs_dir, md)

    fixed_keyword_to_markdown(inp_fname, out_fname, desc)


# Generate fixed format for the h5 file formats
for inp_fname in glob.glob('../file_formats/*.h5.spec'):

    bname = os.path.basename(inp_fname).split('.h5.spec')[0]
    out_fname = os.path.join(mkdocs_docs_dir, bname + '_h5_spec.md')
    desc = '#`{}.h5` file format'.format(bname)

    h5_spec_to_markdown(inp_fname, out_fname, desc)
