# BerkeleyGW input file specification

All input files, with documentation, should be put here!

These files are processed by the `assemble_manual.py` script under the
`../user_manual/` folder.  If you add any new input file here, make sure to
change the `assemble_manual.py` script and put the script in the list
`inputs_markdown_title`, which should hold input files that get processed by
`inread`-type of parsers, or `fixed_inputs_markdown_desc`, which simply copies
the content of the input file without putting any formatting.
