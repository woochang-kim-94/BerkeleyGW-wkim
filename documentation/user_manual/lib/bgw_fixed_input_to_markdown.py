
def fixed_keyword_to_markdown(inp_fname, md_fname, description):
    """Add a description, then insert the indented inp file."""

    md_content = str(description).rstrip() + '\n\n'

    with open(inp_fname, 'r') as f:
        for line in f.readlines():
            md_content += '    ' + line

    with open(md_fname, 'w') as f:
        f.write(md_content)

