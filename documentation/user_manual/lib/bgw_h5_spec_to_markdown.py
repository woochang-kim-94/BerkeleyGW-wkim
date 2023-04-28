
from __future__ import print_function
from collections import OrderedDict
import re
import markdown

def h5_spec_to_markdown_simple(inp_fname, md_fname, description):
    """Add a description, then insert the h5.spec file."""

    md_content = str(description).rstrip() + '\n\n'
    md_content += '```\n'

    with open(inp_fname, 'r') as f:
        for line in f.readlines():
            #md_content += '    ' + line
            md_content += line

    md_content += '```'

    with open(md_fname, 'w') as f:
        f.write(md_content)


re_kw = re.compile(r'([^ ]*?): \s*(.*)')

class Node(object):
    EXTENSIONS = [
        'fenced_code',
        'codehilite',
        'footnotes',
        'admonition',
        'attr_list',
        'pymdownx.arithmatex',
        'pymdownx.details',
        'extra',
        'pymdownx.extra']

    def __init__(self, parent=None):
        self.parent = parent
        self.children = []
        self.keywords = OrderedDict()
        self.raw_content = ''
        if self.parent is None:
            self.level = 0
        else:
            self.parent.get_type() in ('root', 'group')
            self.parent.children.append(self)
            self.level = self.parent.level + 1

    def get_type(self):
        if self.level==0:
            return 'root'
        if 'Group' in self.keywords.keys():
            return 'group'
        if 'Dataset' in self.keywords.keys():
            return 'dataset'
        print(self.keywords)
        raise Exception('Invalid node type: {}'.format(str(self)))

    def get_parent_with_level(self, level):
        '''Return parent node of level `level`'''
        if level > self.level:
            raise ValueError('Requested level larger than node level.')
        node = self
        while node.level > level:
            node = node.parent
        return node

    def get_full_path(self):
        path = []
        node = self
        while node.level > 0:
            if node.get_type()=='dataset':
                s = node.keywords['Dataset']
            else:
                s = node.keywords['Group']
            path.insert(0, s)
            node = node.parent
        return path

    def render_markdown(self):
        _type = self.get_type()
        s = ""

        if _type=='group':
            s += "<div class='h5-spec-group'>\n"

        if len(self.keywords):
            s += "<div class='h5-spec-{}-desc'>\n".format(_type)
            for k, v in self.keywords.items():
                if k in ('Group', 'Dataset'):
                    path = self.get_full_path()
                    parent_path = '/'.join(['']+path[:-1]+[''])
                    my_path = path[-1]
                    s += '<p><span class="h5-keyword">{}</span>: '.format(k)
                    s += '<span class="h5-path-parent">{}</span>'.format(parent_path)
                    s += '<span class="h5-path-mine">{}</span></p>\n'.format(my_path)
                else:
                    s2 = ("<span class='h5-keyword'>{}</span>: {}\n".
                        format(k, v.replace(r'\n','\n')))
                    s += markdown.markdown(s2, extensions=self.EXTENSIONS)
            s += '</div>'

        if len(self.children):
            if _type=='group':
                s += "<div class='h5-spec-group-children'>\n"
            for child in self.children:
                s += child.render_markdown()
            if _type=='group':
                s += '</div>\n'

        if _type=='group':
            s += '</div>\n'

        return s

    def __repr__(self):
        s = '{}<node level={} len(children)={} len(keywords)={}/>'.format(
                '  '*self.level, self.level, len(self.children), len(self.keywords))
        return s


    def __str__(self):
        if len(self.children):
            s = '{}<node level={} len(children)={} len(keywords)={}>\n'.format(
                    '  '*self.level, self.level, len(self.children), len(self.keywords))
            for child in self.children:
                s += str(child)
            s += '{}</node>\n'.format('  '*self.level)
        else:
            s = repr(self) + '\n'
        return s


def h5_spec_to_markdown_fancy(inp_fname, md_fname, description):
    """Add a description, then insert the h5.spec file."""

    # Parse spec file
    root = Node()
    last_group = root
    cur_group = None
    last_kw = None
    with open(inp_fname, 'r') as f:
        for line_orig in f.readlines():
            line_orig = line_orig.rstrip()
            line = line_orig.strip('\t')
            if len(line)==0:
                # Got an empty line: create new group
                if cur_group is not None:
                    last_group = cur_group
                cur_group = None
                last_kw = None
                continue
            if line[0]=='#':
                # Ignore comments
                continue

            level = line_orig.count('\t') + 1
            if cur_group is None:
                # Create a new group. First, find the parent.
                parent = last_group.get_parent_with_level(level-1)
                cur_group = Node(parent)

            cur_group.raw_content += line
            match = re_kw.match(line)
            if match is None:
                # We are continuing the definition of a keyword.
                assert last_kw is not None
                cur_group.keywords[last_kw] += " " + line
            else:
                # New keyword
                last_kw = match.group(1)
                text = match.group(2)
                cur_group.keywords[last_kw] = text

    md_content = str(description).rstrip() + '\n\n'
    md_content += root.render_markdown()

    with open(md_fname, 'w') as f:
        f.write(md_content)


def h5_spec_to_markdown(inp_fname, md_fname, description):
    h5_spec_to_markdown_fancy(inp_fname, md_fname, description)
    #h5_spec_to_markdown_simple(inp_fname, md_fname, description)
