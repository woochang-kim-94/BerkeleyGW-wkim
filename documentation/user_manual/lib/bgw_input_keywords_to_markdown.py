from __future__ import print_function
import os
import re
try:
    from cStringIO import StringIO
except:
    from io import StringIO
from mako.template import Template
from mako.exceptions import RichTraceback

class Metadata(object):
    """Metadata descriptor"""

    _type_toggle = 'toggle'
    _type_integer = 'integer'
    _type_float = 'float'
    _type_array = 'array'
    _type_block = 'block'
    _allowed_types = [_type_toggle, _type_integer, _type_float,
                      _type_array, _type_block]

    def __init__(self, line=None):

        self.unit = None
        self._type = None
        self.required = None
        self.default = None
        self.has_default = False

        self.parse(line)

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, value):
        if value not in self._allowed_types:
            raise TypeError('type must be one of the following:\n{}'.format(
                            ', '.join(self._allowed_types)))
        self._type = value

    def parse(self, line):
        """
        Parse a line of the form
        #[unit=Ry,type=float,required=True,default=None]
        """
        if not line:
            return

        line = line.lstrip('# [').rstrip(']')

        for item in line.split(','):
            item = item.strip()
            k, v = item.split('=', 1)
            if not k in ('unit', 'type', 'required', 'default'):
                raise Exception('Unknown metadata: {}'.format(k))
            setattr(self, k, v)

    def update(self, arg, has_pound):
        pass

    def __repr__(self):
        return ('[unit:{}, type:{}, required:{}, default:{}]'.format(
            self.unit, self.type, self.required, self.default))


class Keyword(object):
    """An input keyword"""
    def __init__(self, name, metadata, doc_lines,
                 n_cols=None, col_labels=None, **kwargs):
        self.name = name
        self.metadata = metadata
        self.doc_lines = doc_lines
        self.related_keywords = list()

        # For keywords of block type, enclosed in begin/end statements
        # (eg: "begin qpoint").
        # Each block has `n_cols` columns, with columns labeled according to the
        # array `col_labels`. For instance, a "qpoint" block has 5 columns,
        # each column with a different label:
        # ['qx', 'qy', 'qz', 'factor', 'is_q0']
        self.n_cols = n_cols
        self.col_labels = col_labels

    @property
    def is_block(self):
        return self.metadata.type == self.metadata._type_block

    @property
    def title_str(self):
        return self.name + '\n' + len(self.name) * '-'

    @property
    def metadata_str(self):
        S = str(self.metadata)
        if self.related_keywords:
            S += '\n' + self.related_str
        return S

    @property
    def related_str(self):
        return 'Related to: ' + str(self.related_keywords)

    @property
    def doc_str(self):
        S = ''
        if self.is_block:
            S += 'Columns: ' + ' '.join(self.col_labels) + '\n'
        S += '\n'.join(self.doc_lines)
        return S

    def __eq__(self, other):
        return isinstance(other, Keyword) and self.name==other.name

    def __hash__(self):
        return hash(self.name)

    def __repr__(self):
        return '<\n{}\n{}\n\n{}\n/>'.format(
                self.title_str, self.metadata_str, self.doc_str)


class KeywordGroup(object):
    """
    A keyword group, which consist of documentation on the keyword, plus
    metadata associated with keywords, and keywords. Groups are separated by
    newlines, documentation starts with a pound sign plus a white space,
    metadata starts with a pound sign plus a square bracket, and keywords
    are lines that optionally start with a pound sign, followed by any
    character other than a white space, plus the keyword, plus a sample
    argument. Example:
    # Energy cutoff for the dielectric matrix, in Ry. The dielectric matrix
    # epsilon_{GG`} will contain all G-vectors with kinetic energy |q+G|^2
    # up to this cutoff.
    #[unit=Ry,type=float,required=True,default=None]
    epsilon_cutoff 35.0

    In this case, the keyword is "epsilon_cutoff" and the sample argument is
    "35.0". The metadata for the keyword epsilon_cutoff is:
    "[unit=Ry,type=float,required=True,default=None]".

    The list of metadata associated with each flag is:
    - unit: defaults to None
    - type: see below
    - required: see below
    - default: see below

    Some defaults:
    - Keywords that don't start with a pound sign are automatically marked as
      required, unless explicitly set in the metadata.
    - If a keyword is required, its default is None, unless set in the metadata.
    - If a keyword is not required, its default value is automatically set to
      the sample argument, unless set in the metadata.
    - If the type of the keyword is not specified, it is guessed by the type
      of the sample argument. In this example, "35.0" is interpreted as a float.
    - If there is no sample argument, the type of the keyword is set to
      "toggle".
    """
    def __init__(self, lines):
        self.raw_lines = []
        self.doc_lines = []
        self.keywords = []

        metadata = Metadata()
        keyword_doc_lines = []
        doc_lines = []

        num_empty = 0

        iter_lines = iter(lines)
        while True:
            try:
                #line = iter_lines.next().strip()
                line = next(iter_lines).strip()
            except StopIteration:
                break

            if 'INTERNAL_ONLY' in line:
                continue

            self.raw_lines.append(line)

            has_pound = line.startswith('#')
            line = line.lstrip('#')

            if line:
                num_empty = 0
            else:

                num_empty += 1

                doc_lines.append('')

                # Two or more empty lines separate the group's documentation
                # from the individual keywords documentation.
                if num_empty > 1 and not self.keywords:
                    self.doc_lines.extend(doc_lines)
                    doc_lines = []

                continue

            # What type of line is it?
            is_keyword = line[0] not in (' ', '[')
            is_metadata = line.startswith('[')
            is_doc = has_pound and not (is_keyword or is_metadata)

            if sum([is_keyword, is_metadata, is_doc]) != 1:
                raise Exception('Invalid line:\n{}'.format(line))

            if is_metadata:
                metadata = Metadata(line)

            elif is_doc:
                cur_line = line[1:]
                if cur_line.startswith('### '):
                    cur_line = cur_line + " {:style='margin:0;'}"
                doc_lines.append(cur_line)

            elif is_keyword:

                kwargs = dict()

                metadata.required = not has_pound

                # Use the doc lines collected so far and reset the list
                keyword_doc_lines.extend(doc_lines)
                doc_lines = []

                tokens = line.split()

                # block keyword
                if tokens[0]=='begin':
                    metadata.type = metadata._type_block

                    if len(tokens) < 2:
                        raise Exception(
                            'Block begins without a keyword after lines:\n{}'.format('\n'.join(raw_lines)))

                    name = tokens[1]

                    # TODO: Remove this line form the doc_lines.
                    #       but first we have to make sure we get it right.
                    col_labels = self.raw_lines[0][1:].split()
                    n_cols = len(col_labels)
                    kwargs['n_cols'] = n_cols
                    kwargs['col_labels'] = col_labels

                    # Skip the content of the block
                    while True:
                        #line = iter_lines.next()
                        line = next(iter_lines)
                        self.raw_lines.append(line)
                        if 'end' in line:
                            break
                    
                # simple keyword
                else:
                    name = tokens[0]

                    # Figure out the type
                    if len(tokens) == 1:
                        metadata.type = metadata._type_toggle
                    elif len(tokens) == 2:
                        value = tokens[1]
                        if is_integer(value):
                            metadata.type = metadata._type_integer
                        elif is_float(value):
                            metadata.type = metadata._type_float
                    else:
                        metadata.type = metadata._type_array

                    # TODO figure out default value

                # Add the keyword
                self.keywords.append(Keyword(name, metadata, keyword_doc_lines, **kwargs))

                # Reset the metadata and doc lines for the next keyword
                metadata = Metadata()
                keyword_doc_lines = []

        # Tag this group of keywords as related
        for i, keyword1 in enumerate(self.keywords):
            for keyword2 in self.keywords[i+1:]:
                keyword1.related_keywords.append(keyword2.name)
                keyword2.related_keywords.append(keyword1.name)
        # FIXME We assume that each keyword is present only once in the file,
        #       which is actually not the case (e.g. broadening)...

        # Add any remaining lines to the doc
        self.doc_lines.extend(doc_lines)


    def __repr__(self):
        lines = ['<Group']
        lines.extend(self.doc_lines)
        for keyword in self.keywords:
            lines.append(str(keyword))
        lines.append('/>')
        return '\n'.join(lines)


def is_integer(S):
    """Check whether a string represents an integer"""
    return S.replace('-','').isdigit()


def is_float(S):
    """Check whether a string represents a float"""
    return S.replace('.','').replace('d','').replace('-','').isdigit()


class BGWInputParser:
    """
    Class that parses a BerkeleyGW input file and return a list
    of keyword groups
    """
    def __init__(self, f):
        '''Given a file `f`, return a list of instances of KeywordGroup'''
        self.groups = list()
        self.keywords = list()
        self.unsorted_lines = list()

        for block in self.iter_blocks(f):

            group = KeywordGroup(block)

            self.unsorted_lines.extend(group.doc_lines)
            self.keywords.extend(group.keywords)

            if not group.keywords:
                continue
            self.groups.append(group)

    @staticmethod
    def iter_blocks(f):
        """Iterate over blocks of non-empty lines."""
        block = list()
        for line in f.readlines():
            line = line.strip()
            if line:
                block.append(line)
            elif block:
                yield block
                block = list()
        if block:
            yield block
        return

    def print_groups(self):
        for group in self.groups:
            print(group)

    def print_keywords(self):
        for keyword in self.keywords:
            print(keyword)

    def print_groups_or_keywords(self):
        for group in self.groups:
            if len(group.keywords) > 1:
                print(group)
            else:
                for keyword in group.keywords:
                    print(keyword)

    #def print_unsorted_lines(self):
    #    if self.unsorted_lines:
    #        print('\n'.join([14*'=','Unsorted lines', 14*'=']))
    #        print('\n'.join(self.unsorted_lines))



class MarkdownParserRenderer(object):
    """
    Create a markdown file (.md) from the input files (.inp) using a parser.
    """
    def __init__(self, parser, title):
        self.parser = parser
        self.title = title
        self.re_pounds = re.compile(r' (#+)')
        self.re_auto_link = re.compile(r'\[\[([a-zA-Z_]+)\]\]')
        self.group_keywords()

    @classmethod
    def from_filtered_file(cls, fname, title):
        '''Creates a MarkdownParserRenderer instance by opening and parsing
        input file `fname`, removing INTERNAL_ONLY blocks, and recursively
        including files.
        '''

        def filter_file(fname_in):
            # GKA: I wanted to avoid using sed because it is machine dependent
            #      and causes problems on my laptop
            tag_start = '#' + 'BEGIN_INTERNAL_ONLY'
            tag_stop = '#' + 'END_INTERNAL_ONLY'

            content = ''
            with open(fname_in, 'r') as f:
                while True:
                    try:
                        #line = f.next()
                        line = next(f)
                    except StopIteration:
                        break

                    if line.startswith(tag_start):
                        while True:
                            try:
                                #line = f.next()
                                line = next(f)
                            except StopIteration:
                                break
                            if line.startswith(tag_stop):
                                break
                        continue

                    content += line

            return content

        def recursive_include(str_in, max_depth=5):
            if max_depth<1:
                return str_in

            re_include = re.compile(r'%include ([a-zA-Z_.]+)')

            def f_include(match):
                s = match.group(1)
                return filter_file(s)

            str_out = re_include.sub(f_include, str_in)
            if str_in==str_out:
                return str_out
            else:
                return recursive_include(str_out, max_depth-1)


        # Because we want to expand '%include' statements within the file
        # we have to do it from inside the subdirectory.

        subdir = os.path.dirname(fname)
        basename = os.path.basename(fname)
        back = os.path.relpath(os.path.curdir, subdir)

        os.chdir(subdir)
        output = recursive_include(filter_file(basename))
        os.chdir(back)

        f = StringIO(output)
        parser = BGWInputParser(f)
        return cls(parser, title)

    def get_kw_repr(self, kw):
        s = '`{}'.format(kw.name)
        if kw.metadata.type=='integer':
            s += ' [integer]'
        elif kw.metadata.type=='float':
            s += ' [float]'
        elif kw.metadata.type=='array':
            s += ' [array of integers]'
        elif kw.metadata.type=='block':
            s = '`begin {} ... end'.format(kw.name)
        elif kw.metadata.type=='toggle':
            pass
        else:
            #raise Exception('Unsupported type: {}'.format(kw.metadata.type))
            pass
        s += '`'
        return s

    def get_doc(self, obj):
        s = self.re_auto_link.sub(r'[`\1`](#\1)',
            '\n'.join(obj.doc_lines))
        s = self.re_pounds.sub(r'\1', s)
        return s

    def group_keywords(self):
        self.groups = []
        self.unstructured_groups = []
        self.structured_groups = []
        self.all_keywords = set(self.parser.keywords)
        self.all_keywords = list(self.all_keywords)
        self.all_keywords.sort(key=lambda kw: kw.name)
        self.required_keywords = []
        self.optional_keywords = []
        for kw in self.all_keywords:
            if kw.metadata.required:
                self.required_keywords.append(kw)
            else:
                self.optional_keywords.append(kw)

        for group in self.parser.groups:
            if len(group.keywords):
                self.groups.append(group)
                if any(['###' in line for line in group.doc_lines]):
                    self.structured_groups.append(group)
                else:
                    self.unstructured_groups.append(group)

    def render(self, fname_out):
        template_fname = os.path.join(os.path.dirname(__file__), 'template.md')
        template = Template(filename=template_fname)
        try:
            s = template.render(renderer=self)
        except:
            traceback = RichTraceback()
            for (filename, lineno, function, line) in traceback.traceback:
                print("File %s, line %s, in %s" % (filename, lineno, function))
                print(line, "\n")
            print("%s: %s" % (str(traceback.error.__class__.__name__), traceback.error))
            raise

        with open(fname_out, 'w') as f:
            f.write(s)
