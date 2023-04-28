#!/usr/bin/env python

# A simple script to quickly enable/disable test files.
# Felipe H. da Jornada <jornada@berkeley.edu> (30 Jan 2012)


from __future__ import print_function
import re
try:
    input = raw_input
except NameError:
    pass


class TestCase(object):
    '''The TestCase object represents a test file stored in disk. It contains
    information about the test title, and whether it's enabled. Supports 
    reading and writing to the file.'''

    enabled_str = {-1:'Unknown', 0:'No', 1:'Yes'}

    def __init__(self, id, file, title='', enabled=-1, auto_read=True):
        '''Create a TestCase Object corresponding to a particular file in disk'''

        self.id = id
        self.file = file
        self.title = title
        self.enabled = enabled
        self.changed = False
        if auto_read:
            self.read_file()

    def read_file(self):
        '''Reads the file to update the "title" and "enabled" fields.
        "enabled" will be -1 if we can't read the state from the file.'''

        f = open(self.file)
        str = f.read(4096)
        res = re.search(r'\s*Test\s*:\s*(.*)', str)
        try:
            self.title = res.group(1).strip()
        except:
            self.title = ''
        res = re.search(r'\s*Enabled\s*:\s*(.*)', str)
        try:
            enabled = res.group(1).lower().strip()
            if enabled=='yes':
                self.enabled = 1
            elif enabled=='no':
                self.enabled = 0
            else:
                self.enabled = -1
        except:
            self.enabled = -1
        f.close()

    def save_state(self):
        '''Update the "enabled" field in the test file.'''

        f = open(self.file)
        str = f.read()
        f.close()
        str = re.sub(r'(\s*Enabled\s*:\s*).*',
            r'\1'+self.enabled_str[self.enabled], str)
        f = open(self.file, 'w')
        f.write(str)
        f.close()

    def set_enabled(self, enabled):
        '''Update the "enabled" field in the TestCase object. Note that no 
        changes will be made if the object was originally in an unkown state.'''

        if self.enabled>-1:
            if self.enabled!=enabled:
                self.enabled = enabled
                self.changed = True

    def get_summary(self):
        '''Returns a tuple of strings with test id, title and enabled.'''

        return '%d'%(self.id), self.title, self.enabled_str[self.enabled]

def print_tests_table(tests):
    '''Prints a table summarizing all test files.'''
    def p_line():
        print('+' + '-'*4 + '+' + '-'*60 + '+' + '-'*9 + '+')
    def p_text(s1,s2,s3):
        print('| %2s | %58s | %7s |'%(s1,s2,s3))
        
    p_line()
    p_text('Id','Title','Enabled')
    p_line()
    for test in tests:
        p_text(*test.get_summary())
    p_line()

def options_loop(ntests):
    while True:
        print('Select an option:')
        print('[e#] enable test #')
        print('[d#] disable test #')
        print('[E] enable all tests')
        print('[D] disable all tests')
        print('[q] quit')
        str = input('> ')
        if len(str)==1:
            if str in 'DEq':
                return str, 0
            else:
                print('Invalid option')
        elif len(str)>1:
            if str[0] in 'de':
                try:
                    id = int(str[1:])
                    if id<1 or id>ntests:
                        raise ValueError()
                    return str[0], id-1
                except:
                    print('Invalid test id')
            else:
                print('Invalid option')
        print('')

def main(test_files):
    '''Load all test files, create the TestCase objects, and run main loop.'''

    ntests = len(test_files)
    tests = []
    id = 0
    for test_file in test_files:
        id += 1
        tests.append(TestCase(id, test_file))

    while True:
        print_tests_table(tests)
        str, id = options_loop(ntests)
        if str=='E':
            for test in tests:
                test.set_enabled(1)
        elif str=='D':
            for test in tests:
                test.set_enabled(0)
        elif str=='e':
            tests[id].set_enabled(1)
        elif str=='d':
            tests[id].set_enabled(0)
        else:
            break

    n_changes = sum([test.changed for test in tests])
    if n_changes>0:
        ans = input('Save %d file(s) before leaving (yes/NO)? '%(n_changes)).lower()
        if ans=='yes' or ans=='y':
            for test in tests:
                test.save_state()
            print('%d test(s) updated.'%(n_changes))
        else:
            print('No changes were made.')

if __name__ == '__main__':
    import glob

    test_files = glob.glob('*/*/*.test') + glob.glob('*/*.test')
    test_files.sort()
    main(test_files)
