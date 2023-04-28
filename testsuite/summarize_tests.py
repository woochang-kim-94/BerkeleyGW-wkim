#!/usr/bin/env python

# Felipe H. da Jornad (2019).

from __future__ import print_function
import sys
import re
import tabulate


if len(sys.argv)!=2:
	sys.exit('Usage: %s test_result'%(sys.argv[0]))

fname = sys.argv[1]

results = []

kw_parsers = {
        'test': re.compile(r'\*\*\*\*\* (.*) \*\*\*\*\*'),
        'input': re.compile(r'Using input file : (.*)'),
        'app': re.compile(r'Executing:.*/bin/([^ ]*) '),
        'elapsed': re.compile(r'Elapsed time:\s*([0-9.]*)'),
        }

keys = list(kw_parsers.keys())
last_kws = dict.fromkeys(keys)

with open(fname, 'r') as f:
    for line in f.readlines():

        for k in kw_parsers.keys():
            m = kw_parsers[k].search(line)

            if m is not None:
                val = m.group(1)

                if k=='elapsed':
                    val = float(val)

                last_kws[k] = val

                if k=='elapsed':
                    results.append(dict(last_kws))

results.sort(key=lambda d: d['elapsed'])
#print(results)
table = [[result[k] for k in keys] for result in results]
t = tabulate.tabulate(table, headers=keys)
print(t)
