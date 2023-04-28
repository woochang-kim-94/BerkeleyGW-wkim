#!/usr/bin/env python

import sys
import base64

with open(sys.argv[1], 'rb') as f:
    s = base64.b64encode(f.read())

print s
