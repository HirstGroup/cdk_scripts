#!/usr/bin/env python

import os
import sys
import time

arg_list = sys.argv

skip = 0

for arg in arg_list:
    if 'dependency' in arg:
        skip += 1
        time.sleep(2)

    if 'parseable' in arg:
        skip += 1
        print(100)

cmd = ' '.join(arg_list[skip+1:])

os.system(cmd)