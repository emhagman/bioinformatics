#!/usr/bin/python

import sys
import re
file_to_convert = ''
if len(sys.argv) != 2:
    print('./fastqtoa.py fastqfile.txt')
    exit(1)
else:
    file_to_convert = sys.argv[1]
file_handle = open(file_to_convert, 'r')
i = 0
output = ''
for line in file_handle:
    if re.match('^@', line) and i == 0:
        output += line.replace('^@', '>')
    elif i == 1:
        output += line
        i = -3
    i += 1
print(output)
