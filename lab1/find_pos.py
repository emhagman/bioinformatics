#!/usr/bin/python

from subprocess import Popen as call
from subprocess import PIPE


def call_program(arr):
    output = call(arr, stdout=PIPE).communicate()[0]
    return output


def get_file(fname):
    f = open(fname, 'r')
    out = f.read()
    f.close()
    return out

chromsome2 = get_file("c2.fa")
chromsome2 = ''.join([line.strip() for line in chromsome2])

end1 = 7366591
start1 = end1 - 100000
start2 = 7368846
end2 = start2 + 100000

output1 = ""
for a, v in enumerate(chromsome2):
    if a >= start1 and a <= end1:
        output1 += chromsome2[a]

output2 = ""
for a, v2 in enumerate(chromsome2):
    if a >= start2 and a <= end2:
        output2 += chromsome2[a]

import re
find_pos = re.compile('agatct')
for m in find_pos.finditer(output1):
    print m.start(), m.end()
find_pos = re.compile('agatct')
for m in find_pos.finditer(output2):
    print m.start(), m.end()