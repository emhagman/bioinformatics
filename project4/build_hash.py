#!/usr/bin/python

from matrix import gen_html_table
import re

# initialize the matrix
aminoacid = {}
submatrix = []

matrix_nums = None
with open('matrix.txt', 'r') as matrix_file:
    i = 0
    for line in matrix_file:
        if i == 0:
            line = line.replace('\n', '')
            matrix_nums = re.split('\s+', line)
            for i in xrange(0, 20):
                aminoacid[matrix_nums[i]] = i
        else:
            line = line.strip('\n')
            temp = re.split('\s+', line)
            submatrix.append(temp)
        i += 1

seq1 = 'EREHSISIVLE'
seq2 = 'QNHKTLGFICN'
len1 = len(seq1)
len2 = len(seq2)

print(aminoacid)
submatrix.insert(0, matrix_nums)
for i, v in enumerate(submatrix):
    submatrix[i].insert(0, matrix_nums[i])
submatrix[0][0] = ''
gen_html_table(submatrix, 'sub.html')

gapscore = -10
matrix = [[0 for y in xrange(len2)] for x in xrange(len1)]
# go through and do the gap value init stufffff
for x, xv in enumerate(matrix):
    for y, yv in enumerate(xv):
        matrix[x][0] = gapscore * x
        matrix[0][y] = gapscore * y
gen_html_table(matrix, 'before.html')

for x in xrange(1, len1):
    for y in xrange(1, len2):
        val1 = matrix[x][y-1] + gapscore
        val2 = matrix[x-1][y] + gapscore
        seq1_nuc = seq1[x-1]
        seq2_nuc = seq2[y-1]
        amino1 = aminoacid[seq1[x-1]]+1
        amino2 = aminoacid[seq2[y-1]]+1
        matrix_val = int(submatrix[amino1][amino2])
        val3 = matrix[x-1][y-1] + matrix_val
        matrix[x][y] = max(val1, val2, val3)
gen_html_table(matrix, 'after.html')

# get strand lengths
dna1_len = len(seq1)
dna2_len = len(seq2)


# score
def score(xind, yind):
    amino1 = aminoacid[seq1[xind-1]]+1
    amino2 = aminoacid[seq2[yind-1]]+1
    matrix_val = int(submatrix[amino1][amino2])
    val3 = matrix_val
    print("AMINO ACID ONE: [%s] :: AMINO ACID 2: [%s] :: Matrix Val: [%d] :: Value 3 : [%d]" % (seq1[xind-1], seq2[yind-1], matrix_val, val3))
    return val3


# check if the values can get us to the path
def got_us_there(x, y, xi, yi):
    if matrix[xi][yi] == matrix[x][y] - score(x, y):
        return True, xi, yi
    elif matrix[xi][yi] == matrix[x][y] - gapscore:
        return True, xi, yi
    else:
        return False, x, y

# lets gooooooooooooooooooooooooo
# start in bottom right corner
x = len(matrix) - 1
y = len(matrix[0]) - 1

# align this son
alignment = []
while x != 1 or y != 1:
    # check diag, cause paul
    match = False
    direction = ""
    if x > 0 and y > 0:
        match, x, y = got_us_there(x, y, x-1, y-1)
        direction = "d"
    if not match and x > 0:
        match, x, y = got_us_there(x, y, x-1, y)
        direction = "u"
    if not match and y > 0:
        match, x, y = got_us_there(x, y, x, y-1)
        direction = "l"
    if not match:
        print('never ever should this happen, ever.')
    else:
        alignment.append((x, y, direction))

# reverse order of our path taken
import os
print(os.linesep + 'Path Taken:')
alignment.reverse()
print(''.join([a[2].upper() for a in alignment]))
