#!/usr/bin/python

from matrix import gen_html_table
import re

# initialize the matrix
aminoacid = {}
submatrix = []

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
        val3 = matrix[x-1][y-1] + int(submatrix[aminoacid[seq1[x]]][aminoacid[seq2[y]]])
        matrix[x][y] = max(val1, val2, val3)
gen_html_table(matrix, 'after.html')

