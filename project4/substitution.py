#!/usr/bin/python

from matrix import gen_matrix
from matrix import gen_html_table
import math

# take in two inputs
seq1_frags = ["FRFRFR", "YFRFR", "YFYFR-F"]
seq2_frags = ["ARFRFR", "YF-FR", "YFRFRYF"]

# header amino acids
aminos = ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'W', 'S', 'T', 'C', 'M', 'N', 'Q', 'K', 'R', 'H', 'D', 'E']
matrix = gen_matrix(20, 20, ivalue=1)

# combine the frags to make it easier
seq1 = ""
seq2 = ""
for a in seq1_frags:
    seq1 += a
for a in seq2_frags:
    seq2 += a


def frequency(seqq, seqq2, aa):
    count = 0
    minnn = min(len(seqq), len(seqq2))
    aa = aa.strip()
    for iv in xrange(minnn):
        if seqq[iv].strip() == aa and seqq2[iv].strip() != '-':
            count += 1
        elif seqq2[iv].strip() == aa and seqq[iv].strip() != '-':
            count += 1
        if seqq[iv].strip() == aa and seqq2[iv].strip() == aa and seqq[iv] != '-':
            count += 1
    return count+1

# get alignment count
aligned = {}
minn = min(len(seq1), len(seq2))
for a in xrange(minn):
    if seq1[a] in aligned:
        if seq2[a] in aligned[seq1[a]]:
            if seq1[a] != '-' and seq2[a] != '-':
                aligned[seq1[a]][seq2[a]] += 1
        else:
            if seq1[a] != '-' and seq2[a] != '-':
                aligned[seq1[a]][seq2[a]] = 1
    else:
        if seq1[a] != '-' and seq2[a] != '-':
            aligned[seq1[a]] = {}
            aligned[seq1[a]][seq2[a]] = 1

# get total aligment positions
aligned_pos = 0
for k in aligned.values():
    for a in k.values():
        aligned_pos += a

# generate the new substitution matrix
for i in xrange(0, 20):
    for j in xrange(0, 20):
        if aminos[i] in aligned:
            if aminos[j] in aligned[aminos[i]]:
                qij = float(aligned[aminos[i]][aminos[j]] + 1) / float(aligned_pos)
                pi = float(frequency(seq1, seq2, aminos[i])) / float(aligned_pos*2)
                pj = float(frequency(seq1, seq2, aminos[j])) / float(aligned_pos*2)
                eij = 2*pi*pj if aminos[i] != aminos[j] else pi*pi
                sij = math.log(qij/eij, 2)
                print('checking for: %s and %s' % (aminos[i], aminos[j]))
                print('number of times i and j are aligned: %d' % (float(aligned[aminos[i]][aminos[j]]) + 1))
                print('qij = %f : pi = %f : pj = %f : eij = %f : sij = %f' % (qij, pi, pj, eij, sij))
                matrix[i][j] = sij

# add headers
matrix.insert(0, aminos)
acount = 0
for a in matrix:
    a.insert(0, aminos[acount])
    acount += 1
matrix[0][0] = ''
gen_html_table(matrix, 'gen_matrix.html')