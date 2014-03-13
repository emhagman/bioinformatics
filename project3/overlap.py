#!/usr/bin/python

import random
from prettytable import PrettyTable


def generate_shit_for_shit(coverage_a, start, end):
    for i in range(start, end):
        coverage_a[i] += 1
    return coverage_a


def generate_fragments(strand, mini=2, maxi=4, coverage=10):
    fragments = []
    rough_length = coverage * len(strand)
    actual_length = 0
    point_of_no_return = len(strand) - maxi
    if point_of_no_return <= 0:
        raise Exception('strand is too short for maximum range')
    if mini > maxi or maxi < mini:
        raise Exception('the minimum or maximum values are incorrect')
    while actual_length < rough_length:
        start_pos = random.randint(0, point_of_no_return)
        fragment_len = random.randint(mini, maxi)
        fragment = strand[start_pos:start_pos+fragment_len]
        actual_length += fragment_len
        fragments.append((fragment, start_pos, start_pos+fragment_len))

    # calc coverage
    coverage_arr = [0 for i in range(0, len(strand))]
    for fragment in fragments:
        coverage_arr = generate_shit_for_shit(coverage_arr, fragment[1], fragment[2])
    return fragments, actual_length, rough_length, coverage_arr

# test out some fragments
frags, actual, rough, coverage = generate_fragments("ATGCATATGCATATGCATATGCAT", coverage=3)
tab2 = PrettyTable(['Fragment', 'Fragment Start Position', 'Fragment Length'])
for frag in frags:
    tab2.add_row([frag[0], frag[1], frag[2]])
print(tab2)

# generate overlap
tab = PrettyTable(['Fragment 1', 'Fragment 2', 'Contig', 'Overlap Length'])
def overlap(fragments):
    num_frags = len(fragments)
    for i in range(0, num_frags-1):
        for j in range(i+1, num_frags):
            frag1_len = len(fragments[i])
            frag2_len = len(fragments[j])
            frag_min = min(frag1_len, frag2_len)
            overlapped = 0
            frag1 = fragments[i]
            frag2 = fragments[j]
            k = frag_min - 1
            while k >= 1 and overlapped == 0:
                if frag1[frag1_len-k:frag1_len] == frag2[0:k]:
                    contig = frag1[0:frag1_len-k] + frag2
                    overlapped = k
                    tab.add_row([frag1, frag2, contig, overlapped])
                elif frag2[frag2_len-k:frag2_len] == frag1[0:k]:
                    contig = frag2[0:frag2_len-k] + frag1
                    overlapped = k
                    tab.add_row([frag1, frag2, contig, overlapped])
                k -= 1
overlap([f[0] for f in frags])
print(tab)