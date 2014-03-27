#!/usr/bin/python
# some random strand
# min # of nucleotides in fragment, 8
# max # of nucleotides in fragment, 12
# random start position, take random # between mix + max
# do this times until coverage is X (# of nucleotides in strand * coverage)

import random
from prettytable import PrettyTable


def determine_indiv_coverage(coverage_a, start, end):
    for i in range(start, end):
        coverage_a[i] += 1
    return coverage_a


def generate_fragments(strand, mini=9, maxi=12, coverage=10):
    fragments = []
    rough_length = coverage * len(strand)
    actual_length = 0
    point_of_no_return = len(strand) - mini - 1
    if point_of_no_return <= 0:
        raise Exception('strand is too short for maximum range')
    if mini > maxi or maxi < mini:
        raise Exception('the minimum or maximum values are incorrect')
    while actual_length < rough_length:
        start_pos = random.randint(0, point_of_no_return)
        fragment_len = random.randint(mini, maxi)
        if start_pos+fragment_len >= len(strand):
            fragment = strand[start_pos:]
            actual_length += len(strand) - start_pos
            fragment_len = len(strand) - start_pos
        else:
            fragment = strand[start_pos:start_pos+fragment_len]
            actual_length += fragment_len
        fragments.append((fragment, start_pos, start_pos+fragment_len))

    # calc coverage
    coverage_arr = [0 for i in range(0, len(strand))]
    for fragment in fragments:
        coverage_arr = determine_indiv_coverage(coverage_arr, fragment[1], fragment[2])
    return fragments, actual_length, rough_length, coverage_arr

strand = "ATGCATATGCAT"
frags, actual, rough, coverage = generate_fragments(strand)

tab = PrettyTable([x+1 for x in range(0, len(strand))])
tab.add_row([x for x in strand])
tab.add_row(coverage)
print(tab)
print('Sequence Reads: ' + str(len(frags)))
