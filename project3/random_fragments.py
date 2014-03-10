#!/usr/bin/python
# some random strand
# min # of nucleotides in fragment, 8
# max # of nucleotides in fragment, 12
# random start position, take random # between mix + max
# do this times until coverage is X (# of nucleotides in strand * coverage)

import random


def generate_fragments(strand, mini=9, maxi=12, coverage=10):
    fragments = []
    rough_length = coverage * len(strand)
    actual_length = 0
    point_of_no_return = len(strand) - maxi - 1
    if point_of_no_return <= 0:
        raise Exception('strand is too short for maximum range')
    if mini > maxi or maxi < mini:
        raise Exception('the minimum or maximum values are incorrect')
    while actual_length < rough_length:
        start_pos = random.randint(0, point_of_no_return)
        fragment_len = random.randint(mini, maxi)
        fragment = strand[start_pos:start_pos+fragment_len]
        actual_length += fragment_len
        fragments.append(fragment)
    return fragments, actual_length, rough_length