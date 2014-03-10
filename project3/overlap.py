#!/usr/bin/python

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

# test out some fragments
frags, actual, rough = generate_fragments("ATGCATATGCATATGCATATGCAT")


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
                    print(frag1, frag2, contig, overlapped)
                elif frag2[frag2_len-k:frag2_len] == frag1[0:k]:
                    contig = frag2[0:frag2_len-k] + frag1
                    overlapped = k
                    print(frag1, frag2, contig, overlapped)
                k -= 1

overlap(['AATGCCTGA', 'TGACGAGTTAATGC'])