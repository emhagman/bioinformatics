#!/usr/bin/python

from subprocess import Popen as call
from subprocess import PIPE
import re
import sqlite3


def roman(vv):
    if vv == 1:
        return 'I'
    if vv == 2:
        return 'II'
    if vv == 3:
        return 'III'
    if vv == 4:
        return 'IV'
    if vv == 5:
        return 'V'
    if vv == 10 or vv == 6:
        return 'X'


def get_chromosome_header(chrr):
    return '>CHROMOSOME_' + roman(chrr)


def get_vcf_filename(chrr):
    return 'HwVarVcf_Chr' + roman(chrr) + '.vcf'


def call_program(arr):
    output = call(arr, stdout=PIPE).communicate()[0]
    return output


def get_file(fname):
    f = open(fname, 'r')
    out = f.read()
    f.close()
    return out


def get_file_array(fname):
    f = open(fname, 'r')
    lines = []
    for l in f:
        lines.append(l)
    f.close()
    return lines


def init_gene_database():
    target_genes = {}
    gene_f = get_file('target_genes.txt')
    lines = gene_f.split('\n')
    count = 0
    for line in lines:
        if count != 0:
            columns = re.split('\s+', line)
            target_genes[columns[0].strip()] = {
                'chromosome': columns[1].strip(),
                'start': int(columns[2].strip()),
                'end': int(columns[3].strip())
            }
        count += 1
    return target_genes


def main():

    chromosome_num = raw_input('Chromosome (#): ')
    chromosome_num = int(chromosome_num)

    # check for create
    import os
    create_db = False
    db_file = 'chr' + roman(chromosome_num) + '.db'
    if not os.path.isfile(db_file):
        create_db = True

    sqlite_conn = sqlite3.connect(db_file)
    sc = sqlite_conn.cursor()

    if create_db:
        sc.execute('CREATE TABLE vcf (pos INT, ref text, alt text, qual int)')
        sqlite_conn.commit()

    # restriction site
    all_restrictions = get_file_array('rebaseN2.txt')
    restrictions = []
    for a in all_restrictions:
        cols = re.split('\s+', a)
        cols[0] = cols[0].strip()
        cols[1] = cols[1].replace('^', '').lower().strip()
        if cols[1] != '':
            restrictions.append((cols[0], cols[1]))
    print('Number of Enzymes: %d' % (len(all_restrictions),))

    # combine this into the two places we want to look
    output1 = ""

    # get positions of our polymorphisms
    restriction_sites = []
    for restriction in restrictions:
        find_pos = re.compile(restriction[1])
        for m in find_pos.finditer(output1):
            restriction_sites.append((m.start(), m.end(), restriction[0], restriction[1]))

    # include the gene
    print('Number of Restriction Sites Found: %d' % (len(restriction_sites),))

    # let's check this against our polymorph file
    if create_db:
        polymorph_file = get_file_array(get_vcf_filename(chromosome_num))
        print('Creating sqlite db for quick access...')
        for l in polymorph_file:
            line = re.split('\s+', l)
            line[1] = line[1].strip()
            line[3] = line[3].strip()
            line[4] = line[4].strip()
            line[5] = line[5].strip()
            into_db = (line[1], line[3], line[4], line[5])
            if line[4] != '.':
                sc.execute('INSERT INTO vcf VALUES (?, ?, ?, ?)', into_db)
        sqlite_conn.commit()
        print('Done creating database.')

    # now that it is cut down, let's go to work
    polymorphisms = []
    for se in restriction_sites:
        for s in sc.execute('SELECT * FROM vcf WHERE pos >= ? AND pos <= ?', (se[0], se[1])):
            polymorphisms.append((se[2], se[3], s))

    print('Searching for creation of restriction sites...')

    # part 2 poly
    kb = [s for s in sc.execute('SELECT * FROM vcf ')]
    print('Total # of Polymorphisms In +/- 100K Range: %d' % (len(kb),))

    # +/- 100k all polymorphisms
    create_restrictions = []
    for p in kb:
        for restriction in restrictions:

            # search width
            search_width = len(restriction[1])

            # get mutated substring pos
            pos = p[0]

            # get each one
            if ',' in p[2]:
                muts = p[2].split(',')
            else:
                muts = p[2]

            # check all possible ones
            if isinstance(muts, list):
                for mu in muts:
                    check = output1[pos - search_width:pos] + mu + output1[pos + len(mu):pos + search_width]
                    #print('----------------------------------------*----------------------------------------')
                    #print(chromosome[pos - search_width:pos + search_width])
                    #print(check)
                    #print('----------------------------------------' + mu)
                    #print('-----')
                    find_pos = re.compile(restriction[1].lower())
                    for m in find_pos.finditer(check):
                        start = (p[0] - search_width) + m.start()
                        end = (p[0] - search_width) + m.end()
                        if start <= p[0] <= end:
                            create_restrictions.append((restriction[0], p[0], p[1].lower(), p[2].lower()))
            else:
                check = output1[pos - search_width:pos] + muts + output1[pos + len(muts):pos + search_width]
                find_pos = re.compile(restriction[1].lower())
                for m in find_pos.finditer(check):
                    start = (p[0] - search_width) + m.start()
                    end = (p[0] - search_width) + m.end()
                    if start <= p[0] <= end:
                        create_restrictions.append((restriction[0], p[0], p[1].lower(), p[2].lower()))

    # print it out
    out = open('test_results.txt', 'w')
    for p in polymorphisms:
        out.write('%s\t%s\t%s\t%s\t%s\n' % ('DESTROY', p[0], p[2][0], p[2][1], p[2][2]))

    for r in create_restrictions:
        out.write('%s\t%s\t%s\t%s\t%s\n' % ('CREATE', r[0], r[1], r[2], r[3]))
    out.close()


if __name__ == '__main__':
    main()