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
                'name': columns[0].strip(),
                'chromosome': columns[1].strip(),
                'start': int(columns[2].strip()),
                'end': int(columns[3].strip())
            }
        count += 1
    return target_genes


def main():

    rrange = 100000
    target_genes = init_gene_database()
    print('Loading target genes into database...')
    print('')

    chromosome_num = raw_input('Chromosome (#): ')
    chromosome_num = int(chromosome_num)

    # check for create
    import os
    db_file = 'chr' + roman(chromosome_num) + '.db'
    sqlite_conn = sqlite3.connect(db_file)
    sc = sqlite_conn.cursor()
    #sc.execute('CREATE TABLE vcf (pos INT, ref text, alt text, qual int)')
    #sqlite_conn.commit()

    # for testing
    genes = []
    for k in target_genes:
        if target_genes[k]['chromosome'] == 'CHROMOSOME_' + roman(chromosome_num):
            genes.append(target_genes[k])

    # get the chromsome
    chr_file = 'chromosome_' + roman(chromosome_num) + '.fa'
    if os.path.isfile(chr_file):
        print('Using pre-existing chromosome file...')
        print('')
        chromosome = get_file(chr_file)
    else:
        print('Scanning genome... may take awhile.')
        print('')
        genome = get_file("full.fa")
        genome = ''.join([line.strip() for line in genome])
        find_chrom_start = re.compile(get_chromosome_header(chromosome_num))
        start = find_chrom_start.search(genome).end()
        find_chrom_end = re.compile(get_chromosome_header(chromosome_num+1))
        end = find_chrom_end.search(genome).start()
        chromosome = genome[start:end]
        print('')
        print('Chromosome %s Start: %d' % (roman(chromosome_num), start))
        print('Chromosome %s End: %d' % (roman(chromosome_num), end))
        print('Cutting from FASTA...')
        print('')
        c_out = open(chr_file, 'w')
        c_out.write(get_chromosome_header(chromosome_num) + '\n')
        c_out.write(chromosome)
        c_out.close()

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

    # let's check this against our polymorph file
    #polymorph_file = get_file_array(get_vcf_filename(chromosome_num))
    #print('Creating sqlite db for quick access...')
    #for l in polymorph_file:
    #    line = re.split('\s+', l)
    #    line[1] = line[1].strip()
    #    line[3] = line[3].strip()
    #    line[4] = line[4].strip()
    #    line[5] = line[5].strip()
    #    into_db = (line[1], line[3], line[4], line[5])
    #    if line[4] != '.':
    #        sc.execute('INSERT INTO vcf VALUES (?, ?, ?, ?)', into_db)
    #sqlite_conn.commit()
    #sqlite_conn.close()
    #print('Done creating database.')

    # let's create a database for each gene, speeds up things significantly
    import shutil
    for g in genes:
        g['database'] = 'databases/' + 'chr' + roman(chromosome_num) + g['name'] + '.db'
    #    shutil.copyfile(db_file, g['database'])

    # set the positions 100k back and front of
    for g in genes:

        # connecto to gene specific db
        gene_conn = sqlite3.connect(g['database'])
        sc_gene = gene_conn.cursor()

        gene_data = g
        start1 = gene_data['start'] - rrange
        end1 = gene_data['start']
        start2 = gene_data['end']
        end2 = gene_data['end'] + rrange

        sc_gene.execute('DELETE FROM vcf WHERE pos <= ? OR pos >= ?', (start1, end2))
        gene_conn.commit()

        # combine this into the two places we want to look
        output1 = ""
        for a, v in enumerate(chromosome):
            if start1 <= a < end2:
                output1 += chromosome[a]

        # get positions of our polymorphisms
        restriction_sites = []
        for restriction in restrictions:
            find_pos = re.compile(restriction[1])
            for m in find_pos.finditer(output1):
                restriction_sites.append((start1 + m.start(), start1 + m.end(), restriction[0], restriction[1]))

        # include the gene
        restriction_sites.append((end1, start2, 'Gene', ''))
        print('Gene: %s' % (g['name']))
        print('Number of Restriction Sites Found: %d' % (len(restriction_sites),))

        # now that it is cut down, let's go to work
        polymorphisms = []
        for se in restriction_sites:
            for s in sc_gene.execute('SELECT * FROM vcf WHERE pos >= ? AND pos <= ?', (se[0], se[1])):
                polymorphisms.append((se[2], se[3], s))

        print('Done finding %s +/- 100K polymorphisms.' % (g['name'],))
        print('Searching for creation of restriction sites...')

        # part 2 poly
        kb = [s for s in sc_gene.execute('SELECT * FROM vcf WHERE pos >= ? AND pos <= ?', (start1, end1))]
        kf = [s for s in sc_gene.execute('SELECT * FROM vcf WHERE pos >= ? AND pos <= ?', (start2, end2))]
        kb.extend(kf)
        print('Total # of Polymorphisms In +/- 100K Range: %d' % (len(kb),))
        gene_conn.close()

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
                        check = chromosome[pos - search_width:pos] + mu + chromosome[pos + len(mu):pos + search_width]
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
                    check = chromosome[pos - search_width:pos] + muts + chromosome[pos + len(muts):pos + search_width]
                    find_pos = re.compile(restriction[1].lower())
                    for m in find_pos.finditer(check):
                        start = (p[0] - search_width) + m.start()
                        end = (p[0] - search_width) + m.end()
                        if start <= p[0] <= end:
                            create_restrictions.append((restriction[0], p[0], p[1].lower(), p[2].lower()))

        print('Done with gene.')
        print('')

        # print it out
        out = open('genes/' + g['name'] + '.txt', 'w')
        for p in polymorphisms:
            if p[0] == 'Gene':
                out.write('%s\t%s\t%s\t%s\t%s\n' % ('GENE', g['name'], p[2][0], p[2][1], p[2][2]))
            else:
                out.write('%s\t%s\t%s\t%s\t%s\n' % ('DESTROY', p[0], p[2][0], p[2][1], p[2][2]))

        for r in create_restrictions:
            out.write('%s\t%s\t%s\t%s\t%s\n' % ('CREATE', r[0], r[1], r[2], r[3]))
        out.close()

if __name__ == '__main__':
    main()