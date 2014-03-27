#!/usr/bin/python

from subprocess import Popen as call
from subprocess import PIPE
import re
import sqlite3
sqlite_conn = sqlite3.connect('chrII.db')
sc = sqlite_conn.cursor()

# create table
#sc.execute('CREATE TABLE vcf (pos INT, ref text, alt text, qual int)')
#sqlite_conn.commit()


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

    target_genes = init_gene_database()
    print('Loading target genes into database...')
    print('')

    chromosome_num = raw_input('Chromosome (#): ')
    chromosome_num = int(chromosome_num)
    print('Scanning genome... may take awhile.')
    print('')

    genome = get_file("full.fa")
    genome = ''.join([line.strip() for line in genome])
    gene = raw_input('Gene of interest: ').lower()

    # for testing
    if gene == '':
        gene = 'F43E2.6'

    # get the gene data
    gene_data = target_genes[gene]

    # get start pos
    find_chrom_start = re.compile(get_chromosome_header(chromosome_num))
    start = find_chrom_start.search(genome).end()

    # get end pos
    find_chrom_end = re.compile(get_chromosome_header(chromosome_num+1))
    end = find_chrom_end.search(genome).start()

    print('')
    print('Chromosome %s Start: %d' % (roman(chromosome_num), start))
    print('Chromosome %s End: %d' % (roman(chromosome_num), end))
    print('Cutting from FASTA...')
    print('')

    # get the chromsome
    chromosome = genome[start:end]
    c_out = open('Chromosome_' + roman(chromosome_num) + '.fa', 'w')
    c_out.write(get_chromosome_header(chromosome_num) + '\n')
    c_out.write(chromosome)
    c_out.close()

    # set the positions 100k back and front of
    start1 = gene_data['start'] - 100000
    end1 = gene_data['start']
    start2 = gene_data['end']
    end2 = gene_data['end'] + 100000

    # print
    chromosome = genome[start:end]
    c_out = open('gene_sequence.txt', 'w')
    c_out.write(chromosome[end1:start2])
    c_out.close()

    print('Searching -100K Back Starting: %d' % (start1,))
    print('Searching -100K Back Ending: %d' % (end1,))
    print('Searching +100K Forward Starting: %d' % (start2,))
    print('Searching +100K Forward Ending: %d' % (end2,))
    print('')

    # restriction site
    all_restrictions = get_file_array('rebaseN.txt')
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
    for a, v in enumerate(chromosome):
        if start1 <= a < end1:
            output1 += chromosome[a]
    output2 = ""
    for a, v2 in enumerate(chromosome):
        if start2 <= a < end2:
            output2 += chromosome[a]

    # get positions of our polymorphisms
    restriction_sites = []
    for restriction in restrictions:
        find_pos = re.compile(restriction[1])
        for m in find_pos.finditer(output1):
            restriction_sites.append((start1 + m.start(), start1 + m.end(), restriction[0], restriction[1]))
        for m in find_pos.finditer(output2):
            restriction_sites.append((start2 + m.start(), start2 + m.end(), restriction[0], restriction[1]))

    # include the gene
    restriction_sites.append((end1, start2, 'Gene', ''))
    print('Number of Restriction Sites Found: %d' % (len(restriction_sites),))

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
    #sc.execute('DELETE FROM vcf WHERE pos < ?', (start1,))
    #sc.execute('DELETE FROM vcf WHERE pos > ?', (end2,))
    #sqlite_conn.commit()
    #print('Done creating database.')

    # now that it is cut down, let's go to work
    polymorphisms = []
    for se in restriction_sites:
        for s in sc.execute('SELECT * FROM vcf WHERE pos >= ? AND pos <= ?', (se[0], se[1])):
            polymorphisms.append((se[2], se[3], s))

    print('Done finding %s +/- 100K polymorphisms.' % (gene,))
    print('')
    print('Searching for creation of restriction sites...')

    # part 2 poly
    kb = [s for s in sc.execute('SELECT * FROM vcf WHERE pos >= ? AND pos <= ?', (start1, end1))]
    kf = [s for s in sc.execute('SELECT * FROM vcf WHERE pos >= ? AND pos <= ?', (start2, end2))]
    kb.extend(kf)
    print('Total # of Polymorphisms In +/- 100K Range: %d' % (len(kb),))

    # +/- 100k all polymorphisms
    create_restrictions = []
    for p in kb:
        for restriction in restrictions:
            # get mutated substring
            pos = p[0]

            # get each one
            if ',' in p[2]:
                muts = p[2].split(',')
            else:
                muts = p[2]

            # check all possible ones
            if isinstance(muts, list):
                for mu in muts:
                    check = chromosome[pos - 40:pos] + mu + chromosome[pos + 1:pos + 40]
                    #print('----------------------------------------*----------------------------------------')
                    #print(chromosome[pos - 40:pos + 40])
                    #print(check)
                    #print('----------------------------------------' + mu)
                    #print('-----')
                    find_pos = re.compile(restriction[1].lower())
                    for m in find_pos.finditer(check):
                        start = (p[0] - 40) + m.start()
                        end = (p[0] - 40) + m.end()
                        if start <= p[0] <= end:
                            create_restrictions.append((restriction[0], p[0], p[1].lower(), p[2].lower()))
            else:
                check = chromosome[pos - 40:pos] + muts + chromosome[pos + 1:pos + 40]
                find_pos = re.compile(restriction[1].lower())
                for m in find_pos.finditer(check):
                    start = (p[0] - 40) + m.start()
                    end = (p[0] - 40) + m.end()
                    if start <= p[0] <= end:
                        create_restrictions.append((restriction[0], p[0], p[1].lower(), p[2].lower()))

    # print it out
    out = open('results.txt', 'w')
    for p in polymorphisms:
        if p[0] == 'Gene':
            out.write('%s\t%s\t%s\t%s\t%s\n' % ('GENE', gene, p[2][0], p[2][1], p[2][2]))
        else:
            out.write('%s\t%s\t%s\t%s\t%s\n' % ('DESTROY', p[0], p[2][0], p[2][1], p[2][2]))

    for r in create_restrictions:
        out.write('%s\t%s\t%s\t%s\t%s\n' % ('CREATE', r[0], r[1], r[2], r[3]))
    out.close()


if __name__ == '__main__':
    main()