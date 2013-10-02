# -*- coding: utf-8 -*-
"""
Created on Mon May 27 17:25:26 2013

@author: Jeff
"""

import cPickle
import gzip
import re

read_dict = cPickle.load(open('42_reads_genomes_mapping.p', 'rb'))
gi_tax = cPickle.load(open('bacterial_genomes_taxonomy.p', 'rb'))
gi_id = cPickle.load(open('bacterial_genomes_ncbiid_gi.p', 'rb'))

## establish the number of genomes top 0.1 % and bottom 0.1 % of reads are assigned to
mean = []
for read in read_dict.keys():
    mean.append(len(read_dict[read]))

mean.sort()    
top = mean[int(len(mean) * 0.999)]
bottom = 1

## identify the top 0.1 % and bottom 0.1 %
top_reads = set()
bottom_reads = set()  
for read in read_dict.keys():
    if len(read_dict[read]) >= top:
        top_reads.add(read)
    elif len(read_dict[read]) <= bottom:
        bottom_reads.add(read)

## create fasta of most broadly mapped reads and least broadly mapped reads
with gzip.open('../42_trim_mates_singles.fasta.gz', 'rb') as fasta_in, open('42_bottom_mapped_reads.fasta', 'w') as bottom, open('42_top_mapped_reads.fasta', 'w') as top:          
    for line in fasta_in:
        line = line.rstrip('\n')
        if line.startswith('>'):
            top_keep = False
            bottom_keep = False
            name = line.strip('>')
            name = re.split(' *', name)
            name = name[0]
            if name in top_reads:
                top_keep = True
                print >> top, line
                print name, 'top'
            elif name in bottom_reads:
                print >> bottom, line
                bottom_keep = True
                print name, 'bottom'
        elif top_keep == True:
            print >> top, line
        elif bottom_keep == True:
            print >> bottom, line
            
## tally the number of different taxonomies present in taxonomic levels for
## most broadly mapped reads 
top_output = open('42_taxonomic_tallies_top_reads.txt', 'w')
error = open('bad_gi.txt', 'w')
bad_gis = set()  
for read in top_reads:
    genus_set = set()
    family_set = set()
    order_set = set()
    class_set = set()
    phylum_set = set()
    ncbis = read_dict[read]
    for ncbi in ncbis:
        try:
            gi = gi_id[ncbi]
            tax = gi_tax[gi]
            genus_set.add(tax[5])
            family_set.add(tax[4])
            order_set.add(tax[3])
            class_set.add(tax[2])
            phylum_set.add(tax[1])
        except KeyError:
            print gi
            print >> error, gi
            bad_gis.add(gi)
    print read, len(genus_set), len(family_set), len(order_set), len(class_set), len(phylum_set)
    print >> top_output, read, len(genus_set), len(family_set), len(order_set), len(class_set), len(phylum_set)
top_output.close()

## bottom should all have only one map.  return the taxonomies.
bottom_output = open('42_taxonomic_tallies_bottom_reads.txt', 'w')  
for read in bottom_reads:
    ncbis = read_dict[read]
    for ncbi in ncbis:
        try:
            gi = gi_id[ncbi]
            tax = gi_tax[gi]
            species = tax[6]
            genus = tax[5]
            family = tax[4]
            order = tax[3]
            clss = tax[2]
            phylum = tax[1]
        except KeyError:
            print gi
            print >> error, gi
            bad_gis.add(gi)
    print read, species, genus, family, order, clss, phylum
    print >> bottom_output, read,'\t',species,'\t',genus,'\t',family,'\t',order,'\t',clss,'\t',phylum
bottom_output.close()
error.close()