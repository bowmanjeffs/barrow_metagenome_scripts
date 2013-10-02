# -*- coding: utf-8 -*-
"""
Created on Sun May 26 09:38:57 2013

@author: Jeff

This script requires gzipped sam files that have been reduced to only those
reads that map to the reference genomes (for example obtained by fgrep 'gi').  
It produces a white space delimited file of taxonomy, number of mapped reads,
and length of the reference genome.

V3 adds an additional dictionary (saved as pickle), with keys as reads and
bacterial species (or taxonomies) as values.  This reduces bias against 
well represented taxa in among the genomes, for which it is difficult to obtain
a unique read.

"""
#### parse sam files for genome length and number of reads ####

import gzip
import os

gis = set()

## from each sam file in genome_alignment return the length of the reference and
## the number of reads that mapped, saving these values to a dictionary with
## gi number as key to allow mapping to taxonomy

genome_dict = {} # key = gi, value = ncbi, genome length, number of mapped reads
gi_id = {} # key = ncbi, value = gi
for filename in os.listdir('.'):
    if filename.endswith('sam.gz'):
        name = filename.split('.')
        name = name[1]
        with gzip.open(filename, 'rb') as sam:
            l = 0
            for line in sam:
                if line.startswith('@'):
                    line = line.split('\t')
                    length = line[2]
                    gi = line[1].split('|')[1]
                    length = length.split(':')
                    length = length[1].rstrip()
                else:
                    l = l + 1
            print gi, l
            gis.add(gi)
            genome_dict[gi] = name, length, l
            gi_id[name] = gi
            
## pickle gi_id for use in downstream script 
import cPickle
cPickle.dump(gi_id, open('bacterial_genomes_ncbiid_gi.p', 'wb'))

## get taxids for each genome, searching by gi.  to save time on repeat runs
## the mapping is saved to a file, which is parsed as a list.  only gi numbers
## not already mapped to a taxid are searched

print 'looking for taxids'

gi_complete = set()
try:
    with open('genome_alignment_gi_taxid.txt', 'r') as id_complete:
        for complete in id_complete:
            gi_complete.add(complete.split('\t')[0])
    bad_gis = cPickle.load(open('genome_gis_no_taxid.p', 'rb'))
    have_bad_gi = True
    
except IOError:
    print 'no genome_alignment_gi_taxid.txt file, making new'
    gi_complete = set()
    bad_gis = set()
    have_bad_gi = False

gi_get = set()       
for gi in gis:
    if gi not in gi_complete:
        if gi not in bad_gis:
            gi_get.add(gi)
        
if len(gi_get) > 0:
    l = 0
    with gzip.open('/volumes/deming/databases/gi_taxid_nucl.dmp.gz', 'rb') as gi_taxid_dmp, open('genome_alignment_gi_taxid.txt', 'a') as gi_taxid_output:
        for line in gi_taxid_dmp:
            l = l + 1
            print l
            line_split = line.split('\t')
            if str(line_split[0]) in gis:
                print >> gi_taxid_output, line,

## look for gis with no taxid and generate a pickle
                
if have_bad_gi == False:
    bad_gi = set()
    with open('genome_alignment_gi_taxid.txt', 'r') as id_complete:
        for complete in id_complete:
            gi_complete.add(complete.split('\t')[0])
    for gi in gis:
        if gi not in gi_complete:
            bad_gi.add(gi)
    cPickle.dump(bad_gi, open('genome_gis_no_taxid.p', 'wb'))
    
print 'obtaining taxonomy'

#### obtain taxonomy ####

## import taxonomy module
from cogent.parse.ncbi_taxonomy import NcbiTaxonomyFromFiles

## open taxonomy files and set desired ranks, takes a long time to load, so see
## if it's already loaded first
try:
    root = tree.Root
except NameError:
    tree = NcbiTaxonomyFromFiles(open('/volumes/deming/databases/nodes.dmp'), open('/volumes/deming/databases/names.dmp'))
    root = tree.Root
    
ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

## define function to get taxonomy from taxid
def get_lineage(node, my_ranks):
    ranks_lookup = dict([(r,idx) for idx,r in enumerate(my_ranks)])
    lineage = [None] * len(my_ranks)
    curr = node
    while curr.Parent is not None:
        if curr.Rank in ranks_lookup:
            lineage[ranks_lookup[curr.Rank]] = curr.Name
        curr = curr.Parent
    return lineage

## create new dictionary to map gi to taxonomy

print 'mapping gi to taxonomy'

tax_dict = {} # key = gi, # value = taxonomy

## map gi to taxonomy  
with open('genome_alignment_gi_taxid.txt', 'r') as gi_taxid:
    for each in gi_taxid:
        each = each.split('\t')
        taxid = int(each[1].rstrip())
        node = tree.ById[taxid]
        tax = get_lineage(node, ranks)
        tax_dict[each[0]] = tax

## pickle the dictionary for use by downstream scripts
cPickle.dump(tax_dict, open('bacterial_genomes_taxonomy.p', 'wb'))

## get read data for each taxonomy, searching both dictionaries by gi number
output = open('mapped_reads_length.txt', 'w')        
for key in tax_dict.keys():
    if key in genome_dict.keys():
        print >> output, tax_dict[key][3], '\t', tax_dict[key][6], '\t', genome_dict[key][0], '\t', genome_dict[key][1], '\t', genome_dict[key][2]

output.close()

## create dictionary with keys as reads and mapped species as values

print 'mapping reads to species'

species_dict = {} # key = species, value = ncbi ids
for gi in tax_dict.keys():
    if gi in genome_dict.keys():
        try:
            name = genome_dict[gi][0]
            sp = tax_dict[gi][6]
            temp = species_dict[sp]
            temp.append(name)
            species_dict[sp] = temp
        except KeyError:
            name = genome_dict[gi][0]
            sp = tax_dict[gi][6]
            species_dict[sp] = [name]
        
read_dict = cPickle.load(open('44_reads_genomes_mapping.p', 'rb'))

species_unique = {} # key = species, value = number of unique reads

i = 0
for read in read_dict.keys():
    i = i + 1
    n = 0
    s = ''
    for name in read_dict[read]:
        for species in species_dict.keys():
            if name in species_dict[species]:
                n = n + 1
                s = species
    if n == 1:
        try:
            temp = species_unique[s]
            temp = temp + 1
            species_unique[s] = temp
        except KeyError:
            species_unique[s] = 1
        print i,s, species_unique[s]
        
with open('num_reads_mapped_2_unique_species.txt', 'w') as output_2:
    for s in species_unique.keys():
        print >> output_2, s,'\t',str(species_unique[s])
        
    