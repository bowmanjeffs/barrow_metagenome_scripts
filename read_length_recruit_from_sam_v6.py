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

V6 combines get_coverage_gaps_v2.py with this script.

"""
read_size = 96 ## mean QC'd read length
bins = 1000 ## size of bin for determing coverage breadth
cutoff = 0.5 ## coverage below this value will not be reported in mapped_genomes_regions.txt.gz
min_gap = 10000

from Bio import SeqIO
import subprocess

## pasting my standard function for finding and translating orfs.  overkill for this...

def find_pros_with_trans(name, seq, trans_table, min_protein_length):
    error_file=open('find_pros_with_trans_errors.txt','a')
    errors = []
    answer = []
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            untrans = str(nuc) #not frame, you only want to count DNA from the true start
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start) #starting from aa_start find first instance of * and return index
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end-aa_start >= min_protein_length:
                    start = frame+aa_start*3 # counting forward from left side of reverse complement!
                    end = frame+aa_end*3 #+3   
                    #print strand, frame, aa_start, aa_end, 'dna start = ', start, 'dna end = ', end
                    pro=trans[aa_start:aa_end]
                    #print 'pro=',pro[:10],'...',pro[-10:]
                    dna=untrans[start:end]
                    dna_trans=str(nuc[start:end].translate(trans_table))
                    #print 'dna_trans=',dna_trans[:10],'...',dna_trans[-10:]
                    if pro!=dna_trans:
                        errors.append((start,end,frame))
                    answer.append((strand, frame, start, end, pro, dna))
                aa_start = aa_end+1
    answer.sort()
    errors.sort()
    if len(errors) != 0:
        for error in errors:
            print 'error!', name, start, end, frame, strand
            print >> error_file, name, start, end, frame, strand
    else:
        #print 'you have no errors!'
        print >> error_file, name, 'has no errors!'
    return answer
    error_file.close()

#### parse sam files for genome length and number of reads ####

import gzip
import os
import re

gis = set()

"""
from each sam file in genome_alignment obtain:
    length of the reference
    number of mapped reads
    coverage
    breadth
    plasmid or chromosome
Save these values to a dictionary with gi number as key to allow mapping to
taxonomy.  For each reference above the specicied coverage, print the start
postition for each mapped read and create a fasta file for each gap region
larger than the specified size.  Set up a blast analysis on Hyak for these
gap regions.
"""

genome_dict = {} # key = gi, value = ncbi, genome length, number of mapped reads, coverage
gi_id = {} # key = ncbi, value = gi

regions_output = gzip.open('mapped_genomes_regions.txt.gz', 'wb')
master_shell = open('gap_blast_submit.sh', 'w')

for filename in os.listdir('genome_alignment'):
    if filename.endswith('sam.gz'):
        name = filename.split('.')
        full_name = name[1]+'.'+name[2]
        name = name[1]
        with gzip.open('/volumes/deming/databases/bwa/'+full_name+'.fasta.gz', 'rb') as ref_fasta:
            for line in ref_fasta:
                if line.startswith('>'):
                    if 'plasmid' in line:
                        typ = 'plasmid'
                    elif 'Plasmid' in line:
                        typ = 'plasmid'
                    else:
                        typ = 'chromosome'
        
        with gzip.open('genome_alignment/'+filename, 'rb') as sam:
            cov_reg = {} ## coverage regions within reference as start, end
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
                    line = line.split('\t')
                    start = int(line[3])
                    end = int(line[3]) + len(line[9])
                    cov_reg[l] = start, end
                                
            gaps = {}            
            cov_set = set() ## unique positions in genome with a read mapped
            cov_dict = {} ## key = reference genome position, value = number of bases mapped to site
            
            for key in cov_reg.keys():
                reg = cov_reg[key]
                for r in range(reg[0], reg[1] + 1):
                    cov_set.add(r)
                    try:
                        t = cov_dict[r]
                        cov_dict[r] = t + 1
                    except KeyError:
                        cov_dict[r] = 1
            
            g = 0            
            cov_set = sorted(cov_set)
            for i,e in enumerate(cov_set):
                try:
                    if cov_set[i + 1] - e > min_gap: ## smallest for fasta creation
                        g = g + 1
                        gaps[g] = e, cov_set[i + 1]
                except IndexError:
                    continue                    
                    
            cov = (read_size * l) / float(length) ## coverage by traditional calculation
            br = len(cov_set) / float(length) ## fraction of positions with a read mapped
            
            print gi, length, len(cov_set), l, cov, br, typ
            gis.add(gi)
            genome_dict[gi] = name, length, l, cov, br, typ, len(cov_set)
            gi_id[name] = gi
            
            ## if coverage is above the cutoff, go ahead and print all the mapped locations to file
            
            if cov >= cutoff:
                for each in cov_dict.keys():
                    print >> regions_output, name+'\t'+str(each)+'\t'+str(cov_dict[each])+'\t'+str(length)
                
                with open(name+'_coverage_gaps.fasta', 'w') as fasta_out:
                    fasta_in = SeqIO.read(gzip.open('/Volumes/deming/databases/bwa/'+full_name+'.fasta.gz', 'rb'), 'fasta')
                    for gap in gaps.keys():
                        start = gaps[gap][0]
                        end = gaps[gap][1]
                        fasta = str(fasta_in.seq)[start:end]
                        print >> fasta_out, '>'+name+'_'+str(start)+'_'+str(end)+'\n'+fasta
                        print 'got gap!', name, start, end
                
                with open(name+'_gap_orfs.fasta','w') as orf_fasta:      
                    for gap in SeqIO.parse(name+'_coverage_gaps.fasta','fasta'):
                        pro_list = find_pros_with_trans(gap.id, gap.seq, 11, 100)
                        for strand, frame, start, end, pro, dna in pro_list:
                            print >> orf_fasta, '>'+gap.id+'_'+str(start)+'_'+str(end)+'_strand='+str(strand)
                            print >> orf_fasta, pro

                blast_out = open(name+'_blast.sh', 'w') 
                
                print >> blast_out, '#!/bin/bash'
                print >> blast_out, '#PBS -N "'+name+'"'
                print >> blast_out, '#PBS -d /gscratch/coenv/bowmanjs/barrow_metagenome_working/gap_blast'
                print >> blast_out, '#PBS -l nodes=1:ppn=12,feature=12core,walltime=4:00:00'
                print >> blast_out, '#PBS -M bowmanjs@u.washington.edu'
                print >> blast_out, '#PBS -m abe'
                print >> blast_out, 'module load epd_7.3_2'
                
                print >> blast_out, 'blastp -num_threads 12 '\
                '-max_target_seqs 1 ' \
                '-evalue 1e-30 ' \
                '-db /gscratch/coenv/bowmanjs/working_database/refseq_protein ' \
                '-outfmt 5 ' \
                '-out '+name+'_gap_orfs_blast.xml ' \
                '-query '+name+'_gap_orfs.fasta;' \
                'gzip -f '+name+'_gap_orfs_blast.xml;' \
                'python parse_blast_xml_gz_v4.py '+name+'_gap_orfs_blast.xml.gz'
            
                blast_out.close()
                print >> master_shell, 'qsub -q bf '+name+'_blast.sh'            

regions_output.close()
master_shell.close()

subprocess.call('rm -r to_hyak;' \
'mkdir to_hyak;' \
'cp *_gap_orfs.fasta to_hyak/;' \
'cp *blast.sh to_hyak/;' \
'cp gap_blast_submit.sh to_hyak/;' \
'cp parse_blast_xml_gz_v4.py to_hyak/;' \
'tar -czvf to_hyak.tgz to_hyak', shell = True)
            
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

## get read data for each taxonomy (strain), searching both dictionaries by gi number
with open('mapped_reads_length.txt', 'w') as output:
    print >> output, 'order'+'\t'+'strain'+'\t'+'id'+'\t'+'length'+'\t'+'mapped'+'\t'+'cov'+'\t'+'breadth'+'\t'+'type'+'\t'+'pos_mapped'       
    for key in tax_dict.keys():
        if key in genome_dict.keys():
            clean_tax = re.sub('\'', '', str(tax_dict[key][6]))
            print >> output, str(tax_dict[key][3])+'\t'+clean_tax+'\t'+str(genome_dict[key][0])+'\t'+str(genome_dict[key][1])+'\t'+str(genome_dict[key][2])+'\t'+str(genome_dict[key][3])+'\t'+str(genome_dict[key][4])+'\t'+str(genome_dict[key][5])+'\t'+str(genome_dict[key][6])

## create dictionary with reads as keys and mapped species as values

print 'mapping reads to species'

species_dict = {} # key = species, value = ncbi ids
length_dict = {} # key = species, value = lengths of all strains in species

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
        try: 
            length = genome_dict[gi][1]
            sp = tax_dict[gi][6]
            temp = length_dict[sp]
            temp.append(length)
            length_dict[sp] = temp 
        except KeyError:
            length = genome_dict[gi][1]
            sp = tax_dict[gi][6]
            length_dict[sp] = [length]
        
read_dict = cPickle.load(open('genome_alignment/44_reads_genomes_mapping.p', 'rb'))

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
        
with open('44_num_reads_mapped_2_unique_species.txt', 'w') as output_2, open('44_cov_reads_mapped_2_unique_species.txt', 'w') as output_3:
    for s in species_unique.keys():
        cov = (int(species_unique[s]) * 96) / (sum(map(float, length_dict[s])))
        print >> output_2, s+'\t'+str(species_unique[s])
        print >> output_3, s+'\t'+str(cov)