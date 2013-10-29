# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 12:55:25 2013

@author: Jeff
"""
min_gap = 5000

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

#genome_dict = {} # key = gi, value = ncbi, genome length, number of mapped reads, coverage
#gi_id = {} # key = ncbi, value = gi

subprocess.call('rm -r to_hyak;mkdir to_hyak', shell = True)

master_shell = open('to_hyak/gap_blast_submit.sh', 'w')

target = []
with open('ncbi_id_top_five_rhizobiales.txt', 'r') as ids: ## coming from get_coverage_gaps.r
    for line in ids:
        line = line.rstrip('\n')
        target.append(line)

for filename in os.listdir('genome_alignment'):
    if filename.endswith('sam.gz'):
        for tar in target:
            if tar in filename:
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
                    size = []
                    
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
                            size.append(int(end) - int(start))
                                        
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
                                    
                    with open('to_hyak/'+name+'_coverage_gaps.fasta', 'w') as fasta_out:
                        fasta_in = SeqIO.read(gzip.open('/Volumes/deming/databases/bwa/'+full_name+'.fasta.gz', 'rb'), 'fasta')
                        for gap in gaps.keys():
                            start = gaps[gap][0]
                            end = gaps[gap][1]
                            fasta = str(fasta_in.seq)[start:end]
                            print >> fasta_out, '>'+name+'_'+str(start)+'_'+str(end)+'\n'+fasta
                            print 'got gap!', name, start, end
                    
                    with open('to_hyak/'+name+'_gap_orfs.fasta','w') as orf_fasta:      
                        for gap in SeqIO.parse(name+'_coverage_gaps.fasta','fasta'):
                            pro_list = find_pros_with_trans(gap.id, gap.seq, 11, 100)
                            for strand, frame, start, end, pro, dna in pro_list:
                                print >> orf_fasta, '>'+gap.id+'_'+str(start)+'_'+str(end)+'_strand='+str(strand)
                                print >> orf_fasta, pro
        
                        blast_out = open('to_hyak/'+name+'_blast.sh', 'w') 
                        
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

master_shell.close()

subprocess.call('cp parse_blast_xml_gz_v4.py to_hyak/;' \
'tar -czvf to_hyak.tgz to_hyak', shell = True)