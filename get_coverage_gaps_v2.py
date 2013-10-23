# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 09:40:25 2013

@author: Jeff
"""
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
                    print strand, frame, aa_start, aa_end, 'dna start = ', start, 'dna end = ', end
                    pro=trans[aa_start:aa_end]
                    print 'pro=',pro[:10],'...',pro[-10:]
                    dna=untrans[start:end]
                    dna_trans=str(nuc[start:end].translate(trans_table))
                    print 'dna_trans=',dna_trans[:10],'...',dna_trans[-10:]
                    if pro!=dna_trans:
                        errors.append((start,end,frame))
                    answer.append((strand, frame, start, end, pro, dna))
                aa_start = aa_end+1
    answer.sort()
    errors.sort()
    if len(errors)!=0:
        print 'you have errors in:'
        for error in errors:
            print name,start,end,frame,strand
            print >> error_file, name, start, end, frame, strand
    else:
        print 'you have no errors!'
        print >> error_file, name, 'has no errors!'
    return answer
    error_file.close()

import os
import gzip
import re

with gzip.open('mapped_genomes_regions.txt.gz', 'wb') as output, open('gap_blast_submit.sh', 'w') as master_shell, open('genome_cov_breadth.txt', 'w') as cov_out:
    print >> master_shell, '#!/bin/bash'    
    for f in os.listdir('.'):
        if f.endswith('fasta.gz.small.sam.gz'):
            temp = [] # recycled list for start of aligned region
            l = 0
            name = re.sub('.small.sam.gz', '', f)
            name = name.lstrip('44.')
            short_name = re.sub('.fasta.gz', '', name)
            with gzip.open(f, 'rb') as sam:
                starts = []
                ends = []
                for line in sam:
                    if line.startswith('@') == False:
                        l = l + 1
                        line = line.split('\t')
                        temp.append(line[3])
                        start = int(line[3])
                        end = int(line[3]) + len(line[9])
                        starts.append(start)
                        ends.append(end)
                    else:
                        length = re.search('LN:[0-9]*', line)
                        length = length.group()
                        length = re.sub('LN:', '', length)
                        length = int(length)
            starts = sorted(starts)
            ends = sorted(ends)                
            gaps = {}
            c = 0
            for i,e in enumerate(ends):
                try:
                    if starts[i + 1] - e > 10000: # lowest size limit for annotating gaps
                        c = c + 1
                        #print e, starts[i + 1]
                        gaps[c] = e, starts[i + 1]
                except IndexError:
                    continue
                
            ## calculate coverage and breadth and set limit for continuing analysis
            
            cov = (96.0 * l) / length
            br = len(gaps.keys()) / float((length / 1000)) # fraction of 1000 bp bins with at least on read mapped
            
            print >> cov_out, short_name, length, l, cov, br
            
            if cov >= 0.5:
                for each in temp:
                    print >> output, short_name+'\t'+str(each)
                
                with open(short_name+'_coverage_gaps.fasta', 'w') as fasta_out:
                    fasta_in = SeqIO.read(gzip.open('/Volumes/deming/databases/bwa/'+name, 'rb'), 'fasta')
                    print 'finding', len(gaps.keys()), 'gaps in', short_name
                    for gap in gaps.keys():
                        start = gaps[gap][0]
                        end = gaps[gap][1]
                        fasta = str(fasta_in.seq)[start:end]
                        print >> fasta_out, '>'+short_name+'_'+str(start)+'_'+str(end)+'\n'+fasta
                        print short_name, start, end
                
                with open(short_name+'_gap_orfs.fasta','w') as orf_fasta:      
                    for gap in SeqIO.parse(short_name+'_coverage_gaps.fasta','fasta'):
                        pro_list = find_pros_with_trans(gap.id, gap.seq, 11, 100)
                        for strand, frame, start, end, pro, dna in pro_list:
                            print >> orf_fasta, '>'+gap.id+'_'+str(start)+'_'+str(end)+'_strand='+str(strand)
                            print >> orf_fasta, pro
        
                blast_out = open(short_name+'_blast.sh', 'w') 
                
                print >> blast_out, '#!/bin/bash'
                print >> blast_out, '#PBS -N "'+short_name+'"'
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
                '-out '+short_name+'_gap_orfs_blast.xml ' \
                '-query '+short_name+'_gap_orfs.fasta;' \
                'gzip -f '+short_name+'_gap_orfs_blast.xml;' \
                'python parse_blast_xml_gz_v4.py '+short_name+'_gap_orfs_blast.xml.gz'
            
                blast_out.close()
                print >> master_shell, 'qsub -q bf '+short_name+'_blast.sh'
       
master_shell.close()
subprocess.call('rm -r to_hyak;' \
'mkdir to_hyak;' \
'cp *_gap_orfs.fasta to_hyak/;' \
'cp *blast.sh to_hyak/;' \
'cp gap_blast_submit.sh to_hyak/;' \
'cp parse_blast_xml_gz_v4.py to_hyak/;' \
'tar -czvf to_hyak.tgz to_hyak', shell = True)
                        