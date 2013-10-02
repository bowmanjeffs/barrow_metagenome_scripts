# -*- coding: utf-8 -*-
"""
Created on Sun May 26 16:38:15 2013

@author: Jeff

This script searches a directory for all sam files (gzipped).  The sam files
should be from one dataset against different references (one ref per sam) and
must be reduced to only mapped reads.  The script tallys the number of times 
each read aligns, and tracks the reference genomes the reads align to.  The 
script produces a dictionary that is saved as a pickle.  The dictionary keys 
are read ids, values are reference names.
"""
import gzip
import os
import cPickle

## initialize dictionary that will hold reads as keys and genome names as 
## values
reads_dict = {}
reads_set = set()

## parse sam.gz files, creating a new dictionary key for new reads or adding
## the genome name as a value to existing keys.
n = 0
for filename in os.listdir('.'):
    if filename[-6:] == 'sam.gz':
        n = n + 1
        name = filename.split('.')
        name = name[1]
        print n, name, len(reads_set)
        with gzip.open(filename, 'rb') as sam:
            for line in sam:
                if line.startswith('@') == False:
                    line = line.split('\t')
                    read = line[0]
                    if read not in reads_set:
                        reads_dict[read] = [name]
                        reads_set.add(read)
                    else:
                        reads_dict[read].append(name)

## print the key and number of values for each read in dictionary
output = open('44_reads_mapped_genomes.txt', 'w')
for read in reads_dict.keys():
    print read, len(reads_dict[read])
    print >> output, read, '\t', len(reads_dict[read])   
output.close()

## this script likely took a long time to run.  save the dictionary as a pickle
## for downstream analysis so that you don't have to run again!
cPickle.dump(reads_dict, open('44_reads_genomes_mapping.p', 'wb'))