# -*- coding: utf-8 -*-
"""
Created on Sat Aug 03 20:30:05 2013

@author: Jeff
"""

import gzip
import sys
import subprocess

#44_dedup.read1.fastq.gz

n = 0
l = 0
with gzip.open(sys.argv[1], 'rb') as fasta, gzip.open('temp.fasta.gz', 'wb') as output:
    for line in fasta:
        line = line.rstrip('\n')
        n = n + 1
        if n == 2:
            strip = False       
            if line.startswith('GGTGTGTTGGGTGTGTTTGGATG'):
                l = l + 1
                print 
                line = line[23:]
                strip = True
        if n == 4:
            if strip == True:
                line = line[23:]
            n = 0
        print >> output, line
    

subprocess.call('mv '+sys.argv[1]+' '+sys.argv[1]+'.original', shell = True)
subprocess.call('mv temp.fasta.gz '+sys.argv[1], shell = True)