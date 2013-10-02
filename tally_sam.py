# -*- coding: utf-8 -*-
"""
Created on Sun May 19 11:58:37 2013

@author: Jeff
"""
import gzip
import sys

tally = {}
with gzip.open(sys.argv[1], 'rb') as sam:
    for line in sam:
        if line.startswith('@') == False:
            line = line.split('\t')
            if line[2] not in tally.keys():
                tally[line[2]] = 1
            else:
                tally[line[2]] = tally[line[2]] + 1

output = open(sys.argv[1]+'.tally.txt', 'w')
for key in tally.keys():
    print key, tally[key]
    print >> output, key, tally[key]
output.close()
    
