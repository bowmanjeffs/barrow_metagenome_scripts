barrow_metagenome_scripts
=========================
These scripts and workflow describe our analysis from Bowman et. al, 2014, FEMS Microbial Ecology 89:2 376-387.  Please cite this publication if you use any portion of this code in your own work.

The python and R scripts included here are used as indicated in the below workflow.  The workflow does not describe a complete pipeline, some modification will be required to work on your system.


For a second round assembly and annotation I'm prepping my sequence files with SEAtAR.  SEAStAR has some advantages over solexaQA + khmer, including the ability to remove technical replicates from the dataset.  This is pretty slow for large datasets, and works best if you have a large insert size.  For my dataset it removes ~ %70 of reads, which is probably excessive, but does make downstream analysis a lot faster and presumably makes quantitative representations of different proteins in the dataset more defendable.  Much of this analysis follows the guidance of the lambda phage vignette available in the SEAStAR package (Google "Armbrust lab seastar").

	#!/bin/bash

	## set memory limit
	ulimit -v 600000000


	## remove technical replicates (reducing PCR bias) and low-complexity sequence
	fastq_nodup -z -l 15 -d 1 -e 2 -v JB44_AGTTCC_L002 44_dedup

	## important new step here... the enzymatic primer removal doesn't get all primers from WGA derived DNA.  Cleave primers from those reads with.  Will not get all, but will get most.
	python2.7 primer_purge.py 44_dedup.read1.fastq.gz
	python2.7 primer_purge.py 44_dedup.read2.fastq.gz

	## trim for quality
	trimfastq -z --mates_file -p 0.5 -l 34 -m 34 --add_len 44_dedup 44_trim

	## assembly round 1
	velveth seastar_44/ 23 -fastq.gz -shortPaired 44_trim.mates.fastq.gz -short 44_trim.single.fastq.gz > seastar_44.velveth.log 2>&1

	velvetg seastar_44/ -scaffolding no -read_trkg no -ins_length auto -ins_length_sd auto -exp_cov 50 -cov_cutoff 5 -min_contig_lgth 100 > seastar_44.velvetg.log 2>&1

	## round 1 graph2 will be huge, remove
	rm seastar_44/graph2

	## assembly round 2

	bwa index -a is seastar_44/contigs.fa

	bwa aln -n 0.001 -l 18 -t 8 seastar_44/contigs.fa 44_trim.read1.fastq.gz > 44_trim.read1.sai

	bwa samse -n 1000000 seastar_44/contigs.fa 44_trim.read1.sai 44_trim.read1.fastq.gz 2>44_trim.read1.samse.log > 44_trim.read1.sam

	bwa aln -n 0.001 -l 18 -t 8 seastar_44/contigs.fa 44_trim.read2.fastq.gz > 44_trim.read2.sai

	bwa samse -n 1000000 seastar_44/contigs.fa 44_trim.read2.sai 44_trim.read2.fastq.gz 2> 44_trim.read2.samse.log > 44_trim.read2.sam

	bwa aln -n 0.001 -l 18 -t 8 seastar_44/contigs.fa 44_trim.single.fastq.gz > 44_trim.single.sai

	bwa samse -n 1000000 seastar_44/contigs.fa 44_trim.single.sai 44_trim.single.fastq.gz 2> 44_trim.single.samse.log > 44_trim.single.sam

	python sam_to_fasta.py 44_trim.read1.sam

	python sam_to_fasta.py 44_trim.read2.sam

	python sam_to_fasta.py 44_trim.single.sam

	velveth seastar_44.2/ 23 -fasta -long seastar_44/contigs.fa -short 44_trim.single.sam.aligned.fasta 44_trim.read1.sam.aligned.fasta 44_trim.read2.sam.aligned.fasta > seastar_44.r2.velveth.log 2>&1

	velvetg seastar_44.2/ -scaffolding no -read_trkg no -ins_length auto -ins_length_sd auto -exp_cov 50 -cov_cutoff 5 -min_contig_lgth 100 > seastar_44.r2.velvetg.log 2>&1

	## repeat assembly round 2 until N50 stabilizes.  If you have a large insert, or succeed in building large contigs, after assembly make use of SEAStAR's graph manipulation tools for tetranucleotide clustering and adding in mate-pairing information.

	#### add annotation ####

follow-on analysis
--------

I'm interested in evaluating the mapped distribution of reads across all available prokaryotic genomes, as a means of identifying both highly conserved and unique reads.  The latter should help in verifying the taxonomic makeup of the community.  I start by downloading the 4,723 genetic elements available as a portions of whole genomes from Genbank, and mapping all reads to them:

	#!/bin/bash

	ulimit -v 20000000
	
	rm master_tally.txt
	
	for f in /volumes/deming/databases/bwa/*.fasta.gz; do
		n=$(basename "$f" .fasta.gz)
		bwa index -a is $f
		bwa aln -n 0.001 -t 8 -l 12 -k 2 $f ../44_trim_mates_singles.fasta.gz > 44.$n.sai
		bwa samse $f 44.$n.sai ../44_trim_mates_singles.fasta.gz > 44.$n.sam
		fgrep 'gi|' 44.$n.sam > 44.$n.small.sam
		gzip 44.$n.small.sam  
		rm 44.$n.sam
		python tally_sam.py 44.$n.small.sam.gz
		grep '>' 44.$n.small.sam.gz.tally.txt >> master_tally.txt
	done

Next I utilize a script tally_reads_in_sam.py which produces a pickled dictionary where the keys are reads and the values are lists of genomes that each key maps to.

Next I implement read_length_recruit_from_sam.py.  This produces a white-space delimited file of taxonomy, number of mapped reads, and length for each reference genome.  This allows me to plot the number of mapped reads against element length (with knowledge of taxonomy).  This script also creates two more pickled dictionaries, one mapping ncbi id to gi number for each reference genome, and another mapping gi number to taxonomy for each reference genome.

The final script in the series (reads_genomes_distribution.py) attempts to do something analytical with this information.  It reads the pickles and produces a file, with taxonomy, of all reads that have mapped to only one genome.  It produces a second file of reads that have mapped to many (top 99.9 %) genomes.  A fasta file of these is two read sets is created. 

The next step is to blast all these reads, and to interpret the pattern of mapping among the single-mapping set to see if a statement can be made about the composition of the community.  For example, if these reads have mapped (uniquely) to a diverse array of bacteria it suggests that the community is truly diverse.  If the opposite were true we would expect to see many reads map to a single genome and no other (if a bacterium was dominant in the environment that was only distantly related to the available genomes). 

Gap regions (regions with no read mapping) are of interest in this analysis.  To identify gap regions from my sam files (from all available genomes) I use get_coverage_gaps.py and get_coverage_gaps.r.  Together these produce histogram recruitment plots, fasta files of the gapped regions, and then blast the gapped regions (blastx against refseq_protein).

