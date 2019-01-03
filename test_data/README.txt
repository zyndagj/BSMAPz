##############################################
Reference
##############################################

Generated using:
   http://www.faculty.ucr.edu/~mmaduro/random.htm

chr1 - 3000bp and 0.40 GC
chr2 - 3000bp and 0.60 GC

and formatted using fasta_formatter from fastx_toolkit

##############################################
Reads
##############################################

https://github.com/FelixKrueger/Sherman

Using Sherman

Single reads:
   perl Sherman -l 100 -n 2000 --genome_folder . -q 40 -I 100 -X 500 -cr 95 -e 3

#Genome folder was specified as ./
#
#Selected general parameters:
#----------------------------------------------------------------------------------------------------
#Single-end read output format selected
#sequence length:	100 bp
#number of sequences being generated:	2000
#
#Possible sources of contamination:
#----------------------------------------------------------------------------------------------------
#overall error rate:	3%
#bisulfite conversion rate:	95%
#
#
#Now reading in and storing sequence information of the genome specified in: ./
#
#chr chr1 (3000 bp)
#chr chr2 (3000 bp)
#
#Generating quality values with a user defined decaying per-bp error rate of 3%
#Starting to work out the slope of the error curve
#Error rates per bp will be modelled according to the formula:
#default base quality - 0.034286*position[bp] + 0.0009263*(position[bp]**2)) - 0.00001*(position[bp]**3)*4.2857)
#
#
#Final report:
#----------------------------------------------------------------------------------------------------
#2000 genomic sequences were successfully generated in total (+ strand: 1002	 - strand: 998)
#Cytosines bisulfite converted in any context: 95.05%
#Random sequencing errors introduced in total: 5294 (of 200000 bp in total) (percentage: 2.65)

Paired reads:
   perl Sherman -l 100 -n 2000 --genome_folder . -q 40 -I 100 -X 500 -cr 95 -e 3 -pe

#Genome folder was specified as ./
#
#Selected general parameters:
#----------------------------------------------------------------------------------------------------
#Paired-end reads selected. Fragment length will be 100-500 bp
#sequence length:	100 bp
#number of sequences being generated:	2000
#
#Possible sources of contamination:
#----------------------------------------------------------------------------------------------------
#overall error rate:	3%
#bisulfite conversion rate:	95%
#
#
#Now reading in and storing sequence information of the genome specified in: ./
#
#chr chr1 (3000 bp)
#chr chr2 (3000 bp)
#
#Generating quality values with a user defined decaying per-bp error rate of 3%
#Starting to work out the slope of the error curve
#Error rates per bp will be modelled according to the formula:
#default base quality - 0.034286*position[bp] + 0.0009263*(position[bp]**2)) - 0.00001*(position[bp]**3)*4.2857)
#
#The length of seq1 or seq2 were not equal to the sequence length of 101
#
#Final report:
#----------------------------------------------------------------------------------------------------
#2000 genomic sequences were successfully generated in total (+ strand: 1012	 - strand: 988)
#Cytosines bisulfite converted in any context: 94.95%
#Random sequencing errors introduced in total: 10583 (of 400000 bp in total) (percentage: 2.65)

##############################################
Alignments
##############################################

./bsmap -a simulated.fastq -z 33 -p 2 -q 20 -d test.fasta -S 77345 -w 1000 -o r1.bam
./bsmap -a simulated_1.fastq -b simulated_2.fastq -z 33 -p 2 -q 20 -d test.fasta -S 77345 -w 1000 -o paired.bam

##############################################
Methylation calls
##############################################

python methratio.py -z -r -d test.fasta -o original.mr r1.bam 
python methratio.py -z -r -d test.fasta -o original_paired.mr paired.bam 
