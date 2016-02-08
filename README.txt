BSMAP 2.73

1. Introduction

BSMAP is a short reads mapping program for bisulfite sequencing in DNA methylation study.  Bisulfite treatment coupled with next generation sequencing could estimate the methylation ratio of every single Cytosine location in the genome by mapping high throughput bisulfite reads to the reference sequences.

Bisulfite mapping is different from usual sequence mapping in two aspects: 1) The additional C/T mapping is asymmetric, a T in the read could be aligned to C in the reference but not vice versa. 2) The Watson and Crick strand are not complimentary after bisulfite treatment.  Each read need to be compared with 4 reference sequences, namely BSW(bisulfite Watson), BSWC(reverse complimentary of BSW), BSC(bisulfite Crick) and BSCC(reverse complimentary of BS).

BSMAP is designed to be a general-purpose mapping program to handle these special characteristics of bisulfite mapping.  It is based on the open source program SOAP (Short Oligo Alignment Program).  

Main features: 
    read length up to 144nt, allow up to 15 mismatches, 1 continous gap up to 3nt
    support pair end mapping, support parallel mapping
    support SAM format input/output, support gzipped FASTA/FASTQ format input  
    support both whole genome (WGBS) and reduced representation bisulfite sequencing (RRBS)
    support trimming adapters and low quality sequences from 3'end of reads
    allow different running modes with flexible speed/memory/sensitivity to run on different hardware configurations
    include script to extract methylation ratios

BSMAP is under GNU Public License (GPL).


2. Installation

BSMAP is designed for linux64 platform. 
First unpackage the source code:
    $ tar zxfv bamsp-X.Y.tgz

Make executable binary:
    $ make

Install the binary into system default path: (optional)
    $ make install
  

3. Usage

bsmap <option>

option:
-a  <str>   query file, FASTA/FASTQ/BAM format. support gzipped FASTA/FASTQ format. The input format will be auto-detected. (required)
-b  <str>   query file b for pair end data, FASTA/FASTQ/BAM format.  support gzipped FASTA/FASTQ format. The input format will be auto-detected. 
            if the input format is in BAM format, it should be the same as the file specified by "-a" option.
            BSMAP will read the two sets of reads w.r.t to the 0x40/0x80 flag in the input BAM file. 
            (required for pair-end mapping)
-d  <str>   reference sequences file, FASTA format. support gzipped FASTA format. (required)
-o  <str>   output alignment file, if filename has .sam suffix, the output will be in SAM format, 
            if the filename has .bam suffix, the output file be in sorted and indexed BAM file and a filename.bai index file will be generated, 
            for other filename suffix the output is in BSP format. 
            if this option is omitted, the output will be written to STDOUT in SAM format to allow pipe into downstream programs.
-s  <int>   seed size, default=16, min=8, max=16. 
            Selecting short seed size may drastically slow down the alignment, usually the default setting is good enough for mammalian genome.
-v  <float> if this value is between 0 and 1, it's interpreted as the mismatch rate w.r.t to the read length.
            otherwise it's interpreted as the maximum number of mismatches allowed on a read, <=15.
            example: -v 5 (max #mismatches = 5), -v 0.1 (max #mismatches = read_length * 10%)
            default=0.08.
-3          using 3-nucleotide mapping approach. (default: off)
-g  <int>   gap size, BSMAP only allow 1 continuous gap (insertion or deletion) with up to 3 nucleotides. default=0
            gaps will not be allowed within 6nt of the read edges.
            the number of mismatches of gapped algnment is calculated as #gap_size+#mismatches+1
-k  <float> set the cut-off ratio for over-represented kmers, default=1e-06
            example: -k 1e-6 means the top 0.0001% over-represented kmer will be skipped in alignment
-w  <int>   max number of equal best hits to count, smaller will be faster, default=MAXHITS in makefile
-r  [0,1,2] how to report repeat hits, 0=none(unique hit/pair only); 1=random one; 2=all(large output file size), default=1.
-q  <int>   quality threshold in trimming 3'end of reads, 0-40, default=0. (no trim)
-z  <int>   base quality, default=33 [Illumina is using 64, Sanger Institute is using 33]
-f  <int>   filter low-quality reads containing >n Ns, default=5
-p  <int>   number of processors to use, default=CPU cores detected (up to 8 threads). 
            The parallel performance scales well with 12 threads or less, no significant speed gain for >12 threads.
-x  <int>   max insertion size for pair end mapping, default=500
-m  <int>   min insertion size for pair end mapping, default=28
-L  <int>   mapping the first N nucleotide of the read, default: 0 (map the whole read).
-I  <int>   index interval (1~16), meaning the reference genome will be indexed every Nbp, default=4. (WGBS mode)
            For RRBS mode, index_interval is fixed to 1bp and this command line option is neglected.
            larger index interval uses memory, and slightly reduces mapping sensitivity. (~0.5% difference) 
            for human genome, -I 16 uses ~5GB, compared with ~9GB at the default -I 4.
-A  <str>   set the adapter sequence(s) and trim from 3'end of reads, default=none, requires at least 4nt matched, no mismatch allowed.
            Multiple -A options could be specified to set more than one adapter sequences, i.e. in pair-end sequencing case. 
            default: none (no adapter trimming)
-R          include the reference sequences as the XR:Z:<string> field in SAM output. default=do not include.
-H          do not print header information in SAM format output
-u          report unmapped reads, default=do not report.
-B  <int>   start from the Nth read or read pair, default: 1.
-E  <int>   end at the Nth read or read pair, default: -1 (map all reads until input ends).
            Using -B and -E options user can specify part of the input file to be mapped, so that the input file 
            could be divided into several parts and mapped parallely over distributed system, without creating temporary files. 
-D  <str>   set restriction enzyme digestion site and activate reduced representation bisulfite mapping mode (RRBS mode), can be used multiple times for multiple enzyme digestion.
            i.e. reads must be mapped to digestion sites, the digestion site must be palindromic, digestion position is marked by '-', 
            for example: '-D C-CGG' for MspI digestion, -D T-CGA for TaqI digestion.
            default: none, meaning whole genome shot gun mapping (WGBS mode).
-S  <int>   seed for random number generation in selecting multiple hits.  default: 0 (seed set from system clock).
            other seed values generate pseudo random number based on read index number, so that mapping results are reproducible. 
-n  [0,1]   set mapping strand information:
            -n 0: only map to 2 forward strands, i.e. BSW(++) and BSC(-+)    (i.e. the "Lister protocol")
            for PE sequencing, map read#1 to ++ and -+, read#2 to +- and --. 
            -n 1: map SE or PE reads to all 4 strands, i.e. ++, +-, -+, --    (i.e. the "Cokus protocol")
            default: -n 0. Most bisulfite sequencing data is generated only from forward strands.
-M  <str>   set the alignment information for the additional nucleotide transition. <str> is in the form of two different nucleotides, 
            the first one in the reads could be mapped to the second one in the reference sequences.
            default: -M TC, corresponds to C=>U(T) transition in bisulfite conversion.
            example: -M GA could be used to detect to A=>I(G) transition in RNA editing. 
-V  [0,1,2] verbose level: 0=no message displayed (quiet mode); 1=major message (default); 2=detailed message.
-h          help


4. Output

4.1 BSP format, includes the following tab delimited fields:

    1) id: read ID
    2) seq: mapped read sequence
    3) qual: quality scores
    4) map_flag: 
        UM: unique map (unique pair for paired mapping).
        MA: multiple map (multiple pair for paired mapping)
        OF: over map (#multiple map exceeds MAXHITS)
        NM: no map
        QC: low quality reads
    5) ref: reference sequence name (chromosome name)
    6) ref_loc: mapping location(1 based, 5'-end coordinates of the mapping region on the Watson strand of reference)
    7) strand: 
        ++: forward strand of Watson of reference (BSW)
        +-: reverse strand of Watson of reference (BSWC)
        -+: forward strand of Crick of reference (BSC)  
        --: reverse strand of Crick of reference (BSCC) 
    8) ins_size: insertion size for pair-end mapping, measured by the total nucleotide of the pair-end segment. 
        (5'end to 3'end length of the DNA fragment). 0 means single-end or unpaired mapping.
    9) refseq: Waston reference sequence at the mapping location, with two
        flanking nucleotides in lower cases on each end.
    10) for ungapped hits:  #mismatches
        for gapped hits:    #mismatches:#gap_size:gap_position
            gap_size > 0: insertion on reads
            gap_size < 0: deletion on reads 
    11) mismatch_info:  #hits of 0 mismatch to #hits of max_mismatches, separated by ':'

4.2 SAM format
    FLAG field: 
        UM: 0x0
        MA: 0x100 (non-unique hits)
        OF: 0x100 (non unique hits)
        NM: 0x4
        QC: 0x204
        for mapping on BSC(-+) or BSWC(+-): FLAG=FLAG+0x10, 
            meaning read sequence is reverse complimentarized of the raw reads
        flag 0x400 is not used.
    for pair-end mapping:
        FLAG=FLAG+0x1
            if it's from read set #1, FLAG=FLAG+0x40
            if it's from read set #2, FLAG=FLAG+0x80
        if mappings are paired, FLAG=FLAG+0x2
            if mate is unmapped, FLAG=FLAG+0x8
            if mate is mapped on BSC(-+) or BSWC(+-), FLAG=FLAG+0x20
    aux field: 
        ZS:Z:<strand info> same as BSP column 6).
        XR:Z:<reference sequence> same as BSP column 8).
        NM:i:<#mismatches> same as BSP column 9).    
        ZP:i:<int> RRBS fragment start location, only for RRBS mode.
        ZL:i:<int> RRBS fragment size, only for RRBS mode.
    

    for more details, please refer to SAM format specification: 
    http://samtools.sourceforge.net/SAM1.pdf

Note: all read sequences are recorded as the corresponding sequence following the reference Watson strand direction.


5. Speed and sensitivity
    
The longer seed size(option -s), the faster speed. With seed size increase every bp, mapping time reduces by ~1.5-fold. 
On the other hand, the max number of mismatches that could be detected with 100% sensitivity is bounded by the 
seed_size.  

    max_mismatches_with_100%_sensitivity = (read_len+1-index_interval) / seed_size - 1

If the -v option set max mismatches larger than this number, those mappings with larger max mismatches may not be 
guaranteed to be detected. 

In case full sensitivity can not be achieved within feasible time, user will need to make a decision on the trade off 
between the speed and sensitivity by setting the optimal seed size.  


6 Examples

    single_end: (allow up to 100 multiple hits [-w 100], map to all 4 possible strands [-n 1])
    $bsmap -a SE_read.fastq.gz -d ~/ref/hg19/hg19.fa -o out.sam -n 1 -w 100 -p 8

    single_end: (allow gap with up to 2 nucleotides)
    $bsmap -a SE_read.fastq.gz -d ~/ref/hg19/hg19.fa -o out.sam -w 1000 -g 2

    pair_end: (set -b option)
    $ bsmap -a read1.fq -b read2.fq -d ~/ref/hg19/hg19.fa -o out_pair.bsp -2 out_upair.bsp -p 8 -w 100
    $ bsmap -a PE_reads.bam -b PE_reads.bam -d ~/ref/hg19/hg19.fa -o out.sam -p 8 -w 100  -v 0.07 -m 50 -x 300

    write to STDOUT and use pipe to convert to BAM format output (recommended way to get BAM output)
    $ bsmap -a PE_reads.bam -b PE_reads.bam -d ~/ref/hg19/hg19.fa -p 4 -v 5 | samtools view -bS - > out.bam

	write to STDOUT and use pipe to get methratio file simultaneously, and save a copy of alignment results.
	$ bsmap -a PE_reads.bam -b PE_reads.bam -d ~/ref/hg19/hg19.fa | methratio -d ~/ref/hg19/hg19.fa -o methratio.txt -O alignment.bam -

    trim adapter sequence from 3'end: (set -A option, can be used more than once)
    $ bsmap -a SE_read.bam -d ~/ref/hg19/hg19.fa -o out.bam -p 8 -w 100 -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA

    set max #mismatch to readlen * 5%: (set -v option between 0 and 1)
    $bsmap -a SE_read.bam -d ~/ref/hg19/hg19.fa -o out.bam -n 1 -w 100 -v 0.05

    set max #mismatch to 5: (set -v option NOT between 0 and 1)
    $bsmap -a SE_read.bam -d ~/ref/hg19/hg19.fa -o out.bam -n 1 -w 100 -v 5

    mapping from read pair #10001 to read pair #20000 in the input file: (set -B and -E option)
    $ bsmap -a PE_reads.bam -b PE_reads.bam -d ~/ref/hg19/hg19.fa -o out.bam -w 100  -v 5 -B 10001 -E 20000

    using Illumina quality: (set -z option)
    $ bsmap -a PE_read1.fq -b PE_read2.fq -d ~/ref/hg19/hg19.fa -o out.bam -z 64
    
    report only uniquely mapped reads/read pairs: (set -r option)
    $ bsmap -a reads.fastq -d hg19.fa -o out.bsp -r 0

    trim low quality 3'end: (set -q option)
    $ bsmap -a PE_reads.bam -b PE_reads.bam -d ~/ref/hg19/hg19.fa -o out.bam -q 2

    detect A=>G editing in RNA_seq instead of C=>T conversion in bisulfite sequencing (set -M option)
    $ bsmap -a reads.bam -d RNA_ref.fa -M GA -o out.bam

	RRBS mode: (set -D option to specify digestion site information, can be used more than once for multiple enzyme digestion)
	$ bsmap -a PE_reads.bam -b PE_reads.bam  -d ~/ref/hg19/hg19.fa -o out.bam -p 8 -w 100 -v 5 -D C-CGG -D G-CWGC 


7. Scripts: 

7.1 methratio.py
python script to extract methylation ratios from BSMAP mapping results. Require python 2.X. 
For human genome, methratio.py needs ~26GB memory.  
For systems with limited memory, user can set the -c/--chr option to process specified chromosomes only,
and combine results for all chromosomes afterwards.

Usage: python methratio.py [options] BSMAP_MAPPING_FILES

BSMAP_MAPPING_FILES could be one or more output files from BSMAP.
The format will be determined by the filename suffix. 
(SAM format for *.sam and *.bam, BSP format for other filenames.)

Options:
  -h, --help            show this help message and exit
  -o FILE, --out=FILE   output file name, if it's omitted, the output will be written to STDOUT for downstream pipes.
  -d FILE, --ref=FILE   reference genome fasta file. (required)
  -c CHR, --chr=CHR     process only specified chromosomes. [default: all]
                        example: --chr=chr1,chr2 (this uses ~4.5GB compared with ~26GB for the whole genome)
  -s PATH, --sam-path=PATH
                        path to samtools. [default: none]
  -u, --unique          process only unique mappings/pairs.
  -p, --pair            process only properly paired mappings.
  -z, --zero-meth       report loci with zero methylation ratios.
  -q, --quiet           don't print progress on stderr.
  -r, --remove-duplicate
                        remove duplicated mappings to reduce PCR bias. (This option should not be used on RRBS data. For WGBS, sometimes 
                        it's hard to tell if duplicates are caused by PCR due to high seqeuncing depth.)
  -t N, --trim-fillin=N
                        trim N fill-in nucleotides in DNA fragment end-repairing. [default:0] 
                        (This option is only for pair-end mapping. For RRBS, N could be detetmined by the distance between
                        cuttings sites on forward and reverse strands. For WGBS, N is usually between 0~3.) 
  -g, --combine-CpG     combine CpG methylaion ratio from both strands. [default: False]
  -m FOLD, --min-depth=FOLD
                        report loci with sequencing depth>=FOLD. [default: 1]
  -n, --no-header       don't print a header line
  -i CT_SNP, --ct-snp=CT_SNP
                        how to handle CT SNP ("no-action", "correct", "skip"),
                        default: "correct".
                        "correct":      correct the methylation ratio according to the C/T SNP information
                        estimated by the G/A counts on reverse strand, see the output format below for details.
                        "skip":         do not report loci with C/T SNP detected (i.e. detected A on reverse strand)
                        "no-action":    do not consider C/T SNP. this uses half the memory compared to the other two options.
  -x TYPE, --context=TYPE
                        methylation pattern type [CG|CHG|CHH], multiple types separated by ','. [default: all]


Output format: tab delimited txt file with the following columns:
    1) chromorome
    2) coordinate (1-based)
    3) strand
    4) sequence context (CG|CHG|CHH)
    5) methylation ratio, calculated as #C_counts / #eff_CT_counts
    6) number of effective total C+T counts on this locus (#eff_CT_counts) 
            CT_SNP="no action", #eff_CT_counts = #CT_counts
            CT_SNP="correct", #eff_CT_counts = #CT_counts * (#rev_G_counts / #rev_GA_counts)
    7) number of total C counts on this locus (#C_counts)
    8) number of total C+T counts on this locuso (#CT_counts)
    9) number of total G counts on this locus of reverse strand (#rev_G_counts)
    10) number of total G+A counts on this locus of reverse strand (#rev_GA_counts)
    11) lower bound of 95% confidence interval of methylation ratio, calculated by Wilson score interval for binomial proportion.
    12) upper bound of 95% confidence interval of methylation ratio, calculated by Wilson score interval for binomial proportion.

Example:
    python methratio.py --chr=chr1,chr2 --ref=hg19.fa --out=methratio.txt rrbsmap_sample*.sam
    python methratio.py -d mm9.fa -o meth.txt -p bsmap_sample1.bsp bsmap_sample2.sam bsmap_sample3.bam 
    python methratio.py -s /home/tools/samtools -t 1 -d arab.fa -o meth.txt bsmap_sample.sam

Note: For overlapping paired hits, nucleotides in the overlapped part should be counted only once instead of twice.
methratio.py can correctly handle such cases for SAM format output, but for BSP format it will still be counted twice,
because the BSP format does not contain mapping information of the mate.

7.2 methdiff.py
python script for differential methylation analysis from methratio.py results. Require python 2.X. 
The differntial test is based on whether the confidence interval of methylation ratio in two groups overlap, this provides a conservative estimation of the p-values.
The C_count/CT_count for replicates within the same sample group will be added up together, 
the methylation ratio can be calculated for a predefined bin size (-b option), single nucleotide methylation ratio is equivalent to -b 1. 

Usage: methdiff.py [options] <GROUP1_SAMPLE1,GROUP1_SAMPLE2,...> <GROUP2_SAMPLE1,GROUP2_SAMPLE2,...>

Options:
  -h, --help            show this help message and exit
  -o FILE, --out=FILE   output differential methylation file name. (required)
  -d FILE, --ref=FILE   reference genome fasta file. (required)
  -b BIN, --bin=BIN     bin size. [default: 100]
  -p PVAL, --pval=PVAL  p-value cut-off. [default: 0.01]
  -r DIFF, --diff=DIFF  minimal abs meth ratio difference. [default: 0.33]
  -x TYPE, --context=TYPE
                        methylation pattern type [CG|CHG|CHH], multiple types separated by ','. [default: all]
  -l LABELS, --labels=LABELS
                        output label for each group, separated by ','. [default: filenames]
  -m FOLD, --min-depth=FOLD
                        minimal average coverage. [default: 1]
  -s STRAND, --strand=STRAND
                        which strand to process, [both|forward|reverse]. [default: both]
  -q, --quiet           don't print progress on stderr.


7.3 sam2bam.sh
Shell script to convert SAM format to sorted and indexed BAM format.
The input SAM file will be deleted if the conversion is successful.
This script is automatically called by BSMAP if the output file has .bam suffix, which requires samtools and sam2bam.sh installed in system default path.  
It can also be used manually.

Usage: ./sam2bam.sh INPUT_SAM_FILE

Example: ./sam2bam.sh sample1.sam 
This will generate sorted BAM file sample1.bam and index file sample1.bam.bai. 


7.4 bsp2sam.py
Python script to convert .BSP format output to SAM format output. 

Usage: bsp2sam.py [options] BSP_FILE

Options:
  -h, --help           show this help message and exit
  -o FILE, --out=FILE  output file name. (required) 
  -d FILE, --ref=FILE  reference genome fasta file. (required)
  -q, --quiet          don't print progress on stderr.

Example: python bsp2sam.py -d hg19.fa -o sample.sam sample.bsp

Note: For pair-end BSP output, this script will conver it into single-end SAM output, i.e. the mapping information remains, but the pairing information lost. 


8. Citation
    Yuanxin Xi and Wei Li, "BSMAP: whole genome bisulfite sequence MAPping program" (2009) BMC Bioinformatics 2009, 10:232


9. Contact
    Yuanxin Xi
    Bioinformatics Division, 
    Dan L. Duncan Cancer Center,
    Baylor College of Medicince, 
    Houston, TX 77030, USA
    713-798-6254
    yxi@bcm.tmc.edu, xiyuanxin@gmail.com
    
