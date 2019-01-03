# BSMAPz

Greg Zynda's modified fork of BSMAP

BSMAPz began as a fork of [BSMAP v2.90](https://sites.google.com/a/brown.edu/bioinformatics-in-biomed/bsmap-for-methylation) at commit 4b164c833fff23d455b4f9b86834f829954fadb4, which is still hosted on the now deprecated [Google Code repository](https://code.google.com/archive/p/bsmap/).
That code base has gone stale and has been incompatible with modern versions of [samtools](http://www.htslib.org/).

To reduce confusion, I have named my modified version BSMAPz, where the Z signifies it is the version modified by Greg Zynda.

## Introduction

BSMAPz is a short reads mapping program for bisulfite sequencing in DNA methylation study.
Bisulfite treatment coupled with next generation sequencing could estimate the methylation ratio of every single Cytosine location in the genome by mapping high throughput bisulfite reads to the reference sequences.

Bisulfite mapping is different from usual sequence mapping in two aspects:

1. The additional C/T mapping is asymmetric, a T in the read could be aligned to C in the reference but not vice versa
2. The Watson and Crick strand are not complimentary after bisulfite treatment.

Each read need to be compared with 4 reference sequences:

1. BSW (bisulfite Watson)
2. BSWC (reverse complimentary of BSW)
3. BSC (bisulfite Crick)
4. BSCC (reverse complimentary of BS)

BSMAPz is designed to be a general-purpose mapping program to handle these special characteristics of bisulfite mapping.
It is based on the open source program [SOAPv1](http://soap.genomics.org.cn/soap1/) (Short Oligo Alignment Program).  

**Main features:**

* read length up to 144nt, allow up to 15 mismatches, 1 continous gap up to 3nt
* support pair end mapping, support parallel mapping
* support SAM format input/output, support gzipped FASTA/FASTQ format input  
* support both whole genome (WGBS) and reduced representation bisulfite sequencing (RRBS)
* support trimming adapters and low quality sequences from 3'end of reads
* allow different running modes with flexible speed/memory/sensitivity to run on different hardware configurations
* include script to extract methylation ratios

BSMAPz is under [GNU Public License (GPL)](GPL_3.0.txt).

## Installation

BSMAPz is designed for linux64 platform, and has been tested on

* Mac OSX
* Linux x86-64

Requirements:

* samtools (on PATH)
* pysam

First, clone the source code:
```
git clone https://github.com/zyndagj/BSMAPz.git
```
and *optionally* check out a specific **commit** or **release**:
```
git checkout [release]
```

Next, compile the binary:
```
make bsmapz
```

Test the code:
```
make test
```
all tests should result in an "OK" status. If not, please submit an issue describing your system.

Install the binary into system default path: (optional)
```
make install
```
or a custom location
```
mkdir -p /opt/bsmapz
make DESTDIR=/opt/bsmapz install
```

## Usage
```
Usage:	bsmapz [options]
       -a  <str>   query a file, FASTA/FASTQ/BAM format
       -d  <str>   reference sequences file, FASTA format
       -o  <str>   output alignment file, BSP/SAM/BAM format, if omitted, the output will be written to STDOUT in SAM format.

  Options for alignment:
       -s  <int>   seed size, default=16(WGBS mode), 12(RRBS mode). min=8, max=16.
       -v  <float> if this value is between 0 and 1, it's interpreted as the mismatch rate w.r.t to the read length.
                   otherwise it's interpreted as the maximum number of mismatches allowed on a read, <=15.
                   example: -v 5 (max #mismatches = 5), -v 0.1 (max #mismatches = read_length * 10%)
                   default=0.08.
       -g  <int>   gap size, BSMAPz only allows 1 continuous gap (insert or deletion) with up to 3 nucleotides
                   default=0
       -w  <int>   maximum number of equal best hits to count, <=10000
       -3          using 3-nucleotide mapping approach
       -B  <int>   start from the Nth read or read pair, default: 1
       -E  <int>   end at the Nth read or read pair, default: 4,294,967,295
       -I  <int>   index interval, default=4
       -k  <float> set the cut-off ratio for over-represented kmers, default=5e-07
                   example: -k 1e-6 means the top 0.0001% over-represented kmer will be skipped in alignment
       -p  <int>   number of processors to use, default=8
       -D  <str>   activating RRBS mapping mode and set restriction enzyme digestion sites.
                   digestion position marked by '-', example: -D C-CGG for MspI digestion.
                   default: none (whole genome shotgun bisulfite mapping mode)
       -S  <int>   seed for random number generation used in selecting multiple hits
                   other seed values generate pseudo random number based on read index number, to allow reproducible mapping results.
                   default=0. (get seed from system clock, mapping results not resproducible.)
       -n  [0,1]   set mapping strand information. default: 0
                   -n 0: only map to 2 forward strands, i.e. BSW(++) and BSC(-+),
                   for PE sequencing, map read#1 to ++ and -+, read#2 to +- and --.
                   -n 1: map SE or PE reads to all 4 strands, i.e. ++, +-, -+, --
       -M  <str>   set alignment information for the additional nucleotide transition.
                   <str> is in the form of two different nucleotides N1N2,
                   indicating N1 in the reads could be mapped to N2 in the reference sequences.
                   default: -M TC, corresponds to C=>U(T) transition in bisulfite conversion.
                   example: -M GA could be used to detect A=>I(G) transition in RNA editing.

  Options for trimming:
       -q  <int>   quality threshold in trimming, 0-40, default=0 (no trim)
       -z  <int>   base quality, default=33 [Illumina is using 64, Sanger Institute is using 33]
       -f  <int>   filter low-quality reads containing >n Ns, default=5
       -A  <str>   3-end adapter sequence, default: none (no trim)
       -L  <int>   map the first N nucleotides of the read, default:160 (map the whole read).

  Options for reporting:
       -r  [0,1,2] how to report repeat hits, 0=none(unique hit/pair); 1=random one; 2=all(slow), default:1.
       -R          print corresponding reference sequences in SAM output, default=off
       -u          report unmapped reads, default=off
       -H          do not print header information in SAM format output
       -V  [0,1,2] verbose level: 0=no message displayed (quiet mode);
                   1=major message (default); 2=detailed message.

  Options for pair-end alignment:
       -b  <str>   query b file
       -m  <int>   minimal insert size allowed, default=28
       -x  <int>   maximal insert size allowed, default=1000

       -h          help
```

## Output

### BSP format

The BSP format includes the following tab delimited fields:

1. **id**: read ID
2. **seq**: mapped read sequence
3. **qual**: quality scores
4. **map_flag**:
   * UM: unique map (unique pair for paired mapping).
   * MA: multiple map (multiple pair for paired mapping)
   * OF: over map (#multiple map exceeds MAXHITS)
   * NM: no map
   * QC: low quality reads
5. **ref**: reference sequence name (chromosome name)
6. **ref_loc**: mapping location(1 based, 5'-end coordinates of the mapping region on the Watson strand of reference)
7. **strand**:
   * `++`: forward strand of Watson of reference (BSW)
   * `+-`: reverse strand of Watson of reference (BSWC)
   * `-+`: forward strand of Crick of reference (BSC)  
   * `--`: reverse strand of Crick of reference (BSCC)
8. **ins_size**: insertion size for pair-end mapping, measured by the total nucleotide of the pair-end segment. (5'end to 3'end length of the DNA fragment). 0 means single-end or unpaired mapping.
9. **refseq**: Waston reference sequence at the mapping location, with two flanking nucleotides in lower cases on each end.
10. **for ungapped hits**:  #mismatches

    **for gapped hits**:    #mismatches:#gap_size:gap_position
    * gap_size > 0: insertion on reads
    * gap_size < 0: deletion on reads
11. **mismatch_info**:  #hits of 0 mismatch to #hits of max_mismatches, separated by `:`

### SAM format

FLAG field:
```
 UM: 0x0
 MA: 0x100 (non-unique hits)
 OF: 0x100 (non unique hits)
 NM: 0x4
 QC: 0x204

for mapping on BSC(-+) or BSWC(+-):
    FLAG=FLAG+0x10, meaning read sequence is the reverse compliment of the raw reads flag 0x400 is not used.

for pair-end mapping:
    FLAG=FLAG+0x1
        if it's from read set #1, FLAG=FLAG+0x40
        if it's from read set #2, FLAG=FLAG+0x80
    if mappings are paired, FLAG=FLAG+0x2
        if mate is unmapped, FLAG=FLAG+0x8
        if mate is mapped on BSC(-+) or BSWC(+-), FLAG=FLAG+0x20
```
aux field:
```
 ZS:Z:<strand info> same as BSP column 6).
 XR:Z:<reference sequence> same as BSP column 8).
 NM:i:<#mismatches> same as BSP column 9).    
 ZP:i:<int> RRBS fragment start location, only for RRBS mode.
 ZL:i:<int> RRBS fragment size, only for RRBS mode.
```

For more details, please refer to [SAM format specification](http://samtools.sourceforge.net/SAM1.pdf)

> Note: all read sequences are recorded as the corresponding sequence following the reference Watson strand direction.

## Speed and sensitivity

The longer seed size (`-s`), the faster speed. With seed size increase every bp, mapping time reduces by ~1.5-fold.
On the other hand, the max number of mismatches that could be detected with 100% sensitivity is bounded by the `seed_size`.

```
    max_mismatches_with_100%_sensitivity = (read_len + 1 - index_interval) / seed_size - 1
```

If the `-v` option set max mismatches larger than this number, those mappings with larger max mismatches may not be guaranteed to be detected.

In case full sensitivity can not be achieved within feasible time, user will need to make a decision on the trade off between the speed and sensitivity by setting the optimal seed size.

## Examples

### Single end input

Allow up to 100 multiple hits [`-w 100`] and map to all 4 possible strands [`-n 1`]
```
bsmapz -a SE_read.fastq.gz -d ~/ref/hg19/hg19.fa -o out.sam -n 1 -w 100 -p 8
```

Allow gaps with up to 2 nucleotides
```
bsmapz -a SE_read.fastq.gz -d ~/ref/hg19/hg19.fa -o out.sam -w 1000 -g 2
```

Trim adapter sequence from 3'end: (set -A option, can be used more than once)
```
bsmapz -a SE_read.bam -d ~/ref/hg19/hg19.fa -o out.bam -p 8 -w 100 -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA
```

Set max #mismatch to readlen * 5%: (set -v option between 0 and 1)
```
bsmapz -a SE_read.bam -d ~/ref/hg19/hg19.fa -o out.bam -n 1 -w 100 -v 0.05
```

Set max #mismatch to 5: (set -v option NOT between 0 and 1)
```
bsmapz -a SE_read.bam -d ~/ref/hg19/hg19.fa -o out.bam -n 1 -w 100 -v 5
```

Report only uniquely mapped reads/read pairs: (set -r option)
```
bsmapz -a reads.fastq -d hg19.fa -o out.bsp -r 0
```

Detect A=>G editing in RNA_seq instead of C=>T conversion in bisulfite sequencing (set -M option)
```
bsmapz -a reads.bam -d RNA_ref.fa -M GA -o out.bam
```

### Paired end input
Use paired fastq inputs
```
bsmapz -a read1.fq -b read2.fq -d ~/ref/hg19/hg19.fa -o out_pair.bsp -2 out_upair.bsp -p 8 -w 100
```

Use paired BAM inputs (same file)
```
bsmapz -a PE_reads.bam -b PE_reads.bam -d ~/ref/hg19/hg19.fa -o out.sam -p 8 -w 100  -v 0.07 -m 50 -x 300
```

Write to STDOUT and use pipe to convert to BAM format output (recommended way to get BAM output)
```
bsmapz -a PE_reads.bam -b PE_reads.bam -d ~/ref/hg19/hg19.fa -p 4 -v 5 | samtools view -bS - > out.bam
```

Write to STDOUT and use pipe to get methratio file simultaneously, and save a copy of alignment results.
```
bsmapz -a PE_reads.bam -b PE_reads.bam -d ~/ref/hg19/hg19.fa | methratio -d ~/ref/hg19/hg19.fa -o methratio.txt -O alignment.bam -
```

Mapping from read pair #10001 to read pair #20000 in the input file: (set -B and -E option)
```
bsmapz -a PE_reads.bam -b PE_reads.bam -d ~/ref/hg19/hg19.fa -o out.bam -w 100  -v 5 -B 10001 -E 20000
```

Use Illumina quality: (set `-z` option)
```
bsmapz -a PE_read1.fq -b PE_read2.fq -d ~/ref/hg19/hg19.fa -o out.bam -z 64
```

Trim low quality 3'end: (set -q option)
```
bsmapz -a PE_reads.bam -b PE_reads.bam -d ~/ref/hg19/hg19.fa -o out.bam -q 2
```

RRBS mode: (set -D option to specify digestion site information, can be used more than once for multiple enzyme digestion)
```
bsmapz -a PE_reads.bam -b PE_reads.bam  -d ~/ref/hg19/hg19.fa -o out.bam -p 8 -w 100 -v 5 -D C-CGG -D G-CWGC
```

## Scripts:

### methratio.py
A python script to extract methylation ratios from BSMAPz mapping results.

Requires:
* python 2.*
* samtools
* pysam

To extract methylation ratios on the human genome, methratio.py needs about 20GB of memory.
While methratio.py will no longer overallocate memory and crash, you can also manually limit memory usage with the `-M` parameter when packing concurrent runs.

```
usage: methratio.py [-h] [-o FILE] [-w FILE] [-b BIN] -d FILE [-c CHR] [-u]
                    [-p] [-z] [-q] [-r] [-t N] [-g] [-m FOLD] [-n] [-i CT_SNP]
                    [-x TYPE] [-f] [-M MB] [-N NP]
                    FILES [FILES ...]

Calls single-base methylation ratios by context.

positional arguments:
  FILES                 Files from BSMAP output [BAM|SAM|BSP]

optional arguments:
  -h, --help            show this help message and exit
  -o FILE, --out FILE   output methylation ratio file name. [default: STDOUT]
  -w FILE, --wig FILE   output methylation ratio wiggle file. [default: none]
  -b BIN, --wig-bin BIN
                        wiggle file bin size. [default: 25]
  -d FILE, --ref FILE   reference genome fasta file. (required)
  -c CHR, --chroms CHR  process only specified chromosomes, separated by ','.
                        [default: all] example: --chroms=chr1,chr2
  -u, --unique          process only unique mappings/pairs.
  -p, --pair            process only properly paired mappings.
  -z, --zero-meth       report loci with zero methylation ratios. (deprecated,
                        -z will be always enabled)
  -q, --quiet           don't print progress on stderr.
  -r, --remove-duplicate
                        remove duplicated reads.
  -t N, --trim-fillin N
                        trim N end-repairing fill-in nucleotides from
                        fragments. [default: 0]
  -g, --combine-CpG     combine CpG methylaion ratios on both strands.
  -m FOLD, --min-depth FOLD
                        report loci with sequencing depth>=FOLD. [default: 1]
  -n, --no-header       don't print a header line
  -i CT_SNP, --ct-snp CT_SNP
                        how to handle CT SNP ("no-action", "correct", "skip"),
                        default: "correct".
  -x TYPE, --context TYPE
                        methylation pattern type [CG|CHG|CHH], multiple types
                        separated by ','. [default: all]
  -f, --full            Report full context (CHG -> CAG)
  -M MB, --mem MB       Maximum memory in megabytes to use [16000]
  -N NP, --np NP        Maximum number of processes to use [-1]
```

#### Output format

A tab delimited txt file with the following columns:

1. **chr** - chromorome
2. **pos** - coordinate (1-based)
3. **strand** - reference strand (+/-)
4. **context** - methylation context (CG|CHG|CHH)
5. **ratio** - methylation ratio, calculated as `#C_counts / #eff_CT_counts`
6. **eff_CT_counts** - number of effective total C+T counts on this locus
   * `CT_SNP="no action"`, then `#eff_CT_counts = #CT_counts`
   * `CT_SNP="correct"`, then `#eff_CT_counts = #CT_counts * (#rev_G_counts / #rev_GA_counts)`
7. **C_counts** - number of total C counts on this locus
8. **CT_counts** - number of total C+T counts on this locus
9. **rev_G_counts** - number of total G counts on the reverse strand at this locus
10. **rev_GA_counts** - number of total G+A counts on this locus of reverse strand
11. **CI_lower** - lower bound of 95% confidence interval of methylation ratio, calculated by Wilson score interval for binomial proportion.
12. **CI_upper** - upper bound of 95% confidence interval of methylation ratio, calculated by Wilson score interval for binomial proportion.

#### Example:

```
python methratio.py --chr=chr1,chr2 --ref=hg19.fa --out=methratio.txt rrbsmap_sample*.sam

python methratio.py -d mm9.fa -o meth.txt -p bsmap_sample1.bsp bsmap_sample2.sam bsmap_sample3.bam

python methratio.py -s /home/tools/samtools -t 1 -d arab.fa -o meth.txt bsmap_sample.sam
```

> Note: The original version of methratio.py counted overlapping regions of reads twice when using BSP input. We recommend re-running all methylation extraction if you utilize the BSP format.

When read pairs overlap, methratio.py will trim off the overlapping region of the read with the lower index.
In the future, we hope to analyze both pairs at once so we can employ the `samtools mpileup` method like [MethylDackel](https://github.com/dpryan79/MethylDackel#a-note-on-overlapping-reads).

### methdiff.py

A python script for differential methylation analysis from `methratio.py` results.

Requires:
* python 2.X.

The differential test is based on whether the confidence interval of methylation ratio in two groups overlap, this provides a conservative estimation of the p-values.
The `C_count/CT_count` for replicates within the same sample group will be added up together, the methylation ratio can be calculated for a predefined bin size (`-b`), single nucleotide methylation ratio is equivalent to `-b 1`.

```
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
```

### sam2bam.sh

Shell script to convert SAM format to sorted and indexed BAM format.
The input SAM file will be deleted if the conversion is successful.
This script is automatically called by BSMAPz if the output file has `.bam` suffix, which requires samtools and sam2bam.sh installed in system default path.

It can also be used manually.
```
Usage: ./sam2bam.sh INPUT_SAM_FILE

Example: ./sam2bam.sh sample1.sam
This will generate sorted BAM file sample1.bam and index file sample1.bam.bai.
```

### bsp2sam.py

A python script to convert .BSP format output to SAM format output.

```
Usage: bsp2sam.py [options] BSP_FILE

Options:
  -h, --help           show this help message and exit
  -o FILE, --out=FILE  output file name. (required)
  -d FILE, --ref=FILE  reference genome fasta file. (required)
  -q, --quiet          don't print progress on stderr.

Example: python bsp2sam.py -d hg19.fa -o sample.sam sample.bsp
```

> Note: For pair-end BSP output, this script will conver it into single-end SAM output, i.e. the mapping information remains, but the pairing information lost.

## Citation

Greg Zynda. "BSMAPz." GitHub Repository. https://github.com/zyndagj/BSMAPz. 2018.

Yuanxin Xi and Wei Li, "BSMAP: whole genome bisulfite sequence MAPping program" (2009) BMC Bioinformatics 2009, 10:232
