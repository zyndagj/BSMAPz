#! /usr/bin/env python

###############################################################################
# Author: Greg Zynda
# Last Modified: 10/05/2018
###############################################################################
# BSD 3-Clause License
# 
# Copyright (c) 2018, Texas Advanced Computing Center
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# 
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###############################################################################

import sys, time, os, array, argparse, re
from itertools import compress, ifilter, izip
import subprocess as sp
import multiprocessing as mp
import multiprocessing.pool
import pysam
from string import maketrans

quiet = False

def main():
	usage = "usage: %prog [options] BSMAP_MAPPING_FILES"
	parser = argparse.ArgumentParser(description="Calls single-base methylation ratios by context.")
	parser.add_argument("-o", "--out", dest="outfile", metavar="FILE", help="output methylation ratio file name. [default: STDOUT]", default="")
	parser.add_argument("-w", "--wig", dest="wigfile", metavar="FILE", help="output methylation ratio wiggle file. [default: none]", default="")
	parser.add_argument("-b", "--wig-bin", dest="wigbin", type=int, metavar='BIN', help="wiggle file bin size. [default: %(default)s]", default=25)
	parser.add_argument("-d", "--ref", dest="reffile", metavar="FILE", help="reference genome fasta file. (required)", required=True)
	parser.add_argument("-c", "--chroms", metavar="CHR", help="process only specified chromosomes, separated by ','. [default: all]\nexample: --chroms=chr1,chr2", default=[])
	parser.add_argument("-u", "--unique", action="store_true", dest="unique", help="process only unique mappings/pairs.", default=False)
	parser.add_argument("-p", "--pair", action="store_true", dest="pair", help="process only properly paired mappings.", default=False)
	parser.add_argument("-z", "--zero-meth", action="store_true", dest="meth0", help="report loci with zero methylation ratios. (deprecated, -z will be always enabled)", default=True)
	parser.add_argument("-q", "--quiet", action="store_true", dest="quiet", help="don't print progress on stderr.", default=False)
	parser.add_argument("-r", "--remove-duplicate", action="store_true", dest="rm_dup", help="remove duplicated reads.", default=False)
	parser.add_argument("-t", "--trim-fillin", dest="trim_fillin", type=int, metavar='N', help="trim N end-repairing fill-in nucleotides. [default: %(default)s]", default=0)
	parser.add_argument("-g", "--combine-CpG", action="store_true", dest="combine_CpG", help="combine CpG methylaion ratios on both strands.", default=False)
	parser.add_argument("-m", "--min-depth", dest="min_depth", type=int, metavar='FOLD', help="report loci with sequencing depth>=FOLD. [default: %(default)s]", default=1)
	parser.add_argument("-n", "--no-header", action="store_true", dest="no_header", help="don't print a header line", default=False)
	parser.add_argument("-i", "--ct-snp", dest="CT_SNP", help='how to handle CT SNP ("no-action", "correct", "skip"), default: "correct".', default="correct")
	parser.add_argument("-x", "--context", dest="context", metavar='TYPE', help="methylation pattern type [CG|CHG|CHH], multiple types separated by ','. [default: all]", default='')
	parser.add_argument("-f", "--full", action="store_true", help="Report full context (CHG -> CAG)")
	parser.add_argument("-M", "--mem", type=int, metavar='MB', help="Maximum memory in megabytes to use [%(default)s]", default=4000)
	parser.add_argument("-N", "--np", type=int, metavar='NP', help="Maximum number of processes to use [%(default)s]", default=-1)
	parser.add_argument("infiles", metavar="FILES", help="Files from BSMAP output [BAM|SAM|BSP]", nargs="+")
	# Parse Options
	options = parser.parse_args()
	# Set up any globals
	quiet = options.quiet
	# Check options
	if len(options.chroms) > 0: options.chroms = options.chroms.split(',')
	CT_SNP_val = {"no-action": 0, "correct": 1, "skip": 2}
	try: options.CT_SNP = CT_SNP_val[options.CT_SNP.lower()]
	except: parser.error('Invalid -i value, select "no-action", "correct" or "skip"')
	if options.min_depth <= 0: parser.error('Invalid -m value, must >= 1')
	if options.trim_fillin < 0: parser.error('Invalid -t value, must >= 0')
	if len(options.context) > 0: options.context = options.context.split(',')
	if len(options.outfile) == 0: disp("Missing output file name, write to STDOUT.")
	if options.mem > memAvail(0.95): sys.exit("Only %i MB available, not %i"%(memAvail(0.95), options.mem))
	if options.np == -1: options.np = min(mp.cpu_count(), 64)
	# Parse fasta.fai for chrom lengths
	if not os.path.exists(options.reffile+'.fai'):
		disp("Indexing reference")
		sp.call(['samtools','faidx',options.reffile])
	chromDict = ParseFai(options.reffile+'.fai', options.chroms)
	sortedChroms = sorted(chromDict.keys())
	# Initialize output
	if options.outfile:
		fout = open(options.outfile, 'w')
	else:
		fout = sys.stdout
	if not options.no_header: 
		fout.write('chr\tpos\tstrand\tcontext\tratio\teff_CT_count\tC_count\tCT_count\trev_G_count\trev_GA_count\tCI_lower\tCI_upper\n')
	if options.outfile: fout.close()
	if options.wigfile: open(options.wigfile, 'w').close()
	# Presort inputs
	disp("Presorting inputs")
	sortedFiles = map(lambda x: sortFile(x, N=options.np, M=options.mem), options.infiles)
	# Calculate memory requirements and limits
	largestChromSize = max(chromDict.values())
	maxMemChrom = 2*largestChromSize*4+largestChromSize
	maxMemChromMB = maxMemChrom/(1024**2)
	if options.mem < maxMemChromMB: sys.exit("methratio.py needs at least %i MB of memory to process this reference\n"%(maxMemChromMB))
	if maxMemChromMB:
		maxChromProcs = min(options.mem/maxMemChromMB, options.np, len(chromDict))
	else:
		maxChromProcs = min(options.np, len(chromDict))
	disp("Processing %i chromosomes at a time"%(maxChromProcs))
	# Create shared array for synchronization
	global syncArray
	syncArray = mp.RawArray('B', [0]*len(chromDict))
	argList = [(chrom, chromDict[chrom], options, sortedFiles, pid) for pid, chrom in enumerate(sortedChroms)]
	if maxChromProcs > 1:
		# Create worker pool
		chromPool = ChromPool(maxChromProcs)
		# Launch workers
		ret = chromPool.map(chromWorker, argList, chunksize=1)
		# Close pool
		chromPool.close()
		chromPool.join()
	else:
		ret = map(chromWorker, argList)
	# Calculate stats
	nmap, nc, nd = (0, 0, 0)
	for sChrom, retList in zip(sortedChroms, ret):
		if not retList:
			sys.exit("%s failed to process\n"%(sChrom))
		cChrom, cNmap, cNc, cNd = retList
		assert(cChrom == sChrom)
		nmap += cNmap
		nc += cNc
		nd += cNd
	# Display stats
	disp('total %i valid mappings, %i covered cytosines, average coverage: %.2f fold.' % (nmap, nc, float(nd)/nc))
	# Delete temporary files
	for f in sortedFiles:
		if 'tmpSrt.bam' in f:
			os.remove(f)
			os.remove(f+'.bai')
		if 'tmpSrt.bsp' in f:
			os.remove(f)

# Can't use daemon processes in the main pool
class NoDaemonProcess(mp.Process):
	def _get_daemon(self):
		return False
	def _set_daemon(self, value):
		pass
	daemon = property(_get_daemon, _set_daemon)
class ChromPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess

def sortFile(infile, N=1, M=1000):
	'''
	Presorts an input file and returns the temprorary file, which should be deleted afterwards
	'''
	fileEXT = infile.split('.')[-1].upper()
	if fileEXT == 'BAM':
		sortedFile = '.'.join(infile.split('.')[:-1]+['tmpSrt','bam'])
		if bamIsSorted(infile):
			disp("%s is already sorted"%(infile))
			return infile
		else:
			samSortMem = int(M/N)
			disp("Calling samtools sort on %s and using %i MB of memory"%(infile, samSortMem))
			sp.check_call('samtools sort -m %iM -@ %i -O bam -o %s -T %s_tmp %s'%(samSortMem, N, sortedFile, sortedFile, infile), shell=True)
			sp.check_call('samtools index %s'%(sortedFile), shell=True)
			return sortedFile
	elif fileEXT == 'SAM':
		sortedFile = '.'.join(infile.split('.')[:-1]+['tmpSrt','bam'])
		samSortMem = int(M/N)
		disp("Calling samtools sort on %s and using %i MB of memory"%(infile, samSortMem))
		sp.check_call('samtools view -uS %s | samtools sort -m %iM -@ %i -O bam -o %s -T %s_tmp'%(infile, samSortMem, N, sortedFile, sortedFile), shell=True)
		sp.check_call('samtools index %s'%(sortedFile), shell=True)
		return sortedFile
	elif fileEXT == 'BSP':
		sortedFile = '.'.join(infile.split('.')[:-1]+['tmpSrt','bsp'])
		disp("Running manual sort on %s using %i MB of memory"%(infile, M))
		try:
			sp.check_call('LC_ALL=C sort --parallel=%i -k5,5 -k6,6n -k1,1n -S %iM %s > %s 2>/dev/null'%(N, M, infile, sortedFile), shell=True)
		except sp.CalledProcessError as e:
			disp("Could not sort %s in parallel. Falling back to single core"%(infile))
			sp.check_call('LC_ALL=C sort -k5,5 -k6,6n -k1,1n -S %iM %s > %s'%(M, infile, sortedFile), shell=True)
		return sortedFile
	else:
		sys.exit("methratio.py does not handle %s files\n"%(fileEXT))

class refcache:
	def __init__(self, FA, chrom, chromLen, cacheSize=50000):
		self.FA = FA
		self.chrom = chrom
		self.chromLen = chromLen
		self.start = 0
		self.cacheSize = cacheSize
		self.end = min(cacheSize, chromLen)
		self.seq = self.FA.fetch(self.chrom, 0, self.end)
	def fetch(self, pos, pos2):
		assert(pos >= self.start)
		if pos2 > self.end:
			assert(pos2 <= self.chromLen)
			self.start = pos
			self.end = pos+self.cacheSize
			self.seq = self.FA.fetch(self.chrom, self.start, self.end)
		sI = pos-self.start
		eI = pos2-self.start
		retseq = self.seq[sI:eI]
		return retseq

def bamIsSorted(inFile):
	if not os.path.exists(inFile+".bai"):
		return False
	if 'coordinate' not in sp.check_output("samtools view -H %s"%(inFile), shell=True):
		return False
	return True

def chromWorker(argList):
	chrom, chromSize, options, sortedFiles, pid = argList
	# load chrom with pysam
	FA=pysam.FastaFile(options.reffile)
	# Create process pool before allocation
	p = mp.Pool(min(4,options.np), initializer=initWorker, initargs=(options, chrom))
	# allocate arrays
	meth = array.array('H',[0])*chromSize
	depth = array.array('H',[0])*chromSize
	# create coverage variable
	if options.rm_dup:
		coverage = array.array('B',[0])*chromSize
	else:
		# breaks if it can't be passsed
		coverage = []
	if options.CT_SNP:
		meth1 = array.array('H',[0])*chromSize
		depth1 = array.array('H',[0])*chromSize
	# Reference cache
	refCache = refcache(FA, chrom, chromSize)
	BS_conversion = {'+': ('C','T','G','A'), '-': ('G','A','C','T')}
	nmap = 0
	# read and filter input files
	for infile in sortedFiles:
		fileEXT = infile.split('.')[-1].upper()
		if fileEXT == 'BAM':
			# Read SAM with samtools
			disp("Reading %s from %s with samtools"%(chrom, infile))
			samflags = "-F 4"
			if options.pair: samflags += " -f 2"
			if options.unique: samflags += " -F 256"
			if bamIsSorted(infile):
				fin = sp.Popen("samtools view %s %s %s"%(samflags, infile, chrom), shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
			else:
				disp("%s should have been presorted\n"%(infile))
				return 0
			get_alignment = get_sam_alignment
		elif fileEXT == 'BSP':
			# Read BSP with awk filter
			disp("Reading %s from %s"%(chrom, infile))
			fin = sp.Popen("awk '$5 == \"%s\"' %s"%(chrom, infile), shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
			get_alignment = get_bsp_alignment
		else:
			disp("%s should be either presorted bsp or bam\n"%(infile))
			return 0
		for line in ifilter(lambda x: x[0] != "@", fin.stdout):
			map_info = get_alignment(line, options.unique, options.pair, \
				options.rm_dup, options.trim_fillin, coverage, chromSize)
			if not map_info: continue
			seq, strand, cr, pos = map_info
			pos2 = pos + len(seq)
			refseq = refCache.fetch(pos, pos2)
			assert(len(refseq) == len(seq))
			nmap += 1
			match, convert, rc_match, rc_convert = BS_conversion[strand]
			searchFunc(refseq, seq, depth, meth, convert, match, pos)
			if options.CT_SNP == 0: continue
			searchFunc(refseq, seq, depth1, meth1, rc_convert, rc_match, pos)
		fin.stdout.close()
		# Check for error messages
		STDERR = fin.stderr.read()
		if STDERR: disp("%s ERROR:\n%s"%(chrom, STDERR))
		fin.stderr.close()
	disp("%s used %i reads"%(chrom, nmap))
	if options.combine_CpG:
		disp('combining CpG methylation from both strands of %s...'%(chrom,))
		# TODO modify this to window across the chromosome instead of loading the whole thing
		for m in re.finditer('CG', FA.fetch(chrom)):
			pos = m.start()
			posp1 = pos+1
			depth[pos] += depth[posp1]
			meth[pos] += meth[posp1]
			depth[posp1] = 0
			meth[posp1] = 0
			if options.CT_SNP > 0:
				depth1[pos] += depth1[posp1]
				meth1[pos] += meth1[posp1]
				depth1[posp1] = 0
				meth1[posp1] = 0
	FA.close()
	# Block until turn to write
	if pid != 0:
		while syncArray[pid-1] == 0:
			time.sleep(5)
	# Open output files
	if not options.outfile:
		fout, outfile = sys.stdout, 'STDOUT'
	else:
		fout = open(options.outfile, 'a', 100000000)
		disp('writing %s to %s ...' % (chrom, options.outfile))
	if options.wigfile: 
		fwig = open(options.wigfile, 'a')
		fwig.write('variableStep chrom=%s span=%d\n' % (chrom, options.wigbin))
		bin = wigd = wigm = 0
	nc, nd, dep0 = 0, 0, options.min_depth
	# Write output
	indexedIter = izip(xrange(len(depth)), depth, meth, depth1, meth1)
	filteredIter = ifilter(lambda x: x[1] >= dep0, indexedIter)
	for ret in p.imap(calcMeth, filteredIter, chunksize=1000):
		if not ret: continue
		retStr, i, d, m = ret
		nc += 1
		nd += d
		if options.wigfile:
			if i / options.wigbin != bin:
				if wigd > 0: fwig.write('%d\t%.3f\n' % (bin*options.wigbin+1, min(wigm/wigd,1)))
				bin = i / options.wigbin
				wigd = wigm = 0.0
			wigd += d
			wigm += m
		fout.write(retStr)
	# Close and join pool
	p.close()
	p.join()
	# Close files
	if options.outfile != 'STDOUT': fout.close()
	if options.wigfile: fwig.close()
	# Allow next process to write
	syncArray[pid] = 1
	return (chrom, nmap, nc, nd)

def calcMeth(argList):
	global options
	global chrom
	i, dd, m, d1, m1 = argList
	if options.CT_SNP > 0:
		if m1 != d1:
			if options.CT_SNP == 2: return 0
			d = float(dd) * m1 / d1
		else: d = float(dd)
	else: d = float(dd)
	strand, seq = getContext(FA, chrom, i, options.full)
	if not strand: return 0
	if len(options.context) > 0:
		if seq not in options.context: return 0
	try: ratio = min(m, d) / d
	except ZeroDivisionError: return 0
	CIl, CIu = wilsonScore(ratio, d)
	if options.CT_SNP:
		vartup = (chrom, i+1, strand, seq, ratio, d, m, dd, m1, d1, CIl, CIu)
		retStr = '%s\t%d\t%c\t%s\t%.3f\t%.2f\t%d\t%d\t%d\t%d\t%.3f\t%.3f\n' % vartup
	else:
		vartup = (chrom, i+1, strand, seq, ratio, d, m, dd, CIl, CIu)
		retStr = '%s\t%d\t%c\t%s\t%.3f\t%.2f\t%d\t%d\tNA\tNA\t%.3f\t%.3f\n' % vartup
	return (retStr, i, d, m)

def initWorker(useOptions, useChrom):
	# Save full
	global options
	options = useOptions
	global chrom
	chrom = useChrom
	# Open FA reader
	global FA
	FA=pysam.FastaFile(options.reffile)
	
z95, z95sq = 1.96, 1.96 * 1.96
def wilsonScore(ratio, d):
	pmid = ratio + z95sq / (2.0 * d)
	sd = z95 * ((ratio*(1.0-ratio)/d + z95sq/(4.0*d*d)) ** 0.5)
	denorminator = 1.0 + z95sq / d
	CIl, CIu = (pmid - sd) / denorminator, (pmid + sd) / denorminator
	return CIl, CIu

def disp(txt, nt=0):
    if not quiet: sys.stderr.write('[methratio] @%s \t%s\n'%(time.asctime(), txt))

samFLAGS = 'pPuUrR12sfdS'
def parseFLAG(intFlag):
	'''
	Parse the samtools integer flag into a string format since
	this feature was removed.

	>>> parseFLAG(3)
	set(('p','P'))
	>>> parseFLAG(67)
	set(('p','P','1'))
	'''
	binFlag = bin(intFlag)[:1:-1]
	return set(compress(samFLAGS,map(int,binFlag)))

cigarRE = re.compile(r'\d+[a-zA-Z]')
def parseCigar(seq, cigar):
	'''
	Delete insertions and dash "-" deletions

	>>> parseCigar('ACTAGAATGGCT','3M1I3M1D5M')
	'ACTGAA-TGGCT'
	>>> parseCigar('ACTG','2M')
	Traceback (most recent call last):
	...
	AssertionError: String length does not match CIGAR
	'''
	index = 0 # index in seq
	originalLen = len(seq)
	cigarMatch = cigarRE.findall(cigar)
	for align in cigarMatch:
		length = int(align[:-1])
		op = align[-1]
		if op == 'M':
			index += length
		elif op == 'I':
			seq = seq[:index]+seq[index+length:]
		elif op == 'D':
			seq = seq[:index]+'-'*length+seq[index:]
			index += length
		else:
			raise ValueError("%c not a valid CIGAR character"%(op))
	assert originalLen == index, "String length does not match CIGAR"
	return seq

def get_sam_alignment(line, unique, pair, rm_dup, trim_fillin, coverage, chromLen):
	col = line.split('\t')
	cr, pos, cigar, seq, strand, insert = col[2], int(col[3])-1, col[5], col[9], '', int(col[8])
	strand_index = line.find('ZS:Z:')
	assert strand_index >= 0, 'missing strand information "ZS:Z:xx"'
	strand = line[strand_index+5:strand_index+7]
	seq = parseCigar(seq, cigar)
	if pos + len(seq) >= chromLen:
		#print("Seq is beyond %s"%(cr,))
		return ()
	if rm_dup:  # remove duplicate hits
	# Duplicates are matched based on the last base and strand of the fragment
		if strand in ('+-','-+'):
			frag_end, direction = pos+len(seq), 2
		else:
			frag_end, direction = pos, 1
		if coverage[frag_end] & direction:
			#print("Duplicate read (%i, %i)"%(pos, frag_end))
			return ()
		coverage[frag_end] |= direction # Record fragment
	if trim_fillin > 0: # trim fill in nucleotides
		if strand in ('+-','-+'):
			seq = seq[:-trim_fillin]
		elif strand in ('++','--'):
			seq, pos = seq[trim_fillin:], pos+trim_fillin
	if insert > 0:
		seq = seq[:int(col[7])-1-pos] # remove overlapped regions in paired hits, SAM format only
	return (seq, strand[0], cr, pos)

def get_bsp_alignment(line, unique, pair, rm_dup, trim_fillin, coverage, chromLen):
	col = line.split('\t')
	flag = col[3][:2]
	if flag == 'NM' or flag == 'QC': return ()
	if unique and flag != 'UM': return ()
	if pair and col[7] == '0': return ()
	seq, strand, cr, pos, insert, mm = col[1], col[6], col[4], int(col[5])-1, int(col[7]), col[9]
	#if cr not in chroms: return [] #shouldn't need this with grep
	if ':' in mm:
		tmp = mm.split(':')
		gap_pos, gap_size = int(tmp[1]), int(tmp[2])
		if gap_size < 0: seq = seq[:gap_pos] + seq[gap_pos-gap_size:]  # insertion on reference
		else: seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]
	if pos + len(seq) >= chromLen: return ()
	if rm_dup:  # remove duplicate hits
	# Duplicates are matched based on the last base and strand of the fragment
		if strand == '+-' or strand == '-+':
			frag_end, direction = pos+len(seq), 2
		else:
			frag_end, direction = pos, 1
		if coverage[frag_end] & direction:
			# This fragment has already been seen
			return ()
		coverage[frag_end] |= direction # Record fragment
	if trim_fillin > 0: # trim fill in nucleotides
		if strand == '+-' or strand == '-+': seq = seq[:-trim_fillin]
		elif strand == '++' or strand == '--': seq, pos = seq[trim_fillin:], pos+trim_fillin
	return (seq, strand[0], cr, pos)

def ParseFai(inFile, useChroms):
	'''
	Parses a fa.fai into a python dictionary
	Paramteters
	================================
	inFile	FILE	fai file
	'''
	chromDict = dict(map(lambda y: (y[0], int(y[1])), map(lambda y: y.split('\t'), open(inFile,'r').readlines())))
	if useChroms:
		return {k:v for k,v in chromDict.iteritems() if k in useChroms}
	else:
		return chromDict

def searchFunc(refseq, seq, depth, meth, convert, match, pos):
	for m in re.finditer(match, refseq):
		index = m.start()
		ip = index+pos
		if seq[index] == convert:
			depth[ip] += 1
		elif seq[index] == match:
			depth[ip] += 1
			meth[ip] += 1

revTab = maketrans('AGCT','TCGA')
def revcomp(seq):
	return seq[::-1].translate(revTab)

def getContext(FA, chrom, index, fine=False):
	# Will break for negative numbers, but not when past chrom
	try:
		refseq = FA.fetch(chrom, index-2, index+3)
		assert(len(refseq) == 5)
	except: return 0,0
	if refseq[2] == 'C':
		strand = '+'
		if refseq[3] == 'G': seq = 'CG'
		elif fine: seq = refseq[2:5]
		elif refseq[4] == 'G': seq = 'CHG'
		else: seq = 'CHH'
	else:
		strand = '-'
		if refseq[1] == 'C': seq = 'CG'
		elif fine: seq = revcomp(refseq[:3])
		elif refseq[0] == 'C': seq = 'CHG'
		else: seq = 'CHH'
	return strand, seq

def memAvail(p=0.8):
	'''
	Returns a fraction (p) of the free memory in megabytes.
	'''
	try:
		mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_AVPHYS_PAGES')  # e.g. 4015976448
	except:
		mem_bytes = int(sp.check_output("vm_stat | grep free | awk '{print $3}'", shell=True).rstrip('.\n'))*4096
	mem_mb = mem_bytes/(1024.**2)  # e.g. 3.74
	#mem_gib = mem_bytes/(1024.**3)  # e.g. 3.74
	return int(mem_mb*p)

if __name__ == "__main__":
	main()
