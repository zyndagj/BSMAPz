#! /usr/bin/env python
import sys, time, os, array, optparse
from itertools import compress
usage = "usage: %prog [options] BSMAP_MAPPING_FILES"
parser = optparse.OptionParser(usage=usage)

parser.add_option("-o", "--out", dest="outfile", metavar="FILE", help="output methylation ratio file name. [default: STDOUT]", default="")
parser.add_option("-O", "--alignment-copy", dest="alignfile", metavar="FILE", help="save a copy of input alignment for BSMAP pipe input. (in BAM format) [default: none]", default="")
parser.add_option("-w", "--wig", dest="wigfile", metavar="FILE", help="output methylation ratio wiggle file. [default: none]", default="")
parser.add_option("-b", "--wig-bin", dest="wigbin", type="int", metavar='BIN', help="wiggle file bin size. [default: 25]", default=25)
parser.add_option("-d", "--ref", dest="reffile", metavar="FILE", help="reference genome fasta file. (required)", default="")
parser.add_option("-c", "--chr", dest="chroms", metavar="CHR", help="process only specified chromosomes, separated by ','. [default: all]\nexample: --chroms=chr1,chr2", default=[])
parser.add_option("-s", "--sam-path", dest="sam_path", metavar="PATH", help="path to samtools. [default: none]", default='')
parser.add_option("-u", "--unique", action="store_true", dest="unique", help="process only unique mappings/pairs.", default=False)                                    
parser.add_option("-p", "--pair", action="store_true", dest="pair", help="process only properly paired mappings.", default=False)
parser.add_option("-z", "--zero-meth", action="store_true", dest="meth0", help="report loci with zero methylation ratios. (depreciated, -z will be always enabled)", default=True)
parser.add_option("-q", "--quiet", action="store_true", dest="quiet", help="don't print progress on stderr.", default=False)
parser.add_option("-r", "--remove-duplicate", action="store_true", dest="rm_dup", help="remove duplicated reads.", default=False)
parser.add_option("-t", "--trim-fillin", dest="trim_fillin", type="int", metavar='N', help="trim N end-repairing fill-in nucleotides. [default: 0]", default=0)
parser.add_option("-g", "--combine-CpG", action="store_true", dest="combine_CpG", help="combine CpG methylaion ratios on both strands.", default=False)
parser.add_option("-m", "--min-depth", dest="min_depth", type="int", metavar='FOLD', help="report loci with sequencing depth>=FOLD. [default: 1]", default=1)
parser.add_option("-n", "--no-header", action="store_true", dest="no_header", help="don't print a header line", default=False)
parser.add_option("-i", "--ct-snp", dest="CT_SNP", help='how to handle CT SNP ("no-action", "correct", "skip"), default: "correct".', default="correct")
parser.add_option("-x", "--context", dest="context", metavar='TYPE', help="methylation pattern type [CG|CHG|CHH], multiple types separated by ','. [default: all]", default='')

options, infiles = parser.parse_args()

if len(options.reffile) == 0: parser.error("Missing reference file, use -d or --ref option.")
if len(infiles) == 0: parser.error("Require at least one BSMAP_MAPPING_FILE.") 
if len(options.chroms) > 0: options.chroms = options.chroms.split(',')
CT_SNP_val = {"no-action": 0, "correct": 1, "skip": 2}
try: options.CT_SNP = CT_SNP_val[options.CT_SNP.lower()]
except: parser.error('Invalid -i value, select "no-action", "correct" or "skip"')
if options.min_depth <= 0: parser.error('Invalid -m value, must >= 1')
if options.trim_fillin < 0: parser.error('Invalid -t value, must >= 0')
if len(options.context) > 0: options.context = options.context.split(',')

if len(options.sam_path) > 0: 
    if options.sam_path[-1] != '/': options.sam_path += '/'

def disp(txt, nt=0):
    if not options.quiet: print >> sys.stderr, '[methratio] @%s \t%s' %(time.asctime(), txt)

samFLAGS = 'pPuUrR12sfdS'
def parseFLAG(intFlag):
	'''
	Parse the samtools integer flag into a string format since
	this feature was removed.
	'''
	binFlag = bin(intFlag)[:1:-1]
	return set(compress(flags,map(int,binFlag)))

if len(options.outfile) == 0: disp("Missing output file name, write to STDOUT.")
def get_alignment(line):
    col = line.split('\t')
    if sam_format:
        if line[0] == '@': return []
        flag = parseFLAG(int(col[1]))
        if 'u' in flag: return []
        if options.unique and 's' in flag: return []
        if options.pair and 'P' not in flag: return []
        cr, pos, cigar, seq, strand, insert = col[2], int(col[3])-1, col[5], col[9], '', int(col[8])
        if cr not in options.chroms: return []
        strand_index = line.find('ZS:Z:')
        assert strand_index >= 0, 'missing strand information "ZS:Z:xx"'
        strand = line[strand_index+5:strand_index+7]
        gap_pos, gap_size = 0, 0
        while 'I' in cigar or 'D' in cigar:
            for sep in 'MID':
                try: gap_size = int(cigar.split(sep, 1)[0])
                except ValueError: continue
                break
            if sep == 'M': gap_pos += gap_size
            elif sep == 'I': seq = seq[:gap_pos] + seq[gap_pos+gap_size:]
            elif sep == 'D': 
                seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]                        
                gap_pos += gap_size
            cigar = cigar[cigar.index(sep)+1:]
    else:
        flag = col[3][:2]
        if flag == 'NM' or flag == 'QC': return []
        if options.unique and flag != 'UM': return []
        if options.pair and col[7] == '0': return []
        seq, strand, cr, pos, insert, mm = col[1], col[6], col[4], int(col[5])-1, int(col[7]), col[9]
        if cr not in options.chroms: return []
        if ':' in mm:
            tmp = mm.split(':')
            gap_pos, gap_size = int(tmp[1]), int(tmp[2])
            if gap_size < 0: seq = seq[:gap_pos] + seq[gap_pos-gap_size:]  # insertion on reference
            else: seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]
    if pos + len(seq) >= len(ref[cr]): return []
    if options.rm_dup:  # remove duplicate hits
        if strand == '+-' or strand == '-+': frag_end, direction = pos+len(seq), 2
        else: frag_end, direction = pos, 1
        if coverage[cr][frag_end] & direction: return []
        coverage[cr][frag_end] |= direction
    if options.trim_fillin > 0: # trim fill in nucleotides
        if strand == '+-' or strand == '-+': seq = seq[:-options.trim_fillin]
        elif strand == '++' or strand == '--': seq, pos = seq[options.trim_fillin:], pos+options.trim_fillin
    if sam_format and insert > 0: seq = seq[:int(col[7])-1-pos] # remove overlapped regions in paired hits, SAM format only
    return (seq, strand[0], cr, pos)

# open pipes to alignment files
pipes = []
for infile in infiles:
    nline = 0
    if infile.strip() == '-': sam_format, fin, infile = True, os.popen('%ssamtools view -Sh -' % options.sam_path), 'STDIN'
    elif infile[-4:].upper() == '.SAM': sam_format, fin = True, os.popen('%ssamtools view -S %s' % (options.sam_path, infile)) 
    elif infile[-4:].upper() == '.BAM': sam_format, fin = True, os.popen('%ssamtools view %s' % (options.sam_path, infile))
    else: sam_format, fin = False, open(infile)
    pipes.append((sam_format,fin))

# Read in chromosomes
ref, cr, seq = {}, '', ''
disp('loading reference file: %s ...' % options.reffile)
for line in open(options.reffile):
    if line[0] == '>': 
        if len(cr) > 0: 
            if len(options.chroms) == 0 or cr in options.chroms: ref[cr] = seq.upper()
        cr, seq = line[1:-1].split()[0], ''
    else: seq += line.strip()

if len(options.chroms) == 0 or cr in options.chroms: ref[cr] = seq.upper()
del seq

meth, depth, coverage, meth1, depth1 = {}, {}, {}, {}, {}
for cr in ref:
    meth[cr] = array.array('H', [0]) * len(ref[cr])
    depth[cr] = array.array('H', [0]) * len(ref[cr])
    if options.rm_dup: coverage[cr] = array.array('B', [0]) * len(ref[cr])
    if options.CT_SNP > 0:
        meth1[cr] = array.array('H', [0]) * len(ref[cr])
        depth1[cr] = array.array('H', [0]) * len(ref[cr])

options.chroms = set(ref.keys())
BS_conversion = {'+': ('C','T','G','A'), '-': ('G','A','C','T')}
nmap = 0
for sam_format, fin in pipes:
    nline = 0
    if len(options.alignfile) > 0: pout = os.popen('%ssamtools view -bS - > %s' % (options.sam_path, options.alignfile), 'w')
    for line in fin:
        if len(options.alignfile) > 0: pout.write(line)
        nline += 1
        if nline % 10000000 == 0: disp('read %d lines' % nline, nt=1)
        map_info = get_alignment(line)
        if len(map_info) == 0: continue
        seq, strand, cr, pos = map_info
        depthcr = depth[cr]
        pos2 = pos + len(seq)
        nmap += 1
        methcr = meth[cr]
        refseq = ref[cr]
        match, convert, rc_match, rc_convert = BS_conversion[strand]
        index = refseq.find(match, pos, pos2)
        while index >= 0:
            if seq[index-pos] == convert: 
                try: depthcr[index] += 1
                except OverflowError: depthcr[index] = 65535
            elif seq[index-pos] == match and depthcr[index] < 65535:
                depthcr[index] += 1
                methcr[index] += 1
            index = refseq.find(match, index+1, pos2)
        if options.CT_SNP == 0: continue
        methcr1 = meth1[cr]
        depthcr1 = depth1[cr]
        index = refseq.find(rc_match, pos, pos2)
        while index >= 0:
            if seq[index-pos] == rc_convert: 
                try: depthcr1[index] += 1
                except OverflowError: depthcr1[index] = 65535
            elif seq[index-pos] == rc_match and depthcr1[index] < 65535: 
                depthcr1[index] += 1
                methcr1[index] += 1
            index = refseq.find(rc_match, index+1, pos2)
    
    fin.close()
    if len(options.alignfile) > 0: pout.close()

disp('read %d lines' % nline, nt=1)
if options.combine_CpG:
    disp('combining CpG methylation from both strands ...')
    for cr in depth:
        depthcr, methcr, refcr = depth[cr], meth[cr], ref[cr]
        if options.CT_SNP > 0: depthcr1, methcr1 = depth1[cr], meth1[cr]
        pos = refcr.find('CG')
        while pos >= 0:
            try: 
                depthcr[pos] += depthcr[pos+1]
                methcr[pos] += methcr[pos+1]
            except OverflowError: 
                depthcr[pos] = (depthcr[pos] + depthcr[pos+1]) / 2
                methcr[pos] = (methcr[pos] + methcr[pos+1]) / 2
            depthcr[pos+1] = 0
            methcr[pos+1] = 0
            if options.CT_SNP > 0:
                try: 
                    depthcr1[pos] += depthcr1[pos+1]
                    methcr1[pos] += methcr1[pos+1]
                except OverflowError: 
                    depthcr1[pos] = (depthcr1[pos] + depthcr1[pos+1]) / 2
                    methcr1[pos] = (methcr1[pos] + methcr1[pos+1]) / 2
            pos = refcr.find('CG', pos+2)

if len(options.outfile) == 0: fout, outfile = sys.stdout, 'STDOUT'
else: fout = open(options.outfile, 'w')
disp('writing %s ...' % options.outfile)
if options.wigfile: 
    fwig = open(options.wigfile, 'w')
    fwig.write('track type=wiggle_0\n')
if not options.no_header: 
    fout.write('chr\tpos\tstrand\tcontext\tratio\teff_CT_count\tC_count\tCT_count\trev_G_count\trev_GA_count\tCI_lower\tCI_upper\n')
z95, z95sq = 1.96, 1.96 * 1.96
nc, nd, dep0 = 0, 0, options.min_depth
for cr in sorted(depth.keys()):
    depthcr, methcr, refcr = depth[cr], meth[cr], ref[cr]
    if options.CT_SNP > 0: depthcr1, methcr1 = depth1[cr], meth1[cr]
    if options.wigfile:
        fwig.write('variableStep chrom=%s span=%d\n' % (cr, options.wigbin))
        bin = wigd = wigm = 0
    for i, dd in enumerate(depthcr):
        if dd < dep0: continue
        if options.CT_SNP > 0: 
            m1, d1 = methcr1[i], depthcr1[i]
            if m1 != d1:
                if options.CT_SNP == 2: continue
                d = float(dd) * m1 / d1
            else: d = float(dd)
        else: d = float(dd)
        if refcr[i] == 'C':
            strand = '+'
            try:
                if refcr[i+1] == 'G': seq = 'CG'
                elif refcr[i+2] == 'G': seq = 'CHG'
                else: seq = 'CHH'
            except IndexError: continue
        else:
            strand = '-'
            if i == 0: continue
            if refcr[i-1] == 'C': seq = 'CG'
            else:
                if i == 1: continue
                if refcr[i-2] == 'C': seq = 'CHG'
                else: seq = 'CHH'
        if len(options.context) > 0:
            if seq not in options.context: continue
        m = methcr[i]
        try: ratio = min(m, d) / d
        except ZeroDivisionError: continue
        nc += 1
        nd += d
        if options.wigfile:
            if i / options.wigbin != bin:
                if wigd > 0: fwig.write('%d\t%.3f\n' % (bin*options.wigbin+1, min(wigm/wigd,1)))
                bin = i / options.wigbin
                wigd = wigm = 0.  
            wigd += d
            wigm += m
        pmid = ratio + z95sq / (2 * d)
        sd = z95 * ((ratio*(1-ratio)/d + z95sq/(4*d*d)) ** 0.5)
        denorminator = 1 + z95sq / d
        CIl, CIu = (pmid - sd) / denorminator, (pmid + sd) / denorminator
        if options.CT_SNP: fout.write('%s\t%d\t%c\t%s\t%.3f\t%.2f\t%d\t%d\t%d\t%d\t%.3f\t%.3f\n' % (cr, i+1, strand, seq, ratio, d, m, dd, m1, d1, CIl, CIu))
        else: fout.write('%s\t%d\t%c\t%s\t%.3f\t%.2f\t%d\t%d\tNA\tNA\t%.3f\t%.3f\n' % (cr, i+1, strand, seq, ratio, d, m, dd, CIl, CIu))
                            
if options.outfile != 'STDOUT': fout.close()
if options.wigfile: fwig.close()
disp('total %d valid mappings, %d covered cytosines, average coverage: %.2f fold.' % (nmap, nc, float(nd)/nc))
