#! /usr/bin/env python
import sys, time, os, array, optparse, math
usage = "usage: %prog [options] <GROUP1_SAMPLE1,GROUP1_SAMPLE2,...> <GROUP2_SAMPLE1,GROUP2_SAMPLE2,...>"
parser = optparse.OptionParser(usage=usage)

parser.add_option("-o", "--out", dest="outfile", metavar="FILE", help="output differential methylation file name. (required)", default="")
parser.add_option("-d", "--ref", dest="reffile", metavar="FILE", help="reference genome fasta file. (required)", default="")
parser.add_option("-b", "--bin", dest="bin", type="int", metavar='BIN', help="bin size. [default: 100]", default=100)
parser.add_option("-p", "--pval", dest="pval", type="float", metavar='PVAL', help="p-value cut-off. [default: 0.01]", default=0.01)
parser.add_option("-r", "--diff", dest="diff", type="float", metavar='DIFF', help="minimal abs meth ratio difference. [default: 0.33]", default=0.33)
parser.add_option("-x", "--context", dest="context", metavar='TYPE', help="methylation pattern type [CG|CHG|CHH], multiple types separated by ','. [default: all]", default='')
parser.add_option("-l", "--labels", dest="labels", metavar='LABELS', help="output label for each group, separated by ','. [default: filenames]", default='')
parser.add_option("-m", "--min-depth", dest="min_depth", type="int", metavar='FOLD', help="minimal average coverage. [default: 1]", default=1.)
parser.add_option("-s", "--strand", dest="strand", metavar='STRAND', help="which strand to process, [both|forward|reverse]. [default: both]", default='both')
parser.add_option("-q", "--quiet", action="store_true", dest="quiet", help="don't print progress on stderr.", default=False)

options, infiles = parser.parse_args()

assert options.reffile, "Require reference file name, use -d/--ref option."
assert options.outfile, "Require output file name, use -o/--out option."
assert len(infiles) >= 2, "Require at least two METHRATIO_FILE groups."
infiles = [x.split(',') for x in infiles] 
if len(options.context) > 0: options.context = options.context.split(',')
if len(options.labels) == 0: labels = ['group%d' % i for i in xrange(len(infiles))]
else: labels = options.labels.split(',')
assert len(labels) == len(infiles), "# labels not equal to # groups."
strand_dict = {'both': ('+-', 'CG'), 'forward': ('+', 'C'), 'reverse': ('-', 'G')}
assert options.strand in strand_dict, "Invalid strand: %s" % options.strand
strands, nts = strand_dict[options.strand]

PVAL = [1,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001,5e-4,2e-4,1e-4,5e-05,2e-05,1e-05,5e-06,1e-06,5e-07,1e-07,5e-08,1e-08,5e-09,1e-09,5e-10,1e-10,1e-11,\
  1e-12,1e-13,1e-14,1e-15,1e-16,1e-17,1e-18,1e-19,1e-20,1e-22,1e-24,1e-26,1e-28,1e-30,1e-32,1e-34,1e-36,1e-38,1e-40,1e-42,1e-44,1e-46,1e-48,1e-50,1e-55,\
  1e-60,1e-65,1e-70,1e-75,1e-80,1e-85,1e-90,1e-95,1e-100,1e-110,1e-120,1e-130,1e-140,1e-150,1e-160,1e-170,1e-180,1e-190,1e-200,1e-220,1e-240,1e-260,1e-280,1e-300,0]
ZVAL = [0.,1.28,1.65,1.96,2.33,2.58,2.81,3.09,3.29,3.48,3.72,3.89,4.06,4.26,4.42,4.56,4.89,5.03,5.33,5.45,5.73,5.85,6.11,6.22,6.47,6.81,7.13,7.44,7.74,\
  8.03,8.30,8.57,8.84,9.09,9.34,9.81,10.27,10.70,11.12,11.52,11.91,12.29,12.66,13.02,13.36,13.70,14.03,14.35,14.67,14.98,15.73,16.44,17.12,17.78,18.41,\
  19.03,19.62,20.20,20.76,21.31,22.36,23.36,24.33,25.25,26.15,27.01,27.85,28.67,29.46,30.23,31.71,33.13,34.49,35.80,37.07,40.]

zleft, zright = 0, len(PVAL)-1
while zleft < zright - 1:
    zmid = (zleft + zright) / 2
    if PVAL[zmid] < options.pval: zright = zmid
    else: zleft = zmid
z0 = ZVAL[zleft]

def disp(txt, nt=0):
    if not options.quiet: print >> sys.stderr, '[methdiff] @%s \t%s' %(time.asctime(), txt)

disp('reading reference file: %s ...' % options.reffile)
ref, cr = {}, ''
for line in open(options.reffile):
    if line[0] == '>': 
        if cr != '': ref[cr] = refcr.upper()
        cr, refcr = line[1:].strip(), ''
    else: refcr += line[:-1]

ref[cr] = refcr.upper()

def get_chrom(group, sample, cr):
    line, fin = lines[group][sample], fins[group][sample]
    mm, dd = meth[group], depth[group]
    while len(line) > 0:  
        col = line.split('\t', 7)
        if col[2] not in strands: 
            line = fin.readline()
            continue
        if col[0] != cr: break
        pos = int(col[1]) - 1
        if len(options.context) > 0:
            if col[3] not in options.context:
                line = fin.readline() 
                continue
	index = int(pos/options.bin)
        mm[index] += float(col[6])
        dd[index] += float(col[5])
        line = fin.readline()
    lines[group][sample] = line
    return

def conf_intv(m, d, z):    
    p, z2 = min(m,d) / d, z * z   
    pmid = p + z2 / (2.0 * d)
    span = z * (p * (1.0 - p) / d + z2 / (4.0 * d * d)) ** 0.5
    denorm = 1.0 + z2 / d
    return ((pmid-span)/denorm, (pmid+span)/denorm)

def get_pval(m0, d0, m1, d1):
    l0, u0 = conf_intv(m0, d0, z0)
    l1, u1 = conf_intv(m1, d1, z0)
    if l0 < u1 and l1 < u0: return 1
    left, right = zleft, len(ZVAL)-1
    while right - left > 1:
        mid = (left + right) / 2
        zm = ZVAL[mid]
        l0, u0 = conf_intv(m0, d0, zm)
        l1, u1 = conf_intv(m1, d1, zm)
        if l0 >= u1 or l1 >= u0: left = mid
        else: right = mid
    return PVAL[left]
    
def cmp_chrom(cr):
    m, d = [0. for i in infiles], [0. for i in infiles]
    refcr = ref[cr]
    for pos in xrange(len(refcr)/options.bin+1):    
        for i in xrange(len(infiles)):
            m[i] = meth[i][pos]
            d[i] = depth[i][pos]
        if d[0] * d[1] == 0.: continue
        r0, r1 = m[0]/d[0], m[1]/d[1]
        if abs(r0-r1) < options.diff: continue
        start, end = pos*options.bin, (pos+1)*options.bin
        nc = sum(refcr.count(nt, start, end) for nt in nts)
        if d[0] / nc < options.min_depth or d[1] / nc < options.min_depth: continue
        pval = get_pval(m[0], d[0], m[1], d[1])
        if pval > options.pval: continue
        fout.write('%s\t%d\t%d\t%.1g\t%.3f\t%.3f\t%.1f\t%.1f\t%.3f\t%.1f\t%.1f\n' % (cr,start+1,end,pval,r0-r1,r0,d[0],m[0],r1,d[1],m[1]))

def xopen(filename):
    if filename[-3:].lower() == '.gz': return os.popen('zcat %s' % filename)
    else: return open(filename)

meth = [0 for infile in infiles]
depth = [0 for infile in infiles]
fins = [[xopen(infile) for infile in group] for group in infiles]
lines = [[fin.readline() for fin in group] for group in fins]
fout = open(options.outfile, 'w')
name0, name1 = labels[0], labels[1]
fout.write('cr\tstart\tend\tp_value\tdiff_ratio\t%s_ratio\t%s_depth\t%s_meth\t%s_ratio\t%s_depth\t%s_meth\n' % (name0,name0,name0,name1,name1,name1))
for cr in sorted(ref.keys()):
    disp('processing %s ...' % cr)
    for g, group in enumerate(infiles): 
        n_bins = int(math.ceil(len(ref[cr])/float(options.bin)))
        meth[g] = array.array('f', [0.]) * n_bins
        depth[g] = array.array('f', [0.]) * n_bins
        for i in xrange(len(group)): get_chrom(g, i, cr)
    cmp_chrom(cr)

fout.close()
disp('done.')
