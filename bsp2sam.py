#! /usr/bin/env python
import sys, optparse, time

usage = "usage: %prog [options] BSMAP_MAPPING_FILE"
parser = optparse.OptionParser(usage=usage)

parser.add_option("-o", "--out", dest="outfile", metavar="FILE", help="output file name. (required)", default="")
parser.add_option("-d", "--ref", dest="reffile", metavar="FILE", help="reference genome fasta file. (required)", default="")
parser.add_option("-q", "--quiet", action="store_true", dest="quiet", help="don't print progress on stderr.", default=False)

options, infile = parser.parse_args()
infile = infile[0]

def disp(txt, nt=0):
    if not options.quiet: print >> sys.stderr, ''.join(['\t' for i in xrange(nt)]+['@ ',time.asctime(),': ',txt])

assert any(options.reffile), "Missing reference file, must set -d/--ref."
assert any(options.outfile), "Missing output file, must set -o/--out."

fout = open(options.outfile, 'w')
disp('reading reference %s ...' % options.reffile)
fout.write('@HD\tVN:1.0\n')
cr, crlen = '', 0
for line in open(options.reffile):
    if line[0] == '>': 
        if any(cr): fout.write('@SQ\tSN:%s\tLN:%d\n' % (cr,crlen))
        cr, crlen = line[1:].split()[0], 0
    else: crlen += len(line) - 1

fout.write('@SQ\tSN:%s\tLN:%d\n@PG\tID:BSMAP_2.43\n' % (cr,crlen))

n = 0
for line in open(infile):
    col = line[:-1].split('\t')
    name, read, qual, flag = col[:4]
    n += 1
    if n % 10000000 == 0: disp('read %d lines' % n, nt=1)
    if flag == 'NM': fout.write('%s\tu\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n' % (name,read,qual))
    elif flag == 'QC': fout.write('%s\tuf\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n' % (name,read,qual))
    else:
        cr, pos, strand, ref, mm, samflag = col[4], col[5], col[6], col[8], col[9], ''
        if strand == '+-' or strand == '-+': samflag += 'r'
        if flag == 'MA' or flag == 'OF': samflag += 's'
        fout.write('%s\t%s\t%s\t%s\t255\t%dM\t*\t0\t0\t%s\t%s\tNM:i:%s\tZS:Z:%s\n' % (name,samflag,cr,pos,len(read),read,qual,mm,strand))

fout.close()
