#!/usr/bin/python
import progressbar
from progressbar import Bar, Percentage, ProgressBar
import argparse
parser = argparse.ArgumentParser(description='Process Filter INDEL by read deep and base INDEL.')
parser.add_argument('-i', '--in', dest='input', type=str, help='INDEL sort file from SAM as Input file')
parser.add_argument('-S', '--sam', dest='sam', type=str, action='store', default=0, help='SAM file')
parser.add_argument('-o', '--output', dest='output', type=str, action='store', default='output', help='Base name for output file')
args = parser.parse_args()
sam=args.sam
indelfile=args.input
_input1=open(sam, 'r')
tem =[x.split('\t') for x in _input1.read().split('\n')[:-1]]
samhead=[x for x in tem  if x[0] in ['@HD', '@SQ', '@RG', '@PG', '@CO']]
samdat=[x for x in tem  if x[0] not in ['@HD', '@SQ', '@RG', '@PG', '@CO']]
tigar=[x for x in samdat if x[5] != '*']
reflist=[x[1][3:] for x in samhead if '@SQ' in x]
cigar=[]
pbar1 = ProgressBar(widgets=['Load SAM file  : ', Percentage(), Bar()], maxval=len(tigar)).start()
i=0
for x in tigar:
    y=[x[0], x[2], x[3], x[5], x[9]]
    cigar.append(y)
    pbar1.update(i+1)
    i=i+1
pbar1.finish()
_input2=open(indelfile, 'r')
indel=[x.split('\t') for x in _input2.read().split('\n')[:-1]]
selindel=[x for x in indel if int(x[5]) >= 10]
countindel=[]
for x in selindel:
    j=0
    for y in cigar:
        if x[0]==y[1] and int(y[2])< int(x[1]) < int(y[2])+len(y[4]):
            j=j+1
    z=[x[0], x[1], x[2], x[3], x[4], x[5], j]
    countindel.append(z)
out=open(args.output+'_indel_count.tsv', 'w')
for x in countindel:
    z=str(x).replace(' ', '').replace('[', '').replace(']', '').replace(',', '\t').replace('\'', '')+'\n'
    out.write(z)
out.close()
print 'Fnish'