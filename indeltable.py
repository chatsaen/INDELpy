#!/usr/bin/python
import progressbar
from progressbar import Bar, Percentage, ProgressBar
import argparse
parser = argparse.ArgumentParser(description='Process Filter compaered table INDEL by read deep and base INDEL.')
parser.add_argument('-i1', '--in1', dest='input1', type=file, help='Full Input file name')
parser.add_argument('-i2', '--in2', dest='input2', type=int, action='store', default=0, help='Number of minimun read deep [default=0]')
parser.add_argument('-s1', '--sam1', dest='sam1', type=int, action='store', default=0, help='Number of minimun base INDEL [default=0]')
parser.add_argument('-s2', '--sam2', dest='sam2', type=int, action='store', default=0, help='Number of minimun base INDEL [default=0]')
parser.add_argument('-o', '--output', dest='output', type=str, action='store', default='output', help='Base name for output file')
args = parser.parse_args()
m=args.input1#'estdura.sam_indel_count.tsv'
n=args.input2#'esttenera.sam_indel_count.tsv'
o=args.sam1#'estdura.sam'
p=args.sam2#'esttenera.sam'
sam1=open(o, 'r')
tem =[x.split('\t') for x in sam1.read().split('\n')[:-1]]
samhead=[x for x in tem  if x[0] in ['@HD', '@SQ', '@RG', '@PG', '@CO']]
samdat=[x for x in tem  if x[0] not in ['@HD', '@SQ', '@RG', '@PG', '@CO']]
tigar=[x for x in samdat if x[5] != '*']
cigar_1=[]
for x in tigar:
    y=[x[0], x[2], x[3], x[5], x[9]]
    cigar_1.append(y)
sam2=open(p, 'r')
tem =[x.split('\t') for x in sam2.read().split('\n')[:-1]]
samhead=[x for x in tem  if x[0] in ['@HD', '@SQ', '@RG', '@PG', '@CO']]
samdat=[x for x in tem  if x[0] not in ['@HD', '@SQ', '@RG', '@PG', '@CO']]
tigar=[x for x in samdat if x[5] != '*']
cigar_2=[]
for x in tigar:
    y=[x[0], x[2], x[3], x[5], x[9]]
    cigar_2.append(y)
indelcount1=open(m, 'r')
indel_1=[x.replace(' ', '').split('\t') for x in indelcount1.read().split('\n')[:-1]]
indelcount2=open(n, 'r')
indel_2=[x.replace(' ', '').split('\t') for x in indelcount2.read().split('\n')[:-1]]
tem=[x[0] for x in indel_1]+[x[0] for x in indel_2]
refname=list(set(tem))
refname=refname.sort()
def uniqset(__iter__):
    check=[]
    res=[]
    for x in __iter__:
        if x not in check:
            res.append(x)
            check.append(x)
    return res
tem3=[]
for x in refname:
    tem2=uniqset([[int(y[1]), int(y[2])] for y in indel_1 if x==y[0]]+[[int(y[1]), int(y[2])] for y in indel_2 if x==y[0]])
    tem2.sort()
    for z in tem2:
        _1=[]
        for a in indel_1:
            if x==a[0] and [int(a[1]), int(a[2])] == z:
                _1=a
        _2=[]
        for a in indel_2:
            if x==a[0] and [int(a[1]), int(a[2])] == z:
                _2=a
        if len(_1)==0:
            j=0
            for b in cigar_1:
                if x==b[1] and int(b[2])< int(z[0]) < int(b[2])+len(b[4]):
                    j=j+1
            c=[x, z[0], z[1], _2[3], _2[4], 0, j]
            _1=c
        if len(_2)==0:
            j=0
            for b in cigar_2:
                if x==b[1] and int(b[2])< int(z[0]) < int(b[2])+len(b[4]):
                    j=j+1
            c=[x, z[0], z[1], _1[3], _1[4], 0, j]
            _2=c
        d=[x, z[0], z[1], _1[3], _1[4], int(_1[5]), int(_1[6]), int(_2[5]), int(_2[6])]
        tem3.append(d)
out=open(args.output+'indeltable.tsv', 'w')
out.write('name\tstart\tend\tindel\tseq\tindel1\tnum.read1\tindel2\tnum.read2\n')
for x in tem3:
    y='\t'.join(map(str, x)).replace(' ', '')
    out.write(y+'\n')
out.close()