#!/usr/bin/python
import argparse
parser = argparse.ArgumentParser(description='Process Filter compaered table INDEL by read deep and base INDEL.')
parser.add_argument('-i', '--in', dest='input', type=file, help='Full Input file name')
parser.add_argument('-R', '--minread', dest='minread', type=int, action='store', default=0, help='Number of minimun read deep [default=0]')
parser.add_argument('-D', '--minindel', dest='minindel', type=int, action='store', default=0, help='Number of minimun base INDEL [default=0]')
parser.add_argument('-o', '--output', dest='output', type=str, action='store', default='output', help='Base name for output file')
args = parser.parse_args()
m=args.input #'est_indeltable.tsv'
n=args.minread #'min. read deep'
o=args.minindel #'min base indel'
indel_table=[x.replace(' ', '').split('\t') for x in m.read().split('\n')[1:-1]]
_1=[x for x in indel_table if sum(int(x[5]), int(x[6]))!=0]
_2=[x for x in _1 if sum(int(x[7]), int(x[8]))!=0]
_3=_3=[x for x in _2 if len(x[4]) >= o]
print 'name\tstart\tend\tindel\tseq\tindel1\tnum.read1\tindel2\tnum.read2\n'
output=open(args.output+'.tsv', 'w')
for x in _3:
    print 'name\tstart\tend\tindel\tseq\tindel1\tnum.read1\tindel2\tnum.read2\n'
    if int(x[6]) > n and int(x[8])> n:
        if int(x[5]) !=0 and int(x[5]) > int(x[6])/2:
            print '\t'.join(x)
            output.wirte('\t'.join(x)+'\n')
        if int(x[7]) !=0 and int(x[7]) > int(x[8])/2:
            print '\t'.join(x)
            output.wirte('\t'.join(x)+'\n')
output.close()
print 'End of Program'
