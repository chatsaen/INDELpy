#!/usr/bin/python
#
#
#python exp.py 'contigreadlist.txt '  'output base name'
#
#
import progressbar
from progressbar import Bar, Percentage, ProgressBar
import argparse
parser = argparse.ArgumentParser(description='Process Sorting  INDEL from SAM file by CIGAR code.')
parser.add_argument('-i', '--in', dest='input', type=file, help='Full Input SAM file name')
args = parser.parse_args()
m=args.input #'est_indeltable.tsv'
_input=open(m, 'r')
tem =[x.split('\t') for x in _input.read().split('\n')[:-1]]
samhead=[x for x in tem  if x[0] in ['@HD', '@SQ', '@RG', '@PG', '@CO']]
samdat=[x for x in tem  if x[0] not in ['@HD', '@SQ', '@RG', '@PG', '@CO']]
tigar=[x for x in samdat if x[5] != '*']
#onmap=[x for x in samdat if x[5] != '*']
#nomap=[x for x in samdat if x[5] != '*']
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
print 'cigar[read, ref, pos, cigar, seq]'
tem3=[]
print "tem3[read, ref, tem2, seq]"
for xx in cigar:
    a=int(xx[2])
    b=[[int(x.split(' ')[0]), x.split(' ')[1]] for x in xx[3].replace('M', ' M,').replace('I', ' I,').replace('D', ' D,').replace('N', ' N,').replace('S', ' S,').replace('H', ' H,').replace('P', ' P,').replace('=', ' =,').replace('X', ' X,').split(',')[:-1]]
    i=0
    tem2=[]
    #tem2[[pos start, pos end], num seq , aling type, cut seq], ...]
    for z in b:
        k=[]
        if z[1] in ['M', 'I', 'X', '=']:
            y=xx[4][i:i+z[0]]
            i=i+z[0]
            if z[1] != 'M':
                k=[a, a]
            else:
                k=[a, a+z[0]]
                a=a+z[0]
        else:
            k=[a, a+z[0]]
            a=a+z[0]
            y='-'*z[0]
        tem2.append([k, z[1], z[0], y])
    tem3.append([xx[0], xx[1], tem2, xx[4]])
indel=[]
pbar2 = ProgressBar(widgets=['Sorting indel : ', Percentage(), Bar()], maxval=len(reflist)).start()
i=0
for x in reflist:
    checkI=[]
    temI=[]
    checkD=[]
    temD=[]
    checkX=[]
    temX=[]
    for y in tem3:
        if x==y[1]:
            for z in y[2]:
                if z[1] =='I' :
                    zz=[z[0], z[3]]
                    if zz not in checkI:
                        checkI.append(zz)
                    temI.append([z[0],z[3]])
                if z[1] =='D' :
                    zz=[z[0], z[3]]
                    if zz not in checkD:
                        checkD.append(zz)
                    temD.append([z[0],z[3]])
                if z[1] =='X' :
                    zz=[z[0], z[3]]
                    if zz not in checkX:
                        checkX.append(zz)
                    temX.append([z[0],z[3]])
    tem5=[]
    for x2 in checkI:
        nI=temI.count(x2)
        tem4=[x, x2[0][0], x2[0][1], 'I', x2[1], nI]
        tem5.append(tem4)
    for x2 in checkD:
        nD=temD.count(x2)
        tem4=[x, x2[0][0], x2[0][1], 'D', x2[1], nD]
        tem5.append(tem4)
    indel.append(tem5)
    pbar2.update(i+1)
    i=i+1
pbar2.finish()
print 'output: '+m.split('.')[0]+'_indel.tsv'
out=open(m.split('.')[0]+'_indel.tsv', 'w')
for x in indel:
    for y in x:
        z=str(y).replace('[', '').replace(']', '').replace(',', '\t').replace('\'', '')+'\n'
        out.write(z)
out.close()
countindel=[]
_input=open(m+'_indel.tsv', 'r')
indel=[x.split('\t') for x in _input.read().split('\n')[:-1]]
selindel=[x for x in indel if int(x[5]) >= 10]
countindel=[]
pbar3 = ProgressBar(widgets=['Count read in indel position  : ', Percentage(), Bar()], maxval=len(selindel)).start()
i=0
for x in selindel:
    j=0
    for y in cigar:
        if x[0]==y[1] and int(y[2])< int(x[1]) < int(y[2])+len(y[4]):
            j=j+1
    z=[x[0], x[1], x[2], x[3], x[4], x[5], j]
    countindel.append(z)
    pbar3.update(i+1)
    i=i+1
pbar3.finish()
out=open(m+'_indel_count.tsv', 'w')
for x in countindel:
    z=str(x).replace('[', '').replace(']', '').replace(',', '\t').replace('\'', '')+'\n'
    out.write(z)
out.close()
print 'Fnish'

_input=open('esttenera.sam_indel.tsv', 'r')
indel=[x.split('\t') for x in _input.read().split('\n')[:-1]]
selindel=[x for x in indel if int(x[5]) >= 10]
countindel=[]
for x in selindel:
    j=0
    for y in cigar:
        if x[0]==y[1] and int(y[2])< int(x[1]) < int(y[2])+len(y[4]):
            j=j+1
    z=[x[0], x[1], x[2], x[3], x[4], x[5], j]
    countindel.append(z)