#!/usr/bin/python
from collections import defaultdict
import glob
import re
from math import log
import random
from sys import argv

cnt=defaultdict(lambda:defaultdict(lambda:0))
srrs=set()
spes=set()

tots=defaultdict(lambda:0);

g=glob.glob(argv[1])
for fn in g:
    try:
        srr=re.match('.*/SRR([0-9]*).results',fn).group(1)
    except:
        print 'ERROR: '+fn
        continue
    srrs.add(srr)
    for line in file(fn):
        try:
            (n,spe)=re.match(" *([0-9]*) (.*)\n",line).groups()
        except:
            print 'ERROR: '+line
        n=int(n)
        if n>1:
            cnt[spe][srr]=n
            spes.add(spe);
        tots[srr]+=n

srrs=list(srrs)
spes=list(spes)
srrs.sort()
spes.sort()

meta=file(argv[2])
sick=defaultdict(lambda:'No')
for line in meta:
    try:
        (srr,s)= re.match('[0-9x]* SRR([0-9]*) ([A-z]*)\n',line).groups();
        sick[srr]=s;
    except:
        print 'ERROR: '+line

o=file('table.csv','w')
o.write('species')
for srr in srrs:
    o.write(','+srr)
o.write('\nsick')
for srr in srrs:
    o.write(','+sick[srr])
for spe in spes:
    o.write('\n'+spe)
    for srr in srrs:
        o.write(','+str(cnt[spe][srr]))
o.close()
