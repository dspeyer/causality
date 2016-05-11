#!/usr/bin/python

from datareader import Data
from utils import *
from scipy.stats import chi2_contingency
from scipy.stats import chi2
from scipy.stats import chisquare
from math import log
from sys import stdout


data = Data()
data.read_from_files()

pl = data.get_pid_list(lambda(p): ('nod2' in p and 
                                   p.health in ['control', 'ileal CD'] and
                                   len(p.fractions)>0))

sick = data.get_data(pl, 'health', bucketizer=lambda(v):v=='ileal CD')
#vals = data.get_data(pl, 'fractions',  'Lactobacillus acidophilus')
#vals = data.get_data(pl, 'fractions',  'Bacteroides xylanolyticus')
#vals = data.get_data(pl, 'fractions',  'Dorea longicatena')
vals = data.get_data(pl, 'fractions',  'Eubacterium rectale')
nod2 = data.get_data(pl, 'nod2', bucketizer=lambda(v):v>0)

z=[0,0]
buckets=[[0]*20, [0]*20]
for s,v in zip(sick,vals):
    if v==0:
        z[s]+=1
    else:
        b=int(round(-3*(log(v)/log(10))))
        buckets[s][b]+=1

print '0, %d, %d' % (z[0], z[1])
for i in range(19,0,-1):
    print '%g, %d, %d' % (10**(-i/3.0), buckets[0][i], buckets[1][i])

co = findcutoff(vals, sick)
print co

vals2 = data.get_data(pl, 'fractions',  'Eubacterium rectale', bucketizer=lambda(v):(v>co.threshold))

print count(zip(sick, vals2, nod2))

cnt=count(zip(nod2, sick)) 
p_nod2 = float(sum(nod2))/len(nod2)
p_cd_given_nod2 = [ float(cnt[0][1]) / sum(cnt[0]),
                    float(cnt[1][1]) / sum(cnt[1])]

print 'p(nod2)=%s' % p_nod2
print 'p(cd|nod2)=%s' % p_cd_given_nod2

print direction(nod2, sick, vals2, len(pl), p_nod2, p_cd_given_nod2)

#(/ 24.0 58)
