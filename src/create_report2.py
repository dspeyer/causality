#!/usr/bin/python

from datareader import Data
from utils import *
from scipy.stats import chi2_contingency
from scipy.stats import chi2
from scipy.stats import chisquare
from math import log
from sys import stdout, argv
import direction, severs
from common_cause import has_common_cause


data = Data()
data.read_from_files()

pl = data.get_pid_list(lambda(p): ('nod2' in p and 
                                   p.health in ['control', 'ileal CD'] and
                                   len(p.fractions)>0))

sick = data.get_data(pl, 'health', bucketizer=lambda(v):v=='ileal CD')
nod2 = data.get_data(pl, 'nod2', bucketizer=lambda(v):v>0)

print '%d patients with both gene and microbiome data' % len(pl)

print 'Total bacterial species: %d' % len(data.bacteria)
interesting = {}
cos = {}
bfs = []
for species in data.bacteria:
    vals = data.get_data(pl, 'fractions', species)
    co = findcutoff(vals, sick)
    if co.sick_when_more==None:
        continue
    boolvals = [(i>co.threshold)==(co.sick_when_more) for i in vals]
    if entropy(boolvals)<.5:
        continue
    interesting[species]=boolvals
    cos[species]=co
    cnt = count(zip(boolvals, nod2))
    sumarr(cnt, 0.1)
    p = chi2_contingency(cnt)[1]
    if p < .01:
        continue
    cnt = count(zip(boolvals, sick))
    sumarr(cnt, 0.1)
    p = chi2_contingency(cnt)[1]
    if p > .01:
        continue
    bf = direction.montecarlo(nod2, sick, boolvals, len(nod2)).bayes_fwd_rev
    bfs.append(bf)
    if bf > 2 and '-print-gene-test' in argv:
        print '%s & %s & %.2g & %.2f' % (species, prettyco(co), p, bf)

print

if '-histogram' in argv:
    maxv = log(max(bfs))
    minv = log(min(bfs))
    n = len(bfs)
    half_n_buckets = int(n/10)
    bs = max(maxv, -minv) / half_n_buckets
    buckets = [0] * (2*half_n_buckets+1)
    for i in bfs:
        buckets[ int((log(i)+.5)/bs)+half_n_buckets ] += 1
    for i,v in enumerate(buckets):
#        print '%.2f, %.2f, %d' % (exp((i-1-half_n_buckets)*bs), exp((i-half_n_buckets)*bs), v)
        print '%.2f, ' % (exp((i-1-half_n_buckets)*bs))
        print '%.2f, %d' % (exp((i-.5-half_n_buckets)*bs), v)

if '-trio' in argv:
    for s1 in interesting:
        v1 = interesting[s1]
        for s2 in interesting:
            v2 = interesting[s2]
            for s3 in interesting:
                v3 = interesting[s3]
                if not has_common_cause([v1,v2,v3],0.01):
                    continue
                if link(v1,sick)<.01 and link(v2,sick)<.01:
                    sev=severs.mle(v2,sick,v1)
                    if sev>2:
                        print '%s & %s & %s & %s & %.2f \\\\' % (s1, prettyco(cos[s1]), s2, s3, sev)

print
                    
print '%d interesting species' % len(interesting)
