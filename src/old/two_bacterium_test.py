#!/usr/bin/python

from datareader import Data
from utils import *
from scipy.stats import chi2_contingency
from scipy.stats import chi2
from scipy.stats import chisquare
from math import log
from sys import stdout
import re
import direction, severs
from common_cause import has_common_cause


data = Data()
data.read_from_files()

pl = data.get_pid_list(lambda(p): (p.health in ['control', 'ileal CD'] and
                                   len(p.fractions)>0))
n = len(pl)

sick = data.get_data(pl, 'health', bucketizer=lambda(v):v=='ileal CD')
nod2 = data.get_data(pl, 'nod2', bucketizer=lambda(v):v>0)

print '%d patients with both health and microbiome data' % n

print 'Total bacterial species: %d' % len(data.bacteria)
interesting = {}
cos = {}
pvs = {}
for species in data.bacteria:
    vals = data.get_data(pl, 'fractions', species)
    co = findcutoff(vals, sick)
    if co.sick_when_more==None:
        continue
    boolvals = [(i>co.threshold)==(co.sick_when_more) for i in vals]
    if entropy(boolvals)<.5:
        continue
    cos[species]=co
    cnt = count(zip(boolvals, sick))
    sumarr(cnt, 0.1)
    p = chi2_contingency(cnt)[1]
    pvs[species] = p
    if p > .01:
        continue    
    interesting[species]=boolvals

print "Interesting bacterial species: %d" % len(interesting)

for s1 in interesting:
    v1 = interesting[s1]
    for s2 in interesting:
        if s2>=s1:
            continue
        v2 = interesting[s2]
        bf = direction.montecarlo(interesting[s1], sick, interesting[s2], n).bayes_fwd_rev
        if bf > 4:
            out = '%s & %s & $%.1g$ & ' % (s2, prettyco(cos[s2]), pvs[s2])
            out += '%s & %s & $%.1g$ & ' % (s1, prettyco(cos[s1]), pvs[s1])
            out += '%.1f \\\\' % bf
            out = out.replace('-0','-')
            print out
