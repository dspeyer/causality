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

prefix = ' '.join(argv[1:])

for species in data.bacteria:
    if prefix not in species:
        continue
    
    vals = data.get_data(pl, 'fractions', species)
    co = findcutoff(vals, sick)
    boolvals = [(i>co.threshold)==(co.sick_when_more) for i in vals]
    cnt = count(zip(boolvals, sick))
    sumarr(cnt, 0.1)
    p = chi2_contingency(cnt)[1]
    bf = direction.montecarlo(nod2, sick, boolvals, len(nod2)).bayes_fwd_rev

    print species
    print 'Sick when: %s' % prettyco(co)
    print 'P.V.: %f' % p
    print 'B.F.: %f' % bf
