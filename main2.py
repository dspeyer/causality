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

pl = data.get_pid_list(lambda(p): (p.health in ['control', 'ileal CD'] and
                                   len(p.fractions)>0))
sick = data.get_data(pl, 'health', bucketizer=lambda(v):v=='ileal CD')

print "patients: %d" % len(pl)

interesting=0
bacteria=[]
for species in data.bacteria:
    vals = data.get_data(pl, 'fractions', species, bucketizer=lambda(x):int(log(x,10)))
    if entropy(vals)<0.5:
        continue
    interesting += 1
    vals = data.get_data(pl, 'fractions', species)
    co = findcutoff(vals, sick)
    if co.sick_when_more==None:
        continue
    boolvals = []
    for i in vals:
        boolvals.append((i>co.threshold)==(co.sick_when_more))
    if link(sick,boolvals)>.01:
        continue
    bacteria.append(struct(species=species,data=boolvals))

print "Total species: %d" % len(data.bacteria)
print "Interesting species: %d" % interesting
print "Correlating species: %d" % len(bacteria)

threshold = .01
for i in range(len(bacteria)):
    for j in range(i):
        for k in range(j):
            a=bacteria[i].data
            b=bacteria[j].data
            c=bacteria[k].data
            try:
                if ( link(a,b) < threshold and
                     link(a,c) < threshold and
                     link(b,c) < threshold and
                     link_despite(a,b,c) < threshold and
                     link_despite(a,c,b) < threshold and
                     link_despite(b,c,a) < threshold):
                    print 'Trio: %s, %s, %s' % (bacteria[i].species, bacteria[j].species, bacteria[k].species)
            except ValueError:
                print 'Error for: %s, %s, %s' % (bacteria[i].species, bacteria[j].species, bacteria[k].species)
