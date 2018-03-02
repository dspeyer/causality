#!/usr/bin/python
import datareader
from scipy.stats.mstats import zscore
import utils
from scipy.stats import chi2_contingency
from numpy import histogram
from math import log

d = datareader.Data()
d.read_from_files()
pl = d.get_pid_list(lambda(p): ('nod2' in p and 
                                   p.health in ['control', 'ileal CD'] and
                                   len(p.fractions)>0))
data=d

sick = data.get_data(pl, 'health', bucketizer=lambda(v):v=='ileal CD')

ent={}
for species in data.bacteria:
    vals = data.get_data(pl, 'fractions', species, bucketizer=lambda(x):int(log(x,10)))
    ent[species] = utils.entropy(vals)

spe=ent.keys()
spe.sort(key=lambda x: ent[x], reverse=True)
species_by_entropy = spe

pvs={}
for species in data.bacteria:
     vals = data.get_data(pl, 'fractions', species)
     co = utils.findcutoff(vals, sick)
     if co.sick_when_more==None:
         continue
     boolvals = [(i>co.threshold)==(co.sick_when_more) for i in vals]
     cnts = utils.count(zip(sick, boolvals))
     try:
         pvs[species]=chi2_contingency(cnts)[1]
     except ValueError:
         pass

spe=pvs.keys()
spe.sort(key=lambda x: pvs[x])

filt2 = [x for x in species_by_entropy if x in pvs and pvs[x]<.1]
print len(species_by_entropy)
print len(pvs)
print len(ent)
print len(filt2)
examples=[filt2[i] for i in range(0,20,3)]

from scipy.stats.mstats import zscore
zss={}
for i in examples:
 v=data.get_data(pl, 'fractions', i)
 v = [ log(max(x, 1e-6)) for x in v ]
 zss[i]=zscore(v)
 print '%.2f .. %.2f' % (min(zss[i]), max(zss[i]))
 
for i in zss:
 print '%s, %s' % (i,
', '.join([str(x) for x in histogram(zss[i], bins=60, range=[-3,3])[0]])
)
