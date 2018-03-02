import utils
dir(utils)
import datareader
dir(datareader)
help(datareader.Data)
d = datareader.Data()
d.read_from_files()
d.read_from_files
d.read_from_files.__help__
help(d.read_from_files)
import os
os.chdir(..)
os.chdir('..')
d.read_from_files()
pl = data.get_pid_list(lambda(p): ('nod2' in p and 
                                   p.health in ['control', 'ileal CD'] and
                                   len(p.fractions)>0))
pl = d.get_pid_list(lambda(p): ('nod2' in p and 
                                   p.health in ['control', 'ileal CD'] and
                                   len(p.fractions)>0))
pl
data.bacteria
d.bacteria
data=d
ent={}
for species in data.bacteria:
    vals = data.get_data(pl, 'fractions', species, bucketizer=lambda(x):int(log(x,10)))
    ent[species] = utils.entropy(vals)
ent
data.get_data(pl, 'fractions', 'Microlunatus phosphovorus', bucketizer=lambda(x):int(log(x,10)))
data.get_data(pl, 'fractions', 'Microlunatus phosphovorus')
log(100,10)
from math import log
data.get_data(pl, 'fractions', 'Microlunatus phosphovorus')
data.get_data(pl, 'fractions', 'Microlunatus phosphovorus', bucketizer=lambda(x):int(log(x,10)))
data.get_data(pl, 'fractions', 'Microlunatus phosphovorus', bucketizer=lambda(x):int(log(x,10)) or 0)
None or 2
None or 0
data.get_data(pl, 'fractions', 'Microlunatus phosphovorus', bucketizer=lambda(x):(int(log(x,10)) or 0))
int(None)
data.get_data(pl, 'fractions', 'Microlunatus phosphovorus', bucketizer=lambda(x):int(log(x,10) or -10))
for species in data.bacteria:
    vals = data.get_data(pl, 'fractions', species, bucketizer=lambda(x):int(log(x,10)))
    ent[species] = utils.entropy(vals)
ent
spe=keys(ent)
spe=ent.keys()
spe.sort(key=lambda x: ent[x], reverse=True
)
for i in spe[:10]:
 print '%s: %.2f' % (i, ent[i])
species_by_entropy = spe
sick = data.get_data(pl, 'health', bucketizer=lambda(v):v=='ileal CD')
sick
len(sick)
for species in data.bacteria:
    vals = data.get_data(pl, 'fractions', species, bucketizer=lambda(x):int(log(x,10)))
for species in data.bacteria:
    vals = data.get_data(pl, 'fractions', species)
    co = findcutoff(vals, sick)
    if co.sick_when_more==None:
        continue
    boolvals = [(i>co.threshold)==(co.sick_when_more) for i in vals]
from utils import *
for species in data.bacteria:
     vals = data.get_data(pl, 'fractions', species)
     co = findcutoff(vals, sick)
     if co.sick_when_more==None:
         continue
     boolvals = [(i>co.threshold)==(co.sick_when_more) for i in vals]
     cnts = count(zip(sick, boolvals))
     print cnts
from scipy.stats import chi2_contingency
for species in data.bacteria:
     vals = data.get_data(pl, 'fractions', species)
     co = findcutoff(vals, sick)
     if co.sick_when_more==None:
         continue
     boolvals = [(i>co.threshold)==(co.sick_when_more) for i in vals]
     cnts = count(zip(sick, boolvals))
     print chi2_contingency(cnts)[1]
pvs={}
for species in data.bacteria:
     vals = data.get_data(pl, 'fractions', species)
     co = findcutoff(vals, sick)
     if co.sick_when_more==None:
         continue
     boolvals = [(i>co.threshold)==(co.sick_when_more) for i in vals]
     cnts = count(zip(sick, boolvals))
     try:
         pvs[species]=chi2_contingency(cnts)[1]
     except ValueError:
         pass
pvs
spe=pvs.keys()
spe.sort(key=lambda x: pvs[x])
for i in spe[:10]:
 print '%s: %.2f' % (i, pvs[i])
for i in spe[:10]:
 print '%s: %.5f' % (i, pvs[i])
sort_ent_filt_p = [x for x in species_by_entropy if spe[x]<.001]
sort_ent_filt_p = [x for x in species_by_entropy if pvs[spe[x]]<.001]
sort_ent_filt_p = [x for x in species_by_entropy if pvs[]<.001]
sort_ent_filt_p = [x for x in species_by_entropy if pvs[x]<.001]
sort_ent_filt_p = [x for x in species_by_entropy if x in pvs and pvs[x]<.001]
sort_ent_filt_p
filt2 = [x for x in species_by_entropy if x in pvs and pvs[x]<.001 and ent[x]>.5]
filt2
ent['Acinetobacter johnsonii']
data.get_data(pl, 'fractions', 'Acinetobacter johnsonii']
data.get_data(pl, 'fractions', 'Acinetobacter johnsonii')
len(filt2)
[filt2[i] for i in range(0,46,3)]
examples=[filt2[i] for i in range(0,46,3)]
from numpy import *
v=data.get_data(pl, 'fractions', 'Acinetobacter johnsonii')
std(v)
mean(v)
About 711,000 results (0.64 seconds) 
Search Results
import scipy.stats.mstats.zscore
from scipy.stats.mstats import zscore
zscore(v)
zss={}
for i in examples:
 v=data.get_data(pl, 'fractions', i)
 zss[i]=zscore(v)
zss
max(zss)
[max(zss[i]) for i in zss]
max([max(zss[i]) for i in zss])
min([min(zss[i]) for i in zss])
histogram
histogram(zss['Sporobacter termitidis'], bins=90, range=[-1,8])
histogram(zss['Sporobacter termitidis'], bins=90, range=[-1,8])[0]
for i in zss:
 print '%s, %s' % (i, 
histogram(zss[i], bins=90, range=[-1,8])[0]
)
for i in zss:
histogram(zss['Sporobacter termitidis'], bins=90, range=[-1,8])
histogram(zss['Sporobacter termitidis'], bins=90, range=[-1,8])[0]
[str(x) for x in histogram(zss['Sporobacter termitidis'], bins=90, range=[-1,8])[0]]
', '.join([str(x) for x in histogram(zss['Sporobacter termitidis'], bins=90, range=[-1,8])[0]])
for i in zss:
 print '%s, %s' % (i,
', '.join([str(x) for x in histogram(zss[i], bins=90, range=[-1,8])[0]])
)
2.3 % 1
import readline
readline.write_history_file('create_nonnormal_graph.py')
