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
nod2 = data.get_data(pl, 'nod2', bucketizer=lambda(v):v>0)

print '%d patients with both gene and microbiome data' % len(pl)

cnt=count(zip(nod2, sick)) 
chi = chi2_contingency(cnt) #.0022

p_nod2 = float(sum(nod2))/len(nod2)
p_cd_given_nod2 = [ float(cnt[0][1]) / sum(cnt[0]),
                    float(cnt[1][1]) / sum(cnt[1])]

print 'p(nod2)=%s' % p_nod2
print 'p(cd|nod2)=%s' % p_cd_given_nod2


print 'P-value rejecting nod2 indep CD: %.4f' % chi[1]

print 'Total bacterial species: %d' % len(data.bacteria)
interesting = 0
match_disease = 0
dirs={}
for species in data.bacteria:
    vals = data.get_data(pl, 'fractions', species, bucketizer=lambda(x):int(log(x,10)))
#    if entropy(vals)<0.5:
#        continue
#    interesting += 1
    vals = data.get_data(pl, 'fractions', species)
    co = findcutoff(vals, sick)
    if co.sick_when_more==None:
        continue
    boolvals = [(i>co.threshold)==(co.sick_when_more) for i in vals]

    try:
        dir=direction(nod2, sick, boolvals, len(pl), p_nod2, p_cd_given_nod2)
        dirs[species]=dir.bayes_fwd_rev
    except ValueError:
        pass
#        linknod2 = chi2_contingency(count(zip(boolvals,nod2)))[1]
#    if dir.reject_indep < .05 and linknod2 < .05:
#        print 'intert %s  %.2e  %.2e  %.1f' % (species, dir.reject_indep, linknod2, severs(nod2,sick,boolvals))
#    elif dir.reject_indep < 1e-2 and dir.bayes_fwd_rev > 2:
#        print 'interact %s  %s  %.2e  %.1f' % (species, prettyco(co), dir.reject_indep, dir.bayes_fwd_rev)
#    if dir.reject_indep < 1e-2:
#        match_disease += 1
        
print 'Interesting bacterial species: %d' % interesting
print 'Species that correlated to CD: %d' % match_disease

print


pl = data.get_pid_list(lambda(p): (p.health in ['control', 'ileal CD'] and
                                   len(p.fractions)>0))
sick = data.get_data(pl, 'health', bucketizer=lambda(v):v=='ileal CD')

print "patients: %d" % len(pl)

interesting=0
bacteria=[]
bv={}
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
    bacteria.append(struct(species=species,data=boolvals,co=co))
    bv[species]=boolvals

print "Total species: %d" % len(data.bacteria)
print "Interesting species: %d" % interesting
print "Correlating species: %d" % len(bacteria)

if False:
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
                    for ii in [i,j,k]:
                        for hi in [i,j,k]:
                            if ii==hi:
                                continue
                            s=set([i,j,k])
                            s.remove(ii)
                            s.remove(hi)
                            oi=s.pop()
                            interest=bacteria[ii].data
                            helper=bacteria[hi].data
                            if ( severs(helper, sick, interest) > 100 and 
                                 severs(interest, sick, helper) < 1):
                                print '%s, %s, %s, %s, %.1f, %.1g, %.1g' % (
                                    bacteria[ii].species, bacteria[hi].species, bacteria[oi].species, prettyco(bacteria[oi].co),
                                    severs(helper, sick, interest),
                                    max( link(a,b), link(a,c), link(b,c), link_despite(a,b,c), link_despite(a,c,b), link_despite(b,c,a)), dirs[bacteria[ii].species])
                                                                            
            except ValueError:
                print 'Error for: %s, %s, %s' % (bacteria[i].species, bacteria[j].species, bacteria[k].species)


print severs(bv['Pseudoflavonifractor capillosus'], sick, bv['Bacteroides xylanolyticus'], verbose=True)
