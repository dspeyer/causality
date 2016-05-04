#!/usr/bin/python

from datareader import Data
from utils import *
from scipy.stats import chi2_contingency
from scipy.stats import chi2
from scipy.stats import chisquare
from math import log
from sys import stdout

def print_222(data):
    txt=['! ','']
    print '\t\tHealthy\tCD'
    for b in range(2):
        stdout.write(txt[b]+'Ba')
        for n in range(2):
            stdout.write('\t'+txt[n]+'Nod2\t')
            for d in range(2):
                stdout.write('%d\t' % data[n][d][b])
            stdout.write('\n')


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
print '%s \t sickwhen \t rej ind \t bayes' % ('species'.ljust(32))
for species in data.bacteria: #['Lactobacillus acidophilus']:#
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
    dir=direction(nod2, sick, boolvals, len(pl), p_nod2, p_cd_given_nod2)
    dirtxt=['<','>']
    if dir.reject_indep < 1e-2 and dir.bayes_fwd_rev > 2:
        print '%s,  %s%.2e,  %.2e,  %.1f,  %.1g,  %.1g' % (species.ljust(32), dirtxt[co.sick_when_more], co.threshold, dir.reject_indep, dir.bayes_fwd_rev, dir.reject_fwd, dir.reject_rev)
#        print '%s \t& %s%.2e \t& %.2e \t& %.1f' % (species.ljust(32), dirtxt[co.sick_when_more], co.threshold, dir.reject_indep, dir.bayes_fwd_rev)
#        print 'p_cd_given_n2_b = %s' % div(count(zip(sick,nod2,boolvals))[1], count(zip(nod2,boolvals)))
#    if dir.reject_indep < 1e-2 and dir.bayes_fwd_rev < 2:
#        print 'p_b_given_cd = %s' % div(count(zip(boolvals,sick))[1], count(zip(sick)))
        
print 'Interesting bacterial species: %d' % interesting
print 'Species that correlated to CD: %d' % match_disease



