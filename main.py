#!/usr/bin/python

from datareader import Data
from utils import *
from scipy.stats import chi2_contingency
from scipy.stats import chisquare
from math import log
from sys import stdout

def p_of_val(p, v):
    if v:
        return p
    else:
        return 1-p

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
for species in ['Lactobacillus acidophilus']: #data.bacteria:
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
    cnt = count(zip(sick, boolvals))
    print cnt
    chi_indep = chi2_contingency(cnt)[1]
    if chi_indep>1e-3:
        continue
    match_disease += 1
    p_bact_given_cd = [ float(cnt[0][1]) / sum(cnt[0]),
                        float(cnt[1][1]) / sum(cnt[1]) ]
#    print 'p(%s|cd)=%s' % (species[:5], p_bact_given_cd)
    # cntall = count(zip(nod2, sick, boolvals))
    # print 'cntall'
    # print_222(cntall)
    # exp=[[[0,0],[0,0]],[[0,0],[0,0]]]
    # for n in range(2):
    #     for s in range(2):
    #         for b in range(2):
    #             exp[n][s][b] = (len(pl) * 
    #                             p_of_val(p_nod2, n) *
    #                             p_of_val(p_cd_given_nod2[n], s) *
    #                             p_of_val(p_bact_given_cd[s], b))
    # print 'exp'
    # print_222(exp)
    # chi_rev = chisquare(cntall, exp, axis=None, ddof=4) ## ddof???
    # print '%s \t %.2e \t %.4f' % (species, chi_indep, chi_rev.pvalue)
    exp=[[0,0],[0,0]]
    for n in range(2):
        for b in range(2):
            for s in range(2):
                exp[n][b] += (len(pl) * 
                              p_of_val(p_nod2, n) *
                              p_of_val(p_cd_given_nod2[n], s) *
                              p_of_val(p_bact_given_cd[s], b))                
    cnt = count(zip(nod2, boolvals))
    chi_rev = chisquare(cnt, exp, axis=None, ddof=2)
    dirtxt=['<','>']
    if chi_rev.pvalue < .05:
        print '%s \t %s%.2e \t %.2e \t %.4f' % (species, dirtxt[co.sick_when_more], co.threshold, chi_indep, chi_rev.pvalue)    

        
print 'Interesting bacterial species: %d' % interesting
print 'Species that correlated to CD: %d' % match_disease



