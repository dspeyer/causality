#!/usr/bin/python

from datareader import Data
from utils import *
from scipy.stats import chi2_contingency
from scipy.stats import chi2
from scipy.stats import chisquare
from math import log
from sys import stdout
from random import random
from simulate_causal_net import simulate
from collections import defaultdict

tot=defaultdict(lambda:0)
hit=defaultdict(lambda:0)

cats=20

def record(dir, truth):
    post = 1.0 / (1 + 1/dir.bayes_fwd_rev) # prior=0.5
    rpost = round(post*cats)
    tot[rpost] += 1
    if truth:
        hit[rpost] += 1

n = 58
p_nod2 = 0.39
p_cd_given_nod2=[0.48, 0.91]

i=0
for _ in range(1000):
    net=struct()
    net.a = p_nod2
    p_bact = (random()*.4)+.4
    net.b_given_a=[p_bact, p_bact]
    net.c_given_ab=[[0,0],[0,0]]
    while min(net.c_given_ab[0][1],net.c_given_ab[1][1])<=0 or max(net.c_given_ab[0][1],net.c_given_ab[1][1])>=1:
        net.c_given_ab[0][0] = random()
        net.c_given_ab[1][0] = (net.c_given_ab[0][0] + (random()*.4+.4)) % 1
        for g in [0,1]:
            net.c_given_ab[g][1] = (p_cd_given_nod2[g]-net.c_given_ab[g][0]*(1-p_bact))/p_bact
    #print net
    data = simulate(net, n)
    cnt = count(zip(*data))
#    if any(cnt, lambda(x):x<=5):
#        continue
    #print cnt
    dir=direction(data[0], data[2], data[1], n, p_nod2, p_cd_given_nod2)
    #print dir.bayes_fwd_rev
    record(dir, True)
    i += 1

print '%d fwd' % i

#print "---"

i=0
for _ in range(1000):
    net=struct()
    net.a = p_nod2
    net.b_given_a = p_cd_given_nod2
    net.c_given_ab=[random(), 0]
    net.c_given_ab[1] = (net.c_given_ab[0] + (random()*.4+.4)) % 1
    net.c_given_ab = [net.c_given_ab]*2
    #print net
    data = simulate(net, n)
    cnt = count(zip(*data))
    print cnt
#    if any(cnt, lambda(x):x<=5):
#        continue
    try:
        dir=direction(data[0], data[1], data[2], n, p_nod2, p_cd_given_nod2)
        #print dir.bayes_fwd_rev
        record(dir, False)
        i += 1
    except ValueError:
        print 'ValueError'

print '%d rev' % i

for i in range(cats+1):
    print '%f: %f (%d)' % (i/float(cats), tot[i] and float(hit[i])/tot[i], tot[i])
