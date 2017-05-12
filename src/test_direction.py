#!/usr/bin/python

from datareader import Data
from utils import *
from scipy.stats import chi2_contingency
from scipy.stats import chi2
from scipy.stats import chisquare
from math import log
from sys import stdout, argv, exit
from random import random
from simulate_causal_net import simulate
from collections import defaultdict
from direction import montecarlo as direction

tot=defaultdict(lambda:0)
hit=defaultdict(lambda:0)

cats=20

logloss = 0
count = 0

def record(dir, truth):
    global logloss
    global count
    post = 1.0 / (1 + 1/dir.bayes_fwd_rev) # prior=0.5
    rpost = round(post*cats)
    tot[rpost] += 1
    count += 1
    if truth:
        hit[rpost] += 1
        logloss += log(post)
    else:
        logloss += log(1-post)

n = 58
p_nod2 = 0.39
p_cd_given_nod2=[0.48, 0.91]

mode=argv[1][1:]
if mode not in ['platonic', 'crohns', 'multi', 'multiplat']:
    print "Bad Mode"
    exit(1)
runs=int(argv[2])
if len(argv)>3:
    n=int(argv[3])

if mode=='multiplat':
    n=1000

### Collider

i=0
tries=0
while i < runs:
    tries+=1
    if tries % (runs/10)==0:
        print 'try #%d / run #%d ' % (tries, i)
    net=struct()

    if mode in ['platonic', 'multiplat']:
        net.a = random()
    else:
        net.a = p_nod2

    if mode in ['platonic', 'multiplat']:
        p_bact = random()
    else:
        p_bact = (random()*.4)+.4
    if mode!='multi':
        net.b_given_a=[p_bact, p_bact]
    else:
        net.b_given_a=[p_bact, p_bact + random()*0.2 - 0.1]

    net.c_given_ab=[[0,0],[0,0]]
    if mode in ['platonic', 'multiplat']:
        for b in range(2):
            r=random()
            for a in range(2):
                net.c_given_ab[a][b]=r
    else:
        while min(net.c_given_ab[0][1],net.c_given_ab[1][1])<=0 or max(net.c_given_ab[0][1],net.c_given_ab[1][1])>=1:
            net.c_given_ab[0][0] = random()
            net.c_given_ab[1][0] = (net.c_given_ab[0][0] + (random()*.4+.4)) % 1
            for g in [0,1]:
                net.c_given_ab[g][1] = (p_cd_given_nod2[g]-net.c_given_ab[g][0]*(1-p_bact))/p_bact

    data = simulate(net, n)

    if mode not in ['platonic', 'multiplat']:
        cnt = count(zip(*data))
        if any(cnt, lambda(x):x==0):
            continue
        if chi2_contingency(count(zip(data[0],data[2])))[1]>.01:
            continue
        if mode=='multi' and chi2_contingency(count(zip(data[0],data[1])))[1]<.1:
            continue
    elif mode=='multiplat':
        cnt = count(zip(*data))
        if any(cnt, lambda(x):x==0):
            continue
        if chi2_contingency(count(zip(data[0],data[1])))[1] < chi2_contingency(count(zip(data[0],data[2])))[1]:
            continue

    dir=direction(data[0], data[2], data[1], n, p_nod2, p_cd_given_nod2)
    record(dir, True)
    i += 1

print '%d fwd (took %d tries)' % (i, tries)

if mode=='multiplat':
    for i in range(cats+1):
        print '%f: %f (%d)' % (i/float(cats), tot[i] and float(hit[i])/tot[i], tot[i])
    tot=defaultdict(lambda:0)
    hit=defaultdict(lambda:0)


### Chain

i=0
tries=0
while i < runs:
    net=struct()
    tries += 1

    if mode in ['platonic', 'multiplat']:
        net.a = random()
    else:
        net.a = p_nod2

    if mode in ['platonic', 'multiplat']:
        net.b_given_a = [random(), random()]
    else:
        net.b_given_a = p_cd_given_nod2

    if mode in ['platonic', 'crohns']:
        net.c_given_ab=[random(), random()]
        net.c_given_ab = [net.c_given_ab]*2
    elif mode=='multi':
        p_bact_given_cd = [ random() * 0.4 + 0.1 ]
        p_bact_given_cd.append( p_bact_given_cd[0] + random()*0.4  )
        net.c_given_ab=[[ p_bact_given_cd[0], p_bact_given_cd[1] ],
                        [ p_bact_given_cd[0] + random()*0.2 - 0.1, p_bact_given_cd[1] + random()*0.2 - 0.1 ]]
    else:
        net.c_given_ab = [[ random(), random() ], [ random(), random() ]]

    data = simulate(net, n)
    if mode not in ['platonic', 'multiplat']:
        cnt = count(zip(*data))
        if any(cnt, lambda(x):x==0):
            continue
        if chi2_contingency(count(zip(data[0],data[1])))[1]>.001:
            continue
        if mode=='multi' and chi2_contingency(count(zip(data[0],data[2])))[1]<.01:
            continue
    elif mode=='multiplat':
        cnt = count(zip(*data))
        if any(cnt, lambda(x):x==0):
            continue
        if chi2_contingency(count(zip(data[0],data[2])))[1] < chi2_contingency(count(zip(data[0],data[1])))[1]:
            continue
    try:
        dir=direction(data[0], data[1], data[2], n, p_nod2, p_cd_given_nod2)
        #print dir.bayes_fwd_rev
        record(dir, False)
        i += 1
    except ValueError:
        print 'ValueError'

print '%d rev (took %d tries)' % (i, tries)

for i in range(cats+1):
    print '%f: %f (%d)' % (i/float(cats), tot[i] and float(hit[i])/tot[i], tot[i])

print 'Average log loss = %.2f' % (logloss/count)
