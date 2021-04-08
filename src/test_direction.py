#!/usr/bin/python

from datareader import Data
from utils import struct
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

def record(dir, truth):
    global logloss
    global count
    post = 1.0 / (1 + 1/dir.bayes_fwd_rev) # prior=0.5
    rpost = round(post*cats)
    tot[rpost] += 1
    if truth:
        hit[rpost] += 1
        logloss += log(post)
    else:
        logloss += log(1-post)

n = 58
p_nod2 = 0.39
p_cd_given_nod2=[0.48, 0.91]

mode=argv[1][1:]
if mode not in ['platonic', 'crohns', 'multi', 'multiplat', 'against-v']:
    print "Bad Mode"
    exit(1)
runs=int(argv[2])
if len(argv)>3:
    n=int(argv[3])

### Collider

i=0
tries=0
while i < runs:
    tries+=1
    if tries % (runs/10)==0:
        print 'try #%d / run #%d ' % (tries, i)
    net=struct()

    if mode in ['platonic', 'multiplat', 'against-v']:
        net.a = random()
    elif mode in ['crohns', 'multi']:
        net.a = p_nod2
    else:
        exit('ERROR: netaf')

    if mode in ['platonic', 'multiplat', 'against-v']:
        p_bact = random()
    elif mode in ['multi', 'crohns']:
        p_bact = (random()*.4)+.4
    else:
        exit('ERROR p(bact) f')
    if mode!='multi':
        net.b_given_a=[p_bact, p_bact]
    else:
        net.b_given_a=[p_bact, p_bact + random()*0.2 - 0.1]

    net.c_given_ab=[[0,0],[0,0]]
    if mode in ['platonic', 'multiplat', 'against-v']:
        for b in range(2):
            for a in range(2):
                net.c_given_ab[a][b]=random()
    elif mode in ['crohns', 'multi']:
        while min(net.c_given_ab[0][1],net.c_given_ab[1][1])<=0 or max(net.c_given_ab[0][1],net.c_given_ab[1][1])>=1:
            net.c_given_ab[0][0] = random()
            net.c_given_ab[1][0] = (net.c_given_ab[0][0] + (random()*.4+.4)) % 1
            for g in [0,1]:
                net.c_given_ab[g][1] = (p_cd_given_nod2[g]-net.c_given_ab[g][0]*(1-p_bact))/p_bact
    else:
        exit('ERROR: netc f')
    data = simulate(net, n)

    if mode in ['crohns', 'multi']:
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

    dir=direction(data[0], data[2], data[1], n)
    record(dir, True)
    i += 1

print '%d fwd (took %d tries)' % (i, tries)

### Chain

i=0
tries=0
while i < runs:
    net=struct()
    tries += 1

    if mode in ['platonic', 'multiplat', 'against-v']:
        net.a = random()
    elif mode in ['crohns', 'multi']:
        net.a = p_nod2
    else:
        exit('ERROR neta ch')

    if mode in ['platonic', 'multiplat', 'against-v']:
        net.b_given_a = [random(), random()]
    elif mode in ['crohns', 'multi']:
        net.b_given_a = p_cd_given_nod2
    else:
        exit('ERROR netb ch')

    if mode in ['platonic', 'crohns']:
        net.c_given_ab=[random(), random()]
        net.c_given_ab = [net.c_given_ab]*2
    elif mode=='multi':
        p_bact_given_cd = [ random() * 0.4 + 0.1 ]
        p_bact_given_cd.append( p_bact_given_cd[0] + random()*0.4  )
        net.c_given_ab=[[ p_bact_given_cd[0], p_bact_given_cd[1] ],
                        [ p_bact_given_cd[0] + random()*0.2 - 0.1, p_bact_given_cd[1] + random()*0.2 - 0.1 ]]
    elif mode=='against-v':
        net.c_given_ab = [ [random()]*2, [random()]*2]
    elif mode=='multiplat':
        net.c_given_ab = [[ random(), random() ], [ random(), random() ]]
    else:
        exit('ERROR netc ch')

    data = simulate(net, n)
    if mode in ['crohns', 'multi']:
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
        if chi2_contingency(count(zip(data[0],data[2])))[1] < chi2_contingency(count(zip(data[0],data[1])))[1]*10:
            continue
    try:
        if mode != 'against-v':
            dir=direction(data[0], data[1], data[2], n)
        else:
            dir=direction(data[1], data[0], data[2], n)
        #print dir.bayes_fwd_rev
        record(dir, False)
        i += 1
    except ValueError:
        print 'ValueError'

print '%d rev (took %d tries)' % (i, tries)

for i in range(cats+1):
    print '%f: %f (%d)' % (i/float(cats), tot[i] and float(hit[i])/tot[i], tot[i])

print 'Average log loss = %.2f' % (logloss/(2*runs))
