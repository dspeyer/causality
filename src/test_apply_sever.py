#!/usr/bin/python

from simulate_causal_net import simulate, randpair
from common_cause import has_common_cause
from utils import struct
from collections import defaultdict
from severs import montecarlo as severs
from random import random

tot=defaultdict(lambda:0)
hit=defaultdict(lambda:0)

cats=20

everything=0
logloss=0

def record(sev, truth):
    global everything
    global logloss
    post = 1.0 / (1 + 1/sev) # prior=1/3
    rpost = round(post*cats)
    tot[rpost] += 1
    everything += 1
    if truth:
        hit[rpost] += 1
        logloss += log(post)
    else:
        logloss += log(1-post)

for run in range(1):
    if run%50 == 0:
        print run
    net=struct()
    net.a = random()/5 + .4
    net.b_given_a = randpair()
    c_given_a = randpair() 
    net.c_given_ab = [c_given_a] * 2
    net.d_given_a = randpair()
    net.e_given_b = randpair()
    data = simulate(net, 500)
    print 'ld=%d'%len(data)
    data3 = [ row[1:4] for row in data ]
    print data3
    if not has_common_cause(data3, 0.001):
        print 'no common cause'
        continue
    datat = zip(*data)
    for b1 in range(1,4):
        for b2 in range(1,4):
            if b1==b2:
                continue
            sev = severs(datat[b2], datat[4], datat[b1])
            should = (b1 == 1)
            record(sev,should)

for i in range(cats+1):
    print '%f: %f n= %d' % (i/float(cats), tot[i] and float(hit[i])/tot[i] or -1, tot[i])

print 'Average log loss = %.2f' % (logloss/(everything))
