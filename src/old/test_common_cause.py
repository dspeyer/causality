#!/usr/bin/python

from random import random
from simulate_causal_net import simulate
from collections import defaultdict
from common_cause import has_common_cause_bf
from utils import struct, count, any#, severs

tot=defaultdict(lambda:0)
hit=defaultdict(lambda:0)

cats=20

def record(sev, truth):
    post = 1.0 / (1 + 1/sev) # prior=1/2
    rpost = round(post*cats)
    tot[rpost] += 1
    if truth:
        hit[rpost] += 1

for i in range(1):
    net=struct()
    net.a = random()
    net.b_given_a=[random(), random()]
    c_given_b = [random(), random()]
    net.c_given_ab = [c_given_b] * 2
    net.d_given_a = [random(), random()]
    data = simulate(net, 1000)
    cnt = count(zip(*data))
    if any(cnt, lambda(x):x==0):
        print 'net: %s' % net
        print 'skipping %s' % cnt
        continue
    bf = has_common_cause_bf(data[0:3])
    record(bf, False)
    bf = has_common_cause_bf(data[1:4])
    record(bf, True)
    
for i in range(cats+1):
    print '%f: %f n= %d' % (i/float(cats), tot[i] and float(hit[i])/tot[i] or -1, tot[i])
