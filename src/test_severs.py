#!/usr/bin/python

from datareader import Data
from scipy.stats import chi2_contingency
from scipy.stats import chi2
from scipy.stats import chisquare
from math import log
from sys import stdout
from random import random
from simulate_causal_net import simulate
from collections import defaultdict
from severs import montecarlo as severs
from utils import struct, count, any#, severs

tot=defaultdict(lambda:0)
hit=defaultdict(lambda:0)

cats=20

def record(sev, truth):
    post = 1.0 / (1 + 1/sev) # prior=1/3
    rpost = round(post*cats)
    tot[rpost] += 1
    if truth:
        hit[rpost] += 1

for i in range(10000):
    net=struct()
#    net.a = (random()*.4)+.4
#    net.b_given_a=[random()*.6+.2]
#    net.b_given_a.append((net.b_given_a[0]+(random()*.4+.4))%1)
#    c_given_b=[random()*.6+.2]
#    c_given_b.append((c_given_b[0]+(random()*.4+.4))%1)
    net.a = random()
    net.b_given_a = [random(), random()]
    c_given_b = [random(), random()]
    net.c_given_ab = [c_given_b] * 2
    data = simulate(net, 100)
    cnt = count(zip(*data))
    if any(cnt, lambda(x):x==0):
        continue
    try:
#        print "--- %s ---" % net
#        sev = severs(data[0], data[1], data[2])
#        record(sev, False)
        #        if sev>10:
        #            print net
#        print "Severs:"
        sev = severs(data[0], data[2], data[1])
        record(sev, True)
#        sev = severs(data[1], data[2], data[0])
#        record(sev, False)
    except (ValueError,ZeroDivisionError):
        print 'skipping %s with count %s' % (net, count(zip(*data)))

    net=struct()
    net.a = random()
    net.b_given_a = [random(), random()]
    net.c_given_ab = [[random()]*2, [random()]*2]
    net.d_given_a = [random(), random()]
    data = simulate(net, 100)
    sev = severs(data[1],data[2],data[3])
    record(sev,False)

for i in range(cats+1):
    print '%f: %f n= %d' % (i/float(cats), tot[i] and float(hit[i])/tot[i] or -1, tot[i])
