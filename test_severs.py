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

def record(sev, truth):
    post = 1.0 / (1 + 2/sev) # prior=1/3
    rpost = round(post*cats)
    tot[rpost] += 1
    if truth:
        hit[rpost] += 1

for i in range(10000):
    net=struct()
    net.a = (random()*.4)+.4
    net.b_given_a=[random()*.6+.2]
    net.b_given_a.append((net.b_given_a[0]+(random()*.4+.4))%1)
    net.c_given_ab=[random()*.6+.2]
    net.c_given_ab.append((net.c_given_ab[0]+(random()*.4+.4))%1)
    net.c_given_ab = [net.c_given_ab] * 2
    data = simulate(net, 100)
    try:
        sev = severs(data[0], data[1], data[2])
        record(sev, False)
        sev = severs(data[0], data[2], data[1])
        record(sev, True)
        sev = severs(data[1], data[2], data[0])
        record(sev, False)
    except (ValueError,ZeroDivisionError):
        pass
        #print 'skipping %s with count %s' % (net, count(zip(*data)))

for i in range(cats+1):
    print '%f: %f (%d)' % (i/float(cats), tot[i] and float(hit[i])/tot[i] or -1, tot[i])
