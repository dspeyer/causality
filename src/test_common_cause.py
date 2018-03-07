#!/usr/bin/python

from simulate_causal_net import simulate, randpair
from random import random
from common_cause import has_common_cause
from utils import struct


for pv in [0.001]:
    for n in [58, 100, 500]:
        print
        print 'pv=%f, n=%d' % (pv,n)

        found = 0
        for run in range(10):
            net=struct()
            net.a = random()/5 + .4
            net.b_given_a = randpair()
            c_given_a = randpair() 
            net.c_given_ab = [c_given_a] * 2
            net.d_given_a = randpair()
            data = simulate(net, n)
            data = data[1:]
            print [len(data),len(data[0])]
            if has_common_cause(data, pv):
                found += 1
        print 'True common cause %d' % found

        found = 0
        for run in range(1000):
            net=struct()
            net.a = random()/5 + .4
            net.b_given_a = randpair()
            c_given_b = randpair()
            net.c_given_ab = [ [x,x] for x in c_given_b ]
            data = simulate(net, n)
            if has_common_cause(data, pv):
                found += 1
        print 'Chain %d' % found
