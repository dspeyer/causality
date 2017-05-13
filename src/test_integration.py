#!/usr/bin/python

from simulate_causal_net import create_big_net, simulate_big_net
from scipy.stats import chi2_contingency
from utils import count
from trio_test import has_common_cause

vs = 6

net = create_big_net(vs, 3)
data = simulate_big_net(net, 1000)

for i in range(vs):
    for j in range(vs):
        if net[i][j]:
            cnt = count(zip(data[i], data[j]))
            p = chi2_contingency(cnt)[1]
            print '%d->%d %.2f p<%.2f' % (i,j,net[i][j], p)

for i in range(vs):
    for j in range(vs):
        if i!=j and not net[i][j]:
            cnt = count(zip(data[i], data[j]))
            p = chi2_contingency(cnt)[1]
            if p<.05:
                print '%d,%d correlate' % (i,j)

            
for i in range(vs):
    for j in range(i):
        for k in range(j):
            #print 'testing %s' % [i,j,k]
            if has_common_cause([data[x] for x in [i,j,k]], .05):
                print 'found %s' % [i,j,k]
