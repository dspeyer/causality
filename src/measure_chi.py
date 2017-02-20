#!/usr/bin/python

from utils import *
from simulate_causal_net import *
from intermed import test_intermediate
from copy import copy
from random import shuffle
from scipy.stats import chi2_contingency

n=200


# net=create_net('a', 'b')
# data=simulate(net, 1000)
# #print data
# a=0
# b=2
# c=1
# cnt = count(zip(data[c],data[a],data[b]), [2,2,2])
# print cnt
# chi = chi2_contingency(cnt[0])                
# print chi
# chi = chi2_contingency(cnt[1])                
# print chi
# cnt = count(zip(data[a],data[b]), [2,2])
# chi = chi2_contingency(cnt)
# print chi

# exit()

for ll in [20,1000]:
 print 'll=%d'%ll
 print 'b_dep \tc_dep \tran \tbc \t bc|a \tbc|~a \t\tac \tac|b \tac|~b \t\tab \tab|c \tab|~c'
 for bdo in ['nothing', 'a']:
    for cdo in ['nothing', 'a', 'b', 'both']:
        notblock=[0,0,0]
        blockt=[0,0,0]
        blockf=[0,0,0]
        fail=0
        for i in range(n):
            net=create_net(bdo, cdo)
            data=simulate(net, ll)
            for c in range(3):
                a=(c+1)%3
                b=(c+2)%3
                cnt = count(zip(data[c],data[a],data[b]), [2,2,2])
                if deepin(cnt, 0):
                   fail+=1
                   continue
                chi = chi2_contingency(cnt[0])                
                blockf[c] += chi[1]
                chi = chi2_contingency(cnt[1])                
                blockt[c] += chi[1]
                cnt = count(zip(data[a],data[b]), [2,2])
                chi = chi2_contingency(cnt)
                notblock[c] += chi[1]
        fail /= 3
        txt= '%s\t%s\t%d/%d' % (bdo, cdo, (n-fail), n)
        for c in range(3):
          txt += '\t%.2f' % (notblock[c] / (n-fail))
          txt += '\t%.2f' % (blockt[c] / (n-fail))
          txt += '\t%.2f' % (blockf[c] / (n-fail))
          txt += '\t'
        print txt
