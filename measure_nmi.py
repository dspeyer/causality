#!/usr/bin/python

from utils import *
from simulate_causal_net import *
from intermed import test_intermediate
from copy import copy
from random import shuffle

n=200

for ll in [50,100,1000]:
 print 'll=%d'%ll
 print 'b_dep \tc_dep \tbc \t bc|a \tbc|s(a) \tac \tac|b \tac|s(b) \tab \tab|c \tab|s(c)'
#print 'b_dep \tc_dep \ta int \t b int \t c int'

 for bdo in ['nothing', 'a']:
    for cdo in ['nothing', 'a', 'b', 'both']:
        notblock=[0,0,0]
        block=[0,0,0]
        sblock=[0,0,0]
#        pv=[0,0,0]
        for i in range(n):
            net=create_net(bdo, cdo)
            data=simulate(net, 100)
            for c in range(3):
                a=(c+1)%3
                b=(c+2)%3
                notblock[c]+=mi(data[a],data[b])
                block[c]+=micond(data[a],data[b],data[c])
                #pv[c]+=test_intermediate(data[a],data[b],data[c])
                cs=copy(data[c])
                shuffle(cs)
                sblock[c]+=micond(data[a], data[b], cs)
        txt=bdo+'\t'+cdo
        for c in range(3):
          txt += '\t%.2f' % notblock[c]
          txt += '\t%.2f' % block[c]
          if block[c] < sblock[c] - 1:
            txt+= ' <'
          if block[c] > sblock[c] + 1:
            txt+= ' >'
          txt += '\t%.2f' % sblock[c]
          #txt += '\t%.2f' % (pv[c]/n)
          txt += '\t'
        print txt
