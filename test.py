#!/usr/bin/python

from datareader import Data
from simulate_causal_net import *
from utils import *

data = Data()
data.read_from_files()

pl = data.get_pid_list(lambda(p):p.sex in ['male','female'])
sexes = data.get_data(pl,'sex')
for i in range(len(pl)):
    if sexes[i] != data.patients[pl[i]].sex:
        print 'ERROR (wrong sex) on %d' % i

pl = data.get_pid_list()
nod2 = data.get_data(pl, 'nod2', bucketizer=lambda(v):v>0)
for pid,n in zip(pl,nod2):
    if n==None and 'nod2' in data.patients[pid]:
        print 'ERROR (present nod2) on %d' % i
    if n==False and data.patients[pid].nod2!=0:
        print 'ERROR (nod2!=0) on %d' % i
    if n==True and data.patients[pid].nod2 not in [1,2]:
        print 'ERROR (bad nod2) on %d' % i
    if n not in [None, False, True]:
        print 'ERROR n=%s on %d' % (n, i)

DDEHA = 'Desulfitobacterium dehalogenans'
ddeha = data.get_data(pl, 'fractions', DDEHA)
for pid,dd in zip(pl,ddeha):
    patient = data.patients[pid]
    if 'fractions' not in patient or DDEHA not in patient.fractions:
        if dd!=None:
            print 'ERROR (dd not none) %d' % i
    elif patient.fractions[DDEHA]!=dd:
        print 'ERROR (wrong ddeha) %d' % i

rows = [[0,0],
        [0,1],
        [0,2],
        [1,1],
        [0,1]]
cnts = count(rows, [2,3])
if (cnts != [[1, 2, 1],
             [0, 1, 0]]):
    print 'Counting Error'


net = create_net('a','b')
print net
print

data = simulate(net, 100)
print 'Really everything: %s' % count(zip(data[0],data[2],data[1]),[2,2,2])
severs(data[0], data[2], data[1])
print
print 'Really everything: %s' % count(zip(data[0],data[1],data[2]),[2,2,2])
severs(data[0], data[1], data[2])
