#!/usr/bin/python

from random import random
from math import log
from utils import p_of_val

def adiff(l):
    return abs(l[0]-l[1])

def pformat(l):
    if type(l)==list:
        return '[' + (', '.join( [pformat(x) for x in l] )) + ']'
    else:
        return '%.2f' % l

def pprint(d):
    k = d.keys()
    k.sort(key=lambda s: (len(s),s) )
    for i in k:
        print '%s: %s' % (i, pformat(d[i]))

def attempt():
    p_a = random()
    p_b_given_a = [random(), random()]
    p_c_given_ab = [[ random(), random() ], [ random(), random() ]]

    p_b = 0
    p_c_given_a = [0, 0]
    p_c_given_b = [0, 0]
    p_ac = [[0,0],[0,0]]
    p_c = 0
    for a in range(2):
        for b in range(2):
            for c in range(2):
                p_abc = ( p_of_val(p_a, a) *
                          p_of_val(p_b_given_a[a], b) *
                          p_of_val(p_c_given_ab[a][b], c) )
                if b:
                    p_b += p_abc
                if c:
                    p_c_given_a[a] += p_abc
                    p_c_given_b[b] += p_abc
                    p_c += p_abc
                p_ac[a][c] += p_abc
    p_c_given_a[0] /= (1-p_a)
    p_c_given_a[1] /= p_a
    p_c_given_b[0] /= (1-p_b)
    p_c_given_b[1] /= p_b

    if adiff(p_c_given_a) > adiff(p_c_given_b):
        #print 'rejecting'
        #pprint(locals())
        global rejects
        rejects += 1
        return 

    p_ac_given_chain = [[0,0],[0,0]]
    for a in range(2):
        for b in range(2):
            for c in range(2):
                p_ac_given_chain[a][c] += ( p_of_val(p_a, a) *
                                            p_of_val(p_b_given_a[a], b) *
                                            p_of_val(p_c_given_b[b], c) )
    p_ac_given_collide = [[0,0],[0,0]]
    for a in range(2):
        for c in range(2):
            p_ac_given_collide[a][c] = p_of_val(p_a, a) * p_of_val(p_c, c)

    score = 0
    for a in range(2):
        for c in range(2):
            score += log(p_ac_given_collide[a][c] / p_ac_given_chain[a][c]) * p_ac[a][c]

    #pprint(locals())
    global falsehoods, still_works
    if score > 0:
        falsehoods += 1
    else:
        still_works += 1


rejects=0
falsehoods=0
still_works=0
for i in range(10000):
    attempt()
print 'rejects: %d' % rejects
print 'falsehoods: %d' % falsehoods
print 'still works: %d' % still_works

