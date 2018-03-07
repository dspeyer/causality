#!/usr/bin/python

from utils import *
from random import random, choice

class RandUniq(object):
    def __init__(self):
        self.already = -1
    def get(self):
        while True:
            attempt = random()*.6+.2
            if abs(attempt-self.already)>.3:
                break
        self.already = attempt
        return attempt

def create_net(b_depends_on, c_depends_on):
    p_a = random()
    if b_depends_on == 'a':
        ru=RandUniq()
        p_b = [ru.get(), ru.get()]
    else:
        p_b = [random()]*2
    if c_depends_on == 'a':
        ru=RandUniq()
        p_c = [[ru.get()]*2, [ru.get()]*2]
    elif c_depends_on == 'b':
        ru=RandUniq()
        p_c = [[ru.get(), ru.get()]]*2
    elif c_depends_on == 'both':
        ru=RandUniq()
        p_c = [[ru.get(), ru.get()], [0, ru.get()]]
        p_c[1][0] = ru.get()
    else:
        p_c = [[random()]*2]*2
    return struct(a=p_a,
                  b_given_a=p_b,
                  c_given_ab=p_c)

def randpair():
    p1 = random()/4
    p2 = random()/4
    if random()>.5:
        p1 = 1-p1
    else:
        p2 = 1-p2
    return [p1,p2]


def simulate(p, n):
    out=[[],[],[]]
    if hasattr(p,'d_given_a'):
        out.append([])
        if hasattr(p,'e_given_b'):
            out.append([])
    for i in range(n):
        a=(random()<p.a)
        b=(random()<p.b_given_a[a])
        c=(random()<p.c_given_ab[a][b])
        out[0].append(a)
        out[1].append(b)
        out[2].append(c)
        if hasattr(p,'d_given_a'):
            out[3].append(random()<p.d_given_a[a])
            if hasattr(p,'e_given_b'):
                out[4].append(random()<p.e_given_b[b])
    return out

def create_big_net(n, outlinks):
    out = []
    for i in range(n):
        out.append([0]*n)
    for i in range(n):
        if i < n-outlinks:
            poss = range(i+1, n)
            for j in range(outlinks):
                ol = choice(poss)
                poss.remove(ol)
                out[i][ol] = (random()+1) * (random()<.5 or -1)
    return out

def simulate_big_net(net, n):
    nvar = len(net)
    out=[]
    for i in range(nvar):
        out.append([])
    for i in range(n):
        for var in range(nvar):
            effects = 0
            for cause in range(var):
                val = out[cause][-1]
                factor = net[cause][var]
                effects += (val or -1) * factor
            p = 1 / (1 + exp(-effects))
            out[var].append(random()<p)
    return out
