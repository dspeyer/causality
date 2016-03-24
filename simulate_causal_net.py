#!/usr/bin/python

from utils import *
from random import random

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

def simulate(p, n):
    out=[[],[],[]]
    for i in range(n):
        a=(random()<p.a)
        b=(random()<p.b_given_a[a])
        c=(random()<p.c_given_ab[a][b])
        out[0].append(a)
        out[1].append(b)
        out[2].append(c)
    return out
