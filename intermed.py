#!/usr/bin/python

from utils import *
from copy import copy
from random import shuffle

def test_intermediate(a,b,c):
    val = micond(a,b,c)
    nulldist=[]
    cs=copy(c)
    for i in range(1000):
        shuffle(cs)
        nulldist.append(micond(a,b,cs))
    nulldist.sort()
    return findp(val, nulldist)
