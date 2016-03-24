#!/usr/bin/python

from copy import copy
from copy import deepcopy
from collections import defaultdict
from math import log
from scipy.stats import chisquare

class struct(object):
    def __init__(self,**kwargs):
        for k in kwargs:
            setattr(self,k,kwargs[k])
    def __eq__(self, other):
        return self.__dict__ == other.__dict__
    def __contains__(self, key):
        return key in self.__dict__
    def __str__(self):
        return str(self.__dict__)

def count(data, limits):
#    print data
    out=0
    for limit in reversed(limits):
        out=[out]
        for i in range(limit-1):
            out.append(deepcopy(out[0]))
#    print 'out[0] is out[1]: %s' % (out[0] is out[1])
#    print 'out[0][0] is out[1][0]: %s' % (out[0][0] is out[1][0])
#    print out
    for row in data:
        tmp=out
        for i,col in enumerate(row):
            if (i<len(row)-1):
                tmp=tmp[col]
            else:
#                print 'incrementing'
                tmp[col]+=1
    return out

def entropy(l):
    d=defaultdict(lambda:0)
    s=0.0
    ent=0.0
    for i in l:
        d[i]+=1
        s+=1
    for i in d:
        p=d[i]/s
        ent-=p*log(p,2)
    return ent

def normmi(a,b):
    ea=entropy(a)
    eb=entropy(b)
    ej=entropy(map(lambda x:'%s/%s'%x,zip(a,b)))
    return (ea+eb-ej)/max(ea,eb)

def normmicond(a,b,c):
    eac=entropy(map(lambda x:'%s/%s'%x,zip(a,c)))
    ebc=entropy(map(lambda x:'%s/%s'%x,zip(b,c)))
    ec=entropy(c)
    ej=entropy(map(lambda x:'%s/%s/%s'%x,zip(a,b,c)))
    return (eac+ebc-ej-ec)/max(eac,ebc)

def mi(a,b):
    ea=entropy(a)
    eb=entropy(b)
    ej=entropy(map(lambda x:'%s/%s'%x,zip(a,b)))
    return ea+eb-ej

def micond(a,b,c):
    eac=entropy(map(lambda x:'%s/%s'%x,zip(a,c)))
    ebc=entropy(map(lambda x:'%s/%s'%x,zip(b,c)))
    ec=entropy(c)
    ej=entropy(map(lambda x:'%s/%s/%s'%x,zip(a,b,c)))
    return eac+ebc-ej-ec

def findp(v,l):
    mi=0
    ma=len(l)-1
    while ma>mi+1:
        t=int((mi+ma)/2)
        if v>l[t]:
            mi=t
        else:
            ma=t
    return (len(l)-float(ma))/len(l) # subtracting before dividing avoids annoying rounding errors

def deepin(l, v):
    if type(l)==type(v):
        return l==v
    for i in l:
        if deepin(i, v):
            return True
    return False

def mularr(arr, f):
    for i in range(len(arr)):
        if type(arr[i])==type([]):
            mularr(arr[i],f)
        else:
            arr[i]*=f

def severs(a, b, mid):
    cntall=count(zip(a,b), [2,2])
    cntsplit=count(zip(mid,a,b), [2,2,2])
    cntmid=count(zip(mid), [2])
    tot = sum(cntmid)
    mularr(cntmid, 1.0 / tot)
    print 'Everything: %s' % cntall
    for mv in range(2):
        exp=deepcopy(cntall)
        mularr(exp, cntmid[mv])
        print 'Expected: %s' % exp
        print 'Actual: %s' % cntsplit[mv]
        chi = chisquare(cntsplit[mv], exp, axis=None)
        print chi
