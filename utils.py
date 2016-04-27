#!/usr/bin/python

from copy import copy
from copy import deepcopy
from collections import defaultdict
from math import log
from scipy.stats import chisquare
from scipy.stats import chi2_contingency
from scipy.stats import chi2

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
    def __eq__(self, other):
        return self.__dict__ == other.__dict__
    def __ne__(self, other):
        return self.__dict__ != other.__dict__

def count(data, limits=None):
#    print data
    if limits==None:
        limits=[2]*len(data[0])
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

# def severs(a, b, mid):
#     cntall=count(zip(a,b), [2,2])
#     cntsplit=count(zip(mid,a,b), [2,2,2])
#     cntmid=count(zip(mid), [2])
#     tot = sum(cntmid)
#     mularr(cntmid, 1.0 / tot)
#     print 'Everything: %s' % cntall
#     for mv in range(2):
#         exp=deepcopy(cntall)
#         mularr(exp, cntmid[mv])
#         print 'Expected: %s' % exp
#         print 'Actual: %s' % cntsplit[mv]
#         chi = chisquare(cntsplit[mv], exp, axis=None)
#         print chi

def findcutoff(bact, dis):
    healthy=[]
    sick=[]
    for i in range(len(bact)):
        if dis[i]:
            sick.append(bact[i])
        else:
            healthy.append(bact[i])
    sick.sort()
    healthy.sort(reverse=True)
#    print 'sick levels %s'%sick
#    print 'healthy levels %s'%healthy
    for i in range(min(len(sick),len(healthy))):
        if healthy[i]<sick[i]:
            return struct(threshold=(healthy[i]+sick[i])/2.0, sick_when_more=True)
    sick.reverse()
    healthy.reverse()
    for i in range(min(len(sick),len(healthy))):
        if healthy[i]>sick[i]:
            return struct(threshold=(healthy[i]+sick[i])/2.0, sick_when_more=False)
    return struct(sick_when_more=None)

def expect(val, ex, desc):
    if type(ex)==type(lambda:0):
        exf = ex
    else:
        exf = (lambda x: x==ex)
    if not exf(val):
        print 'Bad %s: expected %s got %s' % (desc, ex, val)


def p_of_val(p, v):
    if v:
        return p
    else:
        return 1-p

def direction(cause, effect, unknown, n, p_cause, p_effect_given_cause):
    cnt = count(zip(effect, unknown))
    chi_indep = chi2_contingency(cnt)[1]
    p_unknown_given_effect = [ float(cnt[0][1]) / sum(cnt[0]),
                               float(cnt[1][1]) / sum(cnt[1]) ]
    #print 'p(bact|cd)=%s' % p_unknown_given_effect
    exp=[[0,0],[0,0]]
    for c in range(2):
        for e in range(2):
            for u in range(2):
                exp[c][u] += (n * 
                              p_of_val(p_cause, c) *
                              p_of_val(p_effect_given_cause[c], e) *
                              p_of_val(p_unknown_given_effect[e], u))
    cnt = count(zip(cause, unknown))
    #print 'cnt=%s' % cnt
    #print 'expected if cd->bact=%s' % exp
    chi_rev = chisquare(cnt, exp, axis=None, ddof=2)
    chi_fwd = chi2_contingency(cnt)
    #print 'expected if bact->cd=%s' % chi_fwd[3]
    bayes_factor = chi2.pdf(chi_fwd[0],1) / chi2.pdf(chi_rev.statistic,1)
    return struct(reject_indep=chi_indep,
                  bayes_fwd_rev=bayes_factor,
                  reject_fwd=chi_fwd[1],
                  reject_rev=chi_rev.pvalue)


def div(a, b):
    if type(a)==type(0):
        return float(a)/b
    else:
        return [div(ai,bi) for (ai,bi) in zip(a,b)]


def link(a,b):
    cnt = count(zip(a, b))
    return chi2_contingency(cnt)[1]

def link_despite(a,b,despite):
    cnt = count(zip(despite, a, b))
    ps=[0,0]
    for i in [0,1]:
        try:
            ps[i]=chi2_contingency(cnt[i])[1]
        except ValueError:
            ps[i]=1
    return min(ps)

def severs(a,b,cut):
    cntall = count(zip(a,b))
    cntcut = count(zip(cut,a,b))
    pvar = count(zip(cut))
    mularr(pvar, 1.0/sum(pvar))
    expnsev = [deepcopy(cntall), deepcopy(cntall)]
    print expnsev
    print pvar
    for i in [0,1]:
        mularr(expnsev[i], pvar[i])
    bayes_factor = 1
    for i in [0,1]:
        chi_sev = chi2_contingency(cntcut[i])
        chi_nsev = chisquare(cntcut[i], expnsev[i], axis=None, ddof=2)
        peg_sev = chi2.pdf(chi_sev[0],1) 
        peg_nsev = chi2.pdf(chi_nsev.statistic,1)
        print "For sevval %d" % i
        print "Seen     = %s" % cntcut[i]
        print "Exp|Sev = %s" % chi_sev[3]
        print "Exp|~Sev  = %s" % expnsev[i]
        print "bf = %g / %g = %g" % (peg_sev, peg_nsev, peg_sev/peg_nsev)
        bayes_factor *= peg_sev/peg_nsev 
    return bayes_factor
