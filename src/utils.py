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

# def findcutoff(bact, dis):
#     healthy=[]
#     sick=[]
#     for i in range(len(bact)):
#         if dis[i]:
#             sick.append(bact[i])
#         else:
#             healthy.append(bact[i])
#     sick.sort()
#     healthy.sort(reverse=True)
# #    print 'sick levels %s'%sick
# #    print 'healthy levels %s'%healthy
#     for i in range(int((len(sick)+len(healthy))/4)):
#         if healthy[i]<sick[i]:
# #            print 'cutting: %d, h=%g, s=%g' % (i, healthy[i], sick[i])
#             out=struct(threshold=(healthy[i]+sick[i])/2.0, sick_when_more=True)
#             if healthy[i]==0:
#                 out.threshold=0
#             return out
#     sick.reverse()
#     healthy.reverse()
#     for i in range(int((len(sick)+len(healthy))/4)):
#         if healthy[i]>sick[i]:
# #            print 'cutting: %d, %f, %f' % (i, healthy[i], sick[i])
#             out=struct(threshold=(healthy[i]+sick[i])/2.0, sick_when_more=False)
#             if sick[i]==0:
#                 out.threshold=0
#             return out
#     return struct(sick_when_more=None)

def findcutoff(bact, dis):
    ths=[0]
    vals = sorted(bact)
    for i in range(len(vals)-1):
        if vals[i]!=vals[i+1]:
            ths.append((vals[i]+vals[i+1])/2)
    best=None
    bestscore=-1
    for i in ths:
        score=0
        for b,d in zip(bact,dis):
            if (b>i)==d:
                score += 1
        swm=True
        if score < len(bact)/2:
            score = len(bact)-score
            swm=False
        if score > bestscore:
            best = struct(threshold=i, sick_when_more=swm)
            bestscore = score
    return best

def prettyco(co):
    if co.threshold==0:
        return co.sick_when_more and 'present' or 'absent'
    dir = (co.sick_when_more and '>' or '\\leq')
    return '$%s$%.2E' % (dir, co.threshold)

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

def yates_chi_square(obs, exp):
    if type(obs)==type([]):
        sum=0
        for o,e in zip(obs,exp):
            sum += yates_chi_square(o,e)
        return sum
    return ((abs(obs-exp)-0.5)**2)/exp

def direction(cause, effect, unknown, n, p_cause, p_effect_given_cause):
    cnt = count(zip(effect, unknown))
    #print cnt
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
    #print "obs=%s" % cnt
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

def blurred_chi2_pdf(score, n):
    blurred = n * (chi2.cdf(score+1.0/n, 1) - chi2.cdf(score, 1))
    if blurred > 1e-3: # i.e. is floating point subtraction reliable?
        return blurred
    else:
        return chi2.pdf(score, 1)

def sumall(l):
    if type(l)==type([]):
        return sum([sumall(i) for i in l])
    else:
        return l

def severs(a,b,cut,verbose=False):
    cntall = count(zip(a,b))
    cntcut = count(zip(cut,a,b))
    p_b_given_a = [float(x[1])/sum(x) for x in cntall]
    p_a_given_b = [float(x[1])/sum(x) for x in zip(*cntall)]
    if verbose:
        print 'orig=%s' % cntall
        print 'split=%s' % cntcut
        print 'p_a_given_b = %s' % p_a_given_b
        print 'p_b_given_a = %s' % p_b_given_a
    pvar = count(zip(cut))
    mularr(pvar, 1.0/sum(pvar))
    expnsev = [deepcopy(cntall), deepcopy(cntall)]
    for i in [0,1]:
        mularr(expnsev[i], pvar[i])
    exptoucha = deepcopy(expnsev) # We'll overwrite everything
    for i in [0,1]:
        for aval in [0,1]:
            n = cntcut[i][aval][0] + cntcut[i][aval][1]
            for bval in [0,1]:
                p = p_of_val(p_b_given_a[aval], bval)
                exptoucha[i][aval][bval] = n * p
                ##print 'for cut=%d a=%d b=%d, n=%d p=%.1f val=%.1f' % (i, aval, bval, n, p, exptoucha[i][aval][bval])
    exptouchb = deepcopy(expnsev) # We'll overwrite everything
    for i in [0,1]:
        for bval in [0,1]:
            n = cntcut[i][0][bval] + cntcut[i][1][bval]
            for aval in [0,1]:
                exptouchb[i][aval][bval] = n * p_of_val(p_a_given_b[bval], aval)
    if verbose:
        print 'exp|touch a = %s' % exptoucha
        print 'exp|touch b = %s' % exptouchb
    exps = [expnsev, exptoucha, exptouchb]
    bayes_factor = [1, 1, 1]
    for model in [0,1,2]:
        if verbose:
            print 'Model Touches %s' % (['neither', 'a', 'b'])[model]
        for i in [0,1]:
            try:
                chi_sev = chi2_contingency(cntcut[i])
            except (ValueError, ZeroDivisionError) as e:
                continue
            peg_sev = blurred_chi2_pdf(chi_sev[0], sumall(cntcut[i]))
            if verbose:
                print ' chi_sev=%s' % str(chi_sev)
                print ' p(e|sev)=%f' % peg_sev
                print ' Cut=%d' % i
                print ' Actual: %s' % cntcut[i]
                print ' Expected: %s' % exps[model][i]
            try:
                chi_nsev = chisquare(cntcut[i], exps[model][i], axis=None, ddof=2)
            except (ValueError, ZeroDivisionError) as e:
                print 'Failure for model %d cut %d act=%s exp=%s' % (model, i, cntcut, exps[model])
                raise e
            peg_nsev = blurred_chi2_pdf(chi_nsev[0], sumall(cntcut[i]))
            if verbose:
                print ' Chi=%s' % str(chi_nsev)
                print ' p=%s' % peg_nsev
            bayes_factor[model] *= peg_sev/peg_nsev 
    return min(bayes_factor)

def any(l, cond):
    if type(l)==type([]):
        for i in l:
            if any(i, cond):
                return True
        return False
    return cond(l)
