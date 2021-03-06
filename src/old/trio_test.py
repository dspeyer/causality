#!/usr/bin/python

from utils import count
from scipy.stats import chi2_contingency

def rotate(l, n):
    return l[n:] + l[:n]

def has_common_cause(vs,thresh):
    #print 'overall %s' % count(zip(*vs))
    counts = [ count(zip(*rotate(vs, i))) for i in range(3) ]
    counts2 = [ count(zip(vs[i],vs[(i+1)%3])) for i in range(3) ]
    for i in range(3):
        p=chi2_contingency(counts2[i])[1]
        if p >= thresh:
            return False
        #print 'checking corr %s' % counts[i]
        corr = False
        for sv in range(2):
            p=chi2_contingency(counts[i][sv])[1]
            #print p
            if p < thresh:
                corr=True
                break
        if not corr:
            return False
    return True

def p_given(a,b):
    cnt = count(zip(b,a))
    return [ float(x[1]) / sum(x) for x in cnt ]

def has_common_cause_bf(vs):
    dir1 = sum([ v[0]==v[1] for v in zip(*vs) ]) > len(vs[0])/2
    dir2 = sum([ v[0]==v[2] for v in zip(*vs) ]) > len(vs[0])/2
    h = [ (v[0] + (v[1]==dir1) + (v[2]==dir2)) > 1 for v in zip(*vs)]
    p_given_prev = [[],[],[]]
    p_given_h = [[],[],[]]
    for i in range(3):
        p_given_prev[i] = p_given(vs[i], vs[(i-1)%3])
        p_given_h = p_given(vs[i], h)
    logp_ch = [0,0,0]
    logp_cc = 0
    for r in range(len(vs[0])):
        for i in range(3):
            logp_ch += log(p_of_val(p_given_h[i][h[r]], vs[i][r]))
            for hy in range(3):
                logp_ch[hy] += log(p_of_val(p_given_prev[hy][vs[(i-1)%3][r]], vs[i][r]))
    return exp(logp_cc - max(logp_ch))
