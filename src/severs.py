#!/usr/bin/python

from utils import count, p_of_val, deepcopy
from math import log, exp
import numpy as np
from numpy.random import beta


def logp_obs_given(cnt, p):
    out = 0
    for a in range(2):
        for b in range(2):
            out += log(p[a][b]) * cnt[a][b]
    return out

def mle(a,b,cut,verbose=False):
    cntall = count(zip(a,b))
    cntcut = count(zip(cut,a,b))
    sumarr(cntall, 0.1)
    sumarr(cntcut, 0.1)
    p_b_given_a = [float(x[1])/sum(x) for x in cntall]
    p_a_given_b = [float(x[1])/sum(x) for x in zip(*cntall)]
    logbfs=[0,0]
    for cutv in range(2):
        cnt = cntcut[cutv]
        tot = sum([sum(l) for l in cnt])
        p_a = float( sum(cnt[1]) ) / tot
        p_b = float( cnt[0][1] + cnt[1][1] ) / tot
        p_ab_given_cuts = [[0,0],[0,0]]
        p_ab_given_toucha = [[0,0],[0,0]]
        p_ab_given_touchb = [[0,0],[0,0]]
        for av in range(2):
            for bv in range(2):
                p_ab_given_cuts[av][bv] = p_of_val(p_a, av) * p_of_val(p_b, bv)
                p_ab_given_toucha[av][bv] = ( p_of_val(p_a, av) *
                                              p_of_val(p_b_given_a[av], bv) )
                p_ab_given_touchb[av][bv] = ( p_of_val(p_b, bv) *
                                              p_of_val(p_a_given_b[bv], av) )
        logp_obs_given_cuts = logp_obs_given(cnt, p_ab_given_cuts)
        alternatives=[p_ab_given_toucha, p_ab_given_touchb]
        for (i,alternative) in enumerate(alternatives):
            logp_obs_given_alt = logp_obs_given(cnt, alternative)
            logbf = logp_obs_given_cuts - logp_obs_given_alt
            logbfs[i] += logbf
    return exp(min(logbfs))

                

def montecarlo(a,b,cut,verbose=False):
    cntall = count(zip(a,b))
    cntcut = count(zip(cut,a,b))
    cntall += 0.1
    cntcut += 0.1
    logbfs=[0,0]
    runs=10
    for cutv in range(2):
        cnt = cntcut[cutv]
        cnt_a = cnt.sum(axis=1)
        cnt_b = cnt.sum(axis=0)
        tot = cnt.sum()
        p_ab_given_cuts = np.zeros([2,2])
        p_ab_given_toucha = np.zeros([2,2])
        p_ab_given_touchb = np.zeros([2,2])
        for i in range(runs):
            p_a = 1-beta(*cnt_a)
            p_b = 1-beta(*cnt_b)
            p_a_given_b = [ 1-beta(*l) for l in zip(*cntall) ]
            p_b_given_a = [ 1-beta(*l) for l in cntall ]
            for av in range(2):
                for bv in range(2):
                    p_ab_given_cuts[av][bv] += p_of_val(p_a, av) * p_of_val(p_b, bv) / runs
                    p_ab_given_toucha[av][bv] += ( p_of_val(p_a, av) *
                                                   p_of_val(p_b_given_a[av], bv) ) / runs
                    p_ab_given_touchb[av][bv] += ( p_of_val(p_b, bv) *
                                                   p_of_val(p_a_given_b[bv], av) ) / runs
        logp_obs_given_cuts = logp_obs_given(cnt, p_ab_given_cuts)
        alternatives=[p_ab_given_toucha, p_ab_given_touchb]
        for (i,alternative) in enumerate(alternatives):
            logp_obs_given_alt = logp_obs_given(cnt, alternative)
            logbf = logp_obs_given_cuts - logp_obs_given_alt
            logbfs[i] += logbf
    return exp(min(logbfs))
