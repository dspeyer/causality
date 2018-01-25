#!/usr/bin/python

from utils import count, p_of_val, struct, acclarr, mularr, sumarr
from math import log,exp
from scipy.stats import chi2_contingency, chi2, chisquare
from numpy.random import beta

def chi2_dir(cause, effect, unknown, n, p_cause, p_effect_given_cause):
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

def get_joints_by_model(p):
    p.cause_unknown_chain=[[0,0],[0,0]]
    p.cause_unknown_collide=[[0,0],[0,0]]
    for c in range(2):
        for u in range(2):
            for e in range(2):
                p.cause_unknown_chain[c][u] += (
                    p_of_val(p.cause, c) *
                    p_of_val(p.effect_given_cause[c], e) *
                    p_of_val(p.unknown_given_effect[e], u))
            p.cause_unknown_collide[c][u] = (
                p_of_val(p.cause, c) *
                p_of_val(p.unknown, u))
    return p

def get_factor(p, cnt):
    logp_obs_given_chain = 0
    logp_obs_given_collide = 0
    for c in range(2):
        for u in range(2):
            logp_obs_given_chain += log(p.cause_unknown_chain[c][u])*cnt[c][u]
            logp_obs_given_collide += log(p.cause_unknown_collide[c][u])*cnt[c][u]
    return exp(logp_obs_given_collide - logp_obs_given_chain)


def mle(cause, effect, unknown, n, p_cause, p_effect_given_cause):
    p=struct(cause=p_cause, effect_given_cause=p_effect_given_cause)
    cnt = count(zip(effect, unknown))
    chi_indep = chi2_contingency(cnt)
    p.unknown_given_effect = [ float(cnt[0][1]) / sum(cnt[0]),
                               float(cnt[1][1]) / sum(cnt[1]) ]
    cnt = count(zip(unknown))
    p.unknown = float(cnt[1]) / sum(cnt)
    p = get_joints_by_model(p)
    cnt = count(zip(cause, unknown))
    bayes_factor = get_factor(p, cnt)
    return struct(reject_indep=chi_indep,
                  bayes_fwd_rev=bayes_factor)

#def beta(a,b):
#    return 1 - (float(a) / (a+b))

def montecarlo(cause, effect, unknown, n, *ignore):
    cnt_cause = count(zip(cause))
    cnt_unknown = count(zip(unknown))
    cnt_cause_effect = count(zip(cause, effect))
    cnt_effect_unknown = count(zip(effect, unknown))
    sumarr(cnt_cause_effect, 0.1) # make beta dist work with zeros
    sumarr(cnt_cause, 0.1)
    sumarr(cnt_unknown, 0.1)
    sumarr(cnt_effect_unknown, 0.1)
    cnt_cause_unknown = count(zip(cause, unknown))
    rounds = 500
    p_overall = struct(cause_unknown_chain=[[0,0],[0,0]],
                       cause_unknown_collide=[[0,0],[0,0]])
    for i in range(rounds):
        p=struct()
        p.cause = 1-beta(*cnt_cause)
        p.unknown = 1-beta(*cnt_unknown)
        p.effect_given_cause = [1-beta(*cnts) for cnts in cnt_cause_effect]
        p.unknown_given_effect = [1-beta(*cnts) for cnts in cnt_effect_unknown]
        p = get_joints_by_model(p)
        acclarr(p_overall.cause_unknown_chain, p.cause_unknown_chain)
        acclarr(p_overall.cause_unknown_collide, p.cause_unknown_collide)
    mularr(p_overall.cause_unknown_chain, 1.0/rounds)
    mularr(p_overall.cause_unknown_collide, 1.0/rounds)
    try:
        bayes_factor = get_factor(p_overall, cnt_cause_unknown)
    except ValueError:
        print '==ValueError=='
        print p_overall.__dict__
        raise ValueError()
    return struct(bayes_fwd_rev=bayes_factor)
    
