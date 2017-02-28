#!/usr/bin/python

from utils import count, p_of_val, struct
from math import log,exp
from scipy.stats import chi2_contingency, chi2

def chi2(cause, effect, unknown, n, p_cause, p_effect_given_cause):
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

def mle(cause, effect, unknown, n, p_cause, p_effect_given_cause):
    cnt = count(zip(effect, unknown))
    chi_indep = chi2_contingency(cnt)
    p_unknown_given_effect = [ float(cnt[0][1]) / sum(cnt[0]),
                               float(cnt[1][1]) / sum(cnt[1]) ]
    cnt = count(zip(unknown))
    p_unknown = float(cnt[1]) / sum(cnt)
    p_cause_unknown_chain=[[0,0],[0,0]]
    p_cause_unknown_collide=[[0,0],[0,0]]
    for c in range(2):
        for u in range(2):
            for e in range(2):
                p_cause_unknown_chain[c][u] += (
                    p_of_val(p_cause, c) *
                    p_of_val(p_effect_given_cause[c], e) *
                    p_of_val(p_unknown_given_effect[e], u))
            p_cause_unknown_collide[c][u] = (
                p_of_val(p_cause, c) *
                p_of_val(p_unknown, u))
    cnt = count(zip(cause, unknown))
    logp_obs_given_chain = 0
    logp_obs_given_collide = 0
    for c in range(2):
        for u in range(2):
            logp_obs_given_chain += log(p_cause_unknown_chain[c][u])*cnt[c][u]
            logp_obs_given_collide += log(p_cause_unknown_collide[c][u])*cnt[c][u]
    bayes_factor = exp(logp_obs_given_collide - logp_obs_given_chain)
    return struct(reject_indep=chi_indep,
                  bayes_fwd_rev=bayes_factor)
