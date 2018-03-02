#!/usr/bin/python

from datareader import Data
from utils import *
from scipy.stats import chi2_contingency
from scipy.stats import chi2
from scipy.stats import chisquare
from math import log
from sys import stdout, argv
import direction, severs
from common_cause import has_common_cause


data = Data()
data.read_from_files()

pl = data.get_pid_list(lambda(p): (p.health in ['control', 'ileal CD'] and
                                   len(p.fractions)>0))

sick = data.get_data(pl, 'health', bucketizer=lambda(v):v=='ileal CD')

names = argv[1:]
boolvals=[]

for species in names:
    if species not in data.bacteria:
        exit('Unrecognized species "%s"' % species)
    
    vals = data.get_data(pl, 'fractions', species)
    co = findcutoff(vals, sick)
    boolvals.append([(i>co.threshold)==(co.sick_when_more) for i in vals])

bf = direction.montecarlo(boolvals[0], sick, boolvals[1], len(sick)).bayes_fwd_rev
