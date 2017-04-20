#!/usr/bin/python
from collections import defaultdict
import glob
import re
from math import log
import random
import copy
from sys import argv
from libmhchi import *
from numpy import float64
from scipy.stats import chi2_contingency

cnt=defaultdict(lambda:defaultdict(lambda:0))

lines=file('../study2_genome.txt').readlines()
for line in lines:
    if line[0] in ['#', 'd', '\n']:
        continue
    w=line.split('\t')
    cnt['NOD2'][w[0]]=(int(w[2])>0)

lines=file('../study2_sick.txt').readlines()
for line in lines:
    if line[0] in ['#', 'd', '\n']:
        continue
    w=line.split('\t')
    if w[3] in ['control\n', 'ileal CD\n']:
        cnt['sick'][w[0]]=w[3][:-1]

patients=[]
for i in cnt['NOD2']:
    if i in cnt['sick']:
        patients.append(i)

spes=set()
tots=defaultdict(lambda:0)

gsi=file('study2/gsi')
for gsiline in gsi:
    [sample, patient]=gsiline[:-1].split(' ')
    if patient not in patients:
        continue
    sf=file('study2/%s.results'%sample)
    for line in sf:
        (n,spe)=re.match(" *([0-9]*) (.*)\n",line).groups()
        n=int(n)
        if n>1:
            cnt[spe][patient]+=n
            spes.add(spe)
            tots[patient]+=n

ncats=5

for spe in spes:
    for patient in patients:
        if cnt[spe][patient]==0:
            cnt[spe][patient]=ncats
        else:
            cnt[spe][patient]=min(int(round(-log(float(cnt[spe][patient])/tots[patient],10))),ncats-1)

#for patient in patients:
#    print '%s, %s, %s'% (cnt['NOD2'][patient], cnt['Desulfitobacterium dehalogenans'][patient], cnt['sick'][patient])

def d2l(di,idx):
    o=[]
    for i in idx:
        o.append(di[i])
    return o

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


print 'Original Species: %d<br>' % len(spes)

newspes=[]
for spe in spes:
    if entropy(d2l(cnt[spe],patients)) > .5:
        newspes.append(spe)
spes=newspes
spes.sort()

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

nulldist=[]
for i in range(10000):
    l1=copy.copy(d2l(cnt['NOD2'],patients))
    l2=copy.copy(d2l(cnt['sick'],patients))
    rspe = random.choice(spes)
    l3=copy.copy(d2l(cnt[rspe], patients))
    random.shuffle(l1)
    random.shuffle(l2)
    random.shuffle(l3)
    nulldist.append(normmicond(l1,l2,l3))
nulldist.sort()
#print nulldist

# for patient in patients:
#     tmp = 0
#     tmp += (cnt['Anoxybacillus flavithermus'][patient] < ncats-1)
#     tmp += (cnt['Anoxybacillus kamchatkensis'][patient] < ncats)
#     tmp += (cnt['Anoxybacillus thermarum'][patient] < ncats)
#     tmp += (cnt['Desulfitobacterium dehalogenans'][patient] < ncats)
#     tmp += (cnt['Flavobacterium thermophilum'][patient] < ncats-1)
#     tmp += (cnt['Gemella haemolysans'][patient] < ncats)
#     tmp += (cnt['Granulicatella elegans'][patient] < ncats)
#     tmp += (cnt['Roseiflexus sp.'][patient] < ncats)
#     tmp += (cnt['Sphaerobacter thermophilus'][patient] < ncats)
#     tmp += (cnt['Sphingomonas leidyi'][patient] < ncats)
#     tmp += (cnt['Sphingomonas koreensis'][patient] < ncats)
#     tmp += (cnt['Veillonella dispar'][patient] < ncats)
#     tmp += (cnt['Veillonella parvula'][patient] < ncats)
#     #cnt['composite'][patient] = tmp
#     if tmp>6:
#         cnt['composite'][patient] = ncats-1
#     else:
#         cnt['composite'][patient] = ncats
    
# spes.append('composite')

#compcnt = [[0]*16, [0]*16]
#for patient in patients:
#     compcnt[cnt['sick'][patient]=='ileal CD'][cnt['composite'][patient]] += 1
#print compcnt
#exit()

def describe(bactcrohn):
    harmful = (bactcrohn[1][ncats] >= bactcrohn[0][ncats])
    i=ncats-1
    changes=-1
    while i>0:
        if bactcrohn[1][i]==0 and bactcrohn[0][i]==1:
            break
        if changes==-1 and (bactcrohn[1][i] >= bactcrohn[0][i])!=harmful:
            changes = i
        elif changes!=-1 and (bactcrohn[1][i] >= bactcrohn[0][i])==harmful:
            return 'complicated'
        i-=1
    if changes==-1:
        return 'error'
    if changes==ncats-1:
        if harmful:
            return 'absent'
        else:
            return 'present'
    if harmful:
        return '&lt; 10<sup>-%d</sup>' % (changes+1)
    else:
        return '&gt; 10<sup>-%d</sup>' % (changes+1)

savedcnts={}
print 'Total species: %d<br>' % (len(spes)-1)
nod2l=d2l(cnt['NOD2'],patients)
sickl=d2l(cnt['sick'],patients)
print 'Entropy(NOD2): %f<br>' % entropy(nod2l)
print 'Entropy(Crohn): %f<br>' % entropy(sickl)
print '<table border>'
print '<tr><td>Name</td><td>cmi</td><td>p(cmi|h0)</td><td>entropy</td><td>Associated with Disease</td><td>p(o|nnaa)</td><td>p(o|nl)</td><td>p(o|dl)</td><td>p(o|nm)</td><td>p(o|am)</td><td>p(o|m)</td></tr>'
for spe in spes:
    spel=d2l(cnt[spe],patients)
    cmi = normmicond(nod2l, sickl, spel)
    p = findp(cmi,nulldist)
#    if p>.05:
    if True:
        spenulldist=[]
        bactcrohn = [defaultdict(lambda:0), defaultdict(lambda:0)]
        for crohn,bact in zip(sickl,spel):
            bactcrohn[crohn=='ileal CD'][bact] += 1
        specnt=float64([[[0, 0], [0, 0]], [[0, 0], [0, 0]]]) # float64 so 0/0=nan instead of exception
        jsncnt=float64([[0, 0], [0, 0]])
        for nod2,crohn,bact in zip(nod2l,sickl,spel):
            bact_bad_sign = (bactcrohn[1][bact] >= bactcrohn[0][bact])
            specnt[bact_bad_sign][nod2][crohn=='ileal CD']+=1
            jsncnt[bact_bad_sign][nod2]+=1
        if spe in ['Bacteroides nordii', 'Sphingomonas koreensis', 'Veillonella dispar']:
            savedcnts[spe] = specnt
        try:
            ps=odds(specnt)
        except ZeroDivisionError:
            print '<td colspan=9>cannot simulate %s</td></tr>' % spe
            continue
        try:
            pnnaa=chi2_contingency(jsncnt)[1]
        except ValueError:
            pnnaa=float('nan')
        pnl=sim(simnolink,ps)[1][1]
        pdl=sim(simdlinkonly,ps)[1][1]
        pnm=sim(simsep,ps)[1][1]
        pam=sim(simantimed,ps)[1][1]
        pm=sim(simmed,ps)[1][1]
        #        if max(pnl, pdl, pnm, pam) < 1e-1
        #if True:
        if pnnaa < .05:
            print '<tr><td>%s</td><td>%f</td><td>%f</td><td>%.2g</td>' % (spe,cmi,p,entropy(spel))
            print '<td>%s</td>' % describe(bactcrohn)
            for v in [pnnaa, pnl, pdl, pnm, pam, pm]:
                if v<1e-4:
                    print '<td bgcolor=yellow>%.2g</td>' % v
                else:
                    print '<td>%.2g</td>' % v
            print '</tr>'

print '</table>'

#printall('Veillonella dispar', savedcnts['Veillonella dispar'])
#printall('Bacteroides nordii', savedcnts['Bacteroides nordii'])
#printall('Sphingomonas koreensis', savedcnts['Sphingomonas koreensis'])
