#!/usr/bin/python
from collections import defaultdict
import glob
import re
from math import log
import random
import copy
from sys import argv

cnt=defaultdict(lambda:defaultdict(lambda:0))
srrs=set()
spes=set()

tots=defaultdict(lambda:0)

g=glob.glob(argv[1])
for fn in g:
    try:
        srr=re.match('.*/SRR([0-9]*).results',fn).group(1)
    except:
        print 'ERROR: '+fn
        continue
    srrs.add(srr)
    for line in file(fn):
        try:
            (n,spe)=re.match(" *([0-9]*) (.*)\n",line).groups()
        except:
            print 'ERROR: '+line
        n=int(n)
        if n>1:
            cnt[spe][srr]=n
            spes.add(spe)
        tots[srr]+=n

for spe in spes:
    for srr in srrs:
        old=cnt[spe][srr]
        if cnt[spe][srr]==0:
            cnt[spe][srr]='absent'
        else:
            cnt[spe][srr]='rareness'+str(int(round(-log(float(cnt[spe][srr])/tots[srr],10))))

srrs=list(srrs)
spes=list(spes)
srrs.sort()
spes.sort()

meta=file(argv[2])
for line in meta:
    try:
        (srr,s)= re.match('[0-9x]* SRR([0-9]*) ([A-z]*)\n',line).groups()
        cnt['sick'][srr]=s
    except:
        pass

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


def normmi(a,b):
    ea=entropy(a)
    eb=entropy(b)
    ej=entropy(map(lambda x:'%s/%s'%x,zip(a,b)))
    return (ea+eb-ej)/max(ea,eb)

def normmibact(a,b,shuffle=False):
    al=d2l(cnt[a],srrs)
    bl=d2l(cnt[b],srrs)
    if shuffle:
        random.shuffle(al)
        random.shuffle(bl)
    return normmi(al,bl)

gcnt={}
for spe in spes:
    l=d2l(cnt[spe],srrs)
    if entropy(l)>.5:
        gcnt[spe]=l

gcnt['sick']=d2l(cnt['sick'],srrs)

gspes=gcnt.keys()


mis=defaultdict(dict)
for i in range(len(gspes)):
    si=gspes[i]
    for j in range(i+1,len(gspes)):
        sj=gspes[j]
        mis[si][sj]=normmi(gcnt[si],gcnt[sj])
        mis[sj][si]=mis[si][sj]

nulldist=[]
for i in range(1000):
    ri1=int(random.random()*len(gspes))
    ri2=int(random.random()*len(gspes))
    l1=copy.copy(gcnt[gspes[ri1]])
    l2=copy.copy(gcnt[gspes[ri2]])
    random.shuffle(l1)
    random.shuffle(l2)
    nulldist.append(normmi(l1,l2))

nulldist.sort()

thresh=nulldist[990]

print "Using threshold %f"%thresh

trios=[]
for i in range(len(gspes)):
    si=gspes[i]
    if si=='sick':
        continue
    for j in range(i+1,len(gspes)):
        sj=gspes[j]
        if sj=='sick':
            continue
        for k in range(j+1,len(gspes)):
            sk=gspes[k]
            if sk=='sick':
                continue
            if  (mis[si][sj]>thresh and mis[si][sk]>thresh and mis[sj][sk]>thresh and
                 mis[si]['sick']>thresh and mis[sj]['sick']>thresh and mis[sk]['sick']>thresh):
                trios.append([si,sj,sk])

def normmicond(a,b,c):
    eac=entropy(map(lambda x:'%s/%s'%x,zip(a,c)))
    ebc=entropy(map(lambda x:'%s/%s'%x,zip(b,c)))
    ec=entropy(c)
    ej=entropy(map(lambda x:'%s/%s/%s'%x,zip(a,b,c)))
    return (eac+ebc-ej-ec)/max(eac,ebc)

for trio in trios:
    if normmicond(gcnt[trio[0]],gcnt[trio[1]],gcnt[trio[2]])<thresh:
        continue
    if normmicond(gcnt[trio[0]],gcnt[trio[2]],gcnt[trio[1]])<thresh:
        continue
    if normmicond(gcnt[trio[1]],gcnt[trio[2]],gcnt[trio[0]])<thresh:
        continue
    mic0b1=normmicond(gcnt[trio[1]],gcnt['sick'],gcnt[trio[0]])
    mic1b0=normmicond(gcnt[trio[0]],gcnt['sick'],gcnt[trio[1]])
    mic2b1=normmicond(gcnt[trio[1]],gcnt['sick'],gcnt[trio[2]])
    mic0b2=normmicond(gcnt[trio[2]],gcnt['sick'],gcnt[trio[0]])
    mic1b2=normmicond(gcnt[trio[2]],gcnt['sick'],gcnt[trio[1]])
    mic2b0=normmicond(gcnt[trio[0]],gcnt['sick'],gcnt[trio[2]])
    if mic0b1<thresh and mic0b2<thresh and mic1b0>thresh and mic2b0>thresh:
        print "%s cuts %s and %s" % (trio[0],trio[1],trio[2])
    if mic1b2<thresh and mic1b0<thresh and mic0b1>thresh and mic2b1>thresh:
        print "%s cuts %s and %s" % (trio[1],trio[2],trio[0])
    if mic2b1<thresh and mic2b0<thresh and mic0b2>thresh and mic1b2>thresh:
        print "%s cuts %s and %s" % (trio[2],trio[1],trio[0])

