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
        srr=re.match('.*/[A-Z]*([0-9]*).results',fn).group(1)
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
        (srr,s)= re.match('[0-9x]* [A-Z]*([0-9]*) ([A-z]*)\n',line).groups()
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

def normmicond(a,b,c):
    eac=entropy(map(lambda x:'%s/%s'%x,zip(a,c)))
    ebc=entropy(map(lambda x:'%s/%s'%x,zip(b,c)))
    ec=entropy(c)
    ej=entropy(map(lambda x:'%s/%s/%s'%x,zip(a,b,c)))
    return (eac+ebc-ej-ec)/max(eac,ebc)

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
for i in range(10000):
    ri1=int(random.random()*len(gspes))
    ri2=int(random.random()*len(gspes))
    l1=copy.copy(gcnt[gspes[ri1]])
    l2=copy.copy(gcnt[gspes[ri2]])
    random.shuffle(l1)
    random.shuffle(l2)
    nulldist.append(normmi(l1,l2))

nulldist.sort()

thresh=nulldist[9900]

nullcdist=[]
for i in range(10000):
    ri1=int(random.random()*len(gspes))
    ri2=int(random.random()*len(gspes))
    ri3=int(random.random()*len(gspes))
    l1=copy.copy(gcnt[gspes[ri1]])
    l2=copy.copy(gcnt[gspes[ri2]])
    l3=copy.copy(gcnt[gspes[ri3]])
    random.shuffle(l1)
    random.shuffle(l2)
    random.shuffle(l3)
    nullcdist.append(normmicond(l1,l2,l3))

nullcdist.sort()

condthresh=nullcdist[9900]

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

#thresh=.116535
#thresh=.1
#condthresh=.067942

print "Using threshold %f and conditional threshold of %f"%(thresh,condthresh)

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

print "We have %d trios"%len(trios)

def abbr(name):
    [g,s]=name.split(' ')
    return g[0]+'.'+s[0:5]

def printres(a,h1,h2,micabh1,micabh2,mich1ba,mich2ba,acut,h1cut,h2cut):
    aab=abbr(a)
    h1ab=abbr(h1)
    h2ab=abbr(h2)
    print "Evidence that %s(%s) causes Crohn's disease:" % (a,aab)
    print "  Helper bacteria: %s(%s) and %s(%s)" % (h1,h1ab,h2,h2ab)
    print "  Evidence that there's some link:"
    print "    NMI(%s,%s) = %f   p<%.4f" % (aab,h1ab,mis[a][h1],findp(mis[a][h1],nulldist))
    print "    NMI(%s,%s) = %f   p<%.4f" % (aab,h2ab,mis[a][h2],findp(mis[a][h2],nulldist))
    print "    NMI(%s,%s) = %f   p<%.4f" % (h1ab,h2ab,mis[h1][h2],findp(mis[h1][h2],nulldist))
    print "    NMI(%s,Crohn's) = %f   p<%.4f" % (aab,mis[a]['sick'],findp(mis[a]['sick'],nulldist))
    print "    NMI(%s,Crohn's) = %f   p<%.4f" % (h1ab,mis[h1]['sick'],findp(mis[h1]['sick'],nulldist))
    print "    NMI(%s,Crohn's) = %f   p<%.4f" % (h2ab,mis[h2]['sick'],findp(mis[h2]['sick'],nulldist))
    print "  Evidence that some common factor causes all three bacteria:"
    print "    NMI(%s,%s|%s) = %f   p<%.4f" % (aab,h1ab,h2ab,h2cut,findp(h2cut,nullcdist))
    print "    NMI(%s,%s|%s) = %f   p<%.4f" % (aab,h2ab,h1ab,h1cut,findp(h1cut,nullcdist))
    print "    NMI(%s,%s|%s) = %f   p<%.4f" % (h1ab,h2ab,aab,acut,findp(acut,nullcdist))
    print "  Evidence that %s screens the other bacteria from Crohn's disease:" % a
    print "    NMI(%s,Crohn's|%s) = %f   p>%.4f" % (h1ab,aab,micabh1,findp(micabh1,nullcdist))
    print "    NMI(%s,Crohn's|%s) = %f   p>%.4f" % (h2ab,aab,micabh2,findp(micabh2,nullcdist))
    print "  Sanity check, the other bacteria do not screen %s from Crohn's disease:" % a
    print "    NMI(%s,Crohn's|%s) = %f   p<%.4f" % (aab,h1ab,mich1ba,findp(mich1ba,nullcdist))
    print "    NMI(%s,Crohn's|%s) = %f   p<%.4f" % (aab,h2ab,mich2ba,findp(mich2ba,nullcdist))
    print ""


for trio in trios:
    mic012=normmicond(gcnt[trio[0]],gcnt[trio[1]],gcnt[trio[2]])
    if mic012<condthresh:
        continue
    mic201=normmicond(gcnt[trio[2]],gcnt[trio[0]],gcnt[trio[1]])
    if mic201<condthresh:
        continue
    mic120=normmicond(gcnt[trio[1]],gcnt[trio[2]],gcnt[trio[0]])
    if mic120<condthresh:
        continue
    mic0b1=normmicond(gcnt[trio[1]],gcnt['sick'],gcnt[trio[0]])
    mic1b0=normmicond(gcnt[trio[0]],gcnt['sick'],gcnt[trio[1]])
    mic2b1=normmicond(gcnt[trio[1]],gcnt['sick'],gcnt[trio[2]])
    mic0b2=normmicond(gcnt[trio[2]],gcnt['sick'],gcnt[trio[0]])
    mic1b2=normmicond(gcnt[trio[2]],gcnt['sick'],gcnt[trio[1]])
    mic2b0=normmicond(gcnt[trio[0]],gcnt['sick'],gcnt[trio[2]])
    if mic0b1<condthresh and mic0b2<condthresh and mic1b0>condthresh and mic2b0>condthresh:
        printres(trio[0],trio[1],trio[2],mic0b1,mic0b2,mic1b0,mic2b0,mic120,mic201,mic012)
    if mic1b2<condthresh and mic1b0<condthresh and mic0b1>condthresh and mic2b1>condthresh:
        printres(trio[1],trio[0],trio[2],mic1b0,mic0b2,mic0b1,mic2b1,mic201,mic120,mic012)
    if mic2b1<condthresh and mic2b0<condthresh and mic0b2>condthresh and mic1b2>condthresh:
        printres(trio[2],trio[0],trio[1],mic2b0,mic2b1,mic0b2,mic1b2,mic012,mic120,mic201)
