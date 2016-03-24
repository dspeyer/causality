#!/usr/bin/python

from utils import struct
from collections import defaultdict
import re


class Data:
    def __init__(self):
        self.patients={}
        self.bacteria={}

    def read_from_files(self, 
                        sick_fn='../study2_sick.txt',
                        gene_fn='../study2_genome.txt',
                        sample_fn='study2/gsi',
                        sex_fn='study2/sex',
                        bact_dir='study2/'):
        for fn in [sick_fn, gene_fn]:
            for line in file(fn):
                line = line.strip()
                if line=='' or line[0]=='#':
                    continue
                words = line.split('\t')
                if words[0]=='dbGaP SubjID':
                    continue
                if fn==sick_fn:
                    self.patients[words[0]] = struct(id=words[0], health=words[3])
                else:
                    self.patients[words[0]].nod2 = int(words[2])
                    self.patients[words[0]].atg16l1 = int(words[3])
        pid_by_sample = {}
        for line in file(sample_fn):
            line=line.strip()
            [sample, pid]=line.split(' ')
            if 'samples' not in self.patients[pid]:
                self.patients[pid].samples=[]
            self.patients[pid].samples.append(sample)
            pid_by_sample[sample]=pid
        for line in file(sex_fn):
            line=line.strip()
            [x, sample, sex] = line.split(' ')
            patient = self.patients[pid_by_sample[sample]]
            if 'sex' in patient and patient.sex!=sex:
                print 'WARNING: inconsistent sex for %s at %s' % (patient.id, sample)
            if sex!='unknown':
                patient.sex = sex
        for pid in self.patients:
            patient = self.patients[pid]
            if 'samples' not in patient:
                continue
            patient.raw_counts = defaultdict(lambda:0)
            patient.total_count = 0
            for sid in patient.samples:
                for line in file(bact_dir+sid+'.results'):
                    try:
                        (n,species)=re.match(" *([0-9]*) (.*)\n",line).groups()
                    except:
                        print 'ERROR: '+line
                    n=int(n)
                    patient.raw_counts[species] += n
                    patient.total_count += n
                    self.bacteria[species] = True
            patient.fractions = defaultdict(lambda:0.0)
            for species in patient.raw_counts:
                if species > 1:
                    patient.fractions[species] = float(patient.raw_counts[species]) / patient.total_count

    def get_pid_list(self, filter=(lambda(x):True)):
        out=[]
        for pid in self.patients:
            try:
                if filter(self.patients[pid]):
                    out.append(pid)
            except AttributeError:
                pass
        out.sort()
        return out

    def get_data(self, pidlist, attr, key=None, bucketizer=lambda(x):x):
        out=[]
        for pid in pidlist:
            try:
                val = getattr(self.patients[pid], attr)
                if key:
                    val = val[key]
                out.append(bucketizer(val))
            except:
                out.append(None)
        return out
