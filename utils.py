#!/usr/bin/python

from copy import copy

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


def count(data, limits):
    out=0
    for limit in reversed(limits):
        out=[out]
        for i in range(limit-1):
            out.append(copy(out[0]))
    for row in data:
        tmp=out
        for i,col in enumerate(row):
            if (i<len(row)-1):
                tmp=tmp[col]
            else:
                tmp[col]+=1
    return out
