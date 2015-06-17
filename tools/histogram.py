#!/usr/bin/env python

import numpy
import sys

style=1

if len(sys.argv)!=3:
    print("Usage:{0} <input data> <nbins>".format(sys.argv[0]))

nbins=int(sys.argv[2])


data=[]

for l in open(sys.argv[1]):
    data.append(float(l))

hist,edges=numpy.histogram(a=data,bins=nbins)

if style==2:
    print("read serr 1")
for i in range(0,len(hist)):
    if style==1:
        print("{0} {1}".format(edges[i],0))
        print("{0} {1}".format(edges[i],hist[i]))
        print("{0} {1}".format(edges[i+1],hist[i]))
        print("{0} {1}".format(edges[i+1],0))
    elif style==2:
        x=(edges[i+1]+edges[i])/2
        xe=(edges[i+1]-edges[i])/2
        print("{0} {1} {2}".format(x,xe,hist[i]))
    

