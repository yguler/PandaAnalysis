#!/usr/bin/env python

from os import getenv,path,popen
from PandaCore.Tools.job_management import *
import subprocess
import sys
import argparse
from glob import glob
from re import sub as rsub
import cPickle as pickle
from itertools import chain 

logdir = getenv('SUBMIT_LOGDIR')
workdir = getenv('SUBMIT_WORKDIR')
parser = argparse.ArgumentParser(description='analyze log files')
parser.add_argument('--verbose',action='store_true')
args = parser.parse_args()

def sub(x, y, z):
    return rsub(y, z, x)

def clean(x):
    x = x.replace('+', '\+')
    x = sub(x, '\n', '')
    x = sub(x, 'input_[0-9A-Z\-]*.root', 'input.root')
    x = sub(x, '[0-9][0-9][0-9][0-9]+', 'X')
    return x 

cmd = 'grep -i error $SUBMIT_LOGDIR/*err'
errors_by_file = {}
for l in subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.readlines():
    f = sub(l.split(':')[0].split('/')[-1], '.err', '')
    if f not in errors_by_file:
        errors_by_file[f] = []
    msg = clean(l.replace(l.split(':')[0], '')[1:])
    errors_by_file[f].append(msg)


aggregates = {}
for f,errors in errors_by_file.iteritems():
    for e in errors:
        if e not  in aggregates:
            aggregates[e] = set([f])
        else:
            aggregates[e].add(f)


correlations = []
sorted_aggregates = sorted(aggregates)
for i in xrange(len(sorted_aggregates)):
    a = sorted_aggregates[i]
    added = False
    for j in xrange(i):
        b = sorted_aggregates[j]
        if aggregates[a] == aggregates[b]:
            for c in correlations:
                if b in c:
                    c.add(a)
                    added = True
                    break
    if not added:
      correlations.append(set([a]))

cache = pickle.load(open(getenv('SUBMIT_WORKDIR') + '/submission.pkl', 'rb'))


for i,c in enumerate(correlations):
    print 'Failure class %i:'%i
    for msg in c:
        print '   ',msg

for i,c in enumerate(correlations):
    if args.verbose:
        print 'Failed with error class %i:'%i
    files = {}
    for f in sorted(aggregates[list(c)[0]]):
        cmd = 'grep -o "pandaf[^ ]*\.root" $SUBMIT_LOGDIR/%s.err'%f
        for l in subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.readlines():
            ll = l.strip()
            if ll not in files:
                files[ll] = 0
            files[ll] += 1
    if args.verbose:
        for f,n in files.iteritems():
            print '   ',f,n
    else:
        try:
            print 'Failure class %i failed on %3i files, an average of %.1f times'%(i, len(files), float(sum(files.values())) / len(files) / 2) # each file appears in error logs twice
        except ZeroDivisionError:
            pass
