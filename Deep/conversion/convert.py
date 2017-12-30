#!/usr/bin/env python

from sys import argv, exit
import numpy as np 
from os import getenv, system
from PandaCore.Tools.Misc import PInfo 
from glob import glob 

singletons = ['rawpt', 'pt', 'msd', 'eta']
events = ['eventNumber']
fractions = {'train':0.7, 'test':0.15}
fcfg = open(argv[1])
name = '_'.join(argv[1].split('/')[-2:]).replace('.txt','')
outdir = getenv('SUBMIT_NPY')
datadir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/deep/'
me = argv[0]
argv = []

import ROOT as root 
f_pt = root.TFile.Open(datadir + 'flatten.root')
h_pt = f_pt.Get('h_%s'%(name.split('_')[0]))
f_pt_scaled = root.TFile.Open(datadir + 'flatten_scaled.root')
h_pt_scaled = f_pt_scaled.Get('h_%s'%(name.split('_')[0]))

data = {}
for fpath in fcfg.readlines():
    d = np.load(fpath.strip())
    for k,v in d.iteritems():
        if v.shape[0]:
            if k not in data:
                data[k] = []
            data[k].append(v)

if not len(data):
    PInfo(me, 'This was an empty config!')
    exit(0)

for k,v in data.iteritems():
    data[k] = np.concatenate(v)

def reweight(x_pt):
    return h_pt.GetBinContent(h_pt.FindBin(x_pt))
reweight = np.vectorize(reweight)

def reweight_s(x_pt):
    return h_pt_scaled.GetBinContent(h_pt_scaled.FindBin(x_pt))
reweight_s = np.vectorize(reweight_s)


data['ptweight'] = reweight(data['pt'])
data['ptweight_scaled'] = reweight_s(data['pt'])


def dump(idx, partition):
    outpath = 'tmp/' + partition + '/' + name + '_%s.npy'

    # singletons
    d = np.vstack([data[x][idx] for x in singletons]).T 
    np.save(outpath%'singletons', d)

    # events
    d = np.vstack([data[x][idx] for x in events]).T 
    np.save(outpath%'events', d)

    # pf
    d = data['pf'][idx, :, :]
    np.save(outpath%'pf', d)

    # pt weights
    d = [data['ptweight'][idx], data['ptweight_scaled'][idx]]
    d = np.vstack(d).T
    np.save(outpath%'ptweight', d)
    

indices = range(data['eventNumber'].shape[0])
np.random.shuffle(indices)

N = {k:int(len(indices) * v) for k,v in fractions.iteritems()}

for d in ['train', 'test', 'validate']:
    system('mkdir -p tmp/%s'%d)

dump(indices[:N['train']], 'train')
dump(indices[N['train']:N['train']+N['test']], 'test')
dump(indices[N['train']+N['test']:], 'validate')

for d in ['train', 'test', 'validate']:
    for ftmp in glob('tmp/'+d+'/*npy'):
        cmd = 'cp -v %s %s/%s'%(ftmp,outdir,ftmp.replace('tmp/',''))
        # cmd = 'gfal-copy -f file://$PWD/%s srm://t3serv006.mit.edu:8443/srm/v2/server?SFN=%s/%s'%(ftmp,outdir,ftmp.replace('tmp/',''))
        PInfo(me, cmd)
        system(cmd)

# system('rm -rf tmp')
