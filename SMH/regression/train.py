#!/usr/bin/env python

from sys import argv
necf = int(argv[1])
nextra = int(argv[2])
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
from os import getenv, system

import cfgMinimal as cfg
#import cfgAll as cfg

Load('TMVATrainer')

workdir = getenv('PANDA_FLATDIR')
system('mkdir -p %s/training'%workdir)


trainer = root.TMVATrainer('bjet_regression',workdir+'training/')
trainer.treename = 'events'
trainer.sigweight = 'normalizedWeight'
sanitycut='hbbm>0 && jetGenFlavor[hbbjtidx0]==5'
trainer.sigcut =sanitycut

variables = [
    'jetPt[hbbjtidx0]',
    'jetEta[hbbjtidx0]',
    'jetE[hbbjtidx0]',
    'jetEMFrac[hbbjtidx0]',
    'jetHadFrac[hbbjtidx0]',
    'jetNLep[hbbjtidx0]',
    'jetLeadingLepPt[hbbjtidx0]',
    'jetLeadingTrkPt[hbbjtidx0]',
    ]

for v in variables:
  trainer.AddVariable(v, v)


trainer.AddRegressionFile(workdir+'/ggHbb.root')
trainer.BookMLP()
trainer.TrainAll()

