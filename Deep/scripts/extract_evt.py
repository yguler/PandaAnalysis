#!/usr/bin/env python

import numpy as np
from glob import glob
from sys import argv, exit
from PandaCore.Tools.Misc import PInfo, PWarning

flist = glob(argv[2])
target_evt = np.int64(argv[1])

for f in flist:
    arr = np.load(f)
    evt = arr['eventNumber']
    mask = (evt == target_evt)
    if np.sum(mask):
        idx = np.argmax(mask)
        pf = arr['pf'][idx]
        msd = arr['msd'][idx]
        pt = arr['pt'][idx]
        found = evt[idx]
        np.savez('sliced.npz', pf=pf, msd=msd, pt=pt, evt=found)
        PInfo(argv[0], 'Found %i in %s'%(target_evt, f))
        exit(0)


PWarning(argv[0], 'Could not find %i in %s'%(target_evt, argv[2]))
exit(1)
