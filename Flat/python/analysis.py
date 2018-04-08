#!/usr/bin/env python

from PandaCore.Tools.Load import Load 
from PandaCore.Tools.Misc import PInfo, PError
import ROOT as root

Load('PandaAnalyzer')

def _dump(a):
    PInfo('PandaAnalysis.Flat.analysis','Summary of analysis %s:'%(a.name))
    for k in dir(a):
        if k[0] == '_':
            continue
        if type(getattr(a, k)) != int:
            continue
        PInfo('PandaAnalysis.Flat.analysis','    %20s = %s'%(k, 'True' if getattr(a, k) else 'False'))



def _analysis(name, verbose, **kwargs):
    a = root.Analysis(name)
    for k,v in kwargs.iteritems():
        if not hasattr(a, k):
            PError('PandaAnalysis.Flat.analysis','Could not set property %s'%k)
            return None 
        setattr(a, k, bool(v))
    setattr(a, 'dump', lambda : _dump(a))
    if verbose:
        a.dump()
    return a

def analysis(name, **kwargs):
    return _analysis(name, verbose=True, **kwargs)


# predefined!
boosted = lambda v=False : _analysis(
        name = 'boosted',
        verbose = v,
        boosted = True,
        btagWeights = True,
        useCMVA = True,
        hfCounting = True,
        jetFlavorJets = True
    )

resolved = lambda v=False : _analysis(
        name = 'resolved',
        verbose = v,
        boson = True,
        bjetRegression = True,
        btagWeights = True,
        useCMVA = True,
        hfCounting = True,
        jetFlavorJets = True
    )

monojet = lambda v=False : _analysis(
        name = 'monojet',
        verbose = v,
	monojet = True,
        btagWeights = True,
        useCMVA = True,
        hfCounting = True,
        jetFlavorJets = True
    )

lepmonotop = lambda v=False : _analysis(
        name = 'lepmonotop',
        verbose = v,
        lepmonotop = True,
        recoil = True,
        fatjet = False,
        btagSFs = True,
        btagWeights = True,
        useCMVA = True,
        complicatedLeptons = True,
        hfCounting = True,
        reclusterGen = True,
        bjetRegression = False,
        varyJES = True,
        rerunJES = True,
        jetFlavorPartons = False,
        jetFlavorJets = True,
    )


vbf = lambda v=False : _analysis(
        name = 'vbf',
        verbose = v,
        vbf = True,
        fatjet = False,
        btagSFs = False,
        puppi_jets = False
    )

