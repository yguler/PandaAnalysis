from PandaCore.Tools.Misc import *
from re import sub
Data_Trig='( (PD_SingleEle && HLT_EleOR) || (PD_SingleMu && (HLT_MuOR) && !(HLT_EleOR) )) && ( METfilters == 1)'
IsoVeto='!iso_Veto'

cuts = {

    'bjet0': 'nLep ==1 && Lep_pt > 25 && Selected == 1  &&  nVeto == 0 && nJets30Clean >= 1 && nBJet == 0 && MT>=160 && !iso_Veto',
    'bjet1': 'nLep ==1 && Lep_pt > 25 && Selected == 1  &&  nVeto == 0 && nJets30Clean >= 1 && nBJet == 1 && MT>=160 && !iso_Veto',
    'bjet2': 'nLep ==1 && Lep_pt > 25 && Selected == 1  &&  nVeto == 0 && nJets30Clean >= 1 && nBJet == 2 && MT>=160 && !iso_Veto',
    'nlep2bjetgt0': 'nLep ==2 && Lep_pt > 25 && Selected == 1  &&  nVeto == 0 && nJets30Clean >= 1 && nBJet > 0 && MT>=160 && !iso_Veto',
}

weights = {
    'bjet0': 'Xsec*1*btagSF*nISRttweight*puRatio*lepSF',
    'bjet1': 'Xsec*1*btagSF*nISRttweight*puRatio*lepSF',
    'bjet2': 'Xsec*1*btagSF*nISRttweight*puRatio*lepSF',
    'nlep2bjetgt0': 'Xsec*1*btagSF*nISRttweight*puRatio*lepSF',
}
