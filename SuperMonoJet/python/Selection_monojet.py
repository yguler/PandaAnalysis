from PandaCore.Tools.Misc import *
from re import sub

metTrigger='(trigger&1)!=0'
eleTrigger='(trigger&2)!=0'
phoTrigger='(trigger&4)!=0'

cuts = {}
weights = {}
triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
}

baseline = 'metFilter==1  && nJet>0 && nTau==0'
calocutSR = '(abs(calomet-pfmet)/pfUmag)<0.5'
calocutW = '(abs(calomet-pfmet)/pfUWmag)<0.5'
calocutZ = '(abs(calomet-pfmet)/pfUZmag)<0.5'
calocutA = '(abs(calomet-pfmet)/pfUAmag)<0.5'

cuts['SR0'] = tAND(baseline,'nLooseLep==0 && nLooseElectron==0 && nLoosePhoton==0 && jetNMBtags==0 && jet1Pt>100 && abs(jet1Eta)<2.5 && dphipfmet>0.5 && pfUmag>250 && calocutSR<0.5')

  #inclusive
cuts['ZmmINC'] = tAND(baseline,'nLooseMuon==2 && looseLep1IsTight==1 && diLepMass>60 && diLepMass<120')
cuts['ZeeINC'] = tAND(baseline,'nLooseElectron==2 && looseLep1IsTight==1 && diLepMass>60 && diLepMass<120')
cuts['WmnINC'] = tAND(baseline,'nLooseMuon==1 && looseLep1IsTight==1')
cuts['WenINC'] = tAND(baseline,'nLooseElectron==1 && looseLep1IsTight==1')

cuts['ZmmCR'] = tAND(cuts['ZmmINC'],tAND(calocutZ,'(nLooseElectron+nLoosePhoton+nTau)==0 && pfUZmag>250 && dphipfUZ>0.5 && jet1Pt>100 && abs(jet1Eta)<2.5'))
cuts['ZeeCR'] = tAND(cuts['ZeeINC'],tAND(calocutZ,'(nLooseMuon+nLoosePhoton+nTau)==0 && calocutZ<0.5 && pfUZmag>250 && dphipfUZ>0.5 && jet1Pt>100 && abs(jet1Eta)<2.5'))
cuts['WmnCR'] = tAND(cuts['WmnINC'],tAND(calocutZ,'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutW<0.5 && pfUWmag>250 && dphipfUW>0.5 && jet1Pt>100 && abs(jet1Eta)<2.5'))
cuts['WenCR'] = tAND(cuts['WenINC'],tAND(calocutZ,'(nLooseMuon+nLoosePhoton+nTau)==0 && calocutW<0.5 && pfUWmag>250 && dphipfUW>0.5 && jet1Pt>100 && abs(jet1Eta)<2.5'))

##weights for specific regions                                                                   
#base
weights['base'] = 'normalizedWeight*sf_pu*sf_ewkV*sf_qcdV'

#signal region in all b-tag categories
weights['SR0'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_tt*sf_Medbtag0')

#inclusive 
weights['ZmmINC'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt')
weights['ZeeINC'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt')
weights['WmnINC'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt')
weights['WenINC'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt')

#Control Region Specific
#0 b-tag
weights['ZmmCR'] = tTIMES(weights['ZmmINC'],'%f*sf_Medbtag0')
weights['ZeeCR'] = tTIMES(weights['ZeeINC'],'%f*sf_Medbtag0')
weights['WmnCR'] = tTIMES(weights['WmnINC'],'%f*sf_Medbtag0')
weights['WenCR'] = tTIMES(weights['WenINC'],'%f*sf_Medbtag0')

#b-tag shift
for r in ['SR0','ZmmCR','ZeeCR','WmnCR','WenCR']:
    for shift in ['BUp','BDown','MUp','MDown']:
        for cent in ['sf_Medbtag']:
            if 'Medbtag0' in weights[r]:
                weights[r+'_'+cent+shift] = sub(cent+'0',cent+'0'+shift,weights[r])

#scale and pdf shift
for r in ['SR0','ZmmCR','ZeeCR','WmnCR','WenCR']:
    for shift in ['ScaleUp','ScaleDown','PDFUp','PDFDown']:
        if shift == 'ScaleUp':
            weights[r+'_'+shift] = weights[r] + "*(scale[3]+1)"
        if shift == 'ScaleDown':
            weights[r+'_'+shift] = weights[r] + "*(scale[5]+1)"
        if shift == 'PDFUp':
            weights[r+'_'+shift] = weights[r] + "*pdfUp"
        if shift == 'PDFDown':
            weights[r+'_'+shift] = weights[r] + "*pdfDown"
