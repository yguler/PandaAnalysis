from PandaCore.Tools.Misc import *
from re import sub

metTrigger='(trigger&1)!=0'
eleTrigger='(trigger&2)!=0'

cuts = {}
weights = {}
triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
}

baseline = 'metFilter==1 && nJet>0 && jet1Pt>30 && nTau==0'
cat1 = '( ( nJet==1 && jet1CSV>0.800 ) || ( nJet==2 && ( ( jet1CSV>0.800 ) + ( jet2CSV>0.800 ) )==1 ) )'
cat2 = '( jet1Pt>50 && jet2Pt>50 && jet3Pt>30 && ( ( nJet==2 && ( ( jet1CSV>0.800 ) + ( jet2CSV>0.800 ) )==2 ) || ( nJet==3 && ( ( jet1CSV>0.800 ) + ( jet2CSV>0.800 ) + ( jet3CSV>0.800 ) )==2  ) ) )'

calocutSR = '(abs(calomet-pfmet)/pfUmag)'
calocutW = '(abs(calomet-pfmet)/pfUWmag)'
calocutZ = '(abs(calomet-pfmet)/pfUZmag)'
calocutA = '(abs(calomet-pfmet)/pfUAmag)'

cuts['SR0'] = tAND(baseline,'nLooseLep==0 && nLooseElectron==0 && nLoosePhoton==0 && jetNMBtags==0 && jet1Pt>50 && abs(jet1Eta)<2.5 && jet2Pt>30 && abs(jet2Eta)<2.5 && dphipfmet>0.5 && pfUmag>250 && calocutSR<0.5')
cuts['SR1'] = tAND(baseline,'nLooseLep==0 && nLooseElectron==0 && nLoosePhoton==0 && jetNMBtags==1 && jet1Pt>50 && abs(jet1Eta)<2.5 && jet2Pt>30 && abs(jet2Eta)<2.5 && dphipfmet>0.5 && pfUmag>250 && calocutSR<0.5 && cat1')
cuts['SR2'] = tAND(baseline,'nLooseLep==0 && nLooseElectron==0 && nLoosePhoton==0 && jetNMBtags==2 && abs(jet1Eta)<2.5 && abs(jet2Eta)<2.5 && abs(jet3Eta)<2.5 && dphipfmet>0.5 && pfUmag>250 && calocutSR<0.5 && cat2')

  #inclusive
cuts['ZmmINC'] = tAND(baseline,'nLooseMuon==2 && looseLep1IsTight==1 && diLepMass>70 && diLepMass<110')
cuts['ZeeINC'] = tAND(baseline,'nLooseElectron==2 && looseLep1IsTight==1 && diLepMass>70 && diLepMass<110')
cuts['WmnINC'] = tAND(baseline,'nLooseMuon==1 && looseLep1IsTight==1')
cuts['WenINC'] = tAND(baseline,'nLooseElectron==1 && looseLep1IsTight==1')
cuts['TemINC'] = tAND(baseline,'nLooseMuon==1 && nLooseElectron==1 && looseLep1IsTight==1 && looseLep2IsTight==1 && isojetNBtags==1 && (looseLep1PdgId==13 && looseLep2PdgId==-11 || looseLep1PdgId==-13 && looseLep2PdgId==11)    ')

cuts['ZmmCR'] = tAND(cuts['ZmmINC'],'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutZ<0.5 && pfUZmag>250 && dphipfUZ>0.5 && jet1Pt>50 && abs(jet1Eta)<2.5')
cuts['ZeeCR'] = tAND(cuts['ZeeINC'],'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutZ<0.5 && pfUZmag>250 && dphipfUZ>0.5 && jet1Pt>50 && abs(jet1Eta)<2.5')
cuts['WmnCR'] = tAND(cuts['WmnINC'],'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutW<0.5 && pfUWmag>250 && dphipfUW>0.5 && jet1Pt>50 && abs(jet1Eta)<2.5')
cuts['WenCR'] = tAND(cuts['WenINC'],'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutW<0.5 && pfUWmag>250 && dphipfUW>0.5 && jet1Pt>50 && abs(jet1Eta)<2.5')
cuts['TemCR'] = tAND(cuts['TemINC'],'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutW<0.5 && pfUWmag>250 && dphipfUW>0.5 && jet1Pt>50 && abs(jet1Eta)<2.5')

cuts['ZmmbCR'] = tAND(cuts['ZmmINC'],'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutZ<0.5 && pfUZmag>250 && dphipfUZ>0.5 && jet1Pt>50 && abs(jet1Eta)<2.5 && cat1')
cuts['ZeebCR'] = tAND(cuts['ZeeINC'],'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutZ<0.5 && pfUZmag>250 && dphipfUZ>0.5 && jet1Pt>50 && abs(jet1Eta)<2.5 && cat1')
cuts['WmnbCR'] = tAND(cuts['WmnINC'],'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutW<0.5 && pfUWmag>250 && dphipfUW>0.5 && jet1Pt>50 && abs(jet1Eta)<2.5 && cat1')
cuts['WenbCR'] = tAND(cuts['WenINC'],'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutW<0.5 && pfUWmag>250 && dphipfUW>0.5 && jet1Pt>50 && abs(jet1Eta)<2.5 && cat1')
cuts['TembCR'] = tAND(cuts['TemINC'],'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutW<0.5 && pfUWmag>250 && dphipfUW>0.5 && jet1Pt>50 && abs(jet1Eta)<2.5 && cat1')

cuts['ZmmbbCR'] = tAND(cuts['ZmmINC'],'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutZ<0.5 && pfUZmag>250 && dphipfUZ>0.5 && jet1Pt>50 && abs(jet1Eta)<2.5 && cat1')
cuts['ZeebbCR'] = tAND(cuts['ZeeINC'],'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutZ<0.5 && pfUZmag>250 && dphipfUZ>0.5 && jet1Pt>50 && abs(jet1Eta)<2.5 && cat1')
cuts['WmnbbCR'] = tAND(cuts['WmnINC'],'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutW<0.5 && pfUWmag>250 && dphipfUW>0.5 && jet1Pt>50 && abs(jet1Eta)<2.5 && cat1')
cuts['WenbbCR'] = tAND(cuts['WenINC'],'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutW<0.5 && pfUWmag>250 && dphipfUW>0.5 && jet1Pt>50 && abs(jet1Eta)<2.5 && cat1')
cuts['TembbCR'] = tAND(cuts['TemINC'],'(nLooseElectron+nLoosePhoton+nTau)==0 && calocutW<0.5 && pfUWmag>250 && dphipfUW>0.5 && jet1Pt>50 && abs(jet1Eta)<2.5 && cat1')

##weights for specific regions                                                                   
#base
weights['base'] = 'normalizedWeight*sf_pu*sf_ewkV*sf_qcdV'

#signal region in all b-tag categories
weights['SR0'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_tt*sf_Medbtag0')
weights['SR1'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_tt*sf_Medbtag1')
weights['SR2'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_tt*sf_Medbtag2')

#inclusive 
weights['ZmmINC'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt')
weights['ZeeINC'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt')
weights['WmnINC'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt')
weights['WenINC'] = tTIMES(weights['base'],'%f*sf_eleTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt')
weights['TemINC'] = tTIMES(weights['base'],'%f*sf_metTrig*sf_lepID*sf_lepIso*sf_lepTrack*sf_tt*sf_Medbtag1')

#Control Region Specific
#0 b-tag
weights['ZmmCR'] = tTIMES(weights['ZmmINC'],'%f*sf_Medbtag0')
weights['ZeeCR'] = tTIMES(weights['ZeeINC'],'%f*sf_Medbtag0')
weights['WmnCR'] = tTIMES(weights['WmnINC'],'%f*sf_Medbtag0')
weights['WenCR'] = tTIMES(weights['WenINC'],'%f*sf_Medbtag0')
weights['TemCR'] = tTIMES(weights['TemINC'],'%f*sf_Medbtag0')

#1 b-tag
weights['ZmmbCR'] = tTIMES(weights['ZmmINC'],'%f*sf_Medbtag1')
weights['ZeebCR'] = tTIMES(weights['ZeeINC'],'%f*sf_Medbtag1')
weights['WmnbCR'] = tTIMES(weights['WmnINC'],'%f*sf_Medbtag1')
weights['WenbCR'] = tTIMES(weights['WenINC'],'%f*sf_Medbtag1')
weights['TembCR'] = tTIMES(weights['TemINC'],'%f*sf_Medbtag1')

#2 b-tag                                                                                                                                                                         
weights['ZmmbbCR'] = tTIMES(weights['ZmmINC'],'%f*sf_Medbtag2')
weights['ZeebbCR'] = tTIMES(weights['ZeeINC'],'%f*sf_Medbtag2')
weights['WmnbbCR'] = tTIMES(weights['WmnINC'],'%f*sf_Medbtag2')
weights['WenbbCR'] = tTIMES(weights['WenINC'],'%f*sf_Medbtag2')
weights['TembbCR'] = tTIMES(weights['TemINC'],'%f*sf_Medbtag2')

#b-tag shift
for r in ['SR0','SR1','SR2','ZmmCR','ZeeCR','WmnCR','WenCR','TemCR','ZmmbCR','ZeebCR','WmnbCR','WenbCR','TembCR','ZmmbbCR','ZeebbCR','WmnbbCR','WenbbCR','TembbCR']:
    for shift in ['BUp','BDown','MUp','MDown']:
        for cent in ['sf_Medbtag']:
            if 'Medbtag0' in weights[r]:
                weights[r+'_'+cent+shift] = sub(cent+'0',cent+'0'+shift,weights[r])
            if 'Medbtag1' in weights[r]:
                weights[r+'_'+cent+shift] = sub(cent+'1',cent+'1'+shift,weights[r])
            if 'Medbtag2' in weights[r]:
                weights[r+'_'+cent+shift] = sub(cent+'2',cent+'2'+shift,weights[r])

#scale and pdf shift
for r in ['SR0','SR1','SR2','ZmmCR','ZeeCR','WmnCR','WenCR','TemCR','ZmmbCR','ZeebCR','WmnbCR','WenbCR','TembCR','ZmmbbCR','ZeebbCR','WmnbbCR','WenbbCR','TembbCR']:
    for shift in ['ScaleUp','ScaleDown','PDFUp','PDFDown']:
        if shift == 'ScaleUp':
            weights[r+'_'+shift] = weights[r] + "*(scale[3]+1)"
        if shift == 'ScaleDown':
            weights[r+'_'+shift] = weights[r] + "*(scale[5]+1)"
        if shift == 'PDFUp':
            weights[r+'_'+shift] = weights[r] + "*pdfUp"
        if shift == 'PDFDown':
            weights[r+'_'+shift] = weights[r] + "*pdfDown"
