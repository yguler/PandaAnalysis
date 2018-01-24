#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define EGMSCALE 1

using namespace panda;
using namespace std;

void PandaAnalyzer::IncrementGenAuxFile(bool close)
{
  if (fAux) {
    fAux->WriteTObject(tAux, "inputs", "Overwrite");
    fAux->Close();
  }
  if (close)
    return;

  TString path = TString::Format(auxFilePath.Data(),auxCounter++);
  fAux = TFile::Open(path.Data(), "RECREATE");
  if (DEBUG) PDebug("PandaAnalyzer::IncrementAuxFile", "Opening "+path);
  tAux = new TTree("inputs","inputs");
  
  genJetInfo.particles.resize(NMAXPF);
  for (unsigned i = 0; i != NMAXPF; ++i) {
    genJetInfo.particles[i].resize(NGENPROPS);
  }
  tAux->Branch("kinematics",&(genJetInfo.particles));

  tAux->Branch("eventNumber",&(gt->eventNumber),"eventNumber/l");
  tAux->Branch("nprongs",&(genJetInfo.nprongs),"nprongs/I");
  tAux->Branch("partonpt",&(genJetInfo.partonpt),"partonpt/F");
  tAux->Branch("partonm",&(genJetInfo.partonm),"partonm/F");
  tAux->Branch("pt",&(genJetInfo.pt),"pt/F");
  tAux->Branch("msd",&(genJetInfo.msd),"msd/F");
  tAux->Branch("eta",&(genJetInfo.eta),"eta/F");
  tAux->Branch("phi",&(genJetInfo.phi),"phi/F");
  tAux->Branch("m",&(genJetInfo.m),"m/F");
  tAux->Branch("tau3",&(genJetInfo.tau3),"tau3/F");
  tAux->Branch("tau2",&(genJetInfo.tau2),"tau2/F");
  tAux->Branch("tau1",&(genJetInfo.tau1),"tau1/F");
  tAux->Branch("tau3sd",&(genJetInfo.tau3sd),"tau3sd/F");
  tAux->Branch("tau2sd",&(genJetInfo.tau2sd),"tau2sd/F");
  tAux->Branch("tau1sd",&(genJetInfo.tau1sd),"tau1sd/F");

  fOut->cd();

  if (tr)
    tr->TriggerEvent("increment aux file");
}

void PandaAnalyzer::IncrementAuxFile(bool close)
{
  if (fAux) {
    fAux->WriteTObject(tAux, "inputs", "Overwrite");
    fAux->Close();
  }
  if (close)
    return;

  TString path = TString::Format(auxFilePath.Data(),auxCounter++);
  fAux = TFile::Open(path.Data(), "RECREATE");
  if (DEBUG) PDebug("PandaAnalyzer::IncrementAuxFile", "Opening "+path);
  tAux = new TTree("inputs","inputs");
  
  pfInfo.resize(NMAXPF);
  for (unsigned i = 0; i != NMAXPF; ++i) {
    pfInfo[i].resize(NPFPROPS);
  }
  tAux->Branch("kinematics",&pfInfo);
  
  svInfo.resize(NMAXSV);
  for (unsigned i = 0; i != NMAXSV; ++i) {
    svInfo[i].resize(NSVPROPS);
  }
  tAux->Branch("svs",&svInfo);

  tAux->Branch("msd",&fjmsd,"msd/F");
  tAux->Branch("pt",&fjpt,"pt/F");
  tAux->Branch("rawpt",&fjrawpt,"rawpt/F");
  tAux->Branch("eta",&fjeta,"eta/F");
  tAux->Branch("phi",&fjphi,"phi/F");
  tAux->Branch("rho",&(gt->fj1Rho),"rho/f");
  tAux->Branch("rawrho",&(gt->fj1RawRho),"rawrho/f");
  tAux->Branch("rho2",&(gt->fj1Rho2),"rho2/f");
  tAux->Branch("rawrho2",&(gt->fj1RawRho2),"rawrho2/f");
  tAux->Branch("nPartons",&(gt->fj1NPartons),"nPartons/I");
  tAux->Branch("nBPartons",&(gt->fj1NBPartons),"nBPartons/I");
  tAux->Branch("nCPartons",&(gt->fj1NCPartons),"nCPartons/I");
  tAux->Branch("partonM",&(gt->fj1PartonM),"partonM/f");
  tAux->Branch("partonPt",&(gt->fj1PartonPt),"partonPt/f");
  tAux->Branch("partonEta",&(gt->fj1PartonEta),"partonEta/f");
  tAux->Branch("tau32",&(gt->fj1Tau32),"tau32/f");
  tAux->Branch("tau32SD",&(gt->fj1Tau32SD),"tau32SD/f");
  tAux->Branch("tau21",&(gt->fj1Tau21),"tau21/f");
  tAux->Branch("tau21SD",&(gt->fj1Tau21SD),"tau21SD/f");
  tAux->Branch("eventNumber",&(gt->eventNumber),"eventNumber/l");
  tAux->Branch("maxcsv",&(gt->fj1MaxCSV),"maxcsv/f");
  tAux->Branch("mincsv",&(gt->fj1MinCSV),"mincsv/f");
  tAux->Branch("doubleb",&(gt->fj1DoubleCSV),"doubleb/f");

  gt->SetAuxTree(tAux);

  fOut->cd();

  if (tr)
    tr->TriggerEvent("increment aux file");
}

void PandaAnalyzer::RegisterTriggers() 
{
  for (auto &th : triggerHandlers) {
    unsigned N = th.paths.size();
    for (unsigned i = 0; i != N; i++) {
      unsigned panda_idx = event.registerTrigger(th.paths.at(i));
      th.indices[i] = panda_idx;
      if (DEBUG) PDebug("PandaAnalyzer::RegisterTriggers",
        Form("Got index %d for trigger path %s", panda_idx, th.paths.at(i).Data())
      );
    }
  }
}

bool PandaAnalyzer::RecoilPresel() 
{
    if ( (preselBits&kMonotop) || (preselBits&kMonohiggs) || 
         (preselBits&kMonojet) || (preselBits&kRecoil) ) 
    {
       if (event.recoil.max<175)
         return false;
    } 
    return true;
}

void PandaAnalyzer::TriggerEffs()
{

    // trigger efficiencies
    gt->sf_metTrig = GetCorr(cTrigMET,gt->pfmetnomu);
    gt->sf_metTrigZmm = GetCorr(cTrigMETZmm,gt->pfmetnomu);

    if (gt->nLooseElectron>0) {
      panda::Electron *ele1=0, *ele2=0;
      if (gt->nLooseLep>0) ele1 = dynamic_cast<panda::Electron*>(looseLeps[0]);
      if (gt->nLooseLep>1) ele2 = dynamic_cast<panda::Electron*>(looseLeps[1]);
      float eff1=0, eff2=0;
      if (ele1 && ele1->tight) {
        eff1 = GetCorr(cTrigEle, ele1->eta(), ele1->pt());
        if (ele2 && ele2->tight)
          eff2 = GetCorr(cTrigEle, ele2->eta(), ele2->pt());
        gt->sf_eleTrig = 1 - (1-eff1)*(1-eff2);
      }
    } // done with ele trig SF
    if (gt->nLooseMuon>0) {
      panda::Muon *mu1=0, *mu2=0;
      if (gt->nLooseLep>0) mu1 = dynamic_cast<panda::Muon*>(looseLeps[0]);
      if (gt->nLooseLep>1) mu2 = dynamic_cast<panda::Muon*>(looseLeps[1]);
      float eff1=0, eff2=0;
      if (mu1 && mu1->tight) {
        eff1 = GetCorr(
          cTrigMu,
          fabs(mu1->eta()),
          TMath::Max((float)26.,TMath::Min((float)499.99,(float)mu1->pt()))
        );
        if (mu2 && mu2->tight)
          eff2 = GetCorr(
            cTrigMu,
            fabs(mu2->eta()),
            TMath::Max((float)26.,TMath::Min((float)499.99,(float)mu2->pt()))
          );
        gt->sf_muTrig = 1 - (1-eff1)*(1-eff2);
      }
    } // done with mu trig SF

    if (gt->nLoosePhoton>0 && gt->loosePho1IsTight)
      gt->sf_phoTrig = GetCorr(cTrigPho,gt->loosePho1Pt);

    if (analysis->vbf) {
      gt->sf_metTrigVBF = GetCorr(cVBF_TrigMET,gt->barrelHTMiss);
      gt->sf_metTrigZmmVBF = GetCorr(cVBF_TrigMETZmm,gt->barrelHTMiss);
    }
    tr->TriggerEvent("triggers");
}

void PandaAnalyzer::Recoil()
{
    TLorentzVector vpfUp; vpfUp.SetPtEtaPhiM(gt->pfmetUp,0,gt->pfmetphi,0);
    TLorentzVector vpfDown; vpfDown.SetPtEtaPhiM(gt->pfmetDown,0,gt->pfmetphi,0);
    TLorentzVector vpfUWUp, vpfUWDown;
    TLorentzVector vObj1, vObj2;
    gt->whichRecoil = 0; // -1=photon, 0=MET, 1,2=nLep
    if (gt->nLooseLep>0) {
      panda::Lepton *lep1 = looseLeps.at(0);
      vObj1.SetPtEtaPhiM(lep1->pt(),lep1->eta(),lep1->phi(),lep1->m());

      // one lep => W
      vpuppiUW = vPuppiMET+vObj1; gt->puppiUWmag=vpuppiUW.Pt(); gt->puppiUWphi=vpuppiUW.Phi();
      vpfUW = vPFMET+vObj1; gt->pfUWmag=vpfUW.Pt(); gt->pfUWphi=vpfUW.Phi();
      
      if (analysis->varyJES) {
        vpfUWUp = vpfUp+vObj1; gt->pfUWmagUp = vpfUWUp.Pt();
        vpfUWDown = vpfDown+vObj1; gt->pfUWmagDown = vpfUWDown.Pt();
      }

      if (gt->nLooseLep>1 && looseLep1PdgId+looseLep2PdgId==0) {
        // two OS lep => Z
        panda::Lepton *lep2 = looseLeps.at(1);
        vObj2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());

        vpuppiUZ=vpuppiUW+vObj2; gt->puppiUZmag=vpuppiUZ.Pt(); gt->puppiUZphi=vpuppiUZ.Phi();
        vpfUZ=vpfUW+vObj2; gt->pfUZmag=vpfUZ.Pt(); gt->pfUZphi=vpfUZ.Phi();

        if (analysis->varyJES) {
          TLorentzVector vpfUZUp = vpfUWUp+vObj2; gt->pfUZmagUp = vpfUZUp.Pt();
          TLorentzVector vpfUZDown = vpfUWDown+vObj2; gt->pfUZmagDown = vpfUZDown.Pt();
        }

        vpuppiU = vpuppiUZ; vpfU = vpfUZ;
        gt->whichRecoil = 2;
      } else {
        vpuppiU = vpuppiUW; vpfU = vpfUW;
        gt->whichRecoil = 1;
      }
    }
    if (gt->nLoosePhoton>0) {
      panda::Photon *pho = loosePhos.at(0);
      vObj1.SetPtEtaPhiM(pho->pt(),pho->eta(),pho->phi(),0.);

      vpuppiUA=vPuppiMET+vObj1; gt->puppiUAmag=vpuppiUA.Pt(); gt->puppiUAphi=vpuppiUA.Phi();
      vpfUA=vPFMET+vObj1; gt->pfUAmag=vpfUA.Pt(); gt->pfUAphi=vpfUA.Phi();

      if (analysis->varyJES) {
        TLorentzVector vpfUAUp = vpfUp+vObj1; gt->pfUAmagUp = vpfUAUp.Pt();
        TLorentzVector vpfUADown = vpfDown+vObj1; gt->pfUAmagDown = vpfUADown.Pt();
      }

      if (gt->nLooseLep==0) {
        vpuppiU = vpuppiUA; vpfU = vpfUA;
        gt->whichRecoil = -1;
      }
    }
    if (gt->nLooseLep==0 && gt->nLoosePhoton==0) {
      vpuppiU = vPuppiMET;
      vpfU = vPFMET;
      gt->whichRecoil = 0;
    }
    gt->puppiUmag = vpuppiU.Pt();
    gt->puppiUphi = vpuppiU.Phi();
    gt->pfUmag = vpfU.Pt();
    gt->pfUphi = vpfU.Phi();

    tr->TriggerEvent("recoils");
}

void PandaAnalyzer::HeavyFlavorCounting() 
{
  // For now, simple B and C counting
  for (auto& gen : event.genParticles) {
    float pt = gen.pt();
    int pdgid = gen.pdgid;
    if (gen.parent.isValid() && gen.parent->pdgid==gen.pdgid)
      continue;
    //count bs and cs
    int apdgid = abs(pdgid);
    if (apdgid!=5 && apdgid!=4) 
      continue;
    if (gen.pt()>5) {
      gt->nHF++;
      if (apdgid==5)
        gt->nB++;
    }
  }
}

void PandaAnalyzer::GetMETSignificance()
{
  /*
  std::vector<panda::Lepton*> msLeps;
  // Get the lepton collections for the MET significance calculation
  // Ported approximately from RecoMET/METProducers/python/METSignificanceObjects_cfi.py
  for (auto& ele : event.electrons) {
    if(!ele.superCluster.isValid()) continue;
    if(!ele.matchedPF.isValid()) continue;
    float abseta = fabs(ele.eta());
    // WARNING: Missing Track-Cluster Matching
    // barrel abs(deltaEtaSuperClusterTrackAtVtx) < 0.007, abs(deltaPhiSuperClusterTrackAtVtx) < 0.8
    // endcap abs(deltaEtaSuperClusterTrackAtVtx) < 0.009, abs(deltaPhiSuperClusterTrackAtVtx) < 0.10
    panda::SuperCluster *theSC = (panda::SuperCluster*)ele.superCluster.get();
    // panda::PFCand *theTrack = ele.matchedPF.get();
    // float dEtaSCTrack = ?
    // float dPhiSCTrack = ?
    float ooEmooP = fabs(1./theSC->rawPt - 1./ele.trackP);
    if (
      // Common cuts
      ele.pt() > 19.5 && abseta < 2.5 && 
      ele.nMissingHits <= 1 && 
      ele.combIso()/ele.pt() < 0.3 &&  // Warning: Not the Delta-R 0.3 Isolation
      (
        ( 
          // Barrel cuts
          abseta < 1.4442 &&
          ele.sieie < 0.01 &&
          ele.hOverE < 0.15 &&
          ooEmooP < 0.05
        ) || ( 
          // Endcap cuts
          abseta > 1.566 &&
          ele.sieie < 0.03 &&
          ele.hOverE < 0.10 &&
          ooEmooP < 0.05
        )
      )
    ) msLeps.push_back(&ele);
  }
  for (auto& mu : event.muons) {
    float abseta = fabs(mu.eta());
    if (
      mu.tracker==true &&
      mu.pt() > 5 &&
      mu.pf==true &&
      mu.global==true && // Same as globalTrack.isNonnull ?
      mu.pixLayersWithMmt > 0 &&
      mu.normChi2 < 10 &&
      mu.nMatched > 0 &&
      mu.trkLayersWithMmt > 5 &&
      mu.nValidMuon > 0 &&
      mu.r03Iso/mu.pt() < 0.3 &&
      fabs(mu.dxy) < 2.0
    ) msLeps.push_back(&mu);
  }
  // metsig covariance
  double cov_xx=0, cov_yy=0, cov_xy=0;
  
  std::vector<panda::PFCand*> footprint; // jet+leptons footprint
  footprint.reserve(event.pfCandidates.size());
  for (unsigned iL=0; iL<(unsigned)msLeps.size(); iL++) {
    if (msLeps[iL]->pt()<10) continue;
    if (msLeps[iL]->matchedPF.isValid())
      footprint.push_back((panda::PFCand*)msLeps[iL]->matchedPF.get());
  }
  std::vector<bool> jetIsLep; jetIsLep.reserve(event.chsAK4Jets.size());
  for (unsigned iJ=0; iJ<(unsigned)event.chsAK4Jets.size(); iJ++) { 
    panda::Jet jet = event.chsAK4Jets[iJ];
    jetIsLep.push_back(false);
    // jet-lepton cleaning
    for (unsigned iL=0; iL<(unsigned)msLeps.size(); iL++) {
      if (msLeps[iL]->p4().DeltaR( jet.p4()) < 0.4) { 
        jetIsLep[iJ]=true;
        break; 
      }
    }
    if (jetIsLep[iJ]) continue;
    RefVector<PFCand> jetCands = jet.constituents;
    for (UShort_t iJC=0; iJC<jetCands.size(); iJC++) { 
      if (!jetCands.at(iJC).isValid()) continue;
      footprint.push_back( (panda::PFCand*) jetCands.at(iJC).get());
    }
  }
  */
  // Calculate sumPt not including the footprint
  double sumPt=0;
  float puppiEt = 0;
  float pfEt = 0;
  TLorentzVector pfcand(0,0,0,0);
  for (auto& pfCand : event.pfCandidates) {
    pfcand.SetPtEtaPhiM(pfCand.pt(),pfCand.eta(),pfCand.phi(),pfCand.m());
    puppiEt += pfcand.Et()*pfCand.puppiW();
    pfEt += pfcand.Et(); continue;
    bool candIsInFootprint=false;
    for (unsigned iFP=0; iFP<(unsigned)footprint.size(); iFP++) {
      if (footprint[iFP] == &pfCand ) candIsInFootprint=true;
    }
    if(candIsInFootprint) continue;
    sumPt += pfCand.pt();
  }
  gt->puppimetsig = event.puppiMet.pt/sqrt(puppiEt);
  gt->pfmetsig = event.pfMet.pt/sqrt(pfEt);
  tr->TriggerEvent("MET significance");
  return;
  
  /*
  // Add jets to covariance matrix
  JME::JetParameters parameters;
  
  for (unsigned iJ=0; iJ<(unsigned)event.chsAK4Jets.size(); iJ++) { 
    panda::Jet jet = event.chsAK4Jets[iJ];
    if(jetIsLep[iJ]) continue;
    parameters.setJetPt(jet.pt());
    parameters.setJetEta(jet.eta());
    parameters.setRho(event.rho);
    // jet energy resolutions
    double sigmapt = 0.01; //resPtObj.getResolution(parameters);
    double sigmaphi = 0.001; //resPhiObj.getResolution(parameters);
    // split into high-pt and low-pt sector
    if( jet.pt() > 15. ){
      // high-pt jets enter into the covariance matrix via JER
      double scale = 0;
      float abseta=fabs(jet.eta());
      if (isData) {
        if      (abseta<msJetEtas[0]) scale = msJetParametersData[0];
        else if (abseta<msJetEtas[1]) scale = msJetParametersData[1];
        else if (abseta<msJetEtas[2]) scale = msJetParametersData[2];
        else if (abseta<msJetEtas[3]) scale = msJetParametersData[3];
        else                          scale = msJetParametersData[4];
      } else {
        if      (abseta<msJetEtas[0]) scale = msJetParametersMC[0];
        else if (abseta<msJetEtas[1]) scale = msJetParametersMC[1];
        else if (abseta<msJetEtas[2]) scale = msJetParametersMC[2];
        else if (abseta<msJetEtas[3]) scale = msJetParametersMC[3];
        else                          scale = msJetParametersMC[4];
      }
      double dpt = scale*jet.pt()*sigmapt;
      double dph = jet.pt()*sigmaphi;
      double c   = jet.p4().Px()/jet.pt();
      double s   = jet.p4().Px()/jet.pt();
      cov_xx += dpt*dpt*c*c + dph*dph*s*s;
      cov_xy += (dpt*dpt-dph*dph)*c*s;
      cov_yy += dph*dph*c*c + dpt*dpt*s*s;
    } else {
       // add the (corrected) jet to the sumPt
       sumPt += jet.pt();
    }
  }
  //protection against unphysical events
  if(sumPt<0) sumPt=0;

  // add pseudo-jet to metsig covariance matrix
  if(isData) {
    cov_xx += pow(msPseudoJetParametersData[0],2) + pow(msPseudoJetParametersData[1],2)*sumPt;
    cov_yy += pow(msPseudoJetParametersData[0],2) + pow(msPseudoJetParametersData[1],2)*sumPt;
  } else {
    cov_xx += pow(msPseudoJetParametersMC[0],2) + pow(msPseudoJetParametersMC[1],2)*sumPt;
    cov_yy += pow(msPseudoJetParametersMC[0],2) + pow(msPseudoJetParametersMC[1],2)*sumPt;
  }

  // covariance matrix determinant
  double det = cov_xx*cov_yy - cov_xy*cov_xy;
  // invert matrix
  double ncov_xx = cov_yy / det; 
  double ncov_xy = -cov_xy / det;
  double ncov_yy = cov_xx / det; 
  // product of met and inverse of covariance
  double metpx = vPFMET.Px(), metpy = vPFMET.Py();
  double sig = metpx*metpx*ncov_xx + 2*metpx*metpy*ncov_xy + metpy*metpy*ncov_yy;                  
  gt->pfmetsig = sig; 
  tr->TriggerEvent("MET significance");
  */
}

