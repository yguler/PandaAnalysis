#ifndef PandaAnalyzer_h
#define PandaAnalyzer_h

// STL
#include "vector"
#include <unordered_set>
#include "map"
#include <string>
#include <cmath>

// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TRandom3.h>
#include <TLorentzVector.h>

#include "AnalyzerUtilities.h"
#include "GeneralTree.h"

// btag
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
//#include "BTagCalibrationStandalone.h"

// JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "PandaAnalysis/Utilities/interface/RoccoR.h"
#include "PandaAnalysis/Utilities/interface/CSVHelper.h"

// TMVA
#include "TMVA/Reader.h"

/////////////////////////////////////////////////////////////////////////////
// some misc definitions

#define NMAXPF 100
#define NMAXSV 10


/////////////////////////////////////////////////////////////////////////////
// PandaAnalyzer definition
class PandaAnalyzer {
public :
    // configuration enums
    enum PreselectionBit {
     kMonotop    =(1<<0),
     kMonohiggs  =(1<<2),
     kMonojet    =(1<<3),
     kPassTrig   =(1<<4),
     kVBF        =(1<<5),
     kRecoil     =(1<<6),
     kFatjet     =(1<<7),
     kFatjet450  =(1<<8),
     kRecoil50   =(1<<9),
     kGenBosonPt =(1<<10),
     kGenFatJet  =(1<<11),
     kVHBB       =(1<<12),
     kLepton     =(1<<13),
     kLeptonFake =(1<<14)
    };
    
    enum LepSelectionBit {
     kLoose   =(1<<0),
     kFake    =(1<<1),
     kMedium  =(1<<2),
     kTight   =(1<<3),
     kDxyz    =(1<<4)
    };

    enum TriggerBits {
        kMETTrig       = 0,
        kSingleEleTrig,
        kSinglePhoTrig,
        kSingleMuTrig,
        kDoubleMuTrig,
        kDoubleEleTrig,
        kEMuTrig,
        kJetHTTrig,
        kMuFakeTrig,
        kEleFakeTrig,
        kNTrig
    };

    //////////////////////////////////////////////////////////////////////////////////////

    PandaAnalyzer(int debug_=0);
    ~PandaAnalyzer();
    int Init(TTree *tree, TH1D *hweights, TTree *weightNames=0);
    void SetOutputFile(TString fOutName);
    void ResetBranches();
    void Run();
    void Terminate();
    void SetDataDir(const char *s);
    void SetPreselectionBit(PreselectionBit b,bool on=true) {
        if (on) 
            preselBits |= b;
        else 
            preselBits &= ~b;
    }
    void AddGoodLumiRange(int run, int l0, int l1);

    // public configuration
    void SetAnalysis(Analysis *a) { analysis = a; }
    bool isData=false;              // to do gen matching, etc
    int firstEvent=-1;
    int lastEvent=-1;               // max events to process; -1=>all

private:
    enum CorrectionType { //!< enum listing relevant corrections applied to MC
        cNPV=0,       //!< npv weight
        cPU,          //!< true pu weight
        cPUUp,        //!< true pu weight
        cPUDown,      //!< true pu weight
        cEleVeto,     //!< monojet SF, Veto ID for e
        cEleLoose,    //!< monojet SF, Tight ID for e
        cEleMedium,   //!< monojet SF, Tight ID for e
        cEleTight,    //!< monojet SF, Tight ID for e
        cEleReco,     //!< monojet SF, tracking for e
        cZHEwkCorr,     //!< ZH Ewk Corr weight  
        cZHEwkCorrUp,   //!< ZH Ewk Corr weight Up  
        cZHEwkCorrDown, //!< ZH Ewk Corr weight Down  
        cWZEwkCorr,
        cqqZZQcdCorr,
        cMuLooseID,   //!< MUO POG SF, Loose ID for mu 
        cMuMediumID,  //!< MUO POG SF, Tight ID for mu 
        cMuTightID,   //!< MUO POG SF, Tight ID for mu 
        cMuLooseIso,  //!< MUO POG SF, Loose Iso for mu 
        cMuMediumIso, //!< MUO POG SF, Loose Iso for mu 
        cMuTightIso,  //!< MUO POG SF, Tight Iso for mu 
        cMuReco,      //!< MUO POG SF, tracking for mu
        cPho,         //!< EGM POG SF, contains ID for gamma
        cTrigMET,     //!< MET trigger eff        
        cTrigMETZmm,  //!< Zmumu MET trigger eff
        cTrigEle,     //!< Ele trigger eff        
        cTrigMu,      //!< Mu trigger eff        
        cTrigPho,     //!< Pho trigger eff        
        cZNLO,        //!< NLO weights for QCD Z,W,A,A+2j
        cWNLO,
        cANLO,
        cANLO2j,
        cZEWK,        //!< EWK weights for QCD Z,W,A,A+2j
        cWEWK,
        cAEWK,
        cVBF_ZNLO,    //!< NLO weights for QCD Z,W in VBF phase space
        cVBF_ZllNLO,  
        cVBF_WNLO,
        cVBFTight_ZNLO,    //!< NLO weights for QCD Z,W in tight VBF phase space
        cVBFTight_ZllNLO,  
        cVBFTight_WNLO,
        cVBF_EWKZ,    //!< k-factors for EWK Z,W in VBF phase space
        cVBF_EWKW,
        cVBF_TrigMET, //!< MET trigger eff as a f'n of mjj/met 
        cVBF_TrigMETZmm,
        cBadECALJets,  //!< bad ECAL clusters to filter jets
        cN
    };

    enum BTagType {
        bJetL=0,
        bSubJetL,
        bJetM,
        bN
    };

    class btagcand {
        public:
            btagcand(unsigned int i, int f,double e,double cent,double up,double down) {
                idx = i;
                flav = f;
                eff = e;
                sf = cent;
                sfup = up;
                sfdown = down;
            }
            ~btagcand() { }
            int flav, idx;
            double eff, sf, sfup, sfdown;
    };

    struct GenJetInfo {
      float pt=0, eta=0, phi=0, m=0;
      float msd=0;
      float tau3=0, tau2=0, tau1=0;
      float tau3sd=0, tau2sd=0, tau1sd=0;
      int nprongs=0;
      float partonpt=0, partonm=0;
      std::vector<std::vector<float>> particles;
    };

    //////////////////////////////////////////////////////////////////////////////////////

    bool PassGoodLumis(int run, int lumi);
    bool PassPreselection();
    void OpenCorrection(CorrectionType,TString,TString,int);
    double GetCorr(CorrectionType ct,double x, double y=0);
    double GetError(CorrectionType ct,double x, double y=0);
    void RegisterTriggers(); 
    void GetMETSignificance(); 

    // these are functions used for analysis-specific tasks inside Run.
    // ideally the return type is void (e.g. they are stateful functions),
    // but that is not always possible (e.g. RecoilPresel)
    void CalcBJetSFs(BTagType bt, int flavor, double eta, double pt, 
                     double eff, double uncFactor, double &sf, double &sfUp, double &sfDown);
    void ComplicatedLeptons();
    void EvalBTagSF(std::vector<btagcand> &cands, std::vector<double> &sfs,
                    GeneralTree::BTagShift shift,GeneralTree::BTagJet jettype, bool do2=false);
    void IncrementAuxFile(bool close=false);
    void IncrementGenAuxFile(bool close=false);
    void FatjetBasics();
    void FatjetMatching();
    void FatjetPartons();
    void FatjetRecluster();
    template <typename T> void FillGenTree(panda::Collection<T>& genParticles);
    void FillPFTree();
    void GenFatJet();
    void GenJetsNu();
    void GenStudyEWK();
    float GetMSDCorr(float, float); 
    void HeavyFlavorCounting();
    void IsoJet(panda::Jet&);
    void JetBRegressionInfo(panda::Jet&);
    void JetBasics();
    void JetBtagSFs();
    void JetCMVAWeights();
    void JetHbbBasics(panda::Jet&);
    void JetHbbReco();
    void JetHbbSoftActivity();
    void JetVBFBasics(panda::Jet&);
    void JetVBFSystem();
    void JetVaryJES(panda::Jet&);
    void LeptonSFs();
    template <typename T> void MatchGenJets(T& genJets);
    void PhotonSFs();
    void Photons();
    void QCDUncs();
    void Recoil();
    bool RecoilPresel();
    void SaveGenLeptons();
    void SetupJES();
    void SignalInfo();
    void SignalReweights();
    void SimpleLeptons();
    void Taus();
    void TopPTReweight();
    void TriggerEffs();
    void VJetsReweight();
    double WeightEWKCorr(float pt, int type);
    double WeightZHEWKCorr(float baseCorr);

    //////////////////////////////////////////////////////////////////////////////////////

    int DEBUG = 0; //!< debug verbosity level
    Analysis *analysis = 0; //!< configure what to run
    TimeReporter *tr = 0; //!< profile time usage
    float FATJETMATCHDR2 = 2.25;

    //////////////////////////////////////////////////////////////////////////////////////

    // stuff for matching objects
    std::map<panda::GenParticle const*,float> genObjects; 
        //!< particles we want to match the jets to, and the 'size' of the daughters
    panda::GenParticle const* MatchToGen(double eta, double phi, double r2, int pdgid=0);   
        //!< private function to match a jet; returns NULL if not found
    std::map<int,std::vector<LumiRange>> goodLumis;
    std::vector<panda::Particle*> matchPhos, matchEles, matchLeps;
    std::map<int, int> pdgToQ; 
    
    // fastjet reclustering
    fastjet::JetDefinition *jetDef=0;
    fastjet::JetDefinition *jetDefKt=0;
    fastjet::contrib::SoftDrop *softDrop=0;
    fastjet::contrib::Njettiness *tauN=0;
    fastjet::AreaDefinition *areaDef=0;
    fastjet::GhostedAreaSpec *activeArea=0;
    fastjet::JetDefinition *jetDefGen=0;
    fastjet::JetDefinition *softTrackJetDefinition=0;

    //////////////////////////////////////////////////////////////////////////////////////

    // CMSSW-provided utilities
    BTagCalibration *btagCalib=0;
    BTagCalibration *sj_btagCalib=0;
    std::vector<BTagCalibrationReader*> btagReaders = std::vector<BTagCalibrationReader*>(bN,0); //!< maps BTagType to a reader 

    Binner btagpt = Binner({});
    Binner btageta = Binner({});
    std::vector<std::vector<double>> lfeff, ceff, beff;
    TMVA::Reader *bjetreg_reader=0; 

    std::map<TString,JetCorrectionUncertainty*> ak8UncReader; //!< calculate JES unc on the fly
    JERReader *ak8JERReader{0}; //!< fatjet jet energy resolution reader
    std::map<TString,JetCorrectionUncertainty*> ak4UncReader; //!< calculate JES unc on the fly
    std::map<TString,FactorizedJetCorrector*> ak4ScaleReader; //!< calculate JES on the fly
    JERReader *ak4JERReader{0}; //!< fatjet jet energy resolution reader
    JetCorrectionUncertainty *uncReader=0;           
    JetCorrectionUncertainty *uncReaderAK4=0;        
    FactorizedJetCorrector *scaleReaderAK4=0;        

    EraHandler eras = EraHandler(2016); //!< determining data-taking era, to be used for era-dependent JEC
    ParticleGridder *grid = 0;

    //////////////////////////////////////////////////////////////////////////////////////

    // files and histograms containing weights
    std::vector<TFile*>   fCorrs  = std::vector<TFile*>  (cN,0); //!< files containing corrections
    std::vector<THCorr1*> h1Corrs = std::vector<THCorr1*>(cN,0); //!< histograms for binned corrections
    std::vector<THCorr2*> h2Corrs = std::vector<THCorr2*>(cN,0); //!< histograms for binned corrections
    std::vector<TF1Corr*> f1Corrs = std::vector<TF1Corr*>(cN,0); //!< TF1s for continuous corrections
    TFile *MSDcorr=0;
    TF1* puppisd_corrGEN=0;
    TF1* puppisd_corrRECO_cen=0;
    TF1* puppisd_corrRECO_for=0;
    RoccoR *rochesterCorrection=0;
    TRandom3 rng;
    CSVHelper *csvReweighter=0, *cmvaReweighter=0;

    //////////////////////////////////////////////////////////////////////////////////////

    // IO for the analyzer
    TString fOutPath;
    TFile *fOut=0;     // output file is owned by PandaAnalyzer
    TTree *tOut=0;
    GeneralTree *gt=0; // essentially a wrapper around tOut
    TString auxFilePath="";
    unsigned auxCounter=0;
    TFile *fAux=0; // auxillary file
    TTree *tAux=0;
    TH1F *hDTotalMCWeight=0;
    TTree *tIn=0;    // input tree to read
    unsigned int preselBits=0;
    panda::EventAnalysis event;

    //////////////////////////////////////////////////////////////////////////////////////

    // configuration read from output tree
    std::vector<int> ibetas;
    std::vector<int> Ns; 
    std::vector<int> orders;

    //////////////////////////////////////////////////////////////////////////////////////

    // any extra signal weights we want
    // stuff that gets passed between modules
    //
    // NB: ensure that any global vectors/maps that are per-event
    // are reset properly in ResetBranches(), or you can really
    // mess up behavior
    std::vector<TriggerHandler> triggerHandlers = std::vector<TriggerHandler>(kNTrig);
    std::vector<panda::Lepton*> looseLeps, tightLeps;
    std::vector<panda::Photon*> loosePhos;
    TLorentzVector vPFMET, vPuppiMET;
    TVector2 vMETNoMu;
    TLorentzVector vpfUW, vpfUZ, vpfUA, vpfU;
    TLorentzVector vpuppiUW, vpuppiUZ, vpuppiUA, vpuppiU;
    panda::FatJet *fj1 = 0;
    std::vector<panda::Jet*> cleanedJets, isoJets, centralJets, bCandJets;
    std::map<panda::Jet*,int> bCandJetGenFlavor;
    std::map<panda::Jet*,float> bCandJetGenPt;
    TLorentzVector vJet, vBarrelJets;
    panda::FatJetCollection *fatjets = 0;
    panda::JetCollection *jets = 0;
    panda::Jet *jot1 = 0, *jot2 = 0;
    panda::Jet *jotUp1 = 0, *jotUp2 = 0;
    panda::Jet *jotDown1 = 0, *jotDown2 = 0;
    panda::Jet *jetUp1 = 0, *jetUp2 = 0;
    panda::Jet *jetDown1 = 0, *jetDown2 = 0;
    std::vector<panda::GenJet> genJetsNu;
    float genBosonPtMin, genBosonPtMax;
    int looseLep1PdgId, looseLep2PdgId, looseLep3PdgId, looseLep4PdgId;
    std::vector<TString> wIDs;
    float *bjetreg_vars = 0;
    float jetPtThreshold=30;
    float bJetPtThreshold=30;

    std::vector<std::vector<float>> pfInfo;
    std::vector<std::vector<float>> svInfo; 
    float fjmsd, fjpt, fjrawpt, fjeta, fjphi;
    int NPFPROPS = 9, NSVPROPS = 13;

    GenJetInfo genJetInfo;
    int NGENPROPS = 8; 
    
    float minSoftTrackPt=0.3; // 300 MeV
};


/** templated functions **/

template <typename T> 
void PandaAnalyzer::MatchGenJets(T& genJets) 
{
  unsigned N = cleanedJets.size();
  for (unsigned i = 0; i != N; ++i) {
    panda::Jet *reco = cleanedJets.at(i);
    for (auto &gen : genJets) {
      if (DeltaR2(gen.eta(), gen.phi(), reco->eta(), reco->phi()) < 0.09) {
        gt->jetGenPt[i] = gen.pt();
        gt->jetGenFlavor[i] = gen.pdgid;
        break;
      }
    }
  }
  tr->TriggerEvent("match gen jets");
}

template <typename T>
void PandaAnalyzer::FillGenTree(panda::Collection<T>& genParticles)
{

  genJetInfo.pt = -1; genJetInfo.eta = -1; genJetInfo.phi = -1; genJetInfo.m = -1;
  genJetInfo.msd = -1;
  genJetInfo.tau3 = -1; genJetInfo.tau2 = -1; genJetInfo.tau1 = -1;
  genJetInfo.tau3sd = -1; genJetInfo.tau2sd = -1; genJetInfo.tau1sd = -1;
  genJetInfo.nprongs = -1;
  genJetInfo.partonpt = -1; genJetInfo.partonm = -1;
  gt->genFatJetPt = 0;
  for (unsigned i = 0; i != NMAXPF; ++i) {
    for (unsigned j = 0; j != NGENPROPS; ++j) {
      genJetInfo.particles[i][j] = 0;
    }
  }

  if (analysis->deepGenGrid)
    grid->clear();
  std::vector<int> leptonIndices; // hard leptons from t->W->lv decay
  std::vector<fastjet::PseudoJet> finalStates;
  unsigned idx = -1;
  for (auto &p : genParticles) {
    idx++;
    unsigned apdgid = abs(p.pdgid);
    if (!p.finalState)
      continue;
    if (apdgid == 12 ||
        apdgid == 14 ||
        apdgid == 16)
      continue; 
    if (p.pt() > 0.001) {
      if (!analysis->deepGenGrid || (pdgToQ[apdgid] != 0)) { // it's charged, so we have tracking
        finalStates.emplace_back(p.px(), p.py(), p.pz(), p.e());
        finalStates.back().set_user_index(idx);

        if (apdgid == 11 ||
            apdgid == 13 ||
            apdgid == 15) {
          const T *parent = &p;
          bool foundW = false, foundT = false;
          while (parent->parent.isValid()) {
            parent = parent->parent.get();
            unsigned parent_apdgid = abs(parent->pdgid);
            if (!foundW) {
              if (parent_apdgid == 24) {
                foundW = true;
                continue; 
              } else if (parent_apdgid != apdgid) {
                break; // if it's not a W, must be a parent of the particle we care about
              }
            } else {  // foundW = true
              if (parent_apdgid == 6) {
                foundT = true;
                break;
              } else if (parent_apdgid != 24) {
                break; // if it's not a top, must be a parent of the W we found
              }
            }
          }
          if (foundT) {
            leptonIndices.push_back(idx);
          }
        }
      } else {
        grid->add(p);
      }
    }
  }
  if (analysis->deepGenGrid) {
    int user_idx = -2;
    for (auto &v : grid->get()) {
      finalStates.emplace_back(v.Px(), v.Py(), v.Pz(), v.E());
      finalStates.back().set_user_index(user_idx); // not associated with a real particle
      --user_idx;
    } 
  }

  // cluster the  jet 
  fastjet::ClusterSequenceArea seq(finalStates, *jetDef, *areaDef);
  std::vector<fastjet::PseudoJet> allJets(fastjet::sorted_by_pt(seq.inclusive_jets(0.)));


  fastjet::PseudoJet *fullJet = NULL;
  for (auto testJet : allJets) {
    if (testJet.perp() < 450)
      break;
    bool jetOverlapsLepton = false;
    for (auto &c : testJet.constituents()) {
      int idx = c.user_index();
      if (std::find(leptonIndices.begin(),leptonIndices.end(),idx) != leptonIndices.end()) {
        jetOverlapsLepton = true;
        break;
      }
    }
    if (!jetOverlapsLepton) {
      fullJet = &testJet;
      break;
    }
  }

  if (fullJet == NULL) {
    tr->TriggerEvent("fill gen tree");
    return;
  }

  gt->genFatJetPt = fullJet->perp();

  if (gt->genFatJetPt < 450) {
    tr->TriggerEvent("fill gen tree");
    return;
  }

  VPseudoJet allConstituents = fastjet::sorted_by_pt(fullJet->constituents());
  genJetInfo.pt = gt->genFatJetPt;
  genJetInfo.m = fullJet->m();
  genJetInfo.eta = fullJet->eta();
  genJetInfo.phi = fullJet->phi();

  // softdrop the jet
  fastjet::PseudoJet sdJet = (*softDrop)(*fullJet);
  VPseudoJet sdConstituents = fastjet::sorted_by_pt(sdJet.constituents());
  genJetInfo.msd = sdJet.m();
  std::vector<bool> survived(allConstituents.size());
  unsigned nC = allConstituents.size();
  for (unsigned iC = 0; iC != nC; ++iC) {
    int idx = allConstituents.at(iC).user_index();
    survived[iC] = false;
    for (auto &sdc : sdConstituents) {
      if (idx == sdc.user_index()) {
        survived[iC] = true; 
        break;
      }
    }
  }

  // get tau  
  genJetInfo.tau1 = tauN->getTau(1, allConstituents);
  genJetInfo.tau2 = tauN->getTau(2, allConstituents);
  genJetInfo.tau3 = tauN->getTau(3, allConstituents);
  genJetInfo.tau1sd = tauN->getTau(1, sdConstituents);
  genJetInfo.tau2sd = tauN->getTau(2, sdConstituents);
  genJetInfo.tau3sd = tauN->getTau(3, sdConstituents);

  // now we have to count the number of prongs 
  float dR2 = FATJETMATCHDR2;
  float base_eta = genJetInfo.eta, base_phi = genJetInfo.phi;
  auto matchJet = [base_eta, base_phi, dR2](const T &p) -> bool {
    return DeltaR2(base_eta, base_phi, p.eta(), p.phi()) < dR2;
  };
  float threshold = genJetInfo.pt * 0.2;
  std::unordered_set<const T*> partons; 
  for (auto &gen : genParticles) {
    unsigned apdgid = abs(gen.pdgid);
    if (apdgid > 5 && 
        apdgid != 21 &&
        apdgid != 15 &&
        apdgid != 11 && 
        apdgid != 13)
      continue; 

    if (gen.pt() < threshold)
      continue; 

    if (!matchJet(gen))
      continue;

    const T *parent = &gen;
    const T *foundParent = NULL;
    while (parent->parent.isValid()) {
      parent = parent->parent.get();
      if (partons.find(parent) != partons.end()) {
        foundParent = parent;
        break;
      }
    }


    T *dau1 = NULL, *dau2 = NULL;
    for (auto &child : genParticles) {
      if (!(child.parent.isValid() && 
            child.parent.get() == &gen))
        continue; 
      
      unsigned child_apdgid = abs(child.pdgid);
      if (child_apdgid > 5 && 
          child_apdgid != 21 &&
          child_apdgid != 15 &&
          child_apdgid != 11 && 
          child_apdgid != 13)
        continue; 

      if (dau1)
        dau2 = &child;
      else
        dau1 = &child;

      if (dau1 && dau2)
        break;
    }

    if (dau1 && dau2 && 
        dau1->pt() > threshold && dau2->pt() > threshold && 
        matchJet(*dau1) && matchJet(*dau2)) {
      if (foundParent) {
        partons.erase(partons.find(foundParent));
      }
      partons.insert(dau1);
      partons.insert(dau2);
    } else if (foundParent) {
      continue; 
    } else {
      partons.insert(&gen);
    }
  }

  genJetInfo.nprongs = partons.size();
  gt->genFatJetNProngs = genJetInfo.nprongs;

  TLorentzVector vPartonSum;
  TLorentzVector vTmp;
  for (auto *p : partons) {
    vTmp.SetPtEtaPhiM(p->pt(), p->eta(), p->phi(), p->m());
    vPartonSum += vTmp;
  }
  genJetInfo.partonm = vPartonSum.M();
  genJetInfo.partonpt = vPartonSum.Pt();

  std::map<const T*, unsigned> partonToIdx;
  for (auto &parton : partons) 
    partonToIdx[parton] = partonToIdx.size(); // just some arbitrary ordering 

  // get the hardest particle with angle wrt jet axis > 0.1
  fastjet::PseudoJet *axis2 = NULL;
  for (auto &c : allConstituents) {
    if (DeltaR2(c.eta(), c.phi(), genJetInfo.eta, genJetInfo.phi) > 0.01) {
      if (!axis2 || (c.perp() > axis2->perp())) {
        axis2 = &c;
      }
    }
  }

  JetRotation rot(fullJet->px(), fullJet->py(), fullJet->pz(),
                  axis2->px(), axis2->py(), axis2->pz());

  // now we fill the particles
  nC = std::min(nC, (unsigned)NMAXPF);
  for (unsigned iC = 0; iC != nC; ++iC) {
    fastjet::PseudoJet &c = allConstituents.at(iC);

    if (c.perp() < 0.001) // not a real particle
      continue;

    // genJetInfo.particles[iC][0] = c.perp() / fullJet->perp();
    // genJetInfo.particles[iC][1] = c.eta() - fullJet->eta();
    // genJetInfo.particles[iC][2] = SignedDeltaPhi(c.phi(), fullJet->phi());
    // genJetInfo.particles[iC][3] = c.m();
    // genJetInfo.particles[iC][4] = c.e();
    float angle = DeltaR2(c.eta(), c.phi(), genJetInfo.eta, genJetInfo.phi);
    float x=c.px(), y=c.py(), z=c.pz();
    rot.Rotate(x, y, z);  // perform two rotations on the jet 
    genJetInfo.particles[iC][0] = x;
    genJetInfo.particles[iC][1] = y;
    genJetInfo.particles[iC][2] = z;
    genJetInfo.particles[iC][3] = c.e();
    genJetInfo.particles[iC][4] = angle;
    genJetInfo.particles[iC][5] = survived[iC] ? 1 : 0;

    unsigned ptype = 0;
    int parent_idx = -1;
    if (c.user_index() >= 0) {
      T &gen = genParticles.at(c.user_index());
      int pdgid = gen.pdgid;
      unsigned apdgid = abs(pdgid);
      if (apdgid == 11) {
        ptype = 1 * sign(pdgid * -11);
      } else if (apdgid == 13) {
        ptype = 2 * sign(pdgid * -13);
      } else if (apdgid == 22) {
        ptype = 3;
      } else {
        float q = pdgToQ[apdgid];
        if (apdgid != pdgid)
          q *= -1;
        if (q == 0) 
          ptype = 4;
        else if (q > 0) 
          ptype = 5;
        else 
          ptype = 6;
      }

      const T *parent = &gen;
      while (parent->parent.isValid()) {
        parent = parent->parent.get();
        if (partons.find(parent) != partons.end()) {
          parent_idx = partonToIdx[parent];
          break;
        }
      }
    }

    genJetInfo.particles[iC][6] = ptype;
    genJetInfo.particles[iC][7] = parent_idx;
  }

  tr->TriggerEvent("fill gen tree");
}


#endif

