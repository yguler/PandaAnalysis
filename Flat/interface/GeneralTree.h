#ifndef GENERALTREE_H
#define GENERALTREE_H 

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TString.h"
#include "genericTree.h"
#include <map>

#define NJET 20
#define NLEP 4
#define NSUBJET 2

class GeneralTree : public genericTree {
    public:
      // public objects
      struct ECFParams {
        int ibeta;
        int N;
        int order;
        bool operator==(const ECFParams &o) const {
          return ibeta==o.ibeta && N==o.N && order==o.order;
        }
        bool operator<(const ECFParams &o) const {
          return ( N<o.N ||
                   (N==o.N && order<o.order) ||
                   (N==o.N && order==o.order && ibeta<o.ibeta) );
        }
        bool operator>(const ECFParams &o) const {
          return ! operator<(o);
        }
      };

      enum BTagShift {
        bCent=0,
        bBUp,
        bBDown,
        bMUp,
        bMDown,
        bNShift
      };
      enum BTagJet {
        bJet=0,
	bMedJet,
        bSubJet,
        bNJet
      };
      enum BTagTags {
        b0=0,
        b1,
        b2,
        bGT0,
        bNTags
      };
      struct BTagParams {
        BTagJet jet;
        BTagTags tag;
        BTagShift shift=(BTagShift)0;
        bool operator==(const BTagParams &o) const {
          return jet==o.jet && tag==o.tag && shift==o.shift;
        }
        bool operator<(const BTagParams &o) const {
          return ( jet<o.jet ||
                   (jet==o.jet && tag<o.tag) ||
                   (jet==o.jet && tag==o.tag && shift<o.shift) ); 
        }
        bool operator>(const BTagParams &o) const {
          return ! operator<(o);
        }
      };
      enum csvShift {
        csvCent=0,
        csvJESup=7,
        csvJESdown=8,
        csvLFup=9,
        csvLFdown=10,
        csvHFup=11,
        csvHFdown=12,
        csvHFStats1up=13,
        csvHFStats1down=14,
        csvHFStats2up=15,
        csvHFStats2down=16,
        csvLFStats1up=17,
        csvLFStats1down=18,
        csvLFStats2up=19,
        csvLFStats2down=20,
        csvCErr1up=21,
        csvCErr1down=22,
        csvCErr2up=23,
        csvCErr2down=24
      };
      // Array of the CSV/CMVA weight enums that can be looped over
      static const unsigned char nCsvShifts=19;
      csvShift csvShifts[nCsvShifts] = {
        csvCent,
        csvJESup,
        csvJESdown,
        csvLFup,
        csvLFdown,
        csvHFup,
        csvHFdown,
        csvHFStats1up,
        csvHFStats1down,
        csvHFStats2up,
        csvHFStats2down,
        csvLFStats1up,
        csvLFStats1down,
        csvLFStats2up,
        csvLFStats2down,
        csvCErr1up,
        csvCErr1down,
        csvCErr2up,
        csvCErr2down
      };

    
      virtual void SetAuxTree(TTree *t);  

    private:
        std::vector<double> betas = {0.5, 1.0, 2.0, 4.0};
        std::vector<int> ibetas = {0,1,2,3};
        std::vector<int> Ns = {1,2,3,4}; 
        std::vector<int> orders = {1,2,3};
        std::vector<ECFParams> ecfParams;
        std::vector<BTagParams> btagParams;

        TString makeECFString(ECFParams p) {
            return TString::Format("ECFN_%i_%i_%.2i",p.order,p.N,int(10*betas.at(p.ibeta)));
        }
        TString makeBTagSFString(BTagParams p) {
          TString s = "sf_";
          if (p.jet==bSubJet)
            s += "sj";
	  if (p.jet==bMedJet)
	    s += "Med";
          s += "btag";
          switch (p.tag) {
            case b0:
              s += "0"; break;
            case b1:
              s += "1"; break;
            case b2:
              s += "2"; break;
            case bGT0:
              s += "GT0"; break;
            default:
              break;
          }
          if (p.shift==bCent)
            return s;
          switch (p.shift) {
            case bBUp:
              s += "BUp"; break;
            case bBDown:
              s += "BDown"; break;
            case bMUp:
              s += "MUp"; break;
            case bMDown:
              s += "MDown"; break;
            default: break;
          }
          return s;
        }
        TString makeCsvWeightString(csvShift shift, bool isCMVA=false) { 
          TString s = isCMVA? "sf_cmvaWeight_" : "sf_csvWeight_";
          switch (shift) {
            case csvCent         :  s += "Cent"         ; break;
            case csvJESup        :  s += "JESup"        ; break;
            case csvJESdown      :  s += "JESdown"      ; break;
            case csvLFup         :  s += "LFup"         ; break;
            case csvLFdown       :  s += "LFdown"       ; break;
            case csvHFup         :  s += "HFup"         ; break;
            case csvHFdown       :  s += "HFdown"       ; break;
            case csvHFStats1up   :  s += "HFStats1up"   ; break;
            case csvHFStats1down :  s += "HFStats1down" ; break;
            case csvHFStats2up   :  s += "HFStats2up"   ; break;
            case csvHFStats2down :  s += "HFStats2down" ; break;
            case csvLFStats1up   :  s += "LFStats1up"   ; break;
            case csvLFStats1down :  s += "LFStats1down" ; break;
            case csvLFStats2up   :  s += "LFStats2up"   ; break;
            case csvLFStats2down :  s += "LFStats2down" ; break;
            case csvCErr1up      :  s += "CErr1up"      ; break;
            case csvCErr1down    :  s += "CErr1down"    ; break;
            case csvCErr2up      :  s += "CErr2up"      ; break;
            case csvCErr2down    :  s += "CErr2down"    ; break;
            default              :  s += "Unknown"      ; break;
          }
          return s;
        }
    public:
      GeneralTree();
      ~GeneralTree();
      void WriteTree(TTree *t);
      void Fill() { treePtr->Fill(); }
      void SetBranchStatus(const char *bname, bool status, UInt_t *ret=0) 
      { 
        treePtr->SetBranchStatus(bname,status,ret); 
      }
      void Reset();

      std::vector<double> get_betas() const { return betas; }
      std::vector<int> get_ibetas() const { return ibetas; }
      std::vector<int> get_Ns() const { return Ns; }
      std::vector<int> get_orders() const { return orders; }
        
      // public config
      bool monohiggs=false, monojet=false,  vbf=false, fatjet=true, leptonic=false, photonic=false, hfCounting=false;
      bool btagWeights=false, useCMVA=false;

//STARTCUSTOMDEF
      std::map<ECFParams,float> fj1ECFNs;
      std::map<BTagParams,float> sf_btags;
      std::map<BTagParams,float> sf_alt_btags;
      std::map<TString,float> signal_weights;
      std::map<csvShift,float> sf_csvWeights; // this is called csvWeights, but may actually include CMVA weights instead

      float fj1sjPt[NSUBJET];
      float fj1sjEta[NSUBJET];
      float fj1sjPhi[NSUBJET];
      float fj1sjM[NSUBJET];
      float fj1sjCSV[NSUBJET];
      float fj1sjQGL[NSUBJET];

      float jetPt[NJET];
      float jetPtUp[NJET];
      float jetPtDown[NJET];
      float jetEta[NJET];
      float jetPhi[NJET];
      float jetM[NJET];
      float jetE[NJET];
      float jetCSV[NJET];
      float jetCMVA[NJET];
      float jetIso[NJET];
      float jetQGL[NJET];
      float jetLeadingLepPt[NJET];
      float jetLeadingLepPtRel[NJET];
      float jetLeadingLepDeltaR[NJET];
      float jetLeadingTrkPt[NJET];
      float jetvtxPt[NJET];
      float jetvtxMass[NJET];
      float jetvtx3Dval[NJET];
      float jetvtx3Derr[NJET];
      int jetvtxNtrk[NJET];
      float jetEMFrac[NJET];
      float jetHadFrac[NJET];
      int jetNLep[NJET];
      float jetGenPt[NJET];
      int jetGenFlavor[NJET];

      int bosonjtidx[2];
      float jetRegFac[2];

      float scale[6];
      
      float muonPt[NLEP];
      float muonEta[NLEP];
      float muonPhi[NLEP];
      float muonD0[NLEP];
      float muonDZ[NLEP];
      float muonSfLoose[NLEP];
      float muonSfMedium[NLEP];
      float muonSfTight[NLEP];
      float muonSfUnc[NLEP];
      float muonSfReco[NLEP];
      int muonSelBit[NLEP];
      int muonPdgId[NLEP];
      int muonIsSoftMuon[NLEP];
      float muonCombIso[NLEP];
      //int muonIsGlobalMuon[NLEP];
      //int muonIsTrackerMuon[NLEP];
      //int muonNValidMuon[NLEP];
      //int muonNValidPixel[NLEP];
      //int muonTrkLayersWithMmt[NLEP];
      //int muonPixLayersWithMmt[NLEP];
      //int muonNMatched[NLEP];
      //int muonChi2LocalPosition[NLEP];
      //int muonTrkKink[NLEP];
      //float muonValidFraction[NLEP];
      //float muonNormChi2[NLEP];
      //float muonSegmentCompatibility[NLEP];

      float electronPt[NLEP];
      float electronEta[NLEP];
      float electronPhi[NLEP];
      float electronD0[NLEP];
      float electronDZ[NLEP];
      float electronSfLoose[NLEP];
      float electronSfMedium[NLEP];
      float electronSfTight[NLEP];
      float electronSfMvaWP90[NLEP];
      float electronSfMvaWP80[NLEP];
      float electronSfUnc[NLEP];
      float electronSfReco[NLEP];
      int electronSelBit[NLEP];
      int electronPdgId[NLEP];
      //float electronChIsoPh[NLEP];
      //float electronNhIsoPh[NLEP];
      //float electronPhIsoPh[NLEP];
      //float electronEcalIso[NLEP];
      //float electronHcalIso[NLEP];
      //float electronTrackIso[NLEP];
      //float electronIsoPUOffset[NLEP];
      //float electronSieie[NLEP];
      //float electronSipip[NLEP];
      //float electronDEtaInSeed[NLEP];
      //float electronDPhiIn[NLEP];
      //float electronEseed[NLEP];
      //float electronHOverE[NLEP];
      //float electronEcalE[NLEP];
      //float electronTrackP[NLEP];
      int electronNMissingHits[NLEP];
      int electronTripleCharge[NLEP];
      float electronCombIso[NLEP];

//ENDCUSTOMDEF
    float jot1PhiUp = -1;
    float jot1PhiDown = -1;
    float jot2PhiUp = -1;
    float jot2PhiDown = -1;
    int loosePho1SelBit = -1;
    int looseGenPho1PdgId = -1;
    int genFatJetNProngs = -1;
    float genFatJetPt = -1;
    int fj1NBPartons = -1;
    int fj1NCPartons = -1;
    float fj1Rho2 = -1;
    float fj1RawRho2 = -1;
    float fj1Rho = -1;
    float fj1RawRho = -1;
    int fj1NPartons = -1;
    float fj1PartonM = -1;
    float fj1PartonPt = -1;
    float fj1PartonEta = -1;
    float trkmetphi = -1;
    float sf_zzUnc = -1;
    float sf_zz = -1;
    float sf_wz = -1;
    float sf_vh = -1;
    float sf_vhUp = -1;
    float sf_vhDown = -1;
    float genLep1Pt = -1;
    float genLep1Eta = -1;
    float genLep1Phi = -1;
    int genLep1PdgId = -1;
    float genLep2Pt = -1;
    float genLep2Eta = -1;
    float genLep2Phi = -1;
    int genLep2PdgId = -1;
    float genLep3Pt = -1;
    float genLep3Eta = -1;
    float genLep3Phi = -1;
    int genLep3PdgId = -1;
    float genLep4Pt = -1;
    float genLep4Eta = -1;
    float genLep4Phi = -1;
    int genLep4PdgId = -1;
    int looseGenLep1PdgId = -1;
    int looseGenLep2PdgId = -1;
    int looseGenLep3PdgId = -1;
    int looseGenLep4PdgId = -1;
    int looseLep1PdgId = -1;
    int looseLep2PdgId = -1;
    int whichRecoil = -1;
    float genJet1Pt = -1;
    float genJet2Pt = -1;
    float genJet1Eta = -1;
    float genJet2Eta = -1;
    float genMjj = -1;
    float sf_qcdV_VBF2lTight = -1;
    float sf_qcdV_VBF2l = -1;
    float barrelHTMiss = -1;
    float barrelJet12Pt = -1;
    float barrelJet1Pt = -1;
    float barrelJet1Eta = -1;
    float barrelHT = -1;
    float sf_puUp = -1;
    float sf_puDown = -1;
    float genMuonPt = -1;
    float genMuonEta = -1;
    float genElectronPt = -1;
    float genElectronEta = -1;
    float genTauPt = -1;
    float genTauEta = -1;
    int badECALFilter = -1;
    float sf_qcdV_VBFTight = -1;
    float sf_metTrigVBF = -1;
    float sf_metTrigZmmVBF = -1;
    float sumETRaw = -1;
    int jot1VBFID = -1;
    float sf_metTrigZmm = -1;
    float sf_qcdV_VBF = -1;
    int jetNMBtags = -1;
    float pfmetRaw = -1;
    int nAK8jet = -1;
    float ak81Pt = -1;
    float ak81Eta = -1;
    float ak81Phi = -1;
    float ak81MaxCSV = -1;
    int nB = -1;
    int nBGenJets = -1;
    float fj1MSDScaleUp_sj = -1;
    float fj1MSDScaleDown_sj = -1;
    float fj1MSDSmeared_sj = -1;
    float fj1MSDSmearedUp_sj = -1;
    float fj1MSDSmearedDown_sj = -1;
    float fj1PtScaleUp_sj = -1;
    float fj1PtScaleDown_sj = -1;
    float fj1PtSmeared_sj = -1;
    float fj1PtSmearedUp_sj = -1;
    float fj1PtSmearedDown_sj = -1;
    float jot2EtaUp = -1;
    float jot2EtaDown = -1;
    float jot1EtaUp = -1;
    float jot1EtaDown = -1;
    float jot1PtUp = -1;
    float jot1PtDown = -1;
    float jot2PtUp = -1;
    float jot2PtDown = -1;
    float jot12MassUp = -1;
    float jot12DEtaUp = -1;
    float jot12DPhiUp = -1;
    float jot12MassDown = -1;
    float jot12DEtaDown = -1;
    float jot12DPhiDown = -1;
    float pfmetUp = -1;
    float pfmetDown = -1;
    float pfUWmagUp = -1;
    float pfUZmagUp = -1;
    float pfUAmagUp = -1;
    float pfUWWmagUp = -1;
    float pfUmagUp = -1;
    float pfUWmagDown = -1;
    float pfUZmagDown = -1;
    float pfUAmagDown = -1;
    float pfUWWmagDown = -1;
    float pfUmagDown = -1;
    int nJot = -1;
    int nJot_jesUp = -1;
    int nJot_jesDown = -1;
    float jot1Phi = -1;
    float jot1Pt = -1;
    float jot1GenPt = -1;
    float jot1Eta = -1;
    float jot2Phi = -1;
    float jot2Pt = -1;
    float jot2GenPt = -1;
    float jot2Eta = -1;
    float jot12DPhi = -1;
    int isGS = -1;
    float fj1SubMaxCSV = -1;
    int runNumber = -1;
    int lumiNumber = -1;
    ULong64_t eventNumber = -1;
    int npv = -1;
    int pu = -1;
    float mcWeight = -1;
    int trigger = -1;
    int metFilter = -1;
    int egmFilter = -1;
    float filter_maxRecoil = -1;
    float filter_whichRecoil = -1;
    float sf_ewkV = -1;
    float sf_qcdV = -1;
    float sf_ewkV2j = -1;
    float sf_qcdV2j = -1;
    float sf_qcdTT = -1;
    float sf_lepID = -1;
    float sf_lepIso = -1;
    float sf_lepTrack = -1;
    float sf_pho = -1;
    float sf_eleTrig = -1;
    float sf_muTrig = -1;
    float sf_phoTrig = -1;
    float sf_metTrig = -1;
    float sf_pu = -1;
    float sf_npv = -1;
    float sf_tt = -1;
    float sf_tt_ext = -1;
    float sf_tt_bound = -1;
    float sf_tt8TeV = -1;
    float sf_tt8TeV_ext = -1;
    float sf_tt8TeV_bound = -1;
    float sf_phoPurity = -1;
    float pfmet = -1;
    float pfmetphi = -1;
    float pfmetnomu = -1;
    float pfmetsig = -1;
    float puppimet = -1;
    float puppimetphi = -1;
    float puppimetsig = -1;
    float calomet = -1;
    float calometphi = -1;
    float pfcalobalance = -1;
    float sumET = -1;
    float trkmet = -1;
    float puppiUWmag = -1;
    float puppiUWphi = -1;
    float puppiUZmag = -1;
    float puppiUZphi = -1;
    float puppiUAmag = -1;
    float puppiUAphi = -1;
    float puppiUperp = -1;
    float puppiUpara = -1;
    float puppiUmag = -1;
    float puppiUphi = -1;
    float puppiUWWmag = -1;
    float puppiUWWphi = -1;
    float pfUWmag = -1;
    float pfUWphi = -1;
    float pfUZmag = -1;
    float pfUZphi = -1;
    float pfUAmag = -1;
    float pfUAphi = -1;
    float pfUWWmag = -1;
    float pfUWWphi = -1;
    float pfUperp = -1;
    float pfUpara = -1;
    float pfUmag = -1;
    float pfUphi = -1;
    float dphipfmet = -1;
    float dphipuppimet = -1;
    float dphipuppiUW = -1;
    float dphipuppiUWW = -1;
    float dphipuppiUZ = -1;
    float dphipuppiUA = -1;
    float dphipfUW = -1;
    float dphipfUZ = -1;
    float dphipfUA = -1;
    float dphipfUWW = -1;
    float dphipuppiU = -1;
    float dphipfU = -1;
    float trueGenBosonPt = -1;
    float genBosonPt = -1;
    float genBosonEta = -1;
    float genBosonMass = -1;
    float genBosonPhi = -1;

    float genBosonGenPx = -1;
    float genBosonGenPy = -1;
    float genBosonGenPz = -1;
    float genBosonGenEn = -1;
    int genBosonGenFlav = -1;
    float genBosonGenMass = -1;

    float genWPlusPt = -1;
    float genWMinusPt = -1;
    float genWPlusEta = -1;
    float genWMinusEta = -1;
    float genTopPt = -1;
    int genTopIsHad = -1;
    float genTopEta = -1;
    float genAntiTopPt = -1;
    int genAntiTopIsHad = -1;
    float genAntiTopEta = -1;
    float genTTPt = -1;
    float genTTEta = -1;
    int nJet = -1;
    int nJet_jesUp = -1;
    int nJet_jesDown = -1;
    int nIsoJet = -1;
    int jet1Flav = -1;
    float jet1Phi = -1;
    float jet1Pt = -1;
    float jet1GenPt = -1;
    float jet1Eta = -1;
    float jet1CSV = -1;
    float jet1CMVA = -1;
    int jet1IsTight = -1;
    int jet2Flav = -1;
    float jet2Phi = -1;
    float jet2Pt = -1;
    float jet2GenPt = -1;
    float jet2Eta = -1;
    float jet2CSV = -1;
    float jet2CMVA = -1;
    float jet2EtaUp = -1;
    float jet2EtaDown = -1;
    float jet1EtaUp = -1;
    float jet1EtaDown = -1;
    float jet1PtUp = -1;
    float jet1PtDown = -1;
    float jet2PtUp = -1;
    float jet2PtDown = -1;
    int jet3Flav = -1;
    float jet3Phi = -1;
    float jet3Pt = -1;
    float jet3GenPt = -1;
    float jet3Eta = -1;
    float jet3CSV = -1;
    int jet4Flav = -1;
    float jet4Phi = -1;
    float jet4Pt = -1;
    float jet4GenPt = -1;
    float jet4Eta = -1;
    float jet4CSV = -1;
    int jet5Flav = -1;
    float jet5Phi = -1;
    float jet5Pt = -1;
    float jet5GenPt = -1;
    float jet5Eta = -1;
    float jet5CSV = -1;
    float forwjet1Pt = -1;
    float forwjet1Phi = -1;
    float forwjet1Eta = -1;
    //float forwjet1GenPt = -1;
    float forwjet2Pt = -1;
    float forwjet2Phi =-1;
    float forwjet2Eta =-1;
    //float forwjet2GenPt= -1;
    //Gen Info
    float jet1GenPx=0;
    float jet1GenPy=0;
    float jet1GenPz=0;
    float jet1GenEn=0;
    int jet1GenStatus=0;
    int jet1GenFlav=0;
    float jet1GenMass=0;

    float jet2GenPx=0;
    float jet2GenPy=0;
    float jet2GenPz=0;
    float jet2GenEn=0;
    int jet2GenStatus=0;
    int jet2GenFlav=0;
    float jet2GenMass=0;

    float jet3GenPx=0;
    float jet3GenPy=0;
    float jet3GenPz=0;
    float jet3GenEn=0;
    int jet3GenStatus=0;
    int jet3GenFlav=0;
    float jet3GenMass=0;

    float jet4GenPx=0;
    float jet4GenPy=0;
    float jet4GenPz=0;
    float jet4GenEn=0;
    int jet4GenStatus=0;
    int jet4GenFlav=0;
    float jet4GenMass=0;

    float jet5GenPx=0;
    float jet5GenPy=0;
    float jet5GenPz=0;
    float jet5GenEn=0;
    int jet5GenStatus=0;
    int jet5GenFlav=0;
    float jet5GenMass=0;

    float isojet1Pt = -1;
    float isojet1CSV = -1;
    int isojet1Flav = -1;
    float isojet2Pt = -1;
    float isojet2CSV = -1;
    int isojet2Flav = -1;
    float jot12Mass = -1;
    float jot12DEta = -1;
    int jetNBtags = -1;
    int isojetNBtags = -1;
    int nFatjet = -1;
    float fj1Tau32 = -1;
    float fj1Tau21 = -1;
    float fj1Tau32SD = -1;
    float fj1Tau21SD = -1;
    float fj1MSD = -1;
    float fj1MSDScaleUp = -1;
    float fj1MSDScaleDown = -1;
    float fj1MSDSmeared = -1;
    float fj1MSDSmearedUp = -1;
    float fj1MSDSmearedDown = -1;
    float fj1MSD_corr = -1;
    float fj1Pt = -1;
    float fj1PtScaleUp = -1;
    float fj1PtScaleDown = -1;
    float fj1PtSmeared = -1;
    float fj1PtSmearedUp = -1;
    float fj1PtSmearedDown = -1;
    float fj1Phi = -1;
    float fj1Eta = -1;
    float fj1M = -1;
    float fj1MaxCSV = -1;
    float fj1MinCSV = -1;
    float fj1DoubleCSV = -1;
    int fj1Nbs = -1;
    int fj1gbb = -1;
    float fj1GenPt = -1;
    float fj1GenSize = -1;
    int fj1IsMatched = -1;
    float fj1GenWPt = -1;
    float fj1GenWSize = -1;
    int fj1IsWMatched = -1;
    int fj1HighestPtGen = -1;
    float fj1HighestPtGenPt = -1;
    int fj1IsTight = -1;
    int fj1IsLoose = -1;
    float fj1RawPt = -1;
    int fj1NHF = -1;
    float fj1HTTMass = -1;
    float fj1HTTFRec = -1;
    int fj1IsClean = -1;
    int fj1NConst = -1;
    int fj1NSDConst = -1;
    float fj1EFrac100 = -1;
    float fj1SDEFrac100 = -1;
    int nHF = -1;
    int nLoosePhoton = -1;
    int nTightPhoton = -1;
    int loosePho1IsTight = -1;
    float loosePho1Pt = -1;
    float loosePho1Eta = -1;
    float loosePho1Phi = -1;
    int nLooseLep = -1;
    int nLooseElectron = -1;
    int nLooseMuon = -1;
    int nTightLep = -1;
    int nTightElectron = -1;
    int nTightMuon = -1;
    float diLepMass = -1;
    int nTau = -1;
    float mT = -1;
    float bosonpt = -1;
    float bosoneta = -1;
    float bosonphi = -1;
    float bosonm = -1;
    float bosonm_reg = -1;
    float bosonpt_reg = -1;
    float bosonpt_jesUp = -1;
    float bosoneta_jesUp = -1;
    float bosonphi_jesUp = -1;
    float bosonm_jesUp = -1;
    float bosonm_reg_jesUp = -1;
    float bosonpt_reg_jesUp = -1;
    float bosonpt_jesDown = -1;
    float bosoneta_jesDown = -1;
    float bosonphi_jesDown = -1;
    float bosonm_jesDown = -1;
    float bosonm_reg_jesDown = -1;
    float bosonpt_reg_jesDown = -1;
    float bosonCosThetaJJ = -1;
    float bosonCosThetaCSJ1 = -1;
    float topMassLep1Met = -1;
    float topMassLep1Met_jesUp = -1;
    float topMassLep1Met_jesDown = -1;
    float topWBosonCosThetaCS = -1;
    float topWBosonPt = -1;
    float topWBosonEta = -1;
    float topWBosonPhi = -1;
    float sumEtSoft1 = -1;
    int nSoft2 = -1;
    int nSoft5 = -1;
    int nSoft10 = -1;
    float scaleUp = -1;
    float scaleDown = -1;
    float pdfUp = -1;
    float pdfDown = -1;
};

#endif
