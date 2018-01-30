#include "PandaAnalysis/Flat/interface/AnalyzerUtilities.h"

////////////////////////////////////////////////////////////////////////////////////

RotationToZ::RotationToZ(float x, float y, float z)
{
  TVector3 v_(x, y, z);
  TVector3 z_(0, 0, 1);
  r.Rotate(v_.Angle(z_), v_.Cross(z_)); 
}

void RotationToZ::Rotate(float& x, float& y, float& z) 
{
  TVector3 v_(x, y, z);
  v_ = r * v_;
  x = v_.x(); y = v_.y(); z = v_.z();
}

////////////////////////////////////////////////////////////////////////////////////

VPseudoJet ConvertPFCands(std::vector<const panda::PFCand*> &incoll, bool puppi, double minPt) 
{
  VPseudoJet vpj;
  vpj.reserve(incoll.size());
  int idx = -1;
  for (auto *incand : incoll) {
    double factor = puppi ? incand->puppiW() : 1;
    idx++;
    if (factor*incand->pt()<minPt)
      continue;
    vpj.emplace_back(factor*incand->px(),factor*incand->py(),
                     factor*incand->pz(),factor*incand->e());
    vpj.back().set_user_index(idx);
  }
  return vpj;
}

VPseudoJet ConvertPFCands(panda::RefVector<panda::PFCand> &incoll, bool puppi, double minPt) 
{
  std::vector<const panda::PFCand*> outcoll;
  outcoll.reserve(incoll.size());
  for (auto incand : incoll)
    outcoll.push_back(incand.get());

  return ConvertPFCands(outcoll, puppi, minPt);
}

VPseudoJet ConvertPFCands(panda::PFCandCollection &incoll, bool puppi, double minPt) 
{
  std::vector<const panda::PFCand*> outcoll;
  outcoll.reserve(incoll.size());
  for (auto &incand : incoll)
    outcoll.push_back(&incand);

  return ConvertPFCands(outcoll, puppi, minPt);
}

////////////////////////////////////////////////////////////////////////////////////

double TTNLOToNNLO(double pt) 
{
    double a = 0.1102;
    double b = 0.1566;
    double c = -3.685e-4;
    double d = 1.098;

    return TMath::Min(1.25,
                        a*TMath::Exp(-b*pow(pt,2)+1) + c*pt + d);
}

bool ElectronIP(double eta, double dxy, double dz) 
{
  double aeta = fabs(eta);
  if (aeta<1.4442) {
    return (dxy < 0.05 && dz < 0.10) ;
  } else {
    return (dxy < 0.10 && dz < 0.20);
  }
}

bool MuonIP(double dxy, double dz) 
{
  return (dxy < 0.02 && dz < 0.10);
}

bool IsMatched(std::vector<panda::Particle*>*objects,
               double deltaR2, double eta, double phi) 
{
  for (auto *x : *objects) {
    if (x->pt()>0) {
      if ( DeltaR2(x->eta(),x->phi(),eta,phi) < deltaR2 )
        return true;
    }
  }
  return false;
}

