#ifndef Scout4B_KinematicFitResult_h
#define Scout4B_KinematicFitResult_h

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"

std::pair<float, float> getAlpha(const GlobalPoint &vtx_position, const GlobalError &vtx_error,
                                 const GlobalPoint &ip_position, const GlobalError &ip_error,
                                 const GlobalVector &momentum,
                                 bool transverse = true);

std::pair<float, float> getAlphaXY(const GlobalPoint &vtx_position, const GlobalError &vtx_error,
                                   const GlobalPoint &ip_position, const GlobalError &ip_error,
                                   const GlobalVector &momentum,
                                   bool transverse = true);

    typedef reco::Candidate::LorentzVector LorentzVector;

class KinematicFitResult
{
public:
  KinematicFitResult() : charge_(-999), treeIsValid(false),
                         lxy_(-1.0), lxyErr_(-1.0), sigLxy_(-1.0),
                         lz_(-1.0), lzErr_(-1.0), sigLz_(-1.0),
                         l_(-1.0), lErr_(-1.0), sigL_(-1.0),
                         alphaBS_(-999.), alphaBSErr_(-999.),
                         alphaBSXY_(-999.), alphaBSXYErr_(-999.)
  {
  }

  void set_tree(RefCountedKinematicTree tree);

  bool valid() const;

  // compute displacement with respect to the beam spot
  void postprocess(const reco::BeamSpot &beamSpot);

  float mass() const;
  float refit_mass(unsigned int i, unsigned int j) const;
  GlobalVector p3() const;
  LorentzVector p4() const;
  unsigned int number_of_daughters() const
  {
    return refitDaughters.size();
  }
  GlobalVector dau_p3(unsigned int i) const;
  LorentzVector dau_p4(unsigned int i) const;
  float massErr() const;
  float chi2() const;
  float ndof() const;
  float vtxProb() const;
  GlobalPoint vtx_position() const;
  VertexState vtx_state() const;
  GlobalError vtx_error() const;
  reco::Vertex vertex() const;
  float sumPt() const;
  float sumPt2() const;
  float lxy() const { return lxy_; }
  float lxyErr() const { return lxyErr_; }
  float sigLxy() const { return sigLxy_; }
  float lz() const { return lz_; }
  float lzErr() const { return lzErr_; }
  float sigLz() const { return sigLz_; }
  float l() const { return l_; }
  float lErr() const { return lErr_; }
  float sigL() const { return sigL_; }
  float alphaBS() const { return alphaBS_; }
  float alphaBSErr() const { return alphaBSErr_; }
  float alphaBSXY() const { return alphaBSXY_; }
  float alphaBSXYErr() const { return alphaBSXYErr_; }

  RefCountedKinematicTree tree() { return refitTree; }
  const RefCountedKinematicParticle particle() const { return refitMother; }

  int charge_;
  std::vector<int> dau_charge_;

  std::vector<reco::Track *> tracks;
  std::vector<unsigned int> trackIndices;

  bool treeIsValid;
  RefCountedKinematicVertex refitVertex;

private:
  float lxy_, lxyErr_, sigLxy_;
  float lz_, lzErr_, sigLz_;
  float l_, lErr_, sigL_;
  float alphaBS_, alphaBSErr_;
  float alphaBSXY_, alphaBSXYErr_;
  RefCountedKinematicParticle refitMother;
  RefCountedKinematicTree refitTree;
  std::vector<RefCountedKinematicParticle> refitDaughters;
};

#endif
