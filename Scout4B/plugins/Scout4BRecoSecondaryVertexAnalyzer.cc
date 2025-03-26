/*
 * Developed by Yiyang Zhao for Run-3 B scouting
 * 2025-03
 * Scout4BRecoSecondaryVertexAnalyzer can produce secondary vertex any number of muon Track or Tracks
 * and fit multi secondary vertex to one X particle
 */

#include <memory>
#include <vector>
#include <cmath>
#include <utility>
#include <string>

#include "TMath.h"
#include "TVector.h"
#include "TMatrix.h"
#include "algorithm"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2F.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/OrphanHandle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/DeepCopyPointer.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "DataFormats/CLHEP/interface/Migration.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"

#include "Scout4B/Scout4B/interface/Displacement.h"
#include "Scout4B/Scout4B/interface/KinematicFitResult.h"
#include "Scout4B/Scout4B/interface/KinFitUtils.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace muon;
using namespace trigger;

//
// struct definition
//
struct XFitResult
{
   std::vector<std::pair<std::vector<TLorentzVector>, std::vector<int>>> mDaughter;
   std::vector<std::pair<TLorentzVector, int>> mMother;
   std::vector<double> mMassError;
   TLorentzVector xMotherP4;
   int xMotherCharge;
   double xMassError;
   double normChi2;

   double xlxy;
   double xlxyErr;
   double xlz;
   double xlzErr;
   double xl;
   double xlErr;
   double xsigLxy;
   double xsigLz;
   double xsigL;
   double xalphaBS;
   double xalphaBSErr;
   double xalphaBSXY;
   double xalphaBSXYErr;

   std::vector<unsigned int> mIndices;
   std::vector<std::vector<unsigned int>> mTrackIndices;

   XFitResult()
   {
      xMotherP4.SetXYZM(0, 0, 0, 0);
      xMotherCharge = 0;
      xMassError = 0;
      normChi2 = 0;
      xlxy = -999;
      xlxyErr = -999;
      xlz = -999;
      xlzErr = -999;
      xl = -999;
      xlErr = -999;
      xsigLxy = -999;
      xsigLz = -999;
      xsigL = -999;
      xalphaBS = -999;
      xalphaBSErr = -999;
      xalphaBSXY = -999;
      xalphaBSXYErr = -999;
   }
};

///////////////////////////////////////////////////////////////////////////
///                             P L U G I N
///////////////////////////////////////////////////////////////////////////

class Scout4BRecoSecondaryVertexAnalyzer : public edm::one::EDAnalyzer<>
{
public:
   explicit Scout4BRecoSecondaryVertexAnalyzer(const edm::ParameterSet &iConfig);
   ~Scout4BRecoSecondaryVertexAnalyzer() override;

   static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
   void beginJob() override;
   void analyze(edm::Event const &iEvent, edm::EventSetup const &iSetup) override;
   void endJob() override;

   void clearVars();

private:
   UInt_t getTriggerBits(const edm::Event &);

   float distanceOfClosestApproach(const reco::Track *track1, const reco::Track *track2, const MagneticField &bFieldHandle);
   Measurement1D distanceOfClosestApproach(const reco::Track *track, const reco::Vertex &vertex, const MagneticField &bFieldHandle);

   KinematicFitResult KinematicFitter(std::vector<reco::Track *> trks, std::vector<double> masses, const MagneticField &bFieldHandle);
   KinematicFitResult MCKinematicFitter(std::vector<reco::Track *> trks, std::vector<double> masses, double MMass, double MMassErr, const MagneticField &bFieldHandle);

   Scout4B::Displacements compute3dDisplacement(const KinematicFitResult &fit, bool closestIn3D = true);

   // Difine core functions for sv fitting
   std::pair<std::vector<KinematicFitResult>, std::vector<KinematicFitResult>> performVertexFit(
       std::vector<reco::Track> tracks,
       unsigned int NP,
       int MCharge,
       double MMass,
       double MMassErr,
       double MassWin,
       double vProbMin,
       const std::vector<double> &PMass,
       const std::vector<double> &PMassErr,
       const std::vector<int> &PCharge,
       double MMassMin,
       bool doIso,
       double PIso,
       const MagneticField &bFieldHandle,
       const unsigned long long max_loop);

   std::pair<std::vector<XFitResult>, std::vector<KinematicFitResult>> performXVertexFit(
       std::vector<reco::Track> tracks,
       unsigned int NM,
       const std::vector<std::vector<KinematicFitResult>> &vertexs,
       int XCharge,
       double XMass,
       double XMassWin,
       double vProbMin,
       const std::vector<double> &MMass,
       const std::vector<double> &MMassErr,
       double XMassMin,
       bool doIso,
       double PIso,
       const MagneticField &bFieldHandle,
       const unsigned long long max_loop);

   // Define input parameters
   const string treename_;
   const bool multiM_;
   const edm::EDGetTokenT<std::vector<reco::Track>> input_recoTrackMuon_token_;
   const edm::EDGetTokenT<std::vector<reco::Track>> input_recoTrack_token_;
   const edm::EDGetTokenT<std::vector<reco::Vertex>> input_recoVertex_token_;

   edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
   const reco::BeamSpot *beamSpot_;

   edm::EDGetTokenT<edm::TriggerResults> triggerresults_;
   std::vector<std::string> FilterNames_;

   const std::vector<double> MMuMass_;
   const std::vector<double> MMuMassErr_;
   const std::vector<double> MTrkMass_;
   const std::vector<double> MTrkMassErr_;
   const double MuMass_;
   const double MuMassErr_;
   const std::vector<int> MuChargeFlat_;
   const std::vector<double> TrkMassFlat_;
   const std::vector<double> TrkMassErrFlat_;
   const std::vector<int> TrkChargeFlat_;
   const unsigned int NMmu_;
   const unsigned int NMtrk_;
   const std::vector<unsigned int> NPmu_;
   const std::vector<unsigned int> NPtrk_;
   const double XMass_;
   const double XMassWin_;
   const double MMassWin_;
   const int XCharge_;
   const std::vector<int> MChargeMu_;
   const std::vector<int> MChargetrk_;
   const double PIso_;
   const bool doIso_;
   const bool doTrigger_;
   const double vProbMin_;
   const unsigned long long maxLoop_;
   const double MMassMin_;
   const double XMassMin_;

   const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;

   const AnalyticalImpactPointExtrapolator *impactPointExtrapolator;
   std::vector<reco::Vertex> vertices_;

   // Define output parameters
   std::string file_name;
   edm::Service<TFileService> fs;

   TTree *Scout4BTree;
   TTree *NMC_X_Scout4BTree;
   TTree *NMC_Scout4BTree;

   // Define branches
   UInt_t run;
   ULong64_t event;
   UInt_t lumiblock;
   UInt_t trigger;

   UInt_t nXcand;
   UInt_t nXcandNMC;
   UInt_t nMmucand[4];
   UInt_t nMtrkcand[4];

   UInt_t nMmu;
   UInt_t nMtrk;
   UInt_t nDmu[4];
   UInt_t nDtrk[4];

   // X particle branches
   double xMass[5];
   double xMassError[5];
   double xPt[5];
   double xEta[5];
   double xPhi[5];
   double xChi2[5];
   int xCharge[5];

   double xMassNMC[5];
   double xMassErrorNMC[5];
   double xPtNMC[5];
   double xEtaNMC[5];
   double xPhiNMC[5];
   double xChi2NMC[5];
   int xChargeNMC[5];

   double xlxy[5];
   double xlxyErr[5];
   double xlz[5];
   double xlzErr[5];
   double xl[5];
   double xlErr[5];
   double xsigLxy[5];
   double xsigLz[5];
   double xsigL[5];
   double xalphaBS[5];
   double xalphaBSErr[5];
   double xalphaBSXY[5];
   double xalphaBSXYErr[5];

   double xlxyNMC[5];
   double xlxyErrNMC[5];
   double xlzNMC[5];
   double xlzErrNMC[5];
   double xlNMC[5];
   double xlErrNMC[5];
   double xsigLxyNMC[5];
   double xsigLzNMC[5];
   double xsigLNMC[5];
   double xalphaBSNMC[5];
   double xalphaBSErrNMC[5];
   double xalphaBSXYNMC[5];
   double xalphaBSXYErrNMC[5];

   // M particle branches
   double mMass[5][4];
   double mMassError[5][4];
   double mPt[5][4];
   double mEta[5][4];
   double mPhi[5][4];
   int mCharge[5][4];

   double mMassNMC[5][4];
   double mPtNMC[5][4];
   double mEtaNMC[5][4];
   double mPhiNMC[5][4];
   int mChargeNMC[5][4];

   double NMC_mMassMu[5][4];
   double NMC_mMassErrMu[5][4];
   double NMC_mPtMu[5][4];
   double NMC_mEtaMu[5][4];
   double NMC_mPhiMu[5][4];
   double NMC_mChi2Mu[5][4];
   int NMC_mChargeMu[5][4];

   double NMC_mlxyMu[5][4];
   double NMC_mlxyErrMu[5][4];
   double NMC_mlzMu[5][4];
   double NMC_mlzErrMu[5][4];
   double NMC_mlMu[5][4];
   double NMC_mlErrMu[5][4];
   double NMC_msigLxyMu[5][4];
   double NMC_msigLzMu[5][4];
   double NMC_msigLMu[5][4];
   double NMC_malphaBSMu[5][4];
   double NMC_malphaBSErrMu[5][4];
   double NMC_malphaBSXYMu[5][4];
   double NMC_malphaBSXYErrMu[5][4];

   double NMC_mMassTrk[5][4];
   double NMC_mMassErrTrk[5][4];
   double NMC_mPtTrk[5][4];
   double NMC_mEtaTrk[5][4];
   double NMC_mPhiTrk[5][4];
   double NMC_mChi2Trk[5][4];
   int NMC_mChargeTrk[5][4];

   double NMC_mlxyTrk[5][4];
   double NMC_mlxyErrTrk[5][4];
   double NMC_mlzTrk[5][4];
   double NMC_mlzErrTrk[5][4];
   double NMC_mlTrk[5][4];
   double NMC_mlErrTrk[5][4];
   double NMC_msigLxyTrk[5][4];
   double NMC_msigLzTrk[5][4];
   double NMC_msigLTrk[5][4];
   double NMC_malphaBSTrk[5][4];
   double NMC_malphaBSErrTrk[5][4];
   double NMC_malphaBSXYTrk[5][4];
   double NMC_malphaBSXYErrTrk[5][4];

   // Daughter particle branches
   double dMass[5][4][8];
   double dPt[5][4][8];
   double dEta[5][4][8];
   double dPhi[5][4][8];
   int dCharge[5][4][8];

   double NMC_dMassMu[5][4][8];
   double NMC_dPtMu[5][4][8];
   double NMC_dEtaMu[5][4][8];
   double NMC_dPhiMu[5][4][8];
   int NMC_dChargeMu[5][4][8];

   double NMC_dMassTrk[5][4][8];
   double NMC_dPtTrk[5][4][8];
   double NMC_dEtaTrk[5][4][8];
   double NMC_dPhiTrk[5][4][8];
   int NMC_dChargeTrk[5][4][8];
};

//
// constructors and destructor
//
Scout4BRecoSecondaryVertexAnalyzer::Scout4BRecoSecondaryVertexAnalyzer(const edm::ParameterSet &iConfig)
    : treename_(iConfig.getUntrackedParameter<std::string>("treename", "Scout4B")),
      multiM_(iConfig.getUntrackedParameter<bool>("multiM", false)),
      input_recoTrackMuon_token_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("recoTrackMuon"))),
      input_recoTrack_token_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("recoTrack"))),
      input_recoVertex_token_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("recoVertex"))),
      beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
      beamSpot_(nullptr),
      triggerresults_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
      FilterNames_(iConfig.getParameter<std::vector<std::string>>("FilterNames")),
      MMuMass_(iConfig.getParameter<std::vector<double>>("MMuMass")),
      MMuMassErr_(iConfig.getParameter<std::vector<double>>("MMuMassErr")),
      MTrkMass_(iConfig.getParameter<std::vector<double>>("MTrkMass")),
      MTrkMassErr_(iConfig.getParameter<std::vector<double>>("MTrkMassErr")),
      MuMass_(iConfig.getUntrackedParameter<double>("MuMass", 0.1056583755)),
      MuMassErr_(iConfig.getUntrackedParameter<double>("MuMassErr", 1e-7)),
      MuChargeFlat_(iConfig.getParameter<std::vector<int>>("MuChargeFlat")),
      TrkMassFlat_(iConfig.getParameter<std::vector<double>>("TrkMassFlat")),
      TrkMassErrFlat_(iConfig.getParameter<std::vector<double>>("TrkMassErrFlat")),
      TrkChargeFlat_(iConfig.getParameter<std::vector<int>>("TrkChargeFlat")),
      NMmu_(iConfig.getParameter<unsigned int>("NMmu")),
      NMtrk_(iConfig.getParameter<unsigned int>("NMtrk")),
      NPmu_(iConfig.getParameter<std::vector<unsigned int>>("NPmu")),
      NPtrk_(iConfig.getParameter<std::vector<unsigned int>>("NPtrk")),
      XMass_(iConfig.getUntrackedParameter<double>("XMass", 0.0)),
      XMassWin_(iConfig.getUntrackedParameter<double>("XMassWin", 1000.0)),
      MMassWin_(iConfig.getUntrackedParameter<double>("MMassWin", 0.5)),
      XCharge_(iConfig.getParameter<int>("XCharge")),
      MChargeMu_(iConfig.getParameter<std::vector<int>>("MChargeMu")),
      MChargetrk_(iConfig.getParameter<std::vector<int>>("MChargetrk")),
      PIso_(iConfig.getUntrackedParameter<double>("PIso", 0.02)),
      doIso_(iConfig.getUntrackedParameter<bool>("doIso", true)),
      doTrigger_(iConfig.getUntrackedParameter<bool>("doTrigger", false)),
      vProbMin_(iConfig.getUntrackedParameter<double>("vProbMin", 1e-3)),
      maxLoop_(iConfig.getUntrackedParameter<unsigned long long>("maxLoop", 1000000000)),
      MMassMin_(iConfig.getUntrackedParameter<double>("MMassMin", 1e-2)),
      XMassMin_(iConfig.getUntrackedParameter<double>("XMassMin", 1e-2)),
      magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>())
{
}

Scout4BRecoSecondaryVertexAnalyzer::~Scout4BRecoSecondaryVertexAnalyzer() = default;

//
// Begin and end job
//
void Scout4BRecoSecondaryVertexAnalyzer::beginJob()
{
   // Initialize the branches
   clearVars();

   Scout4BTree = fs->make<TTree>(Form("%sTree", treename_.c_str()), Form("Tree of %s", treename_.c_str()));
   NMC_X_Scout4BTree = fs->make<TTree>(Form("NMC_X_%sTree", treename_.c_str()), Form("Tree of %s with no mass constraint", treename_.c_str()));
   NMC_Scout4BTree = fs->make<TTree>(Form("NMC_%sTree", treename_.c_str()), Form("Tree of %s with no mass constraint", treename_.c_str()));

   // Define branches
   Scout4BTree->Branch("run", &run, "run/i");
   Scout4BTree->Branch("event", &event, "event/l");
   Scout4BTree->Branch("lumiblock", &lumiblock, "lumiblock/i");
   Scout4BTree->Branch("trigger", &trigger, "trigger/i");
   Scout4BTree->Branch("nXcand", &nXcand, "nXcand/i");
   Scout4BTree->Branch("nMmu", &nMmu, "nMmu/i");
   Scout4BTree->Branch("nMtrk", &nMtrk, "nMtrk/i");
   Scout4BTree->Branch("nDmu", nDmu, "nDmu[4]/i");
   Scout4BTree->Branch("nDtrk", nDtrk, "nDtrk[4]/i");

   Scout4BTree->Branch("X_Mass", xMass, "X_Mass[5]/D");
   Scout4BTree->Branch("X_MassError", xMassError, "X_MassError[5]/D");
   Scout4BTree->Branch("X_Pt", xPt, "X_Pt[5]/D");
   Scout4BTree->Branch("X_Eta", xEta, "X_Eta[5]/D");
   Scout4BTree->Branch("X_Phi", xPhi, "X_Phi[5]/D");
   Scout4BTree->Branch("X_Chi2", xChi2, "X_Chi2[5]/D");
   Scout4BTree->Branch("X_Charge", xCharge, "X_Charge[5]/I");

   Scout4BTree->Branch("X_lxy", xlxy, "X_lxy[5]/D");
   Scout4BTree->Branch("X_lxyErr", xlxyErr, "X_lxyErr[5]/D");
   Scout4BTree->Branch("X_lz", xlz, "X_lz[5]/D");
   Scout4BTree->Branch("X_lzErr", xlzErr, "X_lzErr[5]/D");
   Scout4BTree->Branch("X_l", xl, "X_l[5]/D");
   Scout4BTree->Branch("X_lErr", xlErr, "X_lErr[5]/D");
   Scout4BTree->Branch("X_sigLxy", xsigLxy, "X_sigLxy[5]/D");
   Scout4BTree->Branch("X_sigLz", xsigLz, "X_sigLz[5]/D");
   Scout4BTree->Branch("X_sigL", xsigL, "X_sigL[5]/D");
   Scout4BTree->Branch("X_alphaBS", xalphaBS, "X_alphaBS[5]/D");
   Scout4BTree->Branch("X_alphaBSErr", xalphaBSErr, "X_alphaBSErr[5]/D");
   Scout4BTree->Branch("X_alphaBSXY", xalphaBSXY, "X_alphaBSXY[5]/D");
   Scout4BTree->Branch("X_alphaBSXYErr", xalphaBSXYErr, "X_alphaBSXYErr[5]/D");

   Scout4BTree->Branch("M_Mass", mMass, "M_Mass[5][4]/D");
   Scout4BTree->Branch("M_MassError", mMassError, "M_MassError[5][4]/D");
   Scout4BTree->Branch("M_Pt", mPt, "M_Pt[5][4]/D");
   Scout4BTree->Branch("M_Eta", mEta, "M_Eta[5][4]/D");
   Scout4BTree->Branch("M_Phi", mPhi, "M_Phi[5][4]/D");
   Scout4BTree->Branch("M_Charge", mCharge, "M_Charge[5][4]/I");

   Scout4BTree->Branch("D_Mass", dMass, "D_Mass[5][4][8]/D");
   Scout4BTree->Branch("D_Pt", dPt, "D_Pt[5][4][8]/D");
   Scout4BTree->Branch("D_Eta", dEta, "D_Eta[5][4][8]/D");
   Scout4BTree->Branch("D_Phi", dPhi, "D_Phi[5][4][8]/D");
   Scout4BTree->Branch("D_Charge", dCharge, "D_Charge[5][4][8]/I");

   NMC_Scout4BTree->Branch("run", &run, "run/i");
   NMC_Scout4BTree->Branch("event", &event, "event/l");
   NMC_Scout4BTree->Branch("lumiblock", &lumiblock, "lumiblock/i");
   NMC_Scout4BTree->Branch("trigger", &trigger, "trigger/i");
   NMC_Scout4BTree->Branch("nMmucand", &nMmucand, "nMmucand[4]/i");
   NMC_Scout4BTree->Branch("nMtrkcand", &nMtrkcand, "nMtrkcand[4]/i");
   NMC_Scout4BTree->Branch("nMmu", &nMmu, "nMmu/i");
   NMC_Scout4BTree->Branch("nMtrk", &nMtrk, "nMtrk/i");
   NMC_Scout4BTree->Branch("nDmu", nDmu, "nDmu[4]/i");
   NMC_Scout4BTree->Branch("nDtrk", nDtrk, "nDtrk[4]/i");

   NMC_Scout4BTree->Branch("NMC_M_MassMu", NMC_mMassMu, "NMC_M_MassMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_MassErrMu", NMC_mMassErrMu, "NMC_M_MassErrMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_PtMu", NMC_mPtMu, "NMC_M_PtMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_EtaMu", NMC_mEtaMu, "NMC_M_EtaMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_PhiMu", NMC_mPhiMu, "NMC_M_PhiMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_ChargeMu", NMC_mChargeMu, "NMC_M_ChargeMu[5][4]/I");
   NMC_Scout4BTree->Branch("NMC_M_Chi2Mu", NMC_mChi2Mu, "NMC_M_Chi2Mu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_MassTrk", NMC_mMassTrk, "NMC_M_MassTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_MassErrTrk", NMC_mMassErrTrk, "NMC_M_MassErrTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_PtTrk", NMC_mPtTrk, "NMC_M_PtTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_EtaTrk", NMC_mEtaTrk, "NMC_M_EtaTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_PhiTrk", NMC_mPhiTrk, "NMC_M_PhiTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_ChargeTrk", NMC_mChargeTrk, "NMC_M_ChargeTrk[5][4]/I");
   NMC_Scout4BTree->Branch("NMC_M_Chi2Trk", NMC_mChi2Trk, "NMC_M_Chi2Trk[5][4]/D");

   NMC_Scout4BTree->Branch("NMC_M_lxyMu", NMC_mlxyMu, "NMC_M_lxyMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_lxyErrMu", NMC_mlxyErrMu, "NMC_M_lxyErrMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_lzMu", NMC_mlzMu, "NMC_M_lzMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_lzErrMu", NMC_mlzErrMu, "NMC_M_lzErrMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_lMu", NMC_mlMu, "NMC_M_lMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_lErrMu", NMC_mlErrMu, "NMC_M_lErrMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_sigLxyMu", NMC_msigLxyMu, "NMC_M_sigLxyMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_sigLzMu", NMC_msigLzMu, "NMC_M_sigLzMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_sigLMu", NMC_msigLMu, "NMC_M_sigLMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_alphaBSMu", NMC_malphaBSMu, "NMC_M_alphaBSMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_alphaBSErrMu", NMC_malphaBSErrMu, "NMC_M_alphaBSErrMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_alphaBSXYMu", NMC_malphaBSXYMu, "NMC_M_alphaBSXYMu[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_alphaBSXYErrMu", NMC_malphaBSXYErrMu, "NMC_M_alphaBSXYErrMu[5][4]/D");

   NMC_Scout4BTree->Branch("NMC_M_lxyTrk", NMC_mlxyTrk, "NMC_M_lxyTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_lxyErrTrk", NMC_mlxyErrTrk, "NMC_M_lxyErrTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_lzTrk", NMC_mlzTrk, "NMC_M_lzTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_lzErrTrk", NMC_mlzErrTrk, "NMC_M_lzErrTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_lTrk", NMC_mlTrk, "NMC_M_lTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_lErrTrk", NMC_mlErrTrk, "NMC_M_lErrTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_sigLxyTrk", NMC_msigLxyTrk, "NMC_M_sigLxyTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_sigLzTrk", NMC_msigLzTrk, "NMC_M_sigLzTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_sigLTrk", NMC_msigLTrk, "NMC_M_sigLTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_alphaBSTrk", NMC_malphaBSTrk, "NMC_M_alphaBSTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_alphaBSErrTrk", NMC_malphaBSErrTrk, "NMC_M_alphaBSErrTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_alphaBSXYTrk", NMC_malphaBSXYTrk, "NMC_M_alphaBSXYTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_alphaBSXYErrTrk", NMC_malphaBSXYErrTrk, "NMC_M_alphaBSXYErrTrk[5][4]/D");

   NMC_Scout4BTree->Branch("NMC_D_MassMu", NMC_dMassMu, "NMC_D_MassMu[5][4][8]/D");
   NMC_Scout4BTree->Branch("NMC_D_PtMu", NMC_dPtMu, "NMC_D_PtMu[5][4][8]/D");
   NMC_Scout4BTree->Branch("NMC_D_EtaMu", NMC_dEtaMu, "NMC_D_EtaMu[5][4][8]/D");
   NMC_Scout4BTree->Branch("NMC_D_PhiMu", NMC_dPhiMu, "NMC_D_PhiMu[5][4][8]/D");
   NMC_Scout4BTree->Branch("NMC_D_ChargeMu", NMC_dChargeMu, "NMC_D_ChargeMu[5][4][8]/I");
   NMC_Scout4BTree->Branch("NMC_D_MassTrk", NMC_dMassTrk, "NMC_D_MassTrk[5][4][8]/D");
   NMC_Scout4BTree->Branch("NMC_D_PtTrk", NMC_dPtTrk, "NMC_D_PtTrk[5][4][8]/D");
   NMC_Scout4BTree->Branch("NMC_D_EtaTrk", NMC_dEtaTrk, "NMC_D_EtaTrk[5][4][8]/D");
   NMC_Scout4BTree->Branch("NMC_D_PhiTrk", NMC_dPhiTrk, "NMC_D_PhiTrk[5][4][8]/D");
   NMC_Scout4BTree->Branch("NMC_D_ChargeTrk", NMC_dChargeTrk, "NMC_D_ChargeTrk[5][4][8]/I");

   NMC_X_Scout4BTree->Branch("run", &run, "run/i");
   NMC_X_Scout4BTree->Branch("event", &event, "event/l");
   NMC_X_Scout4BTree->Branch("lumiblock", &lumiblock, "lumiblock/i");
   NMC_X_Scout4BTree->Branch("trigger", &trigger, "trigger/i");
   NMC_X_Scout4BTree->Branch("nXcandNMC", &nXcandNMC, "nXcandNMC/i");

   NMC_X_Scout4BTree->Branch("X_Mass", xMassNMC, "X_Mass[5]/D");
   NMC_X_Scout4BTree->Branch("X_MassError", xMassErrorNMC, "X_MassError[5]/D");
   NMC_X_Scout4BTree->Branch("X_Pt", xPtNMC, "X_Pt[5]/D");
   NMC_X_Scout4BTree->Branch("X_Eta", xEtaNMC, "X_Eta[5]/D");
   NMC_X_Scout4BTree->Branch("X_Phi", xPhiNMC, "X_Phi[5]/D");
   NMC_X_Scout4BTree->Branch("X_Chi2", xChi2NMC, "X_Chi2[5]/D");
   NMC_X_Scout4BTree->Branch("X_Charge", xChargeNMC, "X_Charge[5]/I");

   NMC_X_Scout4BTree->Branch("X_lxy", xlxyNMC, "X_lxy[5]/D");
   NMC_X_Scout4BTree->Branch("X_lxyErr", xlxyErrNMC, "X_lxyErr[5]/D");
   NMC_X_Scout4BTree->Branch("X_lz", xlzNMC, "X_lz[5]/D");
   NMC_X_Scout4BTree->Branch("X_lzErr", xlzErrNMC, "X_lzErr[5]/D");
   NMC_X_Scout4BTree->Branch("X_l", xlNMC, "X_l[5]/D");
   NMC_X_Scout4BTree->Branch("X_lErr", xlErrNMC, "X_lErr[5]/D");
   NMC_X_Scout4BTree->Branch("X_sigLxy", xsigLxyNMC, "X_sigLxy[5]/D");
   NMC_X_Scout4BTree->Branch("X_sigLz", xsigLzNMC, "X_sigLz[5]/D");
   NMC_X_Scout4BTree->Branch("X_sigL", xsigLNMC, "X_sigL[5]/D");
   NMC_X_Scout4BTree->Branch("X_alphaBS", xalphaBSNMC, "X_alphaBS[5]/D");
   NMC_X_Scout4BTree->Branch("X_alphaBSErr", xalphaBSErrNMC, "X_alphaBSErr[5]/D");
   NMC_X_Scout4BTree->Branch("X_alphaBSXY", xalphaBSXYNMC, "X_alphaBSXY[5]/D");
   NMC_X_Scout4BTree->Branch("X_alphaBSXYErr", xalphaBSXYErrNMC, "X_alphaBSXYErr[5]/D");

   NMC_X_Scout4BTree->Branch("M_Mass", mMassNMC, "M_Mass[5][4]/D");
   NMC_X_Scout4BTree->Branch("M_Pt", mPtNMC, "M_Pt[5][4]/D");
   NMC_X_Scout4BTree->Branch("M_Eta", mEtaNMC, "M_Eta[5][4]/D");
   NMC_X_Scout4BTree->Branch("M_Phi", mPhiNMC, "M_Phi[5][4]/D");
   NMC_X_Scout4BTree->Branch("M_Charge", mChargeNMC, "M_Charge[5][4]/I");
}

void Scout4BRecoSecondaryVertexAnalyzer::endJob()
{
   clearVars();
}

//
// ------------ method called to analyze the data  ------------
//
void Scout4BRecoSecondaryVertexAnalyzer::analyze(edm::Event const &iEvent, edm::EventSetup const &iSetup)
{
   // Init the branches
   clearVars();
   const MagneticField* bFieldHandle = &iSetup.getData(magneticFieldToken_);

   AnalyticalImpactPointExtrapolator extrapolator(bFieldHandle);
   impactPointExtrapolator = &extrapolator;

   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByToken(beamSpotToken_, beamSpotHandle);

   if (!beamSpotHandle.isValid())
   {
      cout << "Invalid beam spot" << endl;
      return;
   }

   beamSpot_ = beamSpotHandle.product();

   // Collect inputs
   edm::Handle<std::vector<reco::Track>> recoTrackMuon;
   iEvent.getByToken(input_recoTrackMuon_token_, recoTrackMuon);
   edm::Handle<std::vector<reco::Track>> recoTrack;
   iEvent.getByToken(input_recoTrack_token_, recoTrack);
   edm::Handle<std::vector<reco::Vertex>> recoVertex;
   iEvent.getByToken(input_recoVertex_token_, recoVertex);

   // get vertices_
   vertices_ = *recoVertex;

   run = iEvent.id().run();
   event = iEvent.id().event();
   lumiblock = iEvent.id().luminosityBlock();
   trigger = getTriggerBits(iEvent);

   // Input filter
   if (NMmu_ + NMtrk_ > 4)
   {
      cout << "Too many particles or mass constraints" << endl;
      return;
   }
   for (unsigned int i = 0; i < NMmu_; i++)
   {
      if (NPmu_[i] > 8)
      {
         cout << "Too many muons in a mass constraint" << endl;
         return;
      }
   }
   for (unsigned int i = 0; i < NMtrk_; i++)
   {
      if (NPtrk_[i] > 8)
      {
         cout << "Too many tracks in a mass constraint" << endl;
         return;
      }
   }
   nMmu = NMmu_;
   nMtrk = NMtrk_;
   for (unsigned int i = 0; i < NMmu_; i++)
   {
      nDmu[i] = NPmu_[i];
   }
   for (unsigned int i = 0; i < NMtrk_; i++)
   {
      nDtrk[i] = NPtrk_[i];
   }

   // Event filter
   unsigned int nMuInNeed = 0;
   unsigned int nTrkInNeed = 0;
   for (unsigned int i = 0; i < NMmu_; i++)
   {
      nMuInNeed += NPmu_[i];
   }
   for (unsigned int i = 0; i < NMtrk_; i++)
   {
      nTrkInNeed += NPtrk_[i];
   }
   if (doTrigger_ && !trigger)
      return;
   if (!recoTrackMuon.isValid() || !recoTrack.isValid())
   {
      // cout << "Invalid input collections" << endl;
      return;
   }
   if (recoTrackMuon->size() < nMuInNeed || recoTrack->size() < nTrkInNeed)
   {
      // cout << "Input collections too small: " << recoTrackMuon->size() << " " << recoTrack->size() << endl;
      return;
   }
   // cout << "Input collections: " << recoTrackMuon->size() << " " << recoTrack->size() << endl;

   // Track filter
   std::vector<reco::Track> goodTracks;
   std::vector<reco::Track> goodTrackMuons;
   for (const auto &trk : *recoTrack)
   {
      if (trk.charge() == 0 || trk.pt() < 0.5)
         continue;
      bool overlap = false;
      if (doIso_)
      {
         for (const auto &mu : *recoTrackMuon)
         {
            if (deltaR(trk, mu) < PIso_)
            {
               overlap = true;
               break;
            }
         }
      }
      if (overlap)
         continue;
      goodTracks.push_back(trk);
   }
   for (const auto &mu : *recoTrackMuon)
   {
      if (mu.charge() == 0 || mu.pt() < 0.5)
         continue;
      goodTrackMuons.push_back(mu);
   }

   if (goodTracks.size() < nTrkInNeed || goodTrackMuons.size() < nMuInNeed)
   {
      // cout << "Good tracks or muons too few: " << goodTracks.size() << " " << goodTrackMuons.size() << endl;
      return;
   }

   // Unflat TrkMassFlat, TrkMassErrFlat and MuChargeFlat
   std::vector<std::vector<double>> TrkMass_v_;
   std::vector<std::vector<double>> TrkMassErr_v_;
   std::vector<std::vector<int>> TrkCharge_v_;
   unsigned int offset = 0;
   for (unsigned int i = 0; i < NMtrk_; i++)
   {
      std::vector<double> TrkMass_v;
      std::vector<double> TrkMassErr_v;
      std::vector<int> TrkCharge_v;

      for (unsigned int j = 0; j < NPtrk_[i]; j++)
      {
         TrkMass_v.push_back(TrkMassFlat_[offset + j]);
         TrkMassErr_v.push_back(TrkMassErrFlat_[offset + j]);
         TrkCharge_v.push_back(TrkChargeFlat_[offset + j]);
      }
      TrkMass_v_.push_back(TrkMass_v);
      TrkMassErr_v_.push_back(TrkMassErr_v);
      TrkCharge_v_.push_back(TrkCharge_v);
      offset += NPtrk_[i];
   }

   std::vector<std::vector<int>> MuCharge_v_;
   offset = 0;
   for (unsigned int i = 0; i < NMmu_; i++)
   {
      std::vector<int> MuCharge;
      for (unsigned int j = 0; j < NPmu_[i]; j++)
      {
         MuCharge.push_back(MuChargeFlat_[j + offset]);
      }
      MuCharge_v_.push_back(MuCharge);
      offset += NPmu_[i];
   }

   // Step 1: Fit MMuMass M particles from muon trks
   vector<vector<KinematicFitResult>> mCandidatesMu;
   vector<vector<KinematicFitResult>> NMC_mCandidatesMu;
   for (unsigned int i = 0; i < NMmu_; i++)
   {
      // 构造用于拟合的 particle 质量及误差向量，长度为 NPmu_[i]，取输入的 MuMass 参数（保持不变）
      vector<double> PMass_mu(NPmu_[i], MuMass_);
      vector<double> PMassErr_mu(NPmu_[i], MuMassErr_);

      std::pair<std::vector<KinematicFitResult>, std::vector<KinematicFitResult>> res = performVertexFit(goodTrackMuons,
                                                                                                         NPmu_[i],
                                                                                                         MChargeMu_[i],
                                                                                                         MMuMass_[i],
                                                                                                         MMuMassErr_[i],
                                                                                                         MMassWin_ * 4.0,
                                                                                                         vProbMin_,
                                                                                                         PMass_mu,
                                                                                                         PMassErr_mu,
                                                                                                         MuCharge_v_[i],
                                                                                                         MMassMin_,
                                                                                                         doIso_,
                                                                                                         PIso_,
                                                                                                         *bFieldHandle,
                                                                                                         maxLoop_);
      if (res.second.empty())
         return;
      mCandidatesMu.push_back(res.second);
      NMC_mCandidatesMu.push_back(res.first);
   }

   if (mCandidatesMu.size() < NMmu_ || NMC_mCandidatesMu.size() < NMmu_)
      return;

   // Step 2: Fit M particles from trks
   vector<vector<KinematicFitResult>> mCandidatesTrk;
   vector<vector<KinematicFitResult>> NMC_mCandidatesTrk;
   for (unsigned int i = 0; i < NMtrk_; i++)
   {

      std::pair<std::vector<KinematicFitResult>, std::vector<KinematicFitResult>> res = performVertexFit(goodTracks,
                                                                                                         NPtrk_[i],
                                                                                                         MChargetrk_[i],
                                                                                                         MTrkMass_[i],
                                                                                                         MTrkMassErr_[i],
                                                                                                         MMassWin_,
                                                                                                         vProbMin_,
                                                                                                         TrkMass_v_[i],
                                                                                                         TrkMassErr_v_[i],
                                                                                                         TrkCharge_v_[i],
                                                                                                         MMassMin_,
                                                                                                         doIso_,
                                                                                                         PIso_,
                                                                                                         *bFieldHandle,
                                                                                                         maxLoop_);
      if (res.second.empty())
         return;
      mCandidatesTrk.push_back(res.second);
      NMC_mCandidatesTrk.push_back(res.first);
   }

   if (mCandidatesTrk.size() < NMtrk_ || NMC_mCandidatesTrk.size() < NMtrk_)
      return;

   // Sort the M candidates by M's pT
   for (auto &vec : NMC_mCandidatesMu)
   {
      std::sort(vec.begin(), vec.end(), [](const KinematicFitResult &a, const KinematicFitResult &b)
                { return a.p4().Pt() > b.p4().Pt(); });
   }
   for (auto &vec : NMC_mCandidatesTrk)
   {
      std::sort(vec.begin(), vec.end(), [](const KinematicFitResult &a, const KinematicFitResult &b)
                { return a.p4().Pt() > b.p4().Pt(); });
   }

   // Fill NMC M muon candidates
   for (unsigned int j = 0; j < NMmu_; j++)
   {
      nMmucand[j] = std::min(int(NMC_mCandidatesMu[j].size()), 5);
      for (unsigned int i = 0; i < nMmucand[j]; i++)
      {
         NMC_mMassMu[i][j] = NMC_mCandidatesMu[j][i].p4().M();
         NMC_mMassErrMu[i][j] = NMC_mCandidatesMu[j][i].massErr();
         NMC_mPtMu[i][j] = NMC_mCandidatesMu[j][i].p4().Pt();
         NMC_mEtaMu[i][j] = NMC_mCandidatesMu[j][i].p4().Eta();
         NMC_mPhiMu[i][j] = NMC_mCandidatesMu[j][i].p4().Phi();
         NMC_mChi2Mu[i][j] = NMC_mCandidatesMu[j][i].chi2() / (double)NMC_mCandidatesMu[j][i].ndof();
         NMC_mChargeMu[i][j] = NMC_mCandidatesMu[j][i].charge_;

         NMC_mlxyMu[i][j] = NMC_mCandidatesMu[j][i].lxy();
         NMC_mlxyErrMu[i][j] = NMC_mCandidatesMu[j][i].lxyErr();
         NMC_mlzMu[i][j] = NMC_mCandidatesMu[j][i].lz();
         NMC_mlzErrMu[i][j] = NMC_mCandidatesMu[j][i].lzErr();
         NMC_mlMu[i][j] = NMC_mCandidatesMu[j][i].l();
         NMC_mlErrMu[i][j] = NMC_mCandidatesMu[j][i].lErr();
         NMC_msigLxyMu[i][j] = NMC_mCandidatesMu[j][i].sigLxy();
         NMC_msigLzMu[i][j] = NMC_mCandidatesMu[j][i].sigLz();
         NMC_msigLMu[i][j] = NMC_mCandidatesMu[j][i].sigL();
         NMC_malphaBSMu[i][j] = NMC_mCandidatesMu[j][i].alphaBS();
         NMC_malphaBSErrMu[i][j] = NMC_mCandidatesMu[j][i].alphaBSErr();
         NMC_malphaBSXYMu[i][j] = NMC_mCandidatesMu[j][i].alphaBSXY();
         NMC_malphaBSXYErrMu[i][j] = NMC_mCandidatesMu[j][i].alphaBSXYErr();

         // Access daughters
         for (unsigned int k = 0; k < NMC_mCandidatesMu[j][i].number_of_daughters(); k++)
         {
            NMC_dMassMu[i][j][k] = NMC_mCandidatesMu[j][i].dau_p4(k).M();
            NMC_dPtMu[i][j][k] = NMC_mCandidatesMu[j][i].dau_p4(k).Pt();
            NMC_dEtaMu[i][j][k] = NMC_mCandidatesMu[j][i].dau_p4(k).Eta();
            NMC_dPhiMu[i][j][k] = NMC_mCandidatesMu[j][i].dau_p4(k).Phi();
            NMC_dChargeMu[i][j][k] = NMC_mCandidatesMu[j][i].dau_charge_[k];
         }
      }
   }

   // Fill NMC M track candidates
   for (unsigned int j = 0; j < NMtrk_; j++)
   {
      nMtrkcand[j] = std::min(int(NMC_mCandidatesTrk[j].size()), 5);
      for (unsigned int i = 0; i < nMtrkcand[j]; i++)
      {
         NMC_mMassTrk[i][j] = NMC_mCandidatesTrk[j][i].p4().M();
         NMC_mMassErrTrk[i][j] = NMC_mCandidatesTrk[j][i].massErr();
         NMC_mPtTrk[i][j] = NMC_mCandidatesTrk[j][i].p4().Pt();
         NMC_mEtaTrk[i][j] = NMC_mCandidatesTrk[j][i].p4().Eta();
         NMC_mPhiTrk[i][j] = NMC_mCandidatesTrk[j][i].p4().Phi();
         NMC_mChi2Trk[i][j] = NMC_mCandidatesTrk[j][i].chi2() / (double)NMC_mCandidatesTrk[j][i].ndof();
         NMC_mChargeTrk[i][j] = NMC_mCandidatesTrk[j][i].charge_;

         NMC_mlxyTrk[i][j] = NMC_mCandidatesTrk[j][i].lxy();
         NMC_mlxyErrTrk[i][j] = NMC_mCandidatesTrk[j][i].lxyErr();
         NMC_mlzTrk[i][j] = NMC_mCandidatesTrk[j][i].lz();
         NMC_mlzErrTrk[i][j] = NMC_mCandidatesTrk[j][i].lzErr();
         NMC_mlTrk[i][j] = NMC_mCandidatesTrk[j][i].l();
         NMC_mlErrTrk[i][j] = NMC_mCandidatesTrk[j][i].lErr();
         NMC_msigLxyTrk[i][j] = NMC_mCandidatesTrk[j][i].sigLxy();
         NMC_msigLzTrk[i][j] = NMC_mCandidatesTrk[j][i].sigLz();
         NMC_msigLTrk[i][j] = NMC_mCandidatesTrk[j][i].sigL();
         NMC_malphaBSTrk[i][j] = NMC_mCandidatesTrk[j][i].alphaBS();
         NMC_malphaBSErrTrk[i][j] = NMC_mCandidatesTrk[j][i].alphaBSErr();
         NMC_malphaBSXYTrk[i][j] = NMC_mCandidatesTrk[j][i].alphaBSXY();
         NMC_malphaBSXYErrTrk[i][j] = NMC_mCandidatesTrk[j][i].alphaBSXYErr();

         // Access daughters
         for (unsigned int k = 0; k < NMC_mCandidatesTrk[j][i].number_of_daughters(); k++)
         {
            NMC_dMassTrk[i][j][k] = NMC_mCandidatesTrk[j][i].dau_p4(k).M();
            NMC_dPtTrk[i][j][k] = NMC_mCandidatesTrk[j][i].dau_p4(k).Pt();
            NMC_dEtaTrk[i][j][k] = NMC_mCandidatesTrk[j][i].dau_p4(k).Eta();
            NMC_dPhiTrk[i][j][k] = NMC_mCandidatesTrk[j][i].dau_p4(k).Phi();
            NMC_dChargeTrk[i][j][k] = NMC_mCandidatesTrk[j][i].dau_charge_[k];
         }
      }
   }

   NMC_Scout4BTree->Fill();

   // Step 3: multiM == false -> directly exit
   if (!multiM_)
   {
      clearVars();
      return;
   }

   // Step 4: Fit X particles from both M particles
   vector<vector<KinematicFitResult>> allMCandidates;
   for (auto &vec : mCandidatesMu)
   {
      allMCandidates.push_back(vec);
   }
   for (auto &vec : mCandidatesTrk)
   {
      allMCandidates.push_back(vec);
   }
   if (allMCandidates.size() < NMmu_ + NMtrk_)
      return;
   unsigned int totalM = NMmu_ + NMtrk_;

   // Add goodTrackMuons to goodTracks end
   // and edit the index in NMC_mCandidatesMu to match the new goodTracks
   for (unsigned int i = 0; i < NMmu_; i++)
   {
      for (auto &v : NMC_mCandidatesMu[i])
      {
         for (unsigned int j = 0; j < v.trackIndices.size(); j++)
         {
            v.trackIndices[j] = v.trackIndices[j] + goodTracks.size();
         }
      }
   }
   for (const auto &mu : goodTrackMuons)
   {
      goodTracks.push_back(mu);
   }

   // 对于 X 拟合，构造每个 M 粒子的质量约束要求，取输入的 MTrkMass、MTrkMassErr、MMuMass、MMuMassErr 参数
   vector<double> mMassForX;
   vector<double> mMassErrForX;
   for (unsigned int i = 0; i < totalM; i++)
   {
      if (i < NMmu_)
      {
         mMassForX.push_back(MMuMass_[i]);
         mMassErrForX.push_back(MMuMassErr_[i]);
      }
      else
      {
         mMassForX.push_back(MTrkMass_[i - NMmu_]);
         mMassErrForX.push_back(MTrkMassErr_[i - NMmu_]);
      }
   }

   std::pair<std::vector<XFitResult>, std::vector<KinematicFitResult>> xCandidatesAll = performXVertexFit(goodTracks, totalM, allMCandidates, XCharge_, XMass_, XMassWin_, vProbMin_, mMassForX, mMassErrForX, XMassMin_, doIso_, PIso_, *bFieldHandle, maxLoop_);

   std::vector<KinematicFitResult> xCandidatesNMC = xCandidatesAll.second;

   // Fill NMC_X_Scout4BTree tree
   if (!xCandidatesNMC.empty())
   {
      // Sort the X candidates by X's pT
      std::sort(xCandidatesNMC.begin(), xCandidatesNMC.end(), [](const KinematicFitResult &a, const KinematicFitResult &b)
                { return a.p4().Pt() > b.p4().Pt(); });

      // Cut first 5 X candidates if exist
      if (xCandidatesNMC.size() > 5)
         xCandidatesNMC.resize(5);
      nXcandNMC = xCandidatesNMC.size();

      for (unsigned int i = 0; i < nXcandNMC; i++)
      {
         xMassNMC[i] = xCandidatesNMC[i].p4().M();
         xMassErrorNMC[i] = xCandidatesNMC[i].massErr();
         xPtNMC[i] = xCandidatesNMC[i].p4().Pt();
         xEtaNMC[i] = xCandidatesNMC[i].p4().Eta();
         xPhiNMC[i] = xCandidatesNMC[i].p4().Phi();
         xChi2NMC[i] = xCandidatesNMC[i].chi2() / (double)xCandidatesNMC[i].ndof();
         xChargeNMC[i] = xCandidatesNMC[i].charge_;

         xlxyNMC[i] = xCandidatesNMC[i].lxy();
         xlxyErrNMC[i] = xCandidatesNMC[i].lxyErr();
         xlzNMC[i] = xCandidatesNMC[i].lz();
         xlzErrNMC[i] = xCandidatesNMC[i].lzErr();
         xlNMC[i] = xCandidatesNMC[i].l();
         xlErrNMC[i] = xCandidatesNMC[i].lErr();
         xsigLxyNMC[i] = xCandidatesNMC[i].sigLxy();
         xsigLzNMC[i] = xCandidatesNMC[i].sigLz();
         xsigLNMC[i] = xCandidatesNMC[i].sigL();
         xalphaBSNMC[i] = xCandidatesNMC[i].alphaBS();
         xalphaBSErrNMC[i] = xCandidatesNMC[i].alphaBSErr();
         xalphaBSXYNMC[i] = xCandidatesNMC[i].alphaBSXY();
         xalphaBSXYErrNMC[i] = xCandidatesNMC[i].alphaBSXYErr();

         // Access mothers
         for (unsigned int j = 0; j < xCandidatesNMC[i].number_of_daughters(); j++)
         {
            mMassNMC[i][j] = xCandidatesNMC[i].dau_p4(j).M();
            mPtNMC[i][j] = xCandidatesNMC[i].dau_p4(j).Pt();
            mEtaNMC[i][j] = xCandidatesNMC[i].dau_p4(j).Eta();
            mPhiNMC[i][j] = xCandidatesNMC[i].dau_p4(j).Phi();
            mChargeNMC[i][j] = xCandidatesNMC[i].dau_charge_[j];
         }
      }

      NMC_X_Scout4BTree->Fill();
   }

   std::vector<XFitResult> xCandidates = xCandidatesAll.first;
   if (!xCandidates.empty())
   {
      // Step 4: Fill the tree

      // Sort the X candidates by X's pT
      std::sort(xCandidates.begin(), xCandidates.end(), [](const XFitResult &a, const XFitResult &b)
                { return a.xMotherP4.Pt() > b.xMotherP4.Pt(); });

      // Cut first 5 X candidates if exist
      if (xCandidates.size() > 5)
         xCandidates.resize(5);
      nXcand = xCandidates.size();

      // Fill the tree
      for (unsigned int i = 0; i < nXcand; i++)
      {
         xMass[i] = xCandidates[i].xMotherP4.M();
         xMassError[i] = xCandidates[i].xMassError;
         xPt[i] = xCandidates[i].xMotherP4.Pt();
         xEta[i] = xCandidates[i].xMotherP4.Eta();
         xPhi[i] = xCandidates[i].xMotherP4.Phi();
         xChi2[i] = xCandidates[i].normChi2;
         xCharge[i] = xCandidates[i].xMotherCharge;

         xlxy[i] = xCandidates[i].xlxy;
         xlxyErr[i] = xCandidates[i].xlxyErr;
         xlz[i] = xCandidates[i].xlz;
         xlzErr[i] = xCandidates[i].xlzErr;
         xl[i] = xCandidates[i].xl;
         xlErr[i] = xCandidates[i].xlErr;
         xsigLxy[i] = xCandidates[i].xsigLxy;
         xsigLz[i] = xCandidates[i].xsigLz;
         xsigL[i] = xCandidates[i].xsigL;
         xalphaBS[i] = xCandidates[i].xalphaBS;
         xalphaBSErr[i] = xCandidates[i].xalphaBSErr;
         xalphaBSXY[i] = xCandidates[i].xalphaBSXY;
         xalphaBSXYErr[i] = xCandidates[i].xalphaBSXYErr;

         for (unsigned int j = 0; j < xCandidates[i].mMother.size(); j++)
         {
            mMass[i][j] = xCandidates[i].mMother[j].first.M();
            mMassError[i][j] = xCandidates[i].mMassError[j];
            mPt[i][j] = xCandidates[i].mMother[j].first.Pt();
            mEta[i][j] = xCandidates[i].mMother[j].first.Eta();
            mPhi[i][j] = xCandidates[i].mMother[j].first.Phi();
            mCharge[i][j] = xCandidates[i].mMother[j].second;

            for (unsigned int k = 0; k < xCandidates[i].mDaughter[j].first.size(); k++)
            {
               dMass[i][j][k] = xCandidates[i].mDaughter[j].first[k].M();
               dPt[i][j][k] = xCandidates[i].mDaughter[j].first[k].Pt();
               dEta[i][j][k] = xCandidates[i].mDaughter[j].first[k].Eta();
               dPhi[i][j][k] = xCandidates[i].mDaughter[j].first[k].Phi();
               dCharge[i][j][k] = xCandidates[i].mDaughter[j].second[k];
            }
         }
      }

      Scout4BTree->Fill();
   }

   // Print xCandidates information for debug
   /* for (unsigned int i = 0; i < xCandidates.size(); i++)
   {
      cout << "X candidate " << i << " : " << endl;
      cout << "X mass: " << xCandidates[i].xMotherP4.M() << endl;
      cout << "X charge: " << xCandidates[i].xMotherCharge << endl;
      cout << "X mass error: " << xCandidates[i].xMassError << endl;
      cout << "X chi2: " << xCandidates[i].normChi2 << endl;
      for (unsigned int j = 0; j < xCandidates[i].mDaughter.size(); j++)
      {
         cout << "M candidate " << j << " : " << endl;
         cout << "M mass: " << xCandidates[i].mMother[j].first.M() << endl;
         cout << "M charge: " << xCandidates[i].mMother[j].second << endl;
         for (unsigned int k = 0; k < xCandidates[i].mDaughter[j].first.size(); k++)
         {
            cout << "Daughter " << k << " : " << endl;
            cout << "Pt: " << xCandidates[i].mDaughter[j].first[k].Pt() << endl;
            cout << "Eta: " << xCandidates[i].mDaughter[j].first[k].Eta() << endl;
            cout << "Phi: " << xCandidates[i].mDaughter[j].first[k].Phi() << endl;
            cout << "Charge: " << xCandidates[i].mDaughter[j].second[k] << endl;
         }
      }
   }*/

   clearVars();
}

//
// Funtions definition
//
UInt_t Scout4BRecoSecondaryVertexAnalyzer::getTriggerBits(const edm::Event &iEvent)
{
   UInt_t trigger = 0;
   edm::Handle<edm::TriggerResults> triggerresults;
   iEvent.getByToken(triggerresults_, triggerresults);
   if (triggerresults.isValid())
   {
      const edm::TriggerNames &TheTriggerNames = iEvent.triggerNames(*triggerresults);
      for (unsigned int i = 0; i < FilterNames_.size(); i++)
      {
         bool matched = false;
         for (unsigned int version = 1; (version < 99 && (!matched)); version++)
         {
            std::stringstream ss;
            ss << FilterNames_[i] << "_v" << version;
            unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
            if (bit < triggerresults->size() && triggerresults->accept(bit) && !triggerresults->error(bit))
               matched = true;
         }
         if (matched)
            trigger += (1 << i);
      }
   }
   // else std::cout << "MMrootupler::getTriggerBits: *** NO triggerResults found *** " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;
   return trigger;
}

float Scout4BRecoSecondaryVertexAnalyzer::distanceOfClosestApproach(const reco::Track *track1, const reco::Track *track2, const MagneticField &bFieldHandle)
{
   TwoTrackMinimumDistance md;
   const reco::TransientTrack tt1(*track1, &bFieldHandle);
   const reco::TransientTrack tt2(*track2, &bFieldHandle);
   if (not md.calculate(tt1.initialFreeState(), tt2.initialFreeState()))
      return -1.0;
   return md.distance();
}

Measurement1D Scout4BRecoSecondaryVertexAnalyzer::distanceOfClosestApproach(const reco::Track *track, const reco::Vertex &vertex, const MagneticField &bFieldHandle)
{
   VertexDistance3D distance3D;
   const reco::TransientTrack tt(*track, &bFieldHandle);
   assert(impactPointExtrapolator);
   auto tsos = impactPointExtrapolator->extrapolate(tt.initialFreeState(), GlobalPoint(Basic3DVector<float>(vertex.position())));
   if (not tsos.isValid())
      return Measurement1D(-1.0, -1.0);
   Measurement1D doca = distance3D.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), vertex);
   return doca;
}

KinematicFitResult Scout4BRecoSecondaryVertexAnalyzer::KinematicFitter(std::vector<reco::Track *> trks, std::vector<double> masses, const MagneticField &bFieldHandle)
{
   // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideKinematicVertexFit
   if (trks.size() != masses.size())
      throw cms::Exception("Error") << "number of tracks and number of masses should match";

   KinematicParticleFactoryFromTransientTrack factory;
   KinematicParticleVertexFitter fitter;

   std::vector<RefCountedKinematicParticle> particles;
   double chi = 0.;
   double ndf = 0.;
   for (unsigned int i = 0; i < trks.size(); ++i)
   {
      const reco::TransientTrack tt(*trks[i], &bFieldHandle);
      float mass = masses[i];
      float massErr = mass / 1000.0;
      particles.push_back(factory.particle(tt, mass, chi, ndf, massErr));
   }

   RefCountedKinematicTree vertexFitTree;
   KinematicFitResult result;
   result.tracks = trks;
   try
   {
      vertexFitTree = fitter.fit(particles);
   }
   catch (const std::exception &e)
   {
      return result;
   }
   result.set_tree(vertexFitTree);
   return result;
}

KinematicFitResult Scout4BRecoSecondaryVertexAnalyzer::MCKinematicFitter(std::vector<reco::Track *> trks, std::vector<double> masses, double MMass, double MMassErr, const MagneticField &bFieldHandle)
{
   // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideKinematicVertexFit
   if (trks.size() != masses.size())
      throw cms::Exception("Error") << "number of tracks and number of masses should match";

   KinematicParticleFactoryFromTransientTrack factory;
   KinematicParticleVertexFitter fitter;
   KinematicParticleFitter MCfitter;
   MassKinematicConstraint massConstraint(MMass, MMassErr);

   std::vector<RefCountedKinematicParticle> particles;
   double chi = 0.;
   double ndf = 0.;
   for (unsigned int i = 0; i < trks.size(); ++i)
   {
      const reco::TransientTrack tt(*trks[i], &bFieldHandle);
      float mass = masses[i];
      float massErr = mass / 1000.0;
      particles.push_back(factory.particle(tt, mass, chi, ndf, massErr));
   }

   RefCountedKinematicTree vertexFitTree;
   RefCountedKinematicTree constrainedTree;
   KinematicFitResult result;
   result.tracks = trks;
   try
   {
      vertexFitTree = fitter.fit(particles);
   }
   catch (const std::exception &e)
   {
      return result;
   }

   if (vertexFitTree->isValid())
   {
      try
      {
         constrainedTree = MCfitter.fit(&massConstraint, vertexFitTree);
      }
      catch (const std::exception &e)
      {
         return result;
      }
   }
   else
   {
      return result;
   }

   result.set_tree(constrainedTree);
   return result;
}

Scout4B::Displacements Scout4BRecoSecondaryVertexAnalyzer::compute3dDisplacement(const KinematicFitResult &fit, bool closestIn3D)
{
   // WARNING: all variables need to be filled for even if the fit is not valid

   Scout4B::Displacements result;

   const reco::Vertex *bestVertex(0);
   int bestVertexIndex(-1);
   const reco::Vertex *bestVertex2(0);
   int bestVertexIndex2(-1);

   if (fit.valid())
   {

      // const auto& vertices = *pvHandle_.product();

      auto candTransientTrack = fit.particle()->refittedTransientTrack();

      // find best matching primary vertex
      double minDistance(999.);
      for (unsigned int i = 0; i < vertices_.size(); ++i)
      {
         const auto &vertex = vertices_.at(i);
         if (closestIn3D)
         {
            auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, vertex);
            if (impactParameter3D.first and impactParameter3D.second.value() < minDistance)
            {
               minDistance = impactParameter3D.second.value();
               bestVertex = &vertex;
               bestVertexIndex = i;
            }
         }
         else
         {
            auto impactParameterZ = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0, 0, 1), vertex);
            double distance = fabs(impactParameterZ.second.value());
            if (impactParameterZ.first and distance < minDistance)
            {
               minDistance = distance;
               bestVertex = &vertex;
               bestVertexIndex = i;
            }
         }
      }

      // find second best vertex
      double minDistance2(999.);
      for (unsigned int i = 0; i < vertices_.size(); ++i)
      {
         const auto &vertex = vertices_.at(i);
         if (closestIn3D)
         {
            auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, vertex);
            if (impactParameter3D.first and impactParameter3D.second.value() < minDistance2 and impactParameter3D.second.value() > minDistance)
            {
               minDistance2 = impactParameter3D.second.value();
               bestVertex2 = &vertex;
               bestVertexIndex2 = i;
            }
         }
         else
         {
            auto impactParameterZ = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0, 0, 1), vertex);
            double distance = fabs(impactParameterZ.second.value());
            if (impactParameterZ.first and distance < minDistance2 and distance > minDistance)
            {
               minDistance2 = distance;
               bestVertex2 = &vertex;
               bestVertexIndex2 = i;
            }
         }
      }
   }

   if (bestVertex)
      result.push_back(Scout4B::Displacement("pv", fit, *bestVertex, bestVertexIndex));
   else
      result.push_back(Scout4B::Displacement("pv"));

   if (bestVertex2)
      result.push_back(Scout4B::Displacement("_pv2", fit, *bestVertex2, bestVertexIndex2));
   else
      result.push_back(Scout4B::Displacement("_pv2"));

   return result;
}

std::pair<std::vector<KinematicFitResult>, std::vector<KinematicFitResult>> Scout4BRecoSecondaryVertexAnalyzer::performVertexFit(
    std::vector<reco::Track> tracks,
    unsigned int NP,
    int MCharge,
    double MMass,
    double MMassErr,
    double MassWin,
    double vProbMin,
    const std::vector<double> &PMass,
    const std::vector<double> &PMassErr,
    const std::vector<int> &PCharge,
    double MMassMin,
    bool doIso,
    double PIso,
    const MagneticField &bFieldHandle,
    const unsigned long long max_loop)
{
   if (PMass.size() != NP || PMassErr.size() != NP || tracks.size() < NP || std::abs(MCharge) > int(NP))
   {
      return std::make_pair(std::vector<KinematicFitResult>(), std::vector<KinematicFitResult>());
      cout << "PMass.size() != NP || PMassErr.size() != NP || tracks.size() < NP || std::abs(MCharge) > int(NP)" << endl;
   }

   std::vector<KinematicFitResult> unconstrainedCandidates;
   std::vector<KinematicFitResult> massConstrainedCandidates;

   std::vector<reco::Track> &trackVec = tracks;
   unsigned int nTracks = trackVec.size();

   std::vector<unsigned int> indices(NP);
   for (unsigned int i = 0; i < NP; ++i)
   {
      indices[i] = i;
   }
   unsigned long long loopCounter = 0;

   bool done = false;
   while (!done)
   {
      loopCounter++;
      if (loopCounter > max_loop)
      {
         cout << "loopCounter > max_loop" << endl;
         break;
      }

      bool process = true;

      int sumCharge = 0;
      for (unsigned int i = 0; i < NP; ++i)
      {
         sumCharge += trackVec[indices[i]].charge();
         if (trackVec[indices[i]].charge() != PCharge[i])
         {
            process = false;
         }
      }
      if (process && sumCharge != MCharge)
      {
         process = false;
      }

      if (process && doIso)
      {
         bool isoFail = false;
         for (unsigned int i = 0; i < NP && !isoFail; ++i)
         {
            for (unsigned int j = i + 1; j < NP; ++j)
            {
               if (deltaR(trackVec[indices[i]], trackVec[indices[j]]) < PIso)
               {
                  isoFail = true;
                  break;
               }
            }
         }
         if (isoFail)
         {
            process = false;
         }
      }

      if (process)
      {
         bool docaFail = false;
         // 遍历所有粒子对
         for (unsigned int i = 0; i < NP && !docaFail; ++i)
         {
            for (unsigned int j = i + 1; j < NP && !docaFail; ++j)
            {
               double doca = distanceOfClosestApproach(
                   &trackVec[indices[i]],
                   &trackVec[indices[j]],
                   bFieldHandle);
               if (doca > 0.05)
                  docaFail = true;
            }
         }
         if (docaFail)
         {
            // cout << "docaFail" << endl;
            // process = false;
         }
      }

      KinematicFitResult resultUnconstrained;
      KinematicFitResult resultConstrained;

      std::vector<reco::Track *> trackPointers;
      for (unsigned int i = 0; i < NP && process; ++i)
      {
         trackPointers.push_back(&trackVec[indices[i]]);
      }

      if (process)
      {
         // Use KinematicFitter for unconstrained fit
         resultUnconstrained = KinematicFitter(trackPointers, PMass, bFieldHandle);

         if (resultUnconstrained.valid())
         {
            resultUnconstrained.charge_ = MCharge;
            resultUnconstrained.dau_charge_ = PCharge;

            if (std::abs(resultUnconstrained.p4().M() - MMass) < MassWin && resultUnconstrained.vtxProb() > vProbMin && resultUnconstrained.p4().M() > MMassMin)
            {
               resultUnconstrained.postprocess(*beamSpot_);
               for (unsigned int i = 0; i < NP; ++i)
               {
                  resultUnconstrained.trackIndices.push_back(indices[i]);
               }
               unconstrainedCandidates.push_back(resultUnconstrained);

               // Use MCKinematicFitter for mass-constrained fit
               resultConstrained = MCKinematicFitter(trackPointers, PMass, MMass, MMassErr, bFieldHandle);

               if (resultConstrained.valid())
               {
                  if (resultConstrained.vtxProb() > vProbMin && resultConstrained.p4().M() > MMassMin && abs(resultConstrained.p4().M() - MMass) < MassWin && abs(resultConstrained.p4().M() - MMass) < 3 * resultConstrained.massErr())
                  {
                     resultConstrained.charge_ = MCharge;
                     resultConstrained.dau_charge_ = PCharge;
                     for (unsigned int i = 0; i < NP; ++i)
                     {
                        resultConstrained.trackIndices.push_back(indices[i]);
                     }
                     massConstrainedCandidates.push_back(resultConstrained);
                  }
               }
            }
         }
      }

      int pos = NP - 1;
      bool found = false;
      while (pos >= 0 && !found)
      {
         // 保存当前值并尝试递增
         unsigned int current = indices[pos];
         current++;

         // 查找下一个可用的唯一值
         while (current < nTracks)
         {
            bool duplicate = false;
            for (int i = 0; i < pos; ++i)
            {
               if (indices[i] == current)
               {
                  duplicate = true;
                  break;
               }
            }
            if (!duplicate)
               break;
            current++;
         }

         if (current < nTracks)
         {
            indices[pos] = current;
            found = true;

            // 填充后续位置为最小可用值
            std::vector<bool> used(nTracks, false);
            for (int i = 0; i <= pos; ++i)
               used[indices[i]] = true;

            for (unsigned int i = pos + 1; i < NP; ++i)
            {
               for (unsigned int j = 0; j < nTracks; ++j)
               {
                  if (!used[j])
                  {
                     indices[i] = j;
                     used[j] = true;
                     break;
                  }
               }
            }
         }
         else
         {
            pos--; // 进位到前一位
         }
      }

      if (!found)
         done = true; // 所有排列生成完毕
   }

   return std::make_pair(unconstrainedCandidates, massConstrainedCandidates);
}

std::pair<std::vector<XFitResult>, std::vector<KinematicFitResult>> Scout4BRecoSecondaryVertexAnalyzer::performXVertexFit(
    std::vector<reco::Track> tracks,
    unsigned int NM,
    const std::vector<std::vector<KinematicFitResult>> &vertexs,
    int XCharge,
    double XMass,
    double XMassWin,
    double vProbMin,
    const std::vector<double> &MMass,
    const std::vector<double> &MMassErr,
    double XMassMin,
    bool doIso,
    double PIso,
    const MagneticField &bFieldHandle,
    const unsigned long long max_loop)
{
   // Input checks
   if (MMass.size() != NM || MMassErr.size() != NM || std::abs(XCharge) > int(NM) || vertexs.size() < NM)
   {
      return std::make_pair(std::vector<XFitResult>(), std::vector<KinematicFitResult>());
   }
   std::vector<reco::Track> &trackVec = tracks;
   if (trackVec.empty())
      return std::make_pair(std::vector<XFitResult>(), std::vector<KinematicFitResult>());

   std::vector<XFitResult> xCandidates;
   std::vector<KinematicFitResult> xNoMassFitResults; // 用于存储不加Mass Constraint的X拟合结果

   // Cartesian product over each vertex group
   std::vector<unsigned int> candIndices(NM, 0);

   unsigned long long candLoopCounter = 0;
   bool doneCandidates = false;
   while (!doneCandidates)
   {
      if (candLoopCounter > max_loop)
      {
         cout << "candLoopCounter > max_loop" << endl;
         break;
      }
      candLoopCounter++;

      bool processCandidate = true;
      std::vector<KinematicFitResult> selectedM;
      // Select one candidate from each group
      for (unsigned int i = 0; i < NM; ++i)
      {
         if (i >= vertexs.size() || vertexs[i].empty() || candIndices[i] >= vertexs[i].size())
         {
            processCandidate = false;
            break;
         }
         selectedM.push_back(vertexs[i][candIndices[i]]);
      }
      if (processCandidate && selectedM.size() != NM)
         processCandidate = false;

      // Check charge sum
      if (processCandidate)
      {
         int sumCharge = 0;
         for (unsigned int i = 0; i < NM; ++i)
         {
            sumCharge += selectedM[i].charge_;
         }
         if (sumCharge != XCharge)
            processCandidate = false;
      }

      // Check for duplicate track indices
      if (processCandidate)
      {
         std::set<unsigned int> allTrackIndices;
         bool duplicateFound = false;
         for (unsigned int i = 0; i < NM && !duplicateFound; ++i)
         {
            for (auto idx : selectedM[i].trackIndices)
            {
               if (allTrackIndices.count(idx) > 0)
               {
                  duplicateFound = true;
                  break;
               }
               allTrackIndices.insert(idx);
            }
         }
         if (duplicateFound)
            processCandidate = false;
      }

      // Check ISO and doca between any 2 tracks
      if (processCandidate)
      {
         bool isoFail = false;
         bool docaFail = false;
         for (unsigned int i = 0; i < NM && !docaFail; ++i)
         {
            for (unsigned int j = i + 1; j < NM && !docaFail; ++j)
            {
               for (auto idx1 : selectedM[i].trackIndices)
               {
                  for (auto idx2 : selectedM[j].trackIndices)
                  {
                     if (deltaR(trackVec[idx1], trackVec[idx2]) < PIso)
                     {
                        isoFail = true;
                        break;
                     }
                     double doca = distanceOfClosestApproach(&trackVec[idx1], &trackVec[idx2], bFieldHandle);
                     if (doca > 0.05)
                     {
                        docaFail = true;
                        break;
                     }
                  }
               }
            }
         }
         if ((isoFail && doIso) || docaFail)
            processCandidate = false;
      }

      /*
      // Fit all M particles together to X particle without Mass Constraint
      // 新增部分：使用KinematicFitter对selectedM所有粒子的D粒子进行X级拟合（不加Mass Constraint）
      std::vector<int> mCharge;
      if (processCandidate)
      {
         KinematicFitResult xKinematicFitResult_noMass;
         RefCountedKinematicTree xKinematicFitTree_noMass;
         bool xFitNoMassSuccess = true;
         try
         {
            vector<RefCountedKinematicParticle> mParticles;
            for (unsigned int i = 0; i < NM; ++i)
            {
               mParticles.push_back(selectedM[i].particle());
               mCharge.push_back(selectedM[i].charge_);
            }

            try
            {
               xKinematicFitTree_noMass = KinematicParticleVertexFitter().fit(mParticles);
            }
            catch (const std::exception &e)
            {
               xFitNoMassSuccess = false;
            }

            if (xFitNoMassSuccess && xKinematicFitTree_noMass->isValid())
            {
               xKinematicFitResult_noMass.set_tree(xKinematicFitTree_noMass);
               xKinematicFitResult_noMass.postprocess(*beamSpot_);
               // 设置track indices为所有M候选的track indices合并结果
               std::vector<unsigned int> combinedIndices;
               for (unsigned int i = 0; i < selectedM.size(); ++i)
               {
                  for (auto idx : selectedM[i].trackIndices)
                     combinedIndices.push_back(idx);
               }
               xKinematicFitResult_noMass.trackIndices = combinedIndices;
               xKinematicFitResult_noMass.charge_ = XCharge;
               xKinematicFitResult_noMass.dau_charge_ = mCharge;
            }
            else
            {
               xFitNoMassSuccess = false;
            }
         }
         catch (const std::exception &e)
         {
            xFitNoMassSuccess = false;
         }
         if (xFitNoMassSuccess)
         {
            xNoMassFitResults.push_back(xKinematicFitResult_noMass);
         }
      }
      */

      // Apply mass constraint fitting for each M candidate
      bool massConstraintSuccess = true;
      std::vector<std::pair<std::vector<TLorentzVector>, std::vector<int>>> mDaughter_result;
      std::vector<std::pair<TLorentzVector, int>> mMother_result;
      std::vector<double> mMassErr_result;
      std::vector<std::vector<unsigned int>> mTrackIndices_result;
      std::vector<RefCountedKinematicParticle> fittedMParticles;
      XFitResult xResult;

      for (unsigned int i = 0; i < NM && massConstraintSuccess && processCandidate; ++i)
      {
         const KinematicFitResult &mCandidate = selectedM[i];
         KinematicFitResult mKinematicFitResult;

         std::vector<double> parMasses;
         std::vector<reco::Track *> parTracks;
         // Construct each daughter particle (initially used for fitting)
         for (unsigned int j = 0; j < mCandidate.trackIndices.size() && massConstraintSuccess; ++j)
         {
            unsigned int trackIdx = mCandidate.trackIndices[j];
            if (trackIdx >= trackVec.size())
            {
               massConstraintSuccess = false;
               break;
            }
            parMasses.push_back(mCandidate.dau_p4(j).M());
            parTracks.push_back(&trackVec[trackIdx]);
         }

         if (!massConstraintSuccess)
         {
            break;
         }

         // Use MCKinematicFitter to apply mass constraint
         RefCountedKinematicTree constrainedTree;
         bool mcProcess = true;

         try
         {
            mKinematicFitResult = MCKinematicFitter(parTracks, parMasses, MMass[i], MMassErr[i], bFieldHandle);
         }
         catch (const std::exception &e)
         {
            mcProcess = false;
         }
         if (mcProcess && mKinematicFitResult.valid())
         {
            constrainedTree = mKinematicFitResult.tree();
            mKinematicFitResult.postprocess(*beamSpot_);
         }
         else
         {
            mcProcess = false;
         }

         if (mcProcess && (!constrainedTree || !constrainedTree->isValid()))
            mcProcess = false;

         if (mcProcess)
         {
            constrainedTree->movePointerToTheTop();
            RefCountedKinematicParticle mMotherParticle = constrainedTree->currentParticle();
         }
         if (!mcProcess)
         {
            massConstraintSuccess = false;
            break;
         }

         // Check the validity of the constrained tree
         double vProbMC = -1;
         constrainedTree->movePointerToTheTop();
         RefCountedKinematicParticle mMotherParticleMC = constrainedTree->currentParticle();
         RefCountedKinematicVertex mVertexMC = constrainedTree->currentDecayVertex();
         if (!mMotherParticleMC->currentState().isValid() || !mVertexMC->vertexIsValid() ||
             mVertexMC->chiSquared() < 0 || mVertexMC->degreesOfFreedom() <= 0 ||
             mVertexMC->chiSquared() > 9999.9)
         {
            massConstraintSuccess = false;
            break;
         }
         else
         {
            vProbMC = ChiSquaredProbability(mVertexMC->chiSquared(), mVertexMC->degreesOfFreedom());
         }
         if (vProbMC < vProbMin)
         {
            massConstraintSuccess = false;
            break;
         }

         // Extract daughter particle information
         std::vector<TLorentzVector> fittedDaughters;
         constrainedTree->movePointerToTheFirstChild();
         for (unsigned int j = 0; j < mCandidate.trackIndices.size(); ++j)
         {
            RefCountedKinematicParticle fittedParticle = constrainedTree->currentParticle();
            if (!fittedParticle->currentState().isValid())
            {
               fittedDaughters.clear();
               massConstraintSuccess = false;
               break;
            }
            auto par = fittedParticle->currentState().kinematicParameters();
            TVector3 mom(par.momentum().x(), par.momentum().y(), par.momentum().z());
            double energy = std::sqrt(mom.Mag2() +
                                      fittedParticle->currentState().mass() * fittedParticle->currentState().mass());
            TLorentzVector p4;
            p4.SetPxPyPzE(mom.x(), mom.y(), mom.z(), energy);
            fittedDaughters.push_back(p4);
            if (j < mCandidate.trackIndices.size() - 1)
               constrainedTree->movePointerToTheNextChild();
         }
         if (fittedDaughters.empty())
         {
            massConstraintSuccess = false;
            break;
         }
         mDaughter_result.push_back(std::make_pair(fittedDaughters, mCandidate.dau_charge_));
         auto mPar = mMotherParticleMC->currentState().kinematicParameters();
         TVector3 mMom(mPar.momentum().x(), mPar.momentum().y(), mPar.momentum().z());
         double mEnergy = std::sqrt(mMom.Mag2() + mMotherParticleMC->currentState().mass() * mMotherParticleMC->currentState().mass());
         TLorentzVector mP4;
         mP4.SetPxPyPzE(mMom.x(), mMom.y(), mMom.z(), mEnergy);
         mMother_result.push_back(std::make_pair(mP4, mCandidate.charge_));
         mMassErr_result.push_back(std::sqrt(mMotherParticleMC->currentState().kinematicParametersError().matrix()(6, 6)));
         mTrackIndices_result.push_back(mCandidate.trackIndices);

         constrainedTree->movePointerToTheTop();
         RefCountedKinematicParticle constrainedM = constrainedTree->currentParticle();
         fittedMParticles.push_back(constrainedM);
      }

      // Use new fitted M particles for X-level vertex fitting with mass constraint
      bool xVertexSuccess = true;
      if (processCandidate && massConstraintSuccess)
      {
         RefCountedKinematicTree xVertexFitTree;
         try
         {
            xVertexFitTree = KinematicParticleVertexFitter().fit(fittedMParticles);
         }
         catch (const std::exception &e)
         {
            xVertexSuccess = false;
         }
         if (xVertexSuccess && xVertexFitTree->isValid() && xVertexFitTree->currentDecayVertex()->vertexIsValid())
         {
            xVertexFitTree->movePointerToTheTop();
            RefCountedKinematicParticle xMotherParticle = xVertexFitTree->currentParticle();
            RefCountedKinematicVertex xVertex = xVertexFitTree->currentDecayVertex();
            double vProbX = ChiSquaredProbability(xVertex->chiSquared(), xVertex->degreesOfFreedom());
            double fittedXMass = xMotherParticle->currentState().mass();
            if (vProbX >= vProbMin && fittedXMass >= XMassMin &&
                std::fabs(fittedXMass - XMass) <= XMassWin &&
                xMotherParticle->currentState().kinematicParametersError().matrix()(6, 6) >= 0)
            {
               auto par = xMotherParticle->currentState().kinematicParameters();
               TVector3 mom(par.momentum().x(), par.momentum().y(), par.momentum().z());
               double energy = std::sqrt(mom.Mag2() + fittedXMass * fittedXMass);
               TLorentzVector xP4;
               xP4.SetPxPyPzE(mom.x(), mom.y(), mom.z(), energy);

               KinematicFitResult xKinematicFitResult;
               vector<reco::Track *> tracks;
               for (unsigned int i = 0; i < NM; i++)
               {
                  for (unsigned int j = 0; j < mTrackIndices_result[i].size(); j++)
                  {
                     tracks.push_back(&trackVec[mTrackIndices_result[i][j]]);
                  }
               }
               xKinematicFitResult.tracks = tracks;
               xKinematicFitResult.set_tree(xVertexFitTree);
               if (xKinematicFitResult.valid())
                  xKinematicFitResult.postprocess(*beamSpot_);

               // Fill XFitResult
               xResult.mDaughter = mDaughter_result;
               xResult.mMother = mMother_result;
               xResult.xMotherP4 = xP4;
               xResult.xMotherCharge = XCharge;
               xResult.mMassError = mMassErr_result;
               xResult.xMassError = std::sqrt(xMotherParticle->currentState().kinematicParametersError().matrix()(6, 6));
               xResult.normChi2 = xVertex->chiSquared() / (double)xVertex->degreesOfFreedom();
               xResult.mIndices = candIndices;
               xResult.mTrackIndices = mTrackIndices_result;
               if (xKinematicFitResult.valid())
               {
                  xResult.xlxy = xKinematicFitResult.lxy();
                  xResult.xlxyErr = xKinematicFitResult.lxyErr();
                  xResult.xlz = xKinematicFitResult.lz();
                  xResult.xlzErr = xKinematicFitResult.lzErr();
                  xResult.xl = xKinematicFitResult.l();
                  xResult.xlErr = xKinematicFitResult.lErr();
                  xResult.xsigLxy = xKinematicFitResult.sigLxy();
                  xResult.xsigLz = xKinematicFitResult.sigLz();
                  xResult.xsigL = xKinematicFitResult.sigL();
                  xResult.xalphaBS = xKinematicFitResult.alphaBS();
                  xResult.xalphaBSErr = xKinematicFitResult.alphaBSErr();
                  xResult.xalphaBSXY = xKinematicFitResult.alphaBSXY();
                  xResult.xalphaBSXYErr = xKinematicFitResult.alphaBSXYErr();
               }
               xVertexSuccess = true;
               xCandidates.push_back(xResult);
            }
         }
      }

      // Generate next combination of candidates
      int posCand = NM - 1;
      while (posCand >= 0 && candIndices[posCand] == vertexs[posCand].size() - 1)
         posCand--;
      if (posCand < 0)
         doneCandidates = true;
      else
      {
         candIndices[posCand]++;
         for (unsigned int j = posCand + 1; j < NM; j++)
            candIndices[j] = 0;
      }
   }

   return std::make_pair(xCandidates, xNoMassFitResults);
}

void Scout4BRecoSecondaryVertexAnalyzer::clearVars()
{
   // Clear All Variables
   nXcand = 0;
   nMmu = 0;
   nMtrk = 0;

   for (unsigned int i = 0; i < 5; i++)
   {
      xMass[i] = 0;
      xMassError[i] = 0;
      xPt[i] = 0;
      xEta[i] = 0;
      xPhi[i] = 0;
      xChi2[i] = 0;
      xCharge[i] = 0;

      xMassNMC[i] = 0;
      xMassErrorNMC[i] = 0;
      xPtNMC[i] = 0;
      xEtaNMC[i] = 0;
      xPhiNMC[i] = 0;
      xChi2NMC[i] = 0;
      xChargeNMC[i] = 0;

      xlxy[i] = 0;
      xlxyErr[i] = 0;
      xlz[i] = 0;
      xlzErr[i] = 0;
      xl[i] = 0;
      xlErr[i] = 0;
      xsigLxy[i] = 0;
      xsigLz[i] = 0;
      xsigL[i] = 0;
      xalphaBS[i] = 0;
      xalphaBSErr[i] = 0;
      xalphaBSXY[i] = 0;
      xalphaBSXYErr[i] = 0;

      xlxyNMC[i] = 0;
      xlxyErrNMC[i] = 0;
      xlzNMC[i] = 0;
      xlzErrNMC[i] = 0;
      xlNMC[i] = 0;
      xlErrNMC[i] = 0;
      xsigLxyNMC[i] = 0;
      xsigLzNMC[i] = 0;
      xsigLNMC[i] = 0;
      xalphaBSNMC[i] = 0;
      xalphaBSErrNMC[i] = 0;
      xalphaBSXYNMC[i] = 0;
      xalphaBSXYErrNMC[i] = 0;

      for (unsigned int j = 0; j < 4; j++)
      {
         mMass[i][j] = 0;
         mMassError[i][j] = 0;
         mPt[i][j] = 0;
         mEta[i][j] = 0;
         mPhi[i][j] = 0;
         mCharge[i][j] = 0;

         mMassNMC[i][j] = 0;
         mPtNMC[i][j] = 0;
         mEtaNMC[i][j] = 0;
         mPhiNMC[i][j] = 0;
         mChargeNMC[i][j] = 0;

         NMC_mMassMu[i][j] = 0;
         NMC_mMassTrk[i][j] = 0;
         NMC_mMassErrMu[i][j] = 0;
         NMC_mMassErrTrk[i][j] = 0;
         NMC_mPtMu[i][j] = 0;
         NMC_mPtTrk[i][j] = 0;
         NMC_mEtaMu[i][j] = 0;
         NMC_mEtaTrk[i][j] = 0;
         NMC_mPhiMu[i][j] = 0;
         NMC_mPhiTrk[i][j] = 0;
         NMC_mChi2Mu[i][j] = 0;
         NMC_mChi2Trk[i][j] = 0;
         NMC_mChargeMu[i][j] = 0;
         NMC_mChargeTrk[i][j] = 0;
         NMC_mChi2Mu[i][j] = 0;
         NMC_mChi2Trk[i][j] = 0;

         NMC_mlxyMu[i][j] = 0;
         NMC_mlxyErrMu[i][j] = 0;
         NMC_mlzMu[i][j] = 0;
         NMC_mlzErrMu[i][j] = 0;
         NMC_mlMu[i][j] = 0;
         NMC_mlErrMu[i][j] = 0;
         NMC_msigLxyMu[i][j] = 0;
         NMC_msigLzMu[i][j] = 0;
         NMC_msigLMu[i][j] = 0;
         NMC_malphaBSMu[i][j] = 0;
         NMC_malphaBSErrMu[i][j] = 0;
         NMC_malphaBSXYMu[i][j] = 0;
         NMC_malphaBSXYErrMu[i][j] = 0;

         NMC_mlxyTrk[i][j] = 0;
         NMC_mlxyErrTrk[i][j] = 0;
         NMC_mlzTrk[i][j] = 0;
         NMC_mlzErrTrk[i][j] = 0;
         NMC_mlTrk[i][j] = 0;
         NMC_mlErrTrk[i][j] = 0;
         NMC_msigLxyTrk[i][j] = 0;
         NMC_msigLzTrk[i][j] = 0;
         NMC_msigLTrk[i][j] = 0;
         NMC_malphaBSTrk[i][j] = 0;
         NMC_malphaBSErrTrk[i][j] = 0;
         NMC_malphaBSXYTrk[i][j] = 0;
         NMC_malphaBSXYErrTrk[i][j] = 0;

         nMmucand[j] = 0;
         nMtrkcand[j] = 0;
         nDmu[j] = 0;
         nDtrk[j] = 0;

         for (unsigned int k = 0; k < 8; k++)
         {
            dMass[i][j][k] = 0;
            dPt[i][j][k] = 0;
            dEta[i][j][k] = 0;
            dPhi[i][j][k] = 0;
            dCharge[i][j][k] = 0;

            NMC_dMassMu[i][j][k] = 0;
            NMC_dMassTrk[i][j][k] = 0;
            NMC_dPtMu[i][j][k] = 0;
            NMC_dPtTrk[i][j][k] = 0;
            NMC_dEtaMu[i][j][k] = 0;
            NMC_dEtaTrk[i][j][k] = 0;
            NMC_dPhiMu[i][j][k] = 0;
            NMC_dPhiTrk[i][j][k] = 0;
            NMC_dChargeMu[i][j][k] = 0;
            NMC_dChargeTrk[i][j][k] = 0;
         }
      }
   }
}

//
// Fill Descriptions
//
void Scout4BRecoSecondaryVertexAnalyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
   edm::ParameterSetDescription desc;

   desc.addUntracked<std::string>("treename", "Scout4B");
   desc.addUntracked<bool>("multiM", false);
   desc.add<edm::InputTag>("recoTrackMuon", edm::InputTag("recoTrackMuon"));
   desc.add<edm::InputTag>("recoTrack", edm::InputTag("recoTrack"));
   desc.add<edm::InputTag>("recoVertex", edm::InputTag("recoVertex"));
   desc.add<edm::InputTag>("beamSpot", edm::InputTag("beamSpot"));
   desc.add<edm::InputTag>("TriggerResults");
   desc.add<std::vector<std::string>>("FilterNames");
   desc.add<std::vector<double>>("MMuMass");
   desc.add<std::vector<double>>("MMuMassErr");
   desc.add<std::vector<double>>("MTrkMass");
   desc.add<std::vector<double>>("MTrkMassErr");
   desc.addUntracked<double>("MuMass", 0.1056583755);
   desc.addUntracked<double>("MuMassErr", 0.0000000023);
   desc.add<std::vector<int>>("MuChargeFlat");
   desc.add<std::vector<double>>("TrkMassFlat");
   desc.add<std::vector<double>>("TrkMassErrFlat");
   desc.add<std::vector<int>>("TrkChargeFlat");
   desc.add<unsigned int>("NMmu");
   desc.add<unsigned int>("NMtrk");
   desc.add<std::vector<unsigned int>>("NPmu");
   desc.add<std::vector<unsigned int>>("NPtrk");
   desc.addUntracked<double>("XMass", 0.0);
   desc.addUntracked<double>("XMassWin", 1000.0);
   desc.addUntracked<double>("MMassWin", 0.5);
   desc.add<int>("XCharge");
   desc.add<std::vector<int>>("MChargeMu");
   desc.add<std::vector<int>>("MChargetrk");
   desc.addUntracked<double>("PIso", 0.02);
   desc.addUntracked<bool>("doIso", true);
   desc.addUntracked<bool>("doTrigger", false);
   desc.addUntracked<double>("vProbMin", 1e-3);
   desc.addUntracked<unsigned long long>("maxLoop", 1000000000);
   desc.addUntracked<double>("MMassMin", 1e-2);
   desc.addUntracked<double>("XMassMin", 1e-2);

   descriptions.add("Scout4BRecoSecondaryVertexAnalyzer", desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Scout4BRecoSecondaryVertexAnalyzer);