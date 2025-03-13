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
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

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

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace muon;
using namespace trigger;

//
// struct definition
//
struct VertexFitResult
{
   std::vector<TLorentzVector> daughterP4; // NP个子粒子的四动量（均来自拟合结果）
   TLorentzVector motherP4;                // 母粒子的四动量（来自顶点拟合后的状态）
   double massError;                       // 拟合得到的母粒子质量误差（sigma）
   double normChi2;                        // 归一化chi2 (chi2/ndf)
   int motherCharge;                       // 母粒子的电荷（来自拟合结果）
   std::vector<int> daughterCharge;        // 每个子粒子的电荷（来自拟合结果）
   std::vector<unsigned int> trackIndices; // 记录参与质量约束拟合的tracks在输入tracks中的索引

   VertexFitResult()
   {
      daughterP4.clear();
      motherP4.SetXYZM(0, 0, 0, 0);
      massError = 0;
      normChi2 = 0;
      motherCharge = 0;
      daughterCharge.clear();
      trackIndices.clear();
   }
};

struct XFitResult
{
   // 每个 M 候选经过质量约束拟合后得到的子粒子四动量（分组保存，组数为 NM）
   std::vector<std::pair<std::vector<TLorentzVector>, std::vector<int>>> mDaughter;
   // 每个 M 候选母粒子的四动量及电荷
   std::vector<std::pair<TLorentzVector, int>> mMother;
   std::vector<double> mMassError;
   // X 粒子的拟合结果：四动量和电荷
   TLorentzVector xMotherP4;
   int xMotherCharge;
   // X 粒子拟合的质量误差和归一化 chi2
   double xMassError;
   double normChi2;

   // 新增信息：记录组成该 X 候选的各个 M 候选在 vertexs 中的位置
   // 对应每个 M 候选：所属组号和该组内候选索引
   std::vector<unsigned int> mGroupIndices;
   std::vector<unsigned int> mCandidateIndices;
   // 对应每个 M 候选，其 particle 在原始 tracks 中的索引（从 VertexFitResult 中获得）
   std::vector<std::vector<unsigned int>> mTrackIndices;

   XFitResult()
   {
      mDaughter.clear();
      mMother.clear();
      xMotherP4.SetXYZM(0, 0, 0, 0);
      xMotherCharge = 0;
      xMassError = 0;
      normChi2 = 0;
      mGroupIndices.clear();
      mCandidateIndices.clear();
      mTrackIndices.clear();
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

   // Difine core functions for sv fitting

   std::pair<std::vector<VertexFitResult>, std::vector<VertexFitResult>> performVertexFit(
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

   std::vector<XFitResult> performXVertexFit(
       std::vector<reco::Track> tracks,
       unsigned int NM,
       const std::vector<std::vector<VertexFitResult>> &vertexs,
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

   // Define output parameters
   std::string file_name;
   edm::Service<TFileService> fs;

   TTree *Scout4BTree;
   TTree *NMC_Scout4BTree;

   // Define branches
   UInt_t run;
   ULong64_t event;
   UInt_t lumiblock;
   UInt_t trigger;

   UInt_t nXcand;
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

   // M particle branches
   double mMass[5][4];
   double mMassError[5][4];
   double mPt[5][4];
   double mEta[5][4];
   double mPhi[5][4];
   double mChi2[5][4];
   int mCharge[5][4];

   double NMC_mMassMu[5][4];
   double NMC_mMassErrMu[5][4];
   double NMC_mPtMu[5][4];
   double NMC_mEtaMu[5][4];
   double NMC_mPhiMu[5][4];
   double NMC_mChi2Mu[5][4];
   int NMC_mChargeMu[5][4];

   double NMC_mMassTrk[5][4];
   double NMC_mMassErrTrk[5][4];
   double NMC_mPtTrk[5][4];
   double NMC_mEtaTrk[5][4];
   double NMC_mPhiTrk[5][4];
   double NMC_mChi2Trk[5][4];
   int NMC_mChargeTrk[5][4];

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
   NMC_Scout4BTree->Branch("NMC_M_MassTrk", NMC_mMassTrk, "NMC_M_MassTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_MassErrTrk", NMC_mMassErrTrk, "NMC_M_MassErrTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_PtTrk", NMC_mPtTrk, "NMC_M_PtTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_EtaTrk", NMC_mEtaTrk, "NMC_M_EtaTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_PhiTrk", NMC_mPhiTrk, "NMC_M_PhiTrk[5][4]/D");
   NMC_Scout4BTree->Branch("NMC_M_ChargeTrk", NMC_mChargeTrk, "NMC_M_ChargeTrk[5][4]/I");

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

   const MagneticField &bFieldHandle = iSetup.getData(magneticFieldToken_);

   // Collect inputs
   Handle<std::vector<reco::Track>> recoTrackMuon;
   iEvent.getByToken(input_recoTrackMuon_token_, recoTrackMuon);
   Handle<std::vector<reco::Track>> recoTrack;
   iEvent.getByToken(input_recoTrack_token_, recoTrack);

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
      //  cout << "Input collections too small: " << recoTrackMuon->size() << " " << recoTrack->size() << endl;
      return;
   }

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
      return;

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
   vector<vector<VertexFitResult>> mCandidatesMu;
   vector<vector<VertexFitResult>> NMC_mCandidatesMu;
   for (unsigned int i = 0; i < NMmu_; i++)
   {
      // 构造用于拟合的 particle 质量及误差向量，长度为 NPmu_[i]，取输入的 MuMass 参数（保持不变）
      vector<double> PMass_mu(NPmu_[i], MuMass_);
      vector<double> PMassErr_mu(NPmu_[i], MuMassErr_);
      std::pair<std::vector<VertexFitResult>, std::vector<VertexFitResult>> res = performVertexFit(goodTrackMuons,
                                                                                                   NPmu_[i],
                                                                                                   MChargeMu_[i],
                                                                                                   MMuMass_[i],
                                                                                                   MMuMassErr_[i],
                                                                                                   MMassWin_,
                                                                                                   vProbMin_,
                                                                                                   PMass_mu,
                                                                                                   PMassErr_mu,
                                                                                                   MuCharge_v_[i],
                                                                                                   MMassMin_,
                                                                                                   doIso_,
                                                                                                   PIso_,
                                                                                                   bFieldHandle,
                                                                                                   maxLoop_);
      if (res.second.empty())
         return;
      mCandidatesMu.push_back(res.second);
      NMC_mCandidatesMu.push_back(res.first);
   }

   if (mCandidatesMu.size() < NMmu_ || NMC_mCandidatesMu.size() < NMmu_)
      return;

   // Step 2: Fit M particles from trks
   vector<vector<VertexFitResult>> mCandidatesTrk;
   vector<vector<VertexFitResult>> NMC_mCandidatesTrk;
   for (unsigned int i = 0; i < NMtrk_; i++)
   {
      std::pair<std::vector<VertexFitResult>, std::vector<VertexFitResult>> res = performVertexFit(goodTracks,
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
                                                                                                   bFieldHandle,
                                                                                                   maxLoop_);
      if (res.second.empty())
         return;
      mCandidatesTrk.push_back(res.second);
      NMC_mCandidatesTrk.push_back(res.first);
   }

   if (mCandidatesTrk.size() < NMtrk_ || NMC_mCandidatesTrk.size() < NMtrk_)
      return;

   // Sort the M candidates by M's pT
   for (auto &vec : mCandidatesMu)
   {
      std::sort(vec.begin(), vec.end(), [](const VertexFitResult &a, const VertexFitResult &b)
                { return a.motherP4.Pt() > b.motherP4.Pt(); });
   }
   for (auto &vec : mCandidatesTrk)
   {
      std::sort(vec.begin(), vec.end(), [](const VertexFitResult &a, const VertexFitResult &b)
                { return a.motherP4.Pt() > b.motherP4.Pt(); });
   }

   // Step 2.5: Fill NMC_Scout4BTree with NMC_mCandidatesMu and NMC_mCandidatesTrk

   // Fill NMC M muon candidates
   for (unsigned int j = 0; j < NMmu_; j++)
   {
      nMmucand[j] = min(int(NMC_mCandidatesMu[j].size()), 5);
      for (unsigned int i = 0; i < nMmucand[j]; i++)
      {
         NMC_mMassMu[i][j] = NMC_mCandidatesMu[j][i].motherP4.M();
         NMC_mMassErrMu[i][j] = NMC_mCandidatesMu[j][i].massError;
         NMC_mPtMu[i][j] = NMC_mCandidatesMu[j][i].motherP4.Pt();
         NMC_mEtaMu[i][j] = NMC_mCandidatesMu[j][i].motherP4.Eta();
         NMC_mPhiMu[i][j] = NMC_mCandidatesMu[j][i].motherP4.Phi();
         NMC_mChargeMu[i][j] = NMC_mCandidatesMu[j][i].motherCharge;

         for (unsigned int k = 0; k < NMC_mCandidatesMu[j][i].daughterP4.size(); k++)
         {
            NMC_dMassMu[i][j][k] = NMC_mCandidatesMu[j][i].daughterP4[k].M();
            NMC_dPtMu[i][j][k] = NMC_mCandidatesMu[j][i].daughterP4[k].Pt();
            NMC_dEtaMu[i][j][k] = NMC_mCandidatesMu[j][i].daughterP4[k].Eta();
            NMC_dPhiMu[i][j][k] = NMC_mCandidatesMu[j][i].daughterP4[k].Phi();
            NMC_dChargeMu[i][j][k] = NMC_mCandidatesMu[j][i].daughterCharge[k];
         }
      }
   }
   // Fill NMC M track candidates
   for (unsigned int j = 0; j < NMtrk_; j++)
   {
      nMtrkcand[j] = min(int(NMC_mCandidatesTrk[j].size()), 5);
      for (unsigned int i = 0; i < nMtrkcand[j]; i++)
      {
         NMC_mMassTrk[i][j] = NMC_mCandidatesTrk[j][i].motherP4.M();
         NMC_mMassErrTrk[i][j] = NMC_mCandidatesTrk[j][i].massError;
         NMC_mPtTrk[i][j] = NMC_mCandidatesTrk[j][i].motherP4.Pt();
         NMC_mEtaTrk[i][j] = NMC_mCandidatesTrk[j][i].motherP4.Eta();
         NMC_mPhiTrk[i][j] = NMC_mCandidatesTrk[j][i].motherP4.Phi();
         NMC_mChargeTrk[i][j] = NMC_mCandidatesTrk[j][i].motherCharge;

         for (unsigned int k = 0; k < NMC_mCandidatesTrk[j][i].daughterP4.size(); k++)
         {
            NMC_dMassTrk[i][j][k] = NMC_mCandidatesTrk[j][i].daughterP4[k].M();
            NMC_dPtTrk[i][j][k] = NMC_mCandidatesTrk[j][i].daughterP4[k].Pt();
            NMC_dEtaTrk[i][j][k] = NMC_mCandidatesTrk[j][i].daughterP4[k].Eta();
            NMC_dPhiTrk[i][j][k] = NMC_mCandidatesTrk[j][i].daughterP4[k].Phi();
            NMC_dChargeTrk[i][j][k] = NMC_mCandidatesTrk[j][i].daughterCharge[k];
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
   vector<vector<VertexFitResult>> allMCandidates;
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
   // and edit the index in mCandidatesMu to match the new goodTracks
   for (unsigned int i = 0; i < NMmu_; i++)
   {
      for (auto &v : mCandidatesMu[i])
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

   vector<XFitResult> xCandidates = performXVertexFit(goodTracks, totalM, allMCandidates, XCharge_, XMass_, XMassWin_, vProbMin_, mMassForX, mMassErrForX, XMassMin_, doIso_, PIso_, bFieldHandle, maxLoop_);
   if (xCandidates.empty())
      return;

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

std::pair<std::vector<VertexFitResult>, std::vector<VertexFitResult>> Scout4BRecoSecondaryVertexAnalyzer::performVertexFit(
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
      return std::make_pair(std::vector<VertexFitResult>(), std::vector<VertexFitResult>());
   }

   std::vector<VertexFitResult> unconstrainedCandidates;
   std::vector<VertexFitResult> massConstrainedCandidates;

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
         break;

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

      std::vector<RefCountedKinematicParticle> particles;
      if (process)
      {
         KinematicParticleFactoryFromTransientTrack particleFactory;
         float chi = 0.0;
         float ndf = 0.0;
         bool particleError = false;
         for (unsigned int i = 0; i < NP && !particleError; ++i)
         {
            const reco::Track &trk = trackVec[indices[i]];
            TransientTrack tTrk(trk, &bFieldHandle);
            try
            {
               ParticleMass PMassTmp = PMass[i];
               float PMassErrTmp = PMassErr[i];
               RefCountedKinematicParticle particle = particleFactory.particle(tTrk, PMassTmp, chi, ndf, PMassErrTmp);
               particles.push_back(particle);
            }
            catch (std::exception &e)
            {
               particleError = true;
            }
         }
         if (particleError || particles.size() != NP)
         {
            process = false;
         }
      }

      RefCountedKinematicTree vertexFitTree;
      if (process)
      {
         KinematicParticleVertexFitter vertexFitter;
         try
         {
            vertexFitTree = vertexFitter.fit(particles);
         }
         catch (std::exception &e)
         {
            process = false;
         }
         if (process)
         {
            if (!vertexFitTree->isValid())
            {
               process = false;
            }
         }
      }

      if (process)
      {
         vertexFitTree->movePointerToTheTop();
         RefCountedKinematicParticle motherParticle = vertexFitTree->currentParticle();
         RefCountedKinematicVertex vertex = vertexFitTree->currentDecayVertex();
         double vProb = -1;
         double fittedMass = -999;

         if (!motherParticle->currentState().isValid() || !vertex->vertexIsValid())
         {
            process = false;
         }

         if (process && (vertex->chiSquared() < 0 || vertex->degreesOfFreedom() <= 0 || vertex->chiSquared() > 9999.9))
         {
            process = false;
         }
         else
         {
            vProb = ChiSquaredProbability(vertex->chiSquared(), vertex->degreesOfFreedom());
         }
         if (process && vProb < vProbMin)
         {
            process = false;
         }
         else
         {
            fittedMass = motherParticle->currentState().mass();
         }
         if (process && (fittedMass < MMassMin || std::fabs(fittedMass - MMass) > MassWin))
         {
            process = false;
         }

         if (process && motherParticle->currentState().kinematicParametersError().matrix()(6, 6) < 0)
         {
            process = false;
         }
         double sigma = process ? std::sqrt(motherParticle->currentState().kinematicParametersError().matrix()(6, 6)) : 999;

         if (process)
         {
            VertexFitResult resultUnconstrained;
            std::vector<TLorentzVector> daughtersP4;
            std::vector<int> daughtersCharge;
            bool daughtersValid = true;

            vertexFitTree->movePointerToTheFirstChild();
            for (unsigned int i = 0; i < NP && daughtersValid; ++i)
            {
               RefCountedKinematicParticle daughter = vertexFitTree->currentParticle();
               if (!daughter->currentState().isValid())
               {
                  daughtersValid = false;
               }
               else
               {
                  auto par = daughter->currentState().kinematicParameters();
                  TVector3 momentum(par.momentum().x(), par.momentum().y(), par.momentum().z());
                  double energy = std::sqrt(momentum.Mag2() + daughter->currentState().mass() * daughter->currentState().mass());
                  TLorentzVector p4;
                  p4.SetPxPyPzE(momentum.x(), momentum.y(), momentum.z(), energy);
                  daughtersP4.push_back(p4);
                  // Store charge from track
                  daughtersCharge.push_back(trackVec[indices[i]].charge());
               }
               if (i < NP - 1)
               {
                  vertexFitTree->movePointerToTheNextChild();
               }
            }

            if (daughtersValid)
            {
               resultUnconstrained.daughterP4 = daughtersP4;
               resultUnconstrained.daughterCharge = daughtersCharge;

               auto par = motherParticle->currentState().kinematicParameters();
               TVector3 momentum(par.momentum().x(), par.momentum().y(), par.momentum().z());
               double energy = std::sqrt(momentum.Mag2() + fittedMass * fittedMass);
               TLorentzVector p4;
               p4.SetPxPyPzE(momentum.x(), momentum.y(), momentum.z(), energy);
               resultUnconstrained.motherP4 = p4;
               // Calculate charge from track
               resultUnconstrained.motherCharge = MCharge;
               resultUnconstrained.massError = sigma;
               resultUnconstrained.normChi2 = vertex->chiSquared() / (double)vertex->degreesOfFreedom();

               unconstrainedCandidates.push_back(resultUnconstrained);

               if (std::fabs(fittedMass - MMass) < 3 * sigma)
               {
                  bool mcProcess = true;
                  MassKinematicConstraint massConstraint(MMass, MMassErr);
                  KinematicParticleVertexFitter massConstraintFitter;
                  RefCountedKinematicTree constrainedTree;

                  if (mcProcess)
                  {
                     try
                     {
                        constrainedTree = massConstraintFitter.fit(particles);
                     }
                     catch (std::exception &e)
                     {
                        mcProcess = false;
                     }
                  }

                  if (mcProcess && !constrainedTree->isValid())
                  {
                     mcProcess = false;
                  }

                  KinematicParticleFitter csFitter;
                  if (mcProcess)
                  {
                     try
                     {
                        constrainedTree = csFitter.fit(&massConstraint, constrainedTree);
                     }
                     catch (std::exception &e)
                     {
                        mcProcess = false;
                     }
                  }

                  if (mcProcess && !constrainedTree->isValid())
                  {
                     mcProcess = false;
                  }

                  if (mcProcess)
                  {
                     constrainedTree->movePointerToTheTop();
                     RefCountedKinematicParticle motherParticleMC = constrainedTree->currentParticle();
                     RefCountedKinematicVertex vertexMC = constrainedTree->currentDecayVertex();

                     if (!motherParticleMC->currentState().isValid() || !vertexMC->vertexIsValid())
                     {
                        mcProcess = false;
                     }

                     if (mcProcess && (vertexMC->chiSquared() < 0 || vertexMC->degreesOfFreedom() <= 0 || vertexMC->chiSquared() > 9999.9))
                     {
                        mcProcess = false;
                     }

                     double vProbMC = ChiSquaredProbability(vertexMC->chiSquared(), vertexMC->degreesOfFreedom());
                     if (mcProcess && vProbMC < vProbMin)
                     {
                        mcProcess = false;
                     }

                     if (mcProcess)
                     {
                        VertexFitResult resultConstrained;
                        std::vector<TLorentzVector> daughtersP4MC;
                        std::vector<int> daughtersChargeMC;
                        bool mcDaughtersValid = true;

                        constrainedTree->movePointerToTheFirstChild();
                        for (unsigned int i = 0; i < NP && mcDaughtersValid; ++i)
                        {
                           RefCountedKinematicParticle daughter = constrainedTree->currentParticle();
                           if (!daughter->currentState().isValid())
                           {
                              mcDaughtersValid = false;
                           }
                           else
                           {
                              auto par = daughter->currentState().kinematicParameters();
                              TVector3 momentum(par.momentum().x(), par.momentum().y(), par.momentum().z());
                              double energy = std::sqrt(momentum.Mag2() + daughter->currentState().mass() * daughter->currentState().mass());
                              TLorentzVector p4;
                              p4.SetPxPyPzE(momentum.x(), momentum.y(), momentum.z(), energy);
                              daughtersP4MC.push_back(p4);
                              daughtersChargeMC.push_back(trackVec[indices[i]].charge());
                           }
                           if (i < NP - 1)
                           {
                              constrainedTree->movePointerToTheNextChild();
                           }
                        }

                        if (mcDaughtersValid)
                        {
                           resultConstrained.daughterP4 = daughtersP4MC;
                           resultConstrained.daughterCharge = daughtersChargeMC;

                           auto parMC = motherParticleMC->currentState().kinematicParameters();
                           TVector3 momentumMC(parMC.momentum().x(), parMC.momentum().y(), parMC.momentum().z());
                           double energyMC = std::sqrt(momentumMC.Mag2() + MMass * MMass);
                           TLorentzVector p4MC;
                           p4MC.SetPxPyPzE(momentumMC.x(), momentumMC.y(), momentumMC.z(), energyMC);
                           resultConstrained.motherP4 = p4MC;
                           resultConstrained.motherCharge = MCharge;
                           resultConstrained.massError = std::sqrt(motherParticleMC->currentState().kinematicParametersError().matrix()(6, 6));
                           resultConstrained.normChi2 = vertexMC->chiSquared() / (double)vertexMC->degreesOfFreedom();
                           resultConstrained.trackIndices = indices;

                           massConstrainedCandidates.push_back(resultConstrained);
                        }
                     }
                  }
               }
            }
         }
      }

      // Generate next combination
      int pos = NP - 1;
      while (pos >= 0 && indices[pos] == nTracks - NP + pos)
         pos--;
      if (pos < 0)
      {
         done = true;
      }
      else
      {
         indices[pos]++;
         for (unsigned int i = pos + 1; i < NP; ++i)
         {
            indices[i] = indices[i - 1] + 1;
         }
      }
   }

   return std::make_pair(unconstrainedCandidates, massConstrainedCandidates);
}

std::vector<XFitResult> Scout4BRecoSecondaryVertexAnalyzer::performXVertexFit(
    std::vector<reco::Track> tracks,
    unsigned int NM,
    const std::vector<std::vector<VertexFitResult>> &vertexs,
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
   // 输入检查：MMass、MMassErr大小、电荷要求，以及vertexs至少包含NM个组
   if (MMass.size() != NM || MMassErr.size() != NM || std::abs(XCharge) > int(NM) || vertexs.size() < NM)
   {
      return std::vector<XFitResult>();
   }
   std::vector<reco::Track> &trackVec = tracks;
   if (trackVec.empty())
      return std::vector<XFitResult>();

   std::vector<XFitResult> xCandidates;

   // 使用 candIndices 实现各组内候选的笛卡尔组合，vertexs[i]对应第 i 个 M组
   std::vector<unsigned int> candIndices(NM, 0);
   unsigned long long candLoopCounter = 0;
   bool doneCandidates = false;
   while (!doneCandidates)
   {
      if (candLoopCounter > max_loop)
         break;
      candLoopCounter++;

      bool processCandidate = true;
      std::vector<VertexFitResult> selectedM;
      // 从每个 M组中选取一个候选
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

      // 检查电荷和
      if (processCandidate)
      {
         int sumCharge = 0;
         for (unsigned int i = 0; i < NM; ++i)
         {
            sumCharge += selectedM[i].motherCharge;
         }
         if (sumCharge != XCharge)
            processCandidate = false;
      }

      // 检查各候选中的 track 索引是否重复
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

      // 隔离检查（若要求隔离）
      if (processCandidate && doIso)
      {
         bool isoFail = false;
         for (unsigned int i = 0; i < NM && !isoFail; ++i)
         {
            for (unsigned int j = i + 1; j < NM && !isoFail; ++j)
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
                  }
               }
            }
         }
         if (isoFail)
            processCandidate = false;
      }

      // 对每个 M候选进行质量约束拟合（利用两步拟合过程）
      bool massConstraintSuccess = true;
      std::vector<std::pair<std::vector<TLorentzVector>, std::vector<int>>> mDaughter_result;
      std::vector<std::pair<TLorentzVector, int>> mMother_result;
      std::vector<double> mMassErr_result;
      std::vector<std::vector<unsigned int>> mTrackIndices_result;
      // 保存经过质量约束拟合得到的 M 母粒子，用于 X 顶点拟合
      std::vector<RefCountedKinematicParticle> fittedMParticles;

      for (unsigned int i = 0; i < NM && massConstraintSuccess && processCandidate; ++i)
      {
         const VertexFitResult &mCandidate = selectedM[i];
         if (mCandidate.daughterP4.size() != mCandidate.trackIndices.size())
         {
            massConstraintSuccess = false;
            break;
         }
         std::vector<RefCountedKinematicParticle> particles;
         KinematicParticleFactoryFromTransientTrack particleFactory;
         float chi = 0.0, ndf = 0.0;
         // 构造每个子粒子（初始构造仅用于拟合，后续参数均取自拟合结果）
         for (unsigned int j = 0; j < mCandidate.daughterP4.size() && massConstraintSuccess; ++j)
         {
            unsigned int trackIdx = mCandidate.trackIndices[j];
            if (trackIdx >= trackVec.size())
            {
               massConstraintSuccess = false;
               break;
            }
            const reco::Track &trk = trackVec[trackIdx];
            TransientTrack tt(trk, &bFieldHandle);
            double massVal = mCandidate.daughterP4[j].M();
            double massErr = massVal / 1000.0;
            try
            {
               ParticleMass PMassTmp = massVal;
               float PMassErrTmp = massErr;
               RefCountedKinematicParticle particle = particleFactory.particle(tt, PMassTmp, chi, ndf, PMassErrTmp);
               particles.push_back(particle);
            }
            catch (const std::exception &e)
            {
               massConstraintSuccess = false;
               break;
            }
         }
         if (!massConstraintSuccess || particles.size() != mCandidate.daughterP4.size())
         {
            massConstraintSuccess = false;
            break;
         }

         // 对该 M候选施加质量约束：先利用 KinematicParticleVertexFitter 拟合，再用 KinematicParticleFitter 加约束
         RefCountedKinematicTree constrainedTree;
         bool mcProcess = true;
         {
            MassKinematicConstraint massConstraint(MMass[i], MMassErr[i]);
            KinematicParticleVertexFitter massConstraintFitter;
            
            try
            {
               constrainedTree = massConstraintFitter.fit(particles);
            }
            catch (const std::exception &e)
            {
               mcProcess = false;
            }
            if (mcProcess && (!constrainedTree || !constrainedTree->isValid()))
               mcProcess = false;

            KinematicParticleFitter csFitter;
            if (mcProcess)
            {
               try
               {
                  constrainedTree = csFitter.fit(&massConstraint, constrainedTree);
               }
               catch (const std::exception &e)
               {
                  mcProcess = false;
               }
            }
            if (mcProcess && (!constrainedTree || !constrainedTree->isValid()))
               mcProcess = false;

            if (mcProcess)
            {
               constrainedTree->movePointerToTheTop();
            }
         }
         if (!mcProcess)
         {
            massConstraintSuccess = false;
            break;
         }

         // 对 constrainedTree 进行合法性检查
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

         // 提取子粒子信息（均采用拟合后新参数）
         std::vector<TLorentzVector> fittedDaughters;
         constrainedTree->movePointerToTheFirstChild();
         for (unsigned int j = 0; j < mCandidate.daughterP4.size(); ++j)
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
            if (j < mCandidate.daughterP4.size() - 1)
               constrainedTree->movePointerToTheNextChild();
         }
         if (fittedDaughters.empty())
         {
            massConstraintSuccess = false;
            break;
         }
         mDaughter_result.push_back(std::make_pair(fittedDaughters, mCandidate.daughterCharge));
         auto mPar = mMotherParticleMC->currentState().kinematicParameters();
         TVector3 mMom(mPar.momentum().x(), mPar.momentum().y(), mPar.momentum().z());
         double mEnergy = std::sqrt(mMom.Mag2() + mMotherParticleMC->currentState().mass() * mMotherParticleMC->currentState().mass());
         TLorentzVector mP4;
         mP4.SetPxPyPzE(mMom.x(), mMom.y(), mMom.z(), mEnergy);
         mMother_result.push_back(std::make_pair(mP4, mCandidate.motherCharge));
         mMassErr_result.push_back(std::sqrt(mMotherParticleMC->currentState().kinematicParametersError().matrix()(6, 6)));
         mTrackIndices_result.push_back(mCandidate.trackIndices);

         constrainedTree->movePointerToTheTop();
         RefCountedKinematicParticle constrainedM = constrainedTree->currentParticle();
         fittedMParticles.push_back(constrainedM);
      } // end for 每个 M候选

      /*if (massConstraintSuccess && processCandidate)
      {
         cout << "7" << endl;
         // Print out the fitted M particles
         for (unsigned int i = 0; i < NM; i++)
         {
            cout << "M" << i << " mass: " << mMother_result[i].first.M() << " charge: " << mMother_result[i].second << endl;
         }
         // Print daughters of the fitted M particles
         for (unsigned int i = 0; i < NM; i++)
         {
            for (unsigned int j = 0; j < mDaughter_result[i].first.size(); j++)
            {
               cout << "M" << i << " daughter " << j << " mass: " << mDaughter_result[i].first[j].M() << " charge: " << mDaughter_result[i].second[j] << endl;
            }
         }
      }*/

      // 第三阶段：如果所有 M候选均成功，则用新拟合得到的 M母粒子进行 X级顶点拟合
      bool xVertexSuccess = true;
      XFitResult xResult;
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

               // Check the M particle mass from the X fit
               /*xVertexFitTree->movePointerToTheFirstChild();
               RefCountedKinematicParticle mParticle = xVertexFitTree->currentParticle();
               cout << "X fit M mass: " << mParticle->currentState().mass() << endl;
               cout << "M mass:       " << mMother_result[0].first.M() << endl;*/

               xResult.mDaughter = mDaughter_result;
               xResult.mMother = mMother_result;
               xResult.xMotherP4 = xP4;
               xResult.xMotherCharge = XCharge;
               xResult.xMassError = std::sqrt(xMotherParticle->currentState().kinematicParametersError().matrix()(6, 6));
               xResult.normChi2 = xVertex->chiSquared() / (double)xVertex->degreesOfFreedom();
               xResult.mCandidateIndices = candIndices;
               xResult.mTrackIndices = mTrackIndices_result;
               xResult.mMassError = mMassErr_result;
               xVertexSuccess = true;
               xCandidates.push_back(xResult);
            }
         }
      }

      // 生成下一个候选组合（多维笛卡尔组合遍历）
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
   } // end while candidates

   return xCandidates;
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
      for (unsigned int j = 0; j < 4; j++)
      {
         mMass[i][j] = 0;
         mMassError[i][j] = 0;
         mPt[i][j] = 0;
         mEta[i][j] = 0;
         mPhi[i][j] = 0;
         mChi2[i][j] = 0;
         mCharge[i][j] = 0;

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