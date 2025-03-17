/*
 * Developed by Yiyang Zhao for Run-3 B scouting
 * 2025-03
 * Scout4BScoutToRecoProducer unpacks Run3ScoutingMuon formats to reco formats
 * It creates reco::Muon and reco::Track objects from Run3ScoutingMuon and Run3ScoutingTrack objects
 * It also creates reco::Track objects from Run3ScoutingTrack objects
 */

#include <memory>

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/EDPutToken.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

// #include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
// #include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
// #include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"

// #include "DataFormats/JetReco/interface/PFJet.h"
// #include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

// #include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
// #include "fastjet/contrib/SoftKiller.hh"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

class Scout4BScoutToRecoProducer : public edm::stream::EDProducer<>
{
public:
   explicit Scout4BScoutToRecoProducer(edm::ParameterSet const &iConfig);
   ~Scout4BScoutToRecoProducer() override;

   static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
   void beginStream(edm::StreamID) override {}
   void produce(edm::Event &iEvent, edm::EventSetup const &setup) override;
   void endStream() override {}

   reco::Track createTrack(Run3ScoutingMuon const &scoutingMuon);
   reco::Track createTrack(Run3ScoutingTrack const &scoutingTrack);
   reco::Vertex createVertex(Run3ScoutingVertex const &scoutingVertex);
   reco::Muon createMuon(Run3ScoutingMuon const &scoutingMuon);

   void createMuons(edm::Handle<std::vector<Run3ScoutingMuon>> scoutingmuonHandle,
                    std::unique_ptr<reco::MuonCollection> &scoutingmuons);
   void createTracks(edm::Handle<std::vector<Run3ScoutingMuon>> scoutingmuonHandle,
                     std::unique_ptr<reco::TrackCollection> &scoutingtracks);
   void createTracks(edm::Handle<std::vector<Run3ScoutingTrack>> scoutingtrackHandle,
                     std::unique_ptr<reco::TrackCollection> &scoutingtracks);
   void createVertexs(edm::Handle<std::vector<Run3ScoutingVertex>> scoutingvertexHandle,
                      std::unique_ptr<reco::VertexCollection> &scoutingvertexs);

   void clearVars();

   bool checkOverlap(Run3ScoutingMuon const &scoutingMuon, Run3ScoutingMuon const &scoutingMuonNoVtx);

private:
   const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> input_scoutingmuon_token_;
   const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> input_scoutingmuonNoVtx_token_;
   const edm::EDGetTokenT<std::vector<Run3ScoutingTrack>> input_scoutingtrack_token_;
   const edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> scoutingPrimaryVertex_collection_token_;
};

//
// constructors and destructor
//
Scout4BScoutToRecoProducer::Scout4BScoutToRecoProducer(edm::ParameterSet const &iConfig)
    : input_scoutingmuon_token_(consumes(iConfig.getParameter<edm::InputTag>("scoutingMuon"))),
      input_scoutingmuonNoVtx_token_(consumes(iConfig.getParameter<edm::InputTag>("scoutingMuonNoVtx"))),
      input_scoutingtrack_token_(consumes(iConfig.getParameter<edm::InputTag>("scoutingTrack"))),
      scoutingPrimaryVertex_collection_token_(consumes(iConfig.getParameter<edm::InputTag>("scoutingPrimaryVertex")))
{
   // register products
   produces<reco::MuonCollection>("recoMuons");
   produces<reco::TrackCollection>("recoTrackMuons");
   produces<reco::TrackCollection>("recoTracks");
   produces<reco::VertexCollection>("recoVertexs");
}

Scout4BScoutToRecoProducer::~Scout4BScoutToRecoProducer() = default;

// ------------ method called to produce the data  ------------
void Scout4BScoutToRecoProducer::produce(edm::Event &iEvent, edm::EventSetup const &setup)
{
   using namespace edm;

   Handle<std::vector<Run3ScoutingMuon>> scoutingmuonHandle;
   Handle<std::vector<Run3ScoutingMuon>> scoutingmuonNoVtxHandle;
   Handle<std::vector<Run3ScoutingTrack>> scoutingtrackHandle;
   Handle<std::vector<Run3ScoutingVertex>> scoutingvertexHandle;

   iEvent.getByToken(input_scoutingmuon_token_, scoutingmuonHandle);
   iEvent.getByToken(input_scoutingmuonNoVtx_token_, scoutingmuonNoVtxHandle);
   iEvent.getByToken(input_scoutingtrack_token_, scoutingtrackHandle);
   iEvent.getByToken(scoutingPrimaryVertex_collection_token_, scoutingvertexHandle);

   if (!scoutingmuonHandle.isValid() || !scoutingtrackHandle.isValid() || !scoutingvertexHandle.isValid())
   {
      return;
   }
   // else std::cout << "ScoutingMuon is valid" << std::endl;

   auto pfcands = std::make_unique<reco::MuonCollection>();     // Store reco::Muons for output
   auto tkcands = std::make_unique<reco::TrackCollection>();    // Store reco::Tracks from Muons for output
   auto tkallcands = std::make_unique<reco::TrackCollection>(); // Store reco::Tracks for output
   auto vtxcands = std::make_unique<reco::VertexCollection>();  // Store reco::Vertex for output

   // Create reco::Muons from scoutingmuon and scoutingmuonNoVtx
   createMuons(scoutingmuonHandle, pfcands);
   // Create reco::Muons from scoutingmuonNoVtx
   if (scoutingmuonNoVtxHandle.isValid())
   {
      for (unsigned int icand = 0; icand < scoutingmuonNoVtxHandle->size(); ++icand)
      {
         auto &scoutingmuonNoVtx = (*scoutingmuonNoVtxHandle)[icand];
         bool found = false;
         for (unsigned int jcand = 0; jcand < scoutingmuonHandle->size(); ++jcand)
         {
            auto &scoutingmuon = (*scoutingmuonHandle)[jcand];
            if (checkOverlap(scoutingmuon, scoutingmuonNoVtx))
            {
               found = true;
               break;
            }
         }
         if (!found)
         {
            auto pfcand = createMuon(scoutingmuonNoVtx);
            if (pfcand.energy() != 0)
               pfcands->push_back(pfcand);
         }
      }
   }

   // Create reco::Tracks from scoutingmuon and scoutingmuonNoVtx
   createTracks(scoutingmuonHandle, tkcands);
   // Create reco::Tracks from scoutingmuonNoVtx
   if (scoutingmuonNoVtxHandle.isValid())
   {
      for (unsigned int icand = 0; icand < scoutingmuonNoVtxHandle->size(); ++icand)
      {
         auto &scoutingmuonNoVtx = (*scoutingmuonNoVtxHandle)[icand];
         bool found = false;
         for (unsigned int jcand = 0; jcand < scoutingmuonHandle->size(); ++jcand)
         {
            auto &scoutingmuon = (*scoutingmuonHandle)[jcand];
            if (checkOverlap(scoutingmuon, scoutingmuonNoVtx))
            {
               found = true;
               break;
            }
         }
         if (!found)
         {
            auto tkcand = createTrack(scoutingmuonNoVtx);
            if (tkcand.pt() != 0)
               tkcands->push_back(tkcand);
         }
      }
   }

   // Create reco::Tracks from scoutingtrack
   createTracks(scoutingtrackHandle, tkallcands);

   // Create reco::Vertex from scoutingvertex
   createVertexs(scoutingvertexHandle, vtxcands);

   // std::cout << "Number of reco::Muons: " << pfcands->size() << std::endl;
   // std::cout << "Number of reco::Tracks from Muons: " << tkcands->size() << std::endl;
   // std::cout << "Number of reco::Tracks: " << tkallcands->size() << std::endl;

   edm::OrphanHandle<reco::MuonCollection> oh = iEvent.put(std::move(pfcands), "recoMuons");
   edm::OrphanHandle<reco::TrackCollection> oh2 = iEvent.put(std::move(tkcands), "recoTrackMuons");
   edm::OrphanHandle<reco::TrackCollection> oh3 = iEvent.put(std::move(tkallcands), "recoTracks");
   edm::OrphanHandle<reco::VertexCollection> oh4 = iEvent.put(std::move(vtxcands), "recoVertexs");

   clearVars();
}

void Scout4BScoutToRecoProducer::clearVars()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Scout4BScoutToRecoProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
   edm::ParameterSetDescription desc;
   desc.add<edm::InputTag>("scoutingMuon", edm::InputTag("hltScoutingMuonPackerVtx"));
   desc.add<edm::InputTag>("scoutingMuonNoVtx", edm::InputTag("hltScoutingMuonPackerNoVtx"));
   desc.add<edm::InputTag>("scoutingTrack", edm::InputTag("hltScoutingTrackPacker"));
   desc.add<edm::InputTag>("scoutingPrimaryVertex", edm::InputTag("hltScoutingPrimaryVertexPacker"));
   
   descriptions.add("Scout4BScoutToRecoProducer", desc);
}

reco::Track Scout4BScoutToRecoProducer::createTrack(Run3ScoutingTrack const &scoutingTrack)
{
   float chi2 = scoutingTrack.tk_chi2();
   float ndof = scoutingTrack.tk_ndof();
   reco::TrackBase::Point referencePoint(scoutingTrack.tk_vx(), scoutingTrack.tk_vy(), scoutingTrack.tk_vz());

   float px = scoutingTrack.tk_pt() * cos(scoutingTrack.tk_phi());
   float py = scoutingTrack.tk_pt() * sin(scoutingTrack.tk_phi());
   float pz = scoutingTrack.tk_pt() * sinh(scoutingTrack.tk_eta());
   reco::TrackBase::Vector momentum(px, py, pz);

   int charge = scoutingTrack.tk_charge();

   std::vector<float> cov_vec(15);                                                 // 5*(5+1)/2 = 15
   cov_vec[0] = scoutingTrack.tk_qoverp_Error() * scoutingTrack.tk_qoverp_Error(); // cov(0, 0)
   cov_vec[1] = scoutingTrack.tk_qoverp_lambda_cov();                              // cov(0, 1)
   cov_vec[3] = scoutingTrack.tk_qoverp_phi_cov();                                 // cov(0, 2)
   cov_vec[6] = scoutingTrack.tk_qoverp_dxy_cov();                                 // cov(0, 3)
   cov_vec[10] = scoutingTrack.tk_qoverp_dsz_cov();                                // cov(0, 4)
   cov_vec[2] = scoutingTrack.tk_lambda_Error() * scoutingTrack.tk_lambda_Error(); // cov(1, 1)
   cov_vec[4] = scoutingTrack.tk_lambda_phi_cov();                                 // cov(1, 2)
   cov_vec[7] = scoutingTrack.tk_lambda_dxy_cov();                                 // cov(1, 3)
   cov_vec[11] = scoutingTrack.tk_lambda_dsz_cov();                                // cov(1, 4)
   cov_vec[5] = scoutingTrack.tk_phi_Error() * scoutingTrack.tk_phi_Error();       // cov(2, 2)
   cov_vec[8] = scoutingTrack.tk_phi_dxy_cov();                                    // cov(2, 3)
   cov_vec[12] = scoutingTrack.tk_phi_dsz_cov();                                   // cov(2, 4)
   cov_vec[9] = scoutingTrack.tk_dxy_Error() * scoutingTrack.tk_dxy_Error();       // cov(3, 3)
   cov_vec[13] = scoutingTrack.tk_dxy_dsz_cov();                                   // cov(3, 4)
   cov_vec[14] = scoutingTrack.tk_dsz_Error() * scoutingTrack.tk_dsz_Error();      // cov(4, 4)
   reco::TrackBase::CovarianceMatrix cov(cov_vec.begin(), cov_vec.end());

   reco::TrackBase::TrackAlgorithm algo(reco::TrackBase::undefAlgorithm); // undefined
   reco::TrackBase::TrackQuality quality(reco::TrackBase::confirmed);     // confirmed

   // the rests are default: t0 = 0, beta = 0, covt0t0 = -1, covbetabeta = -1

   reco::Track recoTrack(chi2, ndof, referencePoint, momentum, charge, cov, algo, quality);

   return recoTrack;
}

reco::Track Scout4BScoutToRecoProducer::createTrack(Run3ScoutingMuon const &scoutingMuon)
{
   float chi2 = scoutingMuon.trk_chi2();
   float ndof = scoutingMuon.trk_ndof();
   reco::TrackBase::Point referencePoint(scoutingMuon.trk_vx(), scoutingMuon.trk_vy(), scoutingMuon.trk_vz());

   float px = scoutingMuon.trk_pt() * cos(scoutingMuon.trk_phi());
   float py = scoutingMuon.trk_pt() * sin(scoutingMuon.trk_phi());
   float pz = scoutingMuon.trk_pt() * sinh(scoutingMuon.trk_eta());
   reco::TrackBase::Vector momentum(px, py, pz);

   int charge = scoutingMuon.charge();

   std::vector<float> cov_vec(15);                                               // 5*(5+1)/2 = 15
   cov_vec[0] = scoutingMuon.trk_qoverpError() * scoutingMuon.trk_qoverpError(); // cov(0, 0)
   cov_vec[1] = scoutingMuon.trk_qoverp_lambda_cov();                            // cov(0, 1)
   cov_vec[3] = scoutingMuon.trk_qoverp_phi_cov();                               // cov(0, 2)
   cov_vec[6] = scoutingMuon.trk_qoverp_dxy_cov();                               // cov(0, 3)
   cov_vec[10] = scoutingMuon.trk_qoverp_dsz_cov();                              // cov(0, 4)
   cov_vec[2] = scoutingMuon.trk_lambdaError() * scoutingMuon.trk_lambdaError(); // cov(1, 1)
   cov_vec[4] = scoutingMuon.trk_lambda_phi_cov();                               // cov(1, 2)
   cov_vec[7] = scoutingMuon.trk_lambda_dxy_cov();                               // cov(1, 3)
   cov_vec[11] = scoutingMuon.trk_lambda_dsz_cov();                              // cov(1, 4)
   cov_vec[5] = scoutingMuon.trk_phiError() * scoutingMuon.trk_phiError();       // cov(2, 2)
   cov_vec[8] = scoutingMuon.trk_phi_dxy_cov();                                  // cov(2, 3)
   cov_vec[12] = scoutingMuon.trk_phi_dsz_cov();                                 // cov(2, 4)
   cov_vec[9] = scoutingMuon.trk_dxyError() * scoutingMuon.trk_dxyError();       // cov(3, 3)
   cov_vec[13] = scoutingMuon.trk_dxy_dsz_cov();                                 // cov(3, 4)
   cov_vec[14] = scoutingMuon.trk_dszError() * scoutingMuon.trk_dszError();      // cov(4, 4)
   reco::TrackBase::CovarianceMatrix cov(cov_vec.begin(), cov_vec.end());

   reco::TrackBase::TrackAlgorithm algo(reco::TrackBase::undefAlgorithm); // undefined
   reco::TrackBase::TrackQuality quality(reco::TrackBase::confirmed);     // confirmed

   // the rests are default: t0 = 0, beta = 0, covt0t0 = -1, covbetabeta = -1

   reco::Track recoTrack(chi2, ndof, referencePoint, momentum, charge, cov, algo, quality);

   return recoTrack;
}

reco::Muon Scout4BScoutToRecoProducer::createMuon(Run3ScoutingMuon const &scoutingMuon)
{
   auto m = 0.1056583755;
   auto q = scoutingMuon.charge();

   float px = scoutingMuon.pt() * cos(scoutingMuon.phi());
   float py = scoutingMuon.pt() * sin(scoutingMuon.phi());
   float pz = scoutingMuon.pt() * sinh(scoutingMuon.eta());
   float p = scoutingMuon.pt() * cosh(scoutingMuon.eta());
   float energy = std::sqrt(p * p + m * m);
   reco::Particle::LorentzVector p4(px, py, pz, energy);

   static const reco::Muon dummy;
   auto recomuon = reco::Muon(q, p4, createTrack(scoutingMuon).vertex());

   return recomuon;
}

reco::Vertex Scout4BScoutToRecoProducer::createVertex(Run3ScoutingVertex const &scoutingVertex)
{
   // fill point coordinate
   reco::Vertex::Point point(scoutingVertex.x(), scoutingVertex.y(), scoutingVertex.z());

   // fill error
   std::vector<float> error_vec(6);
   error_vec[0] = scoutingVertex.xError() * scoutingVertex.xError(); // cov(0, 0)
   error_vec[1] = 0;                                                 // cov(0, 1)
   error_vec[2] = 0;                                                 // cov(0, 2)
   error_vec[3] = scoutingVertex.yError() * scoutingVertex.yError(); // cov(1, 1)
   error_vec[4] = 0;                                                 // cov(1, 2)
   error_vec[5] = scoutingVertex.zError() * scoutingVertex.zError(); // cov(2, 2)

   // off-diagonal errors are added in the begining of 2024
   // see https://github.com/cms-sw/cmssw/pull/43758
   try
   {
      error_vec[1] = scoutingVertex.xyCov();
      error_vec[2] = scoutingVertex.xzCov();
      error_vec[4] = scoutingVertex.yzCov();
   }
   catch (...)
   { // do nothing
   }

   reco::Vertex::Error error(error_vec.begin(), error_vec.end());

   return scoutingVertex.isValidVtx() ? reco::Vertex(point, error, scoutingVertex.chi2(), scoutingVertex.ndof(), scoutingVertex.tracksSize()) : reco::Vertex(point, error);
}

void Scout4BScoutToRecoProducer::createMuons(
    edm::Handle<std::vector<Run3ScoutingMuon>> scoutingmuonHandle,
    std::unique_ptr<reco::MuonCollection> &scoutingmuons)
{
   for (unsigned int icand = 0; icand < scoutingmuonHandle->size(); ++icand)
   {
      auto &scoutingmuon = (*scoutingmuonHandle)[icand];

      auto pfcand = createMuon(scoutingmuon);
      if (pfcand.energy() != 0)
         scoutingmuons->push_back(pfcand);
   }
}

void Scout4BScoutToRecoProducer::createTracks(
    edm::Handle<std::vector<Run3ScoutingMuon>> scoutingmuonHandle,
    std::unique_ptr<reco::TrackCollection> &scoutingtracks)
{
   for (unsigned int icand = 0; icand < scoutingmuonHandle->size(); ++icand)
   {
      auto &scoutingmuon = (*scoutingmuonHandle)[icand];

      auto tkcand = createTrack(scoutingmuon);
      if (tkcand.pt() != 0)
         scoutingtracks->push_back(tkcand);
   }
}

void Scout4BScoutToRecoProducer::createTracks(
    edm::Handle<std::vector<Run3ScoutingTrack>> scoutingtrackHandle,
    std::unique_ptr<reco::TrackCollection> &scoutingtracks)
{
   for (unsigned int icand = 0; icand < scoutingtrackHandle->size(); ++icand)
   {
      auto &scoutingmuon = (*scoutingtrackHandle)[icand];

      auto tkcand = createTrack(scoutingmuon);
      if (tkcand.pt() != 0)
         scoutingtracks->push_back(tkcand);
   }
}

void Scout4BScoutToRecoProducer::createVertexs(
    edm::Handle<std::vector<Run3ScoutingVertex>> scoutingvertexHandle,
    std::unique_ptr<reco::VertexCollection> &scoutingvertexs)
{
   for (unsigned int icand = 0; icand < scoutingvertexHandle->size(); ++icand)
   {
      auto &scoutingvertex = (*scoutingvertexHandle)[icand];

      auto vtxcand = createVertex(scoutingvertex);
      scoutingvertexs->push_back(vtxcand);
   }
}

bool Scout4BScoutToRecoProducer::checkOverlap(Run3ScoutingMuon const &scoutingMuon, Run3ScoutingMuon const &scoutingMuonNoVtx)
{
   if (scoutingMuon.charge() != scoutingMuonNoVtx.charge())
      return false;
   // if (abs(scoutingMuon.pt() - scoutingMuonNoVtx.pt()) > 1e-3) return false;
   if (abs(scoutingMuon.eta() - scoutingMuonNoVtx.eta()) > 2e-2)
      return false;
   if (abs(scoutingMuon.phi() - scoutingMuonNoVtx.phi()) > 2e-2)
      return false;
   return true;
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Scout4BScoutToRecoProducer);