#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/FinalTreeBuilder.h"
#include "RecoVertex/VertexTools/interface/SequentialVertexSmoother.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanSmoothedVertexChi2Estimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanTrackToTrackCovCalculator.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "RecoVertex/LinearizationPointFinders/interface/DefaultLinearizationPointFinder.h"
#include "RecoVertex/VertexTools/interface/SequentialVertexFitter.h"
#include "DataFormats/CLHEP/interface/Migration.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

KinematicParticleVertexFitter::KinematicParticleVertexFitter()
{
  edm::ParameterSet pSet = defaultParameters();
  setup(pSet);
}
// Changed by Suleyman
// KinematicParticleVertexFitter::KinematicParticleVertexFitter(const edm::ParameterSet pSet)
KinematicParticleVertexFitter::KinematicParticleVertexFitter(const edm::ParameterSet &pSet)
// suleyman end
{
  setup(pSet);
}

void
// suleyman
// KinematicParticleVertexFitter::setup(const edm::ParameterSet pSet)
KinematicParticleVertexFitter::setup(const edm::ParameterSet &pSet)
// suleyman
{

  pointFinder = new DefaultLinearizationPointFinder();
  vFactory = new VertexTrackFactory<6>();

  KalmanVertexTrackUpdator<6> vtu;
  KalmanSmoothedVertexChi2Estimator<6> vse;
  KalmanTrackToTrackCovCalculator<6> covCalc;
  SequentialVertexSmoother<6> smoother(vtu, vse, covCalc);
  fitter = new SequentialVertexFitter<6>(pSet, *pointFinder, KalmanVertexUpdator<6>(), smoother, ParticleKinematicLinearizedTrackStateFactory());
}

KinematicParticleVertexFitter::~KinematicParticleVertexFitter()
{
  delete vFactory;
  delete pointFinder;
  delete fitter;
}

edm::ParameterSet KinematicParticleVertexFitter::defaultParameters() const
{
  edm::ParameterSet pSet;
  pSet.addParameter<double>("maxDistance", 0.01);
  pSet.addParameter<int>("maxNbrOfIterations", 100); // 10
  return pSet;
}
// suleyman
// RefCountedKinematicTree KinematicParticleVertexFitter::fit(std::vector<RefCountedKinematicParticle> particles) const
RefCountedKinematicTree KinematicParticleVertexFitter::fit(const std::vector<RefCountedKinematicParticle> &particles) const
// suleyma
{
  typedef ReferenceCountingPointer<VertexTrack<6>> RefCountedVertexTrack;
  // sorting the input
  if (particles.size() < 2)
    throw VertexException("KinematicParticleVertexFitter::input states are less than 2");
  InputSort iSort;
  std::pair<std::vector<RefCountedKinematicParticle>, std::vector<FreeTrajectoryState>> input = iSort.sort(particles);
  std::vector<RefCountedKinematicParticle> &newPart = input.first;
  std::vector<FreeTrajectoryState> &freeStates = input.second;

  GlobalPoint linPoint = pointFinder->getLinearizationPoint(freeStates);

  // cout<<"Linearization point found"<<endl;

  // making initial veretx seed with lin point as position and a fake error
  AlgebraicSymMatrix33 we;
  we(0, 0) = we(1, 1) = we(2, 2) = 10000.;
  GlobalError error(we);
  VertexState state(linPoint, error);

  // vector of Vertex Tracks to fit
  std::vector<RefCountedVertexTrack> ttf;
  // for(std::vector<RefCountedKinematicParticle>::const_iterator i = newPart.begin();i != newPart.end();i++)
  // {ttf.push_back(vFactory->vertexTrack((*i)->particleLinearizedTrackState(linPoint),state,1.));}

  for (std::vector<RefCountedKinematicParticle>::const_iterator i = newPart.begin(); i != newPart.end(); i++)
  {
    TrackKinematicStatePropagator propagator_;
    if (!(*i)->currentState().isValid() || !propagator_.propagateToTheTransversePCA((*i)->currentState(), linPoint).isValid())
    {
      std::cout << "Here's the bad state." << std::endl;
      return ReferenceCountingPointer<KinematicTree>(new KinematicTree()); // return invalid vertex
    }
    ttf.push_back(vFactory->vertexTrack((*i)->particleLinearizedTrackState(linPoint), state, 1.));
  }

  CachingVertex<6> vtx = fitter->vertex(ttf);
  if (!vtx.isValid())
  {
    LogDebug("RecoVertex/KinematicParticleVertexFitter") << "Fitted position is invalid. Returned Tree is invalid\n";
    return ReferenceCountingPointer<KinematicTree>(new KinematicTree()); // return invalid vertex
  }
  FinalTreeBuilder tBuilder;
  return tBuilder.buildTree(vtx, newPart);
}
