import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from Configuration.AlCa.GlobalTag import GlobalTag
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process('BPHScoutNANO')
options = VarParsing('analysis')
options.parseArguments()
options.outputFile = 'Scout4B.root'

assert len(options.inputFiles) == 1, 'Only run interactively with file len=1'
print("Running", options.inputFiles)
print("Output file", options.outputFile)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_PromptAnalysis_v1', '') #2022 ABCDE PromptReco
#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_v15', '') #2022 ABCDE ReReco
#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_PromptAnalysis_v2', '') #2022 FG PromptReco
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_PromptAnalysis_v1', '') #2023 CD ReReco

process.MessageLogger.cerr.FwkSummary.reportEvery = 10000
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
#process.MessageLogger.cerr.threshold = 'ERROR'

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) # -1 
)

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(options.inputFiles),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile),
)

process.options = cms.untracked.PSet(
    TryToContinue = cms.untracked.vstring("ProductNotFound","TooManyProducts","TooFewProducts"),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring("ProductNotFound","TooManyProducts","TooFewProducts"),
    Rethrow = cms.untracked.vstring("ProductNotFound","TooManyProducts","TooFewProducts"),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(2),
    numberOfThreads = cms.untracked.uint32(2),
    printDependencies = cms.untracked.bool(False),
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False),
)

process.Scout4BConverter = cms.EDProducer("Scout4BScoutToRecoProducer",
    scoutingMuon = cms.InputTag("hltScoutingMuonPackerVtx"),
    scoutingMuonNoVtx = cms.InputTag("hltScoutingMuonPackerNoVtx"),
    scoutingTrack = cms.InputTag("hltScoutingTrackPacker"),
    scoutingPrimaryVertex = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx")
)

process.Scout4BScoutMuFilter = cms.EDFilter("TrackCountFilter",
    src       = cms.InputTag("Scout4BConverter", "recoTrackMuons"),
    minNumber = cms.uint32(2)
)


process.Scout4BScoutTrkFilter = cms.EDFilter("TrackCountFilter",
    src       = cms.InputTag("Scout4BConverter", "recoTracks"),
    minNumber = cms.uint32(2)
)

# Bs -> J/psi phi -> mu+ mu- K+ K-
process.Scout4BVertexFinderBsJP = cms.EDAnalyzer("Scout4BRecoSecondaryVertexAnalyzer",
    multiM = cms.untracked.bool(True),
    treename = cms.untracked.string("BsJpsiPhiMuMuKK"),
    recoTrackMuon = cms.InputTag("Scout4BConverter", "recoTrackMuons"),
    recoTrack = cms.InputTag("Scout4BConverter", "recoTracks"),
    recoVertex = cms.InputTag("Scout4BConverter", "recoVertexs"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
    FilterNames     = cms.vstring(
    'DST_PFScouting_DoubleMuon',               #  0=        1
    'DST_PFScouting_DatasetMuon',              #  1=        2
    'DST_PFScouting_SingleMuon',               #  2=        4
    'DST_PFScouting_ZeroBias'                  #  3=        8
    'HLT_DoubleMu4_3_Bs',                      #  4=        16
    'HLT_DoubleMu4_Jpsi_Displaced',            #  5=        32
    'HLT_DoubleMu4_Jpsi_NoVertexing',          #  6=        64
    'HLT_Mu7p5_L2Mu2_Jpsi',                    #  7=        128
    'HLT_Dimuon0_Jpsi'                         #  8=        256
    ),
    MMuMass = cms.vdouble(3.0969),                       # J/psi
    MMuMassErr = cms.vdouble(0.000006),
    MTrkMass = cms.vdouble(1.01946),                     # phi(1020) 
    MTrkMassErr = cms.vdouble(0.000016),
    MuMass = cms.untracked.double(0.1056583755),         # muon
    MuMassErr = cms.untracked.double(1e-7),
    MuChargeFlat = cms.vint32(1, -1),
    TrkMassFlat = cms.vdouble(0.493677, 0.493677),       # Kaon
    TrkMassErrFlat = cms.vdouble(0.000015, 0.000015),
    TrkChargeFlat = cms.vint32(1, -1),
    NMmu = cms.uint32(1),
    NMtrk = cms.uint32(1),
    NPmu = cms.vuint32(2),
    NPtrk = cms.vuint32(2),
    XMass = cms.untracked.double(0.0),
    XMassWin = cms.untracked.double(1000.0),
    MMassWin = cms.untracked.double(0.2),
    XCharge = cms.int32(0),
    MChargeMu = cms.vint32(0),
    MChargetrk = cms.vint32(0),
    PIso = cms.untracked.double(0.01),
    doIso = cms.untracked.bool(True),
    doTrigger = cms.untracked.bool(False),
    vProbMin = cms.untracked.double(1e-3),
    maxLoop = cms.untracked.uint64(100000000),
    MMassMin = cms.untracked.double(1e-2),
    XMassMin = cms.untracked.double(1e-2),
)

process.convert_step = cms.Path(process.Scout4BConverter)
process.muonfilter_step = cms.Path(process.Scout4BScoutMuFilter)
process.trkfilter_step = cms.Path(process.Scout4BScoutTrkFilter)
process.bemaspot_step = cms.Path(process.offlineBeamSpot)
process.vertexFinderBsJP = cms.Path(process.Scout4BVertexFinderBsJP)

process.schedule = cms.Schedule(
    process.convert_step,
    process.muonfilter_step,
    process.trkfilter_step,
    process.bemaspot_step,
    process.vertexFinderBsJP
)