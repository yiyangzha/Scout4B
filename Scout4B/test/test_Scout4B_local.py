import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from Configuration.AlCa.GlobalTag import GlobalTag
import time

process = cms.Process('BPHScoutNANO')
# set file name according time
ouput_filename = 'Scout4B_' + time.strftime("%Y%m%d-%H%M%S") + '.root'

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

process.MessageLogger.cerr.FwkSummary.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.threshold = 'ERROR'



process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) # -1 
)

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch/' + '/store/mc/RunIII2024Summer24MiniAOD/BsToJPsiPhi-JPsiToMuMu-PhiToKK_Par-SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/MINIAODSIM/140X_mcRun3_2024_realistic_v26-v2/110000/03bfca1b-ea3b-43f3-a87b-5466d7a48ae2.root',
        #'root://cms-xrd-global.cern.ch/' + '/store/data/Run2024F/ScoutingPFRun3/HLTSCOUT/v1/000/382/299/00000/4f516e00-43f7-4b65-ade1-490a3590faac.root',
        #'root://cms-xrd-global.cern.ch/' + '/store/data/Run2024F/ScoutingPFRun3/HLTSCOUT/v1/000/382/299/00000/5d476985-1a6d-4bd2-8206-f69c068908ec.root',
        #'root://cms-xrd-global.cern.ch/' + '/store/data/Run2024F/ScoutingPFRun3/HLTSCOUT/v1/000/382/299/00000/6b90b87c-76a8-4dcd-91d7-117150371df1.root',
        #'root://cms-xrd-global.cern.ch/' + '/store/data/Run2024F/ScoutingPFRun3/HLTSCOUT/v1/000/382/299/00000/898e8e03-f206-429f-b4f9-79982a891fa3.root',
        #'root://cms-xrd-global.cern.ch/' + '/store/data/Run2024F/ScoutingPFRun3/HLTSCOUT/v1/000/382/299/00000/b8f216c9-24e9-4f22-81bf-2e0e5e705b9a.root',
        #'root://cms-xrd-global.cern.ch/' + '/store/data/Run2024F/ScoutingPFRun3/HLTSCOUT/v1/000/382/299/00000/cc48f5a8-2a33-4067-ad1c-1c800a1af070.root',
        #'root://cms-xrd-global.cern.ch/' + '/store/data/Run2024F/ScoutingPFRun3/HLTSCOUT/v1/000/382/300/00000/07026d37-c60d-4f33-8b6c-71535dbfd4f0.root',
        #'root://cms-xrd-global.cern.ch/' + '/store/data/Run2024F/ScoutingPFRun3/HLTSCOUT/v1/000/382/300/00000/121bc64b-86d9-4385-a93d-1c35378a7059.root',
        #'root://cms-xrd-global.cern.ch/' + '/store/data/Run2024F/ScoutingPFRun3/HLTSCOUT/v1/000/382/300/00000/2bbc6972-6e79-4b87-aebc-8d8944cd8a8a.root',
        )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(ouput_filename),
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

#

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
    'DST_PFScouting_DoubleMuon',               #  1=        1
    'DST_PFScouting_DatasetMuon',              #  2=        2
    'DST_PFScouting_SingleMuon',               #  3=        4
    'DST_PFScouting_ZeroBias'                  #  4=        8
    'HLT_DoubleMu4_3_Bs',                      #  5=        16
    'HLT_DoubleMu4_Jpsi_Displaced',            #  6=        32
    'HLT_DoubleMu4_Jpsi_NoVertexing',          #  7=        64
    'HLT_Mu7p5_L2Mu2_Jpsi',                    #  8=        128
    'HLT_Dimuon0_Jpsi'                         #  9=        256
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

'''
# Bs -> Phi Phi -> K+ K- K+ K-
process.Scout4BVertexFinderBsPP  = cms.EDAnalyzer("Scout4BRecoSecondaryVertexNMCAnalyzer",
    multiM = cms.untracked.bool(True),
    treename = cms.untracked.string("BsPhiPhiKKKK"),
    recoTrackMuon = cms.InputTag("Scout4BConverter", "recoTrackMuons"),
    recoTrack = cms.InputTag("Scout4BConverter", "recoTracks"),
    recoVertex = cms.InputTag("Scout4BConverter", "recoVertexs"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
    FilterNames     = cms.vstring(
    'DST_PFScouting_DoubleMuon',               #  1=        1
    'DST_PFScouting_DatasetMuon',              #  2=        2
    'DST_PFScouting_SingleMuon',               #  3=        4
    'DST_PFScouting_ZeroBias'                  #  4=        8
    'HLT_DoubleMu4_3_Bs',                      #  5=        16
    'HLT_DoubleMu4_Jpsi_Displaced',            #  6=        32
    'HLT_DoubleMu4_Jpsi_NoVertexing',          #  7=        64
    'HLT_Mu7p5_L2Mu2_Jpsi',                    #  8=        128
    'HLT_Dimuon0_Jpsi'                         #  9=        256
    ),
    MMuMass = cms.vdouble(0),
    MMuMassErr = cms.vdouble(0),
    MTrkMass = cms.vdouble(1.01946, 1.01946),                     # phi(1020) 
    MTrkMassErr = cms.vdouble(0.000016, 0.000016),
    MuMass = cms.untracked.double(0.1056583755),         
    MuMassErr = cms.untracked.double(1e-7),
    MuChargeFlat = cms.vint32(0),
    TrkMassFlat = cms.vdouble(0.493677, 0.493677, 0.493677, 0.493677),       # Kaon
    TrkMassErrFlat = cms.vdouble(0.000015, 0.000015, 0.000015, 0.000015),
    TrkChargeFlat = cms.vint32(1, -1, 1, -1),
    NMmu = cms.uint32(0),
    NMtrk = cms.uint32(2),
    NPmu = cms.vuint32(0),
    NPtrk = cms.vuint32(2, 2),
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


# B0 -> J/psi K*0 -> mu+ mu- K+ pi-
process.Scout4BVertexFinderB0JK = cms.EDAnalyzer("Scout4BRecoSecondaryVertexNMCAnalyzer",
    multiM = cms.untracked.bool(True),
    treename = cms.untracked.string("B0JpsiKstarMuMuKpPim"),
    recoTrackMuon = cms.InputTag("Scout4BConverter", "recoTrackMuons"),
    recoTrack = cms.InputTag("Scout4BConverter", "recoTracks"),
    recoVertex = cms.InputTag("Scout4BConverter", "recoVertexs"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
    FilterNames     = cms.vstring(
    'DST_PFScouting_DoubleMuon',               #  1=        1
    'DST_PFScouting_DatasetMuon',              #  2=        2
    'DST_PFScouting_SingleMuon',               #  3=        4
    'DST_PFScouting_ZeroBias'                  #  4=        8
    'HLT_DoubleMu4_3_Bs',                      #  5=        16
    'HLT_DoubleMu4_Jpsi_Displaced',            #  6=        32
    'HLT_DoubleMu4_Jpsi_NoVertexing',          #  7=        64
    'HLT_Mu7p5_L2Mu2_Jpsi',                    #  8=        128
    'HLT_Dimuon0_Jpsi'                         #  9=        256
    ),
    MMuMass = cms.vdouble(3.0969),                       # J/psi
    MMuMassErr = cms.vdouble(0.000006),
    MTrkMass = cms.vdouble(0.89555),                     # K*0(892)
    MTrkMassErr = cms.vdouble(0.0002),
    MuMass = cms.untracked.double(0.1056583755),         # muon
    MuMassErr = cms.untracked.double(1e-7),
    MuChargeFlat = cms.vint32(1, -1),
    TrkMassFlat = cms.vdouble(0.493677, 0.1395704),      # Kaon + Pion
    TrkMassErrFlat = cms.vdouble(0.000015, 1e-7),
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


# B0 -> J/psi K*0 -> mu+ mu- K- pi+
process.Scout4BVertexFinderB0JKB = cms.EDAnalyzer("Scout4BRecoSecondaryVertexNMCAnalyzer",
    multiM = cms.untracked.bool(True),
    treename = cms.untracked.string("B0JpsiKstarMuMuKmPip"),
    recoTrackMuon = cms.InputTag("Scout4BConverter", "recoTrackMuons"),
    recoTrack = cms.InputTag("Scout4BConverter", "recoTracks"),
    recoVertex = cms.InputTag("Scout4BConverter", "recoVertexs"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
    FilterNames     = cms.vstring(
    'DST_PFScouting_DoubleMuon',               #  1=        1
    'DST_PFScouting_DatasetMuon',              #  2=        2
    'DST_PFScouting_SingleMuon',               #  3=        4
    'DST_PFScouting_ZeroBias'                  #  4=        8
    'HLT_DoubleMu4_3_Bs',                      #  5=        16
    'HLT_DoubleMu4_Jpsi_Displaced',            #  6=        32
    'HLT_DoubleMu4_Jpsi_NoVertexing',          #  7=        64
    'HLT_Mu7p5_L2Mu2_Jpsi',                    #  8=        128
    'HLT_Dimuon0_Jpsi'                         #  9=        256
    ),
    MMuMass = cms.vdouble(3.0969),                       # J/psi
    MMuMassErr = cms.vdouble(0.000006),
    MTrkMass = cms.vdouble(0.89555),                     # K*0(892)
    MTrkMassErr = cms.vdouble(0.0002),
    MuMass = cms.untracked.double(0.1056583755),         # muon
    MuMassErr = cms.untracked.double(1e-7),
    MuChargeFlat = cms.vint32(1, -1),
    TrkMassFlat = cms.vdouble(0.493677, 0.1395704),      # Kaon + Pion
    TrkMassErrFlat = cms.vdouble(0.000015, 1e-7),
    TrkChargeFlat = cms.vint32(-1, 1),
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


# D0 -> K+ pi-
process.Scout4BVertexFinderD0 = cms.EDAnalyzer("Scout4BRecoSecondaryVertexNMCAnalyzer",
    multiM = cms.untracked.bool(False),
    treename = cms.untracked.string("D0KpPim"),
    recoTrackMuon = cms.InputTag("Scout4BConverter", "recoTrackMuons"),
    recoTrack = cms.InputTag("Scout4BConverter", "recoTracks"),
    recoVertex = cms.InputTag("Scout4BConverter", "recoVertexs"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
    FilterNames     = cms.vstring(
    'DST_PFScouting_DoubleMuon',               #  1=        1
    'DST_PFScouting_DatasetMuon',              #  2=        2
    'DST_PFScouting_SingleMuon',               #  3=        4
    'DST_PFScouting_ZeroBias'                  #  4=        8
    'HLT_DoubleMu4_3_Bs',                      #  5=        16
    'HLT_DoubleMu4_Jpsi_Displaced',            #  6=        32
    'HLT_DoubleMu4_Jpsi_NoVertexing',          #  7=        64
    'HLT_Mu7p5_L2Mu2_Jpsi',                    #  8=        128
    'HLT_Dimuon0_Jpsi'                         #  9=        256
    ),
    MMuMass = cms.vdouble(0),
    MMuMassErr = cms.vdouble(0),
    MTrkMass = cms.vdouble(1.86484),                      # D0
    MTrkMassErr = cms.vdouble(0.00005),
    MuMass = cms.untracked.double(0.1056583755),
    MuMassErr = cms.untracked.double(1e-7),
    MuChargeFlat = cms.vint32(0),
    TrkMassFlat = cms.vdouble(0.493677, 0.1395704),       # Kaon + Pion
    TrkMassErrFlat = cms.vdouble(0.000015, 1e-7),
    TrkChargeFlat = cms.vint32(1, -1),
    NMmu = cms.uint32(0),
    NMtrk = cms.uint32(1),
    NPmu = cms.vuint32(0),
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


# D0 -> K+ pi-
process.Scout4BVertexFinderD0B = cms.EDAnalyzer("Scout4BRecoSecondaryVertexNMCAnalyzer",
    multiM = cms.untracked.bool(False),
    treename = cms.untracked.string("D0KmPip"),
    recoTrackMuon = cms.InputTag("Scout4BConverter", "recoTrackMuons"),
    recoTrack = cms.InputTag("Scout4BConverter", "recoTracks"),
    recoVertex = cms.InputTag("Scout4BConverter", "recoVertexs"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
    FilterNames     = cms.vstring(
    'DST_PFScouting_DoubleMuon',               #  1=        1
    'DST_PFScouting_DatasetMuon',              #  2=        2
    'DST_PFScouting_SingleMuon',               #  3=        4
    'DST_PFScouting_ZeroBias'                  #  4=        8
    'HLT_DoubleMu4_3_Bs',                      #  5=        16
    'HLT_DoubleMu4_Jpsi_Displaced',            #  6=        32
    'HLT_DoubleMu4_Jpsi_NoVertexing',          #  7=        64
    'HLT_Mu7p5_L2Mu2_Jpsi',                    #  8=        128
    'HLT_Dimuon0_Jpsi'                         #  9=        256
    ),
    MMuMass = cms.vdouble(0),
    MMuMassErr = cms.vdouble(0),
    MTrkMass = cms.vdouble(1.86484),                      # D0
    MTrkMassErr = cms.vdouble(0.00005),
    MuMass = cms.untracked.double(0.1056583755),
    MuMassErr = cms.untracked.double(1e-7),
    MuChargeFlat = cms.vint32(0),
    TrkMassFlat = cms.vdouble(0.493677, 0.1395704),       # Kaon + Pion
    TrkMassErrFlat = cms.vdouble(0.000015, 1e-7),
    TrkChargeFlat = cms.vint32(-1, 1),
    NMmu = cms.uint32(0),
    NMtrk = cms.uint32(1),
    NPmu = cms.vuint32(0),
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
#process.add_(cms.Service("InitRootHandlers", DebugLevel = cms.untracked.int32(1)))
'''
process.convert_step = cms.Path(process.Scout4BConverter)
process.filter_step = cms.Path(process.Scout4BScoutTrkFilter)
process.vertexFinderBsJP = cms.Path(process.Scout4BVertexFinderBsJP)
'''
process.vertexFinderBsPP = cms.Path(process.Scout4BVertexFinderBsPP)
process.vertexFinderB0JK = cms.Path(process.Scout4BVertexFinderB0JK)
process.vertexFinderB0JKB = cms.Path(process.Scout4BVertexFinderB0JKB)
process.vertexFinderD0 = cms.Path(process.Scout4BVertexFinderD0)
process.vertexFinderD0B = cms.Path(process.Scout4BVertexFinderD0B)
'''

process.schedule = cms.Schedule(
    process.convert_step,
    process.filter_step,
    process.vertexFinderBsJP,
    #process.vertexFinderBsPP,
    #process.vertexFinderB0JK,
    #process.vertexFinderB0JKB,
    #process.vertexFinderD0,
    #process.vertexFinderD0B
)