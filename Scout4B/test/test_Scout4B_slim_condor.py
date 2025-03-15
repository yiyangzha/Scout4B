import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process('BScoutNANO')
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
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
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
    numberOfStreams = cms.untracked.uint32(4),
    numberOfThreads = cms.untracked.uint32(4),
    printDependencies = cms.untracked.bool(False),
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False),
)

process.Scout4BConverter = cms.EDProducer("Scout4BScoutToRecoProducer",
    scoutingMuon = cms.InputTag("hltScoutingMuonPackerVtx"),
    scoutingMuonNoVtx = cms.InputTag("hltScoutingMuonPackerNoVtx"),
    scoutingTrack = cms.InputTag("hltScoutingTrackPacker"),
)

'''
process.Scout4BScoutMuFilter = cms.EDFilter("TrackCountFilter",
    src       = cms.InputTag("Scout4BConverter", "recoTrackMuons"),
    minNumber = cms.uint32(2),
)
'''

process.Scout4BScoutTrkFilter = cms.EDFilter("TrackCountFilter",
    src       = cms.InputTag("Scout4BConverter", "recoTracks"),
    minNumber = cms.uint32(2),
)

# Bs -> J/psi phi -> mu+ mu- K+ K-
process.Scout4BVertexFinderBsJP = cms.EDAnalyzer("Scout4BRecoSecondaryVertexAnalyzer",
    multiM = cms.untracked.bool(True),
    treename = cms.untracked.string("BsJpsiPhiMuMuKK"),
    recoTrackMuon = cms.InputTag("Scout4BConverter", "recoTrackMuons"),
    recoTrack = cms.InputTag("Scout4BConverter", "recoTracks"),
    TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
    FilterNames     = cms.vstring(
    'HLT_DoubleMu2_Jpsi_LowPt',                #  1=        1
    'HLT_Dimuon0_Jpsi3p5_Muon2',               #  2=        2
    'HLT_DoubleMu4_JpsiTrkTrk_Displaced',      #  3=        4
    'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05'    #  4=        8
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
    MMassWin = cms.untracked.double(0.5),
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


# Bs -> Phi Phi -> K+ K- K+ K-
process.Scout4BVertexFinderBsPP  = cms.EDAnalyzer("Scout4BRecoSecondaryVertexAnalyzer",
    multiM = cms.untracked.bool(True),
    treename = cms.untracked.string("BsPhiPhiKKKK"),
    recoTrackMuon = cms.InputTag("Scout4BConverter", "recoTrackMuons"),
    recoTrack = cms.InputTag("Scout4BConverter", "recoTracks"),
    TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
    FilterNames     = cms.vstring(
    'HLT_DoubleMu2_Jpsi_LowPt',                #  1=        1
    'HLT_Dimuon0_Jpsi3p5_Muon2',               #  2=        2
    'HLT_DoubleMu4_JpsiTrkTrk_Displaced',      #  3=        4
    'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05'    #  4=        8
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
    MMassWin = cms.untracked.double(0.3),
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
process.Scout4BVertexFinderB0JK = cms.EDAnalyzer("Scout4BRecoSecondaryVertexAnalyzer",
    multiM = cms.untracked.bool(True),
    treename = cms.untracked.string("B0JpsiKstarMuMuKpPim"),
    recoTrackMuon = cms.InputTag("Scout4BConverter", "recoTrackMuons"),
    recoTrack = cms.InputTag("Scout4BConverter", "recoTracks"),
    TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
    FilterNames     = cms.vstring(
    'HLT_DoubleMu2_Jpsi_LowPt',                #  1=        1
    'HLT_Dimuon0_Jpsi3p5_Muon2',               #  2=        2
    'HLT_DoubleMu4_JpsiTrkTrk_Displaced',      #  3=        4
    'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05'    #  4=        8
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
    MMassWin = cms.untracked.double(0.5),
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
process.Scout4BVertexFinderB0JKB = cms.EDAnalyzer("Scout4BRecoSecondaryVertexAnalyzer",
    multiM = cms.untracked.bool(True),
    treename = cms.untracked.string("B0JpsiKstarMuMuKmPip"),
    recoTrackMuon = cms.InputTag("Scout4BConverter", "recoTrackMuons"),
    recoTrack = cms.InputTag("Scout4BConverter", "recoTracks"),
    TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
    FilterNames     = cms.vstring(
    'HLT_DoubleMu2_Jpsi_LowPt',                #  1=        1
    'HLT_Dimuon0_Jpsi3p5_Muon2',               #  2=        2
    'HLT_DoubleMu4_JpsiTrkTrk_Displaced',      #  3=        4
    'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05'    #  4=        8
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
    MMassWin = cms.untracked.double(0.5),
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
process.filter_step = cms.Path(process.Scout4BScoutTrkFilter)
process.vertexFinderBsJP = cms.Path(process.Scout4BVertexFinderBsJP)
process.vertexFinderBsPP = cms.Path(process.Scout4BVertexFinderBsPP)
process.vertexFinderB0JK = cms.Path(process.Scout4BVertexFinderB0JK)
process.vertexFinderB0JKB = cms.Path(process.Scout4BVertexFinderB0JKB)
process.vertexFinderD0 = cms.Path(process.Scout4BVertexFinderD0)
process.vertexFinderD0B = cms.Path(process.Scout4BVertexFinderD0B)

process.schedule = cms.Schedule(
    process.convert_step,
    process.filter_step,
    process.vertexFinderBsJP,
    process.vertexFinderBsPP,
    process.vertexFinderB0JK,
    process.vertexFinderB0JKB
)