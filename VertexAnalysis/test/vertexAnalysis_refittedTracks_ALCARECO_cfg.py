import FWCore.ParameterSet.Config as cms

process = cms.Process("VertexAnalysis")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.threshold             = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 500
#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cout = cms.untracked.PSet(
#    INFO = cms.untracked.PSet(
#        reportEvery = cms.untracked.int32(1000)
#    )
#)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)


## --- Conditions ------------------------------------------------------
from CondCore.DBCommon.CondDBSetup_cfi import *
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string("FT53_V21A_AN6::All")

# load mis-aligned geometry on top of GT geometry
#import CalibTracker.Configuration.Common.PoolDBESSource_cfi
#process.customTrackerAlignment = CalibTracker.Configuration.Common.PoolDBESSource_cfi.poolDBESSource.clone(
#    connect = cms.string("sqlite_file:/afs/desy.de/user/m/matsch/CMSSW_5_3_23_patch1/src/TrackerAlignmentUserCode/VertexAnalysis/data/geometry_MisalignmentScenario200Mu_fromFT53_V21A_AN6_fromRandomTool.db"),
#    toGet = cms.VPSet(
#        cms.PSet(
#            record = cms.string('TrackerAlignmentRcd'),
#            tag = cms.string('Alignments')
#        )
#    )
#)
process.es_prefer_trackerAlignment = cms.ESPrefer("PoolDBESSource","customTrackerAlignment")




## --- Input file ------------------------------------------------------
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/02B803D1-B770-E211-8287-00266CFFC940.root',
        '/store/data/Run2012A/SingleMu/ALCARECO/TkAlMuonIsolated-22Jan2013-v1/20000/08039FC6-6C70-E211-897C-1CC1DE1CEDB2.root'
    )
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5000)
)
#import FWCore.PythonUtilities.LumiList as LumiList
#goodLumiSecs = LumiList.LumiList(filename = '../data/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt').getVLuminosityBlockRange()
#process.source.lumisToProcess = goodLumiSecs



## --- Output file -----------------------------------------------------
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("VertexAnalysis.root")
)


## --- Tracks and PV re-fitting ----------------------------------------

process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.refittedTracks = process.TrackRefitter.clone(
    src ="ALCARECOTkAlMuonIsolated",
)

process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.refittedOfflinePrimaryVertices = process.offlinePrimaryVertices.clone(
    TrackLabel = cms.InputTag("refittedTracks")
)


## --- Vertex analysis -------------------------------------------------
from TrackerAlignmentUserCode.VertexAnalysis.vertexanalysis_cfi import vertexanalysis
process.VertexAnalysis = vertexanalysis.clone(
    TreeName         = cms.string("VertexAnalysis"),
    MaxNTracks       = cms.int32(200),
    TrackCollection  = cms.InputTag("refittedTracks"),
    MaxNVertices     = cms.int32(50),
    VertexCollection = cms.InputTag("refittedOfflinePrimaryVertices"),
    GenParticlesCollection = cms.InputTag("genParticles"),
)


## --- Run the anlaysis ------------------------------------------------

process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.p = cms.Path(
    process.offlineBeamSpot *
    process.refittedTracks *
    process.refittedOfflinePrimaryVertices *
#    process.dump *
    process.VertexAnalysis 
)
