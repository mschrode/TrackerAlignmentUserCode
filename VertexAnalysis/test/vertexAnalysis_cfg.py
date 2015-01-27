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
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string("FT53_V21A_AN6::All")
#process.GlobalTag.globaltag = cms.string("FT_R_53_V18::All")
#process.GlobalTag.globaltag = cms.string("START53_V23::All")
#process.GlobalTag.globaltag = cms.string("DESIGN53_V18::All")



## --- Input file ------------------------------------------------------
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/00A602AA-BF8B-E211-B161-90B11C186EC4.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/00C65ECC-078C-E211-9D30-80000048FE80.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/02B028E3-318B-E211-B4D1-0002C94CD13C.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/02C7F6BA-BF8B-E211-901C-90B11C18C363.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/043B6D4B-CF8B-E211-A4C9-90B11C18C535.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/047D4D37-078C-E211-8581-0002C94CD310.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/0681F199-188C-E211-A121-80000048FE80.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/0A178F52-378C-E211-A934-842B2B298D0A.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/0A214E78-2B8B-E211-882E-001EC9AAD233.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/0A4FF812-2B8C-E211-B159-782BCB6E134C.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/0AA7065C-2B8C-E211-ACEB-001AA00D8CED.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/0AC295C2-2F8B-E211-B2D4-80000048FE80.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/0AC8FFA2-2A8B-E211-92FF-001EC9AA9FFE.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/0C7F7F9D-298B-E211-ACB1-001EC9AAD459.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/0CB9DA0C-2F8B-E211-9BC7-0002C90A363C.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/0CDAB8A4-BF8B-E211-ACF8-001AA009A50D.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/0EE2E699-BF8B-E211-9EE5-90B11C18C879.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/0EECF9FE-BF8B-E211-B239-90B11C186D00.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/1053A00D-2F8B-E211-9EE5-80000048FE80.root',
        '/store/data/Run2012C/LP_MinBias1/RECO/22Jan2013-v1/10000/105B85F3-BF8B-E211-A758-90B11C1862C2.root',


        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/001C3948-A682-E211-A511-00304867906C.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/003DD113-8682-E211-98AE-0025905938B4.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/0094BB86-A982-E211-B52A-00261894395A.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/00AD589F-9E82-E211-94D4-003048FFD730.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/00BAB721-B082-E211-8768-0026189438B3.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/00DC9D73-A682-E211-AE61-002618943910.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/02119125-A282-E211-B952-00248C65A3EC.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/047AD091-A982-E211-945E-0030486790B8.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/04BB0994-9682-E211-9CCF-00261894390A.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/04C897C5-B082-E211-A438-002590596490.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/04DE01DA-AD82-E211-9242-002618943857.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/060C8DDD-AD82-E211-A5C6-003048F9EB46.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/062D5493-9C82-E211-A378-002618943935.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/0652F652-AC82-E211-8414-0025905964C0.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/06C3158F-AD82-E211-B275-002590596498.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/0825FADE-9282-E211-98D2-003048678AE2.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/083EA0EA-B382-E211-B591-002354EF3BE2.root',
        #'/store/data/Run2012A/MET/RECO/22Jan2013-v1/20000/086CC53E-9E82-E211-8CFC-00261894397F.root',
    )
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50000)
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

#process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.refittedTracks = process.TrackRefitter.clone(
    src ="generalTracks",
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
    #TrackCollection  = cms.InputTag("refittedTracks"),
    TrackCollection  = cms.InputTag("generalTracks"),
    MaxNVertices     = cms.int32(50),
    #VertexCollection = cms.InputTag("refittedOfflinePrimaryVertices"),
    VertexCollection = cms.InputTag("offlinePrimaryVertices"),
    GenParticlesCollection = cms.InputTag("genParticles"),
)


## --- Run the anlaysis ------------------------------------------------

process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.p = cms.Path(
#    process.refittedTracks *
#    process.refittedOfflinePrimaryVertices *
#    process.dump *
    process.VertexAnalysis 
)
