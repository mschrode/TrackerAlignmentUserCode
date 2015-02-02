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




## --- Input file ------------------------------------------------------
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Commissioning12/Cosmics/ALCARECO/TkAlCosmics0T-13Jul2012-v1/00000/0274A00C-76CE-E111-9C8F-003048D15E02.root'
    )
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
#import FWCore.PythonUtilities.LumiList as LumiList
#goodLumiSecs = LumiList.LumiList(filename = '../data/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt').getVLuminosityBlockRange()
#process.source.lumisToProcess = goodLumiSecs



## --- Output file -----------------------------------------------------
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("AlcaRecoValidationPlots.root")
)


## --- Track validation plots ------------------------------------------
from TrackerAlignmentUserCode.AlcaRecoValidation.trackValidationPlots_cfi import trackValidationPlots
process.CosmicsTrackPlots = trackValidationPlots.clone(
    TrackCollection  = cms.InputTag("ALCARECOTkAlCosmicsCTF0T"),
    VertexCollection = cms.InputTag(""),
)


## --- Run the anlaysis ------------------------------------------------

process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.p = cms.Path(
#    process.dump *
    process.CosmicsTrackPlots
)
