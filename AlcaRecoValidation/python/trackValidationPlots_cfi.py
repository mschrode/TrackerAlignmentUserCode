import FWCore.ParameterSet.Config as cms

trackValidationPlots = cms.EDAnalyzer(
    'TrackValidationPlots',
    TrackCollection  = cms.InputTag(""),
    VertexCollection = cms.InputTag(""),
)
