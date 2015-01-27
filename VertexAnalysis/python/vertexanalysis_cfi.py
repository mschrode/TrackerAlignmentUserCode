import FWCore.ParameterSet.Config as cms

vertexanalysis = cms.EDAnalyzer(
    'VertexAnalysis',
    TreeName         = cms.string("VertexAnalysis"),
    MaxNTracks       = cms.int32(200),
    TrackCollection  = cms.InputTag(""),
    MaxNVertices     = cms.int32(50),
    VertexCollection = cms.InputTag(""),
    GenParticlesCollection = cms.InputTag("genParticles"),
)
