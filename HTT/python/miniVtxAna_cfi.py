import FWCore.ParameterSet.Config as cms

miniVtxAna = cms.EDAnalyzer(
    "MiniAODVertexAnalyzer",
    src = cms.InputTag("packedPFCandidates"),
    lostSrc = cms.InputTag("lostTracks"),
    lostEleSrc = cms.InputTag("lostTracks:eleTracks"),
    genSrc = cms.InputTag("prunedGenParticles"),
    #pruned = cms.InputTag("genParticles"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    vertexScores = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    taus = cms.InputTag("slimmedTaus"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    useBeamSpot = cms.bool(True),
    useLostCands = cms.bool(True),
    useTauTracks = cms.untracked.bool(True),
    verbose = cms.untracked.bool(False),
)
