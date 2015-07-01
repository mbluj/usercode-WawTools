import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) 
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	"root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/SUSYGluGluToHToTauTau_M-120_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/2C871A0B-3303-E511-B8D0-0025B3E05D74.root"
	)
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(1),# 1->10
  printVertex = cms.untracked.bool(False),
  printOnlyHardInteraction = cms.untracked.bool(False), # Print only status=3 particles. This will not work for Pythia8, which does not have any such particles.
  #src = cms.InputTag("genParticles")
  src = cms.InputTag("prunedGenParticles")
)

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    # Include a PSet for each module label that needs a
    # random engine.  The name is the module label.
    # You must supply a seed or seeds.
    # Optionally an engine type can be specified
    vtxAna = cms.PSet(
        initialSeed = cms.untracked.uint32(82)
    ),
    # This is optional.  If you want the service to save the state
    # of all engines to a separate text file which is overwritten before
    # modules begin processing on each event. The text file is only
    # needed for one type of replay. The filename can be anything
    # you want but the replay process will need to reference it.
    #saveFileName = cms.untracked.string('RandomEngineStates.txt')
)

process.vtxAna = cms.EDAnalyzer(
    "MiniAODVertexAnalyzer",
    src = cms.InputTag("packedPFCandidates"),
    genSrc = cms.InputTag("prunedGenParticles"),
    #pruned = cms.InputTag("genParticles"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    vertexScores = cms.InputTag("offlineSlimmedPrimaryVertices"),

    verbose = cms.untracked.bool(False),
)

process.p = cms.Path(
    #process.printTree
    process.vtxAna
)

process.TFileService = cms.Service(
    "TFileService", 
    fileName = cms.string("genTauTauAna.root") 
)
