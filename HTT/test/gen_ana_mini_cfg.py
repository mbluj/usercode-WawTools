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
	#"root://cms-xrd-global.cern.ch//store/user/bluj/SUSYGluGluH0ToTauTau_M-125_TuneCUETP8M1_13TeV_pythia8_GEN_RunIIWinter15GS-MCRUN2_71_V1-v2/SUSYGluGluH0ToTauTau_M-125_TuneCUETP8M1_13TeV_pythia8_GEN_RunIIWinter15GS-MCRUN2_71_V1-v2/17aeb1344d1bc358be14f25dc5a8c422/USR-RunIIWinter15GS-GEN_99_1_BY0.root"
	#h0"file:/mnt/home/mbluj/work/CMSSW/CMSSW_7_1_13_prod/src/Prod/Test/test/USR-RunIIWinter15GS-GEN.root" #h0 v2
	#A0"root://cms-xrd-global.cern.ch//store/user/bluj/SUSYGluGluA0ToTauTau_M-125_TuneCUETP8M1_13TeV_pythia8_GEN_RunIIWinter15GS-MCRUN2_71_V1-v1/SUSYGluGluA0ToTauTau_M-125_TuneCUETP8M1_13TeV_pythia8_GEN_RunIIWinter15GS-MCRUN2_71_V1-v1/ca5c4cbeea1a7b632759b4a8c2b92091/USR-RunIIWinter15GS-GEN_99_1_Srj.root"
	#Z0"root://cms-xrd-global.cern.ch//store/user/bluj/DYToTauTau_M-50_TuneCUETP8M1_13TeV_pythia8_GEN_RunIIWinter15GS-MCRUN2_71_V1-v1/DYToTauTau_M-50_TuneCUETP8M1_13TeV_pythia8_GEN_RunIIWinter15GS-MCRUN2_71_V1-v1/76d179419fbff7169973d9fdfd1ce825/USR-RunIIWinter15GS-GEN_99_2_5lw.root"
	"root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/SUSYGluGluToHToTauTau_M-120_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/2C871A0B-3303-E511-B8D0-0025B3E05D74.root"
	)
)

#to reporoduce pruned and packed genParticles as in MiniAOD
#process.load("PhysicsTools.PatAlgos.slimming.genParticles_cff")
#process.packedGenParticles.inputVertices # Producer of packedGenParticles hacked to not call vertices (which are anyhow unused)

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
    genAna = cms.PSet(
        initialSeed = cms.untracked.uint32(81)
    ),
    # This is optional.  If you want the service to save the state
    # of all engines to a separate text file which is overwritten before
    # modules begin processing on each event. The text file is only
    # needed for one type of replay. The filename can be anything
    # you want but the replay process will need to reference it.
    #saveFileName = cms.untracked.string('RandomEngineStates.txt')
)

process.genAna = cms.EDAnalyzer(
    "MiniAODGenTauTauAnalyzer",
    packed = cms.InputTag("packedGenParticles"),
    pruned = cms.InputTag("prunedGenParticles"),
    #pruned = cms.InputTag("genParticles"),
    verbose = cms.untracked.bool(False),
)

process.p = cms.Path(
    #process.prunedGenParticlesWithStatusOne+process.prunedGenParticles +
    #process.packedGenParticles
    #+ process.printTree
    #+ process.genAna
    process.genAna
)

process.TFileService = cms.Service(
    "TFileService", 
    fileName = cms.string("genTauTauAna.root") 
)
