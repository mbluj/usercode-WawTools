import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)


process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_PostLS1_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
#process.GlobalTag.globaltag = 'MCRUN2_71_V1' #as MiniAOD

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")


process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) 
)


DYfileNames = cms.untracked.vstring("file:/home/akalinow/scratch/CMS/HiggsCP/Data/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/009D49A5-7314-E511-84EF-0025905A605E.root",
                                    "file:/home/akalinow/scratch/CMS/HiggsCP/Data/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/00C0BECF-6F14-E511-96F8-0025904B739A.root",
                                    "file:/home/akalinow/scratch/CMS/HiggsCP/Data/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/0260F225-7614-E511-A79F-00A0D1EE8EB4.root",
                                    "file:/home/akalinow/scratch/CMS/HiggsCP/Data/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/02B810EA-7214-E511-BDAB-0025905964C2.root")

HTauTauFileNames = cms.untracked.vstring("file:/home/akalinow/scratch/CMS/HiggsCP/Data/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/242A73B4-0A2F-E511-A2A0-00259073E474.root",
                                         "file:/home/akalinow/scratch/CMS/HiggsCP/Data/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/48DF1EE2-1B2F-E511-9EA7-00259073E4C2.root",
                                         "file:/home/akalinow/scratch/CMS/HiggsCP/Data/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/4AF51D45-6A2F-E511-A56A-002618943868.root",
                                         "file:/home/akalinow/scratch/CMS/HiggsCP/Data/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/7ADC24DF-1B2F-E511-B538-00259073E4F6.root")

ATauTauFileNames = cms.untracked.vstring("file:/home/akalinow/scratch/CMS/HiggsCP/Data/SUSYGluGluToHToTauTau_M-120_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/2C871A0B-3303-E511-B8D0-0025B3E05D74.root")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #"root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/SUSYGluGluToHToTauTau_M-120_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/2C871A0B-3303-E511-B8D0-0025B3E05D74.root"
	)
)

#process.source.fileNames = DYfileNames
process.source.fileNames = HTauTauFileNames
#process.source.fileNames.append(HTauTauFileNames[0])
#process.source.fileNames.append(ATauTauFileNames[0])

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
    lostSrc = cms.InputTag("lostTracks"),
    genSrc = cms.InputTag("prunedGenParticles"),
    #pruned = cms.InputTag("genParticles"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    vertexScores = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    taus = cms.InputTag("slimmedTaus"),
    useBeamSpot = cms.bool(True),
    useLostCands = cms.bool(False),
    useTauTracks = cms.untracked.bool(False),
    verbose = cms.untracked.bool(False),
)

process.p = cms.Path(
    process.printTree+
    process.vtxAna
)

process.TFileService = cms.Service(
    "TFileService", 
    fileName = cms.string("vtxAna.root") 
)
