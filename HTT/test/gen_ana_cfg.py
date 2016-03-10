import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

# switches
runOnSignal = False
#verbose = True
verbose = False
#redoGenParticles = True #True for GEN-SIM with CMSSW<74X?
redoGenParticles = False

if runOnSignal:
    print "Runing on signal, toy-Taus build only for jets matched to genTaus"
else:
    print "Runing on background, toy-Taus build for all jets"

if redoGenParticles:
    print "Redoing genParticles"

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)

process.maxEvents = cms.untracked.PSet( 
    #input = cms.untracked.int32(100) 
    input = cms.untracked.int32(-1) 
)
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff') 
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('74X_mcRun2_asymptotic_v2')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	#h0"root://cms-xrd-global.cern.ch//store/user/bluj/SUSYGluGluH0ToTauTau_M-125_TuneCUETP8M1_13TeV_pythia8_GEN_RunIIWinter15GS-MCRUN2_71_V1-v2/SUSYGluGluH0ToTauTau_M-125_TuneCUETP8M1_13TeV_pythia8_GEN_RunIIWinter15GS-MCRUN2_71_V1-v2/17aeb1344d1bc358be14f25dc5a8c422/USR-RunIIWinter15GS-GEN_99_1_BY0.root"
        #A0
        "file:/mnt/home/mbluj/work/data/SUSYGluGluA0ToTauTau_M-125_TuneCUETP8M1_13TeV_pythia8_GEN_RunIIWinter15GS-MCRUN2_71_V1-v1/USR-RunIIWinter15GS-GEN_9_1_MDz.root"
	#h0"file:/mnt/home/mbluj/work/CMSSW/CMSSW_7_1_13_prod/src/Prod/Test/test/USR-RunIIWinter15GS-GEN.root" #h0 v2
	#A0"root://cms-xrd-global.cern.ch//store/user/bluj/SUSYGluGluA0ToTauTau_M-125_TuneCUETP8M1_13TeV_pythia8_GEN_RunIIWinter15GS-MCRUN2_71_V1-v1/SUSYGluGluA0ToTauTau_M-125_TuneCUETP8M1_13TeV_pythia8_GEN_RunIIWinter15GS-MCRUN2_71_V1-v1/ca5c4cbeea1a7b632759b4a8c2b92091/USR-RunIIWinter15GS-GEN_99_1_Srj.root"
	#Z0"root://cms-xrd-global.cern.ch//store/user/bluj/DYToTauTau_M-50_TuneCUETP8M1_13TeV_pythia8_GEN_RunIIWinter15GS-MCRUN2_71_V1-v1/DYToTauTau_M-50_TuneCUETP8M1_13TeV_pythia8_GEN_RunIIWinter15GS-MCRUN2_71_V1-v1/76d179419fbff7169973d9fdfd1ce825/USR-RunIIWinter15GS-GEN_99_2_5lw.root"
        #WJets_HT100-200"root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/20574402-2801-E511-BBFA-0CC47A4D9A1E.root"
	)
)

#---- genParticles
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")

process.dumpPATTauSequence = cms.Sequence()
#--------------------------------------------------------------------------------
# produce collection of genTauJets
process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
process.dumpPATTauSequence += process.tauGenJets
process.load("PhysicsTools.JetMCAlgos.TauGenJetsDecayModeSelectorAllHadrons_cfi")
process.tauGenJetsSelectorAllHadrons.select = cms.vstring(
    'oneProng0Pi0', 
    'oneProng1Pi0', 
    'oneProng2Pi0', 
    'oneProngOther',
    'threeProng0Pi0', 
    'threeProng1Pi0', 
    'threeProngOther', 
    'rare'
)
process.tauGenJetsSelectorAllHadrons.filter = cms.bool(runOnSignal)
process.dumpPATTauSequence += process.tauGenJetsSelectorAllHadrons

process.selectedTauGenJets = cms.EDFilter("GenJetSelector",
    src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
    #MBcut = cms.string("pt > 20. & abs(eta) < 2.3"),
    cut = cms.string("pt > 0. & abs(eta) < 5"), #MB dummy cuts (in practice PFtaus are reconstructed for Pt>~10 and |eta|<2.4
    filter = cms.bool(runOnSignal)
)
process.dumpPATTauSequence += process.selectedTauGenJets

#--------------------------------------------------------------------------------
# produce collection of "toy" PFCandidates 
process.load("RecoTauTag.GenTauAnalysisTools.toyPFCandidates_cfi")

process.dumpPATTauSequence += process.toyPFCandidates

#--------------------------------------------------------------------------------
# run jet reconstruction on "toy" PFCandidates
from RecoJets.Configuration.RecoPFJets_cff import ak4PFJets

process.ak4ToyPFJets = ak4PFJets.clone(
   src = cms.InputTag('toyPFCandidates', 'pfCandidates')
)
process.dumpPATTauSequence += process.ak4ToyPFJets

#--------------------------------------------------------------------------------
# run tau reconstruction on "toy" PFCandidates
if runOnSignal:
    process.genTauMatchedPFJets = cms.EDFilter(
        "PFJetAntiOverlapSelector",
        src = cms.InputTag('ak4ToyPFJets'),
        srcNotToBeFiltered = cms.VInputTag('selectedTauGenJets'),
        dRmin = cms.double(0.3),
        invert = cms.bool(True),
        filter = cms.bool(False)  
    )
    process.dumpPATTauSequence += process.genTauMatchedPFJets

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag 
massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('particleFlow'), cms.InputTag('toyPFCandidates', 'pfCandidates'))
if runOnSignal:
    massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('ak4PFJets'), cms.InputTag('genTauMatchedPFJets'))
else:
    massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('ak4PFJets'), cms.InputTag('ak4ToyPFJets'))
massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('generalTracks'), cms.InputTag('toyPFCandidates', 'tracks'))
massSearchReplaceAnyInputTag(process.PFTau, cms.InputTag('offlinePrimaryVertices'), cms.InputTag('toyPFCandidates', 'vertices'))

process.dumpPATTauSequence += process.PFTau
#MB Remove anti-e MVA discriminats which need RECO info (e.g. gsfElectrons)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVA5rawElectronRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVA5VLooseElectronRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVA5LooseElectronRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVA5MediumElectronRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVA5TightElectronRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVA5VTightElectronRejection)
#MB  Remove anti-mu MVA discriminats which need RECO info (e.g. reco muons)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByLooseMuonRejection2)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMediumMuonRejection2)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByTightMuonRejection2)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByLooseMuonRejection3)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByTightMuonRejection3)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVArawMuonRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVALooseMuonRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVAMediumMuonRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVATightMuonRejection)

#--------------------------------------------------------------------------------
# produce PAT-tuple
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# switch to HPS PFTaus (and disable all "cleaning" cuts)
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)
#MB Remove anti-e MVA discriminats which need RECO info (e.g. gsfElectrons)
oldTauIds =  process.patTaus.tauIDSources.parameters_()
cleanedTauIDs = {}
for key,value in oldTauIds.iteritems():
    if ( ((key.find("MVA5")<0 or key.find("againstElectron")<0) and (key.find("againstMuon")<0)) ):
        cleanedTauIDs[key]=value
    else:
        print "Removing tau-ID :",key,"=",value
#MB add old anti-mu discriminats (do they work for toy taus?)
cleanedTauIDs["againstMuonLoose"] = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection")
cleanedTauIDs["againstMuonMedium"] = cms.InputTag("hpsPFTauDiscriminationByMediumMuonRejection")
cleanedTauIDs["againstMuonTight"] = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection")
process.patTaus.tauIDSources = cms.PSet(**cleanedTauIDs)

#FIXME: add some tau selection?
# tau-Id makes sense only in tracker coverage: |eta|<~2.4 and Pt>~10GeV
# use only taus with reconstructed DM
# keep dummy Pt
process.selectedPatTaus.cut = "tauID(\"decayModeFindingNewDMs\") && abs(eta)<2.4 && pt>0" 
process.cleanPatTaus.preselection = cms.string('')
process.cleanPatTaus.checkOverlaps = cms.PSet()
process.cleanPatTaus.finalCut = cms.string("")

process.makePatTaus.remove(process.patPFCandidateIsoDepositSelection)
process.makePatTaus.remove(process.patPFTauIsolation)
process.dumpPATTauSequence += process.makePatTaus
process.dumpPATTauSequence += process.selectedPatTaus
#MBprocess.dumpPATTauSequence += process.cleanPatTaus

#--------------------------------------------------------------------------------
#to reporoduce pruned and packed genParticles as in MiniAOD
process.load("PhysicsTools.PatAlgos.slimming.genParticles_cff")
#process.packedGenParticles.inputVertices # Producer of packedGenParticles hacked to not call vertices (which are anyhow unused)
process.packedGenParticles.inputVertices = cms.InputTag('toyPFCandidates', 'vertices') #use toy gen vertices

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
    #pruned = cms.InputTag("prunedGenParticles"),
    pruned = cms.InputTag("genParticles"),
    toyTaus = cms.InputTag("selectedPatTaus"),
    isSignal = cms.untracked.bool(True),
    verbose = cms.untracked.bool(verbose),
)

process.p = cms.Path()
if redoGenParticles:
    process.p += process.genParticles #reproduce genParticles in current release for correct statusFlags 
process.p += (
    process.dumpPATTauSequence +
    process.prunedGenParticlesWithStatusOne+process.prunedGenParticles +
    process.packedGenParticles
    #+ process.printTree
    + process.genAna
)

process.TFileService = cms.Service(
    "TFileService", 
    fileName = cms.string("genTauTauAna.root") 
)
