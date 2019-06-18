import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process("ANA",eras.Run2_2017,eras.run2_nanoAOD_94XMiniAODv2)

#runType="test"
#runType="HPy8"
#runType="H"
runType="A"
#runType="Mix"
#runType="DY"
print "Run type:",runType

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')
process.GlobalTag.globaltag = '102X_mc2017_realistic_v6' #as NanoAOD Fall17

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")


process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) 
)
if runType=="test":
	process.maxEvents.input = 10000

DYfileNames = cms.untracked.vstring( #~250k
    '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/20000/00D13F2E-6F44-E811-923E-001E0BED0560.root',
    '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/20000/00DA5E7F-5044-E811-AE50-001E677925CC.root',
    '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/20000/02890AC5-5D44-E811-9572-E0071B7A9800.root',
)

HTauTauPy8FileNames = cms.untracked.vstring( #~250k
    '/store/mc/RunIIFall17MiniAODv2/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/10000/047009D4-3D21-E911-9727-6C3BE5B541F8.root',
    '/store/mc/RunIIFall17MiniAODv2/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/10000/0604D195-1421-E911-A3BD-6C3BE5B59210.root',
    '/store/mc/RunIIFall17MiniAODv2/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/10000/0A52B0C2-3D21-E911-BC74-D8D385B0EE2E.root',
)

HTauTauFileNames = cms.untracked.vstring( #~200k
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_9.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_8.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_7.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_6.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_5.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_40.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_4.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_39.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_38.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_37.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_36.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_35.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_34.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_33.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_32.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_31.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_30.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_3.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_29.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_28.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_27.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_26.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_25.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_24.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_23.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_22.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_21.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_20.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_2.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_19.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_18.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_17.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_16.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_15.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_14.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_13.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_12.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_11.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_10.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1/180528_091819/0000/RunIIFall17MiniAODv2_1.root',
)

ATauTauFileNames = cms.untracked.vstring(# ~200k
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_9.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_8.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_7.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_6.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_5.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_40.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_4.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_39.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_38.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_37.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_36.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_35.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_34.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_33.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_32.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_31.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_30.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_3.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_29.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_28.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_27.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_26.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_25.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_24.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_23.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_22.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_21.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_20.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_2.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_19.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_18.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_16.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_15.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_14.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_13.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_12.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_11.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_10.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1/180528_091905/0000/RunIIFall17MiniAODv2_1.root',
)

MixCPTauTauFileNames = cms.untracked.vstring( #~200k
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_9.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_8.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_7.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_6.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_5.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_40.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_4.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_39.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_38.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_37.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_36.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_35.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_34.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_33.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_32.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_31.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_30.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_3.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_29.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_28.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_27.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_26.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_25.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_24.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_23.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_22.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_21.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_20.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_2.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_19.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_17.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_16.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_15.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_14.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_13.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_12.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_11.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_10.root',
    '/store/user/bluj/MINIAODSIM/GluGluHToTauTau_M125_13TeV_powheg_pythia8/MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1/180528_082138/0000/RunIIFall17MiniAODv2_1.root',
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/mnt/home/mbluj/work/data/94X/MiniAOD/GluGluHToTauTau_M125_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/047009D4-3D21-E911-9727-6C3BE5B541F8.root'
	)
)


if runType=="DY":
    process.source.fileNames = DYfileNames
elif runType=="HPy8":
    process.source.fileNames = HTauTauPy8FileNames
elif runType=="H":
    process.source.fileNames = HTauTauFileNames
elif runType=="A":
    process.source.fileNames = ATauTauFileNames
elif runType=="Mix":
    process.source.fileNames = MixCPTauTauFileNames
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

process.load("WarsawAnalysis.HTT.miniVtxAna_cfi")
process.load("PhysicsTools.NanoAOD.taus_updatedMVAIds_cff")
process.vtxAna = process.miniVtxAna.clone(taus="slimmedTausUpdated")
process.vtxAnaNoTauTracks = process.vtxAna.clone(
    useTauTracks = cms.untracked.bool(False)
)

process.p = cms.Path(
    process.printTree+
    process.patTauMVAIDsSeq+
    process.vtxAna
    +process.vtxAnaNoTauTracks
)

process.TFileService = cms.Service(
    "TFileService", 
    fileName = cms.string("vtxAna.root") 
)
if runType!="test":
    process.TFileService.fileName = "vtxAna_"+runType+".root"
