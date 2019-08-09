#!/usr/bin/env python

import os, re
import commands
import math
import urllib

from crab3_Ana import *

dataset_name = {
    'default-scalar': ['Scalar','/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM','global'],
    'susy-pseudoscalar120': ['Pseudoscalar120','/SUSYGluGluToHToTauTau_M-120_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM','global'],
    'susy-pseudoscalar130': ['Pseudoscalar130','/SUSYGluGluToHToTauTau_M-130_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM','global'],
    'scalar-ext' : ['ScalarTauola','/GluGluHToTauTau_M125_13TeV_powheg_pythia8/bluj-MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-ext_v1-18783c0a07109245951450a1a4f55409/USER','phys03'],
    'scalar' : ['ScalarTauola','/GluGluHToTauTau_M125_13TeV_powheg_pythia8/bluj-MiniAODv2-HIG-RunIIFall17MiniAODv2-ScalarTauola-PU2017-v1-18783c0a07109245951450a1a4f55409/USER','phys03'],
    'pseudoscalar-ext' : ['PseudoscalarTauola','/GluGluHToTauTau_M125_13TeV_powheg_pythia8/bluj-MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-ext_v1-18783c0a07109245951450a1a4f55409/USER','phys03'],
    'pseudoscalar' : ['PseudoscalarTauola','/GluGluHToTauTau_M125_13TeV_powheg_pythia8/bluj-MiniAODv2-HIG-RunIIFall17MiniAODv2-PseudoscalarTauola-PU2017-v1-18783c0a07109245951450a1a4f55409/USER','phys03'],
    'maxmix-ext' : ['MaxMixCPTauola','/GluGluHToTauTau_M125_13TeV_powheg_pythia8/bluj-MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-ext_v1-18783c0a07109245951450a1a4f55409/USER','phys03'],
    'maxmix' : ['MaxMixCPTauola','/GluGluHToTauTau_M125_13TeV_powheg_pythia8/bluj-MiniAODv2-HIG-RunIIFall17MiniAODv2-MaxMixCPTauola-PU2017-v1-18783c0a07109245951450a1a4f55409/USER','phys03'],
    'nospin' : ['NoSpin','/GluGluHToPseudoscalarTauTau_M125_13TeV_powheg_pythia8_2017-GEN_TEST07Jan19/adow-GluGluToHToTauTauNoSpin_M125_13TeV_pythia8_2017-MINIAOD-5f646ecd4e1c7a39ab0ed099ff55ceb9/USER','phys03'],
    'nospin-filter' : ['NoSpinFilter','/GluGluHToTauTau_M125_13TeV_powheg_pythia8_nospinner-filter-v2/dwinterb-GluGluHToTauTau_M125_13TeV_powheg_pythia8_nospinner-filter-v2-miniAOD-5f646ecd4e1c7a39ab0ed099ff55ceb9/USER','phys03'], #filter eff 0.2537
}
def prepareCrabCfg(dataset,
                   eventsPerJob=1000000,
                   totalUnits=-1,
                   storage_element='T2_PL_Swierk',
                   publish_data_suffix='test'):

    workdir = 'crab3_VtxAna'
    config.Data.outputDatasetTag = 'VtxAna'
    #if dataset != 'default':
    config.Data.outputDatasetTag += '-'+dataset_name[dataset][0]   
    config.General.requestName = 'VtxAna_'+dataset
    if publish_data_suffix != '':
        workdir += '_'+publish_data_suffix
        config.Data.outputDatasetTag += '-'+publish_data_suffix
        config.General.requestName += '_'+publish_data_suffix
    #config.Data.outLFNDirBase = '/store/user/bluj/HTTCP_VtxNTuple/'+config.General.requestName #not needed?
    config.General.workArea = workdir
    config.Data.inputDataset = dataset_name[dataset][1]
    config.Data.inputDBS = dataset_name[dataset][2]

    ##Modify CRAB3 configuration
    config.JobType.psetName = '../vtx_ana_mini_cfg.py'
    config.Data.unitsPerJob = eventsPerJob
    config.Data.totalUnits = totalUnits
    #config.Data.lumiMask = 'HIG-RunIIFall17GS-'+cfg_name[dataset]+'-ext.json'

    #config.JobType.maxJobRuntimeMin = 2700 #2630 # 2630 min = ~44hrs
    config.JobType.maxMemoryMB = 2500 # 2000 -> 2500 which is granted at T2s

    #config.Data.ignoreLocality = True
    #config.Site.whitelist = ["T2_PL_Swierk"]
    #config.Site.whitelist = ["T2_CH_CERN","T2_DE_DESY","T2_PL_Swierk"]
    #config.Site.whitelist = ["T2_CH_CERN","T2_DE_DESY"]

    config.Site.storageSite = storage_element

    out = open('crabAnaTmp.py','w')
    out.write(config.pythonise_())
    out.close()
    os.system("crab submit -c crabAnaTmp.py")
    os.system("rm crabAnaTmp.py crabAnaTmp.pyc")

#################

eventsPerJob = 100000
totalUnits = -1

for dataset in ['default-scalar',
                'susy-pseudoscalar120',
                'susy-pseudoscalar130',
                #'scalar',
                #'scalar-ext',
                #'pseudoscalar',
                #'pseudoscalar-ext',
                #'maxmix',
                #'maxmix-ext',
                'nospin',
                'nospin-filter',
                ]:
    prepareCrabCfg(
        dataset=dataset,
        eventsPerJob=eventsPerJob,
        totalUnits=totalUnits,
        storage_element='T2_PL_Swierk',
        publish_data_suffix='2019-v5')
