from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'ana' #MB: job depended
config.General.workArea = 'crab3_Ana' #MB: job depended
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'pset_Ana.py' #MB: job depended

config.Data.inputDataset = '' #MB: job depended
#config.Data.inputDBS = 'global' #MB: job depended
config.Data.inputDBS = 'phys03' #MB: job depended
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 100000
config.Data.totalUnits = -1
#config.Data.outLFNDirBase = '/store/user/%s/SIM' % (getUsernameFromSiteDB())
config.Data.outLFNDirBase = '/store/user/bluj/HTTCP_VtxNTuple'
config.Data.publication = False
config.Data.outputDatasetTag = '' #MB: job depended, set if publication=True

config.JobType.numCores = 1

#config.Site.whitelist = ["T2_PL_Swierk"]
config.Site.storageSite = 'T2_PL_Swierk'
