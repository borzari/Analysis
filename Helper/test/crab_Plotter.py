#!/usr/bin/env python3

import os
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = ''
config.General.workArea = 'crab'
config.General.transferOutputs = True
config.General.transferLogs = True
# config.General.requestName = '4zToElecProbeTrk' # this is the name of the crab folder; useful to keep track of what is happening // for Pveto
config.General.requestName = 'plotterSingMuonCRNoJetSelSkim' # this is the name of the crab folder; useful to keep track of what is happening // for Poffline, Ptrigger and Nctrl
# config.General.requestName = '3singElecCRNoJetMETTrigSelSkim' # this is the name of the crab folder; useful to keep track of what is happening // for Poffline, Ptrigger and Nctrl

config.JobType.pluginName = 'Analysis'
# config.JobType.psetName = 'zToLepProbeTrk_cfg.py' # for Pveto
config.JobType.psetName = 'plotterSingLepCRNoJetSelSkim_cfg.py' # for Poffline, Ptrigger and Nctrl
config.JobType.allowUndistributedCMSSW = True
# List of files that are in the repos, but are used for the selections
config.JobType.inputFiles = [os.environ['CMSSW_BASE'] + '/src/Analysis/Helper/data/electronFiducialMap_2018_data.root', os.environ['CMSSW_BASE'] + '/src/Analysis/Helper/data/muonFiducialMap_2018_data.root']

config.Data.inputDataset = '/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/borzari-2singMuonCRNoJetSelSkim-a7b3b93af29e296702ddceda65a665b0/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 20 # this is the amount of files processed per output file

config.Data.publication = True
# config.Data.outputDatasetTag = '4zToElecProbeTrk' # this is just an example; it will be part of the name of the output dataset // for Pveto
config.Data.outputDatasetTag = 'plotterSingMuonCRNoJetSelSkim' # this is just an example; it will be part of the name of the output dataset // for Poffline, Ptrigger and Nctrl

# Uncomment one of the following pairs

config.Data.outLFNDirBase = '/store/user/borzari/'
config.Site.storageSite = 'T2_BR_SPRACE'
