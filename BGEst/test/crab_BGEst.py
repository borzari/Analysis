#!/usr/bin/env python3

import os
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = ''
config.General.workArea = 'crab'
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = '2zToTauMProbeTrk' # this is the name of the crab folder; useful to keep track of what is happening // for Pveto
# config.General.requestName = 'singEleCRNoJetMETTrigSelSkim' # this is the name of the crab folder; useful to keep track of what is happening // for Poffline, Ptrigger and Nctrl

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'zToLepProbeTrk_cfg.py' # for Pveto
# config.JobType.psetName = 'singLepCRNoJetSelSkim_cfg.py' # for Poffline, Ptrigger and Nctrl
config.JobType.allowUndistributedCMSSW = True
# List of files that are in the repos, but are used for the selections
config.JobType.inputFiles = [os.environ['CMSSW_BASE'] + '/src/Analysis/Helper/data/electronFiducialMap_2018_data.root', os.environ['CMSSW_BASE'] + '/src/Analysis/Helper/data/muonFiducialMap_2018_data.root']

# config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/borzari-3Skimming_2022EE-3a8953c6b719d6563fe908b4ea7a6705/USER' # DY for Pveto
config.Data.inputDataset = '/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/borzari-Skimming_2022EE-3a8953c6b719d6563fe908b4ea7a6705/USER' # ttbar for Poffline, Ptrigger and Nctrl
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 # this is the amount of files processed per output file

config.Data.publication = True
config.Data.outputDatasetTag = '2zToTauMProbeTrk' # this is just an example; it will be part of the name of the output dataset // for Pveto
# config.Data.outputDatasetTag = 'singEleCRNoJetMETTrigSelSkim' # this is just an example; it will be part of the name of the output dataset // for Poffline, Ptrigger and Nctrl

# Uncomment one of the following pairs

config.Data.outLFNDirBase = '/store/user/borzari/'
config.Site.storageSite = 'T2_BR_SPRACE'
