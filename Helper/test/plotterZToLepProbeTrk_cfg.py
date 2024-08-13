import FWCore.ParameterSet.Config as cms
from Analysis.Helper.plotterZToLepProbeTrk_cfi import *

###########################################################
##### Set up process #####
###########################################################

nEvents = 10000

reportEvery = int(nEvents/10)

isCRAB = False
# isCRAB = True

if isCRAB:
    nEvents = -1
    reportEvery = 1

# lepton = 'electron'
lepton = 'muon'
# lepton = 'tauele'
# lepton = 'taumu'

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process ('PLOTTERZTOLEPPROBETRK', Run3)
process.load ('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = reportEvery

process.maxEvents = cms.untracked.PSet (
    input = cms.untracked.int32 (nEvents)
)

process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring ("file:/home/brenoorzari/CMSSW_13_0_13/src/Analysis/BGEst/test/" + lepton + "_zToLepProbeTrk_selected.root"),
)

process.TFileService = cms.Service ('TFileService',
    fileName = cms.string (lepton + '_plotterZToLepProbeTrk_outputHistograms.root')
)

process.options.SkipEvent = cms.untracked.vstring('ProductNotFound')

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "130X_mcRun3_2022_realistic_postEE_v6", '')

process.plotterZToElecProbeTrkFilter = plotterZToElecProbeTrkFilter_.clone()
process.plotterZToMuonProbeTrkFilter = plotterZToMuonProbeTrkFilter_.clone()
process.plotterZToTauEleProbeTrkFilter = plotterZToTauEleProbeTrkFilter_.clone()
process.plotterZToTauMuProbeTrkFilter = plotterZToTauMuProbeTrkFilter_.clone()

if isCRAB:
    process.plotterZToElecProbeTrkFilter.isCRAB = True
    process.plotterZToMuonProbeTrkFilter.isCRAB = True
    process.plotterZToTauEleProbeTrkFilter.isCRAB = True
    process.plotterZToTauMuProbeTrkFilter.isCRAB = True

if lepton == 'electron': process.filterPath = cms.Path(process.plotterZToElecProbeTrkFilter)
if lepton == 'muon': process.filterPath = cms.Path(process.plotterZToMuonProbeTrkFilter)
if lepton == 'tauele': process.filterPath = cms.Path(process.plotterZToTauEleProbeTrkFilter)
if lepton == 'taumu': process.filterPath = cms.Path(process.plotterZToTauMuProbeTrkFilter)

from Configuration.EventContent.EventContent_cff import MINIAODSIMEventContent
process.EXODisappTrkSkimContent = MINIAODSIMEventContent.clone()
process.EXODisappTrkSkimContent.outputCommands.append('keep *_reducedHcalRecHits_*_*')
process.EXODisappTrkSkimContent.outputCommands.append('keep *_reducedEcalRecHits*_*_*')

# Recommended by https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#Run_3_recommendations
process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')
process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
  "EcalBadCalibFilter",
  EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEBRecHits"),
  ecalMinEt        = cms.double(50.),
  baddetEcal    = cms.vuint32([838871812]), # This is the only badly calibrated crystal in 2022 and 2023
  taggingMode = cms.bool(True),
  debug = cms.bool(False)
)
process.passecalBadCalibFilterUpdatePath = cms.Path (process.ecalBadCalibReducedMINIAODFilter)

process.schedule = cms.Schedule(process.passecalBadCalibFilterUpdatePath,process.filterPath)

