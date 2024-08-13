import FWCore.ParameterSet.Config as cms
from Analysis.Helper.plotterSingLepCRNoJetSelSkim_cfi import *

###########################################################
##### Set up process #####
###########################################################

nEvents = 10000

reportEvery = int(nEvents/10)

isCRAB = False
# isCRAB = True

useMETTriggers = False
# useMETTriggers = True

if isCRAB:
    nEvents = -1
    reportEvery = 1

lepton = 'electron'
# lepton = 'muon'
# lepton = 'tau'

if useMETTriggers: lepton = "METTriggers_" + lepton

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process ('PLOTTERSINGLEPCRNOJETSELSKIM', Run3)
process.load ('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = reportEvery

process.maxEvents = cms.untracked.PSet (
    input = cms.untracked.int32 (nEvents)
)

process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring ("file:/home/brenoorzari/CMSSW_13_0_13/src/Analysis/BGEst/test/" + lepton + "_singLepCRNoJetSelSkim_selected.root"),
)

process.TFileService = cms.Service ('TFileService',
    fileName = cms.string (lepton + '_plotterSingLepCRNoJetSelSkim_outputHistograms.root')
)

process.options.SkipEvent = cms.untracked.vstring('ProductNotFound')

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "130X_mcRun3_2022_realistic_postEE_v6", '')

process.plotterSingElecCRNoJetSelSkimFilter = plotterSingElecCRNoJetSelSkimFilter_.clone()
process.plotterSingMuonCRNoJetSelSkimFilter = plotterSingMuonCRNoJetSelSkimFilter_.clone()
process.plotterSingTauCRNoJetSelSkimFilter = plotterSingTauCRNoJetSelSkimFilter_.clone()

if isCRAB:
    process.plotterSingElecCRNoJetSelSkimFilter.isCRAB = True
    process.plotterSingMuonCRNoJetSelSkimFilter.isCRAB = True
    process.plotterSingTauCRNoJetSelSkimFilter.isCRAB = True

if useMETTriggers:
    process.plotterSingElecCRNoJetSelSkimFilter.isMETTriggers = True
    process.plotterSingMuonCRNoJetSelSkimFilter.isMETTriggers = True
    process.plotterSingTauCRNoJetSelSkimFilter.isMETTriggers = True

if 'electron' in lepton: process.filterPath = cms.Path(process.plotterSingElecCRNoJetSelSkimFilter)
if 'muon' in lepton: process.filterPath = cms.Path(process.plotterSingMuonCRNoJetSelSkimFilter)
if 'tau' in lepton: process.filterPath = cms.Path(process.plotterSingTauCRNoJetSelSkimFilter)

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

