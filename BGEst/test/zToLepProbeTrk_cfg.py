import FWCore.ParameterSet.Config as cms
from Analysis.BGEst.zToLepProbeTrk_cfi import *

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
# lepton = 'muon'
# lepton = 'tauele'
lepton = 'taumu'

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process ('ZTOLEPPROBETRK', Run3)
process.load ('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = reportEvery

process.maxEvents = cms.untracked.PSet (
    input = cms.untracked.int32 (nEvents)
)

process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring ("file:/home/brenoorzari/selected.root"),
)

process.TFileService = cms.Service ('TFileService',
    fileName = cms.string (lepton + '_zToLepProbeTrk_outputHistograms.root')
)

process.options.SkipEvent = cms.untracked.vstring('ProductNotFound')

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "130X_mcRun3_2022_realistic_postEE_v6", '')

process.zToElecProbeTrkFilter = zToElecProbeTrkFilter_.clone()
process.zToMuonProbeTrkFilter = zToMuonProbeTrkFilter_.clone()
process.zToTauEleProbeTrkFilter = zToTauEleProbeTrkFilter_.clone()
process.zToTauMuProbeTrkFilter = zToTauMuProbeTrkFilter_.clone()

if isCRAB:
    process.zToElecProbeTrkFilter.isCRAB = True
    process.zToMuonProbeTrkFilter.isCRAB = True
    process.zToTauEleProbeTrkFilter.isCRAB = True
    process.zToTauMuProbeTrkFilter.isCRAB = True

if lepton == 'electron': process.filterPath = cms.Path(process.zToElecProbeTrkFilter)
if lepton == 'muon': process.filterPath = cms.Path(process.zToMuonProbeTrkFilter)
if lepton == 'tauele': process.filterPath = cms.Path(process.zToTauEleProbeTrkFilter)
if lepton == 'taumu': process.filterPath = cms.Path(process.zToTauMuProbeTrkFilter)

from Configuration.EventContent.EventContent_cff import MINIAODSIMEventContent
process.EXODisappTrkSkimContent = MINIAODSIMEventContent.clone()
process.EXODisappTrkSkimContent.outputCommands.append('keep *_reducedHcalRecHits_*_*')
process.EXODisappTrkSkimContent.outputCommands.append('keep *_reducedEcalRecHits*_*_*')

process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('file:' + lepton + '_zToLepProbeTrk_selected.root'),
    outputCommands = process.EXODisappTrkSkimContent.outputCommands,
    overrideBranchesSplitLevel = cms.untracked.VPSet(
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJets__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVerticesWithBS__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        )
    ),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('filterPath')
    ),
)

process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)

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

process.schedule = cms.Schedule(process.passecalBadCalibFilterUpdatePath,process.filterPath,process.MINIAODSIMoutput_step)

