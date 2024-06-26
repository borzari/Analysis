import FWCore.ParameterSet.Config as cms
from Analysis.BGEst.zToLepProbeTrk_cfi import *

###########################################################
##### Set up process #####
###########################################################

nEvents = 100000

# lepton = 'electron'
lepton = 'muon'
# lepton = 'taue'
# lepton = 'taum'

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process ('MUONTAGSKIM', Run3)
process.load ('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = int(nEvents/10)

process.maxEvents = cms.untracked.PSet (
    input = cms.untracked.int32 (nEvents)
)

process.source = cms.Source ("PoolSource",
    # fileNames = cms.untracked.vstring ("root://osg-se.sprace.org.br:1094//store/user/borzari/WtoLNu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8/Skimming_2022EE/240620_202923/0000/selected_1.root"),
    fileNames = cms.untracked.vstring ("root://osg-se.sprace.org.br:1094//store/user/borzari/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/3Skimming_2022EE/240627_204809/0000/selected_1.root","root://osg-se.sprace.org.br:1094//store/user/borzari/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/3Skimming_2022EE/240627_204809/0000/selected_10.root","root://osg-se.sprace.org.br:1094//store/user/borzari/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/3Skimming_2022EE/240627_204809/0000/selected_100.root","root://osg-se.sprace.org.br:1094//store/user/borzari/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/3Skimming_2022EE/240627_204809/0000/selected_101.root","root://osg-se.sprace.org.br:1094//store/user/borzari/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/3Skimming_2022EE/240627_204809/0000/selected_102.root"),
    # fileNames = cms.untracked.vstring ("root://osg-se.sprace.org.br:1094//store/user/borzari/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/2Skimming_2022EE/240526_160303/0000/skim_Skimming_10.root","root://osg-se.sprace.org.br:1094//store/user/borzari/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/2Skimming_2022EE/240526_160303/0000/skim_Skimming_100.root","root://osg-se.sprace.org.br:1094//store/user/borzari/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/2Skimming_2022EE/240526_160303/0000/skim_Skimming_101.root"), // Can't run over skims produced by original analysis code as it has custom types that don't exist here
)

process.TFileService = cms.Service ('TFileService',
    fileName = cms.string (lepton + '_outputHistograms.root')
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

if lepton == 'electron': process.filterPath = cms.Path(process.zToElecProbeTrkFilter)
if lepton == 'muon': process.filterPath = cms.Path(process.zToMuonProbeTrkFilter)
if lepton == 'taue': process.filterPath = cms.Path(process.zToTauEleProbeTrkFilter)
if lepton == 'taum': process.filterPath = cms.Path(process.zToTauMuProbeTrkFilter)

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
    fileName = cms.untracked.string('file:' + lepton + '_selected.root'),
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

# process.CaloGeometryBuilder.SelectedCalos = cms.vstring(
#     'HCAL',
#     'ZDC',
#     'EcalBarrel',
#     'EcalEndcap',
#     'EcalPreshower',
#     'TOWER'
# )

process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)

process.schedule = cms.Schedule(process.filterPath,process.MINIAODSIMoutput_step)

