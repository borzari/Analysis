import FWCore.ParameterSet.Config as cms
import glob, sys, os
import sigLists

###########################################################
##### Set up process #####
###########################################################

process = cms.Process ('TRIGSIG')
process.load ('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet (
    input = cms.untracked.int32 (-1)
    # input = cms.untracked.int32 (10000)
)

lifetime = "1000"
mass = "900"
datayear = "2023"
dataera = ""

inputList = []
if datayear == "2022":
    if dataera == "":
        if mass == "700":
            if lifetime == "10": inputList = sigLists.AMSB_700_10_2022
            if lifetime == "100": inputList = sigLists.AMSB_700_100_2022
            if lifetime == "1000": inputList = sigLists.AMSB_700_1000_2022
    if dataera == "EE":
        if mass == "700":
            if lifetime == "10": inputList = sigLists.AMSB_700_10_2022EE
            if lifetime == "100": inputList = sigLists.AMSB_700_100_2022EE
            if lifetime == "1000": inputList = sigLists.AMSB_700_1000_2022EE
if datayear == "2023":
    if dataera == "":
        if mass == "200":
            if lifetime == "10": inputList = sigLists.AMSB_200_10_2023
            if lifetime == "100": inputList = sigLists.AMSB_200_100_2023
            if lifetime == "1000": inputList = sigLists.AMSB_200_1000_2023
        if mass == "700":
            if lifetime == "10": inputList = sigLists.AMSB_700_10_2023
            if lifetime == "100": inputList = sigLists.AMSB_700_100_2023
            if lifetime == "1000": inputList = sigLists.AMSB_700_1000_2023
        if mass == "900":
            if lifetime == "10": inputList = sigLists.AMSB_900_10_2023
            if lifetime == "100": inputList = sigLists.AMSB_900_100_2023
            if lifetime == "1000": inputList = sigLists.AMSB_900_1000_2023
        if mass == "1200":
            if lifetime == "10": inputList = sigLists.AMSB_1200_10_2023
            if lifetime == "100": inputList = sigLists.AMSB_1200_100_2023
            if lifetime == "1000": inputList = sigLists.AMSB_1200_1000_2023
    if dataera == "BPix":
        if mass == "700":
            if lifetime == "10": inputList = sigLists.AMSB_700_10_2023BPix
            if lifetime == "100": inputList = sigLists.AMSB_700_100_2023BPix
            if lifetime == "1000": inputList = sigLists.AMSB_700_1000_2023BPix

process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring (inputList),
)

process.TFileService = cms.Service ('TFileService',
    fileName = cms.string ('plots/outputFile_' + mass + '_' + lifetime + '_' + datayear + dataera + '.root')
)

gt = ""
if datayear == "2022":
    if dataera == "": gt = "130X_mcRun3_2022_realistic_v5"
    if dataera == "EE": gt = "130X_mcRun3_2022_realistic_postEE_v6"
if datayear == "2023":
    if dataera == "": gt = "130X_mcRun3_2023_realistic_v14"
    if dataera == "BPix": gt = "130X_mcRun3_2023_realistic_postBPix_v2"

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2022_realistic_postEE_v6', '')
process.GlobalTag = GlobalTag(process.GlobalTag, gt, '')

###########################################################
##### Set up the producer and the end path            #####
###########################################################

process.triggerSignalAnalyzer = cms.EDAnalyzer ("triggerSignalAnalyzer",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices",""),
    tracks = cms.InputTag("isolatedTracks",""),
    mets = cms.InputTag("slimmedMETs"),
    muons = cms.InputTag("slimmedMuons",""),
    jets = cms.InputTag("slimmedJets",""),
    triggersPAT = cms.InputTag('TriggerResults','','PAT'),
    triggersHLT = cms.InputTag('TriggerResults','','HLT'),
    trigobjs = cms.InputTag('slimmedPatTrigger'),
)

process.myPath = cms.Path (process.triggerSignalAnalyzer)
