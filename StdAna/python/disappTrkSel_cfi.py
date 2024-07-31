import FWCore.ParameterSet.Config as cms

disappTrkSelFilter_ = cms.EDFilter ("disappTrkSel",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices",""),
    electrons = cms.InputTag("slimmedElectrons",""),
    muons = cms.InputTag("slimmedMuons",""),
    taus = cms.InputTag("slimmedTaus",""),
    tracks = cms.InputTag("isolatedTracks",""),
    ecalHitsEB = cms.InputTag("reducedEcalRecHitsEB"),
    ecalHitsEE = cms.InputTag("reducedEcalRecHitsEE"),
    hcalHits = cms.InputTag("reducedHcalRecHits", "hbhereco"),
    rhoCentralCalo = cms.InputTag("fixedGridRhoFastjetCentralCalo"),
    pfCandidates = cms.InputTag("packedPFCandidates",""),
    jets = cms.InputTag("slimmedJets",""),
    mets = cms.InputTag("slimmedMETs",""),
    triggersPAT = cms.InputTag('TriggerResults','','PAT'),
    triggersHLT = cms.InputTag('TriggerResults','','HLT'),
    HLTName = cms.string('HLT_Ele32_WPTight_Gsf_v'),
    isMETTriggers = cms.bool(True),
    isCRAB = cms.bool(False),
    trigobjs = cms.InputTag('slimmedPatTrigger'),
    ecalBadCalibReducedMINIAODFilter = cms.InputTag('ecalBadCalibReducedMINIAODFilter'),
)