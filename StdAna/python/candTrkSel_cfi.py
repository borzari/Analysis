import FWCore.ParameterSet.Config as cms

candTrkSelFilter_ = cms.EDFilter ("candTrkSel",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices",""),
    electrons = cms.InputTag("slimmedElectrons",""),
    muons = cms.InputTag("slimmedMuons",""),
    taus = cms.InputTag("slimmedTaus",""),
    tracks = cms.InputTag("isolatedTracks",""),
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