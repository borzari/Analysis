import FWCore.ParameterSet.Config as cms

electronTagSkimFilter_ = cms.EDFilter ("electronTagSkim",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices",""),
    leptons = cms.InputTag("slimmedElectrons",""),
    triggersPAT = cms.InputTag('TriggerResults','','PAT'),
    triggersHLT = cms.InputTag('TriggerResults','','HLT'),
    HLTName = cms.string('HLT_Ele32_WPTight_Gsf_v'),
    trigobjs = cms.InputTag('slimmedPatTrigger'),
)

muonTagSkimFilter_ = cms.EDFilter ("muonTagSkim",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices",""),
    leptons = cms.InputTag("slimmedMuons",""),
    triggersPAT = cms.InputTag('TriggerResults','','PAT'),
    triggersHLT = cms.InputTag('TriggerResults','','HLT'),
    HLTName = cms.string('HLT_IsoMu24_v'),
    trigobjs = cms.InputTag('slimmedPatTrigger'),
)

tauTagSkimFilter_ = cms.EDFilter ("tauTagSkim",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices",""),
    leptons = cms.InputTag("slimmedTaus",""),
    triggersPAT = cms.InputTag('TriggerResults','','PAT'),
    triggersHLT = cms.InputTag('TriggerResults','','HLT'),
    HLTName = cms.string('HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1_v'),
    trigobjs = cms.InputTag('slimmedPatTrigger'),
)