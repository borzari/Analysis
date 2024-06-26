import FWCore.ParameterSet.Config as cms

zToElecProbeTrkFilter_ = cms.EDFilter ("zToElecProbeTrk",
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
    triggersPAT = cms.InputTag('TriggerResults','','PAT'),
    triggersHLT = cms.InputTag('TriggerResults','','HLT'),
    HLTName = cms.string('HLT_Ele32_WPTight_Gsf_v'),
    trigobjs = cms.InputTag('slimmedPatTrigger'),
)

zToMuonProbeTrkFilter_ = cms.EDFilter ("zToMuonProbeTrk",
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
    triggersPAT = cms.InputTag('TriggerResults','','PAT'),
    triggersHLT = cms.InputTag('TriggerResults','','HLT'),
    HLTName = cms.string('HLT_IsoMu24_v'),
    trigobjs = cms.InputTag('slimmedPatTrigger'),
)

zToTauEleProbeTrkFilter_ = cms.EDFilter ("zToTauEleProbeTrk",
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
    triggersPAT = cms.InputTag('TriggerResults','','PAT'),
    triggersHLT = cms.InputTag('TriggerResults','','HLT'),
    HLTName = cms.string('HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1_v'),
    trigobjs = cms.InputTag('slimmedPatTrigger'),
)

zToTauMuProbeTrkFilter_ = cms.EDFilter ("zToTauMuProbeTrk",
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
    triggersPAT = cms.InputTag('TriggerResults','','PAT'),
    triggersHLT = cms.InputTag('TriggerResults','','HLT'),
    HLTName = cms.string('HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1_v'),
    trigobjs = cms.InputTag('slimmedPatTrigger'),
)