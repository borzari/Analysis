// -*- C++ -*-
//
// Package:    TriggerSignalAnalysis/triggerSignalAnalyzer
// Class:      triggerSignalAnalyzer
//
/**\class triggerSignalAnalyzer triggerSignalAnalyzer.cc TriggerSignalAnalysis/triggerSignalAnalyzer/plugins/triggerSignalAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Breno Orzari
//         Created:  Thu, 30 May 2024 22:50:06 GMT
//
//

// system include files
#include <memory>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "TH1D.h"
#include "TLorentzVector.h"
#include "TTree.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

// using reco::TrackCollection;

class triggerSignalAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit triggerSignalAnalyzer(const edm::ParameterSet&);
  ~triggerSignalAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  const double getTrackIsolation(const pat::IsolatedTrack&, const std::vector<pat::IsolatedTrack> &, const bool, const bool, const double, const double = 1.0e-12) const;
  bool isMatchedToTriggerObject (const edm::Event &, const edm::TriggerResults &, const pat::IsolatedTrack &, const std::vector<pat::TriggerObjectStandAlone> &, const std::string &, const std::string &, const double = 0.1);

  Int_t nJet;

  Float_t track_pt;
  Float_t track_eta;
  Float_t track_phi;
  Int_t track_highPurity;
  Int_t track_missInnHits;
  Int_t track_missMidHits;
  Int_t track_missOutHits;
  Int_t track_numValidHits;
  Int_t track_numPixelHits;
  Int_t track_numValidLayers;
  Float_t track_iso;
  Int_t track_isTrackPassingTrackLeg;

  Float_t met_pt;
  Float_t met_ptNoMu;
  Float_t met_phi;
  Float_t met_phiNoMu;

  Float_t jet_pt[1000];
  Float_t jet_eta[1000];
  Float_t jet_phi[1000];

  Int_t pass_HLT_MET105_IsoTrk50;
  Int_t pass_HLT_MET120_IsoTrk50;
  Int_t pass_HLT_PFMET105_IsoTrk50;
  Int_t pass_HLT_PFMET120_PFMHT120_IDTight;
  Int_t pass_HLT_PFMET130_PFMHT130_IDTight;
  Int_t pass_HLT_PFMET140_PFMHT140_IDTight;
  Int_t pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
  Int_t pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF;
  Int_t pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF;
  Int_t pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF;
  Int_t pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF;
  Int_t pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
  Int_t pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight;
  Int_t pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;
  Int_t pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60;
  Int_t pass_hltMET105Filter;
  Int_t pass_hltMET90Filter;
  Int_t pass_hltTrk50Filter;
  Int_t pass_hltPFMETNoMu110Filter;

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  edm::Service<TFileService> fs_;
  std::map<std::string, TH1D *> oneDHists_;
  TTree* tree_;
  double minRange = 0.0;
  double maxRange = 1000.0;
  double bins = 100;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken_;
  edm::EDGetTokenT<std::vector<pat::IsolatedTrack>> tracksToken_;
  edm::EDGetTokenT<std::vector<pat::MET>> metsToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
  edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggersPATToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggersHLTToken_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> trigobjsToken_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
triggerSignalAnalyzer::triggerSignalAnalyzer(const edm::ParameterSet& iConfig)
    : verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
      tracksToken_(consumes<std::vector<pat::IsolatedTrack>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      metsToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
      muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      triggersPATToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggersPAT"))),
      triggersHLTToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggersHLT"))),
      trigobjsToken_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("trigobjs"))) {
  //now do what ever initialization is needed

  oneDHists_["4pfmet105"] = fs_->make<TH1D>("4pfmet105", "", bins, minRange, maxRange);
  oneDHists_["4hltpfmet105"] = fs_->make<TH1D>("4hltpfmet105", "", bins, minRange, maxRange);
  oneDHists_["5pfmet105"] = fs_->make<TH1D>("5pfmet105", "", bins, minRange, maxRange);
  oneDHists_["5hltpfmet105"] = fs_->make<TH1D>("5hltpfmet105", "", bins, minRange, maxRange);
  oneDHists_["6pfmet105"] = fs_->make<TH1D>("6pfmet105", "", bins, minRange, maxRange);
  oneDHists_["6hltpfmet105"] = fs_->make<TH1D>("6hltpfmet105", "", bins, minRange, maxRange);
  oneDHists_["pfmet105"] = fs_->make<TH1D>("pfmet105", "", bins, minRange, maxRange);
  oneDHists_["hltpfmet105"] = fs_->make<TH1D>("hltpfmet105", "", bins, minRange, maxRange);

  oneDHists_["1204pfmet105"] = fs_->make<TH1D>("1204pfmet105", "", bins, minRange, maxRange);
  oneDHists_["1204hltpfmet105"] = fs_->make<TH1D>("1204hltpfmet105", "", bins, minRange, maxRange);
  oneDHists_["1205pfmet105"] = fs_->make<TH1D>("1205pfmet105", "", bins, minRange, maxRange);
  oneDHists_["1205hltpfmet105"] = fs_->make<TH1D>("1205hltpfmet105", "", bins, minRange, maxRange);
  oneDHists_["1206pfmet105"] = fs_->make<TH1D>("1206pfmet105", "", bins, minRange, maxRange);
  oneDHists_["1206hltpfmet105"] = fs_->make<TH1D>("1206hltpfmet105", "", bins, minRange, maxRange);
  oneDHists_["120pfmet105"] = fs_->make<TH1D>("120pfmet105", "", bins, minRange, maxRange);
  oneDHists_["120hltpfmet105"] = fs_->make<TH1D>("120hltpfmet105", "", bins, minRange, maxRange);

  oneDHists_["countTotal"] = fs_->make<TH1D>("countTotal", "", 1, 0.0, 1.0);
  oneDHists_["countJet"] = fs_->make<TH1D>("countJet", "", 1, 0.0, 1.0);
  oneDHists_["countEta2p5"] = fs_->make<TH1D>("countEta2p5", "", 1, 0.0, 1.0);
  oneDHists_["countIsHP"] = fs_->make<TH1D>("countIsHP", "", 1, 0.0, 1.0);
  oneDHists_["countD0"] = fs_->make<TH1D>("countD0", "", 1, 0.0, 1.0);
  oneDHists_["countDz"] = fs_->make<TH1D>("countDz", "", 1, 0.0, 1.0);
  oneDHists_["countPixHit"] = fs_->make<TH1D>("countPixHit", "", 1, 0.0, 1.0);
  oneDHists_["countHit"] = fs_->make<TH1D>("countHit", "", 1, 0.0, 1.0);
  oneDHists_["countMissInn"] = fs_->make<TH1D>("countMissInn", "", 1, 0.0, 1.0);
  oneDHists_["countMissMid"] = fs_->make<TH1D>("countMissMid", "", 1, 0.0, 1.0);
  oneDHists_["countIso"] = fs_->make<TH1D>("countIso", "", 1, 0.0, 1.0);
  oneDHists_["countHLT"] = fs_->make<TH1D>("countHLT", "", 1, 0.0, 1.0);

  nJet = 0;

  tree_ = fs_->make<TTree>("tree","tree");

  tree_->Branch("track_pt",&track_pt,"track_pt/F");
  tree_->Branch("track_eta",&track_eta,"track_eta/F");
  tree_->Branch("track_phi",&track_phi,"track_phi/F");
  tree_->Branch("track_highPurity",&track_highPurity,"track_highPurity/I");
  tree_->Branch("track_missInnHits",&track_missInnHits,"track_missInnHits/I");
  tree_->Branch("track_missMidHits",&track_missMidHits,"track_missMidHits/I");
  tree_->Branch("track_missOutHits",&track_missOutHits,"track_missOutHits/I");
  tree_->Branch("track_numValidHits",&track_numValidHits,"track_numValidHits/I");
  tree_->Branch("track_numPixelHits",&track_numPixelHits,"track_numPixelHits/I");
  tree_->Branch("track_numValidLayers",&track_numValidLayers,"track_numValidLayers/I");
  tree_->Branch("track_iso",&track_iso,"track_iso/F");
  tree_->Branch("track_isTrackPassingTrackLeg",&track_isTrackPassingTrackLeg,"track_isTrackPassingTrackLeg/I");

  tree_->Branch("met_pt",&met_pt,"met_pt/F");
  tree_->Branch("met_ptNoMu",&met_ptNoMu,"met_ptNoMu/F");
  tree_->Branch("met_phi",&met_phi,"met_phi/F");
  tree_->Branch("met_phiNoMu",&met_phiNoMu,"met_phiNoMu/F");

  tree_->Branch("nJet",&nJet,"nJet/I");
  tree_->Branch("jet_pt",jet_pt,"jet_pt/F");
  tree_->Branch("jet_eta",jet_eta,"jet_eta/F");
  tree_->Branch("jet_phi",jet_phi,"jet_phi/F");

  tree_->Branch("hlt_pass_HLT_MET105_IsoTrk50",&pass_HLT_MET105_IsoTrk50,"pass_HLT_MET105_IsoTrk50/I");
  tree_->Branch("hlt_pass_HLT_MET120_IsoTrk50",&pass_HLT_MET120_IsoTrk50,"pass_HLT_MET120_IsoTrk50/I");
  tree_->Branch("hlt_pass_HLT_PFMET105_IsoTrk50",&pass_HLT_PFMET105_IsoTrk50,"pass_HLT_PFMET105_IsoTrk50/I");
  tree_->Branch("hlt_pass_HLT_PFMET120_PFMHT120_IDTight",&pass_HLT_PFMET120_PFMHT120_IDTight,"pass_HLT_PFMET120_PFMHT120_IDTight/I");
  tree_->Branch("hlt_pass_HLT_PFMET130_PFMHT130_IDTight",&pass_HLT_PFMET130_PFMHT130_IDTight,"pass_HLT_PFMET130_PFMHT130_IDTight/I");
  tree_->Branch("hlt_pass_HLT_PFMET140_PFMHT140_IDTight",&pass_HLT_PFMET140_PFMHT140_IDTight,"pass_HLT_PFMET140_PFMHT140_IDTight/I");
  tree_->Branch("hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",&pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,"pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60/I");
  tree_->Branch("hlt_pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF",&pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF,"pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF/I");
  tree_->Branch("hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF",&pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF,"pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF/I");
  tree_->Branch("hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF",&pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF,"pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF/I");
  tree_->Branch("hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF",&pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF,"pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF/I");
  tree_->Branch("hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",&pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight,"pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight/I");
  tree_->Branch("hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight",&pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight,"pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight/I");
  tree_->Branch("hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",&pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight,"pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight/I");
  tree_->Branch("hlt_pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60",&pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60,"pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60/I");
  tree_->Branch("hlt_pass_hltMET105Filter",&pass_hltMET105Filter,"pass_hltMET105Filter/I");
  tree_->Branch("hlt_pass_hltMET90Filter",&pass_hltMET90Filter,"pass_hltMET90Filter/I");
  tree_->Branch("hlt_pass_hltTrk50Filter",&pass_hltTrk50Filter,"pass_hltTrk50Filter/I");
  tree_->Branch("hlt_pass_hltPFMETNoMu110Filter",&pass_hltPFMETNoMu110Filter,"pass_hltPFMETNoMu110Filter/I");

}

triggerSignalAnalyzer::~triggerSignalAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void triggerSignalAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  bool passEta2p5 = false;
  bool passIsHP = false;
  bool passD0 = false;
  bool passDz = false;
  bool passPixHit = false;
  bool passHit = false;
  bool passMissInn = false;
  bool passMissMid = false;
  bool passIso = false;

  track_isTrackPassingTrackLeg = 0;

  met_pt = 0.;
  met_ptNoMu = 0.;
  met_phi = 0.;
  met_phiNoMu = 0.;

  pass_HLT_MET105_IsoTrk50 = 0;
  pass_HLT_MET120_IsoTrk50 = 0;
  pass_HLT_PFMET105_IsoTrk50 = 0;
  pass_HLT_PFMET120_PFMHT120_IDTight = 0;
  pass_HLT_PFMET130_PFMHT130_IDTight = 0;
  pass_HLT_PFMET140_PFMHT140_IDTight = 0;
  pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = 0;
  pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF = 0;
  pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF = 0;
  pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF = 0;
  pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF = 0;
  pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = 0;
  pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight = 0;
  pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight = 0;
  pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60 = 0;
  pass_hltMET105Filter = 0;
  pass_hltMET90Filter = 0;
  pass_hltTrk50Filter = 0;
  pass_hltPFMETNoMu110Filter = 0;

  Int_t pass_Flag_goodVertices = 0;
  Int_t pass_Flag_HBHENoiseFilter = 0;
  Int_t pass_Flag_HBHENoiseIsoFilter = 0;
  Int_t pass_Flag_EcalDeadCellTriggerPrimitiveFilter = 0;
  Int_t pass_Flag_globalSuperTightHalo2016Filter = 0;
  Int_t pass_Flag_BadPFMuonFilter = 0;
  Int_t pass_Flag_BadPFMuonDzFilter = 0;
  Int_t pass_Flag_eeBadScFilter = 0;
  Int_t pass_Flag_ecalBadCalibFilter = 0;

  nJet = 0;

  oneDHists_.at("countTotal")->Fill(0.5);

  edm::Handle<edm::TriggerResults> triggerBitsPAT;
  iEvent.getByToken(triggersPATToken_, triggerBitsPAT);

  edm::Handle<edm::TriggerResults> triggerBitsHLT;
  iEvent.getByToken(triggersHLTToken_, triggerBitsHLT);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjs;
  iEvent.getByToken(trigobjsToken_, triggerObjs);

  for(auto triggerObj : *triggerObjs) {
    triggerObj.unpackNamesAndLabels(iEvent, *triggerBitsHLT);
    for (const auto &filterLabel : triggerObj.filterLabels ()){
      if (filterLabel == "hltPFMETNoMu110") pass_hltPFMETNoMu110Filter = 1;
      if (filterLabel == "hltMET105") pass_hltMET105Filter = 1;
      if (filterLabel == "hltMET90") pass_hltMET90Filter = 1;
      if (filterLabel == "hltTrk50Filter") pass_hltTrk50Filter = 1;
    }
  }

  const edm::TriggerNames &allTriggerNamesHLT = iEvent.triggerNames(*triggerBitsHLT);

  for(unsigned i = 0; i < allTriggerNamesHLT.size(); i++) {
    std::string thisName = allTriggerNamesHLT.triggerName(i);
    if (thisName.find("HLT_MET105_IsoTrk50_v") == 0) pass_HLT_MET105_IsoTrk50 = triggerBitsHLT->accept(i);
    if (thisName.find("HLT_MET120_IsoTrk50_v") == 0) pass_HLT_MET120_IsoTrk50 = triggerBitsHLT->accept(i);
    if (thisName.find("HLT_PFMET105_IsoTrk50_v") == 0) pass_HLT_PFMET105_IsoTrk50 = triggerBitsHLT->accept(i);
    if (thisName.find("HLT_PFMET120_PFMHT120_IDTight_v") == 0) pass_HLT_PFMET120_PFMHT120_IDTight = triggerBitsHLT->accept(i);
    if (thisName.find("HLT_PFMET130_PFMHT130_IDTight_v") == 0) pass_HLT_PFMET130_PFMHT130_IDTight = triggerBitsHLT->accept(i);
    if (thisName.find("HLT_PFMET140_PFMHT140_IDTight_v") == 0) pass_HLT_PFMET140_PFMHT140_IDTight = triggerBitsHLT->accept(i);
    if (thisName.find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v") == 0) pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = triggerBitsHLT->accept(i);
    if (thisName.find("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF_v") == 0) pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF = triggerBitsHLT->accept(i);
    if (thisName.find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF_v") == 0) pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF = triggerBitsHLT->accept(i);
    if (thisName.find("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF_v") == 0) pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF = triggerBitsHLT->accept(i);
    if (thisName.find("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF_v") == 0) pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF = triggerBitsHLT->accept(i);
    if (thisName.find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") == 0) pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = triggerBitsHLT->accept(i);
    if (thisName.find("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v") == 0) pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight = triggerBitsHLT->accept(i);
    if (thisName.find("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v") == 0) pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight = triggerBitsHLT->accept(i);
    if (thisName.find("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v") == 0) pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60 = triggerBitsHLT->accept(i);
  }

  const edm::TriggerNames &allTriggerNamesPAT = iEvent.triggerNames(*triggerBitsPAT);

  for(unsigned i = 0; i < allTriggerNamesPAT.size(); i++) {
    std::string thisName = allTriggerNamesPAT.triggerName(i);
    if (thisName.find("Flag_goodVertices") == 0) pass_Flag_goodVertices = triggerBitsPAT->accept(i);
    if (thisName.find("Flag_HBHENoiseFilter") == 0) pass_Flag_HBHENoiseFilter = triggerBitsPAT->accept(i);
    if (thisName.find("Flag_HBHENoiseIsoFilter") == 0) pass_Flag_HBHENoiseIsoFilter = triggerBitsPAT->accept(i);
    if (thisName.find("Flag_EcalDeadCellTriggerPrimitiveFilter") == 0) pass_Flag_EcalDeadCellTriggerPrimitiveFilter = triggerBitsPAT->accept(i);
    if (thisName.find("Flag_globalSuperTightHalo2016Filter") == 0) pass_Flag_globalSuperTightHalo2016Filter = triggerBitsPAT->accept(i);
    if (thisName.find("Flag_BadPFMuonFilter") == 0) pass_Flag_BadPFMuonFilter = triggerBitsPAT->accept(i);
    if (thisName.find("Flag_BadPFMuonDzFilter") == 0) pass_Flag_BadPFMuonDzFilter = triggerBitsPAT->accept(i);
    if (thisName.find("Flag_eeBadScFilter") == 0) pass_Flag_eeBadScFilter = triggerBitsPAT->accept(i);
    if (thisName.find("Flag_ecalBadCalibFilter") == 0) pass_Flag_ecalBadCalibFilter = triggerBitsPAT->accept(i);
  }

  if(pass_Flag_goodVertices == 1 && pass_Flag_HBHENoiseFilter == 1 && pass_Flag_HBHENoiseIsoFilter == 1 && pass_Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && pass_Flag_globalSuperTightHalo2016Filter == 1 && pass_Flag_BadPFMuonFilter == 1 && pass_Flag_BadPFMuonDzFilter == 1 && pass_Flag_eeBadScFilter == 1 && pass_Flag_ecalBadCalibFilter == 1){

    edm::Handle<std::vector<reco::Vertex> > vertices;
    iEvent.getByToken(verticesToken_, vertices);
    const reco::Vertex &pv = vertices->at(0);

    bool isGoodJet = false;
    bool isGoodTrack = false;

    bool isNLayers4 = false;
    bool isNLayers5 = false;
    bool isNLayers6p = false;

    int nJet_aux = 0;
    for (const auto& jet : iEvent.get(jetsToken_)) {

      if(nJet_aux != 0) break;
      if(abs(jet.eta()) < 2.4 || jet.eta() < -998.9) isGoodJet = true;
      nJet_aux = nJet_aux + 1;

    }

    if(isGoodJet){

      oneDHists_.at("countJet")->Fill(0.5);

      edm::Handle<std::vector<pat::MET> > pfmet;
      iEvent.getByToken(metsToken_, pfmet);

      TLorentzVector pfMetNoMu;
      pfMetNoMu.SetPtEtaPhiM(pfmet->at(0).pt(),0.0,pfmet->at(0).phi(),0.0);

      for (const auto& muon : iEvent.get(muonsToken_)) {

        if(abs(muon.pt()) < 55.) continue;
        if(abs(muon.eta()) > 2.1) continue;
        if(!muon.isTightMuon(pv)) continue;
        if(muon.innerTrack()->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) != 0) continue;
        if(muon.innerTrack()->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS) != 0) continue;
        if(((muon.pfIsolationR04().sumChargedHadronPt + std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt()) > 0.15) continue;

        TLorentzVector mu;
        mu.SetPtEtaPhiM(muon.pt(),muon.eta(),muon.phi(),0.1057128);
        pfMetNoMu = (pfMetNoMu + mu);

      }

      pat::IsolatedTrack passedTrack;

      for (const auto& track : iEvent.get(tracksToken_)) {

        if(abs(track.eta()) > 2.5) continue;
        passEta2p5 = true;
        if(!(track.isHighPurityTrack())) continue;
        passIsHP = true;
        if(abs(track.dxy()) > 0.02) continue;
        passD0 = true;
        if(abs(track.dz()) > 0.5) continue;
        passDz = true;
        if(track.hitPattern().numberOfValidPixelHits() < 4) continue;
        passPixHit = true;
        if(track.hitPattern().numberOfValidHits() < 4) continue;
        passHit = true;
        if(track.lostInnerLayers() != 0) continue;
        passMissInn = true;
        if(track.lostLayers() != 0) continue;
        passMissMid = true;
        if((getTrackIsolation(track, iEvent.get(tracksToken_), true, false, 0.3))/track.pt() > 0.01) continue;
        passIso = true;

        // if(isMatchedToTriggerObject(iEvent, *triggerBitsHLT, track, *triggerObjs, "hltMergedTracks::HLT", "hltTrk50Filter")) passedTrack = track;

        if(track.hitPattern().trackerLayersWithMeasurement() == 4) {
          isNLayers4 = true;
        }
        if(track.hitPattern().trackerLayersWithMeasurement() == 5) {
          isNLayers5 = true;
        }
        if(track.hitPattern().trackerLayersWithMeasurement() >= 6) {
          isNLayers6p = true;
        }

        if(isGoodTrack) continue;

        isGoodTrack = true;
        passedTrack = track;

      }

      if(passEta2p5) oneDHists_.at("countEta2p5")->Fill(0.5);
      if(passIsHP) oneDHists_.at("countIsHP")->Fill(0.5);
      if(passD0) oneDHists_.at("countD0")->Fill(0.5);
      if(passDz) oneDHists_.at("countDz")->Fill(0.5);
      if(passPixHit) oneDHists_.at("countPixHit")->Fill(0.5);
      if(passHit) oneDHists_.at("countHit")->Fill(0.5);
      if(passMissInn) oneDHists_.at("countMissInn")->Fill(0.5);
      if(passMissMid) oneDHists_.at("countMissMid")->Fill(0.5);
      if(passIso) oneDHists_.at("countIso")->Fill(0.5);

      if(isGoodTrack){

        if(isNLayers4) oneDHists_.at("4pfmet105")->Fill(pfMetNoMu.Pt());
        if(isNLayers5) oneDHists_.at("5pfmet105")->Fill(pfMetNoMu.Pt());
        if(isNLayers6p) oneDHists_.at("6pfmet105")->Fill(pfMetNoMu.Pt());
        oneDHists_.at("pfmet105")->Fill(pfMetNoMu.Pt());
        if(pfMetNoMu.Pt() > 120.0){
          if(isNLayers4) oneDHists_.at("1204pfmet105")->Fill(pfMetNoMu.Pt());
          if(isNLayers5) oneDHists_.at("1205pfmet105")->Fill(pfMetNoMu.Pt());
          if(isNLayers6p) oneDHists_.at("1206pfmet105")->Fill(pfMetNoMu.Pt());
          oneDHists_.at("120pfmet105")->Fill(pfMetNoMu.Pt());
        }

        if(pass_HLT_MET105_IsoTrk50 == 1){
          oneDHists_.at("countHLT")->Fill(0.5);
          if(isNLayers4) oneDHists_.at("4hltpfmet105")->Fill(pfMetNoMu.Pt());
          if(isNLayers5) oneDHists_.at("5hltpfmet105")->Fill(pfMetNoMu.Pt());
          if(isNLayers6p) oneDHists_.at("6hltpfmet105")->Fill(pfMetNoMu.Pt());
          oneDHists_.at("hltpfmet105")->Fill(pfMetNoMu.Pt());
          if(pfMetNoMu.Pt() > 120.0){
            if(isNLayers4) oneDHists_.at("1204hltpfmet105")->Fill(pfMetNoMu.Pt());
            if(isNLayers5) oneDHists_.at("1205hltpfmet105")->Fill(pfMetNoMu.Pt());
            if(isNLayers6p) oneDHists_.at("1206hltpfmet105")->Fill(pfMetNoMu.Pt());
            oneDHists_.at("120hltpfmet105")->Fill(pfMetNoMu.Pt());
          }
        }

        met_pt = pfmet->at(0).pt();
        met_ptNoMu = pfMetNoMu.Pt();
        met_phi = pfmet->at(0).phi();
        met_phiNoMu = pfMetNoMu.Phi();

        for (const auto& jet : iEvent.get(jetsToken_)) {
          nJet = nJet + 1;
          jet_pt[nJet - 1] = jet.pt();
          jet_eta[nJet - 1] = jet.eta();
          jet_phi[nJet - 1] = jet.phi();
        }

        track_pt = passedTrack.pt();
        track_eta = passedTrack.eta();
        track_phi = passedTrack.phi();
        track_highPurity = passedTrack.isHighPurityTrack();
        track_missInnHits = passedTrack.hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS);
        track_missMidHits = passedTrack.hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
        track_missOutHits = passedTrack.hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_OUTER_HITS);
        track_numValidHits = passedTrack.hitPattern().numberOfValidTrackerHits();
        track_numPixelHits = passedTrack.hitPattern().numberOfValidPixelHits();
        track_numValidLayers = passedTrack.hitPattern().trackerLayersWithMeasurement();
        track_iso = (getTrackIsolation(passedTrack, iEvent.get(tracksToken_), true, false, 0.3))/passedTrack.pt();
        if(isMatchedToTriggerObject(iEvent, *triggerBitsHLT, passedTrack, *triggerObjs, "hltMergedTracks::HLT", "hltTrk50Filter")) track_isTrackPassingTrackLeg = 1;
        // std::cout << "====================================== End of event ======================================" << std::endl;

        tree_->Fill();

      }

    }

  }
  
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void triggerSignalAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("isolatedTracks"));
  desc.add<edm::InputTag>("mets", edm::InputTag("slimmedMETs"));
  desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("triggersPAT", edm::InputTag("TriggerResults","","PAT"));
  desc.add<edm::InputTag>("triggersHLT", edm::InputTag("TriggerResults","","HLT"));
  desc.add<edm::InputTag>("trigobjs", edm::InputTag("slimmedPatTrigger"));

  descriptions.addWithDefaultLabel(desc);
}

const double triggerSignalAnalyzer::getTrackIsolation (const pat::IsolatedTrack& track, const std::vector<pat::IsolatedTrack> &tracks, const bool noPU, const bool noFakes, const double outerDeltaR, const double innerDeltaR) const
{
  double sumPt = 0.0;
  for (const auto &t : tracks)
    {
      if (noFakes && (t.hitPattern().pixelLayersWithMeasurement() < 2 || t.hitPattern().trackerLayersWithMeasurement() < 5 || fabs(t.dxy() / t.dxyError()) > 5.0)) continue;
      if (noPU && track.dz() > 3.0 * hypot(track.dzError(), t.dzError())) continue;
      double dR = deltaR (track, t);
      if (dR < outerDeltaR && dR > innerDeltaR) sumPt += t.pt ();
    }

  return sumPt;
}

bool triggerSignalAnalyzer::isMatchedToTriggerObject (const edm::Event &event, const edm::TriggerResults &triggers, const pat::IsolatedTrack &obj, const std::vector<pat::TriggerObjectStandAlone> &trigObjs, const std::string &collection, const std::string &filter, const double dR)
{
  // std::cout << "This is my collection: " << collection << std::endl;
  // std::cout << "This is my filter: " << filter << std::endl;
  if(collection == "") return false;
  for(auto trigObj : trigObjs) {
    trigObj.unpackNamesAndLabels(event, triggers);
    // if(filter == "hltTrk50Filter") std::cout << trigObj.collection() << std::endl;
    if(trigObj.collection() != collection) continue;
    if(filter != "") {
      bool flag = false;
      for(const auto &filterLabel : trigObj.filterLabels ()) {
        if(filterLabel == filter) {
          flag = true;
          break;
        }
      }
      if (!flag) continue;
    }
    // if(filter == "hltTrk50Filter") std::cout << trigObj.eta() << std::endl;
    if(deltaR (obj, trigObj) > dR) continue;
    return true;
  }
  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(triggerSignalAnalyzer);
