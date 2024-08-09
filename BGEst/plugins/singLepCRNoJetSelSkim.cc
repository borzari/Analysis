// -*- C++ -*-
//
// Package:    BGEst/BGEst
// Class:      BGEst
//
/**\class BGEst BGEst.cc BGEst/BGEst/plugins/BGEst.cc

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
#include <iomanip>
#include <random>
#include <chrono>
#include <cstdlib>
#include <ctime>

#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector2.h"
#include "TGraphAsymmErrors.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "Analysis/Helper/interface/helperFunctions.h"
#include "Analysis/Helper/interface/plotPrintFunctions.h"
#include "Analysis/Helper/interface/selectingFunctions.h"
//
// class declaration
//

// If the filter does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDFilter<>
// This will improve performance in multithreaded jobs.

template<class T>
class singLepCRNoJetSelSkim : public edm::one::EDFilter<edm::one::SharedResources> {
public:
  explicit singLepCRNoJetSelSkim(const edm::ParameterSet&);
  ~singLepCRNoJetSelSkim() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  std::vector<std::string> commonCuts, electronCuts, muonCuts, tauCuts, trackCuts;

private:
  bool filter(edm::Event&, const edm::EventSetup&) override;
  edm::Service<TFileService> fs_;
  std::map<std::string, TH1D *> oneDHists_;
  std::map<std::string, TH2D *> twoDHists_;
  std::map<std::string, TGraphAsymmErrors *> oneDGraph_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken_;
  edm::EDGetTokenT<std::vector<T>> leptonsToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggersPATToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggersHLTToken_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> trigobjsToken_;
  edm::EDGetTokenT<std::vector<pat::IsolatedTrack>> tracksToken_;
  edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken_;
  edm::EDGetTokenT<std::vector<pat::MET>> metsToken_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> pfCandToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  edm::ESGetToken<EcalChannelStatus, EcalChannelStatusRcd> ecalStatusToken_;
  edm::EDGetTokenT<bool> ecalBadCalibFilterUpdateToken_;
  std::string HLTName_;
  bool isCRAB_;
  bool isMETTriggers_;

  edm::ESHandle<CaloGeometry> caloGeometry;
  edm::ESHandle<EcalChannelStatus> ecalStatus;

  EtaPhiList vetoListElec;
  EtaPhiList vetoListMu;

  std::map<DetId, std::vector<double> > EcalAllDeadChannelsValMap;
  std::map<DetId, std::vector<int> >    EcalAllDeadChannelsBitMap;

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
template<class T>
singLepCRNoJetSelSkim<T>::singLepCRNoJetSelSkim(const edm::ParameterSet& iConfig)
    : verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
      leptonsToken_(consumes<std::vector<T>>(iConfig.getParameter<edm::InputTag>("leptons"))),
      triggersPATToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggersPAT"))),
      triggersHLTToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggersHLT"))),
      trigobjsToken_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("trigobjs"))),
      tracksToken_(consumes<std::vector<pat::IsolatedTrack>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      metsToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
      pfCandToken_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("pfCandidates"))),
      caloGeometryToken_(esConsumes()),
      ecalStatusToken_(esConsumes()),
      ecalBadCalibFilterUpdateToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("ecalBadCalibReducedMINIAODFilter"))) {
  //now do what ever initialization is needed

  HLTName_ = iConfig.getParameter<std::string>("HLTName");
  isCRAB_ = iConfig.getParameter<bool>("isCRAB");
  isMETTriggers_ = iConfig.getParameter<bool>("isMETTriggers");

  std::string HLTCut = HLTName_;

  if(isMETTriggers_) HLTCut = HLTCut + " && MET Triggers";

  std::string elecFile = "";
  std::string muonFile = "";

  if(isCRAB_){
    elecFile = "electronFiducialMap_2018_data.root";
    muonFile = "muonFiducialMap_2018_data.root";
  }
  else{
    elecFile = "/home/brenoorzari/CMSSW_13_0_13/src/Analysis/Helper/data/electronFiducialMap_2018_data.root";
    muonFile = "/home/brenoorzari/CMSSW_13_0_13/src/Analysis/Helper/data/muonFiducialMap_2018_data.root";
  }

  helperFunctions::extractFiducialMap<pat::Electron>(vetoListElec,elecFile);
  helperFunctions::extractFiducialMap<pat::Muon>(vetoListMu,muonFile);

  commonCuts = {
    "Total",
    HLTCut,
    "METFilters",
    "passecalBadCalibFilterUpdate"
  };

  electronCuts = {
    "electron pt > 35 GeV",
    "electron fabs(#eta) < 2.1",
    "electron passesVID_tightID (ID + iso)",
    "electron |d0| < 0.05, 0.10 (EB, EE)",
    "electron |dz| < 0.10, 0.20 (EB, EE)"
  };

  muonCuts = {
    "muon pt > 26 GeV",
    "muon fabs(#eta) < 2.1",
    "muon isTightMuonWRTVtx",
    "muon #Delta#beta-corrected rel. iso. < 0.15"
  };

  tauCuts = {
    "tau pt > 50 GeV",
    "tau fabs(#eta) < 2.1",
    "tau passesDecayModeReconstruction && passesLightFlavorRejection"
  };

  std::string trackPt = "";
  std::string trackLepDeltaR = "";
  std::string trackMatchToLep = "";
  std::string trackRandom = "";

  if(std::is_same<T, pat::Electron>::value) {trackPt = "track pt > 35 GeV"; trackLepDeltaR = "deltaR (track, electron) < 0.1"; trackMatchToLep = "match track to electron"; trackRandom = "pick random electron";}
  if(std::is_same<T, pat::Muon>::value) {trackPt = "track pt > 35 GeV"; trackLepDeltaR = "deltaR (track, muon) < 0.1"; trackMatchToLep = "match track to muon"; trackRandom = "pick random muon";}
  if(std::is_same<T, pat::Tau>::value) {trackPt = "track pt > 55 GeV"; trackLepDeltaR = "deltaR (track, tau) < 0.1"; trackMatchToLep = "match track to tau"; trackRandom = "pick random tau";}

  trackCuts = {
    trackPt,
    trackLepDeltaR,
    "track fabs(#eta) < 2.1",
    "track fabs ( eta ) < 1.42 || fabs ( eta ) > 1.65",
    "track fabs ( eta ) < 0.15 || fabs ( eta ) > 0.35",
    "track fabs ( eta ) < 1.55 || fabs ( eta ) > 1.85",
    "track !inTOBCrack",
    "isFiducialElectronTrack",
    "isFiducialMuonTrack",
    "isFiducialECALTrack",
    "track hitPattern_.numberOfValidPixelHits >= 4",
    "track hitPattern_.numberOfValidHits >= 4",
    "track missingInnerHits == 0",
    "track hitDrop_missingMiddleHits == 0",
    "track ((pfIsolationDR03_.chargedHadronIso + pfIsolationDR03_.puChargedHadronIso) / pt) < 0.05",
    "track |d0| < 0.02",
    "track |dz| < 0.5",
    "track dRMinJet > 0.5",
    trackMatchToLep,
    trackRandom // This might only be applied when selecting the correct event from the skim
  };

  if(std::is_same<T, pat::Tau>::value){
    std::string elementToRemove = "track dRMinJet > 0.5";

    // Remove the element using erase function and iterators
    auto it = std::find(trackCuts.begin(), trackCuts.end(), elementToRemove);

    // If element is found found, erase it
    if (it != trackCuts.end()) trackCuts.erase(it);
  }

  int nBins = int(commonCuts.size() + trackCuts.size());
  double maxRange = 0.0;
  std::vector<std::string> lepCuts;

  if(std::is_same<T, pat::Electron>::value){nBins = nBins + int(electronCuts.size()); maxRange = double(nBins) - 0.5; lepCuts = electronCuts;}
  if(std::is_same<T, pat::Muon>::value){nBins = nBins + int(muonCuts.size()); maxRange = double(nBins) - 0.5; lepCuts = muonCuts;}
  if(std::is_same<T, pat::Tau>::value){nBins = nBins + int(tauCuts.size()); maxRange = double(nBins) - 0.5; lepCuts = tauCuts;}

  oneDHists_["cutflow"] = fs_->make<TH1D>("cutflow", "", nBins, -0.5, maxRange);
  oneDHists_["selection"] = fs_->make<TH1D>("selection", "", nBins, -0.5, maxRange);

  for(int i = 1; i <= int(commonCuts.size() + lepCuts.size() + trackCuts.size()); i++){
    if(i <= int(commonCuts.size())) {
      oneDHists_.at("cutflow")->GetXaxis()->SetBinLabel(i,commonCuts[i-1].c_str());
      oneDHists_.at("selection")->GetXaxis()->SetBinLabel(i,commonCuts[i-1].c_str());
    }
    if(i > int(commonCuts.size()) && i <= int(commonCuts.size() + lepCuts.size())) {
      oneDHists_.at("cutflow")->GetXaxis()->SetBinLabel(i,lepCuts[i-1-int(commonCuts.size())].c_str());
      oneDHists_.at("selection")->GetXaxis()->SetBinLabel(i,lepCuts[i-1-int(commonCuts.size())].c_str());
    }
    if(i > int(commonCuts.size() + lepCuts.size())) {
      oneDHists_.at("cutflow")->GetXaxis()->SetBinLabel(i,trackCuts[i-1-int(lepCuts.size())-int(commonCuts.size())].c_str());
      oneDHists_.at("selection")->GetXaxis()->SetBinLabel(i,trackCuts[i-1-int(lepCuts.size())-int(commonCuts.size())].c_str());
    }
  }

  oneDHists_["hist_metNoMu"] = fs_->make<TH1D>("hist_metNoMu", "", 200, 0.0, 1000.0);
  oneDHists_["hist_metNoMu_passMETTriggers"] = fs_->make<TH1D>("hist_metNoMu_passMETTriggers", "", 200, 0.0, 1000.0);
  twoDHists_["hist_metNoMuvsDeltaPhiJetMetNoMu"] = fs_->make<TH2D>("hist_metNoMuvsDeltaPhiJetMetNoMu", "", 200, 0.0, 1000.0, 64, 0.0, 3.2);

}

template<class T>
singLepCRNoJetSelSkim<T>::~singLepCRNoJetSelSkim() {

  plotPrintFunctions::printCuts(oneDHists_);

}

//
// member functions
//

// ------------ method called for each event  ------------
template<class T>
bool singLepCRNoJetSelSkim<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  std::vector<bool> passSel;
  std::vector<bool> passCut;

  for(int j = 0; j < (int(commonCuts.size()) - 1); ++j) {passSel.push_back(false); passCut.push_back(false);}

  bool isGood = false;

  caloGeometry = iSetup.getHandle(caloGeometryToken_);
  ecalStatus = iSetup.getHandle(ecalStatusToken_);

  helperFunctions::getChannelStatusMaps(ecalStatus,caloGeometry,EcalAllDeadChannelsValMap,EcalAllDeadChannelsBitMap);

  edm::Handle<std::vector<pat::Electron>> electrons;
  if(std::is_same<T, pat::Electron>::value) iEvent.getByToken(leptonsToken_, electrons);  

  edm::Handle<std::vector<pat::Muon>> muons;
  if(std::is_same<T, pat::Muon>::value) iEvent.getByToken(leptonsToken_, muons);  

  edm::Handle<std::vector<pat::Tau>> taus;
  if(std::is_same<T, pat::Tau>::value) iEvent.getByToken(leptonsToken_, taus);

  edm::Handle<std::vector<pat::IsolatedTrack>> tracks;
  iEvent.getByToken(tracksToken_, tracks);

  edm::Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(jetsToken_, jets);
  pat::Jet leadingJet;
  for (const auto &jet : *jets) {
    if(!helperFunctions::IsValidJet(jet)) continue;
    leadingJet = jet;
    break;
  }

  edm::Handle<std::vector<pat::MET> > mets;
  iEvent.getByToken(metsToken_, mets);
  const pat::MET &met = mets->at(0);

  edm::Handle<std::vector<pat::PackedCandidate>> pfCandidates;
  iEvent.getByToken(pfCandToken_, pfCandidates);

  TVector2 metNoMu = helperFunctions::calcMetNoMu(met, pfCandidates);

  edm::Handle<bool> passecalBadCalibFilterUpdate;
  iEvent.getByToken(ecalBadCalibFilterUpdateToken_, passecalBadCalibFilterUpdate);

  oneDHists_.at("selection")->Fill(0);
  oneDHists_.at("cutflow")->Fill(0);

  edm::Handle<edm::TriggerResults> triggerBitsPAT;
  iEvent.getByToken(triggersPATToken_, triggerBitsPAT);

  edm::Handle<edm::TriggerResults> triggerBitsHLT;
  iEvent.getByToken(triggersHLTToken_, triggerBitsHLT);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjs;
  iEvent.getByToken(trigobjsToken_, triggerObjs);

  std::vector<pat::Electron> selElectrons;
  std::vector<pat::Muon> selMuons;
  std::vector<pat::Tau> selTaus;

  std::vector<pat::Electron> tkMatchElectrons;
  std::vector<pat::Muon> tkMatchMuons;
  std::vector<pat::Tau> tkMatchTaus;

  std::vector<pat::IsolatedTrack> selTracks;

  if(!isMETTriggers_) {if(helperFunctions::passLepHLTPath(iEvent,triggerBitsHLT,HLTName_)) {passSel[0] = true; passCut[0] = true;}}

  if(isMETTriggers_) {if(helperFunctions::passLepHLTPath(iEvent,triggerBitsHLT,HLTName_) && helperFunctions::passMETHLTPath(iEvent,triggerBitsHLT)) {passSel[0] = true; passCut[0] = true;}}

  if(helperFunctions::passMETFilters(iEvent,triggerBitsPAT)) {passSel[1] = true; if(passCut[0]) passCut[1] = true;}

  if(*passecalBadCalibFilterUpdate) {passSel[2] = true; if(passCut[1]) passCut[2] = true;}

  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(verticesToken_, vertices);
  const reco::Vertex &pv = vertices->at(0);

  std::vector<bool> auxPassCut;

  if(std::is_same<T, pat::Electron>::value){

    for(int j = 0; j < int(electronCuts.size()); ++j) {passSel.push_back(false); passCut.push_back(false); auxPassCut.push_back(false);}

    int startElecIdx = int(commonCuts.size()) - 1;

    for (const auto& electron : *electrons) {

      selectingFunctions::singElecSel(passSel,passCut,auxPassCut,electron,startElecIdx,pv,selElectrons,35.0);
      
      for(int j = 0; j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }

    }
  
  }

  if(std::is_same<T, pat::Muon>::value){

    for(int j = 0; j < int(muonCuts.size()); ++j) {passSel.push_back(false); passCut.push_back(false); auxPassCut.push_back(false);}

    int startMuonIdx = int(commonCuts.size()) - 1;
    
    for (const auto& muon : *muons) {
      
      selectingFunctions::singMuonSel(passSel,passCut,auxPassCut,muon,startMuonIdx,pv,selMuons);
      
      for(int j = 0; j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }
        
    }
  
  }

  if(std::is_same<T, pat::Tau>::value){

    for(int j = 0; j < int(tauCuts.size()); ++j) {passSel.push_back(false); passCut.push_back(false); auxPassCut.push_back(false);}

    int startTauIdx = int(commonCuts.size()) - 1;
    
    for (const auto& tau : *taus) {
      
      selectingFunctions::singTauSel(passSel,passCut,auxPassCut,tau,startTauIdx,selTaus);

      for(int j = 0; j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }
        
    }
  
  }

  for(int j = 0; j < int(trackCuts.size()); ++j) {passSel.push_back(false); passCut.push_back(false); auxPassCut.push_back(false);}

  int startTrackIdx = 0;
  if(std::is_same<T, pat::Electron>::value) startTrackIdx = int(commonCuts.size() + electronCuts.size()) - 1;
  if(std::is_same<T, pat::Muon>::value) startTrackIdx = int(commonCuts.size() + muonCuts.size()) - 1;
  if(std::is_same<T, pat::Tau>::value) startTrackIdx = int(commonCuts.size() + tauCuts.size()) - 1;

  double cutTrackPt = 0.0;
  if(std::is_same<T, pat::Electron>::value || std::is_same<T, pat::Muon>::value) cutTrackPt = 35.0;
  if(std::is_same<T, pat::Tau>::value) cutTrackPt = 55.0;

  int getStTrkIdxLep = int(commonCuts.size()) - 1;

  for (const auto& track : *tracks) {

    if(std::is_same<T, pat::Electron>::value){

      selectingFunctions::singTrackSel<pat::Electron>(passSel, passCut, auxPassCut, track, startTrackIdx, getStTrkIdxLep, tkMatchElectrons, selElectrons, electrons, cutTrackPt, vetoListElec, vetoListMu, EcalAllDeadChannelsValMap, EcalAllDeadChannelsBitMap, jets, selTracks);

    }

    if(std::is_same<T, pat::Muon>::value){

      selectingFunctions::singTrackSel<pat::Muon>(passSel, passCut, auxPassCut, track, startTrackIdx, getStTrkIdxLep, tkMatchMuons, selMuons, muons, cutTrackPt, vetoListElec, vetoListMu, EcalAllDeadChannelsValMap, EcalAllDeadChannelsBitMap, jets, selTracks);

    }

    if(std::is_same<T, pat::Tau>::value){

      selectingFunctions::singTrackSel<pat::Tau>(passSel, passCut, auxPassCut, track, startTrackIdx, getStTrkIdxLep, tkMatchTaus, selTaus, taus, cutTrackPt, vetoListElec, vetoListMu, EcalAllDeadChannelsValMap, EcalAllDeadChannelsBitMap, jets, selTracks);

    }

    if(std::is_same<T, pat::Electron>::value){
      for(int j = int(electronCuts.size()); j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }
    }

    if(std::is_same<T, pat::Muon>::value){
      for(int j = int(muonCuts.size()); j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }
    }

    if(std::is_same<T, pat::Tau>::value){
      for(int j = int(tauCuts.size()); j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }
    }

  }

  TVector2 metNoMuNoLep;
  double deltaPhiMetJetLeadingVsLeptonMetNoMuMinusOnePt;

  if(passCut[int(passSel.size())-1]){
    std::srand(std::time(0)); // use current time as seed for random generator

    if(std::is_same<T, pat::Electron>::value){
      int random_pos = std::rand() % tkMatchElectrons.size();  // Modulo to restrict the number of random values to be at most A.size()-1
      pat::Electron randomLep = tkMatchElectrons[random_pos];
      metNoMuNoLep = TVector2(metNoMu.Px() + randomLep.px(), metNoMu.Py() + randomLep.py());
      deltaPhiMetJetLeadingVsLeptonMetNoMuMinusOnePt = fabs(deltaPhi(metNoMuNoLep.Phi(),leadingJet.phi()));
      if(int(tkMatchElectrons.size()) > 0) twoDHists_.at("hist_metNoMuvsDeltaPhiJetMetNoMu")->Fill(metNoMuNoLep.Mod(),deltaPhiMetJetLeadingVsLeptonMetNoMuMinusOnePt);
      if(isMETTriggers_){
        oneDHists_.at("hist_metNoMu")->Fill(metNoMuNoLep.Mod());
        std::vector<bool> passes = helperFunctions::passMETTriggers(iEvent,triggerBitsHLT,triggerObjs,electrons,"hltEgammaCandidates::HLT");
        if(passes[0]) oneDHists_.at("hist_metNoMu_passMETTriggers")->Fill(metNoMuNoLep.Mod());
      }
    }

    if(std::is_same<T, pat::Muon>::value){
      int random_pos = std::rand() % tkMatchMuons.size();  // Modulo to restrict the number of random values to be at most A.size()-1
      pat::Muon randomLep = tkMatchMuons[random_pos];
      metNoMuNoLep = TVector2(metNoMu.Px(), metNoMu.Py());
      deltaPhiMetJetLeadingVsLeptonMetNoMuMinusOnePt = fabs(deltaPhi(metNoMuNoLep.Phi(),leadingJet.phi()));
      if(int(tkMatchMuons.size()) > 0) twoDHists_.at("hist_metNoMuvsDeltaPhiJetMetNoMu")->Fill(metNoMuNoLep.Mod(),deltaPhiMetJetLeadingVsLeptonMetNoMuMinusOnePt);
      if(isMETTriggers_){
        oneDHists_.at("hist_metNoMu")->Fill(metNoMuNoLep.Mod());
        std::vector<bool> passes = helperFunctions::passMETTriggers(iEvent,triggerBitsHLT,triggerObjs,muons,"hltL3MuonCandidates::HLT");
        if(passes[0]) oneDHists_.at("hist_metNoMu_passMETTriggers")->Fill(metNoMuNoLep.Mod());
      }
    }

    if(std::is_same<T, pat::Tau>::value){
      int random_pos = std::rand() % tkMatchTaus.size();  // Modulo to restrict the number of random values to be at most A.size()-1
      pat::Tau randomLep = tkMatchTaus[random_pos];
      metNoMuNoLep = TVector2(metNoMu.Px() + randomLep.px(), metNoMu.Py() + randomLep.py());
      deltaPhiMetJetLeadingVsLeptonMetNoMuMinusOnePt = fabs(deltaPhi(metNoMuNoLep.Phi(),leadingJet.phi()));
      if(int(tkMatchTaus.size()) > 0) twoDHists_.at("hist_metNoMuvsDeltaPhiJetMetNoMu")->Fill(metNoMuNoLep.Mod(),deltaPhiMetJetLeadingVsLeptonMetNoMuMinusOnePt);
      if(isMETTriggers_){
        oneDHists_.at("hist_metNoMu")->Fill(metNoMuNoLep.Mod());
        std::vector<bool> passes = helperFunctions::passMETTriggers(iEvent,triggerBitsHLT,triggerObjs,taus,"hltHpsSinglePFTau20MediumDitauWPDeepTauNoMatchForVBFIsoTau::HLT");
        if(passes[0]) oneDHists_.at("hist_metNoMu_passMETTriggers")->Fill(metNoMuNoLep.Mod());
      }
    }
    
  }

  for(int i = 1; i <= int(passSel.size()); i++){
    if(passSel[i-1]) oneDHists_.at("selection")->Fill(i);
    if(passCut[i-1]) oneDHists_.at("cutflow")->Fill(i);
    if(i == int(passSel.size()) && passCut[i-1]) isGood = true;
  }

  return isGood;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<class T>
void singLepCRNoJetSelSkim<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("leptons", edm::InputTag("placeholderLeptons"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("isolatedTracks"));
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("mets", edm::InputTag("slimmedMETs"));
  desc.add<edm::InputTag>("pfCandidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("triggersPAT", edm::InputTag("TriggerResults","","PAT"));
  desc.add<edm::InputTag>("triggersHLT", edm::InputTag("TriggerResults","","HLT"));
  desc.add<std::string>("HLTName", std::string("placeholderHLT"));
  desc.add<bool>("isCRAB", bool(false));
  desc.add<bool>("isMETTriggers", bool(false));
  desc.add<edm::InputTag>("trigobjs", edm::InputTag("slimmedPatTrigger"));
  desc.add<edm::InputTag>("ecalBadCalibReducedMINIAODFilter", edm::InputTag("ecalBadCalibReducedMINIAODFilter"));

  descriptions.addWithDefaultLabel(desc);

}

using singElecCRNoJetSelSkim = singLepCRNoJetSelSkim<pat::Electron>;
using singMuonCRNoJetSelSkim = singLepCRNoJetSelSkim<pat::Muon>;
using singTauCRNoJetSelSkim = singLepCRNoJetSelSkim<pat::Tau>;

//define this as a plug-in
// DEFINE_FWK_MODULE(singLepCRNoJetSelSkim);
DEFINE_FWK_MODULE(singElecCRNoJetSelSkim);
DEFINE_FWK_MODULE(singMuonCRNoJetSelSkim);
DEFINE_FWK_MODULE(singTauCRNoJetSelSkim);
