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

#include "TH1D.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"

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

//
// class declaration
//

// If the filter does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDFilter<>
// This will improve performance in multithreaded jobs.

template<char const *T>
class zToLepProbeTrk : public edm::one::EDFilter<edm::one::SharedResources> {
public:
  explicit zToLepProbeTrk(const edm::ParameterSet&);
  ~zToLepProbeTrk() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  std::vector<std::string> commonCuts, electronCuts, muonCuts, taueCuts, taumCuts, trackCuts;

private:
  bool filter(edm::Event&, const edm::EventSetup&) override;
  edm::Service<TFileService> fs_;
  std::map<std::string, TH1D *> oneDHists_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken_;
  edm::EDGetTokenT<std::vector<pat::Electron>> electronsToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
  edm::EDGetTokenT<std::vector<pat::Tau>> tausToken_;
  edm::EDGetTokenT<std::vector<pat::IsolatedTrack>> tracksToken_;
  edm::EDGetTokenT<EBRecHitCollection> ecalHitsEBToken_;
  edm::EDGetTokenT<EERecHitCollection> ecalHitsEEToken_;
  edm::EDGetTokenT<HBHERecHitCollection> hcalHitsToken_;
  edm::EDGetTokenT<double> rhoCentralCaloToken_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> pfCandToken_;
  edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken_;
  edm::EDGetTokenT<std::vector<pat::MET>> metsToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggersPATToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggersHLTToken_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> trigobjsToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  edm::ESGetToken<EcalChannelStatus, EcalChannelStatusRcd> ecalStatusToken_;
  edm::EDGetTokenT<bool> ecalBadCalibFilterUpdateToken_;
  std::string HLTName_;

  edm::ESHandle<CaloGeometry> caloGeometry;
  edm::ESHandle<EcalChannelStatus> ecalStatus;

  int nTPOS = 0;
  int nTPSS = 0;
  int nTPOS_veto = 0;
  int nTPSS_veto = 0;
  int nTPOSLoose_veto = 0;
  int nTPSSLoose_veto = 0;

  EtaPhiList vetoListElec;
  EtaPhiList vetoListMu;

  std::map<DetId, std::vector<double> > EcalAllDeadChannelsValMap;
  std::map<DetId, std::vector<int> >    EcalAllDeadChannelsBitMap;

};

template<char const *T>
zToLepProbeTrk<T>::zToLepProbeTrk(const edm::ParameterSet& iConfig)
    : verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
      electronsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
      muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      tausToken_(consumes<std::vector<pat::Tau>>(iConfig.getParameter<edm::InputTag>("taus"))),
      tracksToken_(consumes<std::vector<pat::IsolatedTrack>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      ecalHitsEBToken_(consumes<EBRecHitCollection>(iConfig.getParameter<edm::InputTag>("ecalHitsEB"))),
      ecalHitsEEToken_(consumes<EERecHitCollection>(iConfig.getParameter<edm::InputTag>("ecalHitsEE"))),
      hcalHitsToken_(consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hcalHits"))),
      rhoCentralCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCentralCalo"))),
      pfCandToken_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("pfCandidates"))),
      jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      metsToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
      triggersPATToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggersPAT"))),
      triggersHLTToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggersHLT"))),
      trigobjsToken_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("trigobjs"))),
      caloGeometryToken_(esConsumes()),
      ecalStatusToken_(esConsumes()),
      ecalBadCalibFilterUpdateToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("ecalBadCalibReducedMINIAODFilter"))) {
  //now do what ever initialization is needed

  HLTName_ = iConfig.getParameter<std::string>("HLTName");

  helperFunctions::extractFiducialMap<pat::Electron>(vetoListElec);
  helperFunctions::extractFiducialMap<pat::Muon>(vetoListMu);

  commonCuts = {
    "Total",
    HLTName_,
    "METFilters",
    "passecalBadCalibFilterUpdate"
  };

  electronCuts = {
    "electron pt > 32 GeV",
    "electron isMatchedToTriggerObject",
    "electron fabs(#eta) < 2.1",
    "electron passesVID_tightID (ID + iso)",
    "electron |d0| < 0.05, 0.10 (EB, EE)",
    "electron |dz| < 0.10, 0.20 (EB, EE)"
  };

  muonCuts = {
    "muon pt > 26 GeV",
    "muon isMatchedToTriggerObject",
    "muon fabs(#eta) < 2.1",
    "muon isTightMuonWRTVtx",
    "muon #Delta#beta-corrected rel. iso. < 0.15"
  };

  taueCuts = electronCuts;
  taueCuts.push_back("electron-mets with transMass (electron, met) < 40");

  taumCuts = muonCuts;
  taumCuts.push_back("muon-mets with transMass (muon, met) < 40");

  std::string closestLep1 = "";
  std::string closestLep2 = "";

  if(strcmp(T, "electron") == 0){closestLep1 = "track deltaRToClosestMuon > 0.15"; closestLep2 = "track deltaRToClosestTauHad > 0.15";}
  if(strcmp(T, "muon") == 0){closestLep1 = "track deltaRToClosestElectron > 0.15"; closestLep2 = "track deltaRToClosestTauHad > 0.15";}
  if(strcmp(T, "tauele") == 0 || strcmp(T, "taumu") == 0){closestLep1 = "track deltaRToClosestMuon > 0.15"; closestLep2 = "track deltaRToClosestElectron > 0.15";}

  trackCuts = {
    "track pt > 30 GeV",
    "track fabs(#eta) < 2.1",
    "track fabs ( eta ) < 1.42 || fabs ( eta ) > 1.65",
    "track fabs ( eta ) < 0.15 || fabs ( eta ) > 0.35",
    "track fabs ( eta ) < 1.55 || fabs ( eta ) > 1.85",
    "track !inTOBCrack",
    //
    "isFiducialElectronTrack",
    "isFiducialMuonTrack",
    "isFiducialECALTrack",
    //
    "track hitPattern_.numberOfValidPixelHits >= 4",
    "track hitPattern_.numberOfValidHits >= 4",
    "track missingInnerHits == 0",
    "track hitDrop_missingMiddleHits == 0",
    "track ((pfIsolationDR03_.chargedHadronIso + pfIsolationDR03_.puChargedHadronIso) / pt) < 0.05",
    "track |d0| < 0.02","track |dz| < 0.5",
    "track dRMinJet > 0.5",
    closestLep1,
    closestLep2
  };

  if(strcmp(T, "muon") == 0) trackCuts.push_back("track caloNewNoPUDRp5CentralCalo < 10");
  if(strcmp(T, "tauele") == 0 || strcmp(T, "taumu") == 0){
    std::string elementToRemove = "track dRMinJet > 0.5";

    // Remove the element using erase function and iterators
    auto it = std::find(trackCuts.begin(), trackCuts.end(), elementToRemove);

    // If element is found found, erase it
    if (it != trackCuts.end()) trackCuts.erase(it);
  }

  int nBins = int(commonCuts.size() + trackCuts.size());
  double maxRange = 0.0;
  std::vector<std::string> lepCuts;

  if(strcmp(T, "electron") == 0){nBins = nBins + int(electronCuts.size()); maxRange = double(nBins) - 0.5; lepCuts = electronCuts;}
  if(strcmp(T, "muon") == 0){nBins = nBins + int(muonCuts.size()); maxRange = double(nBins) - 0.5; lepCuts = muonCuts;}
  if(strcmp(T, "tauele") == 0){nBins = nBins + int(taueCuts.size()); maxRange = double(nBins) - 0.5; lepCuts = taueCuts;}
  if(strcmp(T, "taumu") == 0){nBins = nBins + int(taumCuts.size()); maxRange = double(nBins) - 0.5; lepCuts = taumCuts;}

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

  oneDHists_["hist_nTPOS"] = fs_->make<TH1D>("hist_nTPOS", "", 10, -0.5, 9.5);
  oneDHists_["hist_nTPSS"] = fs_->make<TH1D>("hist_nTPSS", "", 10, -0.5, 9.5);
  oneDHists_["hist_nTPOS_veto"] = fs_->make<TH1D>("hist_nTPOS_veto", "", 10, -0.5, 9.5);
  oneDHists_["hist_nTPSS_veto"] = fs_->make<TH1D>("hist_nTPSS_veto", "", 10, -0.5, 9.5);
  oneDHists_["hist_nTPOSLoose_veto"] = fs_->make<TH1D>("hist_nTPOSLoose_veto", "", 10, -0.5, 9.5);
  oneDHists_["hist_nTPSSLoose_veto"] = fs_->make<TH1D>("hist_nTPSSLoose_veto", "", 10, -0.5, 9.5);

}

template<char const *T>
zToLepProbeTrk<T>::~zToLepProbeTrk() {

  int maxSize = 0;
  for(int i = 1; i <= int(oneDHists_.at("cutflow")->GetNbinsX()); i++){
    std::string label = oneDHists_.at("cutflow")->GetXaxis()->GetLabels()->At(i-1)->GetName();
    if(int(label.size()) > maxSize) maxSize = int(label.size());
  }

  std::string cut = "Cut Name";
  std::cout << cut << std::string((maxSize - int(cut.size()) + 1),' ') << "\t\tCumul.\t\tIndiv." << std::endl;

  double maxValue = double(oneDHists_.at("selection")->GetBinContent(1));
  for(int i = 1; i <= int(oneDHists_.at("cutflow")->GetNbinsX()); i++){
    std::cout << std::setprecision(3) << std::fixed;
    std::string label = oneDHists_.at("cutflow")->GetXaxis()->GetLabels()->At(i-1)->GetName();
    int size = int(label.size());
    int nBlank = maxSize - size + 1;
    std::cout << oneDHists_.at("cutflow")->GetXaxis()->GetLabels()->At(i-1)->GetName() << std::string(nBlank,' ') << "\t\t" << double(oneDHists_.at("cutflow")->GetBinContent(i)) * 100.0 / maxValue << " % \t" << double(oneDHists_.at("selection")->GetBinContent(i)) * 100.0 / maxValue << " %" << std::endl;
  }

  std::cout << "# OS T&P pairs before veto: " << nTPOS << std::endl;
  std::cout << "# SS T&P pairs before veto: " << nTPSS << std::endl;
  std::cout << "# OS T&P pairs after veto: " << nTPOS_veto << std::endl;
  std::cout << "# SS T&P pairs after veto: " << nTPSS_veto << std::endl;
  std::cout << "# OS T&P pairs after loose veto: " << nTPOSLoose_veto << std::endl;
  std::cout << "# SS T&P pairs after loose veto: " << nTPSSLoose_veto << std::endl;
}

// ------------ method called for each event  ------------
template<char const *T>
bool zToLepProbeTrk<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  std::vector<bool> passSel;
  std::vector<bool> passCut;

  for(int j = 0; j < (int(commonCuts.size()) - 1); ++j) {passSel.push_back(false); passCut.push_back(false);}

  bool isGood = false;

  int hist_nTPOS = 0;
  int hist_nTPSS = 0;
  int hist_nTPOS_veto = 0;
  int hist_nTPSS_veto = 0;
  int hist_nTPOSLoose_veto = 0;
  int hist_nTPSSLoose_veto = 0;

  caloGeometry = iSetup.getHandle(caloGeometryToken_);
  ecalStatus = iSetup.getHandle(ecalStatusToken_);

  helperFunctions::getChannelStatusMaps(ecalStatus,caloGeometry,EcalAllDeadChannelsValMap,EcalAllDeadChannelsBitMap);

  edm::Handle<std::vector<pat::IsolatedTrack>> tracks;
  iEvent.getByToken(tracksToken_, tracks);

  edm::Handle<std::vector<pat::PackedCandidate>> pfCandidates;
  iEvent.getByToken(pfCandToken_, pfCandidates);

  edm::Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(jetsToken_, jets);

  edm::Handle<std::vector<pat::MET> > mets;
  iEvent.getByToken(metsToken_, mets);
  const pat::MET &met = mets->at(0);

  edm::Handle<std::vector<pat::Electron>> electrons;
  iEvent.getByToken(electronsToken_, electrons);  

  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);  

  edm::Handle<std::vector<pat::Tau>> taus;
  iEvent.getByToken(tausToken_, taus);

  edm::Handle<EBRecHitCollection> EBRecHits;
  iEvent.getByToken(ecalHitsEBToken_, EBRecHits);

  edm::Handle<EERecHitCollection> EERecHits;
  iEvent.getByToken(ecalHitsEEToken_, EERecHits);

  edm::Handle<HBHERecHitCollection> HBHERecHits;
  iEvent.getByToken(hcalHitsToken_, HBHERecHits);

  edm::Handle<double> rhoCentralCalo;
  iEvent.getByToken (rhoCentralCaloToken_, rhoCentralCalo);

  edm::Handle<bool> passecalBadCalibFilterUpdate;
  iEvent.getByToken(ecalBadCalibFilterUpdateToken_, passecalBadCalibFilterUpdate);

  oneDHists_.at("selection")->Fill(0);
  oneDHists_.at("cutflow")->Fill(0);

  std::vector<pat::Electron> electronTags;
  std::vector<pat::Muon> muonTags;
  std::vector<pat::Electron> taueTags;
  std::vector<pat::Muon> taumTags;
  std::vector<pat::IsolatedTrack> trackProbes;

  edm::Handle<edm::TriggerResults> triggerBitsPAT;
  iEvent.getByToken(triggersPATToken_, triggerBitsPAT);

  edm::Handle<edm::TriggerResults> triggerBitsHLT;
  iEvent.getByToken(triggersHLTToken_, triggerBitsHLT);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjs;
  iEvent.getByToken(trigobjsToken_, triggerObjs);

  if(helperFunctions::passHLTPath(iEvent,triggerBitsHLT,HLTName_)) {passSel[0] = true; passCut[0] = true;}

  if(helperFunctions::passMETFilters(iEvent,triggerBitsPAT)) {passSel[1] = true; if(passCut[0]) passCut[1] = true;}

  if(*passecalBadCalibFilterUpdate) {passSel[2] = true; if(passCut[1]) passCut[2] = true;}

  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(verticesToken_, vertices);
  const reco::Vertex &pv = vertices->at(0);

  std::vector<bool> auxPassCut;

  if(strcmp(T, "electron") == 0){

    for(int j = 0; j < int(electronCuts.size()); ++j) {passSel.push_back(false); passCut.push_back(false); auxPassCut.push_back(false);}

    int startElecIdx = int(commonCuts.size()) - 1;
    
    for (const auto& electron : *electrons) {
      
      int cutIdxInc = 0;

      if(electron.pt() > 32.)
        {passSel[startElecIdx+cutIdxInc] = true; if(passCut[startElecIdx-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(helperFunctions::isMatchedToTriggerObject<pat::Electron> (iEvent, *triggerBitsHLT, electron, *triggerObjs, "hltEgammaCandidates::HLT", "hltEle32WPTightGsfTrackIsoFilter"))
        {passSel[startElecIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(abs(electron.eta()) < 2.1)
        {passSel[startElecIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-tight"))
        {passSel[startElecIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(helperFunctions::elecD0(electron, pv))
        {passSel[startElecIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(helperFunctions::elecDZ(electron, pv))
        {passSel[startElecIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      if(auxPassCut[cutIdxInc]) electronTags.push_back(electron);
      
      for(int j = 0; j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }

    }
  
  }

  if(strcmp(T, "muon") == 0){
    
    for(int j = 0; j < int(muonCuts.size()); ++j) {passSel.push_back(false); passCut.push_back(false); auxPassCut.push_back(false);}

    int startMuonIdx = int(commonCuts.size()) - 1;

    for (const auto& muon : *muons) {

      int cutIdxInc = 0;
      
      if(muon.pt() > 26.)
        {passSel[startMuonIdx+cutIdxInc] = true; if(passCut[startMuonIdx-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(helperFunctions::isMatchedToTriggerObject<pat::Muon> (iEvent, *triggerBitsHLT, muon, *triggerObjs, "hltIterL3MuonCandidates::HLT", "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered"))
        {passSel[startMuonIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(abs(muon.eta()) < 2.1)
        {passSel[startMuonIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(muon.isTightMuon(pv))
        {passSel[startMuonIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(helperFunctions::muonIso(muon) < 0.15)
        {passSel[startMuonIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      if(auxPassCut[cutIdxInc]) muonTags.push_back(muon);

      for(int j = 0; j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }

    }
  
  }

  if(strcmp(T, "tauele") == 0){

    for(int j = 0; j < int(taueCuts.size()); ++j) {passSel.push_back(false); passCut.push_back(false); auxPassCut.push_back(false);}

    int startTauIdx = int(commonCuts.size()) - 1;
    
    for (const auto& electron : *electrons) {
      
      int cutIdxInc = 0;

      if(electron.pt() > 32.)
        {passSel[startTauIdx+cutIdxInc] = true; if(passCut[startTauIdx-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(helperFunctions::isMatchedToTriggerObject<pat::Electron> (iEvent, *triggerBitsHLT, electron, *triggerObjs, "hltEgammaCandidates::HLT", "hltEle32WPTightGsfTrackIsoFilter"))
        {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(abs(electron.eta()) < 2.1)
        {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-tight"))
        {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(helperFunctions::elecD0(electron, pv))
        {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(helperFunctions::elecDZ(electron, pv))
        {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(helperFunctions::transvMassLepton<pat::Electron>(electron,met) < 40.0)
        {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      if(auxPassCut[cutIdxInc]) taueTags.push_back(electron);
      
      for(int j = 0; j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }
        
    }
  
  }

  if(strcmp(T, "taumu") == 0){

    for(int j = 0; j < int(taumCuts.size()); ++j) {passSel.push_back(false); passCut.push_back(false); auxPassCut.push_back(false);}

    int startTauIdx = int(commonCuts.size()) - 1;
    
    for (const auto& muon : *muons) {

      int cutIdxInc = 0;
      
      if(muon.pt() > 26.)
        {passSel[startTauIdx+cutIdxInc] = true; if(passCut[startTauIdx-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(helperFunctions::isMatchedToTriggerObject<pat::Muon> (iEvent, *triggerBitsHLT, muon, *triggerObjs, "hltIterL3MuonCandidates::HLT", "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered"))
        {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(abs(muon.eta()) < 2.1)
        {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(muon.isTightMuon(pv))
        {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(helperFunctions::muonIso(muon) < 0.15)
        {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      ++cutIdxInc;

      if(helperFunctions::transvMassLepton<pat::Muon>(muon,met) < 40.0)
        {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      if(auxPassCut[cutIdxInc]) taumTags.push_back(muon);

      for(int j = 0; j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }

    }
  
  }

// >= 1 tracks with isFiducialElectronTrack                                                                       1983.0         19.830%        100.000%
// >= 1 tracks with isFiducialMuonTrack                                                                           1983.0         19.830%        100.000%
// >= 1 tracks with isFiducialECALTrack                                                                           1863.0         18.630%         98.170%

  for(int j = 0; j < int(trackCuts.size()); ++j) {passSel.push_back(false); passCut.push_back(false); auxPassCut.push_back(false);}

  int startTrackIdx = 0;
  if(strcmp(T, "electron") == 0) startTrackIdx = int(commonCuts.size() + electronCuts.size()) - 1;
  if(strcmp(T, "muon") == 0) startTrackIdx = int(commonCuts.size() + muonCuts.size()) - 1;
  if(strcmp(T, "tauele") == 0) startTrackIdx = int(commonCuts.size() + taueCuts.size()) - 1;
  if(strcmp(T, "taumu") == 0) startTrackIdx = int(commonCuts.size() + taumCuts.size()) - 1;

  int getStTrkIdxLep = int(commonCuts.size()) - 1;

  for (const auto& track : *tracks) {

    int cutIdxInc = 0;
      
    if(track.pt() > 30.) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(passCut[startTrackIdx+cutIdxInc-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(fabs(track.eta()) < 2.1) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if((fabs(track.eta()) < 1.42) || (fabs(track.eta()) > 1.65)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if((fabs(track.eta()) < 0.15) || (fabs(track.eta()) > 0.35)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if((fabs(track.eta()) < 1.55) || (fabs(track.eta()) > 1.85)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(!helperFunctions::inTOBCrack(track)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    // Need to include fiducial cuts here!!!

    if(helperFunctions::isFiducialTrack(track,vetoListElec,0.05,-1.0))
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(helperFunctions::isFiducialTrack(track,vetoListMu,0.05,-1.0))
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(!helperFunctions::isCloseToBadEcalChannel(track,0.05,EcalAllDeadChannelsValMap,EcalAllDeadChannelsBitMap))
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    // Need to include fiducial cuts here!!!

    if(track.hitPattern().numberOfValidPixelHits() >= 4) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(track.hitPattern().numberOfValidHits() >= 4) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(track.lostInnerLayers() == 0) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(helperFunctions::hitDrop_missingMiddleHits(track) == 0) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(((track.pfIsolationDR03().chargedHadronIso() + track.pfIsolationDR03().puChargedHadronIso()) / track.pt()) < 0.05) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(fabs(track.dxy()) < 0.02) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(fabs(track.dz()) < 0.5) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    if(strcmp(T, "electron") == 0 || strcmp(T, "muon") == 0)
      {
        ++cutIdxInc;

        if(helperFunctions::dRMinJet(track, *jets) > 0.5)
          {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}
      }

    if(strcmp(T, "electron") == 0 || strcmp(T, "tauele") == 0 || strcmp(T, "taumu") == 0)
      {
        ++cutIdxInc;

        if(helperFunctions::deltaRToClosestLepton<pat::Muon>(track, *muons) > 0.15)
          {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}
      }

    if(strcmp(T, "muon") == 0 || strcmp(T, "tauele") == 0 || strcmp(T, "taumu") == 0)
      {
        ++cutIdxInc;

        if(helperFunctions::deltaRToClosestLepton<pat::Electron>(track, *electrons) > 0.15)
          {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}
      }

    if(strcmp(T, "electron") == 0 || strcmp(T, "muon") == 0) 
      {
        ++cutIdxInc;

        if(helperFunctions::deltaRToClosestTauHad(track, *taus) > 0.15)
          {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

      } // The individual efficiency of this is different than the original analysis, because it uses distinct selections that follow the Run 3 recommendations from the Tau POG https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3#Kinematic_tau_selection
    
    if(strcmp(T, "muon") == 0){

      ++cutIdxInc;

      if(helperFunctions::caloNewNoPUDRp5CentralCalo(track, *EBRecHits, *EERecHits, *HBHERecHits, *rhoCentralCalo, caloGeometry) < 10.0) 
        {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    }

    if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep]) trackProbes.push_back(track);

    if(strcmp(T, "electron") == 0){
      for(int j = int(electronCuts.size()); j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }
    }

    if(strcmp(T, "muon") == 0){
      for(int j = int(muonCuts.size()); j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }
    }

    if(strcmp(T, "tauele") == 0){
      for(int j = int(taueCuts.size()); j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }
    }

    if(strcmp(T, "taumu") == 0){
      for(int j = int(taumCuts.size()); j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }
    }

  }

  if(strcmp(T, "electron") == 0){
    for(const auto &tag : electronTags){
      for(const auto &probe : trackProbes){
        if(helperFunctions::goodInvMassLepton<pat::Electron>(tag, probe, false) && (tag.charge()*probe.charge()) < 0.0){
          ++nTPOS;
          ++hist_nTPOS;
          if(helperFunctions::passesVeto<pat::Electron>(probe, *pfCandidates, *jets)){
            ++nTPOS_veto;
            ++hist_nTPOS_veto;
          }
          if(helperFunctions::passesLooseElecVeto(probe, *electrons, pv)){
            ++nTPOSLoose_veto;
            ++hist_nTPOSLoose_veto;
          }
        }
        if(helperFunctions::goodInvMassLepton<pat::Electron>(tag, probe, false) && (tag.charge()*probe.charge()) > 0.0){
          ++nTPSS;
          ++hist_nTPSS;
          if(helperFunctions::passesVeto<pat::Electron>(probe, *pfCandidates, *jets)){
            ++nTPSS_veto;
            ++hist_nTPSS_veto;
          }
          if(helperFunctions::passesLooseElecVeto(probe, *electrons, pv)){
            ++nTPSSLoose_veto;
            ++hist_nTPSSLoose_veto;
          }
        }
      }
    }
  }

  if(strcmp(T, "muon") == 0){
    for(const auto &tag : muonTags){
      for(const auto &probe : trackProbes){
        if(helperFunctions::goodInvMassLepton<pat::Muon>(tag, probe, false) && (tag.charge()*probe.charge()) < 0.0){
          ++nTPOS;
          ++hist_nTPOS;
          if(helperFunctions::passesVeto<pat::Muon>(probe, *pfCandidates, *jets)){
            ++nTPOS_veto;
            ++hist_nTPOS_veto;
          }
          if(helperFunctions::passesLooseMuonVeto(probe, *muons)){
            ++nTPOSLoose_veto;
            ++hist_nTPOSLoose_veto;
          }
        }
        if(helperFunctions::goodInvMassLepton<pat::Muon>(tag, probe, false) && (tag.charge()*probe.charge()) > 0.0){
          ++nTPSS;
          ++hist_nTPSS;
          if(helperFunctions::passesVeto<pat::Muon>(probe, *pfCandidates, *jets)){
            ++nTPSS_veto;
            ++hist_nTPSS_veto;
          }
          if(helperFunctions::passesLooseMuonVeto(probe, *muons)){
            ++nTPSSLoose_veto;
            ++hist_nTPSSLoose_veto;
          }
        }
      }
    }
  }

  if(strcmp(T, "tauele") == 0){
    for(const auto &tag : taueTags){
      for(const auto &probe : trackProbes){
        if(helperFunctions::goodInvMassLepton<pat::Electron>(tag, probe, true) && (tag.charge()*probe.charge()) < 0.0){
          ++nTPOS;
          ++hist_nTPOS;
          if(helperFunctions::passesVeto<pat::Tau>(probe, *pfCandidates, *jets)){
            ++nTPOS_veto;
            ++hist_nTPOS_veto;
          }
          if(helperFunctions::passesLooseElecVeto(probe, *electrons, pv)){
            ++nTPOSLoose_veto;
            ++hist_nTPOSLoose_veto;
          }
        }
        if(helperFunctions::goodInvMassLepton<pat::Electron>(tag, probe, true) && (tag.charge()*probe.charge()) > 0.0){
          ++nTPSS;
          ++hist_nTPSS;
          if(helperFunctions::passesVeto<pat::Tau>(probe, *pfCandidates, *jets)){
            ++nTPSS_veto;
            ++hist_nTPSS_veto;
          }
          if(helperFunctions::passesLooseElecVeto(probe, *electrons, pv)){
            ++nTPSSLoose_veto;
            ++hist_nTPSSLoose_veto;
          }
        }
      }
    }
  }

  if(strcmp(T, "taumu") == 0){
    for(const auto &tag : taumTags){
      for(const auto &probe : trackProbes){
        if(helperFunctions::goodInvMassLepton<pat::Muon>(tag, probe, true) && (tag.charge()*probe.charge()) < 0.0){
          ++nTPOS;
          ++hist_nTPOS;
          if(helperFunctions::passesVeto<pat::Tau>(probe, *pfCandidates, *jets)){
            ++nTPOS_veto;
            ++hist_nTPOS_veto;
          }
          if(helperFunctions::passesLooseMuonVeto(probe, *muons)){
            ++nTPOSLoose_veto;
            ++hist_nTPOSLoose_veto;
          }
        }
        if(helperFunctions::goodInvMassLepton<pat::Muon>(tag, probe, true) && (tag.charge()*probe.charge()) > 0.0){
          ++nTPSS;
          ++hist_nTPSS;
          if(helperFunctions::passesVeto<pat::Tau>(probe, *pfCandidates, *jets)){
            ++nTPSS_veto;
            ++hist_nTPSS_veto;
          }
          if(helperFunctions::passesLooseMuonVeto(probe, *muons)){
            ++nTPSSLoose_veto;
            ++hist_nTPSSLoose_veto;
          }
        }
      }
    }
  }

  if(passCut[int(passSel.size())-1]){
    oneDHists_.at("hist_nTPOS")->Fill(hist_nTPOS);
    oneDHists_.at("hist_nTPSS")->Fill(hist_nTPSS);
    oneDHists_.at("hist_nTPOS_veto")->Fill(hist_nTPOS_veto);
    oneDHists_.at("hist_nTPSS_veto")->Fill(hist_nTPSS_veto);
    oneDHists_.at("hist_nTPOSLoose_veto")->Fill(hist_nTPOSLoose_veto);
    oneDHists_.at("hist_nTPSSLoose_veto")->Fill(hist_nTPSSLoose_veto);
  }

  for(int i = 1; i <= int(passSel.size()); i++){
    if(passSel[i-1]) oneDHists_.at("selection")->Fill(i);
    if(passCut[i-1]) oneDHists_.at("cutflow")->Fill(i);
    if(i == int(passSel.size()) && passCut[i-1]) isGood = true;
  }

  return isGood;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<char const *T>
void zToLepProbeTrk<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("electrons", edm::InputTag("slimmedElectrons"));
  desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("taus", edm::InputTag("slimmedTaus"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("isolatedTracks"));
  desc.add<edm::InputTag>("ecalHitsEB", edm::InputTag("reducedEcalRecHitsEB"));
  desc.add<edm::InputTag>("ecalHitsEE", edm::InputTag("reducedEcalRecHitsEE"));
  desc.add<edm::InputTag>("hcalHits", edm::InputTag("reducedHcalRecHits", "hbhereco"));
  desc.add<edm::InputTag>("rhoCentralCalo", edm::InputTag("fixedGridRhoFastjetCentralCalo"));
  desc.add<edm::InputTag>("pfCandidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("mets", edm::InputTag("slimmedMETs"));
  desc.add<edm::InputTag>("triggersPAT", edm::InputTag("TriggerResults","","PAT"));
  desc.add<edm::InputTag>("triggersHLT", edm::InputTag("TriggerResults","","HLT"));
  desc.add<std::string>("HLTName", std::string("placeholderHLT"));
  desc.add<edm::InputTag>("trigobjs", edm::InputTag("slimmedPatTrigger"));
  desc.add<edm::InputTag>("ecalBadCalibReducedMINIAODFilter", edm::InputTag("ecalBadCalibReducedMINIAODFilter"));

  descriptions.addWithDefaultLabel(desc);

}

extern char const charElectron[] = "electron";
extern char const charMuon[] = "muon";
extern char const charTauele[] = "tauele";
extern char const charTaumu[] = "taumu";

using zToElecProbeTrk = zToLepProbeTrk<charElectron>;
using zToMuonProbeTrk = zToLepProbeTrk<charMuon>;
using zToTauEleProbeTrk = zToLepProbeTrk<charTauele>;
using zToTauMuProbeTrk = zToLepProbeTrk<charTaumu>;

//define this as a plug-in
DEFINE_FWK_MODULE(zToElecProbeTrk);
DEFINE_FWK_MODULE(zToMuonProbeTrk);
DEFINE_FWK_MODULE(zToTauEleProbeTrk);
DEFINE_FWK_MODULE(zToTauMuProbeTrk);
