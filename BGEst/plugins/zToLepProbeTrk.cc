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
  // bool isMatchedToElecTriggerObject (const edm::Event &, const edm::TriggerResults &, const pat::Electron &, const std::vector<pat::TriggerObjectStandAlone> &, const std::string &, const std::string &, const double = 0.1);
  // bool isMatchedToMuonTriggerObject (const edm::Event &, const edm::TriggerResults &, const pat::Muon &, const std::vector<pat::TriggerObjectStandAlone> &, const std::string &, const std::string &, const double = 0.1);
  // bool isMatchedToTauTriggerObject (const edm::Event &, const edm::TriggerResults &, const pat::Tau &, const std::vector<pat::TriggerObjectStandAlone> &, const std::string &, const std::string &, const double = 0.1);
  // bool passesDecayModeReconstruction (const pat::Tau &);
  // bool passesLightFlavorRejection (const pat::Tau &);
  // bool inTOBCrack (const pat::IsolatedTrack &);
  // const int extraMissingMiddleHits (const pat::IsolatedTrack &) const;
  // const int hitDrop_missingMiddleHits (const pat::IsolatedTrack &) const;
  // double dRMinJet (const pat::IsolatedTrack &, const std::vector<pat::Jet> &);
  // double deltaRToClosestElectron(const pat::IsolatedTrack &, const std::vector<pat::Electron> &);
  // double deltaRToClosestMuon(const pat::IsolatedTrack &, const std::vector<pat::Muon> &);
  // double deltaRToClosestTauHad(const pat::IsolatedTrack &, const std::vector<pat::Tau> &);
  // double energyGivenMass (const double, const pat::IsolatedTrack &);
  // bool goodInvMassElec (const pat::Electron &, const pat::IsolatedTrack &);
  // bool goodInvMassMuon (const pat::Muon &, const pat::IsolatedTrack &);
  // double deltaRToClosestPFElectron(const pat::IsolatedTrack &, const std::vector<pat::PackedCandidate> &);
  // double deltaRToClosestPFMuon(const pat::IsolatedTrack &, const std::vector<pat::PackedCandidate> &);
  // double deltaRToClosestVetoElectron (const pat::IsolatedTrack &, const std::vector<pat::Electron> &, const reco::Vertex &);
  // double deltaRToClosestLooseMuon (const pat::IsolatedTrack &, const std::vector<pat::Muon> &);
  // bool passesVeto (const pat::IsolatedTrack &, const std::vector<pat::PackedCandidate> &);
  // bool passesLooseElecVeto (const pat::IsolatedTrack &, const std::vector<pat::Electron> &, const reco::Vertex &);
  // bool passesLooseMuonVeto (const pat::IsolatedTrack &, const std::vector<pat::Muon> &);
  // GlobalPoint getPosition(const DetId&);
  // bool insideCone(const pat::IsolatedTrack &, const DetId&);
  // double calculateCaloE (const pat::IsolatedTrack&, const EBRecHitCollection &, const EERecHitCollection &, const HBHERecHitCollection &);
  // double caloNewNoPUDRp5CentralCalo(const pat::IsolatedTrack&, const EBRecHitCollection &, const EERecHitCollection &, const HBHERecHitCollection &, const double);
  // double transvMassElec(const pat::Electron&, const pat::MET&);
  // double transvMassMuon(const pat::Muon&, const pat::MET&);

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
  edm::EDGetTokenT<bool> ecalBadCalibFilterUpdateToken_;
  std::string HLTName_;

  edm::ESHandle<CaloGeometry> caloGeometry;

  int nTPOS = 0;
  int nTPSS = 0;
  int nTPOS_veto = 0;
  int nTPSS_veto = 0;
  int nTPOSLoose_veto = 0;
  int nTPSSLoose_veto = 0;

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
      ecalBadCalibFilterUpdateToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("ecalBadCalibReducedMINIAODFilter"))) {
  //now do what ever initialization is needed

  HLTName_ = iConfig.getParameter<std::string>("HLTName");

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

  trackCuts = {
    "track pt > 30 GeV",
    "track fabs(#eta) < 2.1",
    "track fabs ( eta ) < 0.15 || fabs ( eta ) > 0.35",
    "track fabs ( eta ) < 1.42 || fabs ( eta ) > 1.65",
    "track fabs ( eta ) < 1.55 || fabs ( eta ) > 1.85",
    "track !inTOBCrack",
    "track hitPattern_.numberOfValidPixelHits >= 4",
    "track hitPattern_.numberOfValidHits >= 4",
    "track missingInnerHits == 0",
    "track hitDrop_missingMiddleHits == 0",
    "track ((pfIsolationDR03_.chargedHadronIso + pfIsolationDR03_.puChargedHadronIso) / pt) < 0.05",
    "track |d0| < 0.02","track |dz| < 0.5",
    "track dRMinJet > 0.5",
    "track deltaRToClosestMuon > 0.15",
    "track deltaRToClosestTauHad > 0.15"
  };

  if(strcmp(T, "muon") == 0) trackCuts.push_back("track caloNewNoPUDRp5CentralCalo < 10");

  int nBins = int(commonCuts.size() + trackCuts.size());
  double maxRange = 0.0;
  std::vector<std::string> lepCuts;

  if(strcmp(T, "electron") == 0){nBins = nBins + int(electronCuts.size()); maxRange = double(nBins) - 0.5; lepCuts = electronCuts;}
  if(strcmp(T, "muon") == 0){nBins = nBins + int(muonCuts.size()); maxRange = double(nBins) - 0.5; lepCuts = muonCuts;}
  if(strcmp(T, "taue") == 0){nBins = nBins + int(taueCuts.size()); maxRange = double(nBins) - 0.5; lepCuts = taueCuts;}
  if(strcmp(T, "taum") == 0){nBins = nBins + int(taumCuts.size()); maxRange = double(nBins) - 0.5; lepCuts = taumCuts;}

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

  caloGeometry = iSetup.getHandle(caloGeometryToken_);

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

// -----------------------------------------------------------------------------------------------------------------------------------------------------
// Cut Name                                                                                                       Events     Cumul. Eff.     Indiv. Eff.
// -----------------------------------------------------------------------------------------------------------------------------------------------------
// total                                                                                                         10000.0        100.000%        100.000%
// trigger                                                                                                        3945.0         39.450%         39.450%
//   HLT_Ele32_WPTight_Gsf_v
// MET filter                                                                                                     3890.0         38.900%         99.310%
//  Flag_goodVertices AND
//  Flag_globalTightHalo2016Filter AND
//  Flag_HBHENoiseFilter AND
//  Flag_HBHENoiseIsoFilter AND
//  Flag_EcalDeadCellTriggerPrimitiveFilter
// >= 1 mets with badPFMuonFilter && badChargedCandidateFilter && passecalBadCalibFilterUpdate                    3890.0         38.900%        100.000%
// >= 1 electrons with pt > 32                                                                                    3855.0         38.550%         39.300%
// >= 1 electrons firing trigger                                                                                  3823.0         38.230%         39.170%
//   HLT_Ele32_WPTight_Gsf_v
// >= 1 electrons with fabs(eta) < 2.1                                                                            3725.0         37.250%         53.100%
// >= 1 electrons with passesVID_tightID (ID + iso)                                                               3457.0         34.570%         37.990%
// >= 1 electrons with |d0| < 0.05, 0.10 (EB, EE)                                                                 3454.0         34.540%         52.660%
// >= 1 electrons with |dz| < 0.10, 0.20 (EB, EE)                                                                 3407.0         34.070%         43.780%
// >= 1 electron-mets with transMass (electron, met) < 40                                                         1858.0         18.580%         47.920%
// >= 1 tracks with pt > 30                                                                                       1858.0         18.580%         93.900%
// >= 1 tracks with fabs ( eta ) < 2.1                                                                            1858.0         18.580%         99.980%
// >= 1 tracks with fabs ( eta ) < 1.42 || fabs ( eta ) > 1.65                                                    1803.0         18.030%         99.440%
// >= 1 tracks with fabs ( eta ) < 0.15 || fabs ( eta ) > 0.35                                                    1712.0         17.120%         99.680%
// >= 1 tracks with fabs ( eta ) < 1.55 || fabs ( eta ) > 1.85                                                    1625.0         16.250%         99.150%
// >= 1 tracks with !inTOBCrack                                                                                   1625.0         16.250%        100.000%
// >= 1 tracks with hitPattern_.numberOfValidPixelHits >= 4                                                       1309.0         13.090%         97.730%
// >= 1 tracks with hitPattern_.numberOfValidHits >= 4                                                            1309.0         13.090%         99.970%
// >= 1 tracks with missingInnerHits == 0                                                                         1292.0         12.920%         99.800%
// >= 1 tracks with hitDrop_missingMiddleHits == 0                                                                1278.0         12.780%         99.640%
// >= 1 tracks with  ((pfIsolationDR03_.chargedHadronIso + pfIsolationDR03_.puChargedHadronIso) / pt) < 0.05       201.0          2.010%         26.360%
// >= 1 tracks with |d0| < 0.02                                                                                    201.0          2.010%         99.960%
// >= 1 tracks with |dz| < 0.5                                                                                     201.0          2.010%         99.460%
// >= 1 tracks with deltaRToClosestElectron > 0.15                                                                   1.0          0.010%         93.240%
// >= 1 tracks with deltaRToClosestMuon > 0.15                                                                       1.0          0.010%         89.630%
// -----------------------------------------------------------------------------------------------------------------------------------------------------

  if(strcmp(T, "taue") == 0){

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

      if(helperFunctions::transvMassElec(electron,met) < 40.0)
        {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      if(auxPassCut[cutIdxInc]) taueTags.push_back(electron);
      
      for(int j = 0; j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }
        
    }
  
  }

  if(strcmp(T, "taum") == 0){

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

      if(helperFunctions::transvMassMuon(muon,met) < 40.0)
        {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

      if(auxPassCut[startTauIdx+cutIdxInc-2]) taumTags.push_back(muon);

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
  if(strcmp(T, "taue") == 0) startTrackIdx = int(commonCuts.size() + taueCuts.size()) - 1;
  if(strcmp(T, "taum") == 0) startTrackIdx = int(commonCuts.size() + taumCuts.size()) - 1;

  int getStTrkIdxLep = int(commonCuts.size()) - 1;

  for (const auto& track : *tracks) {

    int cutIdxInc = 0;
      
    if(track.pt() > 30.) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(passCut[startTrackIdx+cutIdxInc-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(fabs(track.eta()) < 2.1) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if((fabs(track.eta()) < 0.15) || (fabs(track.eta()) > 0.35)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if((fabs(track.eta()) < 1.42) || (fabs(track.eta()) > 1.65)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if((fabs(track.eta()) < 1.55) || (fabs(track.eta()) > 1.85)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    // Need to include fiducial cuts here!!!

    if(!helperFunctions::inTOBCrack(track)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

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

    ++cutIdxInc;

    if(helperFunctions::dRMinJet(track, *jets) > 0.5) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(strcmp(T, "electron") == 0 || strcmp(T, "taue") == 0 || strcmp(T, "taum") == 0)
      {
        if(helperFunctions::deltaRToClosestLepton<pat::Muon>(track, *muons) > 0.15)
          {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}
      }

    if(strcmp(T, "muon") == 0 || strcmp(T, "taue") == 0 || strcmp(T, "taum") == 0)
      {
        if(helperFunctions::deltaRToClosestLepton<pat::Electron>(track, *electrons) > 0.15)
          {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}
      }

    ++cutIdxInc;

    if(strcmp(T, "electron") == 0 || strcmp(T, "muon") == 0) 
      {

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

    if(strcmp(T, "taue") == 0){
      for(int j = int(taueCuts.size()); j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }
    }

    if(strcmp(T, "taum") == 0){
      for(int j = int(taumCuts.size()); j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+int(commonCuts.size())-1] = true;
        auxPassCut[j] = false;
      }
    }

  }

  if(strcmp(T, "electron") == 0){
    for(const auto &tag : electronTags){
      for(const auto &probe : trackProbes){
        if(helperFunctions::goodInvMassLepton<pat::Electron>(tag, probe) && (tag.charge()*probe.charge()) < 0.0){
          ++nTPOS;
          ++hist_nTPOS;
          if(helperFunctions::passesVeto<pat::Electron>(probe, *pfCandidates)){
            ++nTPOS_veto;
            ++hist_nTPOS_veto;
          }
          if(helperFunctions::passesLooseElecVeto(probe, *electrons, pv)){
            ++nTPOSLoose_veto;
            ++hist_nTPOSLoose_veto;
          }
        }
        if(helperFunctions::goodInvMassLepton<pat::Electron>(tag, probe) && (tag.charge()*probe.charge()) > 0.0){
          ++nTPSS;
          ++hist_nTPSS;
          if(helperFunctions::passesVeto<pat::Electron>(probe, *pfCandidates)){
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
        if(helperFunctions::goodInvMassLepton<pat::Muon>(tag, probe) && (tag.charge()*probe.charge()) < 0.0){
          ++nTPOS;
          ++hist_nTPOS;
          if(helperFunctions::passesVeto<pat::Muon>(probe, *pfCandidates)){
            ++nTPOS_veto;
            ++hist_nTPOS_veto;
          }
          if(helperFunctions::passesLooseMuonVeto(probe, *muons)){
            ++nTPOSLoose_veto;
            ++hist_nTPOSLoose_veto;
          }
        }
        if(helperFunctions::goodInvMassLepton<pat::Muon>(tag, probe) && (tag.charge()*probe.charge()) > 0.0){
          ++nTPSS;
          ++hist_nTPSS;
          if(helperFunctions::passesVeto<pat::Muon>(probe, *pfCandidates)){
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

  if(strcmp(T, "taue") == 0){
    for(const auto &tag : taueTags){
      for(const auto &probe : trackProbes){
        if(helperFunctions::goodInvMassLepton<pat::Electron>(tag, probe) && (tag.charge()*probe.charge()) < 0.0){
          ++nTPOS;
          ++hist_nTPOS;
          if(helperFunctions::passesVeto<pat::Electron>(probe, *pfCandidates)){
            ++nTPOS_veto;
            ++hist_nTPOS_veto;
          }
        }
        if(helperFunctions::goodInvMassLepton<pat::Electron>(tag, probe) && (tag.charge()*probe.charge()) > 0.0){
          ++nTPSS;
          ++hist_nTPSS;
          if(helperFunctions::passesVeto<pat::Electron>(probe, *pfCandidates)){
            ++nTPSS_veto;
            ++hist_nTPSS_veto;
          }
        }
      }
    }
  }

  if(strcmp(T, "taum") == 0){
    for(const auto &tag : taumTags){
      for(const auto &probe : trackProbes){
        if(helperFunctions::goodInvMassLepton<pat::Muon>(tag, probe) && (tag.charge()*probe.charge()) < 0.0){
          ++nTPOS;
          ++hist_nTPOS;
          if(helperFunctions::passesVeto<pat::Muon>(probe, *pfCandidates)){
            ++nTPOS_veto;
            ++hist_nTPOS_veto;
          }
        }
        if(helperFunctions::goodInvMassLepton<pat::Muon>(tag, probe) && (tag.charge()*probe.charge()) > 0.0){
          ++nTPSS;
          ++hist_nTPSS;
          if(helperFunctions::passesVeto<pat::Muon>(probe, *pfCandidates)){
            ++nTPSS_veto;
            ++hist_nTPSS_veto;
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



// MEMBER FUNCTIONS ===============================================================



// template<char const *T>
// bool zToLepProbeTrk<T>::isMatchedToElecTriggerObject (const edm::Event &event, const edm::TriggerResults &triggers, const pat::Electron &obj, const std::vector<pat::TriggerObjectStandAlone> &trigObjs, const std::string &collection, const std::string &filter, const double dR)
// {
//   if(collection == "") return false;
//   for(auto trigObj : trigObjs) {
//     trigObj.unpackNamesAndLabels(event, triggers);
//     if(trigObj.collection() != collection) continue;
//     if(filter != "") {
//       bool flag = false;
//       for(const auto &filterLabel : trigObj.filterLabels ())
//         if(filterLabel == filter) {
//           flag = true;
//           break;
//         }
//       if (!flag) continue;
//     }
//     if(deltaR (obj, trigObj) > dR) continue;
//     return true;
//   }
//   return false;
// }

// template<char const *T>
// bool zToLepProbeTrk<T>::isMatchedToMuonTriggerObject (const edm::Event &event, const edm::TriggerResults &triggers, const pat::Muon &obj, const std::vector<pat::TriggerObjectStandAlone> &trigObjs, const std::string &collection, const std::string &filter, const double dR)
// {
//   if(collection == "") return false;
//   for(auto trigObj : trigObjs) {
//     trigObj.unpackNamesAndLabels(event, triggers);
//     if(trigObj.collection() != collection) continue;
//     if(filter != "") {
//       bool flag = false;
//       for(const auto &filterLabel : trigObj.filterLabels ())
//         if(filterLabel == filter) {
//           flag = true;
//           break;
//         }
//       if (!flag) continue;
//     }
//     if(deltaR (obj, trigObj) > dR) continue;
//     return true;
//   }
//   return false;
// }

// template<char const *T>
// bool zToLepProbeTrk<T>::isMatchedToTauTriggerObject (const edm::Event &event, const edm::TriggerResults &triggers, const pat::Tau &obj, const std::vector<pat::TriggerObjectStandAlone> &trigObjs, const std::string &collection, const std::string &filter, const double dR)
// {
//   if(collection == "") return false;
//   for(auto trigObj : trigObjs) {
//     trigObj.unpackNamesAndLabels(event, triggers);
//     if(trigObj.collection() != collection) continue;
//     if(filter != "") {
//       bool flag = false;
//       for(const auto &filterLabel : trigObj.filterLabels ())
//         if(filterLabel == filter) {
//           flag = true;
//           break;
//         }
//       if (!flag) continue;
//     }
//     if(deltaR (obj, trigObj) > dR) continue;
//     return true;
//   }
//   return false;
// }

// template<char const *T>
// bool zToLepProbeTrk<T>::passesDecayModeReconstruction (const pat::Tau &tau){
//   return (tau.tauID("decayModeFindingNewDMs") && (tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") || (tau.tauID("chargedIsoPtSumdR03")+std::max(0.,tau.tauID("neutralIsoPtSumdR03")-0.072*tau.tauID("puCorrPtSum"))<2.5) || tau.tauID("byVVVLooseDeepTau2017v2p1VSjet") || tau.tauID("byVVVLooseDeepTau2018v2p5VSjet")));
// }

// template<char const *T>
// bool zToLepProbeTrk<T>::passesLightFlavorRejection (const pat::Tau &tau){
//   return (tau.tauID("byVVVLooseDeepTau2017v2p1VSe") || tau.tauID("byVVVLooseDeepTau2018v2p5VSe")) && (tau.tauID("byVLooseDeepTau2017v2p1VSmu") || tau.tauID("byVLooseDeepTau2018v2p5VSmu"));
// }

// template<char const *T>
// bool zToLepProbeTrk<T>::inTOBCrack (const pat::IsolatedTrack &track){
//   return (fabs(track.dz()) < 0.5 && fabs(1.57079632679489661923 - track.theta()) < 1.0e-3);
// }

// template<char const *T> 
// const int zToLepProbeTrk<T>::extraMissingMiddleHits (const pat::IsolatedTrack &track) const
// {

//   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//   std::default_random_engine generator (seed);
//   std::uniform_real_distribution<double> distribution (0.0, 1.0);
//   double dropTOBProbability_ = 0.00830971251971; // This is taken from Run 2; needs to be updated once correct value is estimated
//   double hitProbability_ = 0.0175874821487; // This is taken from Run 2; needs to be updated once correct value is estimated
//   // bool dropHits = true;
//   bool dropHits = false;
//   bool dropTOBDecision_ = (dropHits ? distribution (generator) : 1.0e6) < dropTOBProbability_;
//   std::vector<bool> dropMiddleHitDecisions_;
//   for (int i = 0; i < 50; i++)
//     dropMiddleHitDecisions_.push_back ((dropHits ? distribution (generator) : 1.0e6) < hitProbability_);

//   int nHits = 0;
//   bool countMissingMiddleHits = false;
//   for (int i = 0; i < track.hitPattern().stripLayersWithMeasurement() - (dropTOBDecision_ ? track.hitPattern ().stripTOBLayersWithMeasurement () : 0); i++)
//     {
//       bool hit = !dropMiddleHitDecisions_.at(i);
//       if (!hit && countMissingMiddleHits)
//         nHits++;
//       if (hit)
//         countMissingMiddleHits = true;
//     }

//   return nHits;
// }

// template<char const *T>
// const int zToLepProbeTrk<T>::hitDrop_missingMiddleHits (const pat::IsolatedTrack &track) const
// {
//   int nDropHits = helperFunctions::extraMissingMiddleHits(track);
//   return track.hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS) + nDropHits;
// }

// template<char const *T>
// double zToLepProbeTrk<T>::dRMinJet (const pat::IsolatedTrack &track, const std::vector<pat::Jet> &jets)
// {
//   double dRMinJet = 99.0;
//   for(const auto &jet : jets){
//     if(jet.pt() > 30 &&
//     fabs(jet.eta()) < 4.5 &&
//     (((jet.neutralHadronEnergyFraction()<0.90 &&
//     jet.neutralEmEnergyFraction()<0.90 &&
//     (jet.chargedMultiplicity() + jet.neutralMultiplicity())>1 &&
//     jet.muonEnergyFraction()<0.8) &&
//     ((fabs(jet.eta())<=2.4 &&
//     jet.chargedHadronEnergyFraction()>0 &&
//     jet.chargedMultiplicity()>0 &&
//     jet.chargedEmEnergyFraction()<0.90) ||
//     fabs(jet.eta())>2.4) &&
//     fabs(jet.eta())<=3.0) ||
//     (jet.neutralEmEnergyFraction()<0.90 && jet.neutralMultiplicity()>10 && fabs(jet.eta())>3.0)))
//     {
//       double dR = deltaR(track, jet);
//       if(dR < dRMinJet || dRMinJet < 0.0) dRMinJet = dR;
//     }
//   }
//   return dRMinJet;
// }

// template<char const *T>
// double zToLepProbeTrk<T>::deltaRToClosestElectron(const pat::IsolatedTrack &track, const std::vector<pat::Electron> &electrons) 
// {

//   double dR;
//   double deltaRToClosestElectron = 99.0;

//   for(const auto &electron : electrons) {
//     dR = deltaR(track, electron);
//     if(dR < deltaRToClosestElectron || deltaRToClosestElectron < 0.0) deltaRToClosestElectron = dR;
//   }

//   return deltaRToClosestElectron;

// }

// template<char const *T>
// double zToLepProbeTrk<T>::deltaRToClosestMuon(const pat::IsolatedTrack &track, const std::vector<pat::Muon> &muons) 
// {

//   double dR;
//   double deltaRToClosestMuon = 99.0;

//   for(const auto &muon : muons) {
//     dR = deltaR(track, muon);
//     if(dR < deltaRToClosestMuon || deltaRToClosestMuon < 0.0) deltaRToClosestMuon = dR;
//   }

//   return deltaRToClosestMuon;

// }

// template<char const *T>
// double zToLepProbeTrk<T>::deltaRToClosestTauHad(const pat::IsolatedTrack &track, const std::vector<pat::Tau> &taus) 
// {

//   double dR;
//   double deltaRToClosestTauHad = 99.0;
//   bool passesDecayModeReconstruction;
//   bool passesLightFlavorRejection;

//   for(const auto &tau : taus) {
//     dR = deltaR(track, tau);

//     passesDecayModeReconstruction = (tau.tauID("decayModeFindingNewDMs"));
//     passesLightFlavorRejection = (tau.tauID("byVVVLooseDeepTau2017v2p1VSe") || tau.tauID("byVVVLooseDeepTau2018v2p5VSe")) && (tau.tauID("byVLooseDeepTau2017v2p1VSmu") || tau.tauID("byVLooseDeepTau2018v2p5VSmu"));

//     if(passesDecayModeReconstruction && passesLightFlavorRejection && (dR < deltaRToClosestTauHad  || deltaRToClosestTauHad  < 0.0)) {
//       deltaRToClosestTauHad = dR;
//     }
//   }

//   return deltaRToClosestTauHad;
// }

// template<char const *T>
// double zToLepProbeTrk<T>::energyGivenMass (const double mass, const pat::IsolatedTrack &track)
// {
//   return sqrt (track.px () * track.px () + track.py () * track.py () + track.pz () * track.pz () + mass * mass);
// }

// template<char const *T> 
// bool zToLepProbeTrk<T>::goodInvMassElec (const pat::Electron &tag, const pat::IsolatedTrack &probe)
// {
//   TLorentzVector t (tag.px(), tag.py(), tag.pz(), tag.energy()),
//                  p (probe.px(), probe.py(), probe.pz(), helperFunctions::energyGivenMass(0.000510998950, probe)); // Electron mass extracted from PDG on 27/06/2024 https://pdg.lbl.gov/2024/tables/contents_tables.html
//   double m = (t + p).M();
//   return (fabs (m - 91.1880) < 10.0); // Z mass extracted from PDG on 27/06/2024 https://pdg.lbl.gov/2024/tables/contents_tables.html
// }

// template<char const *T> 
// bool zToLepProbeTrk<T>::goodInvMassMuon (const pat::Muon &tag, const pat::IsolatedTrack &probe)
// {
//   TLorentzVector t (tag.px(), tag.py(), tag.pz(), tag.energy()),
//                  p (probe.px(), probe.py(), probe.pz(), helperFunctions::energyGivenMass(0.1056583755, probe)); // Muon mass extracted from PDG on 27/06/2024 https://pdg.lbl.gov/2024/tables/contents_tables.html
//   double m = (t + p).M();
//   return (fabs (m - 91.1880) < 10.0); // Z mass extracted from PDG on 27/06/2024 https://pdg.lbl.gov/2024/tables/contents_tables.html
// }

// template<char const *T>
// double zToLepProbeTrk<T>::deltaRToClosestPFElectron(const pat::IsolatedTrack &track, const std::vector<pat::PackedCandidate> &pfCandidates) 
// {
//   double deltaRToClosestPFElectron = 99.0;
//   for(const auto &pfCandidate : pfCandidates) {
//       int pdgid = abs(pfCandidate.pdgId());
//       if(pdgid != 11) continue;

//       double dR = deltaR(track, pfCandidate);

//       if(pdgid == 11 &&
//          (dR < deltaRToClosestPFElectron || deltaRToClosestPFElectron < 0.0))
//         deltaRToClosestPFElectron = dR;
//   }
//   return deltaRToClosestPFElectron;
// }

// template<char const *T>
// double zToLepProbeTrk<T>::deltaRToClosestPFMuon(const pat::IsolatedTrack &track, const std::vector<pat::PackedCandidate> &pfCandidates) 
// {
//   double deltaRToClosestPFMuon = 99.0;
//   for(const auto &pfCandidate : pfCandidates) {
//       int pdgid = abs(pfCandidate.pdgId());
//       if(pdgid != 13) continue;

//       double dR = deltaR(track, pfCandidate);

//       if(pdgid == 11 &&
//          (dR < deltaRToClosestPFMuon || deltaRToClosestPFMuon < 0.0))
//         deltaRToClosestPFMuon = dR;
//   }
//   return deltaRToClosestPFMuon;
// }

// template<char const *T>
// double zToLepProbeTrk<T>::deltaRToClosestVetoElectron (const pat::IsolatedTrack &track, const std::vector<pat::Electron> &electrons, const reco::Vertex &vertex)
// {
//   double deltaRToClosestVetoElectron = 99.0;
  
//   double dR;

//   for(const auto &electron : electrons) {
//     dR = deltaR(track, electron);

//     bool passesVeto_dxy = false, passesVeto_dz = false;

//     // Note in below, these remain false if |eta| >= 2.5; thus an eta cut is also being applied here as intended
//     double ele_d0 = fabs(electron.gsfTrack()->dxy(vertex.position()));
//     double ele_dz = fabs(electron.gsfTrack()->dz(vertex.position()));

//     if(fabs(electron.superCluster ()->eta()) <= 1.479) {
//       passesVeto_dxy = (ele_d0 < 0.05);
//       passesVeto_dz = (ele_dz < 0.10);
//     }
//     else if(fabs(electron.superCluster()->eta()) < 2.5) {
//       passesVeto_dxy = (ele_d0 < 0.10);
//       passesVeto_dz = (ele_dz < 0.20);
//     }

//     if(electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-tight") &&
//        passesVeto_dxy &&
//        passesVeto_dz &&
//        (dR < deltaRToClosestVetoElectron || deltaRToClosestVetoElectron < 0.0)) {
//       deltaRToClosestVetoElectron = dR;
//     }
//   } // for electrons

//   return deltaRToClosestVetoElectron;
// }

// template<char const *T>
// double zToLepProbeTrk<T>::deltaRToClosestLooseMuon(const pat::IsolatedTrack &track, const std::vector<pat::Muon> &muons) 
// {
//   double deltaRToClosestLooseMuon = 99.0;

//   double dR;

//   for(const auto &muon : muons) {
//     dR = deltaR(track, muon);
//     if(muon.isLooseMuon()  && (dR < deltaRToClosestLooseMuon  || deltaRToClosestLooseMuon  < 0.0)) deltaRToClosestLooseMuon = dR;
//   }

//   return deltaRToClosestLooseMuon;
// }

// template<char const *T> 
// bool zToLepProbeTrk<T>::passesVeto (const pat::IsolatedTrack &probe, const std::vector<pat::PackedCandidate> &pfCandidates)
// {
//   bool passesElec = helperFunctions::deltaRToClosestPFElectron(probe, pfCandidates) > 0.15
//              && (probe.matchedCaloJetEmEnergy() + probe.matchedCaloJetHadEnergy()) < 10.0//;
//              // && probe.hitAndTOBDrop_bestTrackMissingOuterHits () >= 3.0; // This is not applied for BG MC
//              && probe.lostOuterLayers() >= 3.0; // This is not applied for BG MC
//   bool passesMuon = helperFunctions::deltaRToClosestPFMuon(probe, pfCandidates) > 0.15
//              // && probe.hitAndTOBDrop_bestTrackMissingOuterHits () >= 3.0; // This is not applied for BG MC
//              && probe.lostOuterLayers() >= 3.0; // This is not applied for BG MC

//   // if(std::is_same<T, pat::Electron>::value) return passesElec;
//   if(strcmp(T, "electron") == 0) return passesElec;
//   // if(std::is_same<T, pat::Muon>::value) return passesMuon;
//   if(strcmp(T, "muon") == 0) return passesMuon;
//   // if(std::is_same<T, pat::Tau>::value) return passesTau;
// }

// template<char const *T> 
// bool zToLepProbeTrk<T>::passesLooseElecVeto (const pat::IsolatedTrack &probe, const std::vector<pat::Electron> &electrons, const reco::Vertex &vertex)
// {
//   bool passes = helperFunctions::deltaRToClosestVetoElectron(probe,electrons,vertex) > 0.15
//              && (probe.matchedCaloJetEmEnergy() + probe.matchedCaloJetHadEnergy()) < 10.0
//              // && probe.hitAndTOBDrop_bestTrackMissingOuterHits () >= 3.0; // This is not applied for BG MC
//              && probe.lostOuterLayers() >= 3.0; // This is not applied for BG MC
//   return passes;
// }

// template<char const *T> 
// bool zToLepProbeTrk<T>::passesLooseMuonVeto (const pat::IsolatedTrack &probe, const std::vector<pat::Muon> &muons)
// {
//   bool passes = helperFunctions::deltaRToClosestLooseMuon(probe, muons) > 0.15
//              // && probe.hitAndTOBDrop_bestTrackMissingOuterHits () >= 3.0; // This is not applied for BG MC
//              && probe.lostOuterLayers() >= 3.0; // This is not applied for BG MC          
//   return passes;
// }

// template<char const *T>
// GlobalPoint zToLepProbeTrk<T>::getPosition(const DetId& id)
// {
//    if ( ! caloGeometry.isValid() ||
//         ! caloGeometry->getSubdetectorGeometry(id) ||
//         ! caloGeometry->getSubdetectorGeometry(id)->getGeometry(id) ) {
//       throw cms::Exception("FatalError") << "Failed to access geometry for DetId: " << id.rawId();
//       return GlobalPoint(0,0,0);
//    }
//    return caloGeometry->getSubdetectorGeometry(id)->getGeometry(id)->getPosition();
// }

// template<char const *T> 
// bool zToLepProbeTrk<T>::insideCone(const pat::IsolatedTrack &candTrack, const DetId& id)
// {
//    GlobalPoint idPosition = helperFunctions::getPosition(id,caloGeometry);
//    if (idPosition.mag()<0.01) return false;
//    math::XYZVector idPositionRoot( idPosition.x(), idPosition.y(), idPosition.z() );
//    return deltaR(candTrack, idPositionRoot) < 0.5;
// }

// template<char const *T> 
// double zToLepProbeTrk<T>::calculateCaloE(const pat::IsolatedTrack& track, const EBRecHitCollection &EBRecHits, const EERecHitCollection &EERecHits, const HBHERecHitCollection &HBHERecHits)
// { 

//   double caloEnergy = 0.0;

//   for (const auto &hit : EBRecHits) {
//     if (helperFunctions::insideCone(track, hit.detid(),caloGeometry)) {
//       caloEnergy += hit.energy();
//     }
//   }
//   for (const auto &hit : EERecHits) {
//     if (helperFunctions::insideCone(track, hit.detid(),caloGeometry)) {
//       caloEnergy += hit.energy();
//     }
//   }

//   for (const auto &hit : HBHERecHits) {
//     if (helperFunctions::insideCone(track, hit.detid(),caloGeometry)) {
//       caloEnergy += hit.energy();
//     }
//   }

//   return caloEnergy;
// }

// template<char const *T> 
// double zToLepProbeTrk<T>::caloNewNoPUDRp5CentralCalo(const pat::IsolatedTrack& track, const EBRecHitCollection &EBRecHits, const EERecHitCollection &EERecHits, const HBHERecHitCollection &HBHERecHits, const double rhoCentralCalo)
// {
  
//   double rawCaloTot = helperFunctions::calculateCaloE(track, EBRecHits, EERecHits, HBHERecHits, caloGeometry);
//   double caloCorr = rhoCentralCalo * TMath::Pi() * 0.5 * 0.5;  // Define effective area as pi*r^2, where r is radius of DeltaR cone.
//   double caloNewNoPUDRp5CentralCalo = TMath::Max(0., rawCaloTot - caloCorr);

//   return caloNewNoPUDRp5CentralCalo;
// }

// template<char const *T>
// void zToLepProbeTrk<T>::envSet (const edm::EventSetup& iSetup)
// {
//   caloGeometry = iSetup.getHandle(caloGeometryToken_);
//   if(!caloGeometry.isValid()) throw cms::Exception("FatalError") << "Failed to get the caloGeometry!";
// }

// template<char const *T> 
// double zToLepProbeTrk<T>::transvMassElec(const pat::Electron& electron, const pat::MET& met)
// {
  
//   double rawCaloTot = electron.pt() - met.pt();

//   return rawCaloTot;
// }

// template<char const *T> 
// double zToLepProbeTrk<T>::transvMassMuon(const pat::Muon& muon, const pat::MET& met)
// {
  
//   double rawCaloTot = muon.pt() - met.pt();

//   return rawCaloTot;
// }


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
