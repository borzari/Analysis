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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "Analysis/Helper/interface/helperFunctions.h"
#include "Analysis/Helper/interface/plotPrintFunctions.h"
#include "Analysis/Helper/interface/selectingFunctions.h"
#include "Analysis/Helper/interface/sfFunctions.h"

//
// class declaration
//

// If the filter does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDFilter<>
// This will improve performance in multithreaded jobs.

template<char const *T>
class plotterZToLepProbeTrk : public edm::one::EDFilter<edm::one::SharedResources> {
public:
  explicit plotterZToLepProbeTrk(const edm::ParameterSet&);
  ~plotterZToLepProbeTrk() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  std::vector<std::string> commonCuts, electronCuts, muonCuts, taueCuts, taumCuts, trackCuts;

private:
  bool filter(edm::Event&, const edm::EventSetup&) override;
  edm::Service<TFileService> fs_;
  std::map<std::string, TH1D *> oneDHists_;
  std::map<std::string, TH2D *> twoDHists_;

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
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken_;
  edm::EDGetTokenT<std::vector<pat::MET>> metsToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggersPATToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggersHLTToken_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> trigobjsToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  edm::ESGetToken<EcalChannelStatus, EcalChannelStatusRcd> ecalStatusToken_;
  edm::EDGetTokenT<bool> ecalBadCalibFilterUpdateToken_;
  std::string HLTName_;
  bool isCRAB_;

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
plotterZToLepProbeTrk<T>::plotterZToLepProbeTrk(const edm::ParameterSet& iConfig)
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
      generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
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
  isCRAB_ = iConfig.getParameter<bool>("isCRAB");

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

  std::string elecStr = "electron";
  std::string muonStr = "muon";
  std::string trackStr = "track";

  TFileDirectory NLayers4 = fs_->mkdir("NLayers4");
  TFileDirectory NLayers5 = fs_->mkdir("NLayers5");
  TFileDirectory NLayers6Plus = fs_->mkdir("NLayers6Plus");

  plotPrintFunctions::createCommonHists(fs_,oneDHists_,twoDHists_,NLayers4,4);
  plotPrintFunctions::createCommonHists(fs_,oneDHists_,twoDHists_,NLayers5,5);
  plotPrintFunctions::createCommonHists(fs_,oneDHists_,twoDHists_,NLayers6Plus,6);
  if(strcmp(T, "electron") == 0) {
    plotPrintFunctions::createObjHists<pat::Electron>(fs_,oneDHists_,elecStr,NLayers4,4);
    plotPrintFunctions::createObjHists<pat::Electron>(fs_,oneDHists_,elecStr,NLayers5,5);
    plotPrintFunctions::createObjHists<pat::Electron>(fs_,oneDHists_,elecStr,NLayers6Plus,6);
  }
  if(strcmp(T, "muon") == 0) {
    plotPrintFunctions::createObjHists<pat::Muon>(fs_,oneDHists_,muonStr,NLayers4,4);
    plotPrintFunctions::createObjHists<pat::Muon>(fs_,oneDHists_,muonStr,NLayers5,5);
    plotPrintFunctions::createObjHists<pat::Muon>(fs_,oneDHists_,muonStr,NLayers6Plus,6);
  }
  if(strcmp(T, "tauele") == 0) {
    plotPrintFunctions::createObjHists<pat::Electron>(fs_,oneDHists_,elecStr,NLayers4,4);
    plotPrintFunctions::createObjHists<pat::Electron>(fs_,oneDHists_,elecStr,NLayers5,5);
    plotPrintFunctions::createObjHists<pat::Electron>(fs_,oneDHists_,elecStr,NLayers6Plus,6);
  }
  if(strcmp(T, "taumu") == 0) {
    plotPrintFunctions::createObjHists<pat::Muon>(fs_,oneDHists_,muonStr,NLayers4,4);
    plotPrintFunctions::createObjHists<pat::Muon>(fs_,oneDHists_,muonStr,NLayers5,5);
    plotPrintFunctions::createObjHists<pat::Muon>(fs_,oneDHists_,muonStr,NLayers6Plus,6);
  }
  plotPrintFunctions::createObjHists<pat::IsolatedTrack>(fs_,oneDHists_,trackStr,NLayers4,4);
  plotPrintFunctions::createObjHists<pat::IsolatedTrack>(fs_,oneDHists_,trackStr,NLayers5,5);
  plotPrintFunctions::createObjHists<pat::IsolatedTrack>(fs_,oneDHists_,trackStr,NLayers6Plus,6);
  plotPrintFunctions::createTPPairsHists(fs_,oneDHists_,NLayers4,4);
  plotPrintFunctions::createTPPairsHists(fs_,oneDHists_,NLayers5,5);
  plotPrintFunctions::createTPPairsHists(fs_,oneDHists_,NLayers6Plus,6);

}

template<char const *T>
plotterZToLepProbeTrk<T>::~plotterZToLepProbeTrk() {

  std::cout << "# OS T&P pairs before veto: " << nTPOS << std::endl;
  std::cout << "# SS T&P pairs before veto: " << nTPSS << std::endl;
  std::cout << "# OS T&P pairs after veto: " << nTPOS_veto << std::endl;
  std::cout << "# SS T&P pairs after veto: " << nTPSS_veto << std::endl;
  std::cout << "# OS T&P pairs after loose veto: " << nTPOSLoose_veto << std::endl;
  std::cout << "# SS T&P pairs after loose veto: " << nTPSSLoose_veto << std::endl;

}

// ------------ method called for each event  ------------
template<char const *T>
bool plotterZToLepProbeTrk<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  std::vector<bool> passSel;
  std::vector<bool> passCut;

  bool isGood = false;

  VariablesToPlot variables;

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

  TVector2 metNoMu = helperFunctions::calcMetNoMu(met, pfCandidates);

  variables.met = met;
  variables.metNoMu = metNoMu;

  edm::Handle<GenEventInfoProduct> generator;
  iEvent.getByToken(generatorToken_, generator);

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

  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(verticesToken_, vertices);
  const reco::Vertex &pv = vertices->at(0);

  std::vector<bool> auxPassCut;

  // This is just for these vectors to not be empty, for the selections to work without further modifications
  for(int i = 0; i < 50; ++i){
    passSel.push_back(true);
    passCut.push_back(true);
    auxPassCut.push_back(false);
  }

  if(strcmp(T, "electron") == 0){

    int startElecIdx = 1;
    
    for (const auto& electron : *electrons) {
      
      selectingFunctions::zToLepElecSel(passSel,passCut,auxPassCut,electron,startElecIdx,pv,electronTags,iEvent,triggerBitsHLT,triggerObjs);
      
      for(int j = 0; j < int(auxPassCut.size()); ++j) auxPassCut[j] = false;

    }

    variables.electrons = electronTags;
  
  }

  if(strcmp(T, "muon") == 0){
    
    int startMuonIdx = 1;

    for (const auto& muon : *muons) {

      selectingFunctions::zToLepMuonSel(passSel,passCut,auxPassCut,muon,startMuonIdx,pv,muonTags,iEvent,triggerBitsHLT,triggerObjs);

      for(int j = 0; j < int(auxPassCut.size()); ++j) auxPassCut[j] = false;

    }

    variables.muons = muonTags;
  
  }

  if(strcmp(T, "tauele") == 0){

    int startTauIdx = 1;
    
    for (const auto& electron : *electrons) {
      
      selectingFunctions::zToLepTauEleSel(passSel,passCut,auxPassCut,electron,startTauIdx,pv,taueTags,iEvent,triggerBitsHLT,triggerObjs,met);
      
      for(int j = 0; j < int(auxPassCut.size()); ++j) auxPassCut[j] = false;
        
    }

    variables.electrons = taueTags;
  
  }

  if(strcmp(T, "taumu") == 0){

    int startTauIdx = 1;
    
    for (const auto& muon : *muons) {

      selectingFunctions::zToLepTauMuSel(passSel,passCut,auxPassCut,muon,startTauIdx,pv,taumTags,iEvent,triggerBitsHLT,triggerObjs,met);

      for(int j = 0; j < int(auxPassCut.size()); ++j) auxPassCut[j] = false;

    }

    variables.muons = taumTags;
  
  }

  int startTrackIdx = 1;

  int getStTrkIdxLep = 0;

  for (const auto& track : *tracks) {

    if(strcmp(T, "electron") == 0) selectingFunctions::zToLepTrackSel<pat::Electron>(passSel, passCut, auxPassCut, track, startTrackIdx, getStTrkIdxLep, trackProbes, vetoListElec, vetoListMu, EcalAllDeadChannelsValMap, EcalAllDeadChannelsBitMap, jets, electrons, muons, taus, EBRecHits, EERecHits, HBHERecHits, rhoCentralCalo, caloGeometry);

    if(strcmp(T, "muon") == 0) selectingFunctions::zToLepTrackSel<pat::Muon>(passSel, passCut, auxPassCut, track, startTrackIdx, getStTrkIdxLep, trackProbes, vetoListElec, vetoListMu, EcalAllDeadChannelsValMap, EcalAllDeadChannelsBitMap, jets, electrons, muons, taus, EBRecHits, EERecHits, HBHERecHits, rhoCentralCalo, caloGeometry);

    if(strcmp(T, "tauele") == 0 || strcmp(T, "taumu") == 0) selectingFunctions::zToLepTrackSel<pat::Tau>(passSel, passCut, auxPassCut, track, startTrackIdx, getStTrkIdxLep, trackProbes, vetoListElec, vetoListMu, EcalAllDeadChannelsValMap, EcalAllDeadChannelsBitMap, jets, electrons, muons, taus, EBRecHits, EERecHits, HBHERecHits, rhoCentralCalo, caloGeometry);

    for(int j = 0; j < int(auxPassCut.size()); ++j) auxPassCut[j] = false;

  }

  variables.tracks = trackProbes;

  for(const auto &track : trackProbes){
    variables.caloNewNoPUDRp5CentralCalo.push_back(helperFunctions::caloNewNoPUDRp5CentralCalo(track, *EBRecHits, *EERecHits, *HBHERecHits, *rhoCentralCalo, caloGeometry));
  }

  if(strcmp(T, "electron") == 0){
    for(const auto &tag : electronTags){
      for(const auto &probe : trackProbes){
        if(helperFunctions::goodInvMassLepton<pat::Electron>(tag, probe, false) && (tag.charge()*probe.charge()) < 0.0){
          ++nTPOS;
          ++hist_nTPOS;
          double caloNewNoPUDRp5CentralCalo = helperFunctions::caloNewNoPUDRp5CentralCalo(probe, *EBRecHits, *EERecHits, *HBHERecHits, *rhoCentralCalo, caloGeometry);
          if(helperFunctions::passesVeto<pat::Electron>(probe, *pfCandidates, *jets, caloNewNoPUDRp5CentralCalo)){
            ++nTPOS_veto;
            ++hist_nTPOS_veto;
          }
          if(helperFunctions::passesLooseElecVeto(probe, *electrons, pv, caloNewNoPUDRp5CentralCalo)){
            ++nTPOSLoose_veto;
            ++hist_nTPOSLoose_veto;
          }
        }
        if(helperFunctions::goodInvMassLepton<pat::Electron>(tag, probe, false) && (tag.charge()*probe.charge()) > 0.0){
          ++nTPSS;
          ++hist_nTPSS;
          double caloNewNoPUDRp5CentralCalo = helperFunctions::caloNewNoPUDRp5CentralCalo(probe, *EBRecHits, *EERecHits, *HBHERecHits, *rhoCentralCalo, caloGeometry);
          if(helperFunctions::passesVeto<pat::Electron>(probe, *pfCandidates, *jets, caloNewNoPUDRp5CentralCalo)){
            ++nTPSS_veto;
            ++hist_nTPSS_veto;
          }
          if(helperFunctions::passesLooseElecVeto(probe, *electrons, pv, caloNewNoPUDRp5CentralCalo)){
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
          double caloNewNoPUDRp5CentralCalo = helperFunctions::caloNewNoPUDRp5CentralCalo(probe, *EBRecHits, *EERecHits, *HBHERecHits, *rhoCentralCalo, caloGeometry);
          if(helperFunctions::passesVeto<pat::Muon>(probe, *pfCandidates, *jets, caloNewNoPUDRp5CentralCalo)){
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
          double caloNewNoPUDRp5CentralCalo = helperFunctions::caloNewNoPUDRp5CentralCalo(probe, *EBRecHits, *EERecHits, *HBHERecHits, *rhoCentralCalo, caloGeometry);
          if(helperFunctions::passesVeto<pat::Muon>(probe, *pfCandidates, *jets, caloNewNoPUDRp5CentralCalo)){
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
          double caloNewNoPUDRp5CentralCalo = helperFunctions::caloNewNoPUDRp5CentralCalo(probe, *EBRecHits, *EERecHits, *HBHERecHits, *rhoCentralCalo, caloGeometry);
          if(helperFunctions::passesVeto<pat::Tau>(probe, *pfCandidates, *jets, caloNewNoPUDRp5CentralCalo)){
            ++nTPOS_veto;
            ++hist_nTPOS_veto;
          }
          if(helperFunctions::passesLooseElecVeto(probe, *electrons, pv, caloNewNoPUDRp5CentralCalo)){
            ++nTPOSLoose_veto;
            ++hist_nTPOSLoose_veto;
          }
        }
        if(helperFunctions::goodInvMassLepton<pat::Electron>(tag, probe, true) && (tag.charge()*probe.charge()) > 0.0){
          ++nTPSS;
          ++hist_nTPSS;
          double caloNewNoPUDRp5CentralCalo = helperFunctions::caloNewNoPUDRp5CentralCalo(probe, *EBRecHits, *EERecHits, *HBHERecHits, *rhoCentralCalo, caloGeometry);
          if(helperFunctions::passesVeto<pat::Tau>(probe, *pfCandidates, *jets, caloNewNoPUDRp5CentralCalo)){
            ++nTPSS_veto;
            ++hist_nTPSS_veto;
          }
          if(helperFunctions::passesLooseElecVeto(probe, *electrons, pv, caloNewNoPUDRp5CentralCalo)){
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
          double caloNewNoPUDRp5CentralCalo = helperFunctions::caloNewNoPUDRp5CentralCalo(probe, *EBRecHits, *EERecHits, *HBHERecHits, *rhoCentralCalo, caloGeometry);
          if(helperFunctions::passesVeto<pat::Tau>(probe, *pfCandidates, *jets, caloNewNoPUDRp5CentralCalo)){
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
          double caloNewNoPUDRp5CentralCalo = helperFunctions::caloNewNoPUDRp5CentralCalo(probe, *EBRecHits, *EERecHits, *HBHERecHits, *rhoCentralCalo, caloGeometry);
          if(helperFunctions::passesVeto<pat::Tau>(probe, *pfCandidates, *jets, caloNewNoPUDRp5CentralCalo)){
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

  variables.nTP.nTPOS = hist_nTPOS;
  variables.nTP.nTPSS = hist_nTPSS;
  variables.nTP.nTPOS_veto = hist_nTPOS_veto;
  variables.nTP.nTPSS_veto = hist_nTPSS_veto;
  variables.nTP.nTPOS_loose_veto = hist_nTPOSLoose_veto;
  variables.nTP.nTPSS_loose_veto = hist_nTPSSLoose_veto;

  double sf = 1.0;
  std::string elecInputFile = "/home/brenoorzari/CMSSW_13_0_13/src/Analysis/Helper/data/electronSFs.root";
  std::string muonInputFile = "/home/brenoorzari/CMSSW_13_0_13/src/Analysis/Helper/data/muonSFs.root";
  std::string elecReco = "electronReco2022;1";
  std::string elecID = "electronID2022Tight;1";
  std::string muonID = "muonID2022Tight;1";
  std::string muonIso = "muonIso2022TightTightID;1";
  std::string muonTrigger = "muonTrigger2022IsoMu24;1";

  sf *= (*generator).weight() / fabs ((*generator).weight ());

  if(strcmp(T, "electron") == 0 || strcmp(T, "tauele") == 0){
    sfFunctions::sfElectrons(sf, elecInputFile, elecReco, selElectrons);
    sfFunctions::sfElectrons(sf, elecInputFile, elecID, selElectrons);
  }
  if(strcmp(T, "muon") == 0 || strcmp(T, "taumu") == 0){
    sfFunctions::sfMuons(sf, muonInputFile, muonID, selMuons);
    sfFunctions::sfMuons(sf, muonInputFile, muonIso, selMuons);
    sfFunctions::sfMuons(sf, muonInputFile, muonTrigger, selMuons);
  }

  plotPrintFunctions::plotVariables(oneDHists_, twoDHists_, variables, 4, sf);
  plotPrintFunctions::plotVariables(oneDHists_, twoDHists_, variables, 5, sf);
  plotPrintFunctions::plotVariables(oneDHists_, twoDHists_, variables, 6, sf);
  plotPrintFunctions::plotTPPairs(oneDHists_, variables, 4, sf);
  plotPrintFunctions::plotTPPairs(oneDHists_, variables, 5, sf);
  plotPrintFunctions::plotTPPairs(oneDHists_, variables, 6, sf);

  return isGood;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<char const *T>
void plotterZToLepProbeTrk<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

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
  desc.add<edm::InputTag>("generator", edm::InputTag("generator"));
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("mets", edm::InputTag("slimmedMETs"));
  desc.add<edm::InputTag>("triggersPAT", edm::InputTag("TriggerResults","","PAT"));
  desc.add<edm::InputTag>("triggersHLT", edm::InputTag("TriggerResults","","HLT"));
  desc.add<std::string>("HLTName", std::string("placeholderHLT"));
  desc.add<bool>("isCRAB", bool(false));
  desc.add<edm::InputTag>("trigobjs", edm::InputTag("slimmedPatTrigger"));
  desc.add<edm::InputTag>("ecalBadCalibReducedMINIAODFilter", edm::InputTag("ecalBadCalibReducedMINIAODFilter"));

  descriptions.addWithDefaultLabel(desc);

}

extern char const charElectron[] = "electron";
extern char const charMuon[] = "muon";
extern char const charTauele[] = "tauele";
extern char const charTaumu[] = "taumu";

using plotterZToElecProbeTrk = plotterZToLepProbeTrk<charElectron>;
using plotterZToMuonProbeTrk = plotterZToLepProbeTrk<charMuon>;
using plotterZToTauEleProbeTrk = plotterZToLepProbeTrk<charTauele>;
using plotterZToTauMuProbeTrk = plotterZToLepProbeTrk<charTaumu>;

//define this as a plug-in
DEFINE_FWK_MODULE(plotterZToElecProbeTrk);
DEFINE_FWK_MODULE(plotterZToMuonProbeTrk);
DEFINE_FWK_MODULE(plotterZToTauEleProbeTrk);
DEFINE_FWK_MODULE(plotterZToTauMuProbeTrk);
