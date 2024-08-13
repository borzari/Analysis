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

template<class T>
class plotterSingLepCRNoJetSelSkim : public edm::one::EDFilter<edm::one::SharedResources> {
public:
  explicit plotterSingLepCRNoJetSelSkim(const edm::ParameterSet&);
  ~plotterSingLepCRNoJetSelSkim() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  bool filter(edm::Event&, const edm::EventSetup&) override;
  edm::Service<TFileService> fs_;
  std::map<std::string, TH1D *> oneDHists_;
  std::map<std::string, TH2D *> twoDHists_;

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
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
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
plotterSingLepCRNoJetSelSkim<T>::plotterSingLepCRNoJetSelSkim(const edm::ParameterSet& iConfig)
    : verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
      leptonsToken_(consumes<std::vector<T>>(iConfig.getParameter<edm::InputTag>("leptons"))),
      triggersPATToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggersPAT"))),
      triggersHLTToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggersHLT"))),
      trigobjsToken_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("trigobjs"))),
      tracksToken_(consumes<std::vector<pat::IsolatedTrack>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      metsToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
      pfCandToken_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("pfCandidates"))),
      generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
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

  std::string elecStr = "electron";
  std::string muonStr = "muon";
  std::string tauStr = "tau";
  std::string trackStr = "track";

  TFileDirectory NLayers4 = fs_->mkdir("NLayers4");
  TFileDirectory NLayers5 = fs_->mkdir("NLayers5");
  TFileDirectory NLayers6Plus = fs_->mkdir("NLayers6Plus");

  plotPrintFunctions::createCommonHists(fs_,oneDHists_,twoDHists_,NLayers4,4);
  plotPrintFunctions::createCommonHists(fs_,oneDHists_,twoDHists_,NLayers5,5);
  plotPrintFunctions::createCommonHists(fs_,oneDHists_,twoDHists_,NLayers6Plus,6);
  if(std::is_same<T, pat::Electron>::value) {
    plotPrintFunctions::createObjHists<pat::Electron>(fs_,oneDHists_,elecStr,NLayers4,4);
    plotPrintFunctions::createObjHists<pat::Electron>(fs_,oneDHists_,elecStr,NLayers5,5);
    plotPrintFunctions::createObjHists<pat::Electron>(fs_,oneDHists_,elecStr,NLayers6Plus,6);
  }
  if(std::is_same<T, pat::Muon>::value) {
    plotPrintFunctions::createObjHists<pat::Muon>(fs_,oneDHists_,muonStr,NLayers4,4);
    plotPrintFunctions::createObjHists<pat::Muon>(fs_,oneDHists_,muonStr,NLayers5,5);
    plotPrintFunctions::createObjHists<pat::Muon>(fs_,oneDHists_,muonStr,NLayers6Plus,6);
  }
  if(std::is_same<T, pat::Tau>::value) {
    plotPrintFunctions::createObjHists<pat::Tau>(fs_,oneDHists_,tauStr,NLayers4,4);
    plotPrintFunctions::createObjHists<pat::Tau>(fs_,oneDHists_,tauStr,NLayers5,5);
    plotPrintFunctions::createObjHists<pat::Tau>(fs_,oneDHists_,tauStr,NLayers6Plus,6);
  }
  plotPrintFunctions::createObjHists<pat::IsolatedTrack>(fs_,oneDHists_,trackStr,NLayers4,4);
  plotPrintFunctions::createObjHists<pat::IsolatedTrack>(fs_,oneDHists_,trackStr,NLayers5,5);
  plotPrintFunctions::createObjHists<pat::IsolatedTrack>(fs_,oneDHists_,trackStr,NLayers6Plus,6);

}

template<class T>
plotterSingLepCRNoJetSelSkim<T>::~plotterSingLepCRNoJetSelSkim() {

}

//
// member functions
//

// ------------ method called for each event  ------------
template<class T>
bool plotterSingLepCRNoJetSelSkim<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  std::vector<bool> passSel;
  std::vector<bool> passCut;

  VariablesToPlot variables;

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

  variables.met = met;
  variables.metNoMu = metNoMu;

  edm::Handle<GenEventInfoProduct> generator;
  iEvent.getByToken(generatorToken_, generator);

  edm::Handle<bool> passecalBadCalibFilterUpdate;
  iEvent.getByToken(ecalBadCalibFilterUpdateToken_, passecalBadCalibFilterUpdate);

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

  if(std::is_same<T, pat::Electron>::value){

    int startElecIdx = 1;

    for (const auto& electron : *electrons) {

      selectingFunctions::singElecSel(passSel,passCut,auxPassCut,electron,startElecIdx,pv,selElectrons,35.0);

      for(int j = 0; j < int(auxPassCut.size()); ++j) auxPassCut[j] = false;

    }

    variables.electrons = selElectrons;
  
  }

  if(std::is_same<T, pat::Muon>::value){

    int startMuonIdx = 1;
    
    for (const auto& muon : *muons) {
      
      selectingFunctions::singMuonSel(passSel,passCut,auxPassCut,muon,startMuonIdx,pv,selMuons);
      
      for(int j = 0; j < int(auxPassCut.size()); ++j) auxPassCut[j] = false;
        
    }

    variables.muons = selMuons;
  
  }

  if(std::is_same<T, pat::Tau>::value){

    int startTauIdx = 1;
    
    for (const auto& tau : *taus) {
      
      selectingFunctions::singTauSel(passSel,passCut,auxPassCut,tau,startTauIdx,selTaus);

      for(int j = 0; j < int(auxPassCut.size()); ++j) auxPassCut[j] = false;
        
    }

    variables.taus = selTaus;
  
  }

  int startTrackIdx = 1;

  double cutTrackPt = 0.0;
  if(std::is_same<T, pat::Electron>::value || std::is_same<T, pat::Muon>::value) cutTrackPt = 35.0;
  if(std::is_same<T, pat::Tau>::value) cutTrackPt = 55.0;

  int getStTrkIdxLep = 0;

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

    for(int j = 0; j < int(auxPassCut.size()); ++j) auxPassCut[j] = false;

  }

  variables.tracks = selTracks;

  TVector2 metNoMuNoLep;
  double deltaPhiMetJetLeadingVsLeptonMet;
  double deltaPhiMetJetLeadingVsLeptonMetNoMu;
  double deltaPhiMetJetLeadingVsLeptonMetNoMuMinusOnePt;

  std::srand(std::time(0)); // use current time as seed for random generator

  if(std::is_same<T, pat::Electron>::value){
    if(tkMatchElectrons.size() > 0)
    {int random_pos = std::rand() % tkMatchElectrons.size();  // Modulo to restrict the number of random values to be at most A.size()-1
    pat::Electron randomLep = tkMatchElectrons[random_pos];
    metNoMuNoLep = TVector2(metNoMu.Px() + randomLep.px(), metNoMu.Py() + randomLep.py());
    variables.metNoMuNoLep = metNoMuNoLep;
    deltaPhiMetJetLeadingVsLeptonMet = fabs(deltaPhi(met.phi(),leadingJet.phi()));
    deltaPhiMetJetLeadingVsLeptonMetNoMu = fabs(deltaPhi(metNoMu.Phi(),leadingJet.phi()));
    deltaPhiMetJetLeadingVsLeptonMetNoMuMinusOnePt = fabs(deltaPhi(metNoMuNoLep.Phi(),leadingJet.phi()));
    variables.DeltaPhiJetMET = deltaPhiMetJetLeadingVsLeptonMet;
    variables.DeltaPhiJetMETNoMu = deltaPhiMetJetLeadingVsLeptonMetNoMu;
    variables.DeltaPhiJetMETNoMuNoLep = deltaPhiMetJetLeadingVsLeptonMetNoMuMinusOnePt;}
  }

  if(std::is_same<T, pat::Muon>::value){
    if(tkMatchMuons.size() > 0)
    {int random_pos = std::rand() % tkMatchMuons.size();  // Modulo to restrict the number of random values to be at most A.size()-1
    pat::Muon randomLep = tkMatchMuons[random_pos];
    metNoMuNoLep = TVector2(metNoMu.Px(), metNoMu.Py());
    variables.metNoMuNoLep = metNoMuNoLep;
    deltaPhiMetJetLeadingVsLeptonMet = fabs(deltaPhi(met.phi(),leadingJet.phi()));
    deltaPhiMetJetLeadingVsLeptonMetNoMu = fabs(deltaPhi(metNoMu.Phi(),leadingJet.phi()));
    deltaPhiMetJetLeadingVsLeptonMetNoMuMinusOnePt = fabs(deltaPhi(metNoMuNoLep.Phi(),leadingJet.phi()));
    variables.DeltaPhiJetMET = deltaPhiMetJetLeadingVsLeptonMet;
    variables.DeltaPhiJetMETNoMu = deltaPhiMetJetLeadingVsLeptonMetNoMu;
    variables.DeltaPhiJetMETNoMuNoLep = deltaPhiMetJetLeadingVsLeptonMetNoMuMinusOnePt;}
  }

  if(std::is_same<T, pat::Tau>::value){
    if(tkMatchTaus.size() > 0)
    {int random_pos = std::rand() % tkMatchTaus.size();  // Modulo to restrict the number of random values to be at most A.size()-1
    pat::Tau randomLep = tkMatchTaus[random_pos];
    metNoMuNoLep = TVector2(metNoMu.Px() + randomLep.px(), metNoMu.Py() + randomLep.py());
    variables.metNoMuNoLep = metNoMuNoLep;
    deltaPhiMetJetLeadingVsLeptonMet = fabs(deltaPhi(met.phi(),leadingJet.phi()));
    deltaPhiMetJetLeadingVsLeptonMetNoMu = fabs(deltaPhi(metNoMu.Phi(),leadingJet.phi()));
    deltaPhiMetJetLeadingVsLeptonMetNoMuMinusOnePt = fabs(deltaPhi(metNoMuNoLep.Phi(),leadingJet.phi()));
    variables.DeltaPhiJetMET = deltaPhiMetJetLeadingVsLeptonMet;
    variables.DeltaPhiJetMETNoMu = deltaPhiMetJetLeadingVsLeptonMetNoMu;
    variables.DeltaPhiJetMETNoMuNoLep = deltaPhiMetJetLeadingVsLeptonMetNoMuMinusOnePt;}
  }

  double sf = 1.0;
  std::string elecInputFile = "/home/brenoorzari/CMSSW_13_0_13/src/Analysis/Helper/data/electronSFs.root";
  std::string muonInputFile = "/home/brenoorzari/CMSSW_13_0_13/src/Analysis/Helper/data/muonSFs.root";
  std::string elecReco = "electronReco2022;1";
  std::string elecID = "electronID2022Tight;1";
  std::string muonID = "muonID2022Tight;1";
  std::string muonIso = "muonIso2022TightTightID;1";
  std::string muonTrigger = "muonTrigger2022IsoMu24;1";

  sf *= (*generator).weight() / fabs ((*generator).weight ());

  if(std::is_same<T, pat::Electron>::value){
    sfFunctions::sfElectrons(sf, elecInputFile, elecReco, selElectrons);
    sfFunctions::sfElectrons(sf, elecInputFile, elecID, selElectrons);
  }
  if(std::is_same<T, pat::Muon>::value){
    sfFunctions::sfMuons(sf, muonInputFile, muonID, selMuons);
    sfFunctions::sfMuons(sf, muonInputFile, muonIso, selMuons);
    sfFunctions::sfMuons(sf, muonInputFile, muonTrigger, selMuons);
  }
    
  plotPrintFunctions::plotVariables(oneDHists_, twoDHists_, variables, 4, sf);
  plotPrintFunctions::plotVariables(oneDHists_, twoDHists_, variables, 5, sf);
  plotPrintFunctions::plotVariables(oneDHists_, twoDHists_, variables, 6, sf);

  return isGood;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<class T>
void plotterSingLepCRNoJetSelSkim<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("leptons", edm::InputTag("placeholderLeptons"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("isolatedTracks"));
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("mets", edm::InputTag("slimmedMETs"));
  desc.add<edm::InputTag>("pfCandidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("generator", edm::InputTag("generator"));
  desc.add<edm::InputTag>("triggersPAT", edm::InputTag("TriggerResults","","PAT"));
  desc.add<edm::InputTag>("triggersHLT", edm::InputTag("TriggerResults","","HLT"));
  desc.add<std::string>("HLTName", std::string("placeholderHLT"));
  desc.add<bool>("isCRAB", bool(false));
  desc.add<bool>("isMETTriggers", bool(false));
  desc.add<edm::InputTag>("trigobjs", edm::InputTag("slimmedPatTrigger"));
  desc.add<edm::InputTag>("ecalBadCalibReducedMINIAODFilter", edm::InputTag("ecalBadCalibReducedMINIAODFilter"));

  descriptions.addWithDefaultLabel(desc);

}

using plotterSingElecCRNoJetSelSkim = plotterSingLepCRNoJetSelSkim<pat::Electron>;
using plotterSingMuonCRNoJetSelSkim = plotterSingLepCRNoJetSelSkim<pat::Muon>;
using plotterSingTauCRNoJetSelSkim = plotterSingLepCRNoJetSelSkim<pat::Tau>;

//define this as a plug-in
DEFINE_FWK_MODULE(plotterSingElecCRNoJetSelSkim);
DEFINE_FWK_MODULE(plotterSingMuonCRNoJetSelSkim);
DEFINE_FWK_MODULE(plotterSingTauCRNoJetSelSkim);
