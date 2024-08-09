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
//
// class declaration
//

// If the filter does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDFilter<>
// This will improve performance in multithreaded jobs.

class basicSel : public edm::one::EDFilter<edm::one::SharedResources> {
public:
  explicit basicSel(const edm::ParameterSet&);
  ~basicSel() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  std::vector<std::string> commonCuts;

private:
  bool filter(edm::Event&, const edm::EventSetup&) override;
  edm::Service<TFileService> fs_;
  std::map<std::string, TH1D *> oneDHists_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken_;
  edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken_;
  edm::EDGetTokenT<std::vector<pat::MET>> metsToken_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> pfCandToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggersPATToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggersHLTToken_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> trigobjsToken_;
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
basicSel::basicSel(const edm::ParameterSet& iConfig)
    : verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
      jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      metsToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
      pfCandToken_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("pfCandidates"))),
      triggersPATToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggersPAT"))),
      triggersHLTToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggersHLT"))),
      trigobjsToken_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("trigobjs"))),
      caloGeometryToken_(esConsumes()),
      ecalStatusToken_(esConsumes()),
      ecalBadCalibFilterUpdateToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("ecalBadCalibReducedMINIAODFilter"))) {
  //now do what ever initialization is needed

  HLTName_ = iConfig.getParameter<std::string>("HLTName");
  isCRAB_ = iConfig.getParameter<bool>("isCRAB");
  isMETTriggers_ = iConfig.getParameter<bool>("isMETTriggers");

  if(isMETTriggers_) HLTName_ = "MET Triggers";

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
    HLTName_,
    "METFilters",
    "passecalBadCalibFilterUpdate",
    "good primary vertex",
    "met with noMuPt > 120",
    "jet with pt > 110",
    "jet with fabs ( eta ) < 2.4",
    "jet passing TightLepVeto ID",
    "veto pairs of jets with #Delta#phi > 2.5",
    "#Delta#phi(E_{T}^{miss},jet) > 0.5 "
  };

  int nBins = int(commonCuts.size());
  double maxRange = double(nBins) - 0.5;

  oneDHists_["cutflow"] = fs_->make<TH1D>("cutflow", "", nBins, -0.5, maxRange);
  oneDHists_["selection"] = fs_->make<TH1D>("selection", "", nBins, -0.5, maxRange);

  for(int i = 1; i <= int(commonCuts.size()); i++){
    if(i <= int(commonCuts.size())) {
      oneDHists_.at("cutflow")->GetXaxis()->SetBinLabel(i,commonCuts[i-1].c_str());
      oneDHists_.at("selection")->GetXaxis()->SetBinLabel(i,commonCuts[i-1].c_str());
    }
  }

}

basicSel::~basicSel() {

  plotPrintFunctions::printCuts(oneDHists_);

}

//
// member functions
//

// ------------ method called for each event  ------------
bool basicSel::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  std::vector<bool> passSel;
  std::vector<bool> passCut;

  for(int j = 0; j < (int(commonCuts.size()) - 1); ++j) {passSel.push_back(false); passCut.push_back(false);}

  bool isGood = false;

  caloGeometry = iSetup.getHandle(caloGeometryToken_);
  ecalStatus = iSetup.getHandle(ecalStatusToken_);

  helperFunctions::getChannelStatusMaps(ecalStatus,caloGeometry,EcalAllDeadChannelsValMap,EcalAllDeadChannelsBitMap);

  edm::Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(jetsToken_, jets);

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

  if(!isMETTriggers_) {if(helperFunctions::passLepHLTPath(iEvent,triggerBitsHLT,HLTName_)) {passSel[0] = true; passCut[0] = true;}}

  if(isMETTriggers_) {if(helperFunctions::passMETHLTPath(iEvent,triggerBitsHLT)) {passSel[0] = true; passCut[0] = true;}}

  if(helperFunctions::passMETFilters(iEvent,triggerBitsPAT)) {passSel[1] = true; if(passCut[0]) passCut[1] = true;}

  if(*passecalBadCalibFilterUpdate) {passSel[2] = true; if(passCut[1]) passCut[2] = true;}

  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(verticesToken_, vertices);
  // const reco::Vertex &pv = vertices->at(0);

  for (const auto &pv : *vertices) {
    if(helperFunctions::isGoodPV(pv)) {passSel[3] = true; if(passCut[2]) passCut[3] = true;}
  }

  if(metNoMu.Mod() > 120.0) {passSel[4] = true; if(passCut[3]) passCut[4] = true;}

  // Block below is done here because it needs to compare all jets, with the IsValidJet, which would imply in different results in jet cuts
  // if done together
  double dijetMaxDeltaPhi = -999.;
  double phiJetLeading = -999.;
  int idx1 = -1;
  for (const auto &jet1 : *jets) {
    idx1++;
    if (!helperFunctions::IsValidJet(jet1)) continue;
    int idx2 = -1;
    if (phiJetLeading == -999.) phiJetLeading = jet1.phi();
    for (const auto &jet2 : *jets) {
      idx2++;
      if (idx2 <= idx1)              continue;  // Avoid double-counting.
      if (!helperFunctions::IsValidJet(jet2)) continue;
      double dPhi = fabs(deltaPhi(jet1.phi(), jet2.phi()));
      if (dPhi > dijetMaxDeltaPhi) {
        dijetMaxDeltaPhi = dPhi;
      }
    }
  }
  
  std::vector<bool> auxPassCutJet;
  for(int j = 0; j < 3; ++j) auxPassCutJet.push_back(false);

  for (const auto &jet : *jets) {
    if(jet.pt() > 110.0)
      {passSel[5] = true; if(passCut[4]) auxPassCutJet[0] = true;}
    if(fabs(jet.eta()) < 2.4)
      {passSel[6] = true; if(auxPassCutJet[0]) auxPassCutJet[1] = true;}
    if(helperFunctions::jetPassesTightLepVeto(jet))
      {passSel[7] = true; if(auxPassCutJet[1]) auxPassCutJet[2] = true;}
    for(int j = 0; j < 3; ++j){
      if(auxPassCutJet[j]) passCut[j+5] = true;
      auxPassCutJet[j] = false;
    }
  }

  if(dijetMaxDeltaPhi < 2.5)
      {passSel[8] = true; if(passCut[7]) passCut[8] = true;}
  if(fabs(deltaPhi(metNoMu.Phi(),phiJetLeading)) > 0.5)
      {passSel[9] = true; if(passCut[8]) passCut[9] = true;}

  for(int i = 1; i <= int(passSel.size()); i++){
    if(passSel[i-1]) oneDHists_.at("selection")->Fill(i);
    if(passCut[i-1]) oneDHists_.at("cutflow")->Fill(i);
    if(i == int(passSel.size()) && passCut[i-1]) isGood = true;
  }

  return isGood;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void basicSel::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("mets", edm::InputTag("slimmedMETs"));
  desc.add<edm::InputTag>("pfCandidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("triggersPAT", edm::InputTag("TriggerResults","","PAT"));
  desc.add<edm::InputTag>("triggersHLT", edm::InputTag("TriggerResults","","HLT"));
  desc.add<edm::InputTag>("trigobjs", edm::InputTag("slimmedPatTrigger"));
  desc.add<std::string>("HLTName", std::string("placeholderHLT"));
  desc.add<edm::InputTag>("ecalBadCalibReducedMINIAODFilter", edm::InputTag("ecalBadCalibReducedMINIAODFilter"));
  desc.add<bool>("isCRAB", bool(false));
  desc.add<bool>("isMETTriggers", bool(false));

  descriptions.addWithDefaultLabel(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(basicSel);
