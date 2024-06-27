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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
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

// If the filter does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDFilter<>
// This will improve performance in multithreaded jobs.

template<class T>
class zToLepProbeTrk : public edm::one::EDFilter<edm::one::SharedResources> {
public:
  explicit zToLepProbeTrk(const edm::ParameterSet&);
  ~zToLepProbeTrk() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  bool isMatchedToElecTriggerObject (const edm::Event &, const edm::TriggerResults &, const pat::Electron &, const std::vector<pat::TriggerObjectStandAlone> &, const std::string &, const std::string &, const double = 0.1);
  bool isMatchedToMuonTriggerObject (const edm::Event &, const edm::TriggerResults &, const pat::Muon &, const std::vector<pat::TriggerObjectStandAlone> &, const std::string &, const std::string &, const double = 0.1);
  bool isMatchedToTauTriggerObject (const edm::Event &, const edm::TriggerResults &, const pat::Tau &, const std::vector<pat::TriggerObjectStandAlone> &, const std::string &, const std::string &, const double = 0.1);
  bool passesDecayModeReconstruction (const pat::Tau &);
  bool passesLightFlavorRejection (const pat::Tau &);
  bool inTOBCrack (const pat::IsolatedTrack &);
  const int extraMissingMiddleHits (const pat::IsolatedTrack &) const;
  const int hitDrop_missingMiddleHits (const pat::IsolatedTrack &) const;
  double dRMinJet (const pat::IsolatedTrack &, const std::vector<pat::Jet> &);
  double deltaRToClosestMuon(const pat::IsolatedTrack &, const std::vector<pat::Muon> &);
  double deltaRToClosestTauHad(const pat::IsolatedTrack &, const std::vector<pat::Tau> &);


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
  edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggersPATToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggersHLTToken_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> trigobjsToken_;
  std::string HLTName_;

  std::vector<std::string> commonCuts{"Total",HLTName_,"METFilters"};

  std::vector<std::string> electronCuts{"electron pt > 32 GeV","electron isMatchedToTriggerObject","electron fabs(#eta) < 2.1","electron passesVID_tightID (ID + iso)","electron |d0| < 0.05, 0.10 (EB, EE)","electron |dz| < 0.10, 0.20 (EB, EE)"};

  std::vector<std::string> trackCuts{"track pt > 30 GeV","track fabs(#eta) < 2.1","track fabs ( eta ) < 0.15 || fabs ( eta ) > 0.35","track fabs ( eta ) < 1.42 || fabs ( eta ) > 1.65","track fabs ( eta ) < 1.55 || fabs ( eta ) > 1.85","track !inTOBCrack","track hitPattern_.numberOfValidPixelHits >= 4","track hitPattern_.numberOfValidHits >= 4","track missingInnerHits == 0","track hitDrop_missingMiddleHits == 0","track ((pfIsolationDR03_.chargedHadronIso + pfIsolationDR03_.puChargedHadronIso) / pt) < 0.05","track |d0| < 0.02","track |dz| < 0.5","track dRMinJet > 0.5","track deltaRToClosestMuon > 0.15","track deltaRToClosestTauHad > 0.15"};

};

template<class T>
zToLepProbeTrk<T>::zToLepProbeTrk(const edm::ParameterSet& iConfig)
    : verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
      electronsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
      muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      tausToken_(consumes<std::vector<pat::Tau>>(iConfig.getParameter<edm::InputTag>("taus"))),
      tracksToken_(consumes<std::vector<pat::IsolatedTrack>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      triggersPATToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggersPAT"))),
      triggersHLTToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggersHLT"))),
      trigobjsToken_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("trigobjs"))) {
  //now do what ever initialization is needed

  HLTName_ = iConfig.getParameter<std::string>("HLTName");

  commonCuts[1] = HLTName_; // HLTName_ is only defined here, so it needs to be updated in commonCuts

  int nBins = int(commonCuts.size() + trackCuts.size());
  double maxRange = 0.0;
  std::vector<std::string> lepCuts;

  if(std::is_same<T, pat::Electron>::value){nBins = nBins + int(electronCuts.size()); maxRange = double(nBins) - 0.5; lepCuts = electronCuts;}
  if(std::is_same<T, pat::Muon>::value){nBins = 7; maxRange = 6.5;}
  if(std::is_same<T, pat::Tau>::value){nBins = 6; maxRange = 5.5;}

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

}

template<class T>
zToLepProbeTrk<T>::~zToLepProbeTrk() {
  int maxSize = 0;
  for(int i = 1; i <= int(oneDHists_.at("cutflow")->GetNbinsX()); i++){
    std::string label = oneDHists_.at("cutflow")->GetXaxis()->GetLabels()->At(i-1)->GetName();
    if(int(label.size()) > maxSize) maxSize = int(label.size());
  }

  std::string cut = "Cut Value";
  std::cout << cut << std::string((maxSize - int(cut.size()) + 1),' ') << "\t\tCumul.\t\tIndiv." << std::endl;

  double maxValue = double(oneDHists_.at("selection")->GetBinContent(1));
  for(int i = 1; i <= int(oneDHists_.at("cutflow")->GetNbinsX()); i++){
    std::cout << std::setprecision(3) << std::fixed;
    std::string label = oneDHists_.at("cutflow")->GetXaxis()->GetLabels()->At(i-1)->GetName();
    int size = int(label.size());
    int nBlank = maxSize - size + 1;
    std::cout << oneDHists_.at("cutflow")->GetXaxis()->GetLabels()->At(i-1)->GetName() << std::string(nBlank,' ') << "\t\t" << double(oneDHists_.at("cutflow")->GetBinContent(i)) * 100.0 / maxValue << " % \t" << double(oneDHists_.at("selection")->GetBinContent(i)) * 100.0 / maxValue << " %" << std::endl;
  }
}

// ------------ method called for each event  ------------
template<class T>
bool zToLepProbeTrk<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  std::vector<bool> passSel;
  std::vector<bool> passCut;

  for(int j = 0; j < (int(commonCuts.size()) - 1); ++j) {passSel.push_back(false); passCut.push_back(false);}

  bool isGood = false;

  int pass_HLT = 0;

  Int_t pass_Flag_goodVertices = 0;
  Int_t pass_Flag_globalTightHalo2016Filter = 0;
  Int_t pass_Flag_HBHENoiseFilter = 0;
  Int_t pass_Flag_HBHENoiseIsoFilter = 0;
  Int_t pass_Flag_EcalDeadCellTriggerPrimitiveFilter = 0;
  Int_t pass_Flag_BadPFMuonFilter = 0;
  Int_t pass_Flag_BadChargedCandidateFilter = 0;
  Int_t pass_Flag_ecalBadCalibFilter = 0;

  edm::Handle<std::vector<pat::IsolatedTrack>> tracks;
  iEvent.getByToken(tracksToken_, tracks);

  edm::Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(jetsToken_, jets);

  edm::Handle<std::vector<pat::Electron>> electrons;
  iEvent.getByToken(electronsToken_, electrons);  

  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);  

  edm::Handle<std::vector<pat::Tau>> taus;
  iEvent.getByToken(tausToken_, taus);

  oneDHists_.at("selection")->Fill(0);
  oneDHists_.at("cutflow")->Fill(0);

  edm::Handle<edm::TriggerResults> triggerBitsPAT;
  iEvent.getByToken(triggersPATToken_, triggerBitsPAT);

  edm::Handle<edm::TriggerResults> triggerBitsHLT;
  iEvent.getByToken(triggersHLTToken_, triggerBitsHLT);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjs;
  iEvent.getByToken(trigobjsToken_, triggerObjs);

  const edm::TriggerNames &allTriggerNamesHLT = iEvent.triggerNames(*triggerBitsHLT);

  for(unsigned i = 0; i < allTriggerNamesHLT.size(); i++) {
    std::string thisName = allTriggerNamesHLT.triggerName(i);
    if (thisName.find(HLTName_) == 0) pass_HLT = triggerBitsHLT->accept(i);
  }

  if(pass_HLT == 1) {passSel[0] = true; passCut[0] = true;}

  const edm::TriggerNames &allTriggerNamesPAT = iEvent.triggerNames(*triggerBitsPAT);

  for(unsigned i = 0; i < allTriggerNamesPAT.size(); i++) {
    std::string thisName = allTriggerNamesPAT.triggerName(i);
    if (thisName.find("Flag_goodVertices") == 0) pass_Flag_goodVertices = triggerBitsPAT->accept(i);
    if (thisName.find("Flag_globalTightHalo2016Filter") == 0) pass_Flag_globalTightHalo2016Filter = triggerBitsPAT->accept(i);
    if (thisName.find("Flag_HBHENoiseFilter") == 0) pass_Flag_HBHENoiseFilter = triggerBitsPAT->accept(i);
    if (thisName.find("Flag_HBHENoiseIsoFilter") == 0) pass_Flag_HBHENoiseIsoFilter = triggerBitsPAT->accept(i);
    if (thisName.find("Flag_EcalDeadCellTriggerPrimitiveFilter") == 0) pass_Flag_EcalDeadCellTriggerPrimitiveFilter = triggerBitsPAT->accept(i);
    if (thisName.find("Flag_BadPFMuonFilter") == 0) pass_Flag_BadPFMuonFilter = triggerBitsPAT->accept(i);
    if (thisName.find("Flag_BadChargedCandidateFilter") == 0) pass_Flag_BadChargedCandidateFilter = triggerBitsPAT->accept(i);
    if (thisName.find("Flag_ecalBadCalibFilter") == 0) pass_Flag_ecalBadCalibFilter = triggerBitsPAT->accept(i);
  }

  if(pass_Flag_goodVertices == 1 && pass_Flag_globalTightHalo2016Filter == 1 && pass_Flag_HBHENoiseFilter == 1 && pass_Flag_HBHENoiseIsoFilter == 1 && pass_Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && pass_Flag_BadPFMuonFilter == 1 && pass_Flag_BadChargedCandidateFilter == 1 && pass_Flag_ecalBadCalibFilter == 1) {passSel[1] = true; if(passCut[0]) passCut[1] = true;}

  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(verticesToken_, vertices);
  const reco::Vertex &pv = vertices->at(0);

  std::vector<bool> auxPassCut;

  if(std::is_same<T, pat::Electron>::value){

    for(int j = 0; j < int(electronCuts.size()); ++j) {passSel.push_back(false); passCut.push_back(false); auxPassCut.push_back(false);}
    
    for (const auto& electron : *electrons) {
      
      if(electron.pt() > 32.) {passSel[2] = true; if(passCut[1]) auxPassCut[0] = true;}
      if(isMatchedToElecTriggerObject (iEvent, *triggerBitsHLT, electron, *triggerObjs, "hltEgammaCandidates::HLT", "hltEle32WPTightGsfTrackIsoFilter")) {passSel[3] = true; if(auxPassCut[0]) auxPassCut[1] = true;}
      if(abs(electron.eta()) < 2.1) {passSel[4] = true; if(auxPassCut[1]) auxPassCut[2] = true;}
      if(electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-tight")) {passSel[5] = true; if(auxPassCut[2]) auxPassCut[3] = true;}
      if(((fabs (electron.superCluster()->eta()) <= 1.479) && (fabs (((electron.vx() - pv.x()) * electron.py() - (electron.vy() - pv.y()) * electron.px()) / electron.pt()) < 0.05)) || ((fabs (electron.superCluster()->eta()) >  1.479) && (fabs (((electron.vx() - pv.x()) * electron.py() - (electron.vy() - pv.y()) * electron.px()) / electron.pt()) < 0.10))) {passSel[6] = true; if(auxPassCut[3]) auxPassCut[4] = true;}
      if(((fabs (electron.superCluster()->eta()) <= 1.479) && (fabs ((electron.vz() - pv.z()) - ((electron.vx() - pv.x()) * electron.px() + (electron.vy() - pv.y()) * electron.py()) / electron.pt() * electron.pz() / electron.pt()) < 0.10)) || ((fabs (electron.superCluster()->eta()) >  1.479) && (fabs ((electron.vz() - pv.z()) - ((electron.vx() - pv.x()) * electron.px() + (electron.vy() - pv.y()) * electron.py()) / electron.pt() * electron.pz() / electron.pt()) < 0.20))) {passSel[7] = true; if(auxPassCut[4]) auxPassCut[5] = true;}
      
      for(int j = 0; j < int(auxPassCut.size()); ++j){
        if(auxPassCut[j]) passCut[j+2] = true;
        auxPassCut[j] = false;
      }
        
    }
  
  }

  // if(std::is_same<T, pat::Muon>::value){
    
  //   for (const auto& muon : *muons) {
      
  //     if(muon.pt() > 26.) {passPtSel = true; if(passMETFiltersCut) auxPassPtCut = true;}
  //     if(abs(muon.eta()) < 2.1) {passEtaSel = true; if(auxPassPtCut) auxPassEtaCut = true;}
  //     if(muon.isTightMuon(pv)) {passIDSel = true; if(auxPassEtaCut) auxPassIDCut = true;}
  //     if(((muon.pfIsolationR04().sumChargedHadronPt + std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt()) < 0.15) {passPFIsoSel = true; if(auxPassIDCut) auxPassPFIsoCut = true;}
      
  //     if(auxPassPtCut && auxPassEtaCut && auxPassIDCut && auxPassPFIsoCut) passPFIsoCut = true;
  //     if(auxPassPtCut && auxPassEtaCut && auxPassIDCut) passIDCut = true;
  //     if(auxPassPtCut && auxPassEtaCut) passEtaCut = true;
  //     if(auxPassPtCut) passPtCut = true;
          
  //     auxPassPtCut = false;
  //     auxPassEtaCut = false;
  //     auxPassIDCut = false;
  //     auxPassPFIsoCut = false;
        
  //   }
  
  // }

  // if(std::is_same<T, pat::Tau>::value){
    
  //   for (const auto& tau : *taus) {
      
  //     if(tau.pt() > 50.) {passPtSel = true; if(passMETFiltersCut) auxPassPtCut = true;}
  //     if(abs(tau.eta()) < 2.1) {passEtaSel = true; if(auxPassPtCut) auxPassEtaCut = true;}
  //     if(passesDecayModeReconstruction(tau) && passesLightFlavorRejection(tau)) {passIDSel = true; if(auxPassEtaCut) auxPassIDCut = true;}
      
  //     if(auxPassPtCut && auxPassEtaCut && auxPassIDCut) passIDCut = true;
  //     if(auxPassPtCut && auxPassEtaCut) passEtaCut = true;
  //     if(auxPassPtCut) passPtCut = true;
          
  //     auxPassPtCut = false;
  //     auxPassEtaCut = false;
  //     auxPassIDCut = false;
        
  //   }
  
  // }

// >= 1 tracks with isFiducialElectronTrack                                                                       1983.0         19.830%        100.000%
// >= 1 tracks with isFiducialMuonTrack                                                                           1983.0         19.830%        100.000%
// >= 1 tracks with isFiducialECALTrack                                                                           1863.0         18.630%         98.170%

  for(int j = 0; j < int(trackCuts.size()); ++j) {passSel.push_back(false); passCut.push_back(false); auxPassCut.push_back(false);}

  for (const auto& track : *tracks) {
      
    if(std::is_same<T, pat::Electron>::value) {if(track.pt() > 30.) {passSel[8] = true; if(passCut[7]) auxPassCut[6] = true;}}
    if(fabs(track.eta()) < 2.1) {passSel[9] = true; if(auxPassCut[6]) auxPassCut[7] = true;}
    if((fabs(track.eta()) < 0.15) || (fabs(track.eta()) > 0.35)) {passSel[10] = true; if(auxPassCut[7]) auxPassCut[8] = true;}
    if((fabs(track.eta()) < 1.42) || (fabs(track.eta()) > 1.65)) {passSel[11] = true; if(auxPassCut[8]) auxPassCut[9] = true;}
    if((fabs(track.eta()) < 1.55) || (fabs(track.eta()) > 1.85)) {passSel[12] = true; if(auxPassCut[9]) auxPassCut[10] = true;}

    // Need to include fiducial cuts here!!!

    if(!inTOBCrack(track)) {passSel[13] = true; if(auxPassCut[10]) auxPassCut[11] = true;}
    if(track.hitPattern().numberOfValidPixelHits() >= 4) {passSel[14] = true; if(auxPassCut[11]) auxPassCut[12] = true;}
    if(track.hitPattern().numberOfValidHits() >= 4) {passSel[15] = true; if(auxPassCut[12]) auxPassCut[13] = true;}
    if(track.lostInnerLayers() == 0) {passSel[16] = true; if(auxPassCut[13]) auxPassCut[14] = true;}
    if(hitDrop_missingMiddleHits(track) == 0) {passSel[17] = true; if(auxPassCut[14]) auxPassCut[15] = true;}
    if(((track.pfIsolationDR03().chargedHadronIso() + track.pfIsolationDR03().puChargedHadronIso()) / track.pt()) < 0.05) {passSel[18] = true; if(auxPassCut[15]) auxPassCut[16] = true;}
    if(fabs(track.dxy()) < 0.02) {passSel[19] = true; if(auxPassCut[16]) auxPassCut[17] = true;}
    if(fabs(track.dz()) < 0.5) {passSel[20] = true; if(auxPassCut[17]) auxPassCut[18] = true;}
    if(dRMinJet(track, *jets) > 0.5) {passSel[21] = true; if(auxPassCut[18]) auxPassCut[19] = true;}
    if(deltaRToClosestMuon(track, *muons) > 0.15) {passSel[22] = true; if(auxPassCut[19]) auxPassCut[20] = true;}
    if(deltaRToClosestTauHad(track, *taus) > 0.15) {passSel[23] = true; if(auxPassCut[20]) auxPassCut[21] = true;} // The individual efficiency of this is different than the original analysis, because it uses distinct selections that follow the Run 3 recommendations from the Tau POG https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3#Kinematic_tau_selection

    for(int j = int(electronCuts.size()); j < int(auxPassCut.size()); ++j){
      if(auxPassCut[j]) passCut[j+2] = true;
      auxPassCut[j] = false;
    }

  }

  for(int i = 1; i <= int(passSel.size()); i++){
    if(passSel[i-1]) oneDHists_.at("selection")->Fill(i);
    if(passCut[i-1]) oneDHists_.at("cutflow")->Fill(i);
    if(i == int(passSel.size()) && passCut[i-1]) isGood = true;
  }

  return isGood;
}

template<class T>
bool zToLepProbeTrk<T>::isMatchedToElecTriggerObject (const edm::Event &event, const edm::TriggerResults &triggers, const pat::Electron &obj, const std::vector<pat::TriggerObjectStandAlone> &trigObjs, const std::string &collection, const std::string &filter, const double dR)
{
  if(collection == "") return false;
  for(auto trigObj : trigObjs) {
    trigObj.unpackNamesAndLabels(event, triggers);
    if(trigObj.collection() != collection) continue;
    if(filter != "") {
      bool flag = false;
      for(const auto &filterLabel : trigObj.filterLabels ())
        if(filterLabel == filter) {
          flag = true;
          break;
        }
      if (!flag) continue;
    }
    if(deltaR (obj, trigObj) > dR) continue;
    return true;
  }
  return false;
}

template<class T>
bool zToLepProbeTrk<T>::isMatchedToMuonTriggerObject (const edm::Event &event, const edm::TriggerResults &triggers, const pat::Muon &obj, const std::vector<pat::TriggerObjectStandAlone> &trigObjs, const std::string &collection, const std::string &filter, const double dR)
{
  if(collection == "") return false;
  for(auto trigObj : trigObjs) {
    trigObj.unpackNamesAndLabels(event, triggers);
    if(trigObj.collection() != collection) continue;
    if(filter != "") {
      bool flag = false;
      for(const auto &filterLabel : trigObj.filterLabels ())
        if(filterLabel == filter) {
          flag = true;
          break;
        }
      if (!flag) continue;
    }
    if(deltaR (obj, trigObj) > dR) continue;
    return true;
  }
  return false;
}

template<class T>
bool zToLepProbeTrk<T>::isMatchedToTauTriggerObject (const edm::Event &event, const edm::TriggerResults &triggers, const pat::Tau &obj, const std::vector<pat::TriggerObjectStandAlone> &trigObjs, const std::string &collection, const std::string &filter, const double dR)
{
  if(collection == "") return false;
  for(auto trigObj : trigObjs) {
    trigObj.unpackNamesAndLabels(event, triggers);
    if(trigObj.collection() != collection) continue;
    if(filter != "") {
      bool flag = false;
      for(const auto &filterLabel : trigObj.filterLabels ())
        if(filterLabel == filter) {
          flag = true;
          break;
        }
      if (!flag) continue;
    }
    if(deltaR (obj, trigObj) > dR) continue;
    return true;
  }
  return false;
}

template<class T>
bool zToLepProbeTrk<T>::passesDecayModeReconstruction (const pat::Tau &tau){
  return (tau.tauID("decayModeFindingNewDMs") && (tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") || (tau.tauID("chargedIsoPtSumdR03")+std::max(0.,tau.tauID("neutralIsoPtSumdR03")-0.072*tau.tauID("puCorrPtSum"))<2.5) || tau.tauID("byVVVLooseDeepTau2017v2p1VSjet") || tau.tauID("byVVVLooseDeepTau2018v2p5VSjet")));
}

template<class T>
bool zToLepProbeTrk<T>::passesLightFlavorRejection (const pat::Tau &tau){
  return (tau.tauID("byVVVLooseDeepTau2017v2p1VSe") || tau.tauID("byVVVLooseDeepTau2018v2p5VSe")) && (tau.tauID("byVLooseDeepTau2017v2p1VSmu") || tau.tauID("byVLooseDeepTau2018v2p5VSmu"));
}

template<class T>
bool zToLepProbeTrk<T>::inTOBCrack (const pat::IsolatedTrack &track){
  return (fabs(track.dz()) < 0.5 && fabs(1.57079632679489661923 - track.theta()) < 1.0e-3);
}

template<class T> 
const int zToLepProbeTrk<T>::extraMissingMiddleHits (const pat::IsolatedTrack &track) const
{

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> distribution (0.0, 1.0);
  double dropTOBProbability_ = 0.00830971251971; // This is taken from Run 2; needs to be updated once correct value is estimated
  double hitProbability_ = 0.0175874821487; // This is taken from Run 2; needs to be updated once correct value is estimated
  // bool dropHits = true;
  bool dropHits = false;
  bool dropTOBDecision_ = (dropHits ? distribution (generator) : 1.0e6) < dropTOBProbability_;
  std::vector<bool> dropMiddleHitDecisions_;
  for (int i = 0; i < 50; i++)
    dropMiddleHitDecisions_.push_back ((dropHits ? distribution (generator) : 1.0e6) < hitProbability_);

  int nHits = 0;
  bool countMissingMiddleHits = false;
  for (int i = 0; i < track.hitPattern().stripLayersWithMeasurement() - (dropTOBDecision_ ? track.hitPattern ().stripTOBLayersWithMeasurement () : 0); i++)
    {
      bool hit = !dropMiddleHitDecisions_.at(i);
      if (!hit && countMissingMiddleHits)
        nHits++;
      if (hit)
        countMissingMiddleHits = true;
    }

  return nHits;
}

template<class T>
const int zToLepProbeTrk<T>::hitDrop_missingMiddleHits (const pat::IsolatedTrack &track) const
{
  int nDropHits = extraMissingMiddleHits(track);
  return track.hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS) + nDropHits;
}

template<class T>
double zToLepProbeTrk<T>::dRMinJet (const pat::IsolatedTrack &track, const std::vector<pat::Jet> &jets)
{
  double dRMinJet = 99.0;
  for(const auto &jet : jets){
    if(jet.pt() > 30 &&
    fabs(jet.eta()) < 4.5 &&
    (((jet.neutralHadronEnergyFraction()<0.90 &&
    jet.neutralEmEnergyFraction()<0.90 &&
    (jet.chargedMultiplicity() + jet.neutralMultiplicity())>1 &&
    jet.muonEnergyFraction()<0.8) &&
    ((fabs(jet.eta())<=2.4 &&
    jet.chargedHadronEnergyFraction()>0 &&
    jet.chargedMultiplicity()>0 &&
    jet.chargedEmEnergyFraction()<0.90) ||
    fabs(jet.eta())>2.4) &&
    fabs(jet.eta())<=3.0) ||
    (jet.neutralEmEnergyFraction()<0.90 && jet.neutralMultiplicity()>10 && fabs(jet.eta())>3.0)))
    {
      double dR = deltaR(track, jet);
      if(dR < dRMinJet || dRMinJet < 0.0) dRMinJet = dR;
    }
  }
  return dRMinJet;
}

template<class T>
double zToLepProbeTrk<T>::deltaRToClosestMuon(const pat::IsolatedTrack &track, const std::vector<pat::Muon> &muons) 
{

  double dR;
  double deltaRToClosestMuon = 999999.9;

  for(const auto &muon : muons) {
    dR = deltaR(track, muon);
    if(dR < deltaRToClosestMuon || deltaRToClosestMuon < 0.0) deltaRToClosestMuon = dR;
  }

  return deltaRToClosestMuon;

}

template<class T>
double zToLepProbeTrk<T>::deltaRToClosestTauHad(const pat::IsolatedTrack &track, const std::vector<pat::Tau> &taus) 
{

  double dR;
  double deltaRToClosestTauHad = 999999.9;
  bool passesDecayModeReconstruction;
  bool passesLightFlavorRejection;

  for(const auto &tau : taus) {
    dR = deltaR(track, tau);

    passesDecayModeReconstruction = (tau.tauID("decayModeFindingNewDMs"));
    passesLightFlavorRejection = (tau.tauID("byVVVLooseDeepTau2017v2p1VSe") || tau.tauID("byVVVLooseDeepTau2018v2p5VSe")) && (tau.tauID("byVLooseDeepTau2017v2p1VSmu") || tau.tauID("byVLooseDeepTau2018v2p5VSmu"));

    if(passesDecayModeReconstruction && passesLightFlavorRejection && (dR < deltaRToClosestTauHad  || deltaRToClosestTauHad  < 0.0)) {
      deltaRToClosestTauHad = dR;
    }
  }

  return deltaRToClosestTauHad;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<class T>
void zToLepProbeTrk<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("electrons", edm::InputTag("slimmedElectrons"));
  desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("taus", edm::InputTag("slimmedTaus"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("isolatedTracks"));
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("triggersPAT", edm::InputTag("TriggerResults","","PAT"));
  desc.add<edm::InputTag>("triggersHLT", edm::InputTag("TriggerResults","","HLT"));
  desc.add<std::string>("HLTName", std::string("placeholderHLT"));
  desc.add<edm::InputTag>("trigobjs", edm::InputTag("slimmedPatTrigger"));

  descriptions.addWithDefaultLabel(desc);

}

using zToElecProbeTrk = zToLepProbeTrk<pat::Electron>;
using zToMuonProbeTrk = zToLepProbeTrk<pat::Muon>;
using zToTauProbeTrk = zToLepProbeTrk<pat::Tau>;

//define this as a plug-in
DEFINE_FWK_MODULE(zToElecProbeTrk);
DEFINE_FWK_MODULE(zToMuonProbeTrk);
DEFINE_FWK_MODULE(zToTauProbeTrk);
