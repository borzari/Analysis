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
class leptonTagSkim : public edm::one::EDFilter<edm::one::SharedResources> {
public:
  explicit leptonTagSkim(const edm::ParameterSet&);
  ~leptonTagSkim() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  bool isMatchedToTriggerObject (const edm::Event &, const edm::TriggerResults &, const T &, const std::vector<pat::TriggerObjectStandAlone> &, const std::string &, const std::string &, const double = 0.1);
  bool passesDecayModeReconstruction (const pat::Tau &);
  bool passesLightFlavorRejection (const pat::Tau &);


private:
  bool filter(edm::Event&, const edm::EventSetup&) override;
  edm::Service<TFileService> fs_;
  std::map<std::string, TH1D *> oneDHists_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken_;
  edm::EDGetTokenT<std::vector<T>> leptonsToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggersPATToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggersHLTToken_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> trigobjsToken_;
  std::string HLTName_;

  bool passHLTSel;
  bool passMETFiltersSel;
  bool passPtSel;
  bool passEtaSel;
  bool passIDSel;
  bool passD0Sel;
  bool passDzSel;
  bool passPFIsoSel;

  bool passHLTCut;
  bool passMETFiltersCut;
  bool passPtCut;
  bool passEtaCut;
  bool passIDCut;
  bool passD0Cut;
  bool passDzCut;
  bool passPFIsoCut;

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
leptonTagSkim<T>::leptonTagSkim(const edm::ParameterSet& iConfig)
    : verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
      leptonsToken_(consumes<std::vector<T>>(iConfig.getParameter<edm::InputTag>("leptons"))),
      triggersPATToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggersPAT"))),
      triggersHLTToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggersHLT"))),
      trigobjsToken_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("trigobjs"))) {
  //now do what ever initialization is needed

  HLTName_ = iConfig.getParameter<std::string>("HLTName");

  int nBins = 0;
  double maxRange = 0.0;
  std::string ptCut = "";
  std::string idCut = "";

  if(std::is_same<T, pat::Electron>::value){nBins = 8; maxRange = 7.5; ptCut = "pt > 32"; idCut = "passesVID_tightID (ID + iso)";}
  if(std::is_same<T, pat::Muon>::value){nBins = 7; maxRange = 6.5; ptCut = "pt > 26"; idCut = "isTightMuonWRTVtx";}
  if(std::is_same<T, pat::Tau>::value){nBins = 6; maxRange = 5.5; ptCut = "pt > 50"; idCut = "passesDecayModeReconstruction && passesLightFlavorRejection";}

  oneDHists_["cutflow"] = fs_->make<TH1D>("cutflow", "", nBins, -0.5, maxRange);
  oneDHists_["selection"] = fs_->make<TH1D>("selection", "", nBins, -0.5, maxRange);

  oneDHists_.at("cutflow")->GetXaxis()->SetBinLabel(1,"Total");
  oneDHists_.at("cutflow")->GetXaxis()->SetBinLabel(2,HLTName_.c_str());
  oneDHists_.at("cutflow")->GetXaxis()->SetBinLabel(3,"METFilters");
  oneDHists_.at("cutflow")->GetXaxis()->SetBinLabel(4,ptCut.c_str());
  oneDHists_.at("cutflow")->GetXaxis()->SetBinLabel(5,"fabs(#eta) < 2.1");
  oneDHists_.at("cutflow")->GetXaxis()->SetBinLabel(6,idCut.c_str());
  if(std::is_same<T, pat::Muon>::value) oneDHists_.at("cutflow")->GetXaxis()->SetBinLabel(7,"#Delta#beta-corrected rel. iso. < 0.15");
  if(std::is_same<T, pat::Electron>::value){
    oneDHists_.at("cutflow")->GetXaxis()->SetBinLabel(7,"|d0| < 0.05, 0.10 (EB, EE)");
    oneDHists_.at("cutflow")->GetXaxis()->SetBinLabel(8,"|dz| < 0.10, 0.20 (EB, EE)");
  }

  oneDHists_.at("selection")->GetXaxis()->SetBinLabel(1,"Total");
  oneDHists_.at("selection")->GetXaxis()->SetBinLabel(2,HLTName_.c_str());
  oneDHists_.at("selection")->GetXaxis()->SetBinLabel(3,"METFilters");
  oneDHists_.at("selection")->GetXaxis()->SetBinLabel(4,ptCut.c_str());
  oneDHists_.at("selection")->GetXaxis()->SetBinLabel(5,"fabs(#eta) < 2.1");
  oneDHists_.at("selection")->GetXaxis()->SetBinLabel(6,idCut.c_str());
  if(std::is_same<T, pat::Muon>::value) oneDHists_.at("selection")->GetXaxis()->SetBinLabel(7,"#Delta#beta-corrected rel. iso. < 0.15");
  if(std::is_same<T, pat::Electron>::value){
    oneDHists_.at("selection")->GetXaxis()->SetBinLabel(7,"|d0| < 0.05, 0.10 (EB, EE)");
    oneDHists_.at("selection")->GetXaxis()->SetBinLabel(8,"|dz| < 0.10, 0.20 (EB, EE)");
  }

}

template<class T>
leptonTagSkim<T>::~leptonTagSkim() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
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
    std::cout << oneDHists_.at("cutflow")->GetXaxis()->GetLabels()->At(i-1)->GetName() << std::string(nBlank,' ') << "\t\t" << double(oneDHists_.at("cutflow")->GetBinContent(i)) * 100.0 / maxValue << "% \t" << double(oneDHists_.at("selection")->GetBinContent(i)) * 100.0 / maxValue << "%" << std::endl;
  }
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
template<class T>
bool leptonTagSkim<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  passHLTSel = false;
  passMETFiltersSel = false;
  passPtSel = false;
  passEtaSel = false;
  passIDSel = false;
  passD0Sel = false;
  passDzSel = false;
  passPFIsoSel = false;

  passHLTCut = false;
  passMETFiltersCut = false;
  passPtCut = false;
  passEtaCut = false;
  passIDCut = false;
  passD0Cut = false;
  passDzCut = false;
  passPFIsoCut = false;

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

  edm::Handle<std::vector<pat::Electron>> electrons;
  if(std::is_same<T, pat::Electron>::value) iEvent.getByToken(leptonsToken_, electrons);  

  edm::Handle<std::vector<pat::Muon>> muons;
  if(std::is_same<T, pat::Muon>::value) iEvent.getByToken(leptonsToken_, muons);  

  edm::Handle<std::vector<pat::Tau>> taus;
  if(std::is_same<T, pat::Tau>::value) iEvent.getByToken(leptonsToken_, taus);

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

  if(pass_HLT == 1) {passHLTSel = true; passHLTCut = true;}

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

  if(pass_Flag_goodVertices == 1 && pass_Flag_globalTightHalo2016Filter == 1 && pass_Flag_HBHENoiseFilter == 1 && pass_Flag_HBHENoiseIsoFilter == 1 && pass_Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && pass_Flag_BadPFMuonFilter == 1 && pass_Flag_BadChargedCandidateFilter == 1 && pass_Flag_ecalBadCalibFilter == 1) {passMETFiltersSel = true; if(passHLTCut) passMETFiltersCut = true;}

  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(verticesToken_, vertices);
  const reco::Vertex &pv = vertices->at(0);

  bool auxPassPtCut = false;
  bool auxPassEtaCut = false;
  bool auxPassIDCut = false;
  bool auxPassD0Cut = false;
  bool auxPassDzCut = false;
  bool auxPassPFIsoCut = false;
  
  if(std::is_same<T, pat::Electron>::value){
    
    for (const auto& electron : *electrons) {
      
      if(electron.pt() > 32.) {passPtSel = true; if(passMETFiltersCut) auxPassPtCut = true;}
      if(abs(electron.eta()) < 2.1) {passEtaSel = true; if(auxPassPtCut) auxPassEtaCut = true;}
      if(electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-tight")) {passIDSel = true; if(auxPassEtaCut) auxPassIDCut = true;}
      if(((fabs (electron.superCluster()->eta()) <= 1.479) && (fabs (((electron.vx() - pv.x()) * electron.py() - (electron.vy() - pv.y()) * electron.px()) / electron.pt()) < 0.05)) || ((fabs (electron.superCluster()->eta()) >  1.479) && (fabs (((electron.vx() - pv.x()) * electron.py() - (electron.vy() - pv.y()) * electron.px()) / electron.pt()) < 0.10))) {passD0Sel = true; if(auxPassIDCut) auxPassD0Cut = true;}
      if(((fabs (electron.superCluster()->eta()) <= 1.479) && (fabs ((electron.vz() - pv.z()) - ((electron.vx() - pv.x()) * electron.px() + (electron.vy() - pv.y()) * electron.py()) / electron.pt() * electron.pz() / electron.pt()) < 0.10)) || ((fabs (electron.superCluster()->eta()) >  1.479) && (fabs ((electron.vz() - pv.z()) - ((electron.vx() - pv.x()) * electron.px() + (electron.vy() - pv.y()) * electron.py()) / electron.pt() * electron.pz() / electron.pt()) < 0.20))) {passDzSel = true; if(auxPassD0Cut) auxPassDzCut = true;}
      
      if(auxPassDzCut) passDzCut = true;
      if(auxPassD0Cut) passD0Cut = true;
      if(auxPassIDCut) passIDCut = true;
      if(auxPassEtaCut) passEtaCut = true;
      if(auxPassPtCut) passPtCut = true;
          
      auxPassPtCut = false;
      auxPassEtaCut = false;
      auxPassIDCut = false;
      auxPassD0Cut = false;
      auxPassDzCut = false;
        
    }
  
  }

  if(std::is_same<T, pat::Muon>::value){
    
    for (const auto& muon : *muons) {
      
      if(muon.pt() > 26.) {passPtSel = true; if(passMETFiltersCut) auxPassPtCut = true;}
      if(abs(muon.eta()) < 2.1) {passEtaSel = true; if(auxPassPtCut) auxPassEtaCut = true;}
      if(muon.isTightMuon(pv)) {passIDSel = true; if(auxPassEtaCut) auxPassIDCut = true;}
      if(((muon.pfIsolationR04().sumChargedHadronPt + std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt()) < 0.15) {passPFIsoSel = true; if(auxPassIDCut) auxPassPFIsoCut = true;}
      
      if(auxPassPFIsoCut) passPFIsoCut = true;
      if(auxPassIDCut) passIDCut = true;
      if(auxPassEtaCut) passEtaCut = true;
      if(auxPassPtCut) passPtCut = true;
          
      auxPassPtCut = false;
      auxPassEtaCut = false;
      auxPassIDCut = false;
      auxPassPFIsoCut = false;
        
    }
  
  }

  if(std::is_same<T, pat::Tau>::value){
    
    for (const auto& tau : *taus) {
      
      if(tau.pt() > 50.) {passPtSel = true; if(passMETFiltersCut) auxPassPtCut = true;}
      if(abs(tau.eta()) < 2.1) {passEtaSel = true; if(auxPassPtCut) auxPassEtaCut = true;}
      if(passesDecayModeReconstruction(tau) && passesLightFlavorRejection(tau)) {passIDSel = true; if(auxPassEtaCut) auxPassIDCut = true;}
      
      if(auxPassIDCut) passIDCut = true;
      if(auxPassEtaCut) passEtaCut = true;
      if(auxPassPtCut) passPtCut = true;
          
      auxPassPtCut = false;
      auxPassEtaCut = false;
      auxPassIDCut = false;
        
    }
  
  }

  if(passHLTSel) oneDHists_.at("selection")->Fill(1);
  if(passMETFiltersSel) oneDHists_.at("selection")->Fill(2);
  if(passPtSel) oneDHists_.at("selection")->Fill(3);
  if(passEtaSel) oneDHists_.at("selection")->Fill(4);
  if(passIDSel) oneDHists_.at("selection")->Fill(5);
  if(std::is_same<T, pat::Muon>::value){if(passPFIsoSel) oneDHists_.at("selection")->Fill(6);}
  if(std::is_same<T, pat::Electron>::value){
    if(passD0Sel) oneDHists_.at("selection")->Fill(6);
    if(passDzSel) oneDHists_.at("selection")->Fill(7);
  }

  if(passHLTCut) oneDHists_.at("cutflow")->Fill(1);
  if(passMETFiltersCut) oneDHists_.at("cutflow")->Fill(2);
  if(passPtCut) oneDHists_.at("cutflow")->Fill(3);
  if(passEtaCut) oneDHists_.at("cutflow")->Fill(4);
  if(passIDCut) {oneDHists_.at("cutflow")->Fill(5); if(std::is_same<T, pat::Tau>::value) isGood = true;}
  if(std::is_same<T, pat::Muon>::value){if(passPFIsoCut) {oneDHists_.at("cutflow")->Fill(6); isGood = true;}}
  if(std::is_same<T, pat::Electron>::value){
    if(passD0Cut) oneDHists_.at("cutflow")->Fill(6);
    if(passDzCut) {oneDHists_.at("cutflow")->Fill(7); isGood = true;}
  }

  return isGood;
}

template<class T>
bool leptonTagSkim<T>::isMatchedToTriggerObject (const edm::Event &event, const edm::TriggerResults &triggers, const T &obj, const std::vector<pat::TriggerObjectStandAlone> &trigObjs, const std::string &collection, const std::string &filter, const double dR)
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
    if(filter == "hltTrk50Filter") std::cout << trigObj.eta() << std::endl;
    if(deltaR (obj, trigObj) > dR) continue;
    return true;
  }
  return false;
}

template<class T>
bool leptonTagSkim<T>::passesDecayModeReconstruction (const pat::Tau &tau){
  return (tau.tauID("decayModeFindingNewDMs") && (tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") || (tau.tauID("chargedIsoPtSumdR03")+std::max(0.,tau.tauID("neutralIsoPtSumdR03")-0.072*tau.tauID("puCorrPtSum"))<2.5) || tau.tauID("byVVVLooseDeepTau2017v2p1VSjet") || tau.tauID("byVVVLooseDeepTau2018v2p5VSjet")));
}

template<class T>
bool leptonTagSkim<T>::passesLightFlavorRejection (const pat::Tau &tau){
  return (tau.tauID("byVVVLooseDeepTau2017v2p1VSe") || tau.tauID("byVVVLooseDeepTau2018v2p5VSe")) && (tau.tauID("byVLooseDeepTau2017v2p1VSmu") || tau.tauID("byVLooseDeepTau2018v2p5VSmu"));
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<class T>
void leptonTagSkim<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("leptons", edm::InputTag("placeholderLeptons"));
  desc.add<edm::InputTag>("triggersPAT", edm::InputTag("TriggerResults","","PAT"));
  desc.add<edm::InputTag>("triggersHLT", edm::InputTag("TriggerResults","","HLT"));
  desc.add<std::string>("HLTName", std::string("placeholderHLT"));
  desc.add<edm::InputTag>("trigobjs", edm::InputTag("slimmedPatTrigger"));

  descriptions.addWithDefaultLabel(desc);

}

using electronTagSkim = leptonTagSkim<pat::Electron>;
using muonTagSkim = leptonTagSkim<pat::Muon>;
using tauTagSkim = leptonTagSkim<pat::Tau>;

//define this as a plug-in
// DEFINE_FWK_MODULE(leptonTagSkim);
DEFINE_FWK_MODULE(electronTagSkim);
DEFINE_FWK_MODULE(muonTagSkim);
DEFINE_FWK_MODULE(tauTagSkim);
