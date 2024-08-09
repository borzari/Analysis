#ifndef Analysis_Helper_interface_plotPrintFunctions_h
#define Analysis_Helper_interface_plotPrintFunctions_h

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
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector2.h"

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

struct VariablesToPlot {

  std::vector<pat::Electron> electrons;
  std::vector<pat::Muon> muons;
  std::vector<pat::Tau> taus;
  std::vector<pat::IsolatedTrack> tracks;
  std::vector<pat::Jet> jets;
  pat::MET met;
  double DeltaPhiJetMET;
  TVector2 metNoMu;
  TVector2 metNoMuNoLep;

};

class plotPrintFunctions {

  public:

  static void printCuts(std::map<std::string, TH1D *> &oneDHists){

    int maxSizeCutName = 0;
    int maxSizeEvents = 0;
    int maxSizeCutFlow = 0;
    int maxSizeSelection = 0;
    double maxValue = double(oneDHists.at("selection")->GetBinContent(1));

    for(int i = 1; i <= int(oneDHists.at("cutflow")->GetNbinsX()); i++){
      std::string label = oneDHists.at("cutflow")->GetXaxis()->GetLabels()->At(i-1)->GetName();
      std::string numEvents = std::to_string(oneDHists.at("cutflow")->GetBinContent(i));
      std::string numCutflow = std::to_string(double(oneDHists.at("cutflow")->GetBinContent(i)) * 100.0 / maxValue) + "%";
      std::string numSelection = std::to_string(double(oneDHists.at("selection")->GetBinContent(i)) * 100.0 / maxValue) + "%";
      if(int(label.size()) > maxSizeCutName) maxSizeCutName = int(label.size());
      if(int(numEvents.size()) > maxSizeEvents) maxSizeEvents = int(numEvents.size());
      if(int(numCutflow.size()) > maxSizeCutFlow) maxSizeCutFlow = int(numCutflow.size());
      if(int(numSelection.size()) > maxSizeSelection) maxSizeSelection = int(numSelection.size());
    }

    std::string cut = "Cut Name";
    std::string events = "Events";
    std::string cumul = "Cumul.";
    std::string indiv = "Indiv.";
    std::cout << cut << std::string((maxSizeCutName - int(cut.size()) + 3),' ') << events << std::string((maxSizeEvents - int(events.size()) + 3),' ') << cumul << std::string((maxSizeCutFlow - int(cumul.size()) + 3),' ') << indiv << std::endl;

    for(int i = 1; i <= int(oneDHists.at("cutflow")->GetNbinsX()); i++){
      std::cout << std::setprecision(3) << std::fixed;
      std::string label = oneDHists.at("cutflow")->GetXaxis()->GetLabels()->At(i-1)->GetName();
      std::string numEvents = std::to_string(oneDHists.at("cutflow")->GetBinContent(i));
      std::string numCutflow = std::to_string(double(oneDHists.at("cutflow")->GetBinContent(i)) * 100.0 / maxValue) + "%";
      std::string numSelection = std::to_string(double(oneDHists.at("selection")->GetBinContent(i)) * 100.0 / maxValue) + "%";
      int sizeLabel = int(label.size());
      int sizeEvents = int(numEvents.size());
      int sizeCutflow = int(numCutflow.size());
      int sizeSelection = int(numSelection.size());
      int rightPrintEvents = maxSizeEvents - sizeEvents;
      int rightPrintCutflow = maxSizeCutFlow - sizeCutflow;
      int rightPrintSelection = maxSizeSelection - sizeSelection;
      int nBlankLabel = maxSizeCutName - sizeLabel + rightPrintEvents + 3;
      int nBlankEvents = maxSizeEvents - sizeEvents - rightPrintEvents + rightPrintCutflow + 8; // +8 here because 5 places are going to be removed below
      int nBlankCutflow = maxSizeCutFlow - sizeCutflow - rightPrintCutflow + rightPrintSelection + 6; // +6 here because 3 places are going to be removed below
      numEvents.erase(numEvents.size() - 5, 5);
      numCutflow.erase(numCutflow.size() - 4, 3);
      numSelection.erase(numSelection.size() - 4, 3);
      std::cout << label << std::string(nBlankLabel,' ') << numEvents << std::string(nBlankEvents,' ') << numCutflow << std::string(nBlankCutflow,' ') << numSelection << std::endl;
    }

  }

  template<class T>
  static void createObjHists(edm::Service<TFileService> &fs, std::map<std::string, TH1D *> &oneDHists, std::string const &objName){

    std::string pt = objName + "_pt";
    std::string eta = objName + "_eta";
    std::string phi = objName + "_phi";
    std::string charge = objName + "_charge";

    oneDHists[objName + "_pt"] = fs->make<TH1D>(pt.c_str(), "", 200, 0.0, 1000.0);
    oneDHists[objName + "_eta"] = fs->make<TH1D>(eta.c_str(), "", 64, 0.0, 3.2);
    oneDHists[objName + "_phi"] = fs->make<TH1D>(phi.c_str(), "", 64, 0.0, 3.2);
    oneDHists[objName + "_charge"] = fs->make<TH1D>(charge.c_str(), "", 3, -1.5, 1.5);
  
  }

  static void createCommonHists(edm::Service<TFileService> &fs, std::map<std::string, TH1D *> &oneDHists, std::map<std::string, TH2D *> &twoDHists){

    oneDHists["met_pt"] = fs->make<TH1D>("met_pt", "", 200, 0.0, 1000.0);
    oneDHists["met_phi"] = fs->make<TH1D>("met_phi", "", 64, 0.0, 3.2);
    oneDHists["metNoMu_pt"] = fs->make<TH1D>("metNoMu_pt", "", 200, 0.0, 1000.0);
    oneDHists["metNoMu_phi"] = fs->make<TH1D>("metNoMu_phi", "", 64, 0.0, 3.2);
    oneDHists["metNoMuNoLep_pt"] = fs->make<TH1D>("metNoMuNoLep_pt", "", 200, 0.0, 1000.0);
    oneDHists["metNoMuNoLep_phi"] = fs->make<TH1D>("metNoMuNoLep_phi", "", 64, 0.0, 3.2);

    twoDHists["metNoMuvsDeltaPhiJetMetNoMu"] = fs->make<TH2D>("metNoMuvsDeltaPhiJetMetNoMu", "", 200, 0.0, 1000.0, 64, 0.0, 3.2);
  
  }

  template<class T>
  static void plotObjs(std::map<std::string, TH1D *> &oneDHists, std::vector<T> const &objs, std::string objName){

    for (const auto& obj : objs) {
      
      oneDHists.at(objName + "_pt")->Fill(obj.pt());
      oneDHists.at(objName + "_eta")->Fill(obj.eta());
      oneDHists.at(objName + "_phi")->Fill(obj.phi());
      oneDHists.at(objName + "_charge")->Fill(obj.charge());
  
    }

  }

  static void plotVariables(std::map<std::string, TH1D *> oneDHists, std::map<std::string, TH2D *> twoDHists, VariablesToPlot const &objs){

    if(objs.electrons.size() > 0) plotObjs<pat::Electron>(oneDHists,objs.electrons,"electron");
    if(objs.muons.size() > 0) plotObjs<pat::Muon>(oneDHists,objs.muons,"muon");
    if(objs.taus.size() > 0) plotObjs<pat::Tau>(oneDHists,objs.taus,"tau");
    if(objs.tracks.size() > 0) plotObjs<pat::IsolatedTrack>(oneDHists,objs.tracks,"track");
    if(objs.jets.size() > 0) plotObjs<pat::Jet>(oneDHists,objs.jets,"jet");

    oneDHists.at("met_pt")->Fill(objs.met.pt());
    oneDHists.at("met_phi")->Fill(objs.met.phi());
    oneDHists.at("metNoMu_pt")->Fill(objs.metNoMu.Mod());
    oneDHists.at("metNoMu_phi")->Fill(objs.metNoMu.Phi());
    oneDHists.at("metNoMuNoLep_pt")->Fill(objs.metNoMuNoLep.Mod());
    oneDHists.at("metNoMuNoLep_phi")->Fill(objs.metNoMuNoLep.Phi());

    twoDHists.at("metNoMuvsDeltaPhiJetMetNoMu")->Fill(objs.metNoMuNoLep.Mod(),objs.DeltaPhiJetMET);

  }

};

#endif // Analysis_Helper_interface_plotPrintFunctions_h