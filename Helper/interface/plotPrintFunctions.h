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
// #include "TFileDirectory.h"

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

struct nTPPairs {

  int nTPOS;
  int nTPSS;
  int nTPOS_veto;
  int nTPSS_veto;
  int nTPOS_loose_veto;
  int nTPSS_loose_veto;

};

struct VariablesToPlot {

  std::vector<pat::Electron> electrons;
  std::vector<pat::Muon> muons;
  std::vector<pat::Tau> taus;
  std::vector<pat::IsolatedTrack> tracks;
  std::vector<pat::Jet> jets;
  pat::MET met;
  double DeltaPhiJetMET;
  double DeltaPhiJetMETNoMu;
  double DeltaPhiJetMETNoMuNoLep;
  std::vector<double> caloNewNoPUDRp5CentralCalo;
  TVector2 metNoMu;
  TVector2 metNoMuNoLep;
  nTPPairs nTP;

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
  static void createObjHists(edm::Service<TFileService> &fs, std::map<std::string, TH1D *> &oneDHists, std::string const &objName, TFileDirectory &dir, const int NLayers){

    std::string pt = objName + "_pt";
    std::string eta = objName + "_eta";
    std::string phi = objName + "_phi";
    std::string charge = objName + "_charge";
    std::string calo = objName + "_caloNewNoPUDRp5CentralCalo";
    std::string matchedjet = objName + "_matchedJetCaloEnergy";

    oneDHists[std::to_string(NLayers) + objName + "_pt"] = dir.make<TH1D>(pt.c_str(), "", 200, 0.0, 1000.0);
    oneDHists[std::to_string(NLayers) + objName + "_eta"] = dir.make<TH1D>(eta.c_str(), "", 128, -3.2, 3.2);
    oneDHists[std::to_string(NLayers) + objName + "_phi"] = dir.make<TH1D>(phi.c_str(), "", 128, -3.2, 3.2);
    oneDHists[std::to_string(NLayers) + objName + "_charge"] = dir.make<TH1D>(charge.c_str(), "", 3, -1.5, 1.5);
    if(std::is_same<T, pat::IsolatedTrack>::value){
      oneDHists[std::to_string(NLayers) + objName + "_caloNewNoPUDRp5CentralCalo"] = dir.make<TH1D>(calo.c_str(), "", 200, 0.0, 1000.0);
      oneDHists[std::to_string(NLayers) + objName + "_matchedJetCaloEnergy"] = dir.make<TH1D>(matchedjet.c_str(), "", 200, 0.0, 1000.0);
    };
  
  }

  static void createCommonHists(edm::Service<TFileService> &fs, std::map<std::string, TH1D *> &oneDHists, std::map<std::string, TH2D *> &twoDHists, TFileDirectory &dir, const int NLayers){

    oneDHists[std::to_string(NLayers) + "met_pt"] = dir.make<TH1D>("met_pt", "", 200, 0.0, 1000.0);
    oneDHists[std::to_string(NLayers) + "met_phi"] = dir.make<TH1D>("met_phi", "", 128, -3.2, 3.2);
    oneDHists[std::to_string(NLayers) + "metNoMu_pt"] = dir.make<TH1D>("metNoMu_pt", "", 200, 0.0, 1000.0);
    oneDHists[std::to_string(NLayers) + "metNoMu_phi"] = dir.make<TH1D>("metNoMu_phi", "", 128, -3.2, 3.2);
    oneDHists[std::to_string(NLayers) + "metNoMuNoLep_pt"] = dir.make<TH1D>("metNoMuNoLep_pt", "", 200, 0.0, 1000.0);
    oneDHists[std::to_string(NLayers) + "metNoMuNoLep_phi"] = dir.make<TH1D>("metNoMuNoLep_phi", "", 128, -3.2, 3.2);

    twoDHists[std::to_string(NLayers) + "metNoMuvsDeltaPhiJetMet"] = dir.make<TH2D>("metNoMuvsDeltaPhiJetMet", "", 200, 0.0, 1000.0, 128, -3.2, 3.2);
    twoDHists[std::to_string(NLayers) + "metNoMuvsDeltaPhiJetMetNoMu"] = dir.make<TH2D>("metNoMuvsDeltaPhiJetMetNoMu", "", 200, 0.0, 1000.0, 128, -3.2, 3.2);
    twoDHists[std::to_string(NLayers) + "metNoMuvsDeltaPhiJetMetNoMuNoLep"] = dir.make<TH2D>("metNoMuvsDeltaPhiJetMetNoMuNoLep", "", 200, 0.0, 1000.0, 128, -3.2, 3.2);
  
  }

  static void createTPPairsHists(edm::Service<TFileService> &fs, std::map<std::string, TH1D *> &oneDHists, TFileDirectory &dir, const int NLayers){

    oneDHists[std::to_string(NLayers) + "nTPOS"] = dir.make<TH1D>("nTPOS", "", 10, -0.5, 9.5);
    oneDHists[std::to_string(NLayers) + "nTPSS"] = dir.make<TH1D>("nTPSS", "", 10, -0.5, 9.5);
    oneDHists[std::to_string(NLayers) + "nTPOS_veto"] = dir.make<TH1D>("nTPOS_veto", "", 10, -0.5, 9.5);
    oneDHists[std::to_string(NLayers) + "nTPSS_veto"] = dir.make<TH1D>("nTPSS_veto", "", 10, -0.5, 9.5);
    oneDHists[std::to_string(NLayers) + "nTPOSLoose_veto"] = dir.make<TH1D>("nTPOSLoose_veto", "", 10, -0.5, 9.5);
    oneDHists[std::to_string(NLayers) + "nTPSSLoose_veto"] = dir.make<TH1D>("nTPSSLoose_veto", "", 10, -0.5, 9.5);

  }

  template<class T>
  static void plotObjs(std::map<std::string, TH1D *> &oneDHists, std::vector<T> const &objs, std::string objName, const int NLayers, const double sf){

    for (const auto& obj : objs) {
      
      oneDHists.at(std::to_string(NLayers) + objName + "_pt")->Fill(obj.pt(),sf);
      oneDHists.at(std::to_string(NLayers) + objName + "_eta")->Fill(obj.eta(),sf);
      oneDHists.at(std::to_string(NLayers) + objName + "_phi")->Fill(obj.phi(),sf);
      oneDHists.at(std::to_string(NLayers) + objName + "_charge")->Fill(obj.charge(),sf);

    }

  }

  static void plotTrackCaloEnergy(std::map<std::string, TH1D *> &oneDHists, std::vector<pat::IsolatedTrack> const &objs, std::vector<double> const &caloNewNoPUDRp5CentralCalo, std::string objName, const int NLayers, const double sf){
    for (const auto& obj : objs) oneDHists.at(std::to_string(NLayers) + objName + "_matchedJetCaloEnergy")->Fill((obj.matchedCaloJetEmEnergy() + obj.matchedCaloJetHadEnergy()),sf);
    for (const auto& obj : caloNewNoPUDRp5CentralCalo) oneDHists.at(std::to_string(NLayers) + objName + "_caloNewNoPUDRp5CentralCalo")->Fill(obj,sf);
  }

  static void plotVariables(std::map<std::string, TH1D *> oneDHists, std::map<std::string, TH2D *> twoDHists, VariablesToPlot const &objs, const int NLayers, const double sf){

    if(NLayers == 4){
      if(objs.tracks.size() > 0) {
        if(objs.tracks[0].hitPattern().trackerLayersWithMeasurement() != 4) return;
      }
    }
    if(NLayers == 5){
      if(objs.tracks.size() > 0) {
        if(objs.tracks[0].hitPattern().trackerLayersWithMeasurement() != 5) return;
      }
    }
    if(NLayers == 6){
      if(objs.tracks.size() > 0) {
        if(objs.tracks[0].hitPattern().trackerLayersWithMeasurement() < 6) return;
      }
    }

    if(objs.electrons.size() > 0) plotObjs<pat::Electron>(oneDHists,objs.electrons,"electron",NLayers,sf);
    if(objs.muons.size() > 0) plotObjs<pat::Muon>(oneDHists,objs.muons,"muon",NLayers,sf);
    if(objs.taus.size() > 0) plotObjs<pat::Tau>(oneDHists,objs.taus,"tau",NLayers,sf);
    if(objs.tracks.size() > 0) {
      plotObjs<pat::IsolatedTrack>(oneDHists,objs.tracks,"track",NLayers,sf);
      plotTrackCaloEnergy(oneDHists,objs.tracks,objs.caloNewNoPUDRp5CentralCalo,"track",NLayers,sf);
    }
    if(objs.jets.size() > 0) plotObjs<pat::Jet>(oneDHists,objs.jets,"jet",NLayers,sf);

    oneDHists.at(std::to_string(NLayers) + "met_pt")->Fill(objs.met.pt(),sf);
    oneDHists.at(std::to_string(NLayers) + "met_phi")->Fill(objs.met.phi(),sf);
    oneDHists.at(std::to_string(NLayers) + "metNoMu_pt")->Fill(objs.metNoMu.Mod(),sf);
    oneDHists.at(std::to_string(NLayers) + "metNoMu_phi")->Fill(objs.metNoMu.Phi_mpi_pi(objs.metNoMu.Phi()),sf);
    oneDHists.at(std::to_string(NLayers) + "metNoMuNoLep_pt")->Fill(objs.metNoMuNoLep.Mod(),sf);
    oneDHists.at(std::to_string(NLayers) + "metNoMuNoLep_phi")->Fill(objs.metNoMuNoLep.Phi_mpi_pi(objs.metNoMuNoLep.Phi()),sf);

    twoDHists.at(std::to_string(NLayers) + "metNoMuvsDeltaPhiJetMet")->Fill(objs.met.pt(),objs.DeltaPhiJetMET,sf);
    twoDHists.at(std::to_string(NLayers) + "metNoMuvsDeltaPhiJetMetNoMu")->Fill(objs.metNoMu.Mod(),objs.DeltaPhiJetMETNoMu,sf);
    twoDHists.at(std::to_string(NLayers) + "metNoMuvsDeltaPhiJetMetNoMuNoLep")->Fill(objs.metNoMuNoLep.Mod(),objs.DeltaPhiJetMETNoMuNoLep,sf);

  }

  static void plotTPPairs(std::map<std::string, TH1D *> oneDHists, VariablesToPlot const &objs, const int NLayers, const double sf){

    if(NLayers == 4){
      if(objs.tracks.size() > 0) {
        if(objs.tracks[0].hitPattern().trackerLayersWithMeasurement() != 4) return;
      }
    }
    if(NLayers == 5){
      if(objs.tracks.size() > 0) {
        if(objs.tracks[0].hitPattern().trackerLayersWithMeasurement() != 5) return;
      }
    }
    if(NLayers == 6){
      if(objs.tracks.size() > 0) {
        if(objs.tracks[0].hitPattern().trackerLayersWithMeasurement() < 6) return;
      }
    }

    oneDHists.at(std::to_string(NLayers) + "nTPOS")->Fill(objs.nTP.nTPOS,sf);
    oneDHists.at(std::to_string(NLayers) + "nTPSS")->Fill(objs.nTP.nTPSS,sf);
    oneDHists.at(std::to_string(NLayers) + "nTPOS_veto")->Fill(objs.nTP.nTPOS_veto,sf);
    oneDHists.at(std::to_string(NLayers) + "nTPSS_veto")->Fill(objs.nTP.nTPSS_veto,sf);
    oneDHists.at(std::to_string(NLayers) + "nTPOSLoose_veto")->Fill(objs.nTP.nTPOS_loose_veto,sf);
    oneDHists.at(std::to_string(NLayers) + "nTPSSLoose_veto")->Fill(objs.nTP.nTPSS_loose_veto,sf);

  }

};

#endif // Analysis_Helper_interface_plotPrintFunctions_h