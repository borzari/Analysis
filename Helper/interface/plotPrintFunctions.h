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

};

#endif // Analysis_Helper_interface_plotPrintFunctions_h