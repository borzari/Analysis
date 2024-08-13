#ifndef Analysis_Helper_interface_sfFunctions_h
#define Analysis_Helper_interface_sfFunctions_h

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

class sfFunctions {

  public:

  static void sfElectrons(double &sf, std::string const &inputFile, std::string const &hist, std::vector<pat::Electron> const &electrons){
    
    TFile *electronInputFile = TFile::Open (inputFile.c_str ());
    TH2F * plot = (TH2F*)electronInputFile->Get(hist.c_str());

    float xMin = plot->GetXaxis()->GetBinCenter(1);
    float xMax = plot->GetXaxis()->GetBinCenter(plot->GetNbinsX());
    float yMin = plot->GetYaxis()->GetBinCenter(1);
    float yMax = plot->GetYaxis()->GetBinCenter(plot->GetNbinsY());

    for (const auto &electron1 : electrons) {
      float eta = electron1.eta();
      // the 2015 ID plots are in |eta| yet the rest are in eta, so check xMin
      if(xMin >= 0) eta = abs(eta);
      // check for eta outside of plot binning
      if(eta > xMax) eta = xMax;
      if(eta < xMin) eta = xMin;

      float pt = electron1.pt();
      if(pt < yMin) pt = yMin;
      if(pt > yMax) pt = yMax;

      float sfValue = plot->GetBinContent(plot->FindBin(eta, pt));
      // float sfError = plot->GetBinError(plot->FindBin(eta, pt));

      sf *= sfValue;
    //   sfUp *= sfValue + sfError;
    //   sfDown *= sfValue - sfError;
    } // end loop over electrons
    delete plot;

  }

  static void sfMuons(double &sf, std::string const &inputFile, std::string const &hist, std::vector<pat::Muon> const &muons){
    
    TFile *muonInputFile = TFile::Open (inputFile.c_str ());
    TH2F * plot = (TH2F*)muonInputFile->Get(hist.c_str());

    float xMin = plot->GetXaxis()->GetBinCenter(1);
    float xMax = plot->GetXaxis()->GetBinCenter(plot->GetNbinsX());
    float yMax = plot->GetYaxis()->GetBinCenter(plot->GetNbinsY());

    for (const auto &muon1 : muons) {
      float pt = muon1.pt();
      if(pt > xMax) pt = xMax;
      if(pt < xMin) pt = xMin;
      float eta = (abs(muon1.eta()) > yMax) ? yMax : abs(muon1.eta());
      int bin = plot->FindBin(pt, eta);

      float sfValue = plot->GetBinContent(bin);
      // float sfError = plot->GetBinError(plot->FindBin(eta, pt));

      sf *= sfValue;
    //   sfUp *= sfValue + sfError;
    //   sfDown *= sfValue - sfError;
    } // end loop over muons
    delete plot;

  }

};

#endif // Analysis_Helper_interface_sfFunctions_h