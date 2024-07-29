#ifndef Analysis_Helper_interface_helperFunctions_h
#define Analysis_Helper_interface_helperFunctions_h

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

struct EtaPhi {
  double eta;
  double phi;
  double sigma;

  EtaPhi(const double a, const double b, const double sigma = -1.0) :
    eta(a),
    phi(b),
    sigma(sigma)
  {
  }
};

struct EtaPhiList : public std::vector<EtaPhi> {
  double minDeltaR;

  EtaPhiList() :
    minDeltaR(0.0)
  {
  }
};

class helperFunctions {

  public:

  static bool passMETHLTPath(const edm::Event &event, const edm::Handle<edm::TriggerResults> &triggerBits)
  {

    Int_t pass_HLT_MET105_IsoTrk50 = 0;
    Int_t pass_HLT_MET120_IsoTrk50 = 0;
    Int_t pass_HLT_PFMET105_IsoTrk50 = 0;
    Int_t pass_HLT_PFMET120_PFMHT120_IDTight = 0;
    Int_t pass_HLT_PFMET130_PFMHT130_IDTight = 0;
    Int_t pass_HLT_PFMET140_PFMHT140_IDTight = 0;
    Int_t pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = 0;
    Int_t pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF = 0;
    Int_t pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF = 0;
    Int_t pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF = 0;
    Int_t pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF = 0;
    Int_t pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = 0;
    Int_t pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight = 0;
    Int_t pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight = 0;
    Int_t pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60 = 0;

    const edm::TriggerNames &allTriggerNamesHLT = event.triggerNames(*triggerBits);

    for(unsigned i = 0; i < allTriggerNamesHLT.size(); i++) {
      std::string thisName = allTriggerNamesHLT.triggerName(i);
      if (thisName.find("HLT_MET105_IsoTrk50_v") == 0) pass_HLT_MET105_IsoTrk50 = triggerBits->accept(i);
      if (thisName.find("HLT_MET120_IsoTrk50_v") == 0) pass_HLT_MET120_IsoTrk50 = triggerBits->accept(i);
      if (thisName.find("HLT_PFMET105_IsoTrk50_v") == 0) pass_HLT_PFMET105_IsoTrk50 = triggerBits->accept(i);
      if (thisName.find("HLT_PFMET120_PFMHT120_IDTight_v") == 0) pass_HLT_PFMET120_PFMHT120_IDTight = triggerBits->accept(i);
      if (thisName.find("HLT_PFMET130_PFMHT130_IDTight_v") == 0) pass_HLT_PFMET130_PFMHT130_IDTight = triggerBits->accept(i);
      if (thisName.find("HLT_PFMET140_PFMHT140_IDTight_v") == 0) pass_HLT_PFMET140_PFMHT140_IDTight = triggerBits->accept(i);
      if (thisName.find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v") == 0) pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = triggerBits->accept(i);
      if (thisName.find("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF_v") == 0) pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF = triggerBits->accept(i);
      if (thisName.find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF_v") == 0) pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF = triggerBits->accept(i);
      if (thisName.find("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF_v") == 0) pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF = triggerBits->accept(i);
      if (thisName.find("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF_v") == 0) pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF = triggerBits->accept(i);
      if (thisName.find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") == 0) pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = triggerBits->accept(i);
      if (thisName.find("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v") == 0) pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight = triggerBits->accept(i);
      if (thisName.find("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v") == 0) pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight = triggerBits->accept(i);
      if (thisName.find("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v") == 0) pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60 = triggerBits->accept(i);
    }

    return (pass_HLT_MET105_IsoTrk50 == 1 || pass_HLT_MET120_IsoTrk50 == 1 || pass_HLT_PFMET105_IsoTrk50 == 1 || pass_HLT_PFMET120_PFMHT120_IDTight == 1 || pass_HLT_PFMET130_PFMHT130_IDTight == 1 || pass_HLT_PFMET140_PFMHT140_IDTight == 1 || pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 == 1 || pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF == 1 || pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF == 1 || pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF == 1 || pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF == 1 || pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight == 1 || pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight == 1 || pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight == 1 || pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60 == 1);

  }

  static bool passLepHLTPath(const edm::Event &event, const edm::Handle<edm::TriggerResults> &triggerBits, const std::string &HLTName)
  {

    Int_t pass_HLT = 0;

    const edm::TriggerNames &allTriggerNamesHLT = event.triggerNames(*triggerBits);

    for(unsigned i = 0; i < allTriggerNamesHLT.size(); i++) {
      std::string thisName = allTriggerNamesHLT.triggerName(i);
      if (thisName.find(HLTName) == 0) pass_HLT = triggerBits->accept(i);
    }

    return (pass_HLT == 1);

  }

  static bool passMETFilters(const edm::Event &event, const edm::Handle<edm::TriggerResults> &triggerBits)
  {

    Int_t pass_Flag_goodVertices = 0;
    Int_t pass_Flag_globalSuperTightHalo2016Filter = 0;
    Int_t pass_Flag_BadPFMuonDzFilter = 0;
    Int_t pass_Flag_hfNoisyHitsFilter = 0;
    Int_t pass_Flag_EcalDeadCellTriggerPrimitiveFilter = 0;
    Int_t pass_Flag_BadPFMuonFilter = 0;
    Int_t pass_Flag_eeBadScFilter = 0;

    const edm::TriggerNames &allTriggerNamesPAT = event.triggerNames(*triggerBits);

    for(unsigned i = 0; i < allTriggerNamesPAT.size(); i++) {
      std::string thisName = allTriggerNamesPAT.triggerName(i);
      if (thisName.find("Flag_goodVertices") == 0) pass_Flag_goodVertices = triggerBits->accept(i);
      if (thisName.find("Flag_globalSuperTightHalo2016Filter") == 0) pass_Flag_globalSuperTightHalo2016Filter = triggerBits->accept(i);
      if (thisName.find("Flag_BadPFMuonDzFilter") == 0) pass_Flag_BadPFMuonDzFilter = triggerBits->accept(i);
      if (thisName.find("Flag_hfNoisyHitsFilter") == 0) pass_Flag_hfNoisyHitsFilter = triggerBits->accept(i);
      if (thisName.find("Flag_EcalDeadCellTriggerPrimitiveFilter") == 0) pass_Flag_EcalDeadCellTriggerPrimitiveFilter = triggerBits->accept(i);
      if (thisName.find("Flag_BadPFMuonFilter") == 0) pass_Flag_BadPFMuonFilter = triggerBits->accept(i);
      if (thisName.find("Flag_eeBadScFilter") == 0) pass_Flag_eeBadScFilter = triggerBits->accept(i);
    }

    return (pass_Flag_goodVertices == 1 && pass_Flag_globalSuperTightHalo2016Filter == 1 && pass_Flag_BadPFMuonDzFilter == 1 && pass_Flag_hfNoisyHitsFilter == 1 && pass_Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && pass_Flag_BadPFMuonFilter == 1 && pass_Flag_eeBadScFilter == 1);

  }

  template<typename T>
  static bool isMatchedToTriggerObject (const edm::Event &event, const edm::TriggerResults &triggers, const T &obj, const std::vector<pat::TriggerObjectStandAlone> &trigObjs, const std::string &collection, const std::string &filter, const double dR = 0.1)
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

  static bool jetPassesTightLepVeto(const pat::Jet &jet)
  {
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV; using CHS jet as it is the type of slimmedJets
    return((jet.neutralHadronEnergyFraction()<0.99 && jet.neutralEmEnergyFraction()<0.90 && (jet.chargedMultiplicity() + jet.neutralMultiplicity())>1 && jet.muonEnergyFraction()<0.8 && jet.chargedHadronEnergyFraction()>0.01 && jet.chargedMultiplicity()>0 && jet.chargedEmEnergyFraction()<0.80 && fabs(jet.eta())<=2.6) || (jet.neutralHadronEnergyFraction()<0.90 && jet.neutralEmEnergyFraction()<0.99 && jet.muonEnergyFraction()<0.8 && jet.chargedMultiplicity()>0 && jet.chargedEmEnergyFraction()<0.80 && fabs(jet.eta())>2.6 && fabs(jet.eta())<=2.7) || (jet.neutralHadronEnergyFraction()<0.99 && jet.neutralEmEnergyFraction()<0.99 && jet.neutralMultiplicity()>=1 && fabs(jet.eta())>2.7 && fabs(jet.eta())<=3.0) || (jet.neutralEmEnergyFraction()<0.4 && jet.neutralMultiplicity()>10 && fabs(jet.eta())>3.0 && fabs(jet.eta())<=5.0));

  }

  static bool IsValidJet(const pat::Jet &jet) {
    if (!(jet.pt() > 30))         return false;
    if (!(fabs(jet.eta()) < 4.5)) return false;
    if (!jetPassesTightLepVeto(jet)) return false;
    return true;
  }

  static bool elecD0 (const pat::Electron& electron, const reco::Vertex& pv){
    return ((fabs (electron.superCluster()->eta()) <= 1.479) && (fabs (((electron.vx() - pv.x()) * electron.py() - (electron.vy() - pv.y()) * electron.px()) / electron.pt()) < 0.05)) || ((fabs (electron.superCluster()->eta()) >  1.479) && (fabs (((electron.vx() - pv.x()) * electron.py() - (electron.vy() - pv.y()) * electron.px()) / electron.pt()) < 0.10));
  }

  static bool elecDZ (const pat::Electron& electron, const reco::Vertex& pv){
    return ((fabs (electron.superCluster()->eta()) <= 1.479) && (fabs ((electron.vz() - pv.z()) - ((electron.vx() - pv.x()) * electron.px() + (electron.vy() - pv.y()) * electron.py()) / electron.pt() * electron.pz() / electron.pt()) < 0.10)) || ((fabs (electron.superCluster()->eta()) >  1.479) && (fabs ((electron.vz() - pv.z()) - ((electron.vx() - pv.x()) * electron.px() + (electron.vy() - pv.y()) * electron.py()) / electron.pt() * electron.pz() / electron.pt()) < 0.20));
  }

  static double muonIso (const pat::Muon& muon){
    return ((muon.pfIsolationR04().sumChargedHadronPt + std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt());
  }

  static bool passesDecayModeReconstruction (const pat::Tau &tau){
    return (tau.tauID("decayModeFindingNewDMs"));
  }

  static bool passesLightFlavorRejection (const pat::Tau &tau){
    return (tau.tauID("byVVVLooseDeepTau2017v2p1VSe") || tau.tauID("byVVVLooseDeepTau2018v2p5VSe")) && (tau.tauID("byVLooseDeepTau2017v2p1VSmu") || tau.tauID("byVLooseDeepTau2018v2p5VSmu"));
  }

  template<typename T>
  static void extractFiducialMap (EtaPhiList &vetoList, std::string histFile)
  {

    std::string era = "";
    std::string beforeVetoHistName = "beforeVeto";
    std::string afterVetoHistName = "afterVeto";
    double thresholdForVeto = 0.0;

    std::cout << "Attempting to extract \"" << beforeVetoHistName << "\" and \"" << afterVetoHistName << "\" from \"" << histFile << "\"..." << std::endl;
    TFile *fin = TFile::Open (histFile.c_str());
    if (!fin || fin->IsZombie ())
      {
        std::cout << "No file named \"" << histFile << "\" found. Skipping..." << std::endl;
        return;
      }

    TH2D * beforeVetoHist = (TH2D*)fin->Get(beforeVetoHistName.c_str());
    TH2D * afterVetoHist  = (TH2D*)fin->Get(afterVetoHistName.c_str());
    if (!beforeVetoHist) {
      std::cout << "No histogram named \"" << beforeVetoHistName.c_str() << "\" found. Skipping..." << std::endl;;
      return;
    }
    if (!afterVetoHist) {
      std::cout << "No histogram named \"" << afterVetoHistName.c_str() << "\" found. Skipping..." << std::endl;;
      return;
    }

    beforeVetoHist->SetDirectory(0);
    afterVetoHist->SetDirectory(0);
    fin->Close();
    delete fin;

    //////////////////////////////////////////////////////////////////////////////
    // First calculate the mean efficiency.
    //////////////////////////////////////////////////////////////////////////////

    double totalEventsBeforeVeto = 0, totalEventsAfterVeto = 0;
    int nBinsWithTags = 0;
    for (int i = 1; i <= beforeVetoHist->GetXaxis ()->GetNbins (); i++)
      {
        for (int j = 1; j <= beforeVetoHist->GetYaxis ()->GetNbins (); j++)
          {
            double binRadius = hypot (0.5 * beforeVetoHist->GetXaxis ()->GetBinWidth (i), 0.5 * beforeVetoHist->GetYaxis ()->GetBinWidth (j));
            (vetoList.minDeltaR < binRadius) && (vetoList.minDeltaR = binRadius);

            double contentBeforeVeto = beforeVetoHist->GetBinContent (i, j);

            if (!contentBeforeVeto) // skip bins that are empty in the before-veto histogram
              continue;

            nBinsWithTags++;

            totalEventsBeforeVeto += contentBeforeVeto;
            totalEventsAfterVeto  += afterVetoHist->GetBinContent (i, j);
          }
      }
    double meanInefficiency = totalEventsAfterVeto / totalEventsBeforeVeto;

    //////////////////////////////////////////////////////////////////////////////
    // Then calculate the standard deviation of the mean inefficiency.
    //////////////////////////////////////////////////////////////////////////////
    double stdDevInefficiency = 0;
    afterVetoHist->Divide (beforeVetoHist);
    for (int i = 1; i <= beforeVetoHist->GetXaxis ()->GetNbins (); i++)
      {
        for (int j = 1; j <= beforeVetoHist->GetYaxis ()->GetNbins (); j++)
          {
            if(beforeVetoHist->GetBinContent (i, j) == 0) // skip bins that are empty in the before-veto histogram
              continue;

            double thisInefficiency = afterVetoHist->GetBinContent (i, j);

            stdDevInefficiency += (thisInefficiency - meanInefficiency) * (thisInefficiency - meanInefficiency);
          }
      }

    if (nBinsWithTags < 2) stdDevInefficiency = 0.;
    else {
      stdDevInefficiency /= nBinsWithTags - 1;
      stdDevInefficiency = sqrt(stdDevInefficiency);
    }

    //////////////////////////////////////////////////////////////////////////////
    // Then find the bins which are greater than the mean by more than
    // thresholdForVeto sigma. Add the coordinates for these bins to the veto
    // list.
    //////////////////////////////////////////////////////////////////////////////
    for (int i = 1; i <= afterVetoHist->GetXaxis ()->GetNbins (); i++)
      {
        for (int j = 1; j <= afterVetoHist->GetYaxis ()->GetNbins (); j++)
          {
            double content = afterVetoHist->GetBinContent (i, j),
              //error = afterVetoHist->GetBinError (i, j),
                   eta = afterVetoHist->GetXaxis ()->GetBinCenter (i),
                   phi = afterVetoHist->GetYaxis ()->GetBinCenter (j);

            if ((content - meanInefficiency) > thresholdForVeto * stdDevInefficiency)
              {
                vetoList.emplace_back (eta, phi, (content - meanInefficiency) / stdDevInefficiency);
              }
          }
      }
    //////////////////////////////////////////////////////////////////////////////
    delete beforeVetoHist;
    delete afterVetoHist;
  }

  static bool isFiducialTrack (const pat::IsolatedTrack &track, const EtaPhiList &vetoList, const double minDeltaR, double maxSigma)
  {
    const double minDR = std::max(minDeltaR, vetoList.minDeltaR); // use the given parameter unless the bin size from which the veto list is calculated is larger
    bool isFiducial = true;
    maxSigma = 0.0;
    for (const auto &etaPhi : vetoList)
      {
        if (deltaR (track.eta (), track.phi (), etaPhi.eta, etaPhi.phi) < minDR)
          {
            isFiducial = false;
            if (etaPhi.sigma > maxSigma)
              maxSigma = etaPhi.sigma;
          }
      }
    return isFiducial;
  }

  static int getChannelStatusMaps (edm::ESHandle<EcalChannelStatus> &ecalStatus, const edm::ESHandle<CaloGeometry>& caloGeometry, std::map<DetId, std::vector<double> > &EcalAllDeadChannelsValMap, std::map<DetId, std::vector<int> > &EcalAllDeadChannelsBitMap)
  {
    bool makeHist = false;
    EcalAllDeadChannelsValMap.clear(); EcalAllDeadChannelsBitMap.clear();
    TH2D *badChannels = (makeHist ? new TH2D ("badChannels", ";#eta;#phi", 360, -3.0, 3.0, 360, -3.2, 3.2) : NULL);

  // Loop over EB ...
    for( int ieta=-85; ieta<=85; ieta++ ){
       for( int iphi=0; iphi<=360; iphi++ ){
          if(! EBDetId::validDetId( ieta, iphi ) )  continue;

          const EBDetId detid = EBDetId( ieta, iphi, EBDetId::ETAPHIMODE );
          EcalChannelStatus::const_iterator chit = ecalStatus->find( detid );
  // refer https://twiki.cern.ch/twiki/bin/viewauth/CMS/EcalChannelStatus
          int status = ( chit != ecalStatus->end() ) ? chit->getStatusCode() & 0x1F : -1;

          const CaloSubdetectorGeometry*  subGeom = caloGeometry->getSubdetectorGeometry (detid);
          auto cellGeom = subGeom->getGeometry (detid);
          double eta = cellGeom->getPosition ().eta ();
          double phi = cellGeom->getPosition ().phi ();
          double theta = cellGeom->getPosition().theta();

          if(status >= 3){
             std::vector<double> valVec; std::vector<int> bitVec;
             valVec.push_back(eta); valVec.push_back(phi); valVec.push_back(theta);
             bitVec.push_back(1); bitVec.push_back(ieta); bitVec.push_back(iphi); bitVec.push_back(status);
             EcalAllDeadChannelsValMap.insert( std::make_pair(detid, valVec) );
             EcalAllDeadChannelsBitMap.insert( std::make_pair(detid, bitVec) );
             if (makeHist)
               badChannels->Fill (eta, phi);
          }
       } // end loop iphi
    } // end loop ieta

  // Loop over EE detid
    for( int ix=0; ix<=100; ix++ ){
       for( int iy=0; iy<=100; iy++ ){
          for( int iz=-1; iz<=1; iz++ ){
             if(iz==0)  continue;
             if(! EEDetId::validDetId( ix, iy, iz ) )  continue;

             const EEDetId detid = EEDetId( ix, iy, iz, EEDetId::XYMODE );
             EcalChannelStatus::const_iterator chit = ecalStatus->find( detid );
             int status = ( chit != ecalStatus->end() ) ? chit->getStatusCode() & 0x1F : -1;

             const CaloSubdetectorGeometry*  subGeom = caloGeometry->getSubdetectorGeometry (detid);
             auto cellGeom = subGeom->getGeometry (detid);
             double eta = cellGeom->getPosition ().eta () ;
             double phi = cellGeom->getPosition ().phi () ;
             double theta = cellGeom->getPosition().theta();

             if(status >= 3){
                std::vector<double> valVec; std::vector<int> bitVec;
                valVec.push_back(eta); valVec.push_back(phi); valVec.push_back(theta);
                bitVec.push_back(2); bitVec.push_back(ix); bitVec.push_back(iy); bitVec.push_back(iz); bitVec.push_back(status);
                EcalAllDeadChannelsValMap.insert( std::make_pair(detid, valVec) );
                EcalAllDeadChannelsBitMap.insert( std::make_pair(detid, bitVec) );
                 if (makeHist)
                   badChannels->Fill (eta, phi);
             }
          } // end loop iz
       } // end loop iy
    } // end loop ix

    if (makeHist)
      {
        TFile *fout = new TFile ("badEcalChannels.root", "recreate");
        fout->cd ();
        badChannels->Write ();
        fout->Close ();

        delete badChannels;
        delete fout;
      }

    return 1;
  }

  static int isCloseToBadEcalChannel (const pat::IsolatedTrack &track, const double deltaRCut, std::map<DetId, std::vector<double> > &EcalAllDeadChannelsValMap, std::map<DetId, std::vector<int> > &EcalAllDeadChannelsBitMap)
  {
     double trackEta = track.eta(), trackPhi = track.phi();

     double min_dist = 999;
     DetId min_detId;

     std::map<DetId, std::vector<int> >::const_iterator bitItor;
     for(bitItor = EcalAllDeadChannelsBitMap.begin(); bitItor != EcalAllDeadChannelsBitMap.end(); bitItor++){

        DetId maskedDetId = bitItor->first;

        std::map<DetId, std::vector<double> >::const_iterator valItor = EcalAllDeadChannelsValMap.find(maskedDetId);
        if( valItor == EcalAllDeadChannelsValMap.end() ){ std::cout<<"Error cannot find maskedDetId in EcalAllDeadChannelsValMap ?!"<<std::endl; continue; }

        double eta = (valItor->second)[0], phi = (valItor->second)[1];

        double dist = reco::deltaR(eta, phi, trackEta, trackPhi);

        if( min_dist > dist ){ min_dist = dist; min_detId = maskedDetId; }
     }

     if( min_dist > deltaRCut && deltaRCut > 0 ) return 0;

     return 1;
  }

  static bool inTOBCrack (const pat::IsolatedTrack &track){
    return (fabs(track.dz()) < 0.5 && fabs(1.57079632679489661923 - track.theta()) < 1.0e-3);
  }

  static int extraMissingMiddleHits (const pat::IsolatedTrack &track)
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
  
  static int hitDrop_missingMiddleHits (const pat::IsolatedTrack &track)
  {
    int nDropHits = extraMissingMiddleHits(track);
    return track.hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS) + nDropHits;
  }
  
  static double dRMinJet (const pat::IsolatedTrack &track, const std::vector<pat::Jet> &jets)
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
  
  template<typename T>
  static double deltaRToClosestLepton(const pat::IsolatedTrack &track, const std::vector<T> &leptons) 
  {
  
    double dR;
    double deltaRToClosestLepton = 99.0;
  
    for(const auto &lepton : leptons) {
      dR = deltaR(track, lepton);
      if(dR < deltaRToClosestLepton || deltaRToClosestLepton < 0.0) deltaRToClosestLepton = dR;
    }
  
    return deltaRToClosestLepton;
  
  }
  
  static double deltaRToClosestTauHad(const pat::IsolatedTrack &track, const std::vector<pat::Tau> &taus) 
  {
  
    double dR;
    double deltaRToClosestTauHad = 99.0;
  
    for(const auto &tau : taus) {
      dR = deltaR(track, tau);
  
      if(passesDecayModeReconstruction(tau) && passesLightFlavorRejection(tau) && (dR < deltaRToClosestTauHad || deltaRToClosestTauHad < 0.0)) {
        deltaRToClosestTauHad = dR;
      }
    }
  
    return deltaRToClosestTauHad;
  }
  
  static double energyGivenMass (const double mass, const pat::IsolatedTrack &track)
  {
    return sqrt (track.px () * track.px () + track.py () * track.py () + track.pz () * track.pz () + mass * mass);
  }

  template<typename T>
  static bool goodInvMassLepton (const T &tag, const pat::IsolatedTrack &probe, bool isTau)
  {
    double lepMass = 0.0;

    if(std::is_same<T, pat::Electron>::value) lepMass = 0.000510998950; // Electron mass extracted from PDG on 27/06/2024 https://pdg.lbl.gov/2024/tables/contents_tables.html
    if(std::is_same<T, pat::Muon>::value) lepMass = 0.1056583755; // Muon mass extracted from PDG on 27/06/2024 https://pdg.lbl.gov/2024/tables/contents_tables.html
    if(isTau) lepMass = 0.13957039; // Pion mass extracted from PDG on 16/07/2024 https://pdg.lbl.gov/2024/tables/contents_tables.html

    TLorentzVector t (tag.px(), tag.py(), tag.pz(), tag.energy()),
                   p (probe.px(), probe.py(), probe.pz(), energyGivenMass(lepMass, probe));
    double m = (t + p).M();
    if(!isTau) return (fabs (m - 91.1880) < 10.0); // Z mass extracted from PDG on 27/06/2024 https://pdg.lbl.gov/2024/tables/contents_tables.html
    return (15.0 < (91.1880 - m) && (91.1880 - m) < 50.0); // Z mass extracted from PDG on 27/06/2024 https://pdg.lbl.gov/2024/tables/contents_tables.html
  }
  
  template<typename T>
  static double deltaRToClosestPFLepton(const pat::IsolatedTrack &track, const std::vector<pat::PackedCandidate> &pfCandidates) 
  {

    int lepPdgid = 0;

    if(std::is_same<T, pat::Electron>::value) lepPdgid = 11;
    if(std::is_same<T, pat::Muon>::value) lepPdgid = 13;
    if(std::is_same<T, pat::Tau>::value) lepPdgid = 211;

    double deltaRToClosestPFLepton = 99.0;
    for(const auto &pfCandidate : pfCandidates) {
        int pdgid = abs(pfCandidate.pdgId());
        if(pdgid != lepPdgid) continue;
  
        double dR = deltaR(track, pfCandidate);
  
        if(pdgid == lepPdgid &&
           (dR < deltaRToClosestPFLepton || deltaRToClosestPFLepton < 0.0))
          deltaRToClosestPFLepton = dR;
    }
    return deltaRToClosestPFLepton;
  }
  
  static double deltaRToClosestVetoElectron (const pat::IsolatedTrack &track, const std::vector<pat::Electron> &electrons, const reco::Vertex &vertex)
  {
    double deltaRToClosestVetoElectron = 99.0;
    
    double dR;
  
    for(const auto &electron : electrons) {
      dR = deltaR(track, electron);
  
      bool passesVeto_dxy = false, passesVeto_dz = false;
  
      // Note in below, these remain false if |eta| >= 2.5; thus an eta cut is also being applied here as intended
      double ele_d0 = fabs(electron.gsfTrack()->dxy(vertex.position()));
      double ele_dz = fabs(electron.gsfTrack()->dz(vertex.position()));
  
      if(fabs(electron.superCluster ()->eta()) <= 1.479) {
        passesVeto_dxy = (ele_d0 < 0.05);
        passesVeto_dz = (ele_dz < 0.10);
      }
      else if(fabs(electron.superCluster()->eta()) < 2.5) {
        passesVeto_dxy = (ele_d0 < 0.10);
        passesVeto_dz = (ele_dz < 0.20);
      }
  
      if(electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-tight") &&
         passesVeto_dxy &&
         passesVeto_dz &&
         (dR < deltaRToClosestVetoElectron || deltaRToClosestVetoElectron < 0.0)) {
        deltaRToClosestVetoElectron = dR;
      }
    } // for electrons
  
    return deltaRToClosestVetoElectron;
  }
  
  static double deltaRToClosestLooseMuon(const pat::IsolatedTrack &track, const std::vector<pat::Muon> &muons) 
  {
    double deltaRToClosestLooseMuon = 99.0;
  
    double dR;
  
    for(const auto &muon : muons) {
      dR = deltaR(track, muon);
      if(muon.isLooseMuon()  && (dR < deltaRToClosestLooseMuon  || deltaRToClosestLooseMuon  < 0.0)) deltaRToClosestLooseMuon = dR;
    }
  
    return deltaRToClosestLooseMuon;
  }
  
  template<typename T>
  static bool passesVeto (const pat::IsolatedTrack &probe, const std::vector<pat::PackedCandidate> &pfCandidates, const std::vector<pat::Jet> &jets)
  {
    bool passesElec = deltaRToClosestPFLepton<pat::Electron>(probe, pfCandidates) > 0.15
               && (probe.matchedCaloJetEmEnergy() + probe.matchedCaloJetHadEnergy()) < 10.0;
               // && probe.hitAndTOBDrop_bestTrackMissingOuterHits () >= 3.0; // This is not applied for BG MC
              //  && probe.lostOuterLayers() >= 3.0; // This is not applied for BG MC
    bool passesMuon = deltaRToClosestPFLepton<pat::Muon>(probe, pfCandidates) > 0.15;
               // && probe.hitAndTOBDrop_bestTrackMissingOuterHits () >= 3.0; // This is not applied for BG MC
              //  && probe.lostOuterLayers() >= 3.0; // This is not applied for BG MC
    bool passesTau = deltaRToClosestPFLepton<pat::Tau>(probe, pfCandidates) > 0.15
               && dRMinJet (probe, jets) > 0.5
               && (probe.matchedCaloJetEmEnergy() + probe.matchedCaloJetHadEnergy()) < 10.0;
               // && probe.hitAndTOBDrop_bestTrackMissingOuterHits () >= 3.0; // This is not applied for BG MC
              //  && probe.lostOuterLayers() >= 3.0; // This is not applied for BG MC
  
    if(std::is_same<T, pat::Electron>::value) return passesElec;
    if(std::is_same<T, pat::Muon>::value) return passesMuon;
    if(std::is_same<T, pat::Tau>::value) return passesTau;
  }

  static bool passesLooseElecVeto (const pat::IsolatedTrack &probe, const std::vector<pat::Electron> &electrons, const reco::Vertex &vertex)
  {
    bool passes = deltaRToClosestVetoElectron(probe,electrons,vertex) > 0.15
               && (probe.matchedCaloJetEmEnergy() + probe.matchedCaloJetHadEnergy()) < 10.0;
               // && probe.hitAndTOBDrop_bestTrackMissingOuterHits () >= 3.0; // This is not applied for BG MC
              //  && probe.lostOuterLayers() >= 3.0; // This is not applied for BG MC
    return passes;
  }
   
  static bool passesLooseMuonVeto (const pat::IsolatedTrack &probe, const std::vector<pat::Muon> &muons)
  {
    bool passes = deltaRToClosestLooseMuon(probe, muons) > 0.15;
               // && probe.hitAndTOBDrop_bestTrackMissingOuterHits () >= 3.0; // This is not applied for BG MC
              //  && probe.lostOuterLayers() >= 3.0; // This is not applied for BG MC
    return passes;
  }
  
  static GlobalPoint getPosition(const DetId& id, const edm::ESHandle<CaloGeometry>& caloGeometry)
  {
     if ( ! caloGeometry.isValid() ||
          ! caloGeometry->getSubdetectorGeometry(id) ||
          ! caloGeometry->getSubdetectorGeometry(id)->getGeometry(id) ) {
        throw cms::Exception("FatalError") << "Failed to access geometry for DetId: " << id.rawId();
        return GlobalPoint(0,0,0);
     }
     return caloGeometry->getSubdetectorGeometry(id)->getGeometry(id)->getPosition();
  }
   
  static bool insideCone(const pat::IsolatedTrack &candTrack, const DetId& id, const edm::ESHandle<CaloGeometry>& caloGeometry)
  {
     GlobalPoint idPosition = getPosition(id, caloGeometry);
     if (idPosition.mag()<0.01) return false;
     math::XYZVector idPositionRoot( idPosition.x(), idPosition.y(), idPosition.z() );
     return deltaR(candTrack, idPositionRoot) < 0.5;
  }
  
  static double calculateCaloE(const pat::IsolatedTrack& track, const EBRecHitCollection &EBRecHits, const EERecHitCollection &EERecHits, const HBHERecHitCollection &HBHERecHits, const edm::ESHandle<CaloGeometry>& caloGeometry)
  { 
  
    double caloEnergy = 0.0;
  
    for (const auto &hit : EBRecHits) {
      if (insideCone(track, hit.detid(), caloGeometry)) {
        caloEnergy += hit.energy();
      }
    }
    for (const auto &hit : EERecHits) {
      if (insideCone(track, hit.detid(), caloGeometry)) {
        caloEnergy += hit.energy();
      }
    }
  
    for (const auto &hit : HBHERecHits) {
      if (insideCone(track, hit.detid(), caloGeometry)) {
        caloEnergy += hit.energy();
      }
    }
  
    return caloEnergy;
  }
  
  static double caloNewNoPUDRp5CentralCalo(const pat::IsolatedTrack& track, const EBRecHitCollection &EBRecHits, const EERecHitCollection &EERecHits, const HBHERecHitCollection &HBHERecHits, const double rhoCentralCalo, const edm::ESHandle<CaloGeometry>& caloGeometry)
  {
    
    double rawCaloTot = calculateCaloE(track, EBRecHits, EERecHits, HBHERecHits, caloGeometry);
    double caloCorr = rhoCentralCalo * TMath::Pi() * 0.5 * 0.5;  // Define effective area as pi*r^2, where r is radius of DeltaR cone.
    double caloNewNoPUDRp5CentralCalo = TMath::Max(0., rawCaloTot - caloCorr);
  
    return caloNewNoPUDRp5CentralCalo;
  }
  
  template<typename T>
  static double transvMassLepton(const T& lepton, const pat::MET& met)
  {

    double dPhi = deltaPhi (lepton.phi(), met.phi());
	  return sqrt(2.0 * lepton.pt() * met.pt() * (1.0 - cos(dPhi)));
    
  }

  template<class T>
  static const pat::TriggerObjectStandAlone * getMatchedTriggerObject (const edm::Event &event, const edm::TriggerResults &triggers, const T &obj, const std::vector<pat::TriggerObjectStandAlone> &trigObjs, const std::string &collection, const std::string &filter, const double dR = 0.1)
  {
    if (collection == "")
      return NULL;
    double minDR = -1.0;
    int i = -1, iMinDR = -1;
    for (auto trigObj : trigObjs)
      {
        i++;
        trigObj.unpackNamesAndLabels(event, triggers);
        if (trigObj.collection () != collection)
          continue;
        if (filter != "")
          {
            bool flag = false;
            for (const auto &filterLabel : trigObj.filterLabels ())
              if (filterLabel == filter)
                {
                  flag = true;
                  break;
                }
            if (!flag)
              continue;
          }
        double currentDR = deltaR (obj, trigObj);
        if (currentDR > dR)
          continue;

        if (minDR < 0 || currentDR < minDR)
          {
            minDR = currentDR;
            iMinDR = i;
          }
      }
    if (iMinDR >= 0)
      return &trigObjs.at (iMinDR);
    return NULL;
  }

  static bool getTriggerObjects (const edm::Event &event, const edm::TriggerResults &triggers, const std::vector<pat::TriggerObjectStandAlone> &trigObjs, const std::string &collection, const std::string &filter, std::vector<const pat::TriggerObjectStandAlone *> &selectedTrigObjs)
  {
    if (collection == "")
      {
        selectedTrigObjs.push_back (NULL);
        return false;
      }
    std::vector<const pat::TriggerObjectStandAlone *> trigObjsToAdd;
    int i = -1;
    for (auto trigObj : trigObjs)
      {
        i++;
        trigObj.unpackNamesAndLabels(event, triggers);
        if (trigObj.collection () != collection)
          continue;
        if (filter != "")
          {
            bool flag = false;
            for (const auto &filterLabel : trigObj.filterLabels ())
              if (filterLabel == filter)
                {
                  flag = true;
                  break;
                }
            if (!flag)
              continue;
          }

        trigObjsToAdd.push_back (&trigObjs.at (i));
      }
    if (!trigObjsToAdd.empty ())
      selectedTrigObjs.insert (selectedTrigObjs.end (), trigObjsToAdd.begin (), trigObjsToAdd.end ());
    else
      selectedTrigObjs.push_back (NULL);
    return (!trigObjsToAdd.empty ());
  }

  static const double getModifiedMissingEnergy (const TVector2 &x, const TVector2 &y, const bool muonsCountedAsVisible, const double shift)
  {
    TVector2 z;
    z.SetMagPhi (shift, y.Phi ());
    if (muonsCountedAsVisible)
      return (x + y - z).Mod ();
    return x.Mod ();
  }

  static bool triggerObjectExists (const edm::Event &event, const edm::TriggerResults &triggers, const std::vector<pat::TriggerObjectStandAlone> &trigObjs, const std::string &collection, const std::string &filter)
  {
    if (collection == "")
      return false;
    for (auto trigObj : trigObjs)
      {
        trigObj.unpackNamesAndLabels(event, triggers);
        if (trigObj.collection () != collection)
          continue;
        if (filter != "")
          {
            bool flag = false;
            for (const auto &filterLabel : trigObj.filterLabels ())
              if (filterLabel == filter)
                {
                  flag = true;
                  break;
                }
            if (!flag)
              continue;
          }

        return true;
      }
    return false;
  }

  template<class T>
  static std::vector<bool> passMETTriggers (const edm::Event &event, edm::Handle<edm::TriggerResults> &triggers, edm::Handle<std::vector<pat::TriggerObjectStandAlone> > &triggerObjects, edm::Handle<std::vector<T> > &tags, std::string tagCollection)
  {

    std::vector<bool> passesVec;

    if (!triggers.isValid () || !triggerObjects.isValid () || !tags.isValid ()){
      passesVec.push_back(false);
      passesVec.push_back(false);
      return passesVec;
    }

    bool passes = false;
    bool passesUp = false;

    std::vector<std::string> filterCategories = {"met", "metClean", "metCleanBH", "metCleanUsingJetID", "mht", "pfMHTTightID", "pfMET", "pfMHT", "pfMETNoMu", "pfMHTNoMu", "pfMHTNoMuCleaned", "pfHT"};

    std::map<std::string, std::vector<std::string> > trigObjCollections;
    std::map<std::string, std::vector<double> > trigObjThresholds;
    std::map<std::string, std::vector<std::string> > trigObjJetsForTag;
    std::map<std::string, bool> muonsCountedAsVisible;

    std::vector<std::string> additionalCollections = {"hltTrk50Filter::HLT","hltTrk50Filter::HLT","hltTrk50Filter::HLT","","","","","","","","","","","",""};
    std::vector<std::string> additionalFilters = {"hltTrk50Filter","hltTrk50Filter","hltTrk50Filter::HLT","","","","","","","","","","","",""};

    trigObjCollections["met"] = {"hltMet::HLT","hltMet::HLT","hltMet::HLT","hltMet::HLT","hltMet::HLT","hltMet::HLT","hltMet::HLT","hltMet::HLT","hltMet::HLT","hltMet::HLT","hltMet::HLT","hltMet::HLT","hltMet::HLT","hltMet::HLT","hltMet::HLT"};
    trigObjCollections["metClean"] = {"","","","","","","","","","","","","","",""};
    trigObjCollections["metCleanBH"] = {"","","","","","","","","","","","","","",""};
    trigObjCollections["metCleanUsingJetID"] = {"","","","","","","","","","","","","","",""};
    trigObjCollections["mht"] = {"","","","hltMht::HLT","hltMht::HLT","hltMht::HLT","hltMht::HLT","hltMht::HLT","hltMht::HLT","hltMht::HLT","hltMht::HLT","hltMht::HLT","hltMht::HLT","hltMht::HLT","hltMht::HLT"};
    trigObjCollections["pfMHTTightID"] = {"","","","hltPFMHTTightID::HLT","hltPFMHTTightID::HLT","hltPFMHTTightID::HLT","","","","","","","","","hltPFMHTTightID::HLT"};
    trigObjCollections["pfMET"] = {"","","hltPFMETProducer::HLT","hltPFMETProducer::HLT","hltPFMETProducer::HLT","hltPFMETProducer::HLT","","","","","","","","","hltPFMETProducer::HLT"};
    trigObjCollections["pfMHT"] = {"","","","","","","","","","","","","","",""};
    trigObjCollections["pfMETNoMu"] = {"","","","","","","hltPFMETNoMuProducer::HLT","hltPFMETNoMuProducer::HLT","hltPFMETNoMuProducer::HLT","hltPFMETNoMuProducer::HLT","hltPFMETNoMuProducer::HLT","hltPFMETNoMuProducer::HLT","hltPFMETNoMuProducer::HLT","hltPFMETNoMuProducer::HLT",""};
    trigObjCollections["pfMHTNoMu"] = {"","","","","","","hltPFMHTNoMuTightID::HLT","","","","","hltPFMHTNoMuTightID::HLT","hltPFMHTNoMuTightID::HLT","hltPFMHTNoMuTightID::HLT",""};
    trigObjCollections["pfMHTNoMuCleaned"] = {"","","","","","","","hltPFMHTNoMuTightIDHFCleaned::HLT","hltPFMHTNoMuTightIDHFCleaned::HLT","hltPFMHTNoMuTightIDHFCleaned::HLT","hltPFMHTNoMuTightIDHFCleaned::HLT","","","",""};
    trigObjCollections["pfHT"] = {"","","","","","","hltPFHTJet30::HLT","","","","","","","","hltPFHTJet30::HLT"};

    trigObjThresholds["met"] = {105.0, 120.0, 75.0, 90.0, 100.0, 110.0, 90.0, 80.0, 90.0, 100.0, 110.0, 90.0, 100.0, 110.0, 90.0};
    trigObjThresholds["metClean"] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    trigObjThresholds["metCleanBH"] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    trigObjThresholds["metCleanUsingJetID"] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    trigObjThresholds["mht"] = {0.0, 0.0, 0.0, 90.0, 100.0, 110.0, 90.0, 80.0, 90.0, 100.0, 110.0, 90.0, 100.0, 110.0, 90.0};
    trigObjThresholds["pfMHTTightID"] = {0.0, 0.0, 0.0, 120.0, 130.0, 140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 120.0};
    trigObjThresholds["pfMET"] = {0.0, 0.0, 105.0, 120.0, 130.0, 140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 120.0};
    trigObjThresholds["pfMHT"] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    trigObjThresholds["pfMETNoMu"] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 120.0, 110.0, 120.0, 130.0, 140.0, 120.0, 130.0, 140.0, 0.0};
    trigObjThresholds["pfMHTNoMu"] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 120.0, 0.0, 0.0, 0.0, 0.0, 120.0, 130.0, 140.0, 0.0};
    trigObjThresholds["pfMHTNoMuCleaned"] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 110.0, 120.0, 130.0, 140.0, 0.0, 0.0, 0.0, 0.0};
    trigObjThresholds["pfHT"] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 60.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 60.0};

    trigObjJetsForTag["met"] = {"","","","","","","","","","","","","","",""};
    trigObjJetsForTag["metClean"] = {"","","","","","","","","","","","","","",""};
    trigObjJetsForTag["metCleanBH"] = {"","","","","","","","","","","","","","",""};
    trigObjJetsForTag["metCleanUsingJetID"] = {"","","","","","","","","","","","","","",""};
    trigObjJetsForTag["mht"] = {"","","","hltAK4CaloJetsCorrectedIDPassed::HLT","hltAK4CaloJetsCorrectedIDPassed::HLT","hltAK4CaloJetsCorrectedIDPassed::HLT","hltAK4CaloJetsCorrectedIDPassed::HLT","hltAK4CaloJetsCorrectedIDPassed::HLT","hltAK4CaloJetsCorrectedIDPassed::HLT","hltAK4CaloJetsCorrectedIDPassed::HLT","hltAK4CaloJetsCorrectedIDPassed::HLT","hltAK4CaloJetsCorrectedIDPassed::HLT","hltAK4CaloJetsCorrectedIDPassed::HLT","hltAK4CaloJetsCorrectedIDPassed::HLT","hltAK4CaloJetsCorrectedIDPassed::HLT"};
    trigObjJetsForTag["pfMHTTightID"] = {"","","","hltAK4PFJetsTightIDCorrected::HLT","hltAK4PFJetsTightIDCorrected::HLT","hltAK4PFJetsTightIDCorrected::HLT","","","","","","","","","hltAK4PFJetsTightIDCorrected::HLT"};
    trigObjJetsForTag["pfMET"] = {"","","","","","","","","","","","","","",""};
    trigObjJetsForTag["pfMHT"] = {"","","","","","","","","","","","","","",""};
    trigObjJetsForTag["pfMETNoMu"] = {"","","","","","","","","","","","","","",""};
    trigObjJetsForTag["pfMHTNoMu"] = {"","","","","","","hltAK4PFJetsTightIDCorrected::HLT","","","","","hltAK4PFJetsTightIDCorrected::HLT","hltAK4PFJetsTightIDCorrected::HLT","hltAK4PFJetsTightIDCorrected::HLT",""};
    trigObjJetsForTag["pfMHTNoMuCleaned"] = {"","","","","","","","hltAK4PFJetsTightIDCorrectedHFCleaned::HLT","hltAK4PFJetsTightIDCorrectedHFCleaned::HLT","hltAK4PFJetsTightIDCorrectedHFCleaned::HLT","hltAK4PFJetsTightIDCorrectedHFCleaned::HLT","","","",""};
    trigObjJetsForTag["pfHT"] = {"","","","","","","hltAK4PFJetsCorrected::HLT","","","","","","","","hltAK4PFJetsCorrected::HLT"};

    muonsCountedAsVisible["met"] = false;
    muonsCountedAsVisible["metClean"] = false;
    muonsCountedAsVisible["metCleanBH"] = false;
    muonsCountedAsVisible["metCleanUsingJetID"] = false;
    muonsCountedAsVisible["mht"] = false;
    muonsCountedAsVisible["pfMHTTightID"] = true;
    muonsCountedAsVisible["pfMET"] = true;
    muonsCountedAsVisible["pfMHT"] = true;
    muonsCountedAsVisible["pfMETNoMu"] = false;
    muonsCountedAsVisible["pfMHTNoMu"] = false;
    muonsCountedAsVisible["pfMHTNoMuCleaned"] = false;
    muonsCountedAsVisible["pfHT"] = false;

    const pat::TriggerObjectStandAlone *hltTag = getMatchedTriggerObject<T> (event, *triggers, tags->at (0), *triggerObjects, tagCollection, "");
    std::map<std::string, std::vector<const pat::TriggerObjectStandAlone *> > hltFilterObjects, hltJetsForTag;
    int n = -1;
    for (const auto &filterCategory : filterCategories)
      {
        const std::vector<std::string> &collections = trigObjCollections.at (filterCategory);
        const std::vector<std::string> &jetsForTag = trigObjJetsForTag.at (filterCategory);

        if (n < 0)
          n = collections.size ();

        std::vector<const pat::TriggerObjectStandAlone *> &objs = hltFilterObjects[filterCategory];
        std::vector<const pat::TriggerObjectStandAlone *> &jets = hltJetsForTag[filterCategory];

        for (int i = 0; i < n; i++)
          {
            const std::string &collection = collections.at (i);
            const std::string &jetForTag = jetsForTag.at (i);

            getTriggerObjects (event, *triggers, *triggerObjects, collection, "", objs);
            jets.push_back (getMatchedTriggerObject<T> (event, *triggers, tags->at (0), *triggerObjects, jetForTag, ""));
          }
      }

    std::map<std::string, std::vector<bool> > filterDecisions, filterDecisionsUp;
    for (const auto &filterCategory : filterCategories)
      {
        std::vector<bool> &filterDecision = filterDecisions[filterCategory];
        std::vector<bool> &filterDecisionUp = filterDecisionsUp[filterCategory];

        const std::vector<const pat::TriggerObjectStandAlone *> &objs = hltFilterObjects.at (filterCategory);
        const std::vector<const pat::TriggerObjectStandAlone *> &jets = hltJetsForTag.at (filterCategory);
        const std::vector<double> &thresholds = trigObjThresholds.at (filterCategory);

        for (int i = 0; i < n; i++)
          {
            bool flag, flagUp;

            const pat::TriggerObjectStandAlone *tag = (trigObjJetsForTag.at (filterCategory).at (i) == "" ? hltTag : jets.at (i));
            const pat::TriggerObjectStandAlone *obj = objs.at (i);

            TVector2 x, y;
            if (!tag)
              y.Set (0.0, 0.0);
            else
              y.Set (tag->px (), tag->py ());
            if (!obj && trigObjCollections.at (filterCategory).at (i) == "")
              flag = flagUp = true;
            else
              {
                if (!obj)
                  x.Set (0.0, 0.0);
                else
                  x.Set (obj->px (), obj->py ());

                const double threshold = thresholds.at (i);
                const double modifiedMissingEnergy = getModifiedMissingEnergy (x, y, muonsCountedAsVisible.at (filterCategory), 0.0);
                const double modifiedMissingEnergyUp = getModifiedMissingEnergy (x, y, muonsCountedAsVisible.at (filterCategory), 10.0);
                flag = (modifiedMissingEnergy > threshold);
                flagUp = (modifiedMissingEnergyUp > threshold);
              }

            filterDecision.push_back (flag);
            filterDecisionUp.push_back (flagUp);
          }
      }

    for (int i = 0; i < n; i++)
      {
        bool triggerPasses = true, triggerPassesUp = true;
        for (const auto &filterCategory : filterCategories)
          {
            triggerPasses = triggerPasses && filterDecisions.at (filterCategory).at (i);
            triggerPassesUp = triggerPassesUp && filterDecisionsUp.at (filterCategory).at (i);
          }
        if (additionalCollections.at (i) != "" && additionalFilters.at (i) != "")
          {
            bool filterPasses = triggerObjectExists (event, *triggers, *triggerObjects, additionalCollections.at (i), additionalFilters.at (i));
            triggerPasses = triggerPasses && filterPasses;
            triggerPassesUp = triggerPassesUp && filterPasses;
          }
        passes = passes || triggerPasses;
        passesUp = passesUp || triggerPassesUp;
      }

    passesVec.push_back(passes);
    passesVec.push_back(passesUp);

    return passesVec;

  }

};
  
#endif  // Analysis_Helper_interface_helperFunctions_h