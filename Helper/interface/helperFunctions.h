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

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"

class helperFunctions {

  public:

  static bool passHLTPath(const edm::Event &event, const edm::Handle<edm::TriggerResults> &triggerBits, const std::string &HLTName)
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
               && (probe.matchedCaloJetEmEnergy() + probe.matchedCaloJetHadEnergy()) < 10.0//;
               // && probe.hitAndTOBDrop_bestTrackMissingOuterHits () >= 3.0; // This is not applied for BG MC
               && probe.lostOuterLayers() >= 3.0; // This is not applied for BG MC
    bool passesMuon = deltaRToClosestPFLepton<pat::Muon>(probe, pfCandidates) > 0.15
               // && probe.hitAndTOBDrop_bestTrackMissingOuterHits () >= 3.0; // This is not applied for BG MC
               && probe.lostOuterLayers() >= 3.0; // This is not applied for BG MC
    bool passesTau = deltaRToClosestPFLepton<pat::Tau>(probe, pfCandidates) > 0.15
               && dRMinJet (probe, jets) > 0.5
               && (probe.matchedCaloJetEmEnergy() + probe.matchedCaloJetHadEnergy()) < 10.0
               // && probe.hitAndTOBDrop_bestTrackMissingOuterHits () >= 3.0; // This is not applied for BG MC
               && probe.lostOuterLayers() >= 3.0; // This is not applied for BG MC
  
    if(std::is_same<T, pat::Electron>::value) return passesElec;
    if(std::is_same<T, pat::Muon>::value) return passesMuon;
    if(std::is_same<T, pat::Tau>::value) return passesTau;
  }

  static bool passesLooseElecVeto (const pat::IsolatedTrack &probe, const std::vector<pat::Electron> &electrons, const reco::Vertex &vertex)
  {
    bool passes = deltaRToClosestVetoElectron(probe,electrons,vertex) > 0.15
               && (probe.matchedCaloJetEmEnergy() + probe.matchedCaloJetHadEnergy()) < 10.0
               // && probe.hitAndTOBDrop_bestTrackMissingOuterHits () >= 3.0; // This is not applied for BG MC
               && probe.lostOuterLayers() >= 3.0; // This is not applied for BG MC
    return passes;
  }
   
  static bool passesLooseMuonVeto (const pat::IsolatedTrack &probe, const std::vector<pat::Muon> &muons)
  {
    bool passes = deltaRToClosestLooseMuon(probe, muons) > 0.15
               // && probe.hitAndTOBDrop_bestTrackMissingOuterHits () >= 3.0; // This is not applied for BG MC
               && probe.lostOuterLayers() >= 3.0; // This is not applied for BG MC          
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

};
  
#endif  // Analysis_Helper_interface_helperFunctions_h