#ifndef Analysis_Helper_interface_selectingFunctions_h
#define Analysis_Helper_interface_selectingFunctions_h

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

#include "helperFunctions.h"

class selectingFunctions {

  public:

  static void singElecSel(std::vector<bool> &passSel, std::vector<bool> &passCut, std::vector<bool> &auxPassCut, const pat::Electron &electron, const int startElecIdx, reco::Vertex pv, std::vector<pat::Electron> &selElectrons, const double cutElecPt){

    int cutIdxInc = 0;

    if(electron.pt() > cutElecPt)
      {passSel[startElecIdx+cutIdxInc] = true; if(passCut[startElecIdx-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(abs(electron.eta()) < 2.1)
      {passSel[startElecIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-tight"))
      {passSel[startElecIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(helperFunctions::elecD0(electron, pv))
      {passSel[startElecIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(helperFunctions::elecDZ(electron, pv))
      {passSel[startElecIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    if(auxPassCut[cutIdxInc]) selElectrons.push_back(electron);
      
  }

  static void singMuonSel(std::vector<bool> &passSel, std::vector<bool> &passCut, std::vector<bool> &auxPassCut, const pat::Muon &muon, const int startMuonIdx, reco::Vertex pv, std::vector<pat::Muon> &selMuons){

    int cutIdxInc = 0;
      
    if(muon.pt() > 35.)
      {passSel[startMuonIdx+cutIdxInc] = true; if(passCut[startMuonIdx-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(abs(muon.eta()) < 2.1)
      {passSel[startMuonIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(muon.isTightMuon(pv))
      {passSel[startMuonIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(helperFunctions::muonIso(muon) < 0.15)
      {passSel[startMuonIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    if(auxPassCut[cutIdxInc]) selMuons.push_back(muon);
      
  }

  static void singTauSel(std::vector<bool> &passSel, std::vector<bool> &passCut, std::vector<bool> &auxPassCut, const pat::Tau &tau, const int startTauIdx, std::vector<pat::Tau> &selTaus){

    int cutIdxInc = 0;

    if(tau.pt() > 50.)
      {passSel[startTauIdx+cutIdxInc] = true; if(passCut[startTauIdx-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(abs(tau.eta()) < 2.1)
      {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;
      
    if(helperFunctions::passesDecayModeReconstruction(tau) && helperFunctions::passesLightFlavorRejection(tau))
      {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    if(auxPassCut[cutIdxInc]) selTaus.push_back(tau);
      
  }

  template<class T>
  static void singTrackSel(std::vector<bool> &passSel, std::vector<bool> &passCut, std::vector<bool> &auxPassCut, const pat::IsolatedTrack &track, const int startTrackIdx, const int getStTrkIdxLep, std::vector<T> &tkMatchLeptons, std::vector<T> &selLeptons, edm::Handle<std::vector<T>> &leptons, const double cutTrackPt, const EtaPhiList vetoListElec, const EtaPhiList vetoListMu, std::map<DetId, std::vector<double> > &EcalAllDeadChannelsValMap, std::map<DetId, std::vector<int> > &EcalAllDeadChannelsBitMap, edm::Handle<std::vector<pat::Jet>> &jets, std::vector<pat::IsolatedTrack> &selTracks){

    int cutIdxInc = 0;

    if(track.pt() > cutTrackPt) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(passCut[startTrackIdx+cutIdxInc-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    // These two cutas are calculated differently for cumulative and individual selections, since individual takes into account all leptons,
    // and cumulative only takes into account leptons that passed the cuts

    for(const auto& lepton : selLeptons) {
      if(deltaR (track, lepton) < 0.1)
        {if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true; tkMatchLeptons.push_back(lepton);}
    }
    
    for(const auto& lepton : *leptons) {
      if(deltaR (track, lepton) < 0.1) passSel[startTrackIdx+cutIdxInc] = true;
    }

    ++cutIdxInc;

    if(fabs(track.eta()) < 2.1) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if((fabs(track.eta()) < 1.42) || (fabs(track.eta()) > 1.65)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if((fabs(track.eta()) < 0.15) || (fabs(track.eta()) > 0.35)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if((fabs(track.eta()) < 1.55) || (fabs(track.eta()) > 1.85)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(!helperFunctions::inTOBCrack(track)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(helperFunctions::isFiducialTrack(track,vetoListElec,0.05,-1.0))
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(helperFunctions::isFiducialTrack(track,vetoListMu,0.05,-1.0))
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(!helperFunctions::isCloseToBadEcalChannel(track,0.05,EcalAllDeadChannelsValMap,EcalAllDeadChannelsBitMap))
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(track.hitPattern().numberOfValidPixelHits() >= 4) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(track.hitPattern().numberOfValidHits() >= 4) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(track.lostInnerLayers() == 0) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(helperFunctions::hitDrop_missingMiddleHits(track) == 0) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(((track.pfIsolationDR03().chargedHadronIso() + track.pfIsolationDR03().puChargedHadronIso()) / track.pt()) < 0.05) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(fabs(track.dxy()) < 0.02) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(fabs(track.dz()) < 0.5) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    if(std::is_same<T, pat::Electron>::value || std::is_same<T, pat::Muon>::value){

      ++cutIdxInc;

      if(helperFunctions::dRMinJet(track, *jets) > 0.5)
        {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    }

    ++cutIdxInc;

    // Both cuts below need the track to pass the track selection, that is why the individual efficiency has to take into account the previous cuts

    if(int(tkMatchLeptons.size()) > 0)
      {if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]){passSel[startTrackIdx+cutIdxInc] = true; auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}}

    ++cutIdxInc;

    if(int(tkMatchLeptons.size()) > 0)
      {if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]){passSel[startTrackIdx+cutIdxInc] = true; auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}}

    if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep]) selTracks.push_back(track);

  }

  static void zToLepElecSel(std::vector<bool> &passSel, std::vector<bool> &passCut, std::vector<bool> &auxPassCut, const pat::Electron &electron, const int startElecIdx, reco::Vertex pv, std::vector<pat::Electron> &electronTags, const edm::Event &event, const edm::Handle<edm::TriggerResults> &triggers, const edm::Handle<std::vector<pat::TriggerObjectStandAlone> > &trigObjs){

    int cutIdxInc = 0;

    if(electron.pt() > 32.)
      {passSel[startElecIdx+cutIdxInc] = true; if(passCut[startElecIdx-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(helperFunctions::isMatchedToTriggerObject<pat::Electron> (event, *triggers, electron, *trigObjs, "hltEgammaCandidates::HLT", "hltEle32WPTightGsfTrackIsoFilter"))
      {passSel[startElecIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(abs(electron.eta()) < 2.1)
      {passSel[startElecIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-tight"))
      {passSel[startElecIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(helperFunctions::elecD0(electron, pv))
      {passSel[startElecIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(helperFunctions::elecDZ(electron, pv))
      {passSel[startElecIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    if(auxPassCut[cutIdxInc]) electronTags.push_back(electron);
      
  }

  static void zToLepMuonSel(std::vector<bool> &passSel, std::vector<bool> &passCut, std::vector<bool> &auxPassCut, const pat::Muon &muon, const int startMuonIdx, reco::Vertex pv, std::vector<pat::Muon> &muonTags, const edm::Event &event, const edm::Handle<edm::TriggerResults> &triggers, const edm::Handle<std::vector<pat::TriggerObjectStandAlone> > &trigObjs){

    int cutIdxInc = 0;
      
    if(muon.pt() > 26.)
      {passSel[startMuonIdx+cutIdxInc] = true; if(passCut[startMuonIdx-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(helperFunctions::isMatchedToTriggerObject<pat::Muon> (event, *triggers, muon, *trigObjs, "hltIterL3MuonCandidates::HLT", "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered"))
      {passSel[startMuonIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(abs(muon.eta()) < 2.1)
      {passSel[startMuonIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(muon.isTightMuon(pv))
      {passSel[startMuonIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(helperFunctions::muonIso(muon) < 0.15)
      {passSel[startMuonIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    if(auxPassCut[cutIdxInc]) muonTags.push_back(muon);
      
  }

  static void zToLepTauEleSel(std::vector<bool> &passSel, std::vector<bool> &passCut, std::vector<bool> &auxPassCut, const pat::Electron &electron, const int startTauIdx, reco::Vertex pv, std::vector<pat::Electron> &taueTags, const edm::Event &event, const edm::Handle<edm::TriggerResults> &triggers, const edm::Handle<std::vector<pat::TriggerObjectStandAlone> > &trigObjs, const pat::MET &met){

    int cutIdxInc = 0;

    if(electron.pt() > 32.)
      {passSel[startTauIdx+cutIdxInc] = true; if(passCut[startTauIdx-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(helperFunctions::isMatchedToTriggerObject<pat::Electron> (event, *triggers, electron, *trigObjs, "hltEgammaCandidates::HLT", "hltEle32WPTightGsfTrackIsoFilter"))
      {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(abs(electron.eta()) < 2.1)
      {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(electron.electronID("cutBasedElectronID-RunIIIWinter22-V1-tight"))
      {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(helperFunctions::elecD0(electron, pv))
      {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(helperFunctions::elecDZ(electron, pv))
      {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(helperFunctions::transvMassLepton<pat::Electron>(electron,met) < 40.0)
      {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    if(auxPassCut[cutIdxInc]) taueTags.push_back(electron);
      
  }

  static void zToLepTauMuSel(std::vector<bool> &passSel, std::vector<bool> &passCut, std::vector<bool> &auxPassCut, const pat::Muon &muon, const int startTauIdx, reco::Vertex pv, std::vector<pat::Muon> &taumTags, const edm::Event &event, const edm::Handle<edm::TriggerResults> &triggers, const edm::Handle<std::vector<pat::TriggerObjectStandAlone> > &trigObjs, const pat::MET &met){

    int cutIdxInc = 0;
      
    if(muon.pt() > 26.)
      {passSel[startTauIdx+cutIdxInc] = true; if(passCut[startTauIdx-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(helperFunctions::isMatchedToTriggerObject<pat::Muon> (event, *triggers, muon, *trigObjs, "hltIterL3MuonCandidates::HLT", "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered"))
      {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(abs(muon.eta()) < 2.1)
      {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(muon.isTightMuon(pv))
      {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(helperFunctions::muonIso(muon) < 0.15)
      {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    ++cutIdxInc;

    if(helperFunctions::transvMassLepton<pat::Muon>(muon,met) < 40.0)
      {passSel[startTauIdx+cutIdxInc] = true; if(auxPassCut[cutIdxInc-1]) auxPassCut[cutIdxInc] = true;}

    if(auxPassCut[cutIdxInc]) taumTags.push_back(muon);
      
  }

  template<class T>
  static void zToLepTrackSel(std::vector<bool> &passSel, std::vector<bool> &passCut, std::vector<bool> &auxPassCut, const pat::IsolatedTrack &track, const int startTrackIdx, const int getStTrkIdxLep, std::vector<pat::IsolatedTrack> &trackProbes, const EtaPhiList vetoListElec, const EtaPhiList vetoListMu, std::map<DetId, std::vector<double> > &EcalAllDeadChannelsValMap, std::map<DetId, std::vector<int> > &EcalAllDeadChannelsBitMap, edm::Handle<std::vector<pat::Jet>> &jets, const edm::Handle<std::vector<pat::Electron>> &electrons, const edm::Handle<std::vector<pat::Muon>> &muons, const edm::Handle<std::vector<pat::Tau>> &taus, const edm::Handle<EBRecHitCollection> &EBRecHits, const edm::Handle<EERecHitCollection> &EERecHits, const edm::Handle<HBHERecHitCollection> &HBHERecHits, const edm::Handle<double> &rhoCentralCalo, const edm::ESHandle<CaloGeometry> &caloGeometry){

    int cutIdxInc = 0;
      
    if(track.pt() > 30.) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(passCut[startTrackIdx+cutIdxInc-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(fabs(track.eta()) < 2.1) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if((fabs(track.eta()) < 1.42) || (fabs(track.eta()) > 1.65)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if((fabs(track.eta()) < 0.15) || (fabs(track.eta()) > 0.35)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if((fabs(track.eta()) < 1.55) || (fabs(track.eta()) > 1.85)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(!helperFunctions::inTOBCrack(track)) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(helperFunctions::isFiducialTrack(track,vetoListElec,0.05,-1.0))
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(helperFunctions::isFiducialTrack(track,vetoListMu,0.05,-1.0))
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(!helperFunctions::isCloseToBadEcalChannel(track,0.05,EcalAllDeadChannelsValMap,EcalAllDeadChannelsBitMap))
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(track.hitPattern().numberOfValidPixelHits() >= 4) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(track.hitPattern().numberOfValidHits() >= 4) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(track.lostInnerLayers() == 0) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(helperFunctions::hitDrop_missingMiddleHits(track) == 0) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(((track.pfIsolationDR03().chargedHadronIso() + track.pfIsolationDR03().puChargedHadronIso()) / track.pt()) < 0.05) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(fabs(track.dxy()) < 0.02) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    ++cutIdxInc;

    if(fabs(track.dz()) < 0.5) 
      {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    if(std::is_same<T, pat::Electron>::value || std::is_same<T, pat::Muon>::value)
      {
        ++cutIdxInc;

        if(helperFunctions::dRMinJet(track, *jets) > 0.5)
          {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}
      }

    if(std::is_same<T, pat::Electron>::value || std::is_same<T, pat::Tau>::value)
      {
        ++cutIdxInc;

        if(helperFunctions::deltaRToClosestLepton<pat::Muon>(track, *muons) > 0.15)
          {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}
      }

    if(std::is_same<T, pat::Muon>::value || std::is_same<T, pat::Tau>::value)
      {
        ++cutIdxInc;

        if(helperFunctions::deltaRToClosestLepton<pat::Electron>(track, *electrons) > 0.15)
          {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}
      }

    if(std::is_same<T, pat::Electron>::value || std::is_same<T, pat::Muon>::value)
      {
        ++cutIdxInc;

        if(helperFunctions::deltaRToClosestTauHad(track, *taus) > 0.15)
          {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

      } // The individual efficiency of this is different than the original analysis, because it uses distinct selections that follow the Run 3 recommendations from the Tau POG https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3#Kinematic_tau_selection
    
    if(std::is_same<T, pat::Muon>::value){

      ++cutIdxInc;

      if(helperFunctions::caloNewNoPUDRp5CentralCalo(track, *EBRecHits, *EERecHits, *HBHERecHits, *rhoCentralCalo, caloGeometry) < 10.0) 
        {passSel[startTrackIdx+cutIdxInc] = true; if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep-1]) auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep] = true;}

    }

    if(auxPassCut[startTrackIdx+cutIdxInc-getStTrkIdxLep]) trackProbes.push_back(track);

  }

};

#endif  // Analysis_Helper_interface_selectingFunctions_h