#include <vector>
#include <string>

#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH2D.h"

void bkgdEstimation(){

    std::string inpVetoElecFile = "hist_merged_4zToElecProbeTrk.root";
    std::string inpSingElecFile = "hist_merged_plotterSingElecCRNoJetSelSkim.root";
    std::string inpElecMETFile = "hist_merged_plotterSingElecCRNoJetMETTrigSelSkim.root";

    std::string inpVetoMuonFile = "hist_merged_2zToMuonProbeTrk.root";
    std::string inpSingMuonFile = "hist_merged_plotterSingMuonCRNoJetSelSkim.root";
    std::string inpMuonMETFile = "hist_merged_plotterSingMuonCRNoJetMETTrigSelSkim.root";

    std::string inpVetoTauEFile = "hist_merged_2zToTauEProbeTrk.root";
    std::string inpVetoTauMFile = "hist_merged_2zToTauMProbeTrk.root";
    std::string inpSingTauFile = "hist_merged_2singTauCRNoJetSelSkim.root";
    std::string inpTauMETFile = "hist_merged_2singTauCRNoJetMETTrigSelSkim.root";

    TFile *fve = TFile::Open(inpVetoElecFile.c_str());
    TFile *fse = TFile::Open(inpSingElecFile.c_str());
    TFile *fem = TFile::Open(inpElecMETFile.c_str());

    TFile *fvm = TFile::Open(inpVetoMuonFile.c_str());
    TFile *fsm = TFile::Open(inpSingMuonFile.c_str());
    TFile *fmm = TFile::Open(inpMuonMETFile.c_str());

    TFile *fvte = TFile::Open(inpVetoTauEFile.c_str());
    TFile *fvtm = TFile::Open(inpVetoTauMFile.c_str());
    TFile *fst = TFile::Open(inpSingTauFile.c_str());
    TFile *ftm = TFile::Open(inpTauMETFile.c_str());

    TH1D *nTPOS_elec;
    TH1D *nTPSS_elec;
    TH1D *nTPOS_veto_elec;
    TH1D *nTPSS_veto_elec;
    TH1D *metNoMu_elec;
    TH1D *metNoMu_elecMETTrig;
    TH2D *metNoMuNoLep_elec;

    TH1D *nTPOS_muon;
    TH1D *nTPSS_muon;
    TH1D *nTPOS_veto_muon;
    TH1D *nTPSS_veto_muon;
    TH1D *metNoMu_muon;
    TH1D *metNoMu_muonMETTrig;
    TH2D *metNoMuNoLep_muon;

    TH1D *nTPOS_taue;
    TH1D *nTPSS_taue;
    TH1D *nTPOS_veto_taue;
    TH1D *nTPSS_veto_taue;
    TH1D *nTPOS_taum;
    TH1D *nTPSS_taum;
    TH1D *nTPOS_veto_taum;
    TH1D *nTPSS_veto_taum;
    TH2D *ht;
    TH2D *htm;

    fve->GetObject("zToElecProbeTrkFilter/hist_nTPOS;1", nTPOS_elec);
    fve->GetObject("zToElecProbeTrkFilter/hist_nTPSS;1", nTPSS_elec);
    fve->GetObject("zToElecProbeTrkFilter/hist_nTPOS_veto;1", nTPOS_veto_elec);
    fve->GetObject("zToElecProbeTrkFilter/hist_nTPSS_veto;1", nTPSS_veto_elec);
    fse->GetObject("plotterSingElecCRNoJetSelSkimFilter/metNoMu_pt;1", metNoMu_elec);
    fem->GetObject("plotterSingElecCRNoJetSelSkimFilter/metNoMu_pt;1", metNoMu_elecMETTrig);
    fse->GetObject("plotterSingElecCRNoJetSelSkimFilter/metNoMuvsDeltaPhiJetMetNoMuNoLep;1", metNoMuNoLep_elec);

    fvm->GetObject("zToMuonProbeTrkFilter/hist_nTPOS;1", nTPOS_muon);
    fvm->GetObject("zToMuonProbeTrkFilter/hist_nTPSS;1", nTPSS_muon);
    fvm->GetObject("zToMuonProbeTrkFilter/hist_nTPOS_veto;1", nTPOS_veto_muon);
    fvm->GetObject("zToMuonProbeTrkFilter/hist_nTPSS_veto;1", nTPSS_veto_muon);
    fsm->GetObject("plotterSingMuonCRNoJetSelSkimFilter/metNoMu_pt;1", metNoMu_muon);
    fmm->GetObject("plotterSingMuonCRNoJetSelSkimFilter/metNoMu_pt;1", metNoMu_muonMETTrig);
    fsm->GetObject("plotterSingMuonCRNoJetSelSkimFilter/metNoMuvsDeltaPhiJetMetNoMuNoLep;1", metNoMuNoLep_muon);

    fvte->GetObject("zToTauEleProbeTrkFilter/hist_nTPOS;1", nTPOS_taue);
    fvte->GetObject("zToTauEleProbeTrkFilter/hist_nTPSS;1", nTPSS_taue);
    fvte->GetObject("zToTauEleProbeTrkFilter/hist_nTPOS_veto;1", nTPOS_veto_taue);
    fvte->GetObject("zToTauEleProbeTrkFilter/hist_nTPSS_veto;1", nTPSS_veto_taue);
    fvtm->GetObject("zToTauMuProbeTrkFilter/hist_nTPOS;1", nTPOS_taum);
    fvtm->GetObject("zToTauMuProbeTrkFilter/hist_nTPSS;1", nTPSS_taum);
    fvtm->GetObject("zToTauMuProbeTrkFilter/hist_nTPOS_veto;1", nTPOS_veto_taum);
    fvtm->GetObject("zToTauMuProbeTrkFilter/hist_nTPSS_veto;1", nTPSS_veto_taum);
    fst->GetObject("singTauCRNoJetSelSkimFilter/hist_metNoMuvsDeltaPhiJetMetNoMu;1", ht);
    ftm->GetObject("singTauCRNoJetSelSkimFilter/hist_metNoMuvsDeltaPhiJetMetNoMu;1", htm);

    double elec_nTPOS = nTPOS_elec->GetBinContent(2) * 1.0 + nTPOS_elec->GetBinContent(3) * 2.0 + nTPOS_elec->GetBinContent(4) * 3.0;
    double elec_nTPSS = nTPSS_elec->GetBinContent(2) * 1.0 + nTPSS_elec->GetBinContent(3) * 2.0 + nTPSS_elec->GetBinContent(4) * 3.0;
    double elec_nTPOS_veto = nTPOS_veto_elec->GetBinContent(2) * 1.0 + nTPOS_veto_elec->GetBinContent(3) * 2.0 + nTPOS_veto_elec->GetBinContent(4) * 3.0;
    double elec_nTPSS_veto = nTPSS_veto_elec->GetBinContent(2) * 1.0 + nTPSS_veto_elec->GetBinContent(3) * 2.0 + nTPSS_veto_elec->GetBinContent(4) * 3.0;

    double muon_nTPOS = nTPOS_muon->GetBinContent(2) * 1.0 + nTPOS_muon->GetBinContent(3) * 2.0 + nTPOS_muon->GetBinContent(4) * 3.0;
    double muon_nTPSS = nTPSS_muon->GetBinContent(2) * 1.0 + nTPSS_muon->GetBinContent(3) * 2.0 + nTPSS_muon->GetBinContent(4) * 3.0;
    double muon_nTPOS_veto = nTPOS_veto_muon->GetBinContent(2) * 1.0 + nTPOS_veto_muon->GetBinContent(3) * 2.0 + nTPOS_veto_muon->GetBinContent(4) * 3.0;
    double muon_nTPSS_veto = nTPSS_veto_muon->GetBinContent(2) * 1.0 + nTPSS_veto_muon->GetBinContent(3) * 2.0 + nTPSS_veto_muon->GetBinContent(4) * 3.0;

    double tau_nTPOS = nTPOS_taue->GetBinContent(2) * 1.0 + nTPOS_taue->GetBinContent(3) * 2.0 + nTPOS_taue->GetBinContent(4) * 3.0 + nTPOS_taum->GetBinContent(2) * 1.0 + nTPOS_taum->GetBinContent(3) * 2.0 + nTPOS_taum->GetBinContent(4) * 3.0;
    double tau_nTPSS = nTPSS_taue->GetBinContent(2) * 1.0 + nTPSS_taue->GetBinContent(3) * 2.0 + nTPSS_taue->GetBinContent(4) * 3.0 + nTPSS_taum->GetBinContent(2) * 1.0 + nTPSS_taum->GetBinContent(3) * 2.0 + nTPSS_taum->GetBinContent(4) * 3.0;
    double tau_nTPOS_veto = nTPOS_veto_taue->GetBinContent(2) * 1.0 + nTPOS_veto_taue->GetBinContent(3) * 2.0 + nTPOS_veto_taue->GetBinContent(4) * 3.0 + nTPOS_veto_taum->GetBinContent(2) * 1.0 + nTPOS_veto_taum->GetBinContent(3) * 2.0 + nTPOS_veto_taum->GetBinContent(4) * 3.0;
    double tau_nTPSS_veto = nTPSS_veto_taue->GetBinContent(2) * 1.0 + nTPSS_veto_taue->GetBinContent(3) * 2.0 + nTPSS_veto_taue->GetBinContent(4) * 3.0 + nTPSS_veto_taum->GetBinContent(2) * 1.0 + nTPSS_veto_taum->GetBinContent(3) * 2.0 + nTPSS_veto_taum->GetBinContent(4) * 3.0;

    double nCtrl_elec = metNoMuNoLep_elec->Integral();
    double nCtrl_muon = metNoMuNoLep_muon->Integral();
    double nCtrl_tau = ht->Integral();

    std::cout << "Electron Nctrl: " << nCtrl_elec << std::endl;
    std::cout << "Muon Nctrl: " << nCtrl_muon << std::endl;
    std::cout << "Tau Nctrl: " << nCtrl_tau << std::endl;

    std::cout << "Electron Pveto: " << (elec_nTPOS_veto - elec_nTPSS_veto)/(elec_nTPOS - elec_nTPSS) << std::endl;
    std::cout << "Muon Pveto: " << (muon_nTPOS_veto - muon_nTPSS_veto)/(muon_nTPOS - muon_nTPSS) << std::endl;
    std::cout << "Tau Pveto: " << (tau_nTPOS_veto - tau_nTPSS_veto)/(tau_nTPOS - tau_nTPSS) << std::endl;

    double the = metNoMuNoLep_elec->Integral();
    double mche = metNoMuNoLep_elec->Integral(metNoMuNoLep_elec->GetXaxis()->FindBin(120.0), metNoMuNoLep_elec->GetNbinsX() + 1, metNoMuNoLep_elec->GetYaxis()->FindBin(0.5), metNoMuNoLep_elec->GetNbinsY() + 1);

    double thm = metNoMuNoLep_muon->Integral();
    double mchm = metNoMuNoLep_muon->Integral(metNoMuNoLep_muon->GetXaxis()->FindBin(120.0), metNoMuNoLep_muon->GetNbinsX() + 1, metNoMuNoLep_muon->GetYaxis()->FindBin(0.5), metNoMuNoLep_muon->GetNbinsY() + 1);

    double tht = ht->Integral();
    double mcht = ht->Integral(ht->GetXaxis()->FindBin(120.0), ht->GetNbinsX() + 1, ht->GetYaxis()->FindBin(0.5), ht->GetNbinsY() + 1);

    TH1D *pht = ht->ProjectionX("pht");

    std::cout << "Electron Poffline: " << mche/the << std::endl;
    std::cout << "Muon Poffline: " << mchm/thm << std::endl;
    std::cout << "Tau Poffline: " << mcht/tht << std::endl;

    TH1 *rphe = metNoMu_elec->Rebin(4,"rphe");
    TH1 *rphem = metNoMu_elecMETTrig->Rebin(4,"rphem");
    TH1 *rphm = metNoMu_muon->Rebin(4,"rphm");
    TH1 *rphmm = metNoMu_muonMETTrig->Rebin(4,"rphmm");
    TH1 *rpht = pht->Rebin(4,"rpht");

    TH1 *effe = (TH1*)rphem->Clone("effe");
    effe->Divide(rphem,rphe,1,1,"pois");
    effe->SetMarkerStyle(20);

    TH1 *effm = (TH1*)rphmm->Clone("effm");
    effm->Divide(rphmm,rphm,1,1,"pois");
    effm->SetMarkerStyle(20);

    // c1->cd();
    // effe->Draw();
    // c1->Update();
    // std::string c1Name = "test.pdf";
    // c1->Print(c1Name.c_str());

    // c2->cd();
    // effm->Draw();
    // c2->Update();
    // c1Name = "test2.pdf";
    // c2->Print(c1Name.c_str());

    metNoMuNoLep_elec->GetXaxis()->SetRangeUser(120.0, 1000.0);
    metNoMuNoLep_muon->GetXaxis()->SetRangeUser(120.0, 1000.0);
    metNoMuNoLep_elec->GetYaxis()->SetRangeUser(0.5, 3.2);
    metNoMuNoLep_muon->GetYaxis()->SetRangeUser(0.5, 3.2);

    TH1D *phe = metNoMuNoLep_elec->ProjectionX("phe");
    TH1D *phm = metNoMuNoLep_muon->ProjectionX("phm");

    TH1 *rrphe = phe->Rebin(4,"rrphe");
    TH1 *rrphm = phm->Rebin(4,"rrphm");

    effe->GetXaxis()->SetRangeUser(120.0, 1000.0);
    effm->GetXaxis()->SetRangeUser(120.0, 1000.0);

    double elecNum = 0.0;
    double muonNum = 0.0;
    double tauNum = 0.0;
    double elecDen = 0.0;
    double muonDen = 0.0;
    double tauDen = 0.0;

    for(int i = 1; i <= int(rrphm->GetNbinsX()); i++){
        elecNum += rrphe->GetBinContent(i) * effe->GetBinContent(i+6);
        elecDen += rrphe->GetBinContent(i);
        muonNum += rrphm->GetBinContent(i) * effm->GetBinContent(i+6);
        muonDen += rrphm->GetBinContent(i);
        tauNum += rpht->GetBinContent(i) * effe->GetBinContent(i+6);
        tauDen += rpht->GetBinContent(i);
    }

    std::cout << "Electron Ptrigger: " << elecNum/elecDen << std::endl;
    std::cout << "Muon Ptrigger: " << muonNum/muonDen << std::endl;
    std::cout << "Tau Ptrigger: " << tauNum/tauDen << std::endl;

    double nEst_elec = (nCtrl_elec * 407.43 * 1000 * 41.5 / 270699232.0) * ((elec_nTPOS_veto - elec_nTPSS_veto)/(elec_nTPOS - elec_nTPSS)) * (mche/the) * (elecNum/elecDen);
    double nEst_muon = (nCtrl_muon * 407.43 * 1000 * 41.5 / 270699232.0) * ((muon_nTPOS_veto - muon_nTPSS_veto)/(muon_nTPOS - muon_nTPSS)) * (mchm/thm) * (muonNum/muonDen);
    double nEst_tau = (nCtrl_tau * 407.43 * 1000 * 41.5 / 270699232.0) * ((tau_nTPOS_veto - tau_nTPSS_veto)/(tau_nTPOS - tau_nTPSS)) * (mcht/tht) * (tauNum/tauDen);

    std::cout << "Electron Nest: " << nEst_elec << std::endl;
    std::cout << "Muon Nest: " << nEst_muon << std::endl;
    std::cout << "Tau Nest: " << nEst_tau << std::endl;

}