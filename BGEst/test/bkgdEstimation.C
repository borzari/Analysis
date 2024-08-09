#include <vector>
#include <string>

void bkgdEstimation(){

    TCanvas *c1 = new TCanvas("c1", "c1", 700, 700);
    TCanvas *c2 = new TCanvas("c2", "c2", 700, 700);
    TCanvas *c3 = new TCanvas("c3", "c3", 700, 700);
    TCanvas *c4 = new TCanvas("c4", "c4", 700, 700);

    std::string inpElecFile = "hist_merged_2singElecCRNoJetSelSkim.root";
    std::string inpElecMETFile = "hist_merged_2singElecCRNoJetMETTrigSelSkim.root";

    std::string inpMuonFile = "hist_merged_2singMuonCRNoJetSelSkim.root";
    std::string inpMuonMETFile = "hist_merged_2singMuonCRNoJetMETTrigSelSkim.root";

    std::string inpTauFile = "hist_merged_2singTauCRNoJetSelSkim.root";
    std::string inpTauMETFile = "hist_merged_2singTauCRNoJetMETTrigSelSkim.root";

    TFile *fe = TFile::Open(inpElecFile.c_str());
    TFile *fem = TFile::Open(inpElecMETFile.c_str());

    TFile *fm = TFile::Open(inpMuonFile.c_str());
    TFile *fmm = TFile::Open(inpMuonMETFile.c_str());

    TFile *ft = TFile::Open(inpTauFile.c_str());
    TFile *ftm = TFile::Open(inpTauMETFile.c_str());

    TH2D *he;
    TH2D *hem;

    TH2D *hm;
    TH2D *hmm;

    TH2D *ht;
    TH2D *htm;

    fe->GetObject("singElecCRNoJetSelSkimFilter/hist_metNoMuvsDeltaPhiJetMetNoMu;1", he);
    fem->GetObject("singElecCRNoJetSelSkimFilter/hist_metNoMuvsDeltaPhiJetMetNoMu;1", hem);

    fm->GetObject("singMuonCRNoJetSelSkimFilter/hist_metNoMuvsDeltaPhiJetMetNoMu;1", hm);
    fmm->GetObject("singMuonCRNoJetSelSkimFilter/hist_metNoMuvsDeltaPhiJetMetNoMu;1", hmm);

    ft->GetObject("singTauCRNoJetSelSkimFilter/hist_metNoMuvsDeltaPhiJetMetNoMu;1", ht);
    ftm->GetObject("singTauCRNoJetSelSkimFilter/hist_metNoMuvsDeltaPhiJetMetNoMu;1", htm);

    double the = he->Integral();
    double mche = he->Integral(he->GetXaxis()->FindBin(120.0), he->GetNbinsX() + 1, he->GetYaxis()->FindBin(0.5), he->GetNbinsY() + 1);

    double thm = hm->Integral();
    double mchm = hm->Integral(hm->GetXaxis()->FindBin(120.0), hm->GetNbinsX() + 1, hm->GetYaxis()->FindBin(0.5), hm->GetNbinsY() + 1);

    double tht = ht->Integral();
    double mcht = ht->Integral(ht->GetXaxis()->FindBin(120.0), ht->GetNbinsX() + 1, ht->GetYaxis()->FindBin(0.5), ht->GetNbinsY() + 1);

    std::cout << "Electron Poffline: " << mche/the << std::endl;
    std::cout << "Muon Poffline: " << mchm/thm << std::endl;
    std::cout << "Tau Poffline: " << mcht/tht << std::endl;

    he->GetYaxis()->SetRangeUser(0.5, 3.2);
    hem->GetYaxis()->SetRangeUser(0.5, 3.2);
    hm->GetYaxis()->SetRangeUser(0.5, 3.2);
    hmm->GetYaxis()->SetRangeUser(0.5, 3.2);
    ht->GetYaxis()->SetRangeUser(0.5, 3.2);
    htm->GetYaxis()->SetRangeUser(0.5, 3.2);

    TH1D *phe = he->ProjectionX("phe");
    TH1D *phem = hem->ProjectionX("phem");

    TH1D *phm = hm->ProjectionX("phm");
    TH1D *phmm = hmm->ProjectionX("phmm");

    TH1D *pht = ht->ProjectionX("pht");
    TH1D *phtm = htm->ProjectionX("phtm");

    TH1 *rphe = phe->Rebin(4,"rphe");
    TH1 *rphem = phem->Rebin(4,"rphem");
    TH1 *rphm = phm->Rebin(4,"rphm");
    TH1 *rphmm = phmm->Rebin(4,"rphmm");
    TH1 *rpht = pht->Rebin(4,"rpht");
    TH1 *rphtm = phtm->Rebin(4,"rphtm");

    TH1 *effe = (TH1*)rphem->Clone("effe");
    effe->Divide(rphem,rphe,1,1,"pois");
    effe->SetMarkerStyle(20);

    TH1 *effm = (TH1*)rphmm->Clone("effm");
    effm->Divide(rphmm,rphm,1,1,"pois");
    effm->SetMarkerStyle(20);

    c1->cd();
    effe->Draw();
    c1->Update();
    std::string c1Name = "test.pdf";
    c1->Print(c1Name.c_str());

    c2->cd();
    effm->Draw();
    c2->Update();
    c1Name = "test2.pdf";
    c2->Print(c1Name.c_str());

    c3->cd();
    phe->Draw();
    c3->Update();
    c1Name = "test3.pdf";
    c3->Print(c1Name.c_str());

    rphe->GetXaxis()->SetRangeUser(120.0, 1000.0);
    rphm->GetXaxis()->SetRangeUser(120.0, 1000.0);
    rpht->GetXaxis()->SetRangeUser(120.0, 1000.0);

    effe->GetXaxis()->SetRangeUser(120.0, 1000.0);
    effm->GetXaxis()->SetRangeUser(120.0, 1000.0);
    effe->GetXaxis()->SetRangeUser(120.0, 1000.0);

    TH1 *mue = (TH1*)rphe->Clone("mue");
    mue->Multiply(effe,rphe);

    TH1 *mum = (TH1*)rphm->Clone("mum");
    mum->Multiply(effm,rphm);
    
    TH1 *mut = (TH1*)rpht->Clone("mut");
    mut->Multiply(effe,rpht);

    double trphe = rphe->Integral();
    double tmue = mue->Integral();

    double trphm = rphm->Integral();
    double tmum = mum->Integral();

    double trpht = rpht->Integral();
    double tmut = mut->Integral();

    std::cout << "Electron Ptrigger: " << tmue/trphe << std::endl;
    std::cout << "Muon Ptrigger: " << tmum/trphm << std::endl;
    std::cout << "Tau Ptrigger: " << tmut/trpht << std::endl;

}