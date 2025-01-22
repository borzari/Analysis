#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <ctime>
#include <sstream>

void logSpace (const unsigned n, const double a, const double b, std::vector<double> &bins)
{
  double step = (b - a) / ((double) n);

  bins.clear ();
  for (double i = a; i < b + 0.5 * step; i += step)
    bins.push_back (pow (10.0, i));
}

using namespace std;

void analyzeDataBGChain(){

    // std::string year = "";
    // std::string process = "WToLNu";

    std::string year = "2022C";
    std::string process = "Muon";

    std::string inpFile = "/afs/cern.ch/work/b/borzari/mergeHistos/histos_allMuonPt_TrigAnalysis_" + process + year + ".txt";
    cout << inpFile << std::endl;
    std::string outFile = "plots/right_" + process + "_" + year + ".root";

    std::ostringstream os;

    TString Str;
    std::ifstream fpr(inpFile, ios::in);
    if(!fpr.is_open()){
      cout << "List of input files not found!" << endl;
      return;
    }

    std::vector<TString> file_name_vector;
    std::string file_chain;
    while(getline(fpr, file_chain)){
      file_name_vector.push_back(file_chain);
    }
    TChain *t2 = new TChain("triggerDataBGAnalyzer/tree");

    for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++)
    {
      TFile *file = TFile::Open(*listIterator);

      cout << "Adding file:--- " << *listIterator << "--- to the chains" << endl;
      t2->Add(*listIterator);
    }

    TFile *out = new TFile(outFile.c_str(), "RECREATE");

    Int_t hlt_pass_hltMET105Filter;
    Int_t hlt_pass_hltTrk50Filter;

    Int_t hlt_pass_HLT_MET105_IsoTrk50;
    Int_t hlt_pass_HLT_MET120_IsoTrk50;
    Int_t hlt_pass_HLT_PFMET105_IsoTrk50;
    Int_t hlt_pass_HLT_PFMET120_PFMHT120_IDTight;
    Int_t hlt_pass_HLT_PFMET130_PFMHT130_IDTight;
    Int_t hlt_pass_HLT_PFMET140_PFMHT140_IDTight;
    Int_t hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
    Int_t hlt_pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF;
    Int_t hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF;
    Int_t hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF;
    Int_t hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF;
    Int_t hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
    Int_t hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight;
    Int_t hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;
    Int_t hlt_pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60;

    Float_t met_ptNoMu;
    Float_t muon_pt;

    t2->SetBranchAddress("hlt_pass_hltMET105Filter",&hlt_pass_hltMET105Filter);
    t2->SetBranchAddress("hlt_pass_hltTrk50Filter",&hlt_pass_hltTrk50Filter);
    t2->SetBranchAddress("hlt_pass_HLT_MET105_IsoTrk50",&hlt_pass_HLT_MET105_IsoTrk50);

    t2->SetBranchAddress("hlt_pass_HLT_MET120_IsoTrk50",&hlt_pass_HLT_MET120_IsoTrk50);
    t2->SetBranchAddress("hlt_pass_HLT_PFMET105_IsoTrk50",&hlt_pass_HLT_PFMET105_IsoTrk50);
    t2->SetBranchAddress("hlt_pass_HLT_PFMET120_PFMHT120_IDTight",&hlt_pass_HLT_PFMET120_PFMHT120_IDTight);
    t2->SetBranchAddress("hlt_pass_HLT_PFMET130_PFMHT130_IDTight",&hlt_pass_HLT_PFMET130_PFMHT130_IDTight);
    t2->SetBranchAddress("hlt_pass_HLT_PFMET140_PFMHT140_IDTight",&hlt_pass_HLT_PFMET140_PFMHT140_IDTight);
    t2->SetBranchAddress("hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",&hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
    t2->SetBranchAddress("hlt_pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF",&hlt_pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF);
    t2->SetBranchAddress("hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF",&hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF);
    t2->SetBranchAddress("hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF",&hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF);
    t2->SetBranchAddress("hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF",&hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF);
    t2->SetBranchAddress("hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",&hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
    t2->SetBranchAddress("hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight",&hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight);
    t2->SetBranchAddress("hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",&hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight);
    t2->SetBranchAddress("hlt_pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60",&hlt_pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60);

    t2->SetBranchAddress("met_ptNoMu",&met_ptNoMu);
    t2->SetBranchAddress("muon_pt",&muon_pt);

    std::vector<double> bins;
    logSpace (100, 0.0, 3.0, bins);

    TH1F *hmuon_ptlog = new TH1F("hmuon_ptlog","Muon pt",bins.size () - 1, bins.data ());
    TH1F *hmuon_pt = new TH1F("hmuon_pt","Muon pt",100,0.0,1000.0);
    TH1F *hmet_pt = new TH1F("hmet_pt","MET pt",100,0.0,1000.0);
    TH1F *hmuon_pt_hlt = new TH1F("hmuon_pt_hlt","Muon pt",100,0.0,1000.0);
    TH1F *hmet_pt_hlt = new TH1F("hmet_pt_hlt","MET pt",100,0.0,1000.0);
    TH1F *hmet_pt_hltGrandOR = new TH1F("hmet_pt_hltGrandOR","MET pt",100,0.0,1000.0);
    TH1F *hmet_pt_hltGrandOR_except = new TH1F("hmet_pt_hltGrandOR_except","MET pt",100,0.0,1000.0);
    TH1F *hmet_pt_hltGrandOR_exceptIsoTrk = new TH1F("hmet_pt_hltGrandOR_exceptIsoTrk","MET pt",100,0.0,1000.0);
    TH1F *hmuon_pt_hltFilterMet = new TH1F("hmuon_pt_hltFilterMet","Muon pt",100,0.0,1000.0);
    TH1F *hmuon_pt_hltFilterMetlog = new TH1F("hmuon_pt_hltFilterMetlog","Muon pt",bins.size () - 1, bins.data ());
    TH1F *hmet_pt_hltFilterMet = new TH1F("hmet_pt_hltFilterMet","MET pt",100,0.0,1000.0);
    TH1F *hmuon_pt_hltFilterTrklog = new TH1F("hmuon_pt_hltFilterTrklog","Muon pt",bins.size () - 1, bins.data ());
    TH1F *hmuon_pt_hltFilterTrk = new TH1F("hmuon_pt_hltFilterTrk","Muon pt",100,0.0,1000.0);
    TH1F *hmuon_pt_onlyhltFilterTrklog = new TH1F("hmuon_pt_onlyhltFilterTrklog","Muon pt",bins.size () - 1, bins.data ());
    TH1F *hmuon_pt_onlyhltFilterTrk = new TH1F("hmuon_pt_onlyhltFilterTrk","Muon pt",100,0.0,1000.0);
    TH1F *hmuon_pt_hltFilterTrkPassPath = new TH1F("hmuon_pt_hltFilterTrkPassPath","Muon pt",100,0.0,1000.0);
    TH1F *hmet_pt_hltFilterTrk = new TH1F("hmet_pt_hltFilterTrk","MET pt",100,0.0,1000.0);
    TH1F *hmet_pt_onlyhltFilterTrk = new TH1F("hmet_pt_onlyhltFilterTrk","MET pt",100,0.0,1000.0);

    TH1F *hmet_pt_MET120 = new TH1F("hmet_pt_MET120","MET pt",100,0.0,1000.0);
    TH1F *hmet_pt_hlt_MET120 = new TH1F("hmet_pt_hlt_MET120","MET pt",100,0.0,1000.0);

    TH1F *hmet_pt_PFMET105 = new TH1F("hmet_pt_PFMET105","MET pt",100,0.0,1000.0);
    TH1F *hmet_pt_hlt_PFMET105 = new TH1F("hmet_pt_hlt_PFMET105","MET pt",100,0.0,1000.0);

    TH1F *hmet_pt_PFMET110 = new TH1F("hmet_pt_PFMET110","MET pt",100,0.0,1000.0);
    TH1F *hmet_pt_hlt_PFMET110 = new TH1F("hmet_pt_hlt_PFMET110","MET pt",100,0.0,1000.0);

    TH1F *hmet_pt_PFMET120 = new TH1F("hmet_pt_PFMET120","MET pt",100,0.0,1000.0);
    TH1F *hmet_pt_hlt_PFMET120 = new TH1F("hmet_pt_hlt_PFMET120","MET pt",100,0.0,1000.0);

    Int_t nentries = (Int_t)t2->GetEntries();

    std::cout << nentries << std::endl;

    double nGrandOR = 0;
    double nGrandORExcept = 0;
    double nGrandORExceptIsoTrk = 0;

    double nGrandOR120 = 0;
    double nGrandORExcept120 = 0;
    double nGrandORExceptIsoTrk120 = 0;

    for (Int_t i=0;i<nentries;i++) {
        // bhlt_pass_hltMET105Filter->GetEntry(i);
        // bhlt_pass_hltTrk50Filter->GetEntry(i);
        // bhlt_pass_HLT_MET105_IsoTrk50->GetEntry(i);

        // bhlt_pass_HLT_MET120_IsoTrk50->GetEntry(i);
        // bhlt_pass_HLT_PFMET105_IsoTrk50->GetEntry(i);
        // bhlt_pass_HLT_PFMET120_PFMHT120_IDTight->GetEntry(i);
        // bhlt_pass_HLT_PFMET130_PFMHT130_IDTight->GetEntry(i);
        // bhlt_pass_HLT_PFMET140_PFMHT140_IDTight->GetEntry(i);
        // bhlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60->GetEntry(i);
        // bhlt_pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF->GetEntry(i);
        // bhlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF->GetEntry(i);
        // bhlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF->GetEntry(i);
        // bhlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF->GetEntry(i);
        // bhlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight->GetEntry(i);
        // bhlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight->GetEntry(i);
        // bhlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight->GetEntry(i);
        // bhlt_pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60->GetEntry(i);

        // bmet_ptNoMu->GetEntry(i);
        // bmuon_pt->GetEntry(i);

        t2->GetEntry(i);

        hmuon_pt->Fill(muon_pt);
        hmuon_ptlog->Fill(muon_pt);
        hmet_pt->Fill(met_ptNoMu);
        hmet_pt_MET120->Fill(met_ptNoMu);
        hmet_pt_PFMET105->Fill(met_ptNoMu);
        hmet_pt_PFMET110->Fill(met_ptNoMu);
        hmet_pt_PFMET120->Fill(met_ptNoMu);

        if(hlt_pass_HLT_MET105_IsoTrk50 == 1 ||
           hlt_pass_HLT_MET120_IsoTrk50 == 1 || 
           hlt_pass_HLT_PFMET105_IsoTrk50 == 1 || 
           hlt_pass_HLT_PFMET120_PFMHT120_IDTight == 1 || 
           hlt_pass_HLT_PFMET130_PFMHT130_IDTight == 1 || 
           hlt_pass_HLT_PFMET140_PFMHT140_IDTight == 1 || 
           hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 == 1 || 
           hlt_pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF == 1 || 
           hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF == 1 || 
           hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF == 1 || 
           hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF == 1 || 
           hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight == 1 || 
           hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight == 1 || 
           hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight == 1 || 
           hlt_pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60 == 1){hmet_pt_hltGrandOR->Fill(met_ptNoMu); nGrandOR = nGrandOR + 1.0; if(met_ptNoMu > 120.0) nGrandOR120 = nGrandOR120 + 1.0;}

        if(hlt_pass_HLT_MET105_IsoTrk50 == 1 ||
           hlt_pass_HLT_MET120_IsoTrk50 == 1 || 
           hlt_pass_HLT_PFMET105_IsoTrk50 == 1 || 
           hlt_pass_HLT_PFMET120_PFMHT120_IDTight == 1 || 
           hlt_pass_HLT_PFMET130_PFMHT130_IDTight == 1 || 
           hlt_pass_HLT_PFMET140_PFMHT140_IDTight == 1 || 
           hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 == 1 || 
           hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF == 1 || 
           hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF == 1 || 
           hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF == 1 || 
           hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight == 1 || 
           hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight == 1 || 
           hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight == 1 || 
           hlt_pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60 == 1){hmet_pt_hltGrandOR_except->Fill(met_ptNoMu); nGrandORExcept = nGrandORExcept + 1.0; if(met_ptNoMu > 120.0) nGrandORExcept120 = nGrandORExcept120 + 1.0;}

        if(hlt_pass_HLT_PFMET120_PFMHT120_IDTight == 1 || 
           hlt_pass_HLT_PFMET130_PFMHT130_IDTight == 1 || 
           hlt_pass_HLT_PFMET140_PFMHT140_IDTight == 1 || 
           hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 == 1 || 
           hlt_pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF == 1 || 
           hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF == 1 || 
           hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF == 1 || 
           hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF == 1 || 
           hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight == 1 || 
           hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight == 1 || 
           hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight == 1 || 
           hlt_pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60 == 1){hmet_pt_hltGrandOR_exceptIsoTrk->Fill(met_ptNoMu); nGrandORExceptIsoTrk = nGrandORExceptIsoTrk + 1.0; if(met_ptNoMu > 120.0) nGrandORExceptIsoTrk120 = nGrandORExceptIsoTrk120 + 1.0;}

        if(hlt_pass_HLT_MET105_IsoTrk50 == 1) {hmuon_pt_hlt->Fill(muon_pt); hmet_pt_hlt->Fill(met_ptNoMu);}

        if(hlt_pass_HLT_MET120_IsoTrk50 == 1) hmet_pt_hlt_MET120->Fill(met_ptNoMu);

        if(hlt_pass_HLT_PFMET105_IsoTrk50 == 1) hmet_pt_hlt_PFMET105->Fill(met_ptNoMu);

        if(hlt_pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF == 1) hmet_pt_hlt_PFMET110->Fill(met_ptNoMu);

        if(hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF == 1) hmet_pt_hlt_PFMET120->Fill(met_ptNoMu);

        if(hlt_pass_hltMET105Filter == 1) {
            hmuon_pt_hltFilterMet->Fill(muon_pt);
            hmuon_pt_hltFilterMetlog->Fill(muon_pt);
            hmet_pt_hltFilterMet->Fill(met_ptNoMu);
            if(hlt_pass_hltTrk50Filter == 1) {
                hmuon_pt_hltFilterTrk->Fill(muon_pt);
                hmuon_pt_hltFilterTrklog->Fill(muon_pt);
                hmet_pt_hltFilterTrk->Fill(met_ptNoMu);
            }
            if(hlt_pass_HLT_MET105_IsoTrk50 == 1 && muon_pt > 55.0) hmuon_pt_hltFilterTrkPassPath->Fill(muon_pt);
        }

        if(hlt_pass_hltTrk50Filter == 1) {
                hmuon_pt_onlyhltFilterTrk->Fill(muon_pt);
                hmuon_pt_onlyhltFilterTrklog->Fill(muon_pt);
                hmet_pt_onlyhltFilterTrk->Fill(met_ptNoMu);
        }
    }

    int marker = 20;

    // TGraphAsymmErrors *h1 = new TGraphAsymmErrors("h1","MuonPt_HLT_MET105_IsoTrk50");
    // h1->SetMarkerStyle(marker);
    // h1->Divide(hmuon_pt_hlt,hmuon_pt);
    // h1->SetTitle("Muon Pt - HLT_MET105_IsoTrk50");

    // TGraphAsymmErrors *h2 = new TGraphAsymmErrors("h2","MuonPt_hltFilterMet");
    // h2->SetMarkerStyle(marker);
    // h2->Divide(hmuon_pt_hltFilterMet,hmuon_pt);
    // h2->SetTitle("Muon Pt - hltMET105");

    // TGraphAsymmErrors *h3 = new TGraphAsymmErrors("h3","MuonPt_hltFilterTrk");
    // h3->SetMarkerStyle(marker);
    // h3->Divide(hmuon_pt_hltFilterTrk,hmuon_pt_hltFilterMet);
    // h3->SetTitle("Muon Pt - hltTrk50Filter");

    TGraphAsymmErrors *h4 = new TGraphAsymmErrors("h4","MissPt_HLT_MET105_IsoTrk50");
    h4->SetMarkerStyle(marker);
    h4->Divide(hmet_pt_hlt,hmet_pt);
    h4->SetTitle("MET Pt - HLT_MET105_IsoTrk50");

    TGraphAsymmErrors *h5 = new TGraphAsymmErrors("h5","MissPt_hltFilterMet");
    h5->SetMarkerStyle(marker);
    h5->Divide(hmet_pt_hltFilterMet,hmet_pt);
    h5->SetTitle("MET Pt - hltMET105");

    TGraphAsymmErrors *h6 = new TGraphAsymmErrors("h6","MissPt_GrandOR");
    h6->SetMarkerStyle(marker);
    h6->Divide(hmet_pt_hltGrandOR,hmet_pt);
    h6->SetTitle("MET Pt - GrandOR");

    TGraphAsymmErrors *h7 = new TGraphAsymmErrors("h7","MissPt_GrandORExcept");
    h7->SetMarkerStyle(marker);
    h7->Divide(hmet_pt_hltGrandOR_except,hmet_pt);
    h7->SetTitle("MET Pt - GrandORExcept");

    TGraphAsymmErrors *h8 = new TGraphAsymmErrors("h8","MissPt_GrandORExceptIsoTrk");
    h8->SetMarkerStyle(marker);
    h8->Divide(hmet_pt_hltGrandOR_exceptIsoTrk,hmet_pt);
    h8->SetTitle("MET Pt - GrandORExceptIsoTrk");

    // TGraphAsymmErrors *h9 = new TGraphAsymmErrors("h9","MissPt_hltFilterTrk");
    // h9->SetMarkerStyle(marker);
    // h9->Divide(hmet_pt_hltFilterTrk,hmet_pt_hltFilterMet);
    // h9->SetTitle("MET Pt - hltTrk50Filter");

    // TGraphAsymmErrors *h10 = new TGraphAsymmErrors("h10","MuonPt_hltFilterTrklog");
    // h10->SetMarkerStyle(marker);
    // h10->Divide(hmuon_pt_hltFilterTrklog,hmuon_ptlog);
    // h10->SetTitle("Muon Pt - hltTrk50Filter");

    // TGraphAsymmErrors *h11 = new TGraphAsymmErrors("h11","MuonPt_onlyhltFilterTrk");
    // h11->SetMarkerStyle(marker);
    // h11->Divide(hmuon_pt_onlyhltFilterTrk,hmuon_pt);
    // h11->SetTitle("Muon Pt - only hltTrk50Filter");

    // TGraphAsymmErrors *h12 = new TGraphAsymmErrors("h12","MuonPt_onlyhltFilterTrklog");
    // h12->SetMarkerStyle(marker);
    // h12->Divide(hmuon_pt_onlyhltFilterTrklog,hmuon_ptlog);
    // h12->SetTitle("Muon Pt - only hltTrk50Filter");

    // TGraphAsymmErrors *h13 = new TGraphAsymmErrors("h13","MissPt_onlyhltFilterTrklog");
    // h13->SetMarkerStyle(marker);
    // h13->Divide(hmet_pt_onlyhltFilterTrk,hmet_pt);
    // h13->SetTitle("MET Pt - onlyhltTrk50Filter");

    // TGraphAsymmErrors *h14 = new TGraphAsymmErrors("h14","MuonPt_hltFilterTrklog");
    // h14->SetMarkerStyle(marker);
    // h14->Divide(hmuon_pt_hltFilterTrklog,hmuon_pt_hltFilterMetlog);
    // h14->SetTitle("Muon Pt - hltTrk50Filter");

    TGraphAsymmErrors *h15 = new TGraphAsymmErrors("h15","MuonPt_hltFilterTrkPassPath");
    h15->SetMarkerStyle(marker);
    h15->Divide(hmuon_pt_hltFilterTrkPassPath,hmuon_pt_hltFilterMet);
    h15->SetTitle("Muon Pt - hltTrk50Filter");

    TGraphAsymmErrors *h16 = new TGraphAsymmErrors("h16","MissPt_HLT_MET120_IsoTrk50");
    h16->SetMarkerStyle(marker);
    h16->Divide(hmet_pt_hlt_MET120,hmet_pt_MET120);
    h16->SetTitle("MET Pt - HLT_MET120_IsoTrk50");

    TGraphAsymmErrors *h17 = new TGraphAsymmErrors("h16","MissPt_HLT_PFMET105_IsoTrk50");
    h17->SetMarkerStyle(marker);
    h17->Divide(hmet_pt_hlt_PFMET105,hmet_pt_PFMET105);
    h17->SetTitle("MET Pt - HLT_PFMET105_IsoTrk50");

    TGraphAsymmErrors *h18 = new TGraphAsymmErrors("h16","MissPt_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF");
    h18->SetMarkerStyle(marker);
    h18->Divide(hmet_pt_hlt_PFMET110,hmet_pt_PFMET110);
    h18->SetTitle("MET Pt - HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF");

    TGraphAsymmErrors *h19 = new TGraphAsymmErrors("h16","MissPt_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF");
    h19->SetMarkerStyle(marker);
    h19->Divide(hmet_pt_hlt_PFMET120,hmet_pt_PFMET120);
    h19->SetTitle("MET Pt - HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF");

    // h1->Write();
    // h2->Write();
    // h3->Write();
    h4->Write();
    h5->Write();
    h6->Write();
    h7->Write();
    h8->Write();
    // h9->Write();
    // h10->Write();
    // h11->Write();
    // h12->Write();
    // h13->Write();
    // h14->Write();
    h15->Write();
    h16->Write();
    h17->Write();
    h18->Write();
    h19->Write();

    out->Close();

}