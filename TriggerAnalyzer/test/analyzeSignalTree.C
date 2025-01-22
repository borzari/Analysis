#include "tdrstyle.C"
#include "CMS_lumi.C"

double errorCalc(Int_t x, Int_t y){
    double dx = sqrt(x);
    double dy = sqrt(y);
    return (1/(double) y)*(dx + (((double) x/(double) y)*dy));
}

void analyzeSignalTree(){

    std::string mass = "900";
    std::string lifetime = "10";

    std::string inpFile = "plots/outputFile_" + mass + "_" + lifetime + "_2023.root";
    std::string outFile = "plots/AMSB_" + mass + "_" + lifetime + "_2023.root";

    TFile *f = new TFile(inpFile.c_str());
    TTree *t2 = (TTree*)f->Get("triggerSignalAnalyzer/tree");

    TFile *out = new TFile(outFile.c_str(), "RECREATE");
    TCanvas *c1 = new TCanvas("c1", "c1", 700, 700);
    TCanvas *c2 = new TCanvas("c2", "c2", 700, 700);
    TCanvas *c3 = new TCanvas("c3", "c3", 700, 700);
    TCanvas *c4 = new TCanvas("c4", "c4", 700, 700);
    TCanvas *c5 = new TCanvas("c5", "c5", 700, 700);
    TCanvas *c6 = new TCanvas("c6", "c6", 700, 700);
    TCanvas *c7 = new TCanvas("c7", "c7", 700, 700);
    TCanvas *c8 = new TCanvas("c8", "c8", 700, 700);
    TCanvas *cesp = new TCanvas("cesp", "cesp", 700, 700);
    TCanvas *cesp2 = new TCanvas("cesp2", "cesp2", 700, 700);
    TCanvas *cesp3 = new TCanvas("cesp3", "cesp3", 700, 700);

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
    Int_t muon_nLayers;

    TBranch *bhlt_pass_hltMET105Filter = t2->GetBranch("hlt_pass_hltMET105Filter");
    TBranch *bhlt_pass_hltTrk50Filter = t2->GetBranch("hlt_pass_hltTrk50Filter");
    TBranch *bhlt_pass_HLT_MET105_IsoTrk50 = t2->GetBranch("hlt_pass_HLT_MET105_IsoTrk50");

    TBranch *bhlt_pass_HLT_MET120_IsoTrk50 = t2->GetBranch("hlt_pass_HLT_MET120_IsoTrk50");
    TBranch *bhlt_pass_HLT_PFMET105_IsoTrk50 = t2->GetBranch("hlt_pass_HLT_PFMET105_IsoTrk50");
    TBranch *bhlt_pass_HLT_PFMET120_PFMHT120_IDTight = t2->GetBranch("hlt_pass_HLT_PFMET120_PFMHT120_IDTight");
    TBranch *bhlt_pass_HLT_PFMET130_PFMHT130_IDTight = t2->GetBranch("hlt_pass_HLT_PFMET130_PFMHT130_IDTight");
    TBranch *bhlt_pass_HLT_PFMET140_PFMHT140_IDTight = t2->GetBranch("hlt_pass_HLT_PFMET140_PFMHT140_IDTight");
    TBranch *bhlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = t2->GetBranch("hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60");
    TBranch *bhlt_pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF = t2->GetBranch("hlt_pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF");
    TBranch *bhlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF = t2->GetBranch("hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF");
    TBranch *bhlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF = t2->GetBranch("hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF");
    TBranch *bhlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF = t2->GetBranch("hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF");
    TBranch *bhlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = t2->GetBranch("hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight");
    TBranch *bhlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight = t2->GetBranch("hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight");
    TBranch *bhlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight = t2->GetBranch("hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight");
    TBranch *bhlt_pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60 = t2->GetBranch("hlt_pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60");

    TBranch *bmet_ptNoMu = t2->GetBranch("met_ptNoMu");
    TBranch *bmuon_pt = t2->GetBranch("track_pt");
    TBranch *bmuon_nLayers = t2->GetBranch("track_numValidLayers");

    bhlt_pass_hltMET105Filter->SetAddress(&hlt_pass_hltMET105Filter);
    bhlt_pass_hltTrk50Filter->SetAddress(&hlt_pass_hltTrk50Filter);
    bhlt_pass_HLT_MET105_IsoTrk50->SetAddress(&hlt_pass_HLT_MET105_IsoTrk50);

    bhlt_pass_HLT_MET120_IsoTrk50->SetAddress(&hlt_pass_HLT_MET120_IsoTrk50);
    bhlt_pass_HLT_PFMET105_IsoTrk50->SetAddress(&hlt_pass_HLT_PFMET105_IsoTrk50);
    bhlt_pass_HLT_PFMET120_PFMHT120_IDTight->SetAddress(&hlt_pass_HLT_PFMET120_PFMHT120_IDTight);
    bhlt_pass_HLT_PFMET130_PFMHT130_IDTight->SetAddress(&hlt_pass_HLT_PFMET130_PFMHT130_IDTight);
    bhlt_pass_HLT_PFMET140_PFMHT140_IDTight->SetAddress(&hlt_pass_HLT_PFMET140_PFMHT140_IDTight);
    bhlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60->SetAddress(&hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
    bhlt_pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF->SetAddress(&hlt_pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF);
    bhlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF->SetAddress(&hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF);
    bhlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF->SetAddress(&hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF);
    bhlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF->SetAddress(&hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF);
    bhlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight->SetAddress(&hlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
    bhlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight->SetAddress(&hlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight);
    bhlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight->SetAddress(&hlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight);
    bhlt_pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60->SetAddress(&hlt_pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60);

    bmet_ptNoMu->SetAddress(&met_ptNoMu);
    bmuon_pt->SetAddress(&muon_pt); 
    bmuon_nLayers->SetAddress(&muon_nLayers);  

    TH1F *hmuon_pt = new TH1F("hmuon_pt","Muon pt",100,0.0,1000.0);
    TH1F *hmet_pt = new TH1F("hmet_pt","MET pt",100,0.0,1000.0);
    TH1F *hmuon_pt_hlt = new TH1F("hmuon_pt_hlt","Muon pt",100,0.0,1000.0);
    TH1F *hmet_pt_hlt = new TH1F("hmet_pt_hlt","MET pt",100,0.0,1000.0);
    TH1F *hmet_pt_hltGrandOR = new TH1F("hmet_pt_hltGrandOR","MET pt",100,0.0,1000.0);
    TH1F *hmet_pt_hltGrandOR_except = new TH1F("hmet_pt_hltGrandOR_except","MET pt",100,0.0,1000.0);
    TH1F *hmet_pt_hltGrandOR_exceptIsoTrk = new TH1F("hmet_pt_hltGrandOR_exceptIsoTrk","MET pt",100,0.0,1000.0);
    TH1F *hmuon_pt_hltFilterMet = new TH1F("hmuon_pt_hltFilterMet","Muon pt",100,0.0,1000.0);
    TH1F *hmet_pt_hltFilterMet = new TH1F("hmet_pt_hltFilterMet","MET pt",100,0.0,1000.0);
    TH1F *hmuon_pt_hltFilterTrk = new TH1F("hmuon_pt_hltFilterTrk","Muon pt",100,0.0,1000.0);
    Int_t nentries = (Int_t)t2->GetEntries();

    double nGrandOR = 0;
    double nGrandORExcept = 0;
    double nGrandORExceptIsoTrk = 0;

    double nGrandOR120 = 0;
    double nGrandORExcept120 = 0;
    double nGrandORExceptIsoTrk120 = 0;

    double nhlt4Layers = 0;
    double n4Layers = 0;
    double nhlt5Layers = 0;
    double n5Layers = 0;
    double nhlt6Layers = 0;
    double n6Layers = 0;
    double nhltAllLayers = 0;
    double nAllLayers = 0;

    double nhlt4LayersExceptIsoTrk = 0;
    double nhlt5LayersExceptIsoTrk = 0;
    double nhlt6LayersExceptIsoTrk = 0;
    double nhltAllLayersExceptIsoTrk = 0;
    double nhlt4LayersOnlyIsoTrk = 0;
    double nhlt5LayersOnlyIsoTrk = 0;
    double nhlt6LayersOnlyIsoTrk = 0;
    double nhltAllLayersOnlyIsoTrk = 0;

    for (Int_t i=0;i<nentries;i++) {
        bhlt_pass_hltMET105Filter->GetEntry(i);
        bhlt_pass_hltTrk50Filter->GetEntry(i);
        bhlt_pass_HLT_MET105_IsoTrk50->GetEntry(i);

        bhlt_pass_HLT_MET120_IsoTrk50->GetEntry(i);
        bhlt_pass_HLT_PFMET105_IsoTrk50->GetEntry(i);
        bhlt_pass_HLT_PFMET120_PFMHT120_IDTight->GetEntry(i);
        bhlt_pass_HLT_PFMET130_PFMHT130_IDTight->GetEntry(i);
        bhlt_pass_HLT_PFMET140_PFMHT140_IDTight->GetEntry(i);
        bhlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60->GetEntry(i);
        bhlt_pass_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF->GetEntry(i);
        bhlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF->GetEntry(i);
        bhlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF->GetEntry(i);
        bhlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF->GetEntry(i);
        bhlt_pass_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight->GetEntry(i);
        bhlt_pass_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight->GetEntry(i);
        bhlt_pass_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight->GetEntry(i);
        bhlt_pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60->GetEntry(i);

        bmet_ptNoMu->GetEntry(i);
        bmuon_pt->GetEntry(i);
        bmuon_nLayers->GetEntry(i);

        hmuon_pt->Fill(muon_pt);
        hmet_pt->Fill(met_ptNoMu);

        if(met_ptNoMu > 120.0 && muon_nLayers == 4) n4Layers = n4Layers + 1.0;
        if(met_ptNoMu > 120.0 && muon_nLayers == 5) n5Layers = n5Layers + 1.0;
        if(met_ptNoMu > 120.0 && muon_nLayers >= 6) n6Layers = n6Layers + 1.0;
        if(met_ptNoMu > 120.0) nAllLayers = nAllLayers + 1.0;

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
           hlt_pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60 == 1){
            hmet_pt_hltGrandOR->Fill(met_ptNoMu);
            nGrandOR = nGrandOR + 1.0;
            if(met_ptNoMu > 120.0 && muon_nLayers == 4) nhlt4Layers = nhlt4Layers + 1.0;
            if(met_ptNoMu > 120.0 && muon_nLayers == 5) nhlt5Layers = nhlt5Layers + 1.0;
            if(met_ptNoMu > 120.0 && muon_nLayers >= 6) nhlt6Layers = nhlt6Layers + 1.0;
            if(met_ptNoMu > 120.0) nGrandOR120 = nGrandOR120 + 1.0;
            if(met_ptNoMu > 120.0) nhltAllLayers = nhltAllLayers + 1.0;
        }

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
           hlt_pass_HLT_PFMET120_PFMHT120_IDTight_PFHT60 == 1){
            hmet_pt_hltGrandOR_exceptIsoTrk->Fill(met_ptNoMu);
            nGrandORExceptIsoTrk = nGrandORExceptIsoTrk + 1.0;
            if(met_ptNoMu > 120.0) nGrandORExceptIsoTrk120 = nGrandORExceptIsoTrk120 + 1.0;
            if(met_ptNoMu > 120.0 && muon_nLayers == 4) nhlt4LayersExceptIsoTrk = nhlt4LayersExceptIsoTrk + 1.0;
            if(met_ptNoMu > 120.0 && muon_nLayers == 5) nhlt5LayersExceptIsoTrk = nhlt5LayersExceptIsoTrk + 1.0;
            if(met_ptNoMu > 120.0 && muon_nLayers >= 6) nhlt6LayersExceptIsoTrk = nhlt6LayersExceptIsoTrk + 1.0;
            if(met_ptNoMu > 120.0) nhltAllLayersExceptIsoTrk = nhltAllLayersExceptIsoTrk + 1.0;
        }

        if(hlt_pass_HLT_MET105_IsoTrk50 == 1) {
            hmuon_pt_hlt->Fill(muon_pt);
            hmet_pt_hlt->Fill(met_ptNoMu);
            if(met_ptNoMu > 120.0 && muon_nLayers == 4) nhlt4LayersOnlyIsoTrk = nhlt4LayersOnlyIsoTrk + 1.0;
            if(met_ptNoMu > 120.0 && muon_nLayers == 5) nhlt5LayersOnlyIsoTrk = nhlt5LayersOnlyIsoTrk + 1.0;
            if(met_ptNoMu > 120.0 && muon_nLayers >= 6) nhlt6LayersOnlyIsoTrk = nhlt6LayersOnlyIsoTrk + 1.0;
            if(met_ptNoMu > 120.0) nhltAllLayersOnlyIsoTrk = nhltAllLayersOnlyIsoTrk + 1.0;
        }
        if(hlt_pass_hltMET105Filter == 1) {
            hmuon_pt_hltFilterMet->Fill(muon_pt);
            hmet_pt_hltFilterMet->Fill(met_ptNoMu);
            if(hlt_pass_hltTrk50Filter == 1) hmuon_pt_hltFilterTrk->Fill(muon_pt);
        }
    }

    std::cout << "Without 120: " << nGrandOR << " -- " << nGrandORExcept << " -- " << 100.0 - ((nGrandORExcept/nGrandOR)*100.0) <<  std::endl;
    std::cout << "With 120: " << nGrandOR120 << " -- " << nGrandORExcept120 << " -- " << 100.0 - ((nGrandORExcept120/nGrandOR120)*100.0) << std::endl;

    int marker = 20;

    TGraphAsymmErrors *h1 = new TGraphAsymmErrors("h1","MuonPt_HLT_MET105_IsoTrk50");
    h1->SetMarkerStyle(marker);
    h1->Divide(hmuon_pt_hlt,hmuon_pt);
    h1->SetTitle("Muon Pt - HLT_MET105_IsoTrk50");

    TGraphAsymmErrors *h2 = new TGraphAsymmErrors("h2","MuonPt_hltFilterMet");
    h2->SetMarkerStyle(marker);
    h2->Divide(hmuon_pt_hltFilterMet,hmuon_pt);
    h2->SetTitle("Muon Pt - hltMET105");

    TGraphAsymmErrors *h3 = new TGraphAsymmErrors("h3","MuonPt_hltFilterTrk");
    h3->SetMarkerStyle(marker);
    h3->Divide(hmuon_pt_hltFilterTrk,hmuon_pt_hltFilterMet);
    h3->SetTitle("Muon Pt - hltTrk50Filter");

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

    c1->cd();
    h1->Draw();
    c1->Update();
    c1->Write();
    h1->Write();

    c2->cd();
    h2->Draw();
    c2->Update();
    c2->Write();
    h2->Write();

    c3->cd();
    h3->Draw();
    c3->Update();
    c3->Write();
    h3->Write();

    c4->cd();
    h4->Draw();
    c4->Update();
    c4->Write();
    h4->Write();

    c5->cd();
    h5->Draw();
    c5->Update();
    c5->Write();
    h5->Write();

    c6->cd();
    h6->Draw();
    c6->Update();
    c6->Write();
    h6->Write();

    c7->cd();
    h7->Draw();
    c7->Update();
    c7->Write();
    h7->Write();

    c8->cd();
    h8->Draw();
    c8->Update();
    c8->Write();
    h8->Write();

    cesp->cd();
    auto overallEff = new TH1D("overallEff","",4,-0.5,3.5);
    overallEff->Fill(0.,nhlt4Layers/n4Layers);
    overallEff->Fill(1.,nhlt5Layers/n5Layers);
    overallEff->Fill(2.,nhlt6Layers/n6Layers);
    overallEff->Fill(3.,nhltAllLayers/nAllLayers);
    overallEff->GetXaxis()->SetBinLabel(1,"4Layers");
    overallEff->GetXaxis()->SetBinLabel(2,"5Layers");
    overallEff->GetXaxis()->SetBinLabel(3,"6PlusLayers");
    overallEff->GetXaxis()->SetBinLabel(4,"Total");
    overallEff->SetBinError(1,errorCalc(nhlt4Layers,n4Layers));
    overallEff->SetBinError(2,errorCalc(nhlt5Layers,n5Layers));
    overallEff->SetBinError(3,errorCalc(nhlt6Layers,n6Layers));
    overallEff->SetBinError(4,errorCalc(nhltAllLayers,nAllLayers));
    overallEff->GetYaxis()->SetRangeUser(0.0,1.0);
    overallEff->Draw();
    CMS_lumi(cesp, 0, 0);
    cesp->Write();
    overallEff->Write();

    cesp2->cd();
    auto overallEff2 = new TH1D("overallEff2","",4,-0.5,3.5);
    overallEff2->Fill(0.,nhlt4Layers/nhlt4LayersExceptIsoTrk);
    overallEff2->Fill(1.,nhlt5Layers/nhlt5LayersExceptIsoTrk);
    overallEff2->Fill(2.,nhlt6Layers/nhlt6LayersExceptIsoTrk);
    overallEff2->Fill(3.,nhltAllLayers/nhltAllLayersExceptIsoTrk);
    overallEff2->GetXaxis()->SetBinLabel(1,"4Layers");
    overallEff2->GetXaxis()->SetBinLabel(2,"5Layers");
    overallEff2->GetXaxis()->SetBinLabel(3,"6PlusLayers");
    overallEff2->GetXaxis()->SetBinLabel(4,"Total");
    overallEff2->SetBinError(1,errorCalc(nhlt4Layers,nhlt4LayersExceptIsoTrk));
    overallEff2->SetBinError(2,errorCalc(nhlt5Layers,nhlt5LayersExceptIsoTrk));
    overallEff2->SetBinError(3,errorCalc(nhlt6Layers,nhlt6LayersExceptIsoTrk));
    overallEff2->SetBinError(4,errorCalc(nhltAllLayers,nhltAllLayersExceptIsoTrk));
    overallEff2->GetYaxis()->SetRangeUser(0.8,1.2);
    overallEff2->Draw();
    CMS_lumi(cesp2, 0, 0);
    cesp2->Write();
    overallEff2->Write();

    cesp3->cd();
    auto overallEff3 = new TH1D("overallEff3","",4,-0.5,3.5);
    overallEff3->Fill(0.,nhlt4LayersOnlyIsoTrk/nhlt4Layers);
    overallEff3->Fill(1.,nhlt5LayersOnlyIsoTrk/nhlt5Layers);
    overallEff3->Fill(2.,nhlt6LayersOnlyIsoTrk/nhlt6Layers);
    overallEff3->Fill(3.,nhltAllLayersOnlyIsoTrk/nhltAllLayers);
    overallEff3->GetXaxis()->SetBinLabel(1,"4Layers");
    overallEff3->GetXaxis()->SetBinLabel(2,"5Layers");
    overallEff3->GetXaxis()->SetBinLabel(3,"6PlusLayers");
    overallEff3->GetXaxis()->SetBinLabel(4,"Total");
    overallEff3->SetBinError(1,errorCalc(nhlt4LayersOnlyIsoTrk,nhlt4Layers));
    overallEff3->SetBinError(2,errorCalc(nhlt5LayersOnlyIsoTrk,nhlt5Layers));
    overallEff3->SetBinError(3,errorCalc(nhlt6LayersOnlyIsoTrk,nhlt6Layers));
    overallEff3->SetBinError(4,errorCalc(nhltAllLayersOnlyIsoTrk,nhltAllLayers));
    overallEff3->GetYaxis()->SetRangeUser(0.0,1.0);
    overallEff3->Draw();
    CMS_lumi(cesp3, 0, 0);
    cesp3->Write();
    overallEff3->Write();

    out->Close();

}