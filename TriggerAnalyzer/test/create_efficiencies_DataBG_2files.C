#include "tdrstyle.C"
#include "CMS_lumi.C"

void create_efficiencies_DataBG_2files()
{
    setTDRStyle();

    gStyle->SetOptStat(0);

    writeExtraText = true;
    extraText = " Preliminary";
    lumi_sqrtS = "27.72 fb^{-1} (13.6 TeV)";
    // lumi_sqrtS = "(13.6 TeV)";

    std::string inputFile1 = "plots/log_right_Muon_2023.root";
    std::string inputFile2 = "plots/log_right_WToLNu_2023.root";
    TFile *f1 = TFile::Open(inputFile1.c_str());
    TFile *f2 = TFile::Open(inputFile2.c_str());

    TGraph *MET105_file1;
    TGraph *MET105_file2;
    f1->GetObject("Graph;1", MET105_file1);
    f2->GetObject("Graph;1", MET105_file2);

    TGraph *filterMET105_file1;
    TGraph *filterMET105_file2;
    f1->GetObject("Graph;2", filterMET105_file1);
    f2->GetObject("Graph;2", filterMET105_file2);

    TGraph *GrandOR_file1;
    TGraph *GrandOR_file2;
    f1->GetObject("Graph;3", GrandOR_file1);
    f2->GetObject("Graph;3", GrandOR_file2);

    TGraph *filterTrk50_file1;
    TGraph *filterTrk50_file2;
    f1->GetObject("Graph;6", filterTrk50_file1);
    f2->GetObject("Graph;6", filterTrk50_file2);

    float xmin = 0.0;
    float xmax = 1000.0;
    float ymin = 0.0;
    float ymax = 1.2;
    float x0 = 0.1891;
    float x1 = 0.4395;
    float y0 = 0.8292;
    float y1 = 0.9214;
    std::string legend_file1 = "2023 data";
    std::string legend_file2 = "W #rightarrow l#nu MC";

    TLine *l = new TLine(xmin, 1.0, xmax, 1.0);
    l->SetLineStyle(7);
    l->SetLineWidth(2);

    TCanvas *c1 = new TCanvas("c1", "c1", 700, 700);
    TCanvas *c2 = new TCanvas("c2", "c2", 700, 700);
    TCanvas *c3 = new TCanvas("c3", "c3", 700, 700);
    TCanvas *c4 = new TCanvas("c4", "c4", 700, 700);

    TCanvas *c5 = new TCanvas("c5", "c5", 700, 700);
    TCanvas *c6 = new TCanvas("c6", "c6", 700, 700);
    TCanvas *c7 = new TCanvas("c7", "c7", 700, 700);
    TCanvas *c8 = new TCanvas("c8", "c8", 700, 700);

    std::string outFile = "DatavsBG_hlt_effs_log_2023.root";
    TFile *file = new TFile(outFile.c_str(), "RECREATE");

    c1->cd();
    c1->SetLogx();
    MET105_file1->GetYaxis()->SetTitle("Trigger Efficiency");
    MET105_file1->GetXaxis()->SetTitle("PF E_{T}^{miss, no mu} [GeV]");
    MET105_file1->GetXaxis()->SetTitleOffset(1.1);
    MET105_file1->GetYaxis()->SetTitleSize(0.04);
    MET105_file1->GetXaxis()->SetTitleSize(0.04);
    MET105_file1->GetYaxis()->SetLabelSize(0.04);
    MET105_file1->GetXaxis()->SetLabelSize(0.04);
    MET105_file1->SetMinimum(ymin);
    MET105_file1->SetMaximum(ymax);
    MET105_file1->SetMarkerColor(1);
    MET105_file1->SetLineColor(1);
    MET105_file1->SetMarkerStyle(20);
    MET105_file1->GetXaxis()->SetLimits(xmin, xmax);

    MET105_file2->SetMarkerColor(4);
    MET105_file2->SetLineColor(4);
    MET105_file2->SetMarkerStyle(20);
    
    auto legend1 = new TLegend(x0, y0, x1, y1);
    legend1->AddEntry(MET105_file1, legend_file1.c_str(), "P");
    legend1->AddEntry(MET105_file2, legend_file2.c_str(), "P");
    legend1->SetLineWidth(0);
    legend1->SetTextSize(0.03559);
    legend1->SetBorderSize(0);
    legend1->SetTextFont(42);
    // legend1->SetHeader("HLT_MET105_IsoTrk50_v*");
    MET105_file1->Draw("AP");
    MET105_file2->Draw("P,SAME");
    legend1->Draw();
    l->Draw();
    CMS_lumi(c1, 0, 0);
    c1->Update();
    c1->Write();
    std::string c1Name = "plots/AN_24_155/2023vsMC_HLT_MET105_IsoTrk50_v.pdf";
    c1->Print(c1Name.c_str());

    c2->cd();
    c2->SetLogx();
    filterMET105_file1->GetYaxis()->SetTitle("Efficiency");
    filterMET105_file1->GetXaxis()->SetTitle("PF E_{T}^{miss, no mu} [GeV]");
    filterMET105_file1->GetXaxis()->SetTitleOffset(1.1);
    filterMET105_file1->GetYaxis()->SetTitleSize(0.04);
    filterMET105_file1->GetXaxis()->SetTitleSize(0.04);
    filterMET105_file1->GetYaxis()->SetLabelSize(0.04);
    filterMET105_file1->GetXaxis()->SetLabelSize(0.04);
    filterMET105_file1->SetMinimum(ymin);
    filterMET105_file1->SetMaximum(ymax);
    filterMET105_file1->SetMarkerColor(1);
    filterMET105_file1->SetLineColor(1);
    filterMET105_file1->SetMarkerStyle(20);
    filterMET105_file1->GetXaxis()->SetLimits(xmin, xmax);

    filterMET105_file2->SetMarkerColor(4);
    filterMET105_file2->SetLineColor(4);
    filterMET105_file2->SetMarkerStyle(20);
    
    // legend1->SetHeader("hltMET105");
    filterMET105_file1->Draw("AP");
    filterMET105_file2->Draw("P,SAME");
    legend1->Draw();
    l->Draw();
    CMS_lumi(c2, 0, 0);
    c2->Update();
    c2->Write();
    std::string c2Name = "plots/AN_24_155/2023vsMC_filterMET105.pdf";
    c2->Print(c2Name.c_str());

    c3->cd();
    c3->SetLogx();
    GrandOR_file1->GetYaxis()->SetTitle("Efficiency");
    GrandOR_file1->GetXaxis()->SetTitle("PF E_{T}^{miss, no mu} [GeV]");
    GrandOR_file1->GetXaxis()->SetTitleOffset(1.1);
    GrandOR_file1->GetYaxis()->SetTitleSize(0.04);
    GrandOR_file1->GetXaxis()->SetTitleSize(0.04);
    GrandOR_file1->GetYaxis()->SetLabelSize(0.04);
    GrandOR_file1->GetXaxis()->SetLabelSize(0.04);
    GrandOR_file1->SetMinimum(ymin);
    GrandOR_file1->SetMaximum(ymax);
    GrandOR_file1->SetMarkerColor(1);
    GrandOR_file1->SetLineColor(1);
    GrandOR_file1->SetMarkerStyle(20);
    GrandOR_file1->GetXaxis()->SetLimits(xmin, xmax);

    GrandOR_file2->SetMarkerColor(4);
    GrandOR_file2->SetLineColor(4);
    GrandOR_file2->SetMarkerStyle(20);
    
    legend1->SetHeader("OR of Signal Paths");
    GrandOR_file1->Draw("AP");
    GrandOR_file2->Draw("P,SAME");
    legend1->Draw();
    l->Draw();
    CMS_lumi(c3, 0, 0);
    c3->Update();
    c3->Write();
    std::string c3Name = "plots/AN_24_155/2023vsMC_GrandOR.pdf";
    c3->Print(c3Name.c_str());

    legend1->SetHeader("");

    c4->cd();
    c4->SetLogx();
    filterTrk50_file1->GetYaxis()->SetTitle("Efficiency");
    filterTrk50_file1->GetXaxis()->SetTitle("Muon p_{T} [GeV]");
    filterTrk50_file1->GetXaxis()->SetTitleOffset(1.1);
    filterTrk50_file1->GetYaxis()->SetTitleSize(0.04);
    filterTrk50_file1->GetXaxis()->SetTitleSize(0.04);
    filterTrk50_file1->GetYaxis()->SetLabelSize(0.04);
    filterTrk50_file1->GetXaxis()->SetLabelSize(0.04);
    filterTrk50_file1->SetMinimum(ymin);
    filterTrk50_file1->SetMaximum(ymax);
    filterTrk50_file1->SetMarkerColor(1);
    filterTrk50_file1->SetLineColor(1);
    filterTrk50_file1->SetMarkerStyle(20);
    filterTrk50_file1->GetXaxis()->SetLimits(xmin, xmax);

    filterTrk50_file2->SetMarkerColor(4);
    filterTrk50_file2->SetLineColor(4);
    filterTrk50_file2->SetMarkerStyle(20);
    
    // legend1->SetHeader("HLT_MET105_IsoTrk50_v*, hltMET105 applied");
    filterTrk50_file1->Draw("AP");
    filterTrk50_file2->Draw("P,SAME");
    legend1->Draw();
    l->Draw();
    CMS_lumi(c4, 0, 0);
    c4->Update();
    c4->Write();
    std::string c4Name = "plots/AN_24_155/2023vsMC_filterIsoTrk50.pdf";
    c4->Print(c4Name.c_str());

    c5->cd();
    c5->SetLogx();
    filterMET105_file1->GetYaxis()->SetTitle("Efficiency");
    filterMET105_file1->GetXaxis()->SetTitle("PF E_{T}^{miss, no mu} [GeV]");
    filterMET105_file1->GetXaxis()->SetTitleOffset(1.1);
    filterMET105_file1->GetYaxis()->SetTitleSize(0.04);
    filterMET105_file1->GetXaxis()->SetTitleSize(0.04);
    filterMET105_file1->GetYaxis()->SetLabelSize(0.04);
    filterMET105_file1->GetXaxis()->SetLabelSize(0.04);
    filterMET105_file1->SetMinimum(ymin);
    filterMET105_file1->SetMaximum(ymax);
    filterMET105_file1->SetMarkerColor(1);
    filterMET105_file1->SetLineColor(1);
    filterMET105_file1->SetMarkerStyle(20);
    filterMET105_file1->GetXaxis()->SetLimits(xmin, xmax);

    auto legend2 = new TLegend(x0, y0, x1, y1);
    legend2->AddEntry(MET105_file1, legend_file1.c_str(), "P");
    legend2->SetLineWidth(0);
    legend2->SetTextSize(0.03559);
    legend2->SetBorderSize(0);
    legend2->SetTextFont(42);
    // legend2->SetHeader("hltMET105");
    filterMET105_file1->Draw("AP");
    legend2->Draw();
    l->Draw();
    CMS_lumi(c5, 0, 0);
    c5->Update();
    c5->Write();
    std::string c5Name = "plots/AN_24_155/2022_filterMET105.pdf";
    c5->Print(c5Name.c_str());

    c6->cd();
    c6->SetLogx();
    filterMET105_file2->GetYaxis()->SetTitle("Efficiency");
    filterMET105_file2->GetXaxis()->SetTitle("PF E_{T}^{miss, no mu} [GeV]");
    filterMET105_file2->GetXaxis()->SetTitleOffset(1.1);
    filterMET105_file2->GetYaxis()->SetTitleSize(0.04);
    filterMET105_file2->GetXaxis()->SetTitleSize(0.04);
    filterMET105_file2->GetYaxis()->SetLabelSize(0.04);
    filterMET105_file2->GetXaxis()->SetLabelSize(0.04);
    filterMET105_file2->SetMinimum(ymin);
    filterMET105_file2->SetMaximum(ymax);
    filterMET105_file2->SetMarkerColor(1);
    filterMET105_file2->SetLineColor(1);
    filterMET105_file2->SetMarkerStyle(20);
    filterMET105_file2->GetXaxis()->SetLimits(xmin, xmax);

    auto legend3 = new TLegend(x0, y0, x1, y1);
    legend3->AddEntry(MET105_file1, legend_file2.c_str(), "P");
    legend3->SetLineWidth(0);
    legend3->SetTextSize(0.03559);
    legend3->SetBorderSize(0);
    legend3->SetTextFont(42);
    // legend3->SetHeader("hltMET105");
    filterMET105_file2->Draw("AP");
    legend3->Draw();
    l->Draw();
    CMS_lumi(c6, 0, 0);
    c6->Update();
    c6->Write();
    std::string c6Name = "plots/AN_24_155/2023_filterMET105.pdf";
    c6->Print(c6Name.c_str());

    c7->cd();
    c7->SetLogx();
    filterTrk50_file1->GetYaxis()->SetTitle("Efficiency");
    filterTrk50_file1->GetXaxis()->SetTitle("Muon p_{T} [GeV]");
    filterTrk50_file1->GetXaxis()->SetTitleOffset(1.3);
    filterTrk50_file1->GetYaxis()->SetTitleSize(0.04);
    filterTrk50_file1->GetXaxis()->SetTitleSize(0.04);
    filterTrk50_file1->GetYaxis()->SetLabelSize(0.04);
    filterTrk50_file1->GetXaxis()->SetLabelSize(0.04);
    filterTrk50_file1->SetMinimum(ymin);
    filterTrk50_file1->SetMaximum(ymax);
    filterTrk50_file1->SetMarkerColor(1);
    filterTrk50_file1->SetLineColor(1);
    filterTrk50_file1->SetMarkerStyle(20);
    filterTrk50_file1->GetXaxis()->SetLimits(xmin, xmax);

    // legend2->SetHeader("HLT_MET105_IsoTrk50_v*, hltMET105 applied");
    filterTrk50_file1->Draw("AP");
    legend2->Draw();
    l->Draw();
    CMS_lumi(c7, 0, 0);
    c7->Update();
    c7->Write();
    std::string c7Name = "plots/AN_24_155/2022_filterIsoTrk50.pdf";
    c7->Print(c7Name.c_str());

    c8->cd();
    c8->SetLogx();
    filterTrk50_file2->GetYaxis()->SetTitle("Efficiency");
    filterTrk50_file2->GetXaxis()->SetTitle("Muon p_{T} [GeV]");
    filterTrk50_file2->GetXaxis()->SetTitleOffset(1.3);
    filterTrk50_file2->GetYaxis()->SetTitleSize(0.04);
    filterTrk50_file2->GetXaxis()->SetTitleSize(0.04);
    filterTrk50_file2->GetYaxis()->SetLabelSize(0.04);
    filterTrk50_file2->GetXaxis()->SetLabelSize(0.04);
    filterTrk50_file2->SetMinimum(ymin);
    filterTrk50_file2->SetMaximum(ymax);
    filterTrk50_file2->SetMarkerColor(1);
    filterTrk50_file2->SetLineColor(1);
    filterTrk50_file2->SetMarkerStyle(20);
    filterTrk50_file2->GetXaxis()->SetLimits(xmin, xmax);

    // legend3->SetHeader("HLT_MET105_IsoTrk50_v*, hltMET105 applied");
    filterTrk50_file2->Draw("AP");
    legend3->Draw();
    l->Draw();
    CMS_lumi(c8, 0, 0);
    c8->Update();
    c8->Write();
    std::string c8Name = "plots/AN_24_155/2023_filterIsoTrk50.pdf";
    c8->Print(c8Name.c_str());

    file->Close();
}