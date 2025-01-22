#include "tdrstyle.C"
#include "CMS_lumi.C"

void create_efficiencies_3files()
{
    setTDRStyle();

    gStyle->SetOptStat(0);

    writeExtraText = true;
    extraText = " Simulation";
    lumi_sqrtS = "13.6 TeV";

    std::string mass = "900";

    // std::string inputFile1 = "Plots_TriggerEffs_10cm/10_hlt_efficiencies.root";
    // std::string inputFile2 = "Plots_TriggerEffs_100cm/100_hlt_efficiencies.root";
    // std::string inputFile3 = "Plots_TriggerEffs_1000cm/1000_hlt_efficiencies.root";
    std::string inputFile1 = "plots/AMSB_" + mass + "_10_2023.root";
    std::string inputFile2 = "plots/AMSB_" + mass + "_100_2023.root";
    std::string inputFile3 = "plots/AMSB_" + mass + "_1000_2023.root";
    // std::string inputFile1 = "plots/AMSB_200_100_2023.root";
    // std::string inputFile2 = "plots/AMSB_700_100_2023.root";
    // std::string inputFile3 = "plots/AMSB_1200_100_2023.root";
    TFile *f1 = TFile::Open(inputFile1.c_str());
    TFile *f2 = TFile::Open(inputFile2.c_str());
    TFile *f3 = TFile::Open(inputFile3.c_str());

    TH1 *PFMET110_file1;
    TH1 *PFMET110_file2;
    TH1 *PFMET110_file3;
    f1->GetObject("overallEff;1", PFMET110_file1);
    f2->GetObject("overallEff;1", PFMET110_file2);
    f3->GetObject("overallEff;1", PFMET110_file3);

    TH1 *PFMET110_file12;
    TH1 *PFMET110_file22;
    TH1 *PFMET110_file32;
    f1->GetObject("overallEff2;1", PFMET110_file12);
    f2->GetObject("overallEff2;1", PFMET110_file22);
    f3->GetObject("overallEff2;1", PFMET110_file32);

    TH1 *PFMET110_file13;
    TH1 *PFMET110_file23;
    TH1 *PFMET110_file33;
    f1->GetObject("overallEff3;1", PFMET110_file13);
    f2->GetObject("overallEff3;1", PFMET110_file23);
    f3->GetObject("overallEff3;1", PFMET110_file33);
    
    float ymin = 0.0;
    float ymax = 1.0;
    float x0 = 0.1891;
    float x1 = 0.5395;
    float y0 = 0.7742;
    float y1 = 0.9164;
    // std::string legend_file1 = "AMSB - " + mass + " GeV, 10cm";
    // std::string legend_file2 = "AMSB - " + mass + " GeV, 100cm";
    // std::string legend_file3 = "AMSB - " + mass + " GeV, 1000cm";
    std::string legend_file1 = "m_{#tilde{#chi}_{1}^{ #pm}} = 900 GeV, c#tau = 10cm";
    std::string legend_file2 = "m_{#tilde{#chi}_{1}^{ #pm}} = 900 GeV, c#tau = 100cm";
    std::string legend_file3 = "m_{#tilde{#chi}_{1}^{ #pm}} = 900 GeV, c#tau = 1000cm";

    TCanvas *c1 = new TCanvas("c1", "c1", 700, 700);
    TCanvas *c2 = new TCanvas("c2", "c2", 700, 700);
    TCanvas *c3 = new TCanvas("c3", "c3", 700, 700);

    std::string outFile = "plots/2023_" + mass + "_10vs100vs1000_hlt_effs_comp.root";
    TFile *file = new TFile(outFile.c_str(), "RECREATE");

    c1->cd();
    PFMET110_file1->GetYaxis()->SetTitle("Trigger Efficiency");
    // PFMET110_file1->GetXaxis()->SetTitle("100");
    PFMET110_file1->GetYaxis()->SetTitleSize(0.04);
    PFMET110_file1->GetXaxis()->SetTitleSize(0.04);
    PFMET110_file1->GetYaxis()->SetLabelSize(0.04);
    PFMET110_file1->GetXaxis()->SetLabelSize(0.04);
    PFMET110_file1->SetMinimum(ymin);
    PFMET110_file1->SetMaximum(ymax);
    PFMET110_file1->SetMarkerColor(1);
    PFMET110_file1->SetLineColor(1);
    PFMET110_file1->SetLineWidth(2);
    PFMET110_file1->SetLineStyle(1);
    PFMET110_file1->GetXaxis()->SetBinLabel(1,"4 layers");
    PFMET110_file1->GetXaxis()->SetBinLabel(2,"5 layers");
    PFMET110_file1->GetXaxis()->SetBinLabel(3,"6+ layers");
    PFMET110_file1->GetXaxis()->SetBinLabel(4,"Total");

    PFMET110_file2->SetMarkerColor(4);
    PFMET110_file2->SetLineColor(4);
    PFMET110_file2->SetLineWidth(2);
    PFMET110_file2->SetLineStyle(2);
    
    PFMET110_file3->SetMarkerColor(2);
    PFMET110_file3->SetLineColor(2);
    PFMET110_file3->SetLineWidth(2);
    PFMET110_file3->SetLineStyle(3);

    auto legend1 = new TLegend(x0, y0, x1, y1);
    legend1->AddEntry(PFMET110_file1, legend_file1.c_str(), "PL");
    legend1->AddEntry(PFMET110_file2, legend_file2.c_str(), "PL");
    legend1->AddEntry(PFMET110_file3, legend_file3.c_str(), "PL");
    legend1->SetLineWidth(0);
    legend1->SetTextSize(0.02559);
    legend1->SetBorderSize(0);
    legend1->SetTextFont(42);
    legend1->SetHeader("GrandOR");
    PFMET110_file1->Draw();
    PFMET110_file2->Draw("SAME");
    PFMET110_file3->Draw("SAME");
    legend1->Draw();
    CMS_lumi(c1, 0, 0);
    c1->Update();
    c1->Write();
    std::string c1Name = "plots/2023_" + mass + "_10vs100vs1000_HLT_GrandOR.pdf";
    c1->Print(c1Name.c_str());

    c2->cd();
    PFMET110_file12->GetYaxis()->SetTitle("Trigger Efficiency");
    // PFMET110_file1->GetXaxis()->SetTitle("100");
    PFMET110_file12->GetYaxis()->SetTitleSize(0.04);
    PFMET110_file12->GetXaxis()->SetTitleSize(0.04);
    PFMET110_file12->GetYaxis()->SetLabelSize(0.04);
    PFMET110_file12->GetXaxis()->SetLabelSize(0.04);
    PFMET110_file12->SetMinimum(0.8);
    PFMET110_file12->SetMaximum(1.2);
    PFMET110_file12->SetMarkerColor(1);
    PFMET110_file12->SetLineColor(1);
    PFMET110_file12->SetLineWidth(2);
    PFMET110_file12->SetLineStyle(1);

    PFMET110_file22->SetMarkerColor(4);
    PFMET110_file22->SetLineColor(4);
    PFMET110_file22->SetLineWidth(2);
    PFMET110_file22->SetLineStyle(2);
    
    PFMET110_file32->SetMarkerColor(2);
    PFMET110_file32->SetLineColor(2);
    PFMET110_file32->SetLineWidth(2);
    PFMET110_file32->SetLineStyle(3);

    auto legend2 = new TLegend(x0, y0, x1, y1);
    legend2->AddEntry(PFMET110_file12, legend_file1.c_str(), "PL");
    legend2->AddEntry(PFMET110_file22, legend_file2.c_str(), "PL");
    legend2->AddEntry(PFMET110_file32, legend_file3.c_str(), "PL");
    legend2->SetLineWidth(0);
    legend2->SetTextSize(0.02559);
    legend2->SetBorderSize(0);
    legend2->SetTextFont(42);
    legend2->SetHeader("Diff MET+IsoTrk");
    PFMET110_file12->Draw();
    PFMET110_file22->Draw("SAME");
    PFMET110_file32->Draw("SAME");
    legend2->Draw();
    CMS_lumi(c2, 0, 0);
    c2->Update();
    c2->Write();
    std::string c2Name = "plots/2023_" + mass + "_10vs100vs1000_HLT_diffMETIsoTrk.pdf";
    c2->Print(c2Name.c_str());

    c3->cd();
    PFMET110_file13->GetYaxis()->SetTitle("Trigger Efficiency");
    PFMET110_file13->GetXaxis()->SetTitle("100");
    PFMET110_file13->GetYaxis()->SetTitleSize(0.04);
    PFMET110_file13->GetXaxis()->SetTitleSize(0.04);
    PFMET110_file13->GetYaxis()->SetLabelSize(0.04);
    PFMET110_file13->GetXaxis()->SetLabelSize(0.06);
    PFMET110_file13->SetMinimum(0.0);
    PFMET110_file13->SetMaximum(1.0);
    PFMET110_file13->SetMarkerColor(1);
    PFMET110_file13->SetMarkerStyle(20);
    PFMET110_file13->SetLineColor(1);
    PFMET110_file13->SetLineWidth(2);
    PFMET110_file13->SetLineStyle(1);
    PFMET110_file13->GetXaxis()->SetBinLabel(1,"4 layers");
    PFMET110_file13->GetXaxis()->SetBinLabel(2,"5 layers");
    PFMET110_file13->GetXaxis()->SetBinLabel(3,"6+ layers");
    PFMET110_file13->GetXaxis()->SetBinLabel(4,"Total");

    PFMET110_file23->SetMarkerColor(4);
    PFMET110_file23->SetMarkerStyle(20);
    PFMET110_file23->SetLineColor(4);
    PFMET110_file23->SetLineWidth(2);
    PFMET110_file23->SetLineStyle(2);
    
    PFMET110_file33->SetMarkerColor(2);
    PFMET110_file33->SetMarkerStyle(20);
    PFMET110_file33->SetLineColor(2);
    PFMET110_file33->SetLineWidth(2);
    PFMET110_file33->SetLineStyle(3);

    auto legend3 = new TLegend(x0, y0, x1, y1);
    legend3->AddEntry(PFMET110_file13, legend_file1.c_str(), "PL");
    legend3->AddEntry(PFMET110_file23, legend_file2.c_str(), "PL");
    legend3->AddEntry(PFMET110_file33, legend_file3.c_str(), "PL");
    legend3->SetLineWidth(0);
    legend3->SetTextSize(0.03559);
    legend3->SetBorderSize(0);
    legend3->SetTextFont(42);
    legend3->SetHeader("#tilde{#chi}_{1}^{ #pm} #rightarrow #tilde{#chi}_{1}^{ 0}+ X");
    PFMET110_file13->Draw();
    PFMET110_file23->Draw("SAME");
    PFMET110_file33->Draw("SAME");
    legend3->Draw();
    CMS_lumi(c3, 0, 0);
    c3->Update();
    c3->Write();
    std::string c3Name = "plots/EXOLLPTriggerPaper/" + mass + "_10vs100vs1000_HLT_MET105_IsoTrk50.pdf";
    c3->Print(c3Name.c_str());

    file->Close();
}