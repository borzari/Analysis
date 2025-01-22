#include "tdrstyle.C"
#include "CMS_lumi.C"

void create_efficiencies_9files()
{
    setTDRStyle();

    gStyle->SetOptStat(0);

    writeExtraText = true;
    extraText = " Simulation";
    lumi_sqrtS = "13.6 TeV";

    std::string inputFile1 = "plots/AMSB_200_10_2023.root";
    std::string inputFile2 = "plots/AMSB_200_100_2023.root";
    std::string inputFile3 = "plots/AMSB_200_1000_2023.root";
    std::string inputFile4 = "plots/AMSB_700_10_2023.root";
    std::string inputFile5 = "plots/AMSB_700_100_2023.root";
    std::string inputFile6 = "plots/AMSB_700_1000_2023.root";
    std::string inputFile7 = "plots/AMSB_1200_10_2023.root";
    std::string inputFile8 = "plots/AMSB_1200_100_2023.root";
    std::string inputFile9 = "plots/AMSB_1200_1000_2023.root";
    TFile *f1 = TFile::Open(inputFile1.c_str());
    TFile *f2 = TFile::Open(inputFile2.c_str());
    TFile *f3 = TFile::Open(inputFile3.c_str());
    TFile *f4 = TFile::Open(inputFile4.c_str());
    TFile *f5 = TFile::Open(inputFile5.c_str());
    TFile *f6 = TFile::Open(inputFile6.c_str());
    TFile *f7 = TFile::Open(inputFile7.c_str());
    TFile *f8 = TFile::Open(inputFile8.c_str());
    TFile *f9 = TFile::Open(inputFile9.c_str());

    TH1 *ratioGrandOR1;
    TH1 *ratioGrandOR2;
    TH1 *ratioGrandOR3;
    TH1 *ratioGrandOR4;
    TH1 *ratioGrandOR5;
    TH1 *ratioGrandOR6;
    TH1 *ratioGrandOR7;
    TH1 *ratioGrandOR8;
    TH1 *ratioGrandOR9;
    f1->GetObject("overallEff2;1", ratioGrandOR1);
    f2->GetObject("overallEff2;1", ratioGrandOR2);
    f3->GetObject("overallEff2;1", ratioGrandOR3);
    f4->GetObject("overallEff2;1", ratioGrandOR4);
    f5->GetObject("overallEff2;1", ratioGrandOR5);
    f6->GetObject("overallEff2;1", ratioGrandOR6);
    f7->GetObject("overallEff2;1", ratioGrandOR7);
    f8->GetObject("overallEff2;1", ratioGrandOR8);
    f9->GetObject("overallEff2;1", ratioGrandOR9);

    TH2D *ratio = new TH2D("ratio","ratio",3,0,3,3,0,3);
    ratio->GetXaxis()->SetBinLabel(1,"200");
    ratio->GetXaxis()->SetBinLabel(2,"700");
    ratio->GetXaxis()->SetBinLabel(3,"1200");
    ratio->GetYaxis()->SetBinLabel(1,"10");
    ratio->GetYaxis()->SetBinLabel(2,"100");
    ratio->GetYaxis()->SetBinLabel(3,"1000");
    // ratio->SetBinContent(1,1,ratioGrandOR1->GetBinContent(4));
    // ratio->SetBinContent(1,2,ratioGrandOR2->GetBinContent(4));
    // ratio->SetBinContent(1,3,ratioGrandOR3->GetBinContent(4));
    // ratio->SetBinContent(2,1,ratioGrandOR4->GetBinContent(4));
    // ratio->SetBinContent(2,2,ratioGrandOR5->GetBinContent(4));
    // ratio->SetBinContent(2,3,ratioGrandOR6->GetBinContent(4));
    // ratio->SetBinContent(3,1,ratioGrandOR7->GetBinContent(4));
    // ratio->SetBinContent(3,2,ratioGrandOR8->GetBinContent(4));
    // ratio->SetBinContent(3,3,ratioGrandOR9->GetBinContent(4));
    ratio->SetBinContent(1,1,0.53);
    ratio->SetBinContent(1,2,6.46);
    ratio->SetBinContent(1,3,7.36);
    ratio->SetBinContent(2,1,0.23);
    ratio->SetBinContent(2,2,2.21);
    ratio->SetBinContent(2,3,2.92);
    ratio->SetBinContent(3,1,0.28);
    ratio->SetBinContent(3,2,2.28);
    ratio->SetBinContent(3,3,2.80);


    TCanvas *c1 = new TCanvas("c1", "c1", 800, 700);

    std::string outFile = "plots/2023_ratioGrandOR.root";
    TFile *file = new TFile(outFile.c_str(), "RECREATE");

    c1->cd();
    c1->SetRightMargin(0.15);
    c1->SetLeftMargin(0.126);
    ratio->GetXaxis()->SetTitle("Mass [GeV]");
    ratio->GetYaxis()->SetTitle("c#tau [cm]");
    // ratio->GetZaxis()->SetTitle("#frac{Regular PFMET paths + IsoTrk paths eff.}{Regular PFMET paths eff.}");
    ratio->GetZaxis()->SetTitle("Extra events when including IsoTrk paths [%]");
    ratio->GetXaxis()->SetTitleOffset(1.3);
    ratio->GetYaxis()->SetTitleOffset(1.7);
    // ratio->GetZaxis()->SetTitleOffset(2.05);
    ratio->GetZaxis()->SetTitleOffset(1.15);
    ratio->GetYaxis()->SetTitleSize(0.04);
    ratio->GetXaxis()->SetTitleSize(0.04);
    ratio->GetZaxis()->SetTitleSize(0.04);
    ratio->GetYaxis()->SetLabelSize(0.06);
    ratio->GetXaxis()->SetLabelSize(0.06);
    ratio->GetZaxis()->SetRangeUser(0,8);
    ratio->SetMarkerColor(1);
    ratio->SetLineColor(1);
    ratio->SetLineWidth(2);
    ratio->SetLineStyle(1);

    ratio->Draw("colz,text");
    CMS_lumi(c1, 0, 0);
    c1->Update();
    c1->Write();
    std::string c1Name = "plots/2023_ratio_GrandOR.pdf";
    c1->Print(c1Name.c_str());

    file->Close();
}