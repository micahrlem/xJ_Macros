#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

#include "xj_functions.h"
#include "read_binning.h"


void Draw_Herwig_QA(string infile01 = "Hists_xJ_ProjectionTest_TREE_DIJET_v6_1_ana462_2024p010_v001_gl10-00047352-00047733.root" )

{
SetsPhenixStyle();
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

TFile* infile3 = TFile::Open(infile01.c_str());//Unmatched Data Hists 30 zvtx

//Unmatched Data Hists 60 zvtx
TH1D* hUnMatchedDataPtLead_60vtx = (TH1D*)infile3->Get("hDataPtLead_60vtx");TH1D* hUnMatchedDataPtSubLead_60vtx = (TH1D*)infile3->Get("hDataPtSubLead_60vtx");
TH2D* hUnMatchedDataEtaPhiLead_60vtx = (TH2D*)infile3->Get("hDataEtaPhiLead_60vtx");TH2D* hUnMatchedDataEtaPhiSubLead_60vtx = (TH2D*)infile3->Get("hDataEtaPhiSubLead_60vtx");
 TH1D* h_pt1pt2 = (TH1D*)infile3->Get("h_pt1pt2");TH1D* h_xj_classical = (TH1D*)infile3->Get("h_xj_classical");TH1D* h_xj_projected = (TH1D*)infile3->Get("h_xj_projected");


TCanvas* can_pt1pt2 = new TCanvas("can_pt1pt2", "p_{T1} vs p_{T2}", 800, 800);
gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
SetHist(h_pt1pt2, "p_{T1} [GeV]", "p_{T2} [GeV]", 2, 20, 1, 1);
h_pt1pt2->Draw("Colz");
h_pt1pt2->SetAxisRange(0, 200, "X"); // Adjust X range if needed
h_pt1pt2->SetAxisRange(0, 200, "Y"); // Adjust Y range if needed
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
h_pt1pt2->GetZaxis()->SetLabelSize(0.035);

// Legend
TLegend *sleg_pt1pt2 = new TLegend(.14, .60, .44, .88);
sleg_pt1pt2->SetFillStyle(0);
sleg_pt1pt2->SetBorderSize(0);
sleg_pt1pt2->SetTextSize(0.030);
sleg_pt1pt2->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
sleg_pt1pt2->AddEntry("", "p+p #sqrt{s}=200 GeV", "");
sleg_pt1pt2->AddEntry("", "anti-#it{k}_{#it{t}} #it{R} = 0.4", "");
sleg_pt1pt2->AddEntry("", "Leading Jets", "");
sleg_pt1pt2->AddEntry("", "p_{T,1} > 15 GeV , p_{T,2} > 8 GeV", "");
sleg_pt1pt2->AddEntry("", "|z_{vtx}| < 60 cm", "");
sleg_pt1pt2->Draw();



// Plotting Xj Classical vs. Projected
TCanvas* can_xj_comparison = new TCanvas("can_xj_comparison", "Xj Classical vs Projected", 700, 800);
can_xj_comparison->cd(1);
gPad->SetLogy(0); // Remove log scale

// Normalize histograms using the function with 19 bins
 xj_functions::normalize_histo(h_xj_classical, 19);
 xj_functions::normalize_histo(h_xj_projected, 19);

 xj_functions::plotxj_and_ratio(can_xj_comparison, h_xj_classical, h_xj_projected, h_xj_classical,h_xj_classical,false, "x_{J, Data}","Projected/Classical",0, 1.0, 2, 20, 4, 53,  2,  20,  1, 20);
/*
// Set histogram properties
SetHist(h_xj_classical, "x_{J}", "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", 2, 20, 3, 1.5);
 h_xj_classical->SetAxisRange(0.1, 0.99, "X");// Set X range from 0 to 0.99
 h_xj_classical->SetAxisRange(0, 3, "Y");h_xj_classical->GetYaxis()->SetTitleSize(0.04);
h_xj_classical->Draw("P");

SetHist(h_xj_projected, "x_{J}", "Normalized Counts", 4, 24, 3, 1.5);
h_xj_projected->Draw("P SAME");*/

// Legend for Classical vs Projected
TLegend *leg_xj = new TLegend(0.65, 0.45, 0.87, 0.55);
leg_xj->SetBorderSize(1);
leg_xj->SetTextSize(0.026);  // Set text size to 0.026
leg_xj->AddEntry(h_xj_classical, "Xj Classical", "P");
leg_xj->AddEntry(h_xj_projected, "Xj Projected", "P");
leg_xj->Draw();

// General Legend
TLegend *sleg_xj = new TLegend(0.45, 0.68, 0.75, 0.93);
sleg_xj->SetFillStyle(0);
sleg_xj->SetBorderSize(0);
sleg_xj->SetTextSize(0.030);  // Set text size to 0.030
sleg_xj->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
sleg_xj->AddEntry("", "p+p #sqrt{s}=200 GeV", "");
sleg_xj->AddEntry("", "anti-#it{k}_{#it{t}} #it{R} = 0.4", "");
sleg_xj->AddEntry("", "All Dijet Pairs", "");
sleg_xj->AddEntry("", "p_{T,1} > 18.28 GeV , p_{T,2} > 5.5 GeV", "");
sleg_xj->AddEntry("", "|z_{vtx}| < 60 cm, |#eta| < 0.7", "");
sleg_xj->AddEntry("", "|#Delta #phi| > 3#pi/4", "");
sleg_xj->Draw();


//Plotting Pt Together
TCanvas* can_lead_sub_pt = new TCanvas("can_lead_sub_pt","",800,800);
can_lead_sub_pt->cd(1);gPad->SetLogy();
//Data UnMatched
SetHist(hUnMatchedDataPtLead_60vtx,"p_{T} [Gev]", "Counts",2,20,3,1.5);
hUnMatchedDataPtLead_60vtx->Draw("P");
hUnMatchedDataPtLead_60vtx->SetAxisRange(10,54); 
SetHist(hUnMatchedDataPtSubLead_60vtx,"p_{T} [Gev]", "Counts",2,24,3,1.5);
 hUnMatchedDataPtSubLead_60vtx->Draw("same, P");
 TLegend *leg1 = new TLegend (0.65, 0.45, 0.87,0.55);
 leg1->SetTextSize(0.026);
 leg1->AddEntry(hUnMatchedDataPtLead_60vtx, "Leading Jets", "P");
 leg1->AddEntry(hUnMatchedDataPtSubLead_60vtx, "Subleading Jets", "P");
 leg1->SetBorderSize(1);
 leg1->Draw();

 TLegend *sleglead60 = new TLegend(.45,.68,.75,.93);
  sleglead60->SetFillStyle(0);
  sleglead60->SetBorderSize(0);
  sleglead60->SetTextSize(0.030);
  sleglead60->AddEntry("","#it{#bf{sPHENIX}} Internal","");
  sleglead60->AddEntry("","p+p #sqrt{s}=200 GeV ","");
  sleglead60->AddEntry("","anti-#it{k}_{#it{t}} #it{R} = 0.4","");
  sleglead60->AddEntry("","All Dijet Pairs","");
  sleglead60->AddEntry("","p_{T,1} > 15 GeV , p_{T,2} > 8 GeV","");
  sleglead60->AddEntry("","|z_{vtx}| < 60 cm, |#eta| < 0.7","");
  sleglead60->AddEntry("","|#Delta #phi| > 3#pi/4","");
  sleglead60->Draw();


TCanvas* can_lead_sub_etaphi = new TCanvas("can_lead_sub_etaphi","",1600,800);
 can_lead_sub_etaphi->Divide(2,1);can_lead_sub_etaphi->cd(1);gPad->SetRightMargin(0.16);
//Data UnMatched
SetHist(hUnMatchedDataEtaPhiLead_60vtx,"#phi", "#eta",2,20,1,1);
 hUnMatchedDataEtaPhiLead_60vtx->Draw("Colz");hUnMatchedDataEtaPhiLead_60vtx->SetAxisRange(-0.7,0.69,"Y"); hUnMatchedDataEtaPhiLead_60vtx->SetAxisRange(0,375,"Z"); 
 hUnMatchedDataEtaPhiLead_60vtx->GetZaxis()->SetLabelSize(0.035);
TLegend *slegleadeta60 = new TLegend(.14,.60,.44,.88);
  slegleadeta60->SetFillStyle(0);
  slegleadeta60->SetBorderSize(0);
  slegleadeta60->SetTextSize(0.034);
  slegleadeta60->AddEntry("","#it{#bf{sPHENIX}} Internal","");
  slegleadeta60->AddEntry("","p+p #sqrt{s}=200 GeV ","");
  slegleadeta60->AddEntry("","anti-#it{k}_{#it{t}} #it{R} = 0.4","");
  slegleadeta60->AddEntry("","Leading Jets","");
  slegleadeta60->AddEntry("","p_{T,1} > 15 GeV , p_{T,2} > 8 GeV","");
  slegleadeta60->AddEntry("","|z_{vtx}| < 60 cm, |#eta| < 0.7","");
  slegleadeta60->AddEntry("","|#Delta #phi| > 3#pi/4","");
  slegleadeta60->Draw();


 can_lead_sub_etaphi->cd(2);gPad->SetRightMargin(0.14);
 //Data UnMatched
SetHist(hUnMatchedDataEtaPhiSubLead_60vtx,"#phi", "#eta",1,27,3,1.5);
hUnMatchedDataEtaPhiSubLead_60vtx->Draw("Colz");;hUnMatchedDataEtaPhiSubLead_60vtx->SetAxisRange(-0.7,0.69,"Y"); hUnMatchedDataEtaPhiSubLead_60vtx->SetAxisRange(0,375,"Z"); 
 hUnMatchedDataEtaPhiSubLead_60vtx->GetZaxis()->SetLabelSize(0.035);
 TLegend *slegincleta60 = new TLegend(.14,.60,.44,.88);
  slegincleta60->SetFillStyle(0);
  slegincleta60->SetBorderSize(0);
  slegincleta60->SetTextSize(0.034);
  slegincleta60->AddEntry("","#it{#bf{sPHENIX}} Internal","");
  slegincleta60->AddEntry("","p+p #sqrt{s}=200 GeV ","");
  slegincleta60->AddEntry("","anti-#it{k}_{#it{t}} #it{R} = 0.4","");
  slegincleta60->AddEntry("","Subleading Jets","");
  slegincleta60->AddEntry("","p_{T,1} > 15 GeV , p_{T,2} > 8 GeV","");
  slegincleta60->AddEntry("","|z_{vtx}| < 60 cm, |#eta| < 0.7","");
  slegincleta60->AddEntry("","|#Delta #phi| > 3#pi/4","");
  slegincleta60->Draw();

}//close macro

void SetHist(TH1* h, string xt ="", string yt ="",int color = 1, int marker = 20,int width = 3, float size = 1.0)
{
    h->SetLineWidth(width);
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerSize(size);
    h->SetMarkerStyle(marker);
    h->GetYaxis()->SetTitle(yt.c_str());
    h->GetYaxis()->SetTitleOffset(1.6);
    h->GetXaxis()->SetTitle(xt.c_str());
}

void SetHist(TH1* h, int color = 1)
{
    h->SetLineWidth(3);
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerSize(1);
    h->GetYaxis()->SetTitleOffset(1.6);
    
}
 

void SetLeg(TLegend* l,float txtsize=0.03){
    l->SetBorderSize(0);
    l->SetFillColor(0);
    l->SetTextSize(txtsize);
}
