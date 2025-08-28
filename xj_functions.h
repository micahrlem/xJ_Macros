#ifndef XJ_FUNCTIONS_H
#define XJ_FUNCTIONS_H

#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPad.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TColor.h"
namespace xj_functions

{

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

// Function to create and draw a legend with customizable options
 void DrawsPHENIXLegend(TLegend*& legend, const char* legendName, double x1, double y1, double x2, double y2, const char* Runs =  "p+p #sqrt{s}=200 GeV Runs 47352-47733", bool truth = true, bool reco = true, bool data = true, bool xj = false) {
    // Create the legend
    legend = new TLegend(x1, y1, x2, y2);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.028);

    // Add entries to the legend
    legend->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
    legend->AddEntry("",Runs, "");
    legend->AddEntry("", "anti-#it{k}_{#it{t}} #it{R} = 0.4", "");
    if (truth & !xj){
    legend->AddEntry("", "p_{T,1, truth} > 14 GeV , p_{T,2, truth} > 5.5 GeV", "");
    }
    if (reco & !xj){
    legend->AddEntry("", "p_{T,1, reco} > 18.3 GeV , p_{T,2, reco} > 8.2 GeV", "");
    }
    if (data & !xj){
    legend->AddEntry("", "p_{T,1, data} > 18.3 GeV , p_{T,2, data} > 8.2 GeV", "");
    }
    if (xj){
       legend->AddEntry("", "20.9 < p_{T,1} < 31.2 GeV", "");
        legend->AddEntry("", "p_{T,2} > 9.4 GeV", "");

    }
    legend->AddEntry("", "|#eta| < 0.7, |#Delta #phi| > 3#pi/4", "");
    legend->AddEntry("", "|z_{vtx}| < 60 cm", "");  // Customizable z-vtx condition

    // Draw the legend
    legend->Draw();
}


 void DrawsPHENIXLegendNew(TLegend*& legend, const char* legendName, double x1, double y1, double x2, double y2, const char* Runs =  "p+p #sqrt{s}=200 GeV Runs 47352-47733", const char* RValue = "anti-#it{k}_{#it{t}} #it{R} = 0.4 ", const char* xjrange = "20.9 < p_{T,1} < 31.2 GeV", const char* dijets = "Matched Truth Dijets", bool truth = true, bool reco = true, bool data = true, bool xj = false, bool dijet = false) {
    // Create the legend
    legend =new TLegend(x1, y1, x2, y2);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.028);

    // Add entries to the legend
    legend->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
    legend->AddEntry("",Runs, "");
    legend->AddEntry("",RValue, "");
    if (dijet){
      legend->AddEntry("",dijets,"");
    }
    if (truth & !xj){
    legend->AddEntry("", "p_{T,1, truth} > 14 GeV , p_{T,2, truth} > 5.5 GeV", "");
    }
    if (reco & !xj){
    legend->AddEntry("", "p_{T,1, reco} > 18.3 GeV , p_{T,2, reco} > 8.2 GeV", "");
    }
    if (data & !xj){
    legend->AddEntry("", "p_{T,1, data} > 18.3 GeV , p_{T,2, data} > 8.2 GeV", "");
    }
    if (xj){
   
       legend->AddEntry("", xjrange, "");
        legend->AddEntry("", "p_{T,2} > 9.4 GeV", "");

    }
    legend->AddEntry("", "|#eta| < 1.1 - R, |#Delta#phi| > 3#pi/4", "");
    legend->AddEntry("", "|z_{vtx}| < 60 cm", "");  // Customizable z-vtx condition

    // Draw the legend
    legend->Draw();
}


 void plotxj_and_ratio(TCanvas* can, TH1D* h1, TH1D *horig, TH1D *hdiv, TH1D *hextra, bool extra = false, const char* xaxis = "x_{J}", const char* yaxisratio = "Unfolded/Truth", double x1 = 0,double x2 = 1.0, int col1 = 9, int mark1 = 20, int col2 = 1, int mark2 = 24, int col3 = 9, int mark3 = 24, int col4 = 1, int mark4 = 20){
   can->Divide(1,2); can->cd(1);
    gPad->SetPad(x1, 0.3, x2, 1.0);  // Adjust top pad to leave space for the ratio plot
 // Reduce bottom margin to minimize space between plots
    gPad->SetBottomMargin(0.001);
    SetHist(h1, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col1, mark1, 3, 1.5);
    h1->SetAxisRange(0.31, 0.99, "X");// Set X range from 0 to 0.99
    h1->SetAxisRange(0.01, 5, "Y");h1->GetYaxis()->SetTitleSize(0.04);
    h1->Draw("P");
     SetHist(hdiv, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col3, mark3, 3, 1.5);
     hdiv->Draw("P SAME");
    SetHist(horig, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col2, mark2, 3, 1.5);
    horig->Draw("P SAME");

  
      if(extra){
	 SetHist(hextra, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col4, mark4, 3, 1.5);
	 hextra->Draw("P SAME");
    }

      //ratio canvas
        can->cd(2);
    gPad->SetPad(x1, 0.0, x2, 0.3);  // Adjust bottom pad size for ratio plot
    gPad->SetTopMargin(0.001);
  // Reduce top margin to bring plots closer
    gPad->SetBottomMargin(0.25);  // Increase bottom margin to ensure axis labels are readable

    //define ratio and divide
   TH1D* hratio = (TH1D*)horig->Clone("hratio"); hratio->Divide(hdiv);
    hratio->SetTitle("");  // Remove title to keep it clean
    // hratio->GetXaxis()->SetTitle(xTitle);
    hratio->GetYaxis()->SetTitle(yaxisratio);
    hratio->SetAxisRange(0.31,0.99,"X");
    hratio->SetAxisRange(0, 1.99, "Y");  // Set Y axis range for ratio plot
    hratio->Draw("PE");
    hratio->SetLabelSize(.08, "XY");
    hratio->GetXaxis()->SetTitleSize(.08); 
    hratio->GetYaxis()->SetTitleSize(.08);
    hratio->GetYaxis()->SetTitleOffset(0.5);
    hratio->SetMarkerStyle(20);
    hratio->SetMarkerColor(1);hratio->SetLineColor(1);
    hratio->SetMarkerSize(1.5);
    TLine *line = new TLine(0.3,1,1,1);
    line->SetLineStyle(9);line->SetLineColor(2);
    line->Draw("same");

    can->cd(1);

    

 }


  void plotxj_and_tworatios(TCanvas* can, TH1D* h1, TH1D *horig, TH1D *hdiv, TH1D *hextra, bool extra = false, const char* xaxis = "x_{J}", const char* yaxisratio = "Unfolded/Truth", double x1 = 0,double x2 = 1.0, int col1 = 9, int mark1 = 20, int col2 = 1, int mark2 = 24, int col3 = 9, int mark3 = 24, int col4 = 1, int mark4 = 20){
   can->Divide(1,2); can->cd(1);
    gPad->SetPad(x1, 0.3, x2, 1.0);  // Adjust top pad to leave space for the ratio plot
 // Reduce bottom margin to minimize space between plots
    gPad->SetBottomMargin(0.001);
    SetHist(h1, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col1, mark1, 3, 1.5);
    h1->SetAxisRange(0.31, 0.99, "X");// Set X range from 0 to 0.99
    h1->SetAxisRange(0.01, 5, "Y");h1->GetYaxis()->SetTitleSize(0.04);
    h1->Draw("P");
     SetHist(hdiv, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col3, mark3, 3, 1.5);
     hdiv->Draw("P SAME");
    SetHist(horig, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col2, mark2, 3, 1.5);
    horig->Draw("P SAME");

  
      if(extra){
	 SetHist(hextra, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col4, mark4, 3, 1.5);
	 hextra->Draw("P SAME");
    }

      //ratio canvas
        can->cd(2);
    gPad->SetPad(x1, 0.0, x2, 0.3);  // Adjust bottom pad size for ratio plot
    gPad->SetTopMargin(0.001);
  // Reduce top margin to bring plots closer
    gPad->SetBottomMargin(0.25);  // Increase bottom margin to ensure axis labels are readable

    //define ratio and divide
   TH1D* hratio = (TH1D*)horig->Clone("hratio"); hratio->Divide(hdiv);
   TH1D* hratio2 = (TH1D*)horig->Clone("hratio2"); hratio2->Divide(hextra);
    hratio->SetTitle("");  // Remove title to keep it clean
    // hratio->GetXaxis()->SetTitle(xTitle);
    hratio->GetYaxis()->SetTitle(yaxisratio);
    hratio->SetAxisRange(0.31,0.99,"X");
    hratio->SetAxisRange(0, 1.99, "Y");  // Set Y axis range for ratio plot
    hratio->Draw("PE");
    hratio2->Draw("PE Same");
    hratio->SetLabelSize(.08, "XY");
    hratio->GetXaxis()->SetTitleSize(.08); 
    hratio->GetYaxis()->SetTitleSize(.08);
    hratio->GetYaxis()->SetTitleOffset(0.5);
    hratio->SetMarkerStyle(mark3);
    hratio->SetMarkerColor(col3);hratio->SetLineColor(col3);
    hratio2->SetMarkerColor(col4);hratio2->SetLineColor(col4);
    hratio2->SetMarkerStyle(mark4);
    hratio->SetMarkerSize(1.5);
    TLine *line = new TLine(0.3,1,1,1);
    line->SetLineStyle(9);line->SetLineColor(2);
    line->Draw("same");

    can->cd(1);

    

 }

  void plotxj_and_compare_new(TCanvas* can, TH1D* h1, TH1D *hdiv1, TH1D *h2, TH1D *hdiv2, const char* xaxis = "x_{J}", const char* yaxisratio = "Unfolded/Truth", double x1 = 0,double x2 = 1.0, int col1 = 9, int mark1 = 20, int col2 = 1, int mark2 = 24, int col3 = 9, int mark3 = 24, int col4 = 1, int mark4 = 20){
   can->Divide(1,2); can->cd(1);
    gPad->SetPad(x1, 0.3, x2, 1.0);  // Adjust top pad to leave space for the ratio plot
 // Reduce bottom margin to minimize space between plots
    gPad->SetBottomMargin(0.001);
    SetHist(h1, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col1, mark1, 3, 1.5);
    h1->SetAxisRange(0.31, 0.99, "X");// Set X range from 0 to 0.99
    h1->SetAxisRange(0.01, 5, "Y");h1->GetYaxis()->SetTitleSize(0.04);
    h1->Draw("P");
     SetHist(h2, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col3, mark3, 3, 1.5);
     h2->Draw("P SAME");
    SetHist(hdiv1, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col2, mark2, 3, 1.5);
    hdiv1->Draw("P SAME");

  

       SetHist(hdiv2, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col4, mark4, 3, 1.5);
	 hdiv2->Draw("P SAME");
    

      //ratio canvas
        can->cd(2);
    gPad->SetPad(x1, 0.0, x2, 0.3);  // Adjust bottom pad size for ratio plot
    gPad->SetTopMargin(0.001);
  // Reduce top margin to bring plots closer
    gPad->SetBottomMargin(0.25);  // Increase bottom margin to ensure axis labels are readable

    //define ratio and divide
   TH1D* hratio = (TH1D*)h1->Clone("hratio"); hratio->Divide(hdiv1);
   TH1D* hratio2 = (TH1D*)h2->Clone("hratio2"); hratio2->Divide(hdiv2);
    hratio->SetTitle("");  // Remove title to keep it clean
    // hratio->GetXaxis()->SetTitle(xTitle);
    hratio->GetYaxis()->SetTitle(yaxisratio);
    hratio->SetAxisRange(0.31,0.99,"X");
    hratio->SetAxisRange(0, 1.99, "Y");  // Set Y axis range for ratio plot
    hratio->Draw("PE");
    hratio2->Draw("PE Same");
    hratio->SetLabelSize(.08, "XY");
    hratio->GetXaxis()->SetTitleSize(.08); 
    hratio->GetYaxis()->SetTitleSize(.08);
    hratio->GetYaxis()->SetTitleOffset(0.5);
    hratio->SetMarkerStyle(mark2);
    hratio->SetMarkerColor(col2);hratio->SetLineColor(col2);
    hratio2->SetMarkerColor(col4);hratio2->SetLineColor(col4);
    hratio2->SetMarkerStyle(mark4);
    hratio->SetMarkerSize(1.5);
    TLine *line = new TLine(0.3,1,1,1);
    line->SetLineStyle(9);line->SetLineColor(2);
    line->Draw("same");

    can->cd(1);

    

 }
 
void plot_xj_and_ratio_compare(
    TCanvas* can, TH1D* h1a, TH1D* h1b, TH1D* horig1, TH1D* hdiv1, TH1D* horig2, TH1D* hdiv2,TH1D* hextra = nullptr,
    bool extra = false, const char* xaxis = "x_{J}", const char* yaxisratio = "Unfolded/Truth",double x1 = 0, double x2 = 1.0,
    int col1a = 9, int mark1a = 20, int col1b = 8, int mark1b = 21,int col2 = 1, int mark2 = 24, int col3 = 9, int mark3 = 24,int col4 = 1, int mark4 = 20, int col5 = 46, int mark5 = 25, int col6 = 38, int mark6 = 26)
{
    can->Divide(1, 2);

    // === Upper plot ===
    can->cd(1);
    gPad->SetPad(x1, 0.3, x2, 1.0);
    gPad->SetBottomMargin(0.001);

    // First h1 histogram
    SetHist(h1a, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col1a, mark1a, 3, 1.5);
    h1a->SetAxisRange(0.31, 0.99, "X");
    h1a->SetAxisRange(0.01, 5, "Y");
    h1a->GetYaxis()->SetTitleSize(0.04);
    h1a->Draw("P");

    // Second h1 histogram
    SetHist(h1b, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col1b, mark1b, 3, 1.5);
    h1b->Draw("P SAME");

    // First div/orig pair
    SetHist(hdiv1, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col3, mark3, 3, 1.5);
    hdiv1->Draw("P SAME");
    SetHist(horig1, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col2, mark2, 3, 1.5);
    horig1->Draw("P SAME");

    // Second div/orig pair
    SetHist(hdiv2, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col5, mark5, 3, 1.5);
    hdiv2->Draw("P SAME");
    SetHist(horig2, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col4, mark4, 3, 1.5);
    horig2->Draw("P SAME");

    // Optional extra histogram
    if (extra && hextra) {
        SetHist(hextra, xaxis, "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", col6, mark6, 3, 1.5);
        hextra->Draw("P SAME");
    }

    // === Ratio plot ===
    can->cd(2);
    gPad->SetPad(x1, 0.0, x2, 0.3);
    gPad->SetTopMargin(0.001);
    gPad->SetBottomMargin(0.25);

    // First ratio
    TH1D* hratio1 = (TH1D*)horig1->Clone("hratio1");
    hratio1->Divide(hdiv1);
    hratio1->SetTitle("");
    hratio1->GetYaxis()->SetTitle(yaxisratio);
    hratio1->SetAxisRange(0.31, 0.99, "X");
    hratio1->SetAxisRange(0, 1.99, "Y");
    hratio1->SetLabelSize(0.08, "XY");
    hratio1->GetXaxis()->SetTitleSize(0.08);
    hratio1->GetYaxis()->SetTitleSize(0.08);
    hratio1->GetYaxis()->SetTitleOffset(0.5);
    hratio1->SetMarkerStyle(mark2);
    hratio1->SetMarkerColor(col2);
    hratio1->SetLineColor(col2);
    hratio1->SetMarkerSize(1.5);
    hratio1->Draw("PE");

    // Second ratio
    TH1D* hratio2 = (TH1D*)horig2->Clone("hratio2");
    hratio2->Divide(hdiv2);
    hratio2->SetMarkerStyle(mark4);
    hratio2->SetMarkerColor(col4);
    hratio2->SetLineColor(col4);
    hratio2->SetMarkerSize(1.5);
    hratio2->Draw("PE SAME");

    // Reference line at y = 1
    TLine* line = new TLine(0.1, 1, 0.99, 1);
    line->SetLineStyle(9);
    line->SetLineColor(2);
    line->Draw("same");

    can->cd(1);
}

 void makeanddraw_xj_ratio(TCanvas* can, TH1D *horig, TH1D *hdiv, const char* yaxis = "Unfolded/Truth", int can1 = 1, int can2 = 2, double x1 = 0,double x2 = 1.0){
   
    can->cd(can1);
    gPad->SetPad(x1, 0.3, x2, 1.0);  // Adjust top pad to leave space for the ratio plot
 // Reduce bottom margin to minimize space between plots
    gPad->SetBottomMargin(0.001);
    can->cd(can2);
    gPad->SetPad(x1, 0.0, x2, 0.3);  // Adjust bottom pad size for ratio plot
      gPad->SetTopMargin(0.001);
  // Reduce top margin to bring plots closer
    gPad->SetBottomMargin(0.25);  // Increase bottom margin to ensure axis labels are readable

    //define ratio and divide
   TH1D* hratio = (TH1D*)horig->Clone("hratio"); hratio->Divide(hdiv);
    hratio->SetTitle("");  // Remove title to keep it clean
    // hratio->GetXaxis()->SetTitle(xTitle);
    hratio->GetYaxis()->SetTitle(yaxis);
    hratio->SetAxisRange(0.1,0.99,"X");
    hratio->SetAxisRange(0, 1.99, "Y");  // Set Y axis range for ratio plot
    hratio->Draw("PE");
    hratio->SetLabelSize(.08, "XY");
    hratio->GetXaxis()->SetTitleSize(.08); 
    hratio->GetYaxis()->SetTitleSize(.08);
    hratio->GetYaxis()->SetTitleOffset(0.5);
    hratio->SetMarkerStyle(20);
    hratio->SetMarkerColor(1);
    hratio->SetMarkerSize(1.5);
    TLine *line = new TLine(0.1,1,1,1);
    line->SetLineStyle(9);line->SetLineColor(2);
    line->Draw("same");
  

 }

void draw_xj_legend(TLegend*& legend, const char* legendName, double x1, double y1, double x2, double y2, TH1D *h1,TH1D *h2,TH1D *h3,TH1D *h4, const char* label1 = " ", const char* label2 = " ", const char* label3 = " ", const char* label4 = " ", bool extra = false){
     // Create the legend
    legend = new TLegend(x1, y1, x2, y2);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.034);

    legend->AddEntry(h1,label1,"P");
    legend->AddEntry(h2,label2,"P");
    legend->AddEntry(h3,label3,"P");
    if(extra){
    legend->AddEntry(h4,label4,"P");
    }
    legend->Draw();

 }

void draw_xj_legend_compare(
    TLegend*& legend,
    const char* legendName,
    double x1, double y1, double x2, double y2,
    TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4 = nullptr,
    TH1D* h5 = nullptr, TH1D* h6 = nullptr, TH1D* h7 = nullptr,
    const char* label1 = " ", const char* label2 = " ", const char* label3 = " ", const char* label4 = " ",
    const char* label5 = " ", const char* label6 = " ", const char* label7 = " ",
    bool extra = false,
    bool compare = false)
{
    legend = new TLegend(x1, y1, x2, y2);
    legend->SetName(legendName);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.034);

    legend->AddEntry(h1, label1, "P");
    legend->AddEntry(h2, label2, "P");
    legend->AddEntry(h3, label3, "P");

    if (extra && h7)
        legend->AddEntry(h7, label7, "P");

    if (compare) {
        if (h4) legend->AddEntry(h4, label4, "P");
        if (h5) legend->AddEntry(h5, label5, "P");
        if (h6) legend->AddEntry(h6, label6, "P");
    }

    legend->Draw();
}


void get_xj_systematics(TH1D *h1, TH1D *hs, const int nbins)
  {
    for (int i = 0; i < nbins; i++)
      {
	h1->SetBinError(i+1, hs->GetBinContent(i+1)*h1->GetBinContent(i+1));
      }
  }
  void set_xj_errors(TH1D *h1, TProfile *hp, const int nbins)
  {
    for (int i = 0; i < nbins; i++)
      {
	float old_error = h1->GetBinError(i+1);
	float stat_error = hp->GetBinError(i+1);
	float new_error = sqrt(TMath::Power(old_error, 2) + TMath::Power(stat_error, 2));
	h1->SetBinError(i+1, new_error);
      }
  }


  void finalize_xj(TH1D *h1, TH1D *h2, const int nbins, float first_xj)
  {
    for (int i = 0; i < nbins; i++)
      {
	if (h1->GetBinLowEdge(i+1) >= first_xj)
	  {
	    h2->SetBinContent(i+1, h1->GetBinContent(i+1));
	    h2->SetBinError(i+1, h1->GetBinError(i+1));
	  }
      }
  }

void normalize_histo(TH1D* h, const int nbins)
  {
    Double_t bin_contents[30] = {0};
    Double_t bin_errors[30] = {0};

    Double_t integral = 0;

    for (int i = 0; i < nbins; i++)
      {

        Double_t v = h->GetBinContent(i+1);
        Double_t w = h->GetBinWidth(i+1);
        Double_t err = h->GetBinError(i+1);

        Double_t vw = v/w;
        Double_t ew = err/w;

        integral += v;

        bin_contents[i] = v/w;
        bin_errors[i] = ew;

      }

    for (int i = 0; i < nbins; i++)
      {
        h->SetBinContent(i+1, bin_contents[i]/integral);
        h->SetBinError(i+1, bin_errors[i]/integral);
      }
    return;
  }

 void clean_pt1pt2(TH2D* hpt1pt2, const int nbins, const float weight){

    for (int ix = 0; ix < nbins + 1; ix++)
      {
	for (int iy = 0; iy < nbins + 1; iy++)
	  {
	    float content  = hpt1pt2->GetBinContent(ix+1, iy+1);
	    float max = weight*5;

	    if (content <= max)
	      {
		hpt1pt2->SetBinContent(ix+1,iy+1,0);
	
	      }
	    else 
	      {
		continue;
	      }
	
	      
	  }
      }


   }

  void skim_down_histo(TH1D *h_skim, TH1D *h_full, TH1D *h_mapping)
  {

    int nbins = h_mapping->GetXaxis()->GetNbins();

    for (int ib = 0; ib < nbins*nbins; ib++)
      {
	if (!(h_mapping->GetBinContent(ib+1) == 0))
	  {
	    h_skim->SetBinContent(h_mapping->GetBinContent(ib+1), h_full->GetBinContent(ib+1));
	    h_skim->SetBinError(h_mapping->GetBinContent(ib+1), h_full->GetBinError(ib+1));
	  }
      }
  }

   void make_sym_pt1pt2(TH1D *hflat, TH2D* hpt1pt2, const int nbins)
  {

    for (int ib = 0; ib < nbins*nbins; ib++)
      {
	int xbin = ib/nbins;
	int ybin = ib%nbins;
      
	int b = hpt1pt2->GetBin(xbin+1, ybin+1);

	hpt1pt2->SetBinContent(b, hflat->GetBinContent(ib+1));

	hpt1pt2->SetBinError(b, hflat->GetBinError(ib+1));
      }
    return;
  }

   void project_xj(TH2D* hpt1pt2, TH1D* h_xj, const int nbins, const int start_leading_bin, const int end_leading_bin, const int start_subleading_bin, const int end_subleading_bin)
  {

    TH1D *h_unc = (TH1D*) h_xj->Clone();
    h_unc->Reset();

    TH2D *h_asym_pt1pt2 = (TH2D*) hpt1pt2->Clone();

    std::cout << "Leading pT bin range: "
          << hpt1pt2->GetXaxis()->GetBinLowEdge(start_leading_bin+1)
          << " – "
          << hpt1pt2->GetXaxis()->GetBinUpEdge(end_leading_bin)
          << std::endl;
  
    for (int ix = 0; ix < nbins; ix++)
      {
	for (int iy = 0; iy < nbins; iy++)
	  {
	    int bin = h_asym_pt1pt2->GetBin(ix+1, iy+1);

	    if (ix > iy)
	      {
		h_asym_pt1pt2->SetBinContent(bin, h_asym_pt1pt2->GetBinContent(bin)*2.);
		h_asym_pt1pt2->SetBinError(bin, h_asym_pt1pt2->GetBinError(bin)*2);
	      }
	    else if (ix < iy)
	      {
		h_asym_pt1pt2->SetBinContent(bin, 0);
		h_asym_pt1pt2->SetBinError(bin, 0);
	      }
	
	      
	  }
      }

    for (int ix = 0; ix < nbins; ix++)
      {
	for (int iy = 0; iy <= ix; iy++)
	  {
	    int low =  iy - ix - 1;
	    int high = iy - ix + 1;

	    int xjbin_low = nbins + low + 1;
	    int xjbin_high = nbins + low + 2;
	    //std::cout << ix << " / " << iy << " -- > " << xjbin_low << "--"<<xjbin_high<<std::endl;
	    int bin = h_asym_pt1pt2->GetBin(ix+1, iy+1);

	    if (ix < start_leading_bin) continue;
	    if (iy < start_subleading_bin) continue;
	    if (ix >= end_leading_bin) continue;
	    if (iy >= end_subleading_bin) continue;

	    if (ix == iy)
	      {

		h_xj->Fill(h_xj->GetBinCenter(xjbin_low), h_asym_pt1pt2->GetBinContent(bin));
		h_unc->Fill(h_xj->GetBinCenter(xjbin_low),TMath::Power( h_asym_pt1pt2->GetBinError(bin), 2));
	      }
	    else
	      {
		h_xj->Fill(h_xj->GetBinCenter(xjbin_high), h_asym_pt1pt2->GetBinContent(bin)/2.);
		h_xj->Fill(h_xj->GetBinCenter(xjbin_low), h_asym_pt1pt2->GetBinContent(bin)/2.);
		h_unc->Fill(h_xj->GetBinCenter(xjbin_high),TMath::Power( h_asym_pt1pt2->GetBinError(bin)/sqrt(2), 2));
		h_unc->Fill(h_xj->GetBinCenter(xjbin_low),TMath::Power( h_asym_pt1pt2->GetBinError(bin)/sqrt(2), 2));
	      }
	  }
      }

    for (int ix = 0; ix < nbins; ix++)
      {
	h_xj->SetBinError(ix+1, sqrt(h_unc->GetBinContent(ix+1)));
      }

    return;
  }

   
   
}
#endif
