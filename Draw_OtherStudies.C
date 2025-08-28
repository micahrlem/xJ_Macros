#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

#include "xj_functions.h"
#include "read_binning.h"

#include <TFile.h>
#include <TH1D.h>
#include <TKey.h>
#include <TClass.h>
#include <map>
#include <vector>
#include <string>
#include <iostream>

#include <algorithm>
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TBox.h"
#include "TColor.h"

// Blend a base ROOT color toward white by factor w in [0,1].
// w=0 -> original color; w=1 -> pure white.
// Replace your BlendWithWhite with this version:
#include "TColor.h"
#include "TROOT.h"

// Lighten a ROOT color by blending toward white.
// whiten in [0,1]: 0 -> ignore whitening for this base; 1 -> fully allow whitening
// alpha  in [0,1]: 0 -> no lightening now; 1 -> maximum white now
static int BlendWithWhite(int baseColor, double whiten, double alpha)
{
  // clamp inputs
  whiten = std::max(0.0, std::min(1.0, whiten));
  alpha  = std::max(0.0, std::min(1.0, alpha));
  const double w = whiten * alpha;

  TColor* c = gROOT->GetColor(baseColor);
  double r = 0.0, g = 0.0, b = 0.0;
  if (c) { r = c->GetRed(); g = c->GetGreen(); b = c->GetBlue(); }

  // linear blend toward white by w
  const double R = r*(1.0 - w) + 1.0*w;
  const double G = g*(1.0 - w) + 1.0*w;
  const double B = b*(1.0 - w) + 1.0*w;

  // disambiguate overload (Float_t)
  return TColor::GetColor((Float_t)R, (Float_t)G, (Float_t)B);
}


TCanvas* DrawArraysAndHistsTogether(
	      // ----- array inputs -----
    const std::vector<double>& edges,              // size = nA+1
    const std::vector<double>& yCentroidArr,       // size = nA
    const std::vector<double>& yLowArr,            // size = nA
    const std::vector<double>& yHighArr,           // size = nA
    // ----- histogram inputs (binning may differ) -----
    const TH1* hCent, const TH1* hLow, const TH1* hHigh,
    // ----- legend labels -----
    const char* arrLegendLabel,                    // e.g. "Arrays (20–30%)"
    const char* histLegendLabel,                   // e.g. "Hists (20–30%)"
    // ----- cosmetics -----
    const char* xTitle = "x_{J}",
    const char* yTitle = "Yield",
    const char* tag    = nullptr,                  // unique suffix
    double gapFrac = 1e-4                          // tiny gap to avoid touching
){
  // Basic checks for arrays
  const int nA = (int)edges.size() - 1;
  if ((int)yCentroidArr.size()!=nA || (int)yLowArr.size()!=nA || (int)yHighArr.size()!=nA) {
    ::Error("DrawArraysAndHistsTogether", "Array sizes mismatch with edges");
    return nullptr;
  }
  // Basic checks for hists
  if (!hCent || !hLow || !hHigh) {
    ::Error("DrawArraysAndHistsTogether", "Null histogram pointer(s)");
    return nullptr;
  }
  const int nH = hCent->GetNbinsX();
  if (hLow->GetNbinsX()!=nH || hHigh->GetNbinsX()!=nH) {
    ::Error("DrawArraysAndHistsTogether", "Histogram bin count mismatch among hist inputs");
    return nullptr;
  }

  // Unique suffix
  static int call_id = 0;
  TString suf = (tag && *tag) ? TString::Format("_%s_%d", tag, call_id++) : TString::Format("_%d", call_id++);

  // Fixed visible ranges
  const double XMIN_VIS = 0.31;
  const double XMAX_VIS = 0.99;
  const double YMIN_VIS = 0.01;
  const double YMAX_VIS = 5.0;

  // Precompute array bin centers
  std::vector<double> xcentArr(nA);
  for (int i=0; i<nA; ++i) xcentArr[i] = 0.5*(edges[i]+edges[i+1]);

  // Precompute hist bin centers
  std::vector<double> xcentHist(nH), ycentHist(nH);
  for (int i=1;i<=nH;++i){
    xcentHist[i-1] = hCent->GetXaxis()->GetBinCenter(i);
    ycentHist[i-1] = hCent->GetBinContent(i);
  }
  gStyle->SetCanvasPreferGL(kTRUE);

  // Canvas & frame (fresh TH1D to avoid Sumw2 warnings)
  TCanvas* c = new TCanvas(TString("c_xj_combo")+suf, "", 800, 800);
  gStyle->SetOptStat(0);

  TH1D* frame = new TH1D(TString("frame")+suf, "", 1, XMIN_VIS, XMAX_VIS);
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle(xTitle);
  frame->GetYaxis()->SetTitle(yTitle);
  frame->SetMinimum(YMIN_VIS);
  frame->SetMaximum(YMAX_VIS);
 
  frame->Draw();

   frame->GetXaxis()->SetTitleSize(0.035);
  frame->GetYaxis()->SetTitleSize(0.035);
  frame->GetXaxis()->SetLabelSize(0.030);
  frame->GetYaxis()->SetLabelSize(0.030);

  frame->GetYaxis()->SetTitleOffset(1.8);

  // Right-side ticks
  gPad->SetTicks(0,1);

  // Gaps based on each object’s *own* full range
  const double arrFullW = (edges.back() - edges.front());
  const double GAP_arr  = gapFrac * (arrFullW > 0 ? arrFullW : (XMAX_VIS - XMIN_VIS));

  const double hxmin = hCent->GetXaxis()->GetXmin();
  const double hxmax = hCent->GetXaxis()->GetXmax();
  const double histFullW = (hxmax - hxmin);
  const double GAP_hist  = gapFrac * (histFullW > 0 ? histFullW : (XMAX_VIS - XMIN_VIS));

  // =========================
  // 1) ARRAY band (solid pale) + markers, clipped to visible X
  // =========================
 // once, outside the loop
  const int arrFillColor = kAzure-5;  // 15% opacity
 int arrLineColor  = kAzure-5;



  std::vector<TBox*> arrBoxes; arrBoxes.reserve(nA);
  for (int i=0;i<nA;++i){
    double xl = edges[i]   + GAP_arr;
    double xr = edges[i+1] - GAP_arr;

    // clamp to frame X range
    if (xr < XMIN_VIS || xl > XMAX_VIS) continue;
    xl = std::max(xl, XMIN_VIS);
    xr = std::min(xr, XMAX_VIS);
    if (xr <= xl) continue;

    double yl = yLowArr[i];
    double yh = yHighArr[i];

    TBox* b = new TBox(xl, yl, xr, yh);
    b->SetFillStyle(1001);         // solid
    b->SetFillColor(arrFillColor);
    b->SetFillColorAlpha(arrFillColor, 0.20);
    b->SetLineColor(arrLineColor);
    b->SetLineWidth(1);
    b->Draw("F");
    arrBoxes.push_back(b);
  }

 

  TGraphAsymmErrors* gArr = new TGraphAsymmErrors(nA);
for (int i = 0; i < nA; ++i) {
  double xl  = edges[i];
  double xr  = edges[i+1];
  double xc  = 0.5 * (xl + xr);
  double exl = xc - xl;  // left half-width
  double exh = xr - xc;  // right half-width

  gArr->SetPoint(i, xc, yCentroidArr[i]);
  gArr->SetPointEXlow (i, exl);
  gArr->SetPointEXhigh(i, exh);
  gArr->SetPointEYlow (i, 0.0);
  gArr->SetPointEYhigh(i, 0.0);
}

gArr->SetName(TString("gArrayCentroid") + suf);
gArr->SetMarkerStyle(20);
gArr->SetMarkerSize(1.5);
gArr->SetMarkerColor(kAzure-5);
gArr->SetLineColor(kAzure-5);
gArr->SetLineWidth(2);

// Remove little "hats" on the ends of the error bars
gStyle->SetEndErrorSize(3);

// Draw points with error bars, no caps
gArr->Draw("PE");

  // =========================
  // 2) HIST band (hatched) + markers, clipped to visible X
  // =========================
  const int hatchFillStyle = 3335;
 // Hists: base red, allow full whitening, apply alpha = 0.50 (less pale)
int histFillColor = kBlack;
int histLineColor = kBlack;

  std::vector<TBox*> histBoxes; histBoxes.reserve(nH);
  for (int i=1;i<=nH;++i){
    double xl = hCent->GetXaxis()->GetBinLowEdge(i) + GAP_hist;
    double xr = hCent->GetXaxis()->GetBinUpEdge(i)  - GAP_hist;

    // clamp to frame X range
    if (xr < XMIN_VIS || xl > XMAX_VIS) continue;
    xl = std::max(xl, XMIN_VIS);
    xr = std::min(xr, XMAX_VIS);
    if (xr <= xl) continue;

    double yl = hLow ->GetBinContent(i);
    double yh = hHigh->GetBinContent(i);

    TBox* b = new TBox(xl, yl, xr, yh);
    b->SetFillStyle(hatchFillStyle);  // pattern (robust in batch)
    b->SetFillColor(histFillColor);
    b->SetLineColor(histLineColor);
     b->SetFillColorAlpha(histFillColor, 0.9);
    b->SetLineWidth(1);
    b->Draw("F");
    histBoxes.push_back(b);
  }

TGraphAsymmErrors* gHist = new TGraphAsymmErrors(nH);
for (int i = 0; i < nH; ++i) {
  double xl  = hCent->GetXaxis()->GetBinLowEdge(i+1);
  double xr  = hCent->GetXaxis()->GetBinUpEdge(i+1);
  double xc  = 0.5 * (xl + xr);
  double exl = xc - xl;  // left half-width
  double exh = xr - xc;  // right half-width

  gHist->SetPoint(i, xc, ycentHist[i]);
  gHist->SetPointEXlow (i, exl);
  gHist->SetPointEXhigh(i, exh);
  gHist->SetPointEYlow (i, 0.0);
  gHist->SetPointEYhigh(i, 0.0);
}

gHist->SetName(TString("gHistCentroid") + suf);
gHist->SetMarkerStyle(20);
gHist->SetMarkerSize(1.5);
gHist->SetMarkerColor(kBlack);
gHist->SetLineColor(kBlack);
gHist->SetLineWidth(2);

// Caps for the ends of horizontal error bars
gStyle->SetEndErrorSize(3); // tweak to taste

gHist->Draw("PE");

  // =========================
  // Legend (two custom entries)
  // =========================
  TLegend* leg = new TLegend(0.20, 0.58, 0.50, 0.63);
  leg->SetName(TString("leg_combo")+suf);
  leg->SetTextSize(0.023);
  leg->SetBorderSize(0);
  /*if (!arrBoxes.empty()) leg->AddEntry(arrBoxes.front(),  arrLegendLabel,  "f");
  else                   leg->AddEntry(gArr,               arrLegendLabel,  "p");
  if (!histBoxes.empty()) leg->AddEntry(histBoxes.front(), histLegendLabel, "f");
  else                    leg->AddEntry(gHist,             histLegendLabel, "p");*/
  // Optional: include marker keys too
  leg->AddEntry(gHist, "Unfolded Data",  "p");
  leg->AddEntry(gArr,  "Unfolded Data QM", "p");
  
  leg->Draw();

  c->RedrawAxis();
  return c;
}


// Pretty pT labels by centrality key
static const std::map<std::string, std::string> pt_labels = {
  {"20_30", "20.9 < p_{T,1} < 31.2 GeV"},
  {"30_40", "31.2 < p_{T,1} < 40.7 GeV"},
  {"40_60", "40.7 < p_{T,1} < 60.8 GeV"}
};

inline void Make_sPHENIX_Legend_for(const std::string& pt_range, TLegend*& leg_out)
{
  // Resolve the pT label from the map (with a safe default)
  auto it = pt_labels.find(pt_range);
  const std::string pt_range_str = (it != pt_labels.end())
      ? it->second
      : "20.9 < p_{T,1} < 60.8 GeV";

  // Unique legend key/name per centrality
  const std::string leg_key = "leg_sphenix_compare_" + pt_range;

  // Draw the legend (positions are your original; tweak if needed)
  xj_functions::DrawsPHENIXLegendNew(
      leg_out,
      leg_key.c_str(),
      /*x1=*/0.15, /*y1=*/0.65, /*x2=*/0.50, /*y2=*/0.90,
      "p+p #sqrt{s}=200 GeV",
      "anti-#it{k}_{#it{t}} #it{R} = 0.4",           // R4 histograms
      pt_range_str.c_str(),
      "Matched Truth Dijets",
      /*drawBox=*/true, /*drawBorder=*/true,
      /*drawShadow=*/false, /*drawFill=*/true,
      /*transparent=*/false
  );
}

void Draw_OtherStudies()

{

   const char* fname10 = "Xj_2D_Response_Bins19_R2_Play_Smear_Flatten_pythia8-Jet10-Run21-multiR.root";
    const char* fname20 = "Xj_2D_Response_Bins19_R2_Play_Smear_Flatten_pythia8-Jet20-Run21-multiR.root";
    const char* fname30 = "Xj_2D_Response_Bins19_R2_Play_Smear_Flatten_pythia8-Jet30-Run21-multiR.root";

    // Open files
    TFile *f10 = TFile::Open(fname10, "READ");
    TFile *f20 = TFile::Open(fname20, "READ");
    TFile *f30 = TFile::Open(fname30, "READ");

    if (!f10 || f10->IsZombie() || !f20 || f20->IsZombie() || !f30 || f30->IsZombie()) {
        Error("add_truth_histos", "One or more files could not be opened.");
        return;
    }

    // Get histograms and clone for accumulation
    TH1D *h_Aj = (TH1D*)f10->Get("h_linear_truth_Aj");
    TH1D *h_xj = (TH1D*)f10->Get("h_linear_truth_xj");

    if (!h_Aj || !h_xj) {
        Error("add_truth_histos", "Could not find h_linear_truth_Aj or h_linear_truth_xj in first file.");
        return;
    }

    TH1D *h_Aj_sum = (TH1D*)h_Aj->Clone("h_Aj_sum");
    TH1D *h_xj_sum = (TH1D*)h_xj->Clone("h_xj_sum");

    // Add from second file
    h_Aj_sum->Add((TH1D*)f20->Get("h_linear_truth_Aj"));
    h_xj_sum->Add((TH1D*)f20->Get("h_linear_truth_xj"));

    // Add from third file
    h_Aj_sum->Add((TH1D*)f30->Get("h_linear_truth_Aj"));
    h_xj_sum->Add((TH1D*)f30->Get("h_linear_truth_xj"));

     gStyle->SetOptStat(0); // fixed from gPad->SetOptStyle(0)

// ---- Normalize helper: PDF (unit area) ----
// Normalize to a *probability density*:
// y_i -> counts_i / (N * bin_width_i), so sum_i y_i * bin_width_i = 1
auto normalize_pdf = [](TH1* h){
    if (!h) return;
    h->Sumw2(kTRUE);                 // ensure errors scale properly
    const double N = h->Integral();  // total counts (no width)
    if (N <= 0) return;
    h->Scale(1.0 / N);               // now contents are per *count* (prob per bin)
    h->Scale(1.0, "width");          // divide each bin by its width => per unit x
    // Now Integral("width") == 1 exactly (within float precision)
};

// ---- Aj canvas ----
TCanvas *c_Aj_sum = new TCanvas("c_Aj_sum", "Summed Truth Histogram: Aj", 800, 800);

// Style pad
auto style_pad = [](TPad* p){
  p->SetLeftMargin(0.16);
  p->SetRightMargin(0.06);
  p->SetBottomMargin(0.12);
  p->SetTopMargin(0.06);
};
style_pad(c_Aj_sum);

// Style histogram helper
auto style_h = [](TH1* h){
  h->SetAxisRange(0.01,5,"Y");// peak at ~60% frame height
  h->GetYaxis()->SetTitleSize(0.032);
  h->GetYaxis()->SetLabelSize(0.028);
  h->GetYaxis()->SetTitleOffset(1.4);
  h->GetXaxis()->SetTitleSize(0.032);
  h->GetXaxis()->SetLabelSize(0.028);
};

// PDF-normalize Aj
normalize_pdf(h_Aj_sum);

// Aj plot (red line)
h_Aj_sum->SetLineColor(kRed);
h_Aj_sum->SetLineWidth(2);
style_h(h_Aj_sum);
h_Aj_sum->GetYaxis()->SetTitle("P(A_{J})");
h_Aj_sum->Draw("HIST");

// Physics block legend
TLegend* Aj_leg = new TLegend(0.43,0.70, 0.73, 0.90);
Aj_leg->SetFillStyle(0);
Aj_leg->SetBorderSize(0);
Aj_leg->SetTextSize(0.028);
Aj_leg->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
Aj_leg->AddEntry("", "#it{R} = 0.2", "");
Aj_leg->AddEntry("", "p_{T,1, truth} > 20 GeV , p_{T,2, truth} > 5.5 GeV", "");
Aj_leg->AddEntry("", "|#eta| < 0.9, |#Delta#phi| > 3#pi/4", "");
Aj_leg->AddEntry("", "|z_{vtx}| < 60 cm", "");
Aj_leg->Draw();

// Label legend (top-left)
TLegend* Aj_leg2 = new TLegend(0.25, 0.85, 0.40, 0.90);
Aj_leg2->SetBorderSize(1);
Aj_leg2->SetFillStyle(0);
Aj_leg2->SetTextSize(0.028);
Aj_leg2->AddEntry(h_Aj_sum, "PYTHIA8", "L");
Aj_leg2->Draw();


// --------------------

// xJ canvas
TCanvas *c_xJ_sum = new TCanvas("c_xJ_sum", "Summed Truth Histogram: xJ", 800, 800);
style_pad(c_xJ_sum);

// PDF-normalize xJ
normalize_pdf(h_xj_sum);

// xJ plot (kGreen+3 line)
h_xj_sum->SetLineColor(kGreen+3);
h_xj_sum->SetLineWidth(2);
style_h(h_xj_sum);
h_xj_sum->GetYaxis()->SetTitle("P(x_{J})");
h_xj_sum->Draw("HIST");

// Physics block legend
TLegend* xJ_sum_leg = new TLegend(0.43,0.70, 0.73, 0.90);
xJ_sum_leg->SetFillStyle(0);
xJ_sum_leg->SetBorderSize(0);
xJ_sum_leg->SetTextSize(0.028);
xJ_sum_leg->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
xJ_sum_leg->AddEntry("", "#it{R} = 0.2", "");
xJ_sum_leg->AddEntry("", "p_{T,1, truth} > 20 GeV , p_{T,2, truth} > 5.5 GeV", "");
xJ_sum_leg->AddEntry("", "|#eta| < 0.9, |#Delta#phi| > 3#pi/4", "");
xJ_sum_leg->AddEntry("", "|z_{vtx}| < 60 cm", "");
xJ_sum_leg->Draw();

// Label legend (top-left)
TLegend* xJ_sum_leg2 = new TLegend(0.25, 0.85, 0.40, 0.90);
xJ_sum_leg2->SetBorderSize(1);
xJ_sum_leg2->SetFillStyle(0);
xJ_sum_leg2->SetTextSize(0.028);
xJ_sum_leg2->AddEntry(h_xj_sum, "PYTHIA8", "L");
xJ_sum_leg2->Draw();

  //Plots with new Axis
  // Open files
TFile *f_p8     = TFile::Open("UnfoldedxJHists_Play_pythia8_R2.root");
TFile *f_hwg    = TFile::Open("UnfoldedxJHists_Play_Herwig_R2.root");
TFile *f_p8pos  = TFile::Open("UnfoldedxJHists_Play_pythia8_posJER_R2.root");
TFile *f_p8neg  = TFile::Open("UnfoldedxJHists_Play_pythia8_negJER_R2.root");

// Grab histograms
TH1D *h_xj_pythia8        = (TH1D*)f_p8    ->Get("h_xj_unfold_data_30_40");
TH1D *h_xj_herwig_rew     = (TH1D*)f_hwg   ->Get("h_xj_unfold_data_30_40");
TH1D *h_xj_pythia8_posJER = (TH1D*)f_p8pos ->Get("h_xj_unfold_data_30_40");
TH1D *h_xj_pythia8_negJER = (TH1D*)f_p8neg ->Get("h_xj_unfold_data_30_40");


 


// Assume you already have:
// TH1D* h_xj_pythia8;
// TH1D* h_xj_herwig_rew;
TCanvas* c_xj = new TCanvas("c_xj","xJ: Pythia8 vs Herwig",800,800);


// Common margin & offset settings
double leftMargin = 0.22;   // increased a bit
double yTitleOffset = 1.8;

// -------- Top pad (spectra) --------
c_xj->cd(1);
 //gPad->SetPad(0.0, 0.30, 1.0, 1.0);
gPad->SetBottomMargin(0.001);
gPad->SetLeftMargin(leftMargin);

// Style & axes
auto style_line = [yTitleOffset](TH1D* h, const char* xtitle, const char* ytitle, int col, int mstyle){
  h->SetTitle("");
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetLabelSize(0.035);
  h->GetYaxis()->SetLabelSize(0.035);
  h->GetYaxis()->SetTitleOffset(yTitleOffset);
  h->SetMarkerStyle(mstyle);
  h->SetMarkerSize(1.5);
  h->SetMarkerColor(col);
  h->SetLineColor(col);
  h->SetLineWidth(3);
};

style_line(h_xj_pythia8, "x_{J}", "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", 6, 20);
style_line(h_xj_herwig_rew, "x_{J}", "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", 9, 53);

h_xj_pythia8->SetAxisRange(0.31, 0.99, "X");
h_xj_pythia8->SetAxisRange(0.01, 2.5, "Y");

h_xj_pythia8->Draw("P");
h_xj_herwig_rew->Draw("P SAME");

// Legend
TLegend* leg_p_h = new TLegend(0.42, 0.78, 0.68, 0.84);
leg_p_h->SetBorderSize(0);
leg_p_h->SetFillStyle(0);
leg_p_h->SetTextSize(0.030);
leg_p_h->SetMargin(0.28);
leg_p_h->AddEntry(h_xj_pythia8, "R=0.2 Data Unfolded w/ PYTHIA8", "p");
leg_p_h->AddEntry(h_xj_herwig_rew, "R=0.2 Data Unfolded w/ HERWIG7", "p");
leg_p_h->Draw();
/*
// -------- Bottom pad (ratio) --------
c_xj->cd(2);
gPad->SetPad(0.0, 0.0, 1.0, 0.30);
gPad->SetTopMargin(0.02);
gPad->SetBottomMargin(0.28);
gPad->SetLeftMargin(leftMargin);

TH1D* h_ratio = (TH1D*)h_xj_pythia8->Clone("h_ratio_p8_over_hwg");
h_ratio->Divide(h_xj_herwig_rew);
h_ratio->GetXaxis()->SetTitle("x_{J}");
h_ratio->GetYaxis()->SetTitle("Pythia8/Herwig");
h_ratio->SetAxisRange(0.31, 0.99, "X");
h_ratio->SetAxisRange(0.0, 1.99, "Y");

h_ratio->GetXaxis()->SetTitleSize(0.08);
h_ratio->GetXaxis()->SetLabelSize(0.08);
h_ratio->GetYaxis()->SetTitleSize(0.08);
h_ratio->GetYaxis()->SetLabelSize(0.08);
h_ratio->GetYaxis()->SetTitleOffset(0.4); // reduced a bit for small pad
h_ratio->GetYaxis()->SetNdivisions(505);

h_ratio->SetMarkerStyle(20);
h_ratio->SetMarkerColor(6);
h_ratio->SetLineColor(6);
h_ratio->SetMarkerSize(1.5);

h_ratio->Draw("PE");

TLine *line = new TLine(0.3,1.0,0.99,1.0);
line->SetLineStyle(9);
line->SetLineColor(2);
line->Draw("SAME");

c_xj->cd(1);
c_xj->Update();

*/

 TLegend* leg_sphenix_xj = nullptr;
 xj_functions::DrawsPHENIXLegendNew(leg_sphenix_xj, "leg_sphenix_xj",
        0.15, 0.52, 0.45, 0.90, "p+p #sqrt{s}=200 GeV",
        "anti-#it{k}_{#it{t}} #it{R} = 0.2","31.2 < p_{T,1} < 40.7 GeV", "Matched Truth Dijets",
        true, true, false, true, false);

TCanvas* c_xj_JER = new TCanvas("c_xj_JER","xJ: Pythia8 vs JER variations",800,800);
c_xj_JER->Divide(1,2);

// ---------- Top pad (spectra) ----------
c_xj_JER->cd(1);
gPad->SetPad(0.0, 0.30, 1.0, 1.0);
gPad->SetBottomMargin(0.001);
gPad->SetLeftMargin(leftMargin);

TH1D* h_xj_pythia8_base = (TH1D*)h_xj_pythia8->Clone("h_xj_pythia8_base");

// colors/markers: central=6•20, +JER=2•24, -JER=4•25
style_line(h_xj_pythia8_base, "x_{J}", "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", 1, 20);
style_line(h_xj_pythia8_posJER, "x_{J}", "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", kMagenta - 2, 20);
style_line(h_xj_pythia8_negJER, "x_{J}", "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", 2, 20);

// Override y-title offset to match bottom pad
h_xj_pythia8_base->GetYaxis()->SetTitleOffset(yTitleOffset);
h_xj_pythia8_posJER->GetYaxis()->SetTitleOffset(yTitleOffset);
h_xj_pythia8_negJER->GetYaxis()->SetTitleOffset(yTitleOffset);

// Ranges
h_xj_pythia8_base->SetAxisRange(0.31, 0.99, "X");
h_xj_pythia8_base->SetAxisRange(0.01, 2.5, "Y");

// Draw (central first sets axes)
h_xj_pythia8_base->Draw("P");
h_xj_pythia8_posJER->Draw("P SAME");
h_xj_pythia8_negJER->Draw("P SAME");

 TLegend* leg_p_h_jer = new TLegend(0.46, 0.76, 0.70, 0.88);
leg_p_h_jer->SetBorderSize(0);
leg_p_h_jer->SetFillStyle(0);
leg_p_h_jer->SetTextSize(0.030);
leg_p_h_jer->SetMargin(0.28);
leg_p_h_jer->AddEntry(h_xj_pythia8_base, "Unfold w/ PYTHIA8", "p");
leg_p_h_jer->AddEntry(h_xj_pythia8_posJER, "Unfold + 5% JER", "p");
 leg_p_h_jer->AddEntry(h_xj_pythia8_negJER, "Unfold - 5% JER", "p");
leg_p_h_jer->Draw();

TLegend* leg_sphenix_xj_jer = nullptr;
xj_functions::DrawsPHENIXLegendNew(leg_sphenix_xj_jer, "leg_sphenix_xj_jer",
    0.15, 0.60, 0.45, 0.88, "p+p #sqrt{s}=200 GeV",
    "anti-#it{k}_{#it{t}} #it{R} = 0.2","20.9 < p_{T,1} < 31.2 GeV", "Matched Truth Dijets",
    true, true, false, true, false);

// ---------- Bottom pad (ratios) ----------
c_xj_JER->cd(2);
gPad->SetPad(0.0, 0.0, 1.0, 0.30);
gPad->SetTopMargin(0.001);
gPad->SetBottomMargin(0.25);
gPad->SetLeftMargin(leftMargin);

// Build ratios: central / (±JER)
TH1D* h_ratio_pos = (TH1D*)h_xj_pythia8_base->Clone("h_ratio_p8_over_posJER");
h_ratio_pos->Divide(h_xj_pythia8_posJER);

TH1D* h_ratio_neg = (TH1D*)h_xj_pythia8_base->Clone("h_ratio_p8_over_negJER");
h_ratio_neg->Divide(h_xj_pythia8_negJER);

// Axis & style for ratios
auto style_ratio = [yTitleOffset](TH1D* r, int col, int mstyle){
  r->SetTitle("");
  r->GetXaxis()->SetTitle("x_{J}");
  r->GetYaxis()->SetTitle("Unfolded/Truth");
  r->SetAxisRange(0.31, 0.99, "X");
  r->SetAxisRange(0.0, 1.99, "Y");
  r->GetXaxis()->SetTitleSize(0.08);
  r->GetXaxis()->SetLabelSize(0.08);
  r->GetYaxis()->SetTitleSize(0.08);
  r->GetYaxis()->SetLabelSize(0.08);
  r->GetYaxis()->SetTitleOffset(0.4); // match top pad
  r->SetMarkerStyle(mstyle);
  r->SetMarkerColor(col);
  r->SetLineColor(col);
  r->SetMarkerSize(1.5);
};

// Match top-pad colors/markers for clarity
style_ratio(h_ratio_pos, kMagenta - 2, 20); // central / +JER
style_ratio(h_ratio_neg, 2, 20); // central / -JER

h_ratio_pos->Draw("PE");
h_ratio_neg->Draw("PE SAME");

// Reference line at 1
TLine *lineJER = new TLine(0.3, 1.0, 0.99, 1.0);
lineJER->SetLineStyle(9);
lineJER->SetLineColor(2);
lineJER->Draw("SAME");




// Conference Note Comparison

//Define arrays for note
 // ===== First dataset: 20-30 =====


std::vector<double> edges = {
  0.300817, 0.34377, 0.392857, 0.448953, 0.513059,
  0.586319, 0.67004, 0.765715, 0.875051, 1.0
}; // size = 10  -> 9 bins


 // ===== First dataset: 20-30 =====
double negJER_20_30[]  = {0.009, 0.086, 0.437, 0.789, 0.989, 1.247, 1.699, 2.000, 1.900};
double Centroid_20_30[] = {0.014, 0.129, 0.581, 0.996, 1.183, 1.398, 1.756, 2.222, 2.237};
double posJER_20_30[]  = {0.027, 0.229, 0.896, 1.434, 1.541, 1.627, 1.806, 2.430, 2.616};

// ===== Second dataset: 30-40 =====
double negJER_30_40[]  = {0.086, 0.151, 0.251, 0.366, 0.509, 0.760, 1.305, 2.208, 2.581};
double Centroid_30_40[] = {0.229, 0.308, 0.423, 0.573, 0.738, 0.975, 1.455, 2.437, 2.925};
double posJER_30_40[]  = {0.530, 0.667, 0.803, 0.925, 1.018, 1.168, 1.591, 2.789, 3.505};

// ===== Third dataset: 40-60 =====
double negJER_40_60[]  = {0.014, 0.057, 0.115, 0.208, 0.308, 0.495, 1.039, 2.330, 2.896};
double Centroid_40_60[] = {0.158, 0.237, 0.330, 0.452, 0.595, 0.789, 1.290, 2.552, 3.269};
double posJER_40_60[]  = {0.416, 0.638, 0.724, 0.717, 0.853, 1.054, 1.548, 3.075, 4.057};


   // Open files
TFile *f_p8_R4     = TFile::Open("UnfoldedxJHists_Play_pythia8_R4.root");
TFile *f_p8pos_R4  = TFile::Open("UnfoldedxJHists_Play_pythia8_posJER_R4.root");
TFile *f_p8neg_R4  = TFile::Open("UnfoldedxJHists_Play_pythia8_negJER_R4.root");

// Grab histograms
// ===== Grab histograms: 30_40 =====
TH1D *h_xj_pythia8_R4_30_40        = (TH1D*)f_p8_R4    ->Get("h_xj_unfold_data_30_40");
TH1D *h_xj_pythia8_posJER_R4_30_40 = (TH1D*)f_p8pos_R4 ->Get("h_xj_unfold_data_30_40");
TH1D *h_xj_pythia8_negJER_R4_30_40 = (TH1D*)f_p8neg_R4 ->Get("h_xj_unfold_data_30_40");

// ===== Grab histograms: 20_30 =====
TH1D *h_xj_pythia8_R4_20_30        = (TH1D*)f_p8_R4    ->Get("h_xj_unfold_data_20_30");
TH1D *h_xj_pythia8_posJER_R4_20_30 = (TH1D*)f_p8pos_R4 ->Get("h_xj_unfold_data_20_30");
TH1D *h_xj_pythia8_negJER_R4_20_30 = (TH1D*)f_p8neg_R4 ->Get("h_xj_unfold_data_20_30");

// ===== Grab histograms: 40_60 =====
TH1D *h_xj_pythia8_R4_40_60        = (TH1D*)f_p8_R4    ->Get("h_xj_unfold_data_40_60");
TH1D *h_xj_pythia8_posJER_R4_40_60 = (TH1D*)f_p8pos_R4 ->Get("h_xj_unfold_data_40_60");
TH1D *h_xj_pythia8_negJER_R4_40_60 = (TH1D*)f_p8neg_R4 ->Get("h_xj_unfold_data_40_60");



// Convert C arrays to std::vector<double>
std::vector<double> yLow_20_30(negJER_20_30,  negJER_20_30  + 9);
std::vector<double> yCent_20_30(Centroid_20_30, Centroid_20_30 + 9);
std::vector<double> yHigh_20_30(posJER_20_30, posJER_20_30 + 9);


 // Convert C arrays to std::vector<double>
std::vector<double> yLow_30_40(negJER_30_40,  negJER_30_40  + 9);
std::vector<double> yCent_30_40(Centroid_30_40, Centroid_30_40 + 9);
std::vector<double> yHigh_30_40(posJER_30_40, posJER_30_40 + 9);


// Convert C arrays to std::vector<double>
std::vector<double> yLow_40_60(negJER_40_60,  negJER_40_60  + 9);
std::vector<double> yCent_40_60(Centroid_40_60, Centroid_40_60 + 9);
std::vector<double> yHigh_40_60(posJER_40_60, posJER_40_60 + 9);

 TLegend* leg_sphenix = nullptr;







// 20–30
auto c20_30 = DrawArraysAndHistsTogether(
  edges, yCent_20_30, yLow_20_30, yHigh_20_30,
  h_xj_pythia8_R4_20_30,        // centroid hist
  h_xj_pythia8_negJER_R4_20_30, // lower hist
  h_xj_pythia8_posJER_R4_20_30, // upper hist
  "Error (Conference Note)", "JER Error (This Analysis)",
  "x_{J}", "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", "20_30");
gPad->SetLeftMargin(leftMargin);
 Make_sPHENIX_Legend_for("20_30", leg_sphenix);  // for the 20–30% canvas
// ... draw or SaveAs


// 30–40 (build yLow_30_40, yCent_30_40, yHigh_30_40 similarly)
auto c30_40 = DrawArraysAndHistsTogether(
  edges, yCent_30_40, yLow_30_40, yHigh_30_40,
  h_xj_pythia8_R4_30_40,
  h_xj_pythia8_negJER_R4_30_40,
  h_xj_pythia8_posJER_R4_30_40,
   "Error (Conference Note)", "JER Error (This Analysis)",
  "x_{J}", "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", "30_40");
gPad->SetLeftMargin(leftMargin);
 Make_sPHENIX_Legend_for("30_40", leg_sphenix);  // for the 30–40% canvas
// ... draw or SaveAs

// 40–60
auto c40_60 = DrawArraysAndHistsTogether(
  edges, yCent_40_60, yLow_40_60, yHigh_40_60,
  h_xj_pythia8_R4_40_60,
  h_xj_pythia8_negJER_R4_40_60,
  h_xj_pythia8_posJER_R4_40_60,
   "Error (Conference Note)", "JER Error (This Analysis)",
  "x_{J}", "#frac{1 dN_{pair}}{N_{pairs} dx_{J}}", "40_60");

gPad->SetLeftMargin(leftMargin);
 
Make_sPHENIX_Legend_for("40_60", leg_sphenix);  // for the 40–60% canvas
// ... draw or SaveAs
}

