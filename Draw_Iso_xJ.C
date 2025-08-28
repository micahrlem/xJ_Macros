#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

#include "xj_functions.h"
#include "read_binning.h"

 //xj histograms
 



void Draw_Iso_xJ(string infile01 = "Hists_Iso_Split_unmatchtruth_R2_pythia8-splice.root",string infile02 = "Hists_Iso_Split_unmatchtruth_R4_pythia8-splice.root", const std::string configfile = "binning_original.config")

{

  read_binning rb(configfile.c_str());
  Int_t read_nbins = rb.get_nbins();
  const int nbins = read_nbins;
  float ipt_bins[nbins+1];
  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
    }

   float truth_leading_cut = rb.get_truth_leading_cut();
  float truth_subleading_cut = rb.get_truth_subleading_cut();

  float reco_leading_cut = rb.get_reco_leading_cut();
  float reco_subleading_cut = rb.get_reco_subleading_cut();

  float measure_leading_cut = rb.get_measure_leading_cut();
  float measure_subleading_cut = rb.get_measure_subleading_cut();

  int truth_leading_bin = rb.get_truth_leading_bin();
  int truth_subleading_bin = rb.get_truth_subleading_bin();

  int reco_leading_bin = rb.get_reco_leading_bin();
  int reco_subleading_bin = rb.get_reco_subleading_bin();

  int measure_leading_bin = rb.get_measure_leading_bin();
  int measure_subleading_bin = rb.get_measure_subleading_bin();

   std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;
  
SetsPhenixStyle();
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
TFile* infile1 = TFile::Open(infile01.c_str());  // R = 0.2
TFile* infile2 = TFile::Open(infile02.c_str());  // R = 0.4

// R2 histograms
TH1D* hTruthPtLead_all_R2 = (TH1D*)infile1->Get("hTruthPtLead_all"); TH1D* hTruthPtSubLead_all_R2 = (TH1D*)infile1->Get("hTruthPtSubLead_all");
TH1D* hTruthPtLead_iso_R2 = (TH1D*)infile1->Get("hTruthPtLead_iso"); TH1D* hTruthPtSubLead_iso_R2 = (TH1D*)infile1->Get("hTruthPtSubLead_iso");
TH1D* hTruthEtaLead_all_R2 = (TH1D*)infile1->Get("hTruthEtaLead_all"); TH1D* hTruthEtaSubLead_all_R2 = (TH1D*)infile1->Get("hTruthEtaSubLead_all");
TH1D* hTruthEtaLead_iso_R2 = (TH1D*)infile1->Get("hTruthEtaLead_iso"); TH1D* hTruthEtaSubLead_iso_R2 = (TH1D*)infile1->Get("hTruthEtaSubLead_iso");
TH1D* hTruthPhiLead_all_R2 = (TH1D*)infile1->Get("hTruthPhiLead_all"); TH1D* hTruthPhiSubLead_all_R2 = (TH1D*)infile1->Get("hTruthPhiSubLead_all");
TH1D* hTruthPhiLead_iso_R2 = (TH1D*)infile1->Get("hTruthPhiLead_iso"); TH1D* hTruthPhiSubLead_iso_R2 = (TH1D*)infile1->Get("hTruthPhiSubLead_iso");
TH1D* hTruthDeltaPhi_all_R2 = (TH1D*)infile1->Get("hTruthDeltaPhi_all"); TH1D* hTruthDeltaEta_all_R2 = (TH1D*)infile1->Get("hTruthDeltaEta_all");
TH1D* hTruthDeltaPhi_iso_R2 = (TH1D*)infile1->Get("hTruthDeltaPhi_iso"); TH1D* hTruthDeltaEta_iso_R2 = (TH1D*)infile1->Get("hTruthDeltaEta_iso");
TH2D* hTruthPt_PhiLead_all_R2 = (TH2D*)infile1->Get("hTruthPt_PhiLead_all"); TH2D* hTruthPt_PhiSubLead_all_R2 = (TH2D*)infile1->Get("hTruthPt_PhiSubLead_all");
TH2D* hTruthPt_PhiLead_iso_R2 = (TH2D*)infile1->Get("hTruthPt_PhiLead_iso"); TH2D* hTruthPt_PhiSubLead_iso_R2 = (TH2D*)infile1->Get("hTruthPt_PhiSubLead_iso");
TH2D* hTruthPt_EtaLead_all_R2 = (TH2D*)infile1->Get("hTruthPt_EtaLead_all"); TH2D* hTruthPt_EtaSubLead_all_R2 = (TH2D*)infile1->Get("hTruthPt_EtaSubLead_all");
TH2D* hTruthPt_EtaLead_iso_R2 = (TH2D*)infile1->Get("hTruthPt_EtaLead_iso"); TH2D* hTruthPt_EtaSubLead_iso_R2 = (TH2D*)infile1->Get("hTruthPt_EtaSubLead_iso");
TH2D* h_pt1pt2_truth_all_R2 = (TH2D*)infile1->Get("h_pt1pt2_truth_all"); TH2D* h_pt1pt2_truth_iso_R2 = (TH2D*)infile1->Get("h_pt1pt2_truth_iso");
TH1D* h_xj_classical_truth_all_R2 = (TH1D*)infile1->Get("h_xj_classical_truth_all"); TH1D* h_xj_classical_truth_iso_R2 = (TH1D*)infile1->Get("h_xj_classical_truth_iso");

// R4 histograms
TH1D* hTruthPtLead_all_R4 = (TH1D*)infile2->Get("hTruthPtLead_all"); TH1D* hTruthPtSubLead_all_R4 = (TH1D*)infile2->Get("hTruthPtSubLead_all");
TH1D* hTruthPtLead_iso_R4 = (TH1D*)infile2->Get("hTruthPtLead_iso"); TH1D* hTruthPtSubLead_iso_R4 = (TH1D*)infile2->Get("hTruthPtSubLead_iso");
TH1D* hTruthEtaLead_all_R4 = (TH1D*)infile2->Get("hTruthEtaLead_all"); TH1D* hTruthEtaSubLead_all_R4 = (TH1D*)infile2->Get("hTruthEtaSubLead_all");
TH1D* hTruthEtaLead_iso_R4 = (TH1D*)infile2->Get("hTruthEtaLead_iso"); TH1D* hTruthEtaSubLead_iso_R4 = (TH1D*)infile2->Get("hTruthEtaSubLead_iso");
TH1D* hTruthPhiLead_all_R4 = (TH1D*)infile2->Get("hTruthPhiLead_all"); TH1D* hTruthPhiSubLead_all_R4 = (TH1D*)infile2->Get("hTruthPhiSubLead_all");
TH1D* hTruthPhiLead_iso_R4 = (TH1D*)infile2->Get("hTruthPhiLead_iso"); TH1D* hTruthPhiSubLead_iso_R4 = (TH1D*)infile2->Get("hTruthPhiSubLead_iso");
TH1D* hTruthDeltaPhi_all_R4 = (TH1D*)infile2->Get("hTruthDeltaPhi_all"); TH1D* hTruthDeltaEta_all_R4 = (TH1D*)infile2->Get("hTruthDeltaEta_all");
TH1D* hTruthDeltaPhi_iso_R4 = (TH1D*)infile2->Get("hTruthDeltaPhi_iso"); TH1D* hTruthDeltaEta_iso_R4 = (TH1D*)infile2->Get("hTruthDeltaEta_iso");
TH2D* hTruthPt_PhiLead_all_R4 = (TH2D*)infile2->Get("hTruthPt_PhiLead_all"); TH2D* hTruthPt_PhiSubLead_all_R4 = (TH2D*)infile2->Get("hTruthPt_PhiSubLead_all");
TH2D* hTruthPt_PhiLead_iso_R4 = (TH2D*)infile2->Get("hTruthPt_PhiLead_iso"); TH2D* hTruthPt_PhiSubLead_iso_R4 = (TH2D*)infile2->Get("hTruthPt_PhiSubLead_iso");
TH2D* hTruthPt_EtaLead_all_R4 = (TH2D*)infile2->Get("hTruthPt_EtaLead_all"); TH2D* hTruthPt_EtaSubLead_all_R4 = (TH2D*)infile2->Get("hTruthPt_EtaSubLead_all");
TH2D* hTruthPt_EtaLead_iso_R4 = (TH2D*)infile2->Get("hTruthPt_EtaLead_iso"); TH2D* hTruthPt_EtaSubLead_iso_R4 = (TH2D*)infile2->Get("hTruthPt_EtaSubLead_iso");
TH2D* h_pt1pt2_truth_all_R4 = (TH2D*)infile2->Get("h_pt1pt2_truth_all"); TH2D* h_pt1pt2_truth_iso_R4 = (TH2D*)infile2->Get("h_pt1pt2_truth_iso");
TH1D* h_xj_classical_truth_all_R4 = (TH1D*)infile2->Get("h_xj_classical_truth_all"); TH1D* h_xj_classical_truth_iso_R4 = (TH1D*)infile2->Get("h_xj_classical_truth_iso");

 // Reco R2 histograms
TH1D* hRecoPtLead_all_R2 = (TH1D*)infile1->Get("hRecoPtLead_all"); TH1D* hRecoPtSubLead_all_R2 = (TH1D*)infile1->Get("hRecoPtSubLead_all");
TH1D* hRecoPtLead_iso_R2 = (TH1D*)infile1->Get("hRecoPtLead_iso"); TH1D* hRecoPtSubLead_iso_R2 = (TH1D*)infile1->Get("hRecoPtSubLead_iso");
TH1D* hRecoEtaLead_all_R2 = (TH1D*)infile1->Get("hRecoEtaLead_all"); TH1D* hRecoEtaSubLead_all_R2 = (TH1D*)infile1->Get("hRecoEtaSubLead_all");
TH1D* hRecoEtaLead_iso_R2 = (TH1D*)infile1->Get("hRecoEtaLead_iso"); TH1D* hRecoEtaSubLead_iso_R2 = (TH1D*)infile1->Get("hRecoEtaSubLead_iso");
TH1D* hRecoPhiLead_all_R2 = (TH1D*)infile1->Get("hRecoPhiLead_all"); TH1D* hRecoPhiSubLead_all_R2 = (TH1D*)infile1->Get("hRecoPhiSubLead_all");
TH1D* hRecoPhiLead_iso_R2 = (TH1D*)infile1->Get("hRecoPhiLead_iso"); TH1D* hRecoPhiSubLead_iso_R2 = (TH1D*)infile1->Get("hRecoPhiSubLead_iso");
TH1D* hRecoDeltaPhi_all_R2 = (TH1D*)infile1->Get("hRecoDeltaPhi_all"); TH1D* hRecoDeltaEta_all_R2 = (TH1D*)infile1->Get("hRecoDeltaEta_all");
TH1D* hRecoDeltaPhi_iso_R2 = (TH1D*)infile1->Get("hRecoDeltaPhi_iso"); TH1D* hRecoDeltaEta_iso_R2 = (TH1D*)infile1->Get("hRecoDeltaEta_iso");
TH2D* hRecoPt_PhiLead_all_R2 = (TH2D*)infile1->Get("hRecoPt_PhiLead_all"); TH2D* hRecoPt_PhiSubLead_all_R2 = (TH2D*)infile1->Get("hRecoPt_PhiSubLead_all");
TH2D* hRecoPt_PhiLead_iso_R2 = (TH2D*)infile1->Get("hRecoPt_PhiLead_iso"); TH2D* hRecoPt_PhiSubLead_iso_R2 = (TH2D*)infile1->Get("hRecoPt_PhiSubLead_iso");
TH2D* hRecoPt_EtaLead_all_R2 = (TH2D*)infile1->Get("hRecoPt_EtaLead_all"); TH2D* hRecoPt_EtaSubLead_all_R2 = (TH2D*)infile1->Get("hRecoPt_EtaSubLead_all");
TH2D* hRecoPt_EtaLead_iso_R2 = (TH2D*)infile1->Get("hRecoPt_EtaLead_iso"); TH2D* hRecoPt_EtaSubLead_iso_R2 = (TH2D*)infile1->Get("hRecoPt_EtaSubLead_iso");
TH2D* h_pt1pt2_reco_all_R2 = (TH2D*)infile1->Get("h_pt1pt2_reco_all"); TH2D* h_pt1pt2_reco_iso_R2 = (TH2D*)infile1->Get("h_pt1pt2_reco_iso");
TH1D* h_xj_classical_reco_all_R2 = (TH1D*)infile1->Get("h_xj_classical_reco_all"); TH1D* h_xj_classical_reco_iso_R2 = (TH1D*)infile1->Get("h_xj_classical_reco_iso");

// Reco R4 histograms
TH1D* hRecoPtLead_all_R4 = (TH1D*)infile2->Get("hRecoPtLead_all"); TH1D* hRecoPtSubLead_all_R4 = (TH1D*)infile2->Get("hRecoPtSubLead_all");
TH1D* hRecoPtLead_iso_R4 = (TH1D*)infile2->Get("hRecoPtLead_iso"); TH1D* hRecoPtSubLead_iso_R4 = (TH1D*)infile2->Get("hRecoPtSubLead_iso");
TH1D* hRecoEtaLead_all_R4 = (TH1D*)infile2->Get("hRecoEtaLead_all"); TH1D* hRecoEtaSubLead_all_R4 = (TH1D*)infile2->Get("hRecoEtaSubLead_all");
TH1D* hRecoEtaLead_iso_R4 = (TH1D*)infile2->Get("hRecoEtaLead_iso"); TH1D* hRecoEtaSubLead_iso_R4 = (TH1D*)infile2->Get("hRecoEtaSubLead_iso");
TH1D* hRecoPhiLead_all_R4 = (TH1D*)infile2->Get("hRecoPhiLead_all"); TH1D* hRecoPhiSubLead_all_R4 = (TH1D*)infile2->Get("hRecoPhiSubLead_all");
TH1D* hRecoPhiLead_iso_R4 = (TH1D*)infile2->Get("hRecoPhiLead_iso"); TH1D* hRecoPhiSubLead_iso_R4 = (TH1D*)infile2->Get("hRecoPhiSubLead_iso");
TH1D* hRecoDeltaPhi_all_R4 = (TH1D*)infile2->Get("hRecoDeltaPhi_all"); TH1D* hRecoDeltaEta_all_R4 = (TH1D*)infile2->Get("hRecoDeltaEta_all");
TH1D* hRecoDeltaPhi_iso_R4 = (TH1D*)infile2->Get("hRecoDeltaPhi_iso"); TH1D* hRecoDeltaEta_iso_R4 = (TH1D*)infile2->Get("hRecoDeltaEta_iso");
TH2D* hRecoPt_PhiLead_all_R4 = (TH2D*)infile2->Get("hRecoPt_PhiLead_all"); TH2D* hRecoPt_PhiSubLead_all_R4 = (TH2D*)infile2->Get("hRecoPt_PhiSubLead_all");
TH2D* hRecoPt_PhiLead_iso_R4 = (TH2D*)infile2->Get("hRecoPt_PhiLead_iso"); TH2D* hRecoPt_PhiSubLead_iso_R4 = (TH2D*)infile2->Get("hRecoPt_PhiSubLead_iso");
TH2D* hRecoPt_EtaLead_all_R4 = (TH2D*)infile2->Get("hRecoPt_EtaLead_all"); TH2D* hRecoPt_EtaSubLead_all_R4 = (TH2D*)infile2->Get("hRecoPt_EtaSubLead_all");
TH2D* hRecoPt_EtaLead_iso_R4 = (TH2D*)infile2->Get("hRecoPt_EtaLead_iso"); TH2D* hRecoPt_EtaSubLead_iso_R4 = (TH2D*)infile2->Get("hRecoPt_EtaSubLead_iso");
TH2D* h_pt1pt2_reco_all_R4 = (TH2D*)infile2->Get("h_pt1pt2_reco_all"); TH2D* h_pt1pt2_reco_iso_R4 = (TH2D*)infile2->Get("h_pt1pt2_reco_iso");
TH1D* h_xj_classical_reco_all_R4 = (TH1D*)infile2->Get("h_xj_classical_reco_all"); TH1D* h_xj_classical_reco_iso_R4 = (TH1D*)infile2->Get("h_xj_classical_reco_iso");

 

// Plotting Xj Iso vs All, R = 0.2 vs R = 0.4
TCanvas* can_xj_comparison_truth = new TCanvas("can_xj_comparison_truth", "Xj Classical vs Projected", 700, 800);
can_xj_comparison_truth->cd(1);
gPad->SetLogy(0); // Remove log scale
 TH1D *h_xj_projected_truth_all_R2 = new TH1D("h_xj_projected_truth_all_R2", ";x_{J};",nbins, ixj_bins);
 TH1D *h_xj_projected_truth_iso_R2 = new TH1D("h_xj_projected_truth_iso_R2", ";x_{J};",nbins, ixj_bins);
 TH1D *h_xj_projected_truth_all_R4 = new TH1D("h_xj_projected_truth_all_R4", ";x_{J};",nbins, ixj_bins);
 TH1D *h_xj_projected_truth_iso_R4 = new TH1D("h_xj_projected_truth_iso_R4", ";x_{J};",nbins, ixj_bins);

  xj_functions::project_xj(h_pt1pt2_truth_all_R2, h_xj_projected_truth_all_R2, nbins, truth_leading_bin, nbins -1, truth_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_truth_iso_R2, h_xj_projected_truth_iso_R2, nbins, truth_leading_bin, nbins -1, truth_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_truth_all_R4, h_xj_projected_truth_all_R4, nbins, truth_leading_bin, nbins - 1, truth_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_truth_iso_R4, h_xj_projected_truth_iso_R4, nbins, truth_leading_bin, nbins - 1, truth_subleading_bin, nbins - 1);
// Normalize histograms using the function with 19 bins
 xj_functions::normalize_histo(h_xj_projected_truth_iso_R2, 19);
 xj_functions::normalize_histo(h_xj_projected_truth_all_R2, 19);
 xj_functions::normalize_histo(h_xj_projected_truth_iso_R4, 19);
 xj_functions::normalize_histo(h_xj_projected_truth_all_R4, 19);


 xj_functions::plot_xj_and_ratio_compare(can_xj_comparison_truth, h_xj_projected_truth_iso_R4, h_xj_projected_truth_iso_R2, h_xj_projected_truth_iso_R4, h_xj_projected_truth_all_R4, h_xj_projected_truth_iso_R2, h_xj_projected_truth_all_R2, nullptr, false, "x_{J}",  "Iso/All",0, 1.0, 1,  53, 2,  20, 1, 53, 2, 20, 4,  54, kMagenta - 2,  21, 38,  26);



// Legend for Classical vs Projected
TLegend *leg_xj = new TLegend(0.65, 0.45, 0.87, 0.55);
leg_xj->SetBorderSize(0);
leg_xj->SetTextSize(0.030);  // Set text size to 0.026
leg_xj->AddEntry(h_xj_projected_truth_iso_R4, "Xj Iso R = 0.4", "P");
leg_xj->AddEntry(h_xj_projected_truth_all_R4, "Xj All R = 0.4", "P");
 leg_xj->AddEntry(h_xj_projected_truth_iso_R2, "Xj Iso R = 0.2", "P");
leg_xj->AddEntry(h_xj_projected_truth_all_R2, "Xj All R = 0.2", "P");
leg_xj->Draw();

// General Legend

  TLegend *sleg_xj_truth; 
 xj_functions::DrawsPHENIXLegendNew(sleg_xj_truth, "slegxjtruth",0.45,0.68,0.75,0.93, "p+p #sqrt{s}=200 GeV", "anti-#it{k}_{#it{t}} #it{R} = 0.2, 0.4", "16 < p_{T,1} < 60.77 GeV", "Matched Truth Dijets",true,false, false, true, true);


 TCanvas* can_xj_comparison_reco = new TCanvas("can_xj_comparison_reco", "Xj Classical vs Projected", 700, 800);
can_xj_comparison_reco->cd(1);
gPad->SetLogy(0); // Remove log scale
 TH1D *h_xj_projected_reco_all_R2 = new TH1D("h_xj_projected_reco_all_R2", ";x_{J};",nbins, ixj_bins);
 TH1D *h_xj_projected_reco_iso_R2 = new TH1D("h_xj_projected_reco_iso_R2", ";x_{J};",nbins, ixj_bins);
 TH1D *h_xj_projected_reco_all_R4 = new TH1D("h_xj_projected_reco_all_R4", ";x_{J};",nbins, ixj_bins);
 TH1D *h_xj_projected_reco_iso_R4 = new TH1D("h_xj_projected_reco_iso_R4", ";x_{J};",nbins, ixj_bins);

  xj_functions::project_xj(h_pt1pt2_reco_all_R2, h_xj_projected_reco_all_R2, nbins, reco_leading_bin, nbins -1, reco_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_reco_iso_R2, h_xj_projected_reco_iso_R2, nbins, reco_leading_bin, nbins -1, reco_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_reco_all_R4, h_xj_projected_reco_all_R4, nbins, reco_leading_bin, nbins - 1, reco_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_reco_iso_R4, h_xj_projected_reco_iso_R4, nbins, reco_leading_bin, nbins - 1, reco_subleading_bin, nbins - 1);
// Normalize histograms using the function with 19 bins
 xj_functions::normalize_histo(h_xj_projected_reco_iso_R2, 19);
 xj_functions::normalize_histo(h_xj_projected_reco_all_R2, 19);
 xj_functions::normalize_histo(h_xj_projected_reco_iso_R4, 19);
 xj_functions::normalize_histo(h_xj_projected_reco_all_R4, 19);


 xj_functions::plot_xj_and_ratio_compare(can_xj_comparison_reco, h_xj_projected_reco_iso_R4, h_xj_projected_reco_iso_R2, h_xj_projected_reco_iso_R4, h_xj_projected_reco_all_R4, h_xj_projected_reco_iso_R2, h_xj_projected_reco_all_R2, nullptr, false, "x_{J}",  "Iso/All",0, 1.0, 1,  53, 2,  20, 1, 53, 2, 20, 4,  54, kMagenta-2,  21, 38,  26);



// Legend for Classical vs Projected
TLegend *leg_xj_reco = new TLegend(0.65, 0.45, 0.87, 0.55);
leg_xj_reco->SetBorderSize(0);
leg_xj_reco->SetTextSize(0.030);  // Set text size to 0.026
leg_xj_reco->AddEntry(h_xj_projected_reco_iso_R4, "Xj Iso R = 0.4", "P");
leg_xj_reco->AddEntry(h_xj_projected_reco_all_R4, "Xj All R = 0.4", "P");
 leg_xj_reco->AddEntry(h_xj_projected_reco_iso_R2, "Xj Iso R = 0.2", "P");
leg_xj_reco->AddEntry(h_xj_projected_reco_all_R2, "Xj All R = 0.2", "P");
leg_xj_reco->Draw();

// General Legend
 TLegend *sleg_xj_reco; 
 xj_functions::DrawsPHENIXLegendNew(sleg_xj_reco, "slegxjreco",0.45,0.68,0.75,0.93, "p+p #sqrt{s}=200 GeV", "anti-#it{k}_{#it{t}} #it{R} = 0.2, 0.4", "20.9 < p_{T,1} < 60.77 GeV", "Matched Reco Dijets",false, true, false, true, true);


//Plotting Pt Together
TCanvas* can_lead_sub_pt_truth = new TCanvas("can_lead_sub_pt_truth","",1600,800);
 can_lead_sub_pt_truth->Divide(2,1);
can_lead_sub_pt_truth->cd(1);gPad->SetLogy();
//Truth UnMatched
xj_functions::SetHist(hTruthPtLead_all_R4,"Leading Truth p_{T} [Gev]", "Integrated Counts",2,20,3,1.5);
hTruthPtLead_all_R4->Draw("P");
 hTruthPtLead_all_R4->SetAxisRange(10,64); hTruthPtLead_all_R4->Scale(1/hTruthPtLead_all_R4->Integral());
xj_functions::SetHist(hTruthPtLead_iso_R4,"p_{T} [Gev]", "Counts",2,24,3,1.5);
  hTruthPtLead_iso_R4->SetAxisRange(10,54); hTruthPtLead_iso_R4->Scale(1/hTruthPtLead_iso_R4->Integral());hTruthPtLead_iso_R4->Draw("same, P");
  xj_functions::SetHist(hTruthPtLead_iso_R2,"p_{T} [Gev]", "Counts",1,24,3,1.5);
  hTruthPtLead_iso_R2->SetAxisRange(10,54); hTruthPtLead_iso_R2->Scale(1/hTruthPtLead_iso_R2->Integral());hTruthPtLead_iso_R2->Draw("same, P");
  xj_functions::SetHist(hTruthPtLead_all_R2,"p_{T} [Gev]", "Counts",1,20,3,1.5);
  hTruthPtLead_all_R2->SetAxisRange(10,54); hTruthPtLead_all_R2->Scale(1/hTruthPtLead_all_R2->Integral());hTruthPtLead_all_R2->Draw("same, P");
 TLegend *legleadpt = new TLegend (0.65, 0.70, 0.87,0.80);
 legleadpt->SetTextSize(0.026);
 legleadpt->AddEntry(hTruthPtLead_all_R4, "R = 0.4 All", "P");
 legleadpt->AddEntry(hTruthPtLead_iso_R4, "R = 0.4 Isolated", "P");
 legleadpt->AddEntry(hTruthPtLead_all_R2, "R = 0.2 All", "P");
 legleadpt->AddEntry(hTruthPtLead_iso_R2, "R = 0.2 Isolated", "P");
 legleadpt->SetBorderSize(0);
 legleadpt->Draw();

TLegend *sleg_pt_lead;
xj_functions::DrawsPHENIXLegendNew(sleg_pt_lead, "sleg_pt_lead",0.15,0.20,0.65,0.50, "p+p #sqrt{s}=200 GeV", "anti-#it{k}_{#it{t}} #it{R} = 0.2, 0.4", "20.9 < p_{T,1} < 60.77 GeV", "Matched Truth Leading Dijets",true, false, false, false, true);

 can_lead_sub_pt_truth->cd(2);gPad->SetLogy();
//Truth UnMatched
xj_functions::SetHist(hTruthPtSubLead_all_R4,"Subleading Truth p_{T} [Gev]", "Integrated Counts",2,20,3,1.5);
hTruthPtSubLead_all_R4->Draw("P");hTruthPtSubLead_all_R4->SetAxisRange(10,64); hTruthPtSubLead_all_R4->Scale(1/hTruthPtSubLead_all_R4->Integral());
xj_functions::SetHist(hTruthPtSubLead_iso_R4,"p_{T} [Gev]", "Counts",2,24,3,1.5);
  hTruthPtSubLead_iso_R4->SetAxisRange(10,54); hTruthPtSubLead_iso_R4->Scale(1/hTruthPtSubLead_iso_R4->Integral());hTruthPtSubLead_iso_R4->Draw("same, P");
  xj_functions::SetHist(hTruthPtSubLead_iso_R2,"p_{T} [Gev]", "Counts",1,24,3,1.5);
  hTruthPtSubLead_iso_R2->SetAxisRange(10,54); hTruthPtSubLead_iso_R2->Scale(1/hTruthPtSubLead_iso_R2->Integral());hTruthPtSubLead_iso_R2->Draw("same, P");
  xj_functions::SetHist(hTruthPtSubLead_all_R2,"p_{T} [Gev]", "Counts",1,20,3,1.5);
  hTruthPtSubLead_all_R2->SetAxisRange(10,54); hTruthPtSubLead_all_R2->Scale(1/hTruthPtSubLead_all_R2->Integral());hTruthPtSubLead_all_R2->Draw("same, P");
 TLegend *legsubleadpt = new TLegend (0.65, 0.70, 0.87,0.80);
 legsubleadpt->SetTextSize(0.026);
 legsubleadpt->AddEntry(hTruthPtSubLead_all_R4, "R = 0.4 All", "P");
 legsubleadpt->AddEntry(hTruthPtSubLead_iso_R4, "R = 0.4 Isolated", "P");
 legsubleadpt->AddEntry(hTruthPtSubLead_all_R2, "R = 0.2 All", "P");
 legsubleadpt->AddEntry(hTruthPtSubLead_iso_R2, "R = 0.2 Isolated", "P");
 legsubleadpt->SetBorderSize(0);
 legsubleadpt->Draw();

TLegend *sleg_pt_sublead;
xj_functions::DrawsPHENIXLegendNew(sleg_pt_sublead, "sleg_pt_sublead",0.15,0.20,0.65,0.50, "p+p #sqrt{s}=200 GeV", "anti-#it{k}_{#it{t}} #it{R} = 0.2, 0.4", "20.9 < p_{T,1} < 60.77 GeV", "Matched Truth Subleading Dijets",true, false, false, false, true);

 TCanvas* can_lead_sub_pt_reco = new TCanvas("can_lead_sub_pt_reco","",1600,800);
 can_lead_sub_pt_reco->Divide(2,1);
can_lead_sub_pt_reco->cd(1);gPad->SetLogy();
//Reco UnMatched
xj_functions::SetHist(hRecoPtLead_all_R4,"Leading Reco p_{T} [Gev]", "Integrated Counts",2,20,3,1.5);
hRecoPtLead_all_R4->Draw("P");
 hRecoPtLead_all_R4->SetAxisRange(10,64); hRecoPtLead_all_R4->Scale(1/hRecoPtLead_all_R4->Integral());
xj_functions::SetHist(hRecoPtLead_iso_R4,"p_{T} [Gev]", "Counts",2,24,3,1.5);
  hRecoPtLead_iso_R4->SetAxisRange(10,54); hRecoPtLead_iso_R4->Scale(1/hRecoPtLead_iso_R4->Integral());hRecoPtLead_iso_R4->Draw("same, P");
  xj_functions::SetHist(hRecoPtLead_iso_R2,"p_{T} [Gev]", "Counts",1,24,3,1.5);
  hRecoPtLead_iso_R2->SetAxisRange(10,54); hRecoPtLead_iso_R2->Scale(1/hRecoPtLead_iso_R2->Integral());hRecoPtLead_iso_R2->Draw("same, P");
  xj_functions::SetHist(hRecoPtLead_all_R2,"p_{T} [Gev]", "Counts",1,20,3,1.5);
  hRecoPtLead_all_R2->SetAxisRange(10,54); hRecoPtLead_all_R2->Scale(1/hRecoPtLead_all_R2->Integral());hRecoPtLead_all_R2->Draw("same, P");
 TLegend *legleadptreco = new TLegend (0.65, 0.70, 0.87,0.80);
 legleadptreco->SetTextSize(0.026);
 legleadptreco->AddEntry(hRecoPtLead_all_R4, "R = 0.4 All", "P");
 legleadptreco->AddEntry(hRecoPtLead_iso_R4, "R = 0.4 Isolated", "P");
 legleadptreco->AddEntry(hRecoPtLead_all_R2, "R = 0.2 All", "P");
 legleadptreco->AddEntry(hRecoPtLead_iso_R2, "R = 0.2 Isolated", "P");
 legleadptreco->SetBorderSize(0);
 legleadptreco->Draw();
TLegend *sleg_pt_lead_reco;
xj_functions::DrawsPHENIXLegendNew(sleg_pt_lead_reco, "sleg_pt_lead_reco",0.15,0.20,0.65,0.50, "p+p #sqrt{s}=200 GeV", "anti-#it{k}_{#it{t}} #it{R} = 0.2, 0.4", "20.9 < p_{T,1} < 60.77 GeV", "Matched Reco Leading Dijets",false, true, false, false, true);

 can_lead_sub_pt_reco->cd(2);gPad->SetLogy();
//Reco UnMatched
xj_functions::SetHist(hRecoPtSubLead_all_R4,"Subleading Reco p_{T} [Gev]", "Integrated Counts",2,20,3,1.5);
hRecoPtSubLead_all_R4->Draw("P");hRecoPtSubLead_all_R4->SetAxisRange(10,64); hRecoPtSubLead_all_R4->Scale(1/hRecoPtSubLead_all_R4->Integral());
xj_functions::SetHist(hRecoPtSubLead_iso_R4,"p_{T} [Gev]", "Counts",2,24,3,1.5);
  hRecoPtSubLead_iso_R4->SetAxisRange(10,54); hRecoPtSubLead_iso_R4->Scale(1/hRecoPtSubLead_iso_R4->Integral());hRecoPtSubLead_iso_R4->Draw("same, P");
  xj_functions::SetHist(hRecoPtSubLead_iso_R2,"p_{T} [Gev]", "Counts",1,24,3,1.5);
  hRecoPtSubLead_iso_R2->SetAxisRange(10,54); hRecoPtSubLead_iso_R2->Scale(1/hRecoPtSubLead_iso_R2->Integral());hRecoPtSubLead_iso_R2->Draw("same, P");
  xj_functions::SetHist(hRecoPtSubLead_all_R2,"p_{T} [Gev]", "Counts",1,20,3,1.5);
  hRecoPtSubLead_all_R2->SetAxisRange(10,54); hRecoPtSubLead_all_R2->Scale(1/hRecoPtSubLead_all_R2->Integral());hRecoPtSubLead_all_R2->Draw("same, P");
 TLegend *legsubleadptreco = new TLegend (0.65, 0.70, 0.87,0.80);
 legsubleadptreco->SetTextSize(0.026);
 legsubleadptreco->AddEntry(hRecoPtSubLead_all_R4, "R = 0.4 All", "P");
 legsubleadptreco->AddEntry(hRecoPtSubLead_iso_R4, "R = 0.4 Isolated", "P");
 legsubleadptreco->AddEntry(hRecoPtSubLead_all_R2, "R = 0.2 All", "P");
 legsubleadptreco->AddEntry(hRecoPtSubLead_iso_R2, "R = 0.2 Isolated", "P");
 legsubleadptreco->SetBorderSize(0);
 legsubleadptreco->Draw();
TLegend *sleg_pt_sublead_reco;
xj_functions::DrawsPHENIXLegendNew(sleg_pt_sublead_reco, "sleg_pt_sublead_reco",0.15,0.20,0.65,0.50, "p+p #sqrt{s}=200 GeV", "anti-#it{k}_{#it{t}} #it{R} = 0.2, 0.4", "20.9 < p_{T,1} < 60.77 GeV", "Matched Reco Subleading Dijets",false, true, false, false, true);
 // === R = 0.4 Eta Plots ===
TCanvas* can_eta_R4 = new TCanvas("can_eta_R4","",1600,800);
can_eta_R4->Divide(2,1);
can_eta_R4->cd(1); gPad->SetLogy();
// Truth + Reco Lead, R = 0.4
xj_functions::SetHist(hTruthEtaLead_all_R4,"Leading Truth #eta", "Integrated Counts",1,20,3,1.5);
hTruthEtaLead_all_R4->SetAxisRange(-0.7, 0.7); hTruthEtaLead_all_R4->Scale(1/hTruthEtaLead_all_R4->Integral()); hTruthEtaLead_all_R4->Draw("P");

xj_functions::SetHist(hTruthEtaLead_iso_R4,"#eta", "Counts",1,24,3,1.5);
hTruthEtaLead_iso_R4->SetAxisRange(-0.7, 0.7); hTruthEtaLead_iso_R4->Scale(1/hTruthEtaLead_iso_R4->Integral()); hTruthEtaLead_iso_R4->Draw("same, P");

xj_functions::SetHist(hRecoEtaLead_all_R4,"#eta", "Counts",2,20,3,1.5);
hRecoEtaLead_all_R4->SetAxisRange(-0.7, 0.7); hRecoEtaLead_all_R4->Scale(1/hRecoEtaLead_all_R4->Integral()); hRecoEtaLead_all_R4->Draw("same, P");

xj_functions::SetHist(hRecoEtaLead_iso_R4,"#eta", "Counts",2,24,3,1.5);
hRecoEtaLead_iso_R4->SetAxisRange(-0.7, 0.7); hRecoEtaLead_iso_R4->Scale(1/hRecoEtaLead_iso_R4->Integral()); hRecoEtaLead_iso_R4->Draw("same, P");

TLegend *leg_eta_lead_R4 = new TLegend (0.25, 0.80, 0.75,0.90);
leg_eta_lead_R4->SetTextSize(0.026);
leg_eta_lead_R4->AddEntry(hTruthEtaLead_all_R4, "Truth R = 0.4 All", "P");
leg_eta_lead_R4->AddEntry(hTruthEtaLead_iso_R4, "Truth R = 0.4 Isolated", "P");
leg_eta_lead_R4->AddEntry(hRecoEtaLead_all_R4, "Reco R = 0.4 All", "P");
leg_eta_lead_R4->AddEntry(hRecoEtaLead_iso_R4, "Reco R = 0.4 Isolated", "P");
leg_eta_lead_R4->SetBorderSize(0);
leg_eta_lead_R4->Draw();

TLegend *sleg_eta_lead_R4;
 xj_functions::DrawsPHENIXLegendNew(sleg_eta_lead_R4, "sleg_eta_lead_R4",0.25,0.25,0.75,0.55, "p+p #sqrt{s}=200 GeV", "anti-#it{k}_{#it{t}} #it{R} = 0.4", "20.9 < p_{T,1} < 60.77 GeV", "Matched Truth and Reco Leading Dijets",true, true, false, false, true);

can_eta_R4->cd(2); gPad->SetLogy();
// Truth + Reco Sublead, R = 0.4
xj_functions::SetHist(hTruthEtaSubLead_all_R4,"Subleading Truth #eta", "Integrated Counts",1,20,3,1.5);
hTruthEtaSubLead_all_R4->SetAxisRange(-0.7, 0.7); hTruthEtaSubLead_all_R4->Scale(1/hTruthEtaSubLead_all_R4->Integral()); hTruthEtaSubLead_all_R4->Draw("P");

xj_functions::SetHist(hTruthEtaSubLead_iso_R4,"#eta", "Counts",1,24,3,1.5);
hTruthEtaSubLead_iso_R4->SetAxisRange(-0.7, 0.7); hTruthEtaSubLead_iso_R4->Scale(1/hTruthEtaSubLead_iso_R4->Integral()); hTruthEtaSubLead_iso_R4->Draw("same, P");

xj_functions::SetHist(hRecoEtaSubLead_all_R4,"#eta", "Counts",2,20,3,1.5);
hRecoEtaSubLead_all_R4->SetAxisRange(-0.7, 0.7); hRecoEtaSubLead_all_R4->Scale(1/hRecoEtaSubLead_all_R4->Integral()); hRecoEtaSubLead_all_R4->Draw("same, P");

xj_functions::SetHist(hRecoEtaSubLead_iso_R4,"#eta", "Counts",2,24,3,1.5);
hRecoEtaSubLead_iso_R4->SetAxisRange(-0.7, 0.7); hRecoEtaSubLead_iso_R4->Scale(1/hRecoEtaSubLead_iso_R4->Integral()); hRecoEtaSubLead_iso_R4->Draw("same, P");

TLegend *leg_eta_sublead_R4 = new TLegend (0.25, 0.80, 0.75,0.90);
leg_eta_sublead_R4->SetTextSize(0.026);
leg_eta_sublead_R4->AddEntry(hTruthEtaSubLead_all_R4, "Truth R = 0.4 All", "P");
leg_eta_sublead_R4->AddEntry(hTruthEtaSubLead_iso_R4, "Truth R = 0.4 Isolated", "P");
leg_eta_sublead_R4->AddEntry(hRecoEtaSubLead_all_R4, "Reco R = 0.4 All", "P");
leg_eta_sublead_R4->AddEntry(hRecoEtaSubLead_iso_R4, "Reco R = 0.4 Isolated", "P");
leg_eta_sublead_R4->SetBorderSize(0);
leg_eta_sublead_R4->Draw();

TLegend *sleg_eta_sublead_R4;
xj_functions::DrawsPHENIXLegendNew(sleg_eta_sublead_R4, "sleg_eta_sublead_R4",0.25,0.25,0.75,0.55, "p+p #sqrt{s}=200 GeV", "anti-#it{k}_{#it{t}} #it{R} = 0.4", "20.9 < p_{T,1} < 60.77 GeV", "Matched Truth and Reco Subleading Dijets",true, true, false, false, true);


// === R = 0.2 Eta Plots ===
TCanvas* can_eta_R2 = new TCanvas("can_eta_R2","",1600,800);
can_eta_R2->Divide(2,1);
can_eta_R2->cd(1); gPad->SetLogy();
// Truth + Reco Lead, R = 0.2
xj_functions::SetHist(hTruthEtaLead_all_R2,"Leading Truth #eta", "Integrated Counts",1,20,3,1.5);
hTruthEtaLead_all_R2->SetAxisRange(-0.9, 0.9); hTruthEtaLead_all_R2->Scale(1/hTruthEtaLead_all_R2->Integral()); hTruthEtaLead_all_R2->Draw("P");

xj_functions::SetHist(hTruthEtaLead_iso_R2,"#eta", "Counts",1,24,3,1.5);
hTruthEtaLead_iso_R2->SetAxisRange(-0.9, 0.9); hTruthEtaLead_iso_R2->Scale(1/hTruthEtaLead_iso_R2->Integral()); hTruthEtaLead_iso_R2->Draw("same, P");

xj_functions::SetHist(hRecoEtaLead_all_R2,"#eta", "Counts",2,20,3,1.5);
hRecoEtaLead_all_R2->SetAxisRange(-0.9, 0.9); hRecoEtaLead_all_R2->Scale(1/hRecoEtaLead_all_R2->Integral()); hRecoEtaLead_all_R2->Draw("same, P");

xj_functions::SetHist(hRecoEtaLead_iso_R2,"#eta", "Counts",2,24,3,1.5);
hRecoEtaLead_iso_R2->SetAxisRange(-0.9, 0.9); hRecoEtaLead_iso_R2->Scale(1/hRecoEtaLead_iso_R2->Integral()); hRecoEtaLead_iso_R2->Draw("same, P");

TLegend *leg_eta_lead_R2 = new TLegend (0.25, 0.80, 0.75,0.90);
leg_eta_lead_R2->SetTextSize(0.026);
leg_eta_lead_R2->AddEntry(hTruthEtaLead_all_R2, "Truth R = 0.2 All", "P");
leg_eta_lead_R2->AddEntry(hTruthEtaLead_iso_R2, "Truth R = 0.2 Isolated", "P");
leg_eta_lead_R2->AddEntry(hRecoEtaLead_all_R2, "Reco R = 0.2 All", "P");
leg_eta_lead_R2->AddEntry(hRecoEtaLead_iso_R2, "Reco R = 0.2 Isolated", "P");
leg_eta_lead_R2->SetBorderSize(0);
leg_eta_lead_R2->Draw();

TLegend *sleg_eta_lead_R2;
xj_functions::DrawsPHENIXLegendNew(sleg_eta_lead_R2, "sleg_eta_lead_R2",0.25,0.25,0.75,0.55, "p+p #sqrt{s}=200 GeV", "anti-#it{k}_{#it{t}} #it{R} = 0.2", "20.9 < p_{T,1} < 60.77 GeV", "Matched Truth and Reco Leading Dijets",true, true, false, false, true);

can_eta_R2->cd(2); gPad->SetLogy();
// Truth + Reco Sublead, R = 0.2
xj_functions::SetHist(hTruthEtaSubLead_all_R2,"Subleading Truth #eta", "Integrated Counts",1,20,3,1.5);
hTruthEtaSubLead_all_R2->SetAxisRange(-0.9, 0.9); hTruthEtaSubLead_all_R2->Scale(1/hTruthEtaSubLead_all_R2->Integral()); hTruthEtaSubLead_all_R2->Draw("P");

xj_functions::SetHist(hTruthEtaSubLead_iso_R2,"#eta", "Counts",1,24,3,1.5);
hTruthEtaSubLead_iso_R2->SetAxisRange(-0.9, 0.9); hTruthEtaSubLead_iso_R2->Scale(1/hTruthEtaSubLead_iso_R2->Integral()); hTruthEtaSubLead_iso_R2->Draw("same, P");

xj_functions::SetHist(hRecoEtaSubLead_all_R2,"#eta", "Counts",2,20,3,1.5);
hRecoEtaSubLead_all_R2->SetAxisRange(-0.9, 0.9); hRecoEtaSubLead_all_R2->Scale(1/hRecoEtaSubLead_all_R2->Integral()); hRecoEtaSubLead_all_R2->Draw("same, P");

xj_functions::SetHist(hRecoEtaSubLead_iso_R2,"#eta", "Counts",2,24,3,1.5);
hRecoEtaSubLead_iso_R2->SetAxisRange(-0.9, 0.9); hRecoEtaSubLead_iso_R2->Scale(1/hRecoEtaSubLead_iso_R2->Integral()); hRecoEtaSubLead_iso_R2->Draw("same, P");

TLegend *leg_eta_sublead_R2 = new TLegend (0.25, 0.80, 0.75,0.90);
leg_eta_sublead_R2->SetTextSize(0.026);
leg_eta_sublead_R2->AddEntry(hTruthEtaSubLead_all_R2, "Truth R = 0.2 All", "P");
leg_eta_sublead_R2->AddEntry(hTruthEtaSubLead_iso_R2, "Truth R = 0.2 Isolated", "P");
leg_eta_sublead_R2->AddEntry(hRecoEtaSubLead_all_R2, "Reco R = 0.2 All", "P");
leg_eta_sublead_R2->AddEntry(hRecoEtaSubLead_iso_R2, "Reco R = 0.2 Isolated", "P");
leg_eta_sublead_R2->SetBorderSize(0);
leg_eta_sublead_R2->Draw();

TLegend *sleg_eta_sublead_R2;
xj_functions::DrawsPHENIXLegendNew(sleg_eta_sublead_R2, "sleg_eta_sublead_R2",0.25,0.25,0.75,0.55, "p+p #sqrt{s}=200 GeV", "anti-#it{k}_{#it{t}} #it{R} = 0.2", "20.9 < p_{T,1} < 60.77 GeV", "Matched Truth and Reco Subleading Dijets",true, true, false, false, true);

 // === R = 0.4 Phi Plots ===
TCanvas* can_phi_R4 = new TCanvas("can_phi_R4","",1600,800);
can_phi_R4->Divide(2,1);
can_phi_R4->cd(1); gPad->SetLogy();
// Truth + Reco Lead, R = 0.4
xj_functions::SetHist(hTruthPhiLead_all_R4,"Leading Truth #phi", "Integrated Counts",1,20,3,1.5);
hTruthPhiLead_all_R4->Scale(1/hTruthPhiLead_all_R4->Integral()); hTruthPhiLead_all_R4->Draw("P");

xj_functions::SetHist(hTruthPhiLead_iso_R4,"#phi", "Counts",1,24,3,1.5);
hTruthPhiLead_iso_R4->Scale(1/hTruthPhiLead_iso_R4->Integral()); hTruthPhiLead_iso_R4->Draw("same, P");

xj_functions::SetHist(hRecoPhiLead_all_R4,"#phi", "Counts",2,20,3,1.5);
hRecoPhiLead_all_R4->Scale(1/hRecoPhiLead_all_R4->Integral()); hRecoPhiLead_all_R4->Draw("same, P");

xj_functions::SetHist(hRecoPhiLead_iso_R4,"#phi", "Counts",2,24,3,1.5);
hRecoPhiLead_iso_R4->Scale(1/hRecoPhiLead_iso_R4->Integral()); hRecoPhiLead_iso_R4->Draw("same, P");

TLegend *leg_phi_lead_R4 = new TLegend (0.25, 0.80, 0.75,0.90);
leg_phi_lead_R4->SetTextSize(0.026);
leg_phi_lead_R4->AddEntry(hTruthPhiLead_all_R4, "Truth R = 0.4 All", "P");
leg_phi_lead_R4->AddEntry(hTruthPhiLead_iso_R4, "Truth R = 0.4 Isolated", "P");
leg_phi_lead_R4->AddEntry(hRecoPhiLead_all_R4, "Reco R = 0.4 All", "P");
leg_phi_lead_R4->AddEntry(hRecoPhiLead_iso_R4, "Reco R = 0.4 Isolated", "P");
leg_phi_lead_R4->SetBorderSize(0);
leg_phi_lead_R4->Draw();
TLegend *sleg_phi_lead_R4;
xj_functions::DrawsPHENIXLegendNew(sleg_phi_lead_R4, "sleg_phi_lead_R4",0.25,0.25,0.75,0.55, "p+p #sqrt{s}=200 GeV", "anti-#it{k}_{#it{t}} #it{R} = 0.4", "20.9 < p_{T,1} < 60.77 GeV", "Matched Truth and Reco Leading Dijets",true, true, false, false, true);

can_phi_R4->cd(2); gPad->SetLogy();
// Truth + Reco Sublead, R = 0.4
xj_functions::SetHist(hTruthPhiSubLead_all_R4,"Subleading Truth #phi", "Integrated Counts",1,20,3,1.5);
hTruthPhiSubLead_all_R4->Scale(1/hTruthPhiSubLead_all_R4->Integral()); hTruthPhiSubLead_all_R4->Draw("P");

xj_functions::SetHist(hTruthPhiSubLead_iso_R4,"#phi", "Counts",1,24,3,1.5);
hTruthPhiSubLead_iso_R4->Scale(1/hTruthPhiSubLead_iso_R4->Integral()); hTruthPhiSubLead_iso_R4->Draw("same, P");

xj_functions::SetHist(hRecoPhiSubLead_all_R4,"#phi", "Counts",2,20,3,1.5);
hRecoPhiSubLead_all_R4->Scale(1/hRecoPhiSubLead_all_R4->Integral()); hRecoPhiSubLead_all_R4->Draw("same, P");

xj_functions::SetHist(hRecoPhiSubLead_iso_R4,"#phi", "Counts",2,24,3,1.5);
hRecoPhiSubLead_iso_R4->Scale(1/hRecoPhiSubLead_iso_R4->Integral()); hRecoPhiSubLead_iso_R4->Draw("same, P");

TLegend *leg_phi_sublead_R4 = new TLegend (0.25, 0.80, 0.75,0.90);
leg_phi_sublead_R4->SetTextSize(0.026);
leg_phi_sublead_R4->AddEntry(hTruthPhiSubLead_all_R4, "Truth R = 0.4 All", "P");
leg_phi_sublead_R4->AddEntry(hTruthPhiSubLead_iso_R4, "Truth R = 0.4 Isolated", "P");
leg_phi_sublead_R4->AddEntry(hRecoPhiSubLead_all_R4, "Reco R = 0.4 All", "P");
leg_phi_sublead_R4->AddEntry(hRecoPhiSubLead_iso_R4, "Reco R = 0.4 Isolated", "P");
leg_phi_sublead_R4->SetBorderSize(0);
leg_phi_sublead_R4->Draw();

TLegend *sleg_phi_sublead_R4;
xj_functions::DrawsPHENIXLegendNew(sleg_phi_sublead_R4, "sleg_phi_sublead_R4",0.25,0.25,0.75,0.55, "p+p #sqrt{s}=200 GeV", "anti-#it{k}_{#it{t}} #it{R} = 0.4", "20.9 < p_{T,1} < 60.77 GeV", "Matched Truth and Reco Subleading Dijets",true, true, false, false, true);


// === R = 0.2 Phi Plots ===
TCanvas* can_phi_R2 = new TCanvas("can_phi_R2","",1600,800);
can_phi_R2->Divide(2,1);
can_phi_R2->cd(1); gPad->SetLogy();
// Truth + Reco Lead, R = 0.2
xj_functions::SetHist(hTruthPhiLead_all_R2,"Leading Truth #phi", "Integrated Counts",1,20,3,1.5);
hTruthPhiLead_all_R2->Scale(1/hTruthPhiLead_all_R2->Integral()); hTruthPhiLead_all_R2->Draw("P");

xj_functions::SetHist(hTruthPhiLead_iso_R2,"#phi", "Counts",1,24,3,1.5);
hTruthPhiLead_iso_R2->Scale(1/hTruthPhiLead_iso_R2->Integral()); hTruthPhiLead_iso_R2->Draw("same, P");

xj_functions::SetHist(hRecoPhiLead_all_R2,"#phi", "Counts",2,20,3,1.5);
hRecoPhiLead_all_R2->Scale(1/hRecoPhiLead_all_R2->Integral()); hRecoPhiLead_all_R2->Draw("same, P");

xj_functions::SetHist(hRecoPhiLead_iso_R2,"#phi", "Counts",2,24,3,1.5);
hRecoPhiLead_iso_R2->Scale(1/hRecoPhiLead_iso_R2->Integral()); hRecoPhiLead_iso_R2->Draw("same, P");

TLegend *leg_phi_lead_R2 = new TLegend (0.25, 0.80, 0.75,0.90);
leg_phi_lead_R2->SetTextSize(0.026);
leg_phi_lead_R2->AddEntry(hTruthPhiLead_all_R2, "Truth R = 0.2 All", "P");
leg_phi_lead_R2->AddEntry(hTruthPhiLead_iso_R2, "Truth R = 0.2 Isolated", "P");
leg_phi_lead_R2->AddEntry(hRecoPhiLead_all_R2, "Reco R = 0.2 All", "P");
leg_phi_lead_R2->AddEntry(hRecoPhiLead_iso_R2, "Reco R = 0.2 Isolated", "P");
leg_phi_lead_R2->SetBorderSize(0);
leg_phi_lead_R2->Draw();

TLegend *sleg_phi_lead_R2;
xj_functions::DrawsPHENIXLegendNew(sleg_phi_lead_R2, "sleg_phi_lead_R2",0.25,0.25,0.75,0.55, "p+p #sqrt{s}=200 GeV", "anti-#it{k}_{#it{t}} #it{R} = 0.2", "20.9 < p_{T,1} < 60.77 GeV", "Matched Truth and Reco Leading Dijets",true, true, false, false, true);

can_phi_R2->cd(2); gPad->SetLogy();
// Truth + Reco Sublead, R = 0.2
xj_functions::SetHist(hTruthPhiSubLead_all_R2,"Subleading Truth #phi", "Integrated Counts",1,20,3,1.5);
hTruthPhiSubLead_all_R2->Scale(1/hTruthPhiSubLead_all_R2->Integral()); hTruthPhiSubLead_all_R2->Draw("P");

xj_functions::SetHist(hTruthPhiSubLead_iso_R2,"#phi", "Counts",1,24,3,1.5);
hTruthPhiSubLead_iso_R2->Scale(1/hTruthPhiSubLead_iso_R2->Integral()); hTruthPhiSubLead_iso_R2->Draw("same, P");

xj_functions::SetHist(hRecoPhiSubLead_all_R2,"#phi", "Counts",2,20,3,1.5);
hRecoPhiSubLead_all_R2->Scale(1/hRecoPhiSubLead_all_R2->Integral()); hRecoPhiSubLead_all_R2->Draw("same, P");

xj_functions::SetHist(hRecoPhiSubLead_iso_R2,"#phi", "Counts",2,24,3,1.5);
hRecoPhiSubLead_iso_R2->Scale(1/hRecoPhiSubLead_iso_R2->Integral()); hRecoPhiSubLead_iso_R2->Draw("same, P");

TLegend *leg_phi_sublead_R2 = new TLegend (0.25, 0.80, 0.75,0.90);
leg_phi_sublead_R2->SetTextSize(0.026);
leg_phi_sublead_R2->AddEntry(hTruthPhiSubLead_all_R2, "Truth R = 0.2 All", "P");
leg_phi_sublead_R2->AddEntry(hTruthPhiSubLead_iso_R2, "Truth R = 0.2 Isolated", "P");
leg_phi_sublead_R2->AddEntry(hRecoPhiSubLead_all_R2, "Reco R = 0.2 All", "P");
leg_phi_sublead_R2->AddEntry(hRecoPhiSubLead_iso_R2, "Reco R = 0.2 Isolated", "P");
leg_phi_sublead_R2->SetBorderSize(0);
leg_phi_sublead_R2->Draw();

TLegend *sleg_phi_sublead_R2;
xj_functions::DrawsPHENIXLegendNew(sleg_phi_sublead_R2, "sleg_phi_sublead_R2",0.25,0.25,0.75,0.55, "p+p #sqrt{s}=200 GeV", "anti-#it{k}_{#it{t}} #it{R} = 0.2", "20.9 < p_{T,1} < 60.77 GeV", "Matched Truth and Reco Subleading Dijets",true, true, false, false, true);
 
}//close macro
