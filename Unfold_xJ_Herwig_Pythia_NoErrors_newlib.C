
#include "xj_functions.h"
#include "read_binning.h"

#include "sPhenixStyle.h"
#include "sPhenixStyle.C"


#include "RooUnfoldResponse.h"
#include "RooUnfold.h"
#include "RooUnfoldBayes.h"

#define _USE_MATH_DEFINES

#include <math.h> 
#include <cmath>
#include <iostream>

#include<vector>
#include<array>


void Unfold_xJ_Herwig_Pythia_NoErrors_newlib(
    string infilepythia1 = "Xj_2D_Response_Base_Bins19_TREE_JET_v5_10_new_ProdA_2024-00000021.root",
    string infilepythia2 = "Xj_2D_Response_Base_Bins19_TREE_JET_v5_30_new_ProdA_2024-00000021.root",
    string infileherwig1 = "Xj_2D_Response_Base_Bins19_All-Herwig-Jet10-R04-Run21.root",
    string infileherwig2 = "Xj_2D_Response_Base_Bins19_All-Herwig-Jet30-R04-Run21.root",
    const std::string configfile = "binning_original.config",
    int niterations = 6)
{
    // === Pythia Histograms ===
    TFile* infilepythia01 = TFile::Open(infilepythia1.c_str());
    RooUnfoldResponse *response_pt1pt2_full_pythia = (RooUnfoldResponse*)infilepythia01->Get("response_pt1pt2_full");
    RooUnfoldResponse *response_pt1pt2_half_fill_pythia = (RooUnfoldResponse*)infilepythia01->Get("response_pt1pt2_half_fill");
    TH2D *h_pt1pt2_truth_pythia = (TH2D*)infilepythia01->Get("h_pt1pt2_truth");
    TH2D *h_pt1pt2_reco_pythia = (TH2D*)infilepythia01->Get("h_pt1pt2_reco");
    TH2D *h_pt1pt2_truth_half_pythia = (TH2D*)infilepythia01->Get("h_pt1pt2_truth_half");
    TH2D *h_pt1pt2_reco_half_test_pythia = (TH2D*)infilepythia01->Get("h_pt1pt2_reco_half_test");
    TH1D* h_pt_lead_truth_matched_pythia = (TH1D*)infilepythia01->Get("h_pt_lead_truth_matched");
    TH1D* h_pt_lead_reco_matched_pythia = (TH1D*)infilepythia01->Get("h_pt_lead_reco_matched");

    TFile* infilepythia02 = TFile::Open(infilepythia2.c_str());
    RooUnfoldResponse *response_pt1pt2_full_30_pythia = (RooUnfoldResponse*)infilepythia02->Get("response_pt1pt2_full");
    RooUnfoldResponse *response_pt1pt2_half_fill_30_pythia = (RooUnfoldResponse*)infilepythia02->Get("response_pt1pt2_half_fill");
    TH2D *h_pt1pt2_truth_30_pythia = (TH2D*)infilepythia02->Get("h_pt1pt2_truth");
    TH2D *h_pt1pt2_reco_30_pythia = (TH2D*)infilepythia02->Get("h_pt1pt2_reco");
    TH2D *h_pt1pt2_truth_half_30_pythia = (TH2D*)infilepythia02->Get("h_pt1pt2_truth_half");
    TH2D *h_pt1pt2_reco_half_test_30_pythia = (TH2D*)infilepythia02->Get("h_pt1pt2_reco_half_test");
    TH1D* h_pt_lead_truth_matched_30_pythia = (TH1D*)infilepythia02->Get("h_pt_lead_truth_matched");
    TH1D* h_pt_lead_reco_matched_30_pythia = (TH1D*)infilepythia02->Get("h_pt_lead_reco_matched");

    response_pt1pt2_full_pythia->Add(*response_pt1pt2_full_30_pythia);
    response_pt1pt2_half_fill_pythia->Add(*response_pt1pt2_half_fill_30_pythia);
    h_pt1pt2_truth_pythia->Add(h_pt1pt2_truth_30_pythia);
    h_pt1pt2_reco_pythia->Add(h_pt1pt2_reco_30_pythia);
    h_pt1pt2_truth_half_pythia->Add(h_pt1pt2_truth_half_30_pythia);
    h_pt1pt2_reco_half_test_pythia->Add(h_pt1pt2_reco_half_test_30_pythia);
    h_pt_lead_truth_matched_pythia->Add(h_pt_lead_truth_matched_30_pythia);
    h_pt_lead_reco_matched_pythia->Add(h_pt_lead_reco_matched_30_pythia);

    // === Herwig Histograms ===
    TFile* infileherwig01 = TFile::Open(infileherwig1.c_str());
    RooUnfoldResponse *response_pt1pt2_full_herwig = (RooUnfoldResponse*)infileherwig01->Get("response_pt1pt2_full");
    RooUnfoldResponse *response_pt1pt2_half_fill_herwig = (RooUnfoldResponse*)infileherwig01->Get("response_pt1pt2_half_fill");
    TH2D *h_pt1pt2_truth_herwig = (TH2D*)infileherwig01->Get("h_pt1pt2_truth");
    TH2D *h_pt1pt2_reco_herwig = (TH2D*)infileherwig01->Get("h_pt1pt2_reco");
    TH2D *h_pt1pt2_truth_half_herwig = (TH2D*)infileherwig01->Get("h_pt1pt2_truth_half");
    TH2D *h_pt1pt2_reco_half_test_herwig = (TH2D*)infileherwig01->Get("h_pt1pt2_reco_half_test");
    TH1D* h_pt_lead_truth_matched_herwig = (TH1D*)infileherwig01->Get("h_pt_lead_truth_matched");
    TH1D* h_pt_lead_reco_matched_herwig = (TH1D*)infileherwig01->Get("h_pt_lead_reco_matched");

    TFile* infileherwig02 = TFile::Open(infileherwig2.c_str());
    RooUnfoldResponse *response_pt1pt2_full_30_herwig = (RooUnfoldResponse*)infileherwig02->Get("response_pt1pt2_full");
    RooUnfoldResponse *response_pt1pt2_half_fill_30_herwig = (RooUnfoldResponse*)infileherwig02->Get("response_pt1pt2_half_fill");
    TH2D *h_pt1pt2_truth_30_herwig = (TH2D*)infileherwig02->Get("h_pt1pt2_truth");
    TH2D *h_pt1pt2_reco_30_herwig = (TH2D*)infileherwig02->Get("h_pt1pt2_reco");
    TH2D *h_pt1pt2_truth_half_30_herwig = (TH2D*)infileherwig02->Get("h_pt1pt2_truth_half");
    TH2D *h_pt1pt2_reco_half_test_30_herwig = (TH2D*)infileherwig02->Get("h_pt1pt2_reco_half_test");
    TH1D* h_pt_lead_truth_matched_30_herwig = (TH1D*)infileherwig02->Get("h_pt_lead_truth_matched");
    TH1D* h_pt_lead_reco_matched_30_herwig = (TH1D*)infileherwig02->Get("h_pt_lead_reco_matched");

    response_pt1pt2_full_herwig->Add(*response_pt1pt2_full_30_herwig);
    response_pt1pt2_half_fill_herwig->Add(*response_pt1pt2_half_fill_30_herwig);
    h_pt1pt2_truth_herwig->Add(h_pt1pt2_truth_30_herwig);
    h_pt1pt2_reco_herwig->Add(h_pt1pt2_reco_30_herwig);
    h_pt1pt2_truth_half_herwig->Add(h_pt1pt2_truth_half_30_herwig);
    h_pt1pt2_reco_half_test_herwig->Add(h_pt1pt2_reco_half_test_30_herwig);
    h_pt_lead_truth_matched_herwig->Add(h_pt_lead_truth_matched_30_herwig);
    h_pt_lead_reco_matched_herwig->Add(h_pt_lead_reco_matched_30_herwig);

    // Now both _pythia and _herwig datasets are fully loaded and merged for 10+30 jet cuts.

// === Plot pT lead comparison: Pythia + Herwig ===
TCanvas* can_ptqa = new TCanvas("can_ptqa", "p_{T1} vs p_{T2}", 800, 800);
can_ptqa->cd(1);
gPad->SetRightMargin(0.16);
gPad->SetLogy();

xj_functions::SetHist(h_pt_lead_truth_matched_pythia, "p_{T1, truth} [GeV]", "Counts", 1, 20, 1, 1);
h_pt_lead_truth_matched_pythia->Draw("P");

xj_functions::SetHist(h_pt_lead_reco_matched_pythia, "p_{T1, truth} [GeV]", "Counts", 2, 20, 1, 1);
h_pt_lead_reco_matched_pythia->Draw("same");

xj_functions::SetHist(h_pt_lead_truth_matched_herwig, "p_{T1, truth} [GeV]", "Counts", 4, 20, 1, 1);
h_pt_lead_truth_matched_herwig->Draw("same");

xj_functions::SetHist(h_pt_lead_reco_matched_herwig, "p_{T1, truth} [GeV]", "Counts", 6, 20, 1, 1);
h_pt_lead_reco_matched_herwig->Draw("same");

// xj histograms
read_binning rb(configfile.c_str());
Int_t read_nbins = rb.get_nbins();
const int nbins = read_nbins;
float ipt_bins[nbins + 1];
float ixj_bins[nbins + 1];

rb.get_pt_bins(ipt_bins);
rb.get_xj_bins(ixj_bins);
for (int i = 0; i < nbins + 1; i++) {
  std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
}

TH1D *h_xj_reco_full_pythia   = new TH1D("h_xj_reco_full_pythia",   ";x_{J};", nbins, ixj_bins);
TH1D *h_xj_truth_full_pythia  = new TH1D("h_xj_truth_full_pythia",  ";x_{J};", nbins, ixj_bins);
TH1D *h_xj_unfold_full_pythia = new TH1D("h_xj_unfold_full_pythia", ";x_{J};", nbins, ixj_bins);
TH1D *h_xj_reco_half_pythia   = new TH1D("h_xj_reco_half_pythia",   ";x_{J};", nbins, ixj_bins);
TH1D *h_xj_truth_half_pythia  = new TH1D("h_xj_truth_half_pythia",  ";x_{J};", nbins, ixj_bins);
TH1D *h_xj_unfold_half_pythia = new TH1D("h_xj_unfold_half_pythia", ";x_{J};", nbins, ixj_bins);

TH1D* h_xj_classical_truth_pythia     = (TH1D*)infilepythia01->Get("h_xj_classical_truth");
TH1D* h_xj_classical_truth_30_pythia  = (TH1D*)infilepythia02->Get("h_xj_classical_truth");
h_xj_classical_truth_pythia->Add(h_xj_classical_truth_30_pythia);

TH1D* h_xj_classical_reco_pythia      = (TH1D*)infilepythia01->Get("h_xj_classical_reco");
TH1D* h_xj_classical_reco_30_pythia   = (TH1D*)infilepythia02->Get("h_xj_classical_reco");
h_xj_classical_reco_pythia->Add(h_xj_classical_reco_30_pythia);


TH1D *h_xj_reco_full_herwig   = new TH1D("h_xj_reco_full_herwig",   ";x_{J};", nbins, ixj_bins);
TH1D *h_xj_truth_full_herwig  = new TH1D("h_xj_truth_full_herwig",  ";x_{J};", nbins, ixj_bins);
TH1D *h_xj_unfold_full_herwig = new TH1D("h_xj_unfold_full_herwig", ";x_{J};", nbins, ixj_bins);
TH1D *h_xj_reco_half_herwig   = new TH1D("h_xj_reco_half_herwig",   ";x_{J};", nbins, ixj_bins);
TH1D *h_xj_truth_half_herwig  = new TH1D("h_xj_truth_half_herwig",  ";x_{J};", nbins, ixj_bins);
TH1D *h_xj_unfold_half_herwig = new TH1D("h_xj_unfold_half_herwig", ";x_{J};", nbins, ixj_bins);

TH1D* h_xj_classical_truth_herwig     = (TH1D*)infileherwig01->Get("h_xj_classical_truth");
TH1D* h_xj_classical_truth_30_herwig  = (TH1D*)infileherwig02->Get("h_xj_classical_truth");
h_xj_classical_truth_herwig->Add(h_xj_classical_truth_30_herwig);

TH1D* h_xj_classical_reco_herwig      = (TH1D*)infileherwig01->Get("h_xj_classical_reco");
TH1D* h_xj_classical_reco_30_herwig   = (TH1D*)infileherwig02->Get("h_xj_classical_reco");
h_xj_classical_reco_herwig->Add(h_xj_classical_reco_30_herwig);



  //Read in Data File and Set up Histos

  TFile* infiledata = TFile::Open("Hists_xJ_ProjectionTest_TREE_DIJET_v6_1_ana462_2024p010_v001_gl10-00047352-00047733.root");
  TH2D *h_pt1pt2_data = (TH2D*)infiledata->Get("h_pt1pt2");
  TH1D *h_xj_data = new TH1D("h_xj_data", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold_data_pythia = new TH1D("h_xj_unfold_data_pythia", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold_data_herwig = new TH1D("h_xj_unfold_data_herwig", ";x_{J};",nbins, ixj_bins);

  //other functions

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

  //clean data 
  xj_functions::clean_pt1pt2(h_pt1pt2_data,nbins,1);
 
// Unfolding Process

// Pythia
TH2D* h_unfolded_pt1pt2_full_pythia;
TH2D* h_unfolded_pt1pt2_half_pythia;
TH2D* h_unfolded_pt1pt2_data_pythia;

RooUnfoldBayes* unfolded_pt1pt2_full_pythia = new RooUnfoldBayes(response_pt1pt2_full_pythia, h_pt1pt2_reco_pythia, niterations); 
RooUnfoldBayes* unfolded_pt1pt2_half_pythia = new RooUnfoldBayes(response_pt1pt2_half_fill_pythia, h_pt1pt2_reco_half_test_pythia, niterations); 
h_unfolded_pt1pt2_full_pythia = (TH2D*)unfolded_pt1pt2_full_pythia->Hunfold(RooUnfold::ErrorTreatment::kCovariance);
h_unfolded_pt1pt2_half_pythia = (TH2D*)unfolded_pt1pt2_half_pythia->Hunfold(RooUnfold::ErrorTreatment::kCovariance);

RooUnfoldBayes* unfolded_pt1pt2_data_pythia = new RooUnfoldBayes(response_pt1pt2_full_pythia, h_pt1pt2_data, niterations);
h_unfolded_pt1pt2_data_pythia = (TH2D*)unfolded_pt1pt2_data_pythia->Hunfold(RooUnfold::ErrorTreatment::kCovariance);

// Herwig
TH2D* h_unfolded_pt1pt2_full_herwig;
TH2D* h_unfolded_pt1pt2_half_herwig;
TH2D* h_unfolded_pt1pt2_data_herwig;

RooUnfoldBayes* unfolded_pt1pt2_full_herwig = new RooUnfoldBayes(response_pt1pt2_full_herwig, h_pt1pt2_reco_herwig, niterations); 
RooUnfoldBayes* unfolded_pt1pt2_half_herwig = new RooUnfoldBayes(response_pt1pt2_half_fill_herwig, h_pt1pt2_reco_half_test_herwig, niterations); 
h_unfolded_pt1pt2_full_herwig = (TH2D*)unfolded_pt1pt2_full_herwig->Hunfold(RooUnfold::ErrorTreatment::kCovariance);
h_unfolded_pt1pt2_half_herwig = (TH2D*)unfolded_pt1pt2_half_herwig->Hunfold(RooUnfold::ErrorTreatment::kCovariance);

RooUnfoldBayes* unfolded_pt1pt2_data_herwig = new RooUnfoldBayes(response_pt1pt2_full_herwig, h_pt1pt2_data, niterations);
h_unfolded_pt1pt2_data_herwig = (TH2D*)unfolded_pt1pt2_data_herwig->Hunfold(RooUnfold::ErrorTreatment::kCovariance);

// Clean data unfoldings
xj_functions::clean_pt1pt2(h_unfolded_pt1pt2_data_pythia, nbins, 1);
xj_functions::clean_pt1pt2(h_unfolded_pt1pt2_data_herwig, nbins, 1);

// Project xJ (Pythia)
xj_functions::project_xj(h_pt1pt2_reco_pythia,         h_xj_reco_full_pythia,    nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
xj_functions::project_xj(h_pt1pt2_truth_pythia,        h_xj_truth_full_pythia,   nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
xj_functions::project_xj(h_unfolded_pt1pt2_full_pythia,h_xj_unfold_full_pythia,  nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
xj_functions::project_xj(h_pt1pt2_reco_half_test_pythia,h_xj_reco_half_pythia,   nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
xj_functions::project_xj(h_pt1pt2_truth_half_pythia,   h_xj_truth_half_pythia,   nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
xj_functions::project_xj(h_unfolded_pt1pt2_half_pythia,h_xj_unfold_half_pythia,  nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
xj_functions::project_xj(h_unfolded_pt1pt2_data_pythia,h_xj_unfold_data_pythia,  nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
xj_functions::project_xj(h_pt1pt2_data,                h_xj_data,                nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);

// Project xJ (Herwig)
xj_functions::project_xj(h_pt1pt2_reco_herwig,         h_xj_reco_full_herwig,    nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
xj_functions::project_xj(h_pt1pt2_truth_herwig,        h_xj_truth_full_herwig,   nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
xj_functions::project_xj(h_unfolded_pt1pt2_full_herwig,h_xj_unfold_full_herwig,  nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
xj_functions::project_xj(h_pt1pt2_reco_half_test_herwig,h_xj_reco_half_herwig,   nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
xj_functions::project_xj(h_pt1pt2_truth_half_herwig,   h_xj_truth_half_herwig,   nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
xj_functions::project_xj(h_unfolded_pt1pt2_half_herwig,h_xj_unfold_half_herwig,  nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
xj_functions::project_xj(h_unfolded_pt1pt2_data_herwig,h_xj_unfold_data_herwig,  nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);

// Normalize
xj_functions::normalize_histo(h_xj_data, nbins);
xj_functions::normalize_histo(h_xj_unfold_data_pythia, nbins);
xj_functions::normalize_histo(h_xj_unfold_data_herwig, nbins);

xj_functions::normalize_histo(h_xj_truth_full_pythia, nbins);
xj_functions::normalize_histo(h_xj_reco_full_pythia, nbins);
xj_functions::normalize_histo(h_xj_truth_half_pythia, nbins);
xj_functions::normalize_histo(h_xj_reco_half_pythia, nbins);
xj_functions::normalize_histo(h_xj_unfold_full_pythia, nbins);
xj_functions::normalize_histo(h_xj_unfold_half_pythia, nbins);
xj_functions::normalize_histo(h_xj_classical_truth_pythia, nbins);
xj_functions::normalize_histo(h_xj_classical_reco_pythia, nbins);

xj_functions::normalize_histo(h_xj_truth_full_herwig, nbins);
xj_functions::normalize_histo(h_xj_reco_full_herwig, nbins);
xj_functions::normalize_histo(h_xj_truth_half_herwig, nbins);
xj_functions::normalize_histo(h_xj_reco_half_herwig, nbins);
xj_functions::normalize_histo(h_xj_unfold_full_herwig, nbins);
xj_functions::normalize_histo(h_xj_unfold_half_herwig, nbins);

  //draw histograms
SetsPhenixStyle();
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
 // 2D histograms
TCanvas* can_pt1pt2 = new TCanvas("can_pt1pt2", "p_{T1} vs p_{T2}", 2400, 1600);
can_pt1pt2->Divide(3, 2);

// Truth
can_pt1pt2->cd(1);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_pt1pt2_truth_pythia, "p_{T1, truth} [GeV]", "p_{T2, truth} [GeV]", 2, 20, 1, 1);
h_pt1pt2_truth_pythia->Draw("Colz");
h_pt1pt2_truth_pythia->SetAxisRange(0, 89, "X");
h_pt1pt2_truth_pythia->SetAxisRange(0, 89, "Y");
h_pt1pt2_truth_pythia->SetAxisRange(pow(10, -9), 1, "Z");
h_pt1pt2_truth_pythia->GetZaxis()->SetLabelSize(0.035);

// Reco
can_pt1pt2->cd(2);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_pt1pt2_reco_pythia, "p_{T1, reco} [GeV]", "p_{T2, reco} [GeV]", 2, 20, 1, 1);
h_pt1pt2_reco_pythia->Draw("Colz");
h_pt1pt2_reco_pythia->SetAxisRange(0, 89, "X");
h_pt1pt2_reco_pythia->SetAxisRange(0, 89, "Y");
h_pt1pt2_reco_pythia->SetAxisRange(pow(10, -9), 1, "Z");
h_pt1pt2_reco_pythia->GetZaxis()->SetLabelSize(0.035);

// Unfolded Full
can_pt1pt2->cd(3);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_unfolded_pt1pt2_full_pythia, "p_{T1, unfold} [GeV]", "p_{T2, unfold} [GeV]", 2, 20, 1, 1);
h_unfolded_pt1pt2_full_pythia->Draw("Colz");
h_unfolded_pt1pt2_full_pythia->SetAxisRange(0, 89, "X");
h_unfolded_pt1pt2_full_pythia->SetAxisRange(0, 89, "Y");
h_unfolded_pt1pt2_full_pythia->SetAxisRange(pow(10, -9), 1, "Z");
h_unfolded_pt1pt2_full_pythia->GetZaxis()->SetLabelSize(0.035);

// Truth Half
can_pt1pt2->cd(4);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_pt1pt2_truth_half_pythia, "p_{T1, truth_half} [GeV]", "p_{T2, truth_half} [GeV]", 2, 20, 1, 1);
h_pt1pt2_truth_half_pythia->Draw("Colz");
h_pt1pt2_truth_half_pythia->SetAxisRange(0, 89, "X");
h_pt1pt2_truth_half_pythia->SetAxisRange(0, 89, "Y");
h_pt1pt2_truth_half_pythia->SetAxisRange(pow(10, -9), 1, "Z");
h_pt1pt2_truth_half_pythia->GetZaxis()->SetLabelSize(0.035);

// Reco Half
can_pt1pt2->cd(5);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_pt1pt2_reco_half_test_pythia, "p_{T1, reco_half_test} [GeV]", "p_{T2, reco_half_test} [GeV]", 2, 20, 1, 1);
h_pt1pt2_reco_half_test_pythia->Draw("Colz");
h_pt1pt2_reco_half_test_pythia->SetAxisRange(0, 89, "X");
h_pt1pt2_reco_half_test_pythia->SetAxisRange(0, 89, "Y");
h_pt1pt2_reco_half_test_pythia->SetAxisRange(pow(10, -9), 1, "Z");
h_pt1pt2_reco_half_test_pythia->GetZaxis()->SetLabelSize(0.035);

// Unfold Half
can_pt1pt2->cd(6);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_unfolded_pt1pt2_half_pythia, "p_{T1, unfold half} [GeV]", "p_{T2, unfold half} [GeV]", 2, 20, 1, 1);
h_unfolded_pt1pt2_half_pythia->Draw("Colz");
h_unfolded_pt1pt2_half_pythia->SetAxisRange(0, 89, "X");
h_unfolded_pt1pt2_half_pythia->SetAxisRange(0, 89, "Y");
h_unfolded_pt1pt2_half_pythia->SetAxisRange(pow(10, -9), 1, "Z");
h_unfolded_pt1pt2_half_pythia->GetZaxis()->SetLabelSize(0.035);

// Legends for can_pt1pt2
can_pt1pt2->cd(1); TLegend* full_truth;
xj_functions::DrawsPHENIXLegend(full_truth, "fulltruth", 0.10, 0.68, 0.5, 0.93,"p+p #sqrt{s}=200 GeV (Pythia)", true,false,false,false);
can_pt1pt2->cd(2); TLegend* full_reco;
xj_functions::DrawsPHENIXLegend(full_reco, "fullreco", 0.10, 0.68, 0.5, 0.93,"p+p #sqrt{s}=200 GeV  (Pythia)", false,true,false,false);
can_pt1pt2->cd(3); TLegend* full_unfold;
xj_functions::DrawsPHENIXLegend(full_unfold, "fullunfold", 0.10, 0.68, 0.5, 0.93,"p+p #sqrt{s}=200 GeV (Pythia)", true,true,false,false);
can_pt1pt2->cd(4); TLegend* half_truth;
xj_functions::DrawsPHENIXLegend(half_truth, "halftruth", 0.10, 0.68, 0.5, 0.93,"p+p #sqrt{s}=200 GeV  (Pythia)", true,false,false,false);
can_pt1pt2->cd(5); TLegend* half_reco;
xj_functions::DrawsPHENIXLegend(half_reco, "halfreco", 0.10, 0.68, 0.5, 0.93,"p+p #sqrt{s}=200 GeV (Pythia)", false,true,false,false);
can_pt1pt2->cd(6); TLegend* half_unfold;
xj_functions::DrawsPHENIXLegend(half_unfold, "halfunfold", 0.10, 0.68, 0.5, 0.93,"p+p #sqrt{s}=200 GeV (Pythia)", true,true,false,false);

// 2D histograms for HERWIG
TCanvas* can_pt1pt2_herwig = new TCanvas("can_pt1pt2_herwig", "p_{T1} vs p_{T2} (Herwig)", 2400, 1600);
can_pt1pt2_herwig->Divide(3, 2);

// Truth
can_pt1pt2_herwig->cd(1);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_pt1pt2_truth_herwig, "p_{T1, truth} [GeV]", "p_{T2, truth} [GeV]", 2, 20, 1, 1);
h_pt1pt2_truth_herwig->Draw("Colz");
h_pt1pt2_truth_herwig->SetAxisRange(0, 89, "X");
h_pt1pt2_truth_herwig->SetAxisRange(0, 89, "Y");
h_pt1pt2_truth_herwig->SetAxisRange(pow(10, -9), 1, "Z");
h_pt1pt2_truth_herwig->GetZaxis()->SetLabelSize(0.035);

// Reco
can_pt1pt2_herwig->cd(2);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_pt1pt2_reco_herwig, "p_{T1, reco} [GeV]", "p_{T2, reco} [GeV]", 2, 20, 1, 1);
h_pt1pt2_reco_herwig->Draw("Colz");
h_pt1pt2_reco_herwig->SetAxisRange(0, 89, "X");
h_pt1pt2_reco_herwig->SetAxisRange(0, 89, "Y");
h_pt1pt2_reco_herwig->SetAxisRange(pow(10, -9), 1, "Z");
h_pt1pt2_reco_herwig->GetZaxis()->SetLabelSize(0.035);

// Unfolded Full
can_pt1pt2_herwig->cd(3);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_unfolded_pt1pt2_full_herwig, "p_{T1, unfold} [GeV]", "p_{T2, unfold} [GeV]", 2, 20, 1, 1);
h_unfolded_pt1pt2_full_herwig->Draw("Colz");
h_unfolded_pt1pt2_full_herwig->SetAxisRange(0, 89, "X");
h_unfolded_pt1pt2_full_herwig->SetAxisRange(0, 89, "Y");
h_unfolded_pt1pt2_full_herwig->SetAxisRange(pow(10, -9), 1, "Z");
h_unfolded_pt1pt2_full_herwig->GetZaxis()->SetLabelSize(0.035);

// Truth Half
can_pt1pt2_herwig->cd(4);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_pt1pt2_truth_half_herwig, "p_{T1, truth_half} [GeV]", "p_{T2, truth_half} [GeV]", 2, 20, 1, 1);
h_pt1pt2_truth_half_herwig->Draw("Colz");
h_pt1pt2_truth_half_herwig->SetAxisRange(0, 89, "X");
h_pt1pt2_truth_half_herwig->SetAxisRange(0, 89, "Y");
h_pt1pt2_truth_half_herwig->SetAxisRange(pow(10, -9), 1, "Z");
h_pt1pt2_truth_half_herwig->GetZaxis()->SetLabelSize(0.035);

// Reco Half
can_pt1pt2_herwig->cd(5);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_pt1pt2_reco_half_test_herwig, "p_{T1, reco_half_test} [GeV]", "p_{T2, reco_half_test} [GeV]", 2, 20, 1, 1);
h_pt1pt2_reco_half_test_herwig->Draw("Colz");
h_pt1pt2_reco_half_test_herwig->SetAxisRange(0, 89, "X");
h_pt1pt2_reco_half_test_herwig->SetAxisRange(0, 89, "Y");
h_pt1pt2_reco_half_test_herwig->SetAxisRange(pow(10, -9), 1, "Z");
h_pt1pt2_reco_half_test_herwig->GetZaxis()->SetLabelSize(0.035);

// Unfold Half
can_pt1pt2_herwig->cd(6);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_unfolded_pt1pt2_half_herwig, "p_{T1, unfold half} [GeV]", "p_{T2, unfold half} [GeV]", 2, 20, 1, 1);
h_unfolded_pt1pt2_half_herwig->Draw("Colz");
h_unfolded_pt1pt2_half_herwig->SetAxisRange(0, 89, "X");
h_unfolded_pt1pt2_half_herwig->SetAxisRange(0, 89, "Y");
h_unfolded_pt1pt2_half_herwig->SetAxisRange(pow(10, -9), 1, "Z");
h_unfolded_pt1pt2_half_herwig->GetZaxis()->SetLabelSize(0.035);

// Legends for can_pt1pt2_herwig
can_pt1pt2_herwig->cd(1); TLegend* full_truth_herwig;
xj_functions::DrawsPHENIXLegend(full_truth_herwig, "fulltruth_herwig", 0.10, 0.68, 0.5, 0.93, "p+p #sqrt{s}=200 GeV (Herwig)", true,false,false,false);
can_pt1pt2_herwig->cd(2); TLegend* full_reco_herwig;
xj_functions::DrawsPHENIXLegend(full_reco_herwig, "fullreco_herwig", 0.10, 0.68, 0.5, 0.93, "p+p #sqrt{s}=200 GeV (Herwig)", false,true,false,false);
can_pt1pt2_herwig->cd(3); TLegend* full_unfold_herwig;
xj_functions::DrawsPHENIXLegend(full_unfold_herwig, "fullunfold_herwig", 0.10, 0.68, 0.5, 0.93, "p+p #sqrt{s}=200 GeV (Herwig)", true,true,false,false);
can_pt1pt2_herwig->cd(4); TLegend* half_truth_herwig;
xj_functions::DrawsPHENIXLegend(half_truth_herwig, "halftruth_herwig", 0.10, 0.68, 0.5, 0.93, "p+p #sqrt{s}=200 GeV (Herwig)", true,false,false,false);
can_pt1pt2_herwig->cd(5); TLegend* half_reco_herwig;
xj_functions::DrawsPHENIXLegend(half_reco_herwig, "halfreco_herwig", 0.10, 0.68, 0.5, 0.93, "p+p #sqrt{s}=200 GeV (Herwig)", false,true,false,false);
can_pt1pt2_herwig->cd(6); TLegend* half_unfold_herwig;
xj_functions::DrawsPHENIXLegend(half_unfold_herwig, "halfunfold_herwig", 0.10, 0.68, 0.5, 0.93, "p+p #sqrt{s}=200 GeV (Herwig)", true,true,false,false);


// XJ: Classical vs Truth
TCanvas* can_xj_projected = new TCanvas("can_xj_projected", "p_{T1} vs p_{T2}", 700, 800);
xj_functions::plotxj_and_ratio(can_xj_projected, h_xj_classical_reco_pythia, h_xj_classical_truth_pythia, h_xj_truth_full_pythia, h_xj_reco_full_pythia, true, "x_{J}", "Classical/Truth", 0, 1.0, 4, 53, 2, 20, kGreen+1, 20, 8, 53);

// XJ: Full
TCanvas* can_xj_full = new TCanvas("can_xj_full", "p_{T1} vs p_{T2}", 700, 800);
xj_functions::plotxj_and_ratio(can_xj_full, h_xj_reco_full_pythia, h_xj_unfold_full_pythia, h_xj_truth_full_pythia, h_xj_reco_full_pythia, false, "x_{J}", "Unfold/Truth", 0, 1.0, 4, 20, 4, 53, kGreen+1, 20, 8, 53);
TLegend* xj_full_unfold;
xj_functions::draw_xj_legend(xj_full_unfold, "xj_full_unfold", 0.17, 0.45, 0.5, 0.59, h_xj_truth_full_pythia, h_xj_reco_full_pythia, h_xj_unfold_full_pythia, h_xj_unfold_full_pythia, "PYTHIA 8", "PTYHIA 8 Reco", "Reco Unfolded", " ");
TLegend* xj_full_2;
xj_functions::DrawsPHENIXLegend(xj_full_2, "xj_full_2", 0.10, 0.60, 0.5, 0.93, "p+p #sqrt{s}=200 GeV", true, true, false, true);
xj_full_2->AddEntry("", "Full Closure Test", "");

// XJ: Half
TCanvas* can_xj_half = new TCanvas("can_xj_half", "p_{T1} vs p_{T2}", 700, 800);
xj_functions::plotxj_and_ratio(can_xj_half, h_xj_reco_half_pythia, h_xj_unfold_half_pythia, h_xj_truth_half_pythia, h_xj_reco_full_pythia, false, "x_{J}", "Unfold/Truth", 0, 1.0, 4, 20, 4, 53, kGreen+1, 20, 8, 53);
TLegend* xj_half_unfold;
xj_functions::draw_xj_legend(xj_half_unfold, "xj_half_unfold", 0.17, 0.45, 0.50, 0.59, h_xj_truth_half_pythia, h_xj_reco_half_pythia, h_xj_unfold_half_pythia, h_xj_unfold_half_pythia, "PYTHIA 8", "PTYHIA 8 Reco", "Reco Unfolded", " ");
TLegend* xj_half_2;
xj_functions::DrawsPHENIXLegend(xj_half_2, "xj_half_2", 0.10, 0.60, 0.5, 0.93, "p+p #sqrt{s}=200 GeV", true, true, false, true);
xj_half_2->AddEntry("", "Half Closure Test", "");

// XJ: Full (HERWIG)
TCanvas* can_xj_full_herwig = new TCanvas("can_xj_full_herwig", "p_{T1} vs p_{T2} (Herwig)", 700, 800);
xj_functions::plotxj_and_ratio(can_xj_full_herwig, h_xj_reco_full_herwig, h_xj_unfold_full_herwig, h_xj_truth_full_herwig, h_xj_reco_full_herwig, false, "x_{J}", "Unfold/Truth", 0, 1.0, 4, 20, 4, 53, kGreen+1, 20, 8, 53);
TLegend* xj_full_unfold_herwig;
xj_functions::draw_xj_legend(xj_full_unfold_herwig, "xj_full_unfold_herwig", 0.17, 0.45, 0.5, 0.59, h_xj_truth_full_herwig, h_xj_reco_full_herwig, h_xj_unfold_full_herwig, h_xj_unfold_full_herwig, "HERWIG", "HERWIG Reco", "Reco Unfolded", " ");
TLegend* xj_full_2_herwig;
xj_functions::DrawsPHENIXLegend(xj_full_2_herwig, "xj_full_2_herwig", 0.10, 0.60, 0.5, 0.93, "p+p #sqrt{s}=200 GeV", true, true, false, true);
xj_full_2_herwig->AddEntry("", "Full Closure Test", "");

// XJ: Half (HERWIG)
TCanvas* can_xj_half_herwig = new TCanvas("can_xj_half_herwig", "p_{T1} vs p_{T2} (Herwig)", 700, 800);
xj_functions::plotxj_and_ratio(can_xj_half_herwig, h_xj_reco_half_herwig, h_xj_unfold_half_herwig, h_xj_truth_half_herwig, h_xj_reco_full_herwig, false, "x_{J}", "Unfold/Truth", 0, 1.0, 4, 20, 4, 53, kGreen+1, 20, 8, 53);
TLegend* xj_half_unfold_herwig;
xj_functions::draw_xj_legend(xj_half_unfold_herwig, "xj_half_unfold_herwig", 0.17, 0.45, 0.50, 0.59, h_xj_truth_half_herwig, h_xj_reco_half_herwig, h_xj_unfold_half_herwig, h_xj_unfold_half_herwig, "HERWIG", "HERWIG Reco", "Reco Unfolded", " ");
TLegend* xj_half_2_herwig;
xj_functions::DrawsPHENIXLegend(xj_half_2_herwig, "xj_half_2_herwig", 0.10, 0.60, 0.5, 0.93, "p+p #sqrt{s}=200 GeV", true, true, false, true);
xj_half_2_herwig->AddEntry("", "Half Closure Test", "");

// Data Histogram
TCanvas* can_xj_data = new TCanvas("can_xj_data", "p_{T1} vs p_{T2}", 700, 800);

// Reduce bottom margin to minimize space between plots
xj_functions::plotxj_and_ratio(can_xj_data, h_xj_data, h_xj_unfold_data_pythia, h_xj_truth_full_pythia, h_xj_data,
                               false, "x_{J, Data}", "Unfold/Truth", 0, 1.0, 1, 20, 1, 53, kGreen+1, 20, 1, 20);

can_xj_data->cd(1);
TLegend* xj_data_unfold;
xj_functions::draw_xj_legend(xj_data_unfold, "xj_data_unfold", 0.17, 0.45, 0.50, 0.59,
                             h_xj_truth_full_pythia, h_xj_data, h_xj_unfold_data_pythia, h_xj_unfold_half_pythia,
                             "PYTHIA 8", "Data", "Data Unfolded", " ");
TLegend* xj_data_2;
xj_functions::DrawsPHENIXLegend(xj_data_2, "xj_data_2", 0.10, 0.60, 0.5, 0.93,
                                "p+p #sqrt{s}=200 GeV  Runs 47352-47733", true, false, true, true);

// 2D histogram Data
TCanvas* can_pt1pt2_data = new TCanvas("can_pt1pt2_data", "p_{T1} vs p_{T2}", 2400, 600);
can_pt1pt2_data->Divide(3, 1);

// Truth
can_pt1pt2_data->cd(1);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_pt1pt2_truth_pythia, "p_{T1, truth} [GeV]", "p_{T2, truth} [GeV]", 2, 20, 1, 1);
h_pt1pt2_truth_pythia->Draw("Colz");
h_pt1pt2_truth_pythia->SetAxisRange(0, 89, "X");
h_pt1pt2_truth_pythia->SetAxisRange(0, 89, "Y");
h_pt1pt2_truth_pythia->SetAxisRange(pow(10, -9), 1, "Z");
h_pt1pt2_truth_pythia->GetZaxis()->SetLabelSize(0.035);

// Data
can_pt1pt2_data->cd(2);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_pt1pt2_data, "p_{T1, data} [GeV]", "p_{T2, data} [GeV]", 2, 20, 1, 1);
h_pt1pt2_data->Draw("Colz");
h_pt1pt2_data->SetAxisRange(0, 89, "X");
h_pt1pt2_data->SetAxisRange(0, 89, "Y");
h_pt1pt2_data->SetAxisRange(1, pow(10, 7), "Z");
h_pt1pt2_data->GetZaxis()->SetLabelSize(0.035);

// Unfolded data
can_pt1pt2_data->cd(3);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_unfolded_pt1pt2_data_pythia, "p_{T1, unfold} [GeV]", "p_{T2, unfold} [GeV]", 2, 20, 1, 1);
h_unfolded_pt1pt2_data_pythia->Draw("Colz");
h_unfolded_pt1pt2_data_pythia->SetAxisRange(0, 89, "X");
h_unfolded_pt1pt2_data_pythia->SetAxisRange(0, 89, "Y");
h_unfolded_pt1pt2_data_pythia->SetAxisRange(1, pow(10, 7), "Z");
h_unfolded_pt1pt2_data_pythia->GetZaxis()->SetLabelSize(0.035);

// Legends
can_pt1pt2_data->cd(1); TLegend* data_truth;
xj_functions::DrawsPHENIXLegend(data_truth, "datatruth", 0.10, 0.68, 0.5, 0.93,
                                "p+p #sqrt{s}=200 GeV", true, false, false);
can_pt1pt2_data->cd(2); TLegend* data_reco;
xj_functions::DrawsPHENIXLegend(data_reco, "datareco", 0.10, 0.68, 0.5, 0.93,
                                "p+p #sqrt{s}=200 GeV", false, false, true);
can_pt1pt2_data->cd(3); TLegend* data_unfold;
xj_functions::DrawsPHENIXLegend(data_unfold, "dataunfold", 0.10, 0.68, 0.5, 0.93,
                                "p+p #sqrt{s}=200 GeV", true, false, true);

// Reco comparison
TCanvas* can_xj_data_reco = new TCanvas("can_xj_data_reco", "p_{T1} vs p_{T2}", 700, 800);
xj_functions::plotxj_and_ratio(can_xj_data_reco, h_xj_reco_half_pythia, h_xj_unfold_data_pythia, h_xj_unfold_half_pythia,
                               h_xj_data, true, "x_{J, Data}", "Data/Reco", 0, 1.0, 4, 20, 1, 53, 4, 53, 1, 20);

can_xj_data_reco->cd(1);
TLegend* xj_data_unfold_2;
xj_functions::draw_xj_legend(xj_data_unfold_2, "xj_data_unfold_2", 0.17, 0.45, 0.50, 0.59,
                             h_xj_reco_full_pythia, h_xj_unfold_half_pythia, h_xj_data, h_xj_unfold_data_pythia,
                             "PYTHIA 8 Reco", "PYTHIA 8 Reco Unfold", "Data", "Data Unfolded", true);
TLegend* xj_data_21;
xj_functions::DrawsPHENIXLegend(xj_data_21, "xj_data_21", 0.10, 0.60, 0.5, 0.93,
                                "p+p #sqrt{s}=200 GeV  Runs 47352-47733", false, true, true, true);


// Data Histogram (HERWIG)
TCanvas* can_xj_data_herwig = new TCanvas("can_xj_data_herwig", "p_{T1} vs p_{T2}", 700, 800);

xj_functions::plotxj_and_ratio(can_xj_data_herwig, h_xj_data, h_xj_unfold_data_herwig, h_xj_truth_full_herwig, h_xj_data,
                               false, "x_{J, Data}", "Unfold/Truth", 0, 1.0, 1, 20, 1, 53, kBlue+1, 20, 1, 20);

can_xj_data_herwig->cd(1);
TLegend* xj_data_unfold_herwig;
xj_functions::draw_xj_legend(xj_data_unfold_herwig, "xj_data_unfold_herwig", 0.17, 0.45, 0.50, 0.59,
                             h_xj_truth_full_herwig, h_xj_data, h_xj_unfold_data_herwig, h_xj_unfold_half_herwig,
                             "HERWIG", "Data", "Data Unfolded", " ");
TLegend* xj_data_2_herwig;
xj_functions::DrawsPHENIXLegend(xj_data_2_herwig, "xj_data_2_herwig", 0.10, 0.60, 0.5, 0.93,
                                "p+p #sqrt{s}=200 GeV  Runs 47352-47733", true, false, true, true);

// 2D histogram Data (HERWIG)
TCanvas* can_pt1pt2_data_herwig = new TCanvas("can_pt1pt2_data_herwig", "p_{T1} vs p_{T2}", 2400, 600);
can_pt1pt2_data_herwig->Divide(3, 1);

// Truth
can_pt1pt2_data_herwig->cd(1);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_pt1pt2_truth_herwig, "p_{T1, truth} [GeV]", "p_{T2, truth} [GeV]", 2, 20, 1, 1);
h_pt1pt2_truth_herwig->Draw("Colz");
h_pt1pt2_truth_herwig->SetAxisRange(0, 89, "X");
h_pt1pt2_truth_herwig->SetAxisRange(0, 89, "Y");
h_pt1pt2_truth_herwig->SetAxisRange(pow(10, -9), 1, "Z");
h_pt1pt2_truth_herwig->GetZaxis()->SetLabelSize(0.035);

// Data
can_pt1pt2_data_herwig->cd(2);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_pt1pt2_data, "p_{T1, data} [GeV]", "p_{T2, data} [GeV]", 2, 20, 1, 1);
h_pt1pt2_data->Draw("Colz");
h_pt1pt2_data->SetAxisRange(0, 89, "X");
h_pt1pt2_data->SetAxisRange(0, 89, "Y");
h_pt1pt2_data->SetAxisRange(1, pow(10, 7), "Z");
h_pt1pt2_data->GetZaxis()->SetLabelSize(0.035);

// Unfolded data
can_pt1pt2_data_herwig->cd(3);
gPad->SetRightMargin(0.16);
gPad->SetLogz();
xj_functions::SetHist(h_unfolded_pt1pt2_data_herwig, "p_{T1, unfold} [GeV]", "p_{T2, unfold} [GeV]", 2, 20, 1, 1);
h_unfolded_pt1pt2_data_herwig->Draw("Colz");
h_unfolded_pt1pt2_data_herwig->SetAxisRange(0, 89, "X");
h_unfolded_pt1pt2_data_herwig->SetAxisRange(0, 89, "Y");
h_unfolded_pt1pt2_data_herwig->SetAxisRange(1, pow(10, 7), "Z");
h_unfolded_pt1pt2_data_herwig->GetZaxis()->SetLabelSize(0.035);

// Legends
can_pt1pt2_data_herwig->cd(1); TLegend* data_truth_herwig;
xj_functions::DrawsPHENIXLegend(data_truth_herwig, "datatruth_herwig", 0.10, 0.68, 0.5, 0.93,
                                "p+p #sqrt{s}=200 GeV", true, false, false);
can_pt1pt2_data_herwig->cd(2); TLegend* data_reco_herwig;
xj_functions::DrawsPHENIXLegend(data_reco_herwig, "datareco_herwig", 0.10, 0.68, 0.5, 0.93,
                                "p+p #sqrt{s}=200 GeV", false, false, true);
can_pt1pt2_data_herwig->cd(3); TLegend* data_unfold_herwig;
xj_functions::DrawsPHENIXLegend(data_unfold_herwig, "dataunfold_herwig", 0.10, 0.68, 0.5, 0.93,
                                "p+p #sqrt{s}=200 GeV", true, false, true);

// Reco comparison
TCanvas* can_xj_data_reco_herwig = new TCanvas("can_xj_data_reco_herwig", "p_{T1} vs p_{T2}", 700, 800);
xj_functions::plotxj_and_ratio(can_xj_data_reco_herwig, h_xj_reco_half_herwig, h_xj_unfold_data_herwig, h_xj_unfold_half_herwig,
                               h_xj_data, true, "x_{J, Data}", "Data/Reco", 0, 1.0, 4, 20, 1, 53, 4, 53, 1, 20);

can_xj_data_reco_herwig->cd(1);
TLegend* xj_data_unfold_2_herwig;
xj_functions::draw_xj_legend(xj_data_unfold_2_herwig, "xj_data_unfold_2_herwig", 0.17, 0.45, 0.50, 0.59,
                             h_xj_reco_full_herwig, h_xj_unfold_half_herwig, h_xj_data, h_xj_unfold_data_herwig,
                             "HERWIG Reco", "HERWIG Reco Unfold", "Data", "Data Unfolded", true);
TLegend* xj_data_21_herwig;
xj_functions::DrawsPHENIXLegend(xj_data_21_herwig, "xj_data_21_herwig", 0.10, 0.60, 0.5, 0.93,
                                "p+p #sqrt{s}=200 GeV  Runs 47352-47733", false, true, true, true);


// === Canvas: Half Closure Comparison (PYTHIA vs HERWIG) ===
TCanvas* can_xj_half_compare = new TCanvas("can_xj_half_compare", "x_{J} Half Closure Compare", 700, 800);
xj_functions::plot_xj_and_ratio_compare(
    can_xj_half_compare,
    h_xj_reco_half_pythia, h_xj_reco_half_herwig,
    h_xj_unfold_half_pythia, h_xj_truth_half_pythia,
    h_xj_unfold_half_herwig, h_xj_truth_half_herwig,
    nullptr, false,
    "x_{J}", "Unfolded/Truth", 0, 1.0,
    kBlue+1, 20, kBlue+1, 21,
    kBlack, 24, kGreen+3, 25,
    kViolet+1, 20, kGreen+3, 26, 38, 26
);

// Legend
can_xj_half_compare->cd(1);
TLegend* leg_xj_half_compare;
xj_functions::draw_xj_legend_compare(
    leg_xj_half_compare, "leg_xj_half_compare",
    0.17, 0.45, 0.50, 0.63,
    h_xj_unfold_half_pythia, h_xj_reco_half_pythia, h_xj_truth_half_pythia,
    h_xj_unfold_half_herwig, h_xj_reco_half_herwig, h_xj_truth_half_herwig,nullptr,
    "PYTHIA 8 Unfolded", "PYTHIA 8 Reco", "PYTHIA 8 Truth",
    "HERWIG Unfolded", "HERWIG Reco", "HERWIG Truth"," ",
    false, true
);
TLegend* leg_half_closure_phx;
xj_functions::DrawsPHENIXLegend(leg_half_closure_phx, "leg_half_closure_phx", 0.10, 0.64, 0.5, 0.93,
    "p+p #sqrt{s}=200 GeV", true, true, false, true);
leg_half_closure_phx->AddEntry("", "Half Closure Test", "");


// === Canvas: Data Unfolding Comparison (PYTHIA vs HERWIG) ===
TCanvas* can_xj_data_compare = new TCanvas("can_xj_data_compare", "x_{J} Data Unfolding Compare", 700, 800);
xj_functions::plot_xj_and_ratio_compare(
    can_xj_data_compare,
    h_xj_unfold_data_pythia, h_xj_unfold_data_herwig,
    h_xj_unfold_data_pythia, h_xj_unfold_half_pythia,
    h_xj_unfold_data_herwig, h_xj_unfold_half_herwig,
    h_xj_data, true,
    "x_{J, Data}", "Unfolded/DataReco", 0, 1.0,
    kGreen+1, 20, kBlue+1, 21,
    kBlack, 24, kBlack, 24,
    kViolet+1, 20, kOrange+1, 26
);

// Legend
can_xj_data_compare->cd(1);
TLegend* leg_xj_data_compare;
xj_functions::draw_xj_legend_compare(
    leg_xj_data_compare, "leg_xj_data_compare",
    0.17, 0.45, 0.50, 0.63,
    h_xj_unfold_data_pythia, h_xj_unfold_half_pythia, h_xj_data,
    h_xj_unfold_data_herwig, h_xj_unfold_half_herwig, nullptr,nullptr,
    "PYTHIA 8 Unfolded", "PYTHIA 8 Reco Unfold", "Data",
    " ",
    "HERWIG Unfolded", "HERWIG Reco Unfold", " ",
    false, true
);
TLegend* leg_data_compare_phx;
xj_functions::DrawsPHENIXLegend(leg_data_compare_phx, "leg_data_compare_phx", 0.10, 0.64, 0.5, 0.93,
    "p+p #sqrt{s}=200 GeV  Runs 47352-47733", false, true, true, true);
leg_data_compare_phx->AddEntry("", "Data Unfolding Comparison", "");



// Simple Compare (HERWIG and Pythia)
TCanvas* can_xj_data_herwig_pythia = new TCanvas("can_xj_data_herwig_pythia", "p_{T1} vs p_{T2}", 700, 800);

xj_functions::plotxj_and_ratio(can_xj_data_herwig_pythia, h_xj_data, h_xj_unfold_data_herwig, h_xj_unfold_data_pythia, h_xj_data,
                               false, "x_{J, Data}", "Herwig Unfold/Pythia Unfold", 0, 1.0, 1, 20, 1, 53, kBlue+1, 20, 1, 20);

can_xj_data_herwig_pythia->cd(1);
TLegend* xj_data_unfold_herwig_pythia;
xj_functions::draw_xj_legend(xj_data_unfold_herwig_pythia, "xj_data_unfold_herwig", 0.17, 0.45, 0.50, 0.59,
                             h_xj_data, h_xj_unfold_data_pythia, h_xj_unfold_data_herwig, h_xj_unfold_half_herwig,
			     "Data Not Unfolded",  "Data Unfold w/ Pythia", "Data Unfolded w/ Herwig", " ");
TLegend* xj_data_2_herwig_pythia;
xj_functions::DrawsPHENIXLegend(xj_data_2_herwig_pythia, "xj_data_2_herwig", 0.10, 0.60, 0.5, 0.93,
                                "p+p #sqrt{s}=200 GeV  Runs 47352-47733", true, false, true, true);



}//end function
