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

void nIteration_Unfold_Plot(string infile1 = "Xj_2D_Response_Bins19_R4_Smear_Flatten_Herwig_Reweight_Herwig-Jet10-Run21-multiR.root", string infile2 = "Xj_2D_Response_Bins19_R4_Smear_Flatten_Herwig_Reweight_Herwig-Jet30-Run21-multiR.root", string infile3 ="Xj_2D_Response_Bins19_R4_Smear_Flatten_pythia8-Jet20-Run21-multiR.root", string infile_data = "Xj_2D_Response_Bins19_R4_Data_Flatten_output_run2pp_ana468_2024p012_v001_data_All.Root", const std::string configfile = "binning_original.config", int niterations = 6, bool Pythia = false){

  // --- Load input files ---
  TFile* infile01 = TFile::Open(infile1.c_str());
  TFile* infile02 = TFile::Open(infile2.c_str());
  TFile* infile03 = TFile::Open(infile3.c_str());

  TFile* infiledata = TFile::Open(infile_data.c_str());

  // --- Get histograms from base file (Jet10) ---

  TH1D *h_flat_truth_pt1pt2_full        = (TH1D*)infile01->Get("h_truth_flat_pt1pt2_full");
  TH1D *h_truth_count_flat_pt1pt2_full  = (TH1D*)infile01->Get("h_truth_count_flat_pt1pt2_full");
  TH1D *h_flat_reco_pt1pt2_full         = (TH1D*)infile01->Get("h_reco_flat_pt1pt2_full");
  TH1D *h_count_flat_reco_pt1pt2_full   = (TH1D*)infile01->Get("h_count_reco_flat_pt1pt2_full");
  TH2D *h_flat_response_pt1pt2_full     = (TH2D*)infile01->Get("h_flat_response_pt1pt2_full");
  RooUnfoldResponse *rooResponse_dense_pt1pt2_full = (RooUnfoldResponse*)infile01->Get("rooResponse_dense_pt1pt2_full");

  TH1D *h_flat_truth_pt1pt2_half_fill        = (TH1D*)infile01->Get("h_truth_flat_pt1pt2_half_fill");
  TH1D *h_truth_count_flat_pt1pt2_half_fill  = (TH1D*)infile01->Get("h_truth_count_flat_pt1pt2_half_fill");
  TH1D *h_flat_reco_pt1pt2_half_fill         = (TH1D*)infile01->Get("h_reco_flat_pt1pt2_half_fill");
  TH1D *h_count_flat_reco_pt1pt2_half_fill   = (TH1D*)infile01->Get("h_count_reco_flat_pt1pt2_half_fill");
  TH2D *h_flat_response_pt1pt2_half_fill     = (TH2D*)infile01->Get("h_flat_response_pt1pt2_half_fill");
  RooUnfoldResponse *rooResponse_dense_pt1pt2_half_fill = (RooUnfoldResponse*)infile01->Get("rooResponse_dense_pt1pt2_half_fill");

  TH1D *h_flat_truth_pt1pt2_half_test        = (TH1D*)infile01->Get("h_truth_flat_pt1pt2_half_test");
  TH1D *h_truth_count_flat_pt1pt2_half_test  = (TH1D*)infile01->Get("h_truth_count_flat_pt1pt2_half_test");
  TH1D *h_flat_reco_pt1pt2_half_test         = (TH1D*)infile01->Get("h_reco_flat_pt1pt2_half_test");
  TH1D *h_count_flat_reco_pt1pt2_half_test   = (TH1D*)infile01->Get("h_count_reco_flat_pt1pt2_half_test");
  TH2D *h_flat_response_pt1pt2_half_test     = (TH2D*)infile01->Get("h_flat_response_pt1pt2_half_test");

  TH1D *h_pt_lead_truth_matched = (TH1D*)infile01->Get("h_pt_lead_truth_matched");
  TH1D *h_pt_lead_reco_matched  = (TH1D*)infile01->Get("h_pt_lead_reco_matched");

  // --- Add Jet30 histograms ---
  h_flat_truth_pt1pt2_full       ->Add((TH1D*)infile02->Get("h_truth_flat_pt1pt2_full"));
  h_truth_count_flat_pt1pt2_full ->Add((TH1D*)infile02->Get("h_truth_count_flat_pt1pt2_full"));
  h_flat_reco_pt1pt2_full        ->Add((TH1D*)infile02->Get("h_reco_flat_pt1pt2_full"));
  h_count_flat_reco_pt1pt2_full  ->Add((TH1D*)infile02->Get("h_count_reco_flat_pt1pt2_full"));
  h_flat_response_pt1pt2_full    ->Add((TH2D*)infile02->Get("h_flat_response_pt1pt2_full"));
  rooResponse_dense_pt1pt2_full  ->Add(*(RooUnfoldResponse*)infile02->Get("rooResponse_dense_pt1pt2_full"));

  h_flat_truth_pt1pt2_half_fill       ->Add((TH1D*)infile02->Get("h_truth_flat_pt1pt2_half_fill"));
  h_truth_count_flat_pt1pt2_half_fill ->Add((TH1D*)infile02->Get("h_truth_count_flat_pt1pt2_half_fill"));
  h_flat_reco_pt1pt2_half_fill        ->Add((TH1D*)infile02->Get("h_reco_flat_pt1pt2_half_fill"));
  h_count_flat_reco_pt1pt2_half_fill  ->Add((TH1D*)infile02->Get("h_count_reco_flat_pt1pt2_half_fill"));
  h_flat_response_pt1pt2_half_fill    ->Add((TH2D*)infile02->Get("h_flat_response_pt1pt2_half_fill"));
  rooResponse_dense_pt1pt2_half_fill  ->Add(*(RooUnfoldResponse*)infile02->Get("rooResponse_dense_pt1pt2_half_fill"));

  h_flat_truth_pt1pt2_half_test       ->Add((TH1D*)infile02->Get("h_truth_flat_pt1pt2_half_test"));
  h_truth_count_flat_pt1pt2_half_test ->Add((TH1D*)infile02->Get("h_truth_count_flat_pt1pt2_half_test"));
  h_flat_reco_pt1pt2_half_test        ->Add((TH1D*)infile02->Get("h_reco_flat_pt1pt2_half_test"));
  h_count_flat_reco_pt1pt2_half_test  ->Add((TH1D*)infile02->Get("h_count_reco_flat_pt1pt2_half_test"));
  h_flat_response_pt1pt2_half_test    ->Add((TH2D*)infile02->Get("h_flat_response_pt1pt2_half_test"));

  h_pt_lead_truth_matched ->Add((TH1D*)infile02->Get("h_pt_lead_truth_matched"));
  h_pt_lead_reco_matched  ->Add((TH1D*)infile02->Get("h_pt_lead_reco_matched"));

  if(Pythia){
  // --- Add Jet20 histograms ---
  h_flat_truth_pt1pt2_full       ->Add((TH1D*)infile03->Get("h_truth_flat_pt1pt2_full"));
  h_truth_count_flat_pt1pt2_full ->Add((TH1D*)infile03->Get("h_truth_count_flat_pt1pt2_full"));
  h_flat_reco_pt1pt2_full        ->Add((TH1D*)infile03->Get("h_reco_flat_pt1pt2_full"));
  h_count_flat_reco_pt1pt2_full  ->Add((TH1D*)infile03->Get("h_count_reco_flat_pt1pt2_full"));
  h_flat_response_pt1pt2_full    ->Add((TH2D*)infile03->Get("h_flat_response_pt1pt2_full"));
  rooResponse_dense_pt1pt2_full  ->Add(*(RooUnfoldResponse*)infile03->Get("rooResponse_dense_pt1pt2_full"));

  h_flat_truth_pt1pt2_half_fill       ->Add((TH1D*)infile03->Get("h_truth_flat_pt1pt2_half_fill"));
  h_truth_count_flat_pt1pt2_half_fill ->Add((TH1D*)infile03->Get("h_truth_count_flat_pt1pt2_half_fill"));
  h_flat_reco_pt1pt2_half_fill        ->Add((TH1D*)infile03->Get("h_reco_flat_pt1pt2_half_fill"));
  h_count_flat_reco_pt1pt2_half_fill  ->Add((TH1D*)infile03->Get("h_count_reco_flat_pt1pt2_half_fill"));
  h_flat_response_pt1pt2_half_fill    ->Add((TH2D*)infile03->Get("h_flat_response_pt1pt2_half_fill"));
  rooResponse_dense_pt1pt2_half_fill  ->Add(*(RooUnfoldResponse*)infile03->Get("rooResponse_dense_pt1pt2_half_fill"));

  h_flat_truth_pt1pt2_half_test       ->Add((TH1D*)infile03->Get("h_truth_flat_pt1pt2_half_test"));
  h_truth_count_flat_pt1pt2_half_test ->Add((TH1D*)infile03->Get("h_truth_count_flat_pt1pt2_half_test"));
  h_flat_reco_pt1pt2_half_test        ->Add((TH1D*)infile03->Get("h_reco_flat_pt1pt2_half_test"));
  h_count_flat_reco_pt1pt2_half_test  ->Add((TH1D*)infile03->Get("h_count_reco_flat_pt1pt2_half_test"));
  h_flat_response_pt1pt2_half_test    ->Add((TH2D*)infile03->Get("h_flat_response_pt1pt2_half_test"));

  h_pt_lead_truth_matched ->Add((TH1D*)infile03->Get("h_pt_lead_truth_matched"));
  h_pt_lead_reco_matched  ->Add((TH1D*)infile03->Get("h_pt_lead_reco_matched"));
  }
  //Draw pt

  
 TCanvas* can_ptqa = new TCanvas("can_ptqa", "p_{T1} vs p_{T2}", 800, 800);
  can_ptqa->cd(1);
gPad->SetRightMargin(0.16);
 gPad->SetLogy();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt_lead_truth_matched, "p_{T1, truth} [GeV]", "Counts", 1, 20, 1, 1);
 h_pt_lead_truth_matched->Draw("P");
 h_pt_lead_truth_matched->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 //h_pt_lead_truth_matched->SetAxisRange(pow(10,-9), 1, "Y"); // Adjust Y range if needed
 //h_pt_lead_truth_matched->SetAxisRange(pow(10,-9), 1, "Z"); // Adjust color scale range
 xj_functions::SetHist(h_pt_lead_reco_matched, "p_{T1, truth} [GeV]", "Counts", 2, 20, 1, 1);
 h_pt_lead_reco_matched->Draw("same");
 h_pt_lead_truth_matched->GetZaxis()->SetLabelSize(0.035);
 
 
  //xj histograms
  read_binning rb(configfile.c_str());
  Int_t read_nbins = rb.get_nbins();
  Int_t minentries = rb.get_minentries();
  const int nbins = read_nbins;
  float ipt_bins[nbins+1];
  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
    }



 

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

  //trim matrices

    // FULL skim setup
    int nbins_total_full = nbins * nbins;
    int nrecobins_full = 0, ntruthbins_full = 0;
    std::vector<int> reco_map_full(nbins_total_full), truth_map_full(nbins_total_full);

    TH1D* h_truth_map_full = (TH1D*) h_flat_truth_pt1pt2_full->Clone("h_truth_mapping_full");
    TH1D* h_reco_map_full = (TH1D*) h_flat_reco_pt1pt2_full->Clone("h_reco_mapping_full");
    h_truth_map_full->Reset(); h_reco_map_full->Reset();

    for (int ib = 0; ib < nbins_total_full; ib++) {
      if (h_truth_count_flat_pt1pt2_full->GetBinContent(ib+1) >= minentries)
	truth_map_full[ib] = ++ntruthbins_full;
      else
	truth_map_full[ib] = 0;

      if (h_count_flat_reco_pt1pt2_full->GetBinContent(ib+1) >= minentries)
	reco_map_full[ib] = ++nrecobins_full;
      else
	reco_map_full[ib] = 0;
    }

    TH1D* h_flat_truth_skim_full = new TH1D("h_flat_truth_skim_full", "", ntruthbins_full, 0, ntruthbins_full);
    TH1D* h_flat_reco_skim_full = new TH1D("h_flat_reco_skim_full", "", nrecobins_full, 0, nrecobins_full);
    TH2D* h_flat_response_skim_full = new TH2D("h_flat_response_skim_full", "", nrecobins_full, 0, nrecobins_full, ntruthbins_full, 0, ntruthbins_full);

    for (int ib = 0; ib < nbins_total_full; ib++) {
      int tbin = truth_map_full[ib];
      int rbin = reco_map_full[ib];
      if (tbin) {
	h_flat_truth_skim_full->SetBinContent(tbin, h_flat_truth_pt1pt2_full->GetBinContent(ib+1));
	h_flat_truth_skim_full->SetBinError(tbin, h_flat_truth_pt1pt2_full->GetBinError(ib+1));
      }
      if (rbin) {
	h_flat_reco_skim_full->SetBinContent(rbin, h_flat_reco_pt1pt2_full->GetBinContent(ib+1));
	h_flat_reco_skim_full->SetBinError(rbin, h_flat_reco_pt1pt2_full->GetBinError(ib+1));
      }
      if (tbin) {
	for (int ibr = 0; ibr < nbins_total_full; ibr++) {
	  int rbin2 = reco_map_full[ibr];
	  if (!rbin2) continue;
	  int rbin_orig = h_flat_response_pt1pt2_full->GetBin(ibr+1, ib+1);
	  int rbin_new = h_flat_response_skim_full->GetBin(rbin2, tbin);
	  h_flat_response_skim_full->SetBinContent(rbin_new, h_flat_response_pt1pt2_full->GetBinContent(rbin_orig));
	  h_flat_response_skim_full->SetBinError(rbin_new, h_flat_response_pt1pt2_full->GetBinError(rbin_orig));
	}
      }
    }
    RooUnfoldResponse rooResponse_skim_full(h_flat_reco_skim_full, h_flat_truth_skim_full, h_flat_response_skim_full);

    // COMBINED skim setup
    int nbins_total = nbins * nbins;
    std::vector<int> truth_map(nbins_total, 0);
    std::vector<int> reco_map(nbins_total, 0);
    int ntruthbins = 0, nrecobins = 0;

    // Create combined maps by taking the union of non-empty bins across half_fill and half_test
    for (int ib = 0; ib < nbins_total; ib++) {
      bool truth_filled = (h_truth_count_flat_pt1pt2_half_fill->GetBinContent(ib+1) >= minentries) ||
	(h_truth_count_flat_pt1pt2_half_test->GetBinContent(ib+1) >= minentries);
      bool reco_filled  = (h_count_flat_reco_pt1pt2_half_fill->GetBinContent(ib+1) >= minentries) ||
	(h_count_flat_reco_pt1pt2_half_test->GetBinContent(ib+1) >= minentries);
  
      if (truth_filled) truth_map[ib] = ++ntruthbins;
      if (reco_filled)  reco_map[ib]  = ++nrecobins;
    }

    // Allocate common histograms
    TH1D* h_flat_truth_skim_half_fill = new TH1D("h_flat_truth_skim_half_fill", "", ntruthbins, 0, ntruthbins);
    TH1D* h_flat_reco_skim_half_fill  = new TH1D("h_flat_reco_skim_half_fill",  "", nrecobins, 0, nrecobins);
    TH2D* h_flat_response_skim_half_fill = new TH2D("h_flat_response_skim_half_fill", "", nrecobins, 0, nrecobins, ntruthbins, 0, ntruthbins);

    TH1D* h_flat_truth_skim_half_test = new TH1D("h_flat_truth_skim_half_test", "", ntruthbins, 0, ntruthbins);
    TH1D* h_flat_reco_skim_half_test  = new TH1D("h_flat_reco_skim_half_test",  "", nrecobins, 0, nrecobins);

    // Fill skimmed histograms for half_fill
    for (int ib = 0; ib < nbins_total; ib++) {
      int tbin = truth_map[ib];
      int rbin = reco_map[ib];
      if (tbin) {
	h_flat_truth_skim_half_fill->SetBinContent(tbin, h_flat_truth_pt1pt2_half_fill->GetBinContent(ib+1));
	h_flat_truth_skim_half_fill->SetBinError(tbin, h_flat_truth_pt1pt2_half_fill->GetBinError(ib+1));
      }
      if (rbin) {
	h_flat_reco_skim_half_fill->SetBinContent(rbin, h_flat_reco_pt1pt2_half_fill->GetBinContent(ib+1));
	h_flat_reco_skim_half_fill->SetBinError(rbin, h_flat_reco_pt1pt2_half_fill->GetBinError(ib+1));
      }

      if (tbin) {
	for (int ibr = 0; ibr < nbins_total; ibr++) {
	  int rbin2 = reco_map[ibr];
	  if (!rbin2) continue;
	  int orig_bin = h_flat_response_pt1pt2_half_fill->GetBin(ibr+1, ib+1);
	  int skim_bin = h_flat_response_skim_half_fill->GetBin(rbin2, tbin);
	  h_flat_response_skim_half_fill->SetBinContent(skim_bin, h_flat_response_pt1pt2_half_fill->GetBinContent(orig_bin));
	  h_flat_response_skim_half_fill->SetBinError(skim_bin, h_flat_response_pt1pt2_half_fill->GetBinError(orig_bin));
	}
      }
    }

    // Fill skimmed histograms for half_test using SAME bin maps
    for (int ib = 0; ib < nbins_total; ib++) {
      int tbin = truth_map[ib];
      int rbin = reco_map[ib];
      if (tbin) {
	h_flat_truth_skim_half_test->SetBinContent(tbin, h_flat_truth_pt1pt2_half_test->GetBinContent(ib+1));
	h_flat_truth_skim_half_test->SetBinError(tbin, h_flat_truth_pt1pt2_half_test->GetBinError(ib+1));
      }
      if (rbin) {
	h_flat_reco_skim_half_test->SetBinContent(rbin, h_flat_reco_pt1pt2_half_test->GetBinContent(ib+1));
	h_flat_reco_skim_half_test->SetBinError(rbin, h_flat_reco_pt1pt2_half_test->GetBinError(ib+1));
      }
    }


     //Read in Data File and Set up Histos
 
  TH1D *h_data_flat_pt1pt2  = (TH1D*)infiledata->Get("h_data_flat_pt1pt2");
  //TH1D *h_count_data_flat_pt1pt2_full  = (TH1D*)infiledata->Get("h_count_data_flat_pt1pt2_full");
  
 

  h_data_flat_pt1pt2->Scale(.5);
  TH1D* h_flat_data_skim = (TH1D*)h_flat_reco_skim_full->Clone();
  h_flat_data_skim->Reset();

  int nbins_total_data = reco_map_full.size();  // should be nbins * nbins

for (int ib = 0; ib < nbins_total_data; ib++) {
  int new_bin = reco_map_full[ib];
  if (new_bin > 0) {
    h_flat_data_skim->SetBinContent(new_bin, h_data_flat_pt1pt2->GetBinContent(ib + 1));
    h_flat_data_skim->SetBinError(new_bin, h_data_flat_pt1pt2->GetBinError(ib + 1));
  }
}


  //skim Data
  


    // FULL set
    std::cout << "FULL SET:" << std::endl;
    std::cout << "  Skimmed Reco Bins:  " << h_flat_reco_skim_full->GetNbinsX() << std::endl;
    std::cout << "  Skimmed Truth Bins: " << h_flat_truth_skim_full->GetNbinsX() << std::endl;
    std::cout << "  Response Matrix Bins: " 
	      << h_flat_response_skim_full->GetXaxis()->GetNbins() << " (Reco) × "
	      << h_flat_response_skim_full->GetYaxis()->GetNbins() << " (Truth)" << std::endl;

    // DATA histogram dimension check
    if (h_flat_data_skim->GetNbinsX() != h_flat_reco_skim_full->GetNbinsX()) {
      std::cerr << "  [WARNING] h_flat_data_skim has " << h_flat_data_skim->GetNbinsX()
		<< " bins, but reco skim has " << h_flat_reco_skim_full->GetNbinsX() << " bins!" << std::endl;
    } else {
      std::cout << "  Data Histogram: h_flat_data_skim has correct dimensions (" 
		<< h_flat_data_skim->GetNbinsX() << " bins)" << std::endl;
    }


    // HALF-FILL set
    std::cout << "\nHALF-FILL SET:" << std::endl;
    std::cout << "  Skimmed Reco Bins:  " << h_flat_reco_skim_half_fill->GetNbinsX() << std::endl;
    std::cout << "  Skimmed Truth Bins: " << h_flat_truth_skim_half_fill->GetNbinsX() << std::endl;
    std::cout << "  Response Matrix Bins: " 
	      << h_flat_response_skim_half_fill->GetXaxis()->GetNbins() << " (Reco) × "
	      << h_flat_response_skim_half_fill->GetYaxis()->GetNbins() << " (Truth)" << std::endl;

    // HALF-TEST set (no response matrix)
    std::cout << "\nHALF-TEST SET:" << std::endl;
    std::cout << "  Skimmed Reco Bins:  " << h_flat_reco_skim_half_test->GetNbinsX() << std::endl;
    std::cout << "  Skimmed Truth Bins: " << h_flat_truth_skim_half_test->GetNbinsX() << std::endl;
    std::cout << "  Response Matrix: Not constructed for half-test set." << std::endl;

 
    RooUnfoldResponse roo_response_flat_skim_full(h_flat_reco_skim_full,h_flat_truth_skim_full,h_flat_response_skim_full); RooUnfoldResponse roo_response_flat_skim_half_fill(h_flat_reco_skim_half_fill,h_flat_truth_skim_half_fill,h_flat_response_skim_half_fill);

     TH1D* h_flat_unfold_skim_full; TH1D* h_flat_unfold_skim_half;  TH1D* h_flat_unfold_skim_data;
      TH1D* h_flat_unfold_pt1pt2_full; TH1D* h_flat_unfold_pt1pt2_half; TH1D* h_flat_unfold_pt1pt2_data;

    if (minentries){//unfold skimmed histos
      h_flat_reco_skim_full->Scale(.5);h_flat_reco_skim_half_test->Scale(.5);
     

      RooUnfoldBayes *unfoldfull = new RooUnfoldBayes(&roo_response_flat_skim_full,h_flat_reco_skim_full,niterations);
      RooUnfoldBayes *unfoldhalf = new RooUnfoldBayes(&roo_response_flat_skim_half_fill, h_flat_reco_skim_half_test, niterations);
      RooUnfoldBayes *unfolddata = new RooUnfoldBayes(&roo_response_flat_skim_full, h_flat_data_skim, niterations);
      h_flat_unfold_skim_full = (TH1D*) unfoldfull->Hunfold(); h_flat_unfold_skim_half = (TH1D*) unfoldhalf->Hunfold();
      h_flat_unfold_skim_data = (TH1D*) unfolddata->Hunfold();

      h_flat_unfold_pt1pt2_full = (TH1D*)h_flat_truth_pt1pt2_full->Clone();
      h_flat_unfold_pt1pt2_data = (TH1D*)h_data_flat_pt1pt2->Clone();
      
      h_flat_unfold_pt1pt2_half = (TH1D*)h_flat_truth_pt1pt2_half_test->Clone();
      h_flat_unfold_pt1pt2_full->Reset();h_flat_unfold_pt1pt2_half->Reset(); h_flat_unfold_pt1pt2_data->Reset();

      for (int ib = 0; ib < nbins*nbins; ib++){

	int fullbin = truth_map_full[ib];
	int halfbin = truth_map[ib];

	if (fullbin){

	  h_flat_unfold_pt1pt2_full->SetBinContent(ib+1, h_flat_unfold_skim_full->GetBinContent(fullbin));
	  h_flat_unfold_pt1pt2_full->SetBinError(ib+1, h_flat_unfold_skim_full->GetBinError(fullbin));
	  h_flat_unfold_pt1pt2_data->SetBinContent(ib+1, h_flat_unfold_skim_data->GetBinContent(fullbin));
	  h_flat_unfold_pt1pt2_data->SetBinError(ib+1, h_flat_unfold_skim_data->GetBinError(fullbin));
	}

       if (halfbin){

	  h_flat_unfold_pt1pt2_half->SetBinContent(ib+1, h_flat_unfold_skim_half->GetBinContent(halfbin));
	  h_flat_unfold_pt1pt2_half->SetBinError(ib+1, h_flat_unfold_skim_half->GetBinError(halfbin));
	}

      }


    }//


    h_flat_truth_pt1pt2_full->Scale(.5); h_flat_truth_pt1pt2_half_test->Scale(.5);
    h_flat_reco_pt1pt2_full->Scale(.5); h_flat_reco_pt1pt2_half_test->Scale(.5);

    
    
    /*  TCanvas* test2 = new TCanvas("test2", "p_{T1} vs p_{T2}", 800, 800);
     h_flat_unfold_pt1pt2_half->Draw();*/

    // create  pt1pt2 hisograms
  

   TH2D *h_pt1pt2_reco_full = new TH2D("h_pt1pt2_reco_full", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
   TH2D *h_pt1pt2_truth_full = new TH2D("h_pt1pt2_truth_full", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
   TH2D *h_pt1pt2_unfold_full = new TH2D("h_pt1pt2_unfold_full", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
   TH2D *h_pt1pt2_reco_half = new TH2D("h_pt1pt2_reco_half", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
   TH2D *h_pt1pt2_truth_half = new TH2D("h_pt1pt2_truth_half", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
   TH2D *h_pt1pt2_unfold_half = new TH2D("h_pt1pt2_unfold_half", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);

   TH2D *h_pt1pt2_data = new TH2D("h_pt1pt2_data", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
   TH2D *h_pt1pt2_unfold_data = new TH2D("h_pt1pt2_unfold_data", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);

  //unflatten histograms
   xj_functions::make_sym_pt1pt2(h_flat_truth_pt1pt2_full,h_pt1pt2_truth_full,nbins);
   xj_functions::make_sym_pt1pt2(h_flat_reco_pt1pt2_full,h_pt1pt2_reco_full,nbins);
   xj_functions::make_sym_pt1pt2(h_flat_unfold_pt1pt2_full,h_pt1pt2_unfold_full,nbins);
   xj_functions::make_sym_pt1pt2(h_flat_truth_pt1pt2_half_test,h_pt1pt2_truth_half,nbins);
   xj_functions::make_sym_pt1pt2(h_flat_reco_pt1pt2_half_test,h_pt1pt2_reco_half,nbins);
   xj_functions::make_sym_pt1pt2(h_flat_unfold_pt1pt2_half,h_pt1pt2_unfold_half,nbins);

   xj_functions::make_sym_pt1pt2(h_data_flat_pt1pt2,h_pt1pt2_data,nbins);
   xj_functions::make_sym_pt1pt2(h_flat_unfold_pt1pt2_data,h_pt1pt2_unfold_data,nbins);
   /*  TCanvas* test3          = new TCanvas("test3", "p_{T1} vs p_{T2}", 800, 800);
    gPad->SetLogz();
    h_pt1pt2_truth_full->Draw("Colz");
    h_pt1pt2_truth_full->SetAxisRange(1,pow(10,8), "Z"); h_pt1pt2_truth_full->SetAxisRange(0,70, "Y");
    h_pt1pt2_truth_full->SetAxisRange(0,70, "X");
    

     TCanvas* test4 = new TCanvas("test4", "p_{T1} vs p_{T2}", 800, 800);
         gPad->SetLogz();
    h_pt1pt2_unfold_full->Draw("Colz");h_pt1pt2_unfold_full->SetAxisRange(1,pow(10,8), "Z");
    h_pt1pt2_unfold_full->SetAxisRange(0,70, "Y");
    h_pt1pt2_unfold_full->SetAxisRange(0,70, "X");*/

    SetsPhenixStyle();
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

   TCanvas* test          = new TCanvas("test", "p_{T1} vs p_{T2}", 800, 800);
    gPad->SetLogz();
    gPad->SetRightMargin(0.16);
    h_flat_response_skim_full->Draw("ColZ");
    h_flat_response_skim_full->GetXaxis()->SetTitle("p_{T,1}^{reco,bin} x nbins + p_{T,2}^{reco,bin}");
     h_flat_response_skim_full->GetYaxis()->SetTitle("p_{T,1}^{truth,bin} x nbins + p_{T,2}^{truth,bin}");


  //2D histograms
  TCanvas* can_pt1pt2 = new TCanvas("can_pt1pt2", "p_{T1} vs p_{T2}", 2400, 1600);
  can_pt1pt2->Divide(3,2);
  can_pt1pt2->cd(1);
gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt1pt2_truth_full, "p_{T1, truth} [GeV]", "p_{T2, truth} [GeV]", 2, 20, 1, 1);
 h_pt1pt2_truth_full->Draw("Colz");
 h_pt1pt2_truth_full->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_pt1pt2_truth_full->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
 h_pt1pt2_truth_full->SetAxisRange(1, pow(10,8), "Z"); // Adjust color scale range
 h_pt1pt2_truth_full->GetZaxis()->SetLabelSize(0.035);
 

 can_pt1pt2->cd(2);
gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt1pt2_reco_full, "p_{T1, reco} [GeV]", "p_{T2, reco} [GeV]", 2, 20, 1, 1);
 h_pt1pt2_reco_full->Draw("Colz");
 h_pt1pt2_reco_full->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_pt1pt2_reco_full->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
 h_pt1pt2_reco_full->SetAxisRange(1, pow(10,8), "Z");
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
 h_pt1pt2_reco_full->GetZaxis()->SetLabelSize(0.035);

 can_pt1pt2->cd(3);
 gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt1pt2_unfold_full, "p_{T1, unfold} [GeV]", "p_{T2, unfold} [GeV]", 2, 20, 1, 1);
 h_pt1pt2_unfold_full->Draw("Colz");
 h_pt1pt2_unfold_full->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_pt1pt2_unfold_full->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
 h_pt1pt2_unfold_full->SetAxisRange(1, pow(10,8), "Z");
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
 h_pt1pt2_unfold_full->GetZaxis()->SetLabelSize(0.035);

  cout << "made it here 2" << endl;

 can_pt1pt2->cd(4);
gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt1pt2_truth_half, "p_{T1, truth_half} [GeV]", "p_{T2, truth_half} [GeV]", 2, 20, 1, 1);
 h_pt1pt2_truth_half->Draw("Colz");
 h_pt1pt2_truth_half->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_pt1pt2_truth_half->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
  h_pt1pt2_truth_half->SetAxisRange(1, pow(10,8), "Z");
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
 h_pt1pt2_truth_half->GetZaxis()->SetLabelSize(0.035);

 can_pt1pt2->cd(5);
gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt1pt2_reco_half, "p_{T1, reco_half_test} [GeV]", "p_{T2, reco_half_test} [GeV]", 2, 20, 1, 1);
 h_pt1pt2_reco_half->Draw("Colz");
 h_pt1pt2_reco_half->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_pt1pt2_reco_half->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
  h_pt1pt2_reco_half->SetAxisRange(1, pow(10,8), "Z");
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
 h_pt1pt2_reco_half->GetZaxis()->SetLabelSize(0.035);

 can_pt1pt2->cd(6);
 gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt1pt2_unfold_half, "p_{T1, unfold half} [GeV]", "p_{T2, unfold half} [GeV]", 2, 20, 1, 1);
 h_pt1pt2_unfold_half->Draw("Colz");
 h_pt1pt2_unfold_half->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_pt1pt2_unfold_half->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
 h_pt1pt2_unfold_half->SetAxisRange(1, pow(10,8), "Z");
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
 h_pt1pt2_unfold_half->GetZaxis()->SetLabelSize(0.035);


 //2D Hists legends
 can_pt1pt2->cd(1); TLegend* full_truth;
 xj_functions::DrawsPHENIXLegend(full_truth, "fulltruth", 0.10, 0.68, 0.5, 0.93,"p+p #sqrt{s}=200 GeV", true,false,false,false);
 can_pt1pt2->cd(2); TLegend* full_reco;
 xj_functions::DrawsPHENIXLegend(full_reco, "fullreco", 0.10, 0.68, 0.5, 0.93,"p+p #sqrt{s}=200 GeV", false,true,false,false);
 can_pt1pt2->cd(3);TLegend* full_unfold;
 xj_functions::DrawsPHENIXLegend(full_unfold, "fullunfold", 0.10, 0.68, 0.5, 0.93,"p+p #sqrt{s}=200 GeV", true,true,false,false);
can_pt1pt2->cd(4); TLegend* half_truth;
 xj_functions::DrawsPHENIXLegend(half_truth, "halftruth", 0.10, 0.68, 0.5, 0.93,"p+p #sqrt{s}=200 GeV", true,false,false,false);
 can_pt1pt2->cd(5); TLegend* half_reco;
 xj_functions::DrawsPHENIXLegend(half_reco, "halfreco", 0.10, 0.68, 0.5, 0.93,"p+p #sqrt{s}=200 GeV", false,true,false,false);
 can_pt1pt2->cd(6);TLegend* half_unfold;
 xj_functions::DrawsPHENIXLegend(half_unfold, "halfunfold", 0.10, 0.68, 0.5, 0.93,"p+p #sqrt{s}=200 GeV", true,true,false,false);



 // create xJ hists



 TH1D* h_xj_reco_full    = new TH1D("h_xj_reco_full_20_60", ";x_{J};", nbins, ixj_bins);
TH1D* h_xj_truth_full   = new TH1D("h_xj_truth_full_20_60", ";x_{J};", nbins, ixj_bins);
TH1D* h_xj_unfold_full  = new TH1D("h_xj_unfold_full_20_60", ";x_{J};", nbins, ixj_bins);
TH1D* h_xj_reco_half    = new TH1D("h_xj_reco_half_20_60", ";x_{J};", nbins, ixj_bins);
TH1D* h_xj_truth_half   = new TH1D("h_xj_truth_half_20_60", ";x_{J};", nbins, ixj_bins);
TH1D* h_xj_unfold_half  = new TH1D("h_xj_unfold_half_20_60", ";x_{J};", nbins, ixj_bins);
TH1D* h_xj_data         = new TH1D("h_xj_data_20_60", ";x_{J};", nbins, ixj_bins);
TH1D* h_xj_unfold_data  = new TH1D("h_xj_unfold_data_20_60", ";x_{J};", nbins, ixj_bins);

 TH1D* h_xj_reco_full_20_30   = new TH1D("h_xj_reco_full_20_30", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_truth_full_20_30  = new TH1D("h_xj_truth_full_20_30", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_unfold_full_20_30 = new TH1D("h_xj_unfold_full_20_30", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_reco_half_20_30   = new TH1D("h_xj_reco_half_20_30", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_truth_half_20_30  = new TH1D("h_xj_truth_half_20_30", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_unfold_half_20_30 = new TH1D("h_xj_unfold_half_20_30", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_data_20_30        = new TH1D("h_xj_data_20_30", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_unfold_data_20_30 = new TH1D("h_xj_unfold_data_20_30", ";x_{J};", nbins, ixj_bins);

 TH1D* h_xj_reco_full_30_40   = new TH1D("h_xj_reco_full_30_40", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_truth_full_30_40  = new TH1D("h_xj_truth_full_30_40", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_unfold_full_30_40 = new TH1D("h_xj_unfold_full_30_40", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_reco_half_30_40   = new TH1D("h_xj_reco_half_30_40", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_truth_half_30_40  = new TH1D("h_xj_truth_half_30_40", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_unfold_half_30_40 = new TH1D("h_xj_unfold_half_30_40", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_data_30_40        = new TH1D("h_xj_data_30_40", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_unfold_data_30_40 = new TH1D("h_xj_unfold_data_30_40", ";x_{J};", nbins, ixj_bins);

 TH1D* h_xj_reco_full_40_60   = new TH1D("h_xj_reco_full_40_60", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_truth_full_40_60  = new TH1D("h_xj_truth_full_40_60", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_unfold_full_40_60 = new TH1D("h_xj_unfold_full_40_60", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_reco_half_40_60   = new TH1D("h_xj_reco_half_40_60", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_truth_half_40_60  = new TH1D("h_xj_truth_half_40_60", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_unfold_half_40_60 = new TH1D("h_xj_unfold_half_40_60", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_data_40_60        = new TH1D("h_xj_data_40_60", ";x_{J};", nbins, ixj_bins);
 TH1D* h_xj_unfold_data_40_60 = new TH1D("h_xj_unfold_data_40_60", ";x_{J};", nbins, ixj_bins);
 

  // bin check base
 cout << "base leading bin range: " <<  ipt_bins[measure_leading_bin + 0] << " - " <<  ipt_bins[nbins - 1] << endl;
  cout << "base subleading bin range: " <<  ipt_bins[measure_subleading_bin] << " - " <<  ipt_bins[nbins - 1] << endl;
 
  //project xj
  xj_functions::project_xj(h_pt1pt2_reco_full, h_xj_reco_full, nbins, measure_leading_bin + 0, nbins - 1, measure_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_truth_full, h_xj_truth_full, nbins, measure_leading_bin  + 0, nbins - 1, measure_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_unfold_full, h_xj_unfold_full, nbins, measure_leading_bin  + 0, nbins - 1, measure_subleading_bin, nbins - 1);
   xj_functions::project_xj(h_pt1pt2_reco_half, h_xj_reco_half, nbins, measure_leading_bin  + 0, nbins - 1, measure_subleading_bin, nbins - 1);
   xj_functions::project_xj(h_pt1pt2_truth_half, h_xj_truth_half, nbins, measure_leading_bin  + 0, nbins - 1, measure_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_unfold_half, h_xj_unfold_half, nbins, measure_leading_bin  + 0, nbins - 1, measure_subleading_bin, nbins - 1);

  xj_functions::project_xj(h_pt1pt2_unfold_data, h_xj_unfold_data, nbins, measure_leading_bin + 0, nbins - 1, measure_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_data, h_xj_data, nbins, measure_leading_bin + 0, nbins - 1, measure_subleading_bin, nbins - 1);

   //bins check
  cout << "20 - 30 leading bin range: " <<  ipt_bins[measure_leading_bin + 0] << " - " <<  ipt_bins[nbins - 6] << endl;
  cout << "20 - 30 subleading bin range: " <<  ipt_bins[measure_subleading_bin] << " - " <<  ipt_bins[nbins - 6] << endl;

// Projection calls
  xj_functions::project_xj(h_pt1pt2_reco_full,   h_xj_reco_full_20_30,   nbins, measure_leading_bin + 0, nbins - 6, measure_subleading_bin, nbins - 6);
  xj_functions::project_xj(h_pt1pt2_truth_full,  h_xj_truth_full_20_30,  nbins, measure_leading_bin + 0, nbins - 6, measure_subleading_bin, nbins - 6);
  xj_functions::project_xj(h_pt1pt2_unfold_full, h_xj_unfold_full_20_30, nbins, measure_leading_bin + 0, nbins - 6, measure_subleading_bin, nbins - 6);
  xj_functions::project_xj(h_pt1pt2_reco_half,   h_xj_reco_half_20_30,   nbins, measure_leading_bin + 0, nbins - 6, measure_subleading_bin, nbins - 6);
  xj_functions::project_xj(h_pt1pt2_truth_half,  h_xj_truth_half_20_30,  nbins, measure_leading_bin + 0, nbins - 6, measure_subleading_bin, nbins - 6);
  xj_functions::project_xj(h_pt1pt2_unfold_half, h_xj_unfold_half_20_30, nbins, measure_leading_bin + 0, nbins - 6, measure_subleading_bin, nbins - 6);
  xj_functions::project_xj(h_pt1pt2_unfold_data, h_xj_unfold_data_20_30, nbins, measure_leading_bin + 0, nbins - 6, measure_subleading_bin, nbins - 6);
  xj_functions::project_xj(h_pt1pt2_data,        h_xj_data_20_30,        nbins, measure_leading_bin + 0, nbins - 6, measure_subleading_bin, nbins - 6);

  //bins check
  cout << "30 - 40 leading bin range: " <<  ipt_bins[measure_leading_bin + 3] << " - " <<  ipt_bins[nbins - 4] << endl;
  cout << "30 - 40 subleading bin range: " <<  ipt_bins[measure_subleading_bin] << " - " <<  ipt_bins[nbins - 4] << endl;

  // Projection calls
  xj_functions::project_xj(h_pt1pt2_reco_full,   h_xj_reco_full_30_40,   nbins, measure_leading_bin + 3, nbins - 4, measure_subleading_bin, nbins - 4);
  xj_functions::project_xj(h_pt1pt2_truth_full,  h_xj_truth_full_30_40,  nbins, measure_leading_bin + 3, nbins - 4, measure_subleading_bin, nbins - 4);
  xj_functions::project_xj(h_pt1pt2_unfold_full, h_xj_unfold_full_30_40, nbins, measure_leading_bin + 3, nbins - 4, measure_subleading_bin, nbins - 4);
  xj_functions::project_xj(h_pt1pt2_reco_half,   h_xj_reco_half_30_40,   nbins, measure_leading_bin + 3, nbins - 4, measure_subleading_bin, nbins - 4);
  xj_functions::project_xj(h_pt1pt2_truth_half,  h_xj_truth_half_30_40,  nbins, measure_leading_bin + 3, nbins - 4, measure_subleading_bin, nbins - 4);
  xj_functions::project_xj(h_pt1pt2_unfold_half, h_xj_unfold_half_30_40, nbins, measure_leading_bin + 3, nbins - 4, measure_subleading_bin, nbins - 4);
  xj_functions::project_xj(h_pt1pt2_unfold_data, h_xj_unfold_data_30_40, nbins, measure_leading_bin + 3, nbins - 4, measure_subleading_bin, nbins - 4);
  xj_functions::project_xj(h_pt1pt2_data,        h_xj_data_30_40,        nbins, measure_leading_bin + 3, nbins - 4, measure_subleading_bin, nbins - 4);

  //bins check
  cout << "40 - 60 leading bin range: " <<  ipt_bins[measure_leading_bin + 5] << " - " <<  ipt_bins[nbins - 1] << endl;
  cout << "40 - 60 subleading bin range: " <<  ipt_bins[measure_subleading_bin] << " - " <<  ipt_bins[nbins - 1] << endl;

  // Projection calls
  xj_functions::project_xj(h_pt1pt2_reco_full,   h_xj_reco_full_40_60,   nbins, measure_leading_bin + 5, nbins - 1, measure_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_truth_full,  h_xj_truth_full_40_60,  nbins, measure_leading_bin + 5, nbins - 1, measure_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_unfold_full, h_xj_unfold_full_40_60, nbins, measure_leading_bin + 5, nbins - 1, measure_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_reco_half,   h_xj_reco_half_40_60,   nbins, measure_leading_bin + 5, nbins - 1, measure_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_truth_half,  h_xj_truth_half_40_60,  nbins, measure_leading_bin + 5, nbins - 1, measure_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_unfold_half, h_xj_unfold_half_40_60, nbins, measure_leading_bin + 5, nbins - 1, measure_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_unfold_data, h_xj_unfold_data_40_60, nbins, measure_leading_bin + 5, nbins - 1, measure_subleading_bin, nbins - 1);
  xj_functions::project_xj(h_pt1pt2_data,        h_xj_data_40_60,        nbins, measure_leading_bin + 5, nbins - 1, measure_subleading_bin, nbins - 1);
  
 
  //normalize xj histos
  xj_functions::normalize_histo(h_xj_truth_full, nbins);
  xj_functions::normalize_histo(h_xj_reco_full, nbins);
  xj_functions::normalize_histo(h_xj_truth_half, nbins);
  xj_functions::normalize_histo(h_xj_reco_half, nbins);
  xj_functions::normalize_histo(h_xj_unfold_full, nbins);
  xj_functions::normalize_histo(h_xj_unfold_half, nbins);
   xj_functions::normalize_histo(h_xj_unfold_data, nbins);
  xj_functions::normalize_histo(h_xj_data, nbins);

    xj_functions::normalize_histo(h_xj_data_20_30, nbins);
  xj_functions::normalize_histo(h_xj_unfold_data_20_30, nbins);
  xj_functions::normalize_histo(h_xj_truth_full_20_30, nbins);
  xj_functions::normalize_histo(h_xj_reco_full_20_30, nbins);
  xj_functions::normalize_histo(h_xj_truth_half_20_30, nbins);
  xj_functions::normalize_histo(h_xj_reco_half_20_30, nbins);
  xj_functions::normalize_histo(h_xj_unfold_full_20_30, nbins);
  xj_functions::normalize_histo(h_xj_unfold_half_20_30, nbins);

  xj_functions::normalize_histo(h_xj_data_30_40, nbins);
  xj_functions::normalize_histo(h_xj_unfold_data_30_40, nbins);
  xj_functions::normalize_histo(h_xj_truth_full_30_40, nbins);
  xj_functions::normalize_histo(h_xj_reco_full_30_40, nbins);
  xj_functions::normalize_histo(h_xj_truth_half_30_40, nbins);
  xj_functions::normalize_histo(h_xj_reco_half_30_40, nbins);
  xj_functions::normalize_histo(h_xj_unfold_full_30_40, nbins);
  xj_functions::normalize_histo(h_xj_unfold_half_30_40, nbins);

  xj_functions::normalize_histo(h_xj_data_40_60, nbins);
  xj_functions::normalize_histo(h_xj_unfold_data_40_60, nbins);
  xj_functions::normalize_histo(h_xj_truth_full_40_60, nbins);
  xj_functions::normalize_histo(h_xj_reco_full_40_60, nbins);
  xj_functions::normalize_histo(h_xj_truth_half_40_60, nbins);
  xj_functions::normalize_histo(h_xj_reco_half_40_60, nbins);
  xj_functions::normalize_histo(h_xj_unfold_full_40_60, nbins);
  xj_functions::normalize_histo(h_xj_unfold_half_40_60, nbins);

  // xj_functions::normalize_histo(h_xj_classical_truth, nbins);
  //xj_functions::normalize_histo(h_xj_classical_reco, nbins);

  
  //draw histograms

  /*
 TCanvas* can_xj_projected = new TCanvas("can_xj_projected", "p_{T1} vs p_{T2}", 700, 800);
 xj_functions::plotxj_and_ratio(can_xj_projected,h_xj_classical_reco, h_xj_classical_truth,h_xj_truth_full, h_xj_reco_full,true, "x_{J}","Classical/Truth",0, 1.0, 4, 53, 2, 20,  kGreen+1,  20,  8, 53);
  */
  const char* SimType = "PYTHIA 8"; const char* SimTypeReco = "PYTHIA 8 Reco"; const char* SimTypeRecoUnfolded = "PYTHIA 8 Reco Unfold";
  if(!Pythia){
    SimType = "HERWIG 7"; SimTypeReco = "HERWIG 7 Reco"; SimTypeRecoUnfolded = "Herwig 7 Reco Unfold";
  }
 TCanvas* can_xj_full = new TCanvas("can_xj_full", "p_{T1} vs p_{T2}", 700, 800);
 xj_functions::plotxj_and_ratio(can_xj_full,h_xj_reco_full, h_xj_unfold_full,h_xj_truth_full, h_xj_reco_full,false, "x_{J}","Unfold/Truth",0, 1.0, 4, 20, 4, 53,  kGreen+1,  20,  8, 53);
 TLegend* xj_full_unfold;
 xj_functions::draw_xj_legend(xj_full_unfold,"xj_full_unfold", 0.17, 0.45, 0.5, 0.59, h_xj_truth_full,h_xj_reco_full, h_xj_unfold_full, h_xj_unfold_full,SimType,SimTypeReco, "Reco Unfolded", " ");
 TLegend* xj_full_2;
 xj_functions::DrawsPHENIXLegend(xj_full_2, "xj_full_2", 0.10, 0.60, 0.5, 0.93,"p+p #sqrt{s}=200 GeV", true,true,false, true);xj_full_2->AddEntry("","Full Closure Test","");

TCanvas* can_xj_half = new TCanvas("can_xj_half", "p_{T1} vs p_{T2}", 700, 800);
 xj_functions::plotxj_and_ratio(can_xj_half,h_xj_reco_half, h_xj_unfold_half,h_xj_truth_half, h_xj_reco_full,false, "x_{J}","Unfold/Truth",0, 1.0, 4, 20, 4, 53,  kGreen+1,  20,  8, 53);
 TLegend* xj_half_unfold;
  xj_functions::draw_xj_legend(xj_half_unfold,"xj_half_unfold", 0.17, 0.45, 0.50, 0.59, h_xj_truth_half, h_xj_reco_half, h_xj_unfold_half, h_xj_unfold_half,SimType,SimTypeReco, "Reco Unfolded", " ");
TLegend* xj_half_2;
 xj_functions::DrawsPHENIXLegend(xj_half_2, "xj_half_2", 0.10, 0.60, 0.5, 0.93,"p+p #sqrt{s}=200 GeV", true,true,false,true);xj_half_2->AddEntry("","Half Closure Test","");


 
   //Data Histogram
  TCanvas* can_xj_data = new TCanvas("can_xj_data", "p_{T1} vs p_{T2}", 700, 800);

   // Reduce bottom margin to minimize space between plots
  xj_functions::plotxj_and_ratio(can_xj_data, h_xj_data, h_xj_unfold_data, h_xj_truth_full,h_xj_data,false, "x_{J, Data}","Unfold/Truth",0, 1.0,1 , 20, 1, 53,  kGreen+1,  20,  1, 20);

  can_xj_data->cd(1);
   TLegend* xj_data_unfold;
  xj_functions::draw_xj_legend(xj_data_unfold,"xj_data_unfold", 0.17, 0.45, 0.50, 0.59, h_xj_truth_full, h_xj_data, h_xj_unfold_data, h_xj_unfold_half,SimType,"Data", "Data Unfolded", " ");
TLegend* xj_data_2;
 xj_functions::DrawsPHENIXLegend(xj_data_2, "xj_data_2", 0.10, 0.60, 0.5, 0.93,"p+p #sqrt{s}=200 GeV  Runs 47352-47733", true,false,true,true);

   //2D histogram Data

  TCanvas* can_pt1pt2_data = new TCanvas("can_pt1pt2_data", "p_{T1} vs p_{T2}", 2400, 600);
  can_pt1pt2_data->Divide(3,1);
  can_pt1pt2_data->cd(1);
  gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt1pt2_truth_full, "p_{T1, truth} [GeV]", "p_{T2, truth} [GeV]", 2, 20, 1, 1);
 h_pt1pt2_truth_full->Draw("Colz");
 h_pt1pt2_truth_full->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_pt1pt2_truth_full->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
 h_pt1pt2_truth_full->SetAxisRange( 1,pow(10,7), "Z");
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
 h_pt1pt2_truth_full->GetZaxis()->SetLabelSize(0.035);

 can_pt1pt2_data->cd(2);
gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt1pt2_data, "p_{T1, data} [GeV]", "p_{T2, data} [GeV]", 2, 20, 1, 1);
 h_pt1pt2_data->Draw("Colz");
 h_pt1pt2_data->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_pt1pt2_data->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
 h_pt1pt2_data->SetAxisRange(1, pow(10,7), "Z");
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
 h_pt1pt2_data->GetZaxis()->SetLabelSize(0.035);

 can_pt1pt2_data->cd(3);
 gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt1pt2_unfold_data, "p_{T1, unfold} [GeV]", "p_{T2, unfold} [GeV]", 2, 20, 1, 1);
 h_pt1pt2_unfold_data->Draw("Colz");
 h_pt1pt2_unfold_data->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_pt1pt2_unfold_data->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
 h_pt1pt2_unfold_data->SetAxisRange(1, pow(10,7), "Z");
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
 h_pt1pt2_unfold_data->GetZaxis()->SetLabelSize(0.035);

 can_pt1pt2_data->cd(1); TLegend* data_truth;
 xj_functions::DrawsPHENIXLegend(data_truth, "datatruth", 0.10, 0.68, 0.5, 0.93,"p+p #sqrt{s}=200 GeV", true,false,false);
 can_pt1pt2_data->cd(2); TLegend* data_reco;
 xj_functions::DrawsPHENIXLegend(data_reco, "datareco", 0.10, 0.68, 0.5, 0.93,"p+p #sqrt{s}=200 GeV", false,false,true);
 can_pt1pt2_data->cd(3);TLegend* data_unfold;
 xj_functions::DrawsPHENIXLegend(data_unfold, "dataunfold", 0.10, 0.68, 0.5, 0.93,"p+p #sqrt{s}=200 GeV", true,false,true);

 
  TCanvas* can_xj_data_reco = new TCanvas("can_xj_data_reco", "p_{T1} vs p_{T2}", 700, 800);
  xj_functions::plotxj_and_ratio(can_xj_data_reco, h_xj_reco_half, h_xj_unfold_data, h_xj_unfold_half,h_xj_data,true, "x_{J, Data}","Data/Reco",0, 1.0, 4, 20, 1, 53,  4,  53,  1, 20);
  can_xj_data_reco->cd(1);
   TLegend* xj_data_unfold_2;
   xj_functions::draw_xj_legend(xj_data_unfold_2,"xj_data_unfold_2", 0.17, 0.45, 0.50, 0.59, h_xj_reco_full, h_xj_unfold_half, h_xj_data, h_xj_unfold_data,SimTypeReco,SimTypeRecoUnfolded, "Data", "Data Unfolded", true);
TLegend* xj_data_21;
 xj_functions::DrawsPHENIXLegend(xj_data_21, "xj_data_21", 0.10, 0.60, 0.5, 0.93,"p+p #sqrt{s}=200 GeV  Runs 47352-47733", false,true,true,true);



 //Name and Write histograms

 
	  TString outpath = "UnfoldedxJHists_test_" ;
	  if (Pythia){

	    outpath += "pythia8_";
	  }
	  else if (!Pythia){

	    outpath += "Herwig_";
	  }
	  if (infile1.find("Reweight") != std::string::npos)
	    {
	      outpath += "Reweight_";
	    }

	   if (infile1.find("posJES") != std::string::npos)
	    {
	      outpath += "posJES_";
	    }
	  else if (infile1.find("negJES") != std::string::npos)
	    {
	      outpath += "negJES_";
	    }
	  else if (infile1.find("posJER") != std::string::npos)
	    {
	      outpath += "posJER_";
	    }
	  else if (infile1.find("negJER") != std::string::npos)
	    {
	      outpath += "negJER_";
	    }

	 if (infile1.find("_R2_") != std::string::npos){
	   outpath += "R2";

		}

	 else if (infile1.find("_R4_") != std::string::npos){
	   outpath += "R4";

		}
		
	 
	
	  outpath += ".root";

	  TFile *outfile = TFile::Open(outpath,"RECREATE");
	  h_xj_data->Write();
	  h_xj_unfold_data->Write();
	  h_xj_truth_full->Write();
	  h_xj_reco_full->Write();
	  h_xj_truth_half->Write();
	  h_xj_reco_half->Write();
	  h_xj_unfold_full->Write();
	  h_xj_unfold_half->Write();

	  h_xj_data_20_30->Write();
	  h_xj_unfold_data_20_30->Write();
	  h_xj_truth_full_20_30->Write();
	  h_xj_reco_full_20_30->Write();
	  h_xj_truth_half_20_30->Write();
	  h_xj_reco_half_20_30->Write();
	  h_xj_unfold_full_20_30->Write();
	  h_xj_unfold_half_20_30->Write();

	  h_xj_data_30_40->Write();
	  h_xj_unfold_data_30_40->Write();
	  h_xj_truth_full_30_40->Write();
	  h_xj_reco_full_30_40->Write();
	  h_xj_truth_half_30_40->Write();
	  h_xj_reco_half_30_40->Write();
	  h_xj_unfold_full_30_40->Write();
	  h_xj_unfold_half_30_40->Write();

	  h_xj_data_40_60->Write();
	  h_xj_unfold_data_40_60->Write();
	  h_xj_truth_full_40_60->Write();
	  h_xj_reco_full_40_60->Write();
	  h_xj_truth_half_40_60->Write();
	  h_xj_reco_half_40_60->Write();
	  h_xj_unfold_full_40_60->Write();
	  h_xj_unfold_half_40_60->Write();

	  



}//end function
