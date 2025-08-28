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

void Unfold_xJ_Base_newlib(string infile1 = "Xj_2D_Response_Newlib_Bins19_All-Herwig-Jet10-R04-Run21.root", string infile2 = "Xj_2D_Response_Newlib_Bins19_All-Herwig-Jet30-R04-Run21.root",const std::string configfile = "binning_original.config", int niterations = 6){

  // grab matrices and histograms from files
  TFile* infile01 = TFile::Open(infile1.c_str());
  RooUnfoldResponse *response_pt1pt2_full = (RooUnfoldResponse*)infile01->Get("response_pt1pt2_full");
  RooUnfoldResponse *response_pt1pt2_half_fill = (RooUnfoldResponse*)infile01->Get("response_pt1pt2_half_fill");
  // RooUnfoldResponse *resp_test= (RooUnfoldResponse*)infile01->Get("resp_test");
  TH2D *h_pt1pt2_truth = (TH2D*)infile01->Get("h_pt1pt2_truth");
  TH2D *h_pt1pt2_reco = (TH2D*)infile01->Get("h_pt1pt2_reco");
  TH2D *h_pt1pt2_truth_half = (TH2D*)infile01->Get("h_pt1pt2_truth_half");
  TH2D *h_pt1pt2_reco_half_test = (TH2D*)infile01->Get("h_pt1pt2_reco_half_test");
  
  TH1D* h_pt_lead_truth_matched = (TH1D*)infile01->Get("h_pt_lead_truth_matched");
  TH1D* h_pt_lead_reco_matched = (TH1D*)infile01->Get("h_pt_lead_reco_matched");

  TFile* infile02 = TFile::Open(infile2.c_str());
  RooUnfoldResponse *response_pt1pt2_full_30 = (RooUnfoldResponse*)infile02->Get("response_pt1pt2_full");
  RooUnfoldResponse *response_pt1pt2_half_fill_30 = (RooUnfoldResponse*)infile02->Get("response_pt1pt2_half_fill");
  // RooUnfoldResponse *resp_test= (RooUnfoldResponse*)infile02->Get("resp_test");
  TH2D *h_pt1pt2_truth_30 = (TH2D*)infile02->Get("h_pt1pt2_truth");
  TH2D *h_pt1pt2_reco_30 = (TH2D*)infile02->Get("h_pt1pt2_reco");
  TH2D *h_pt1pt2_truth_half_30 = (TH2D*)infile02->Get("h_pt1pt2_truth_half");
  TH2D *h_pt1pt2_reco_half_test_30 = (TH2D*)infile02->Get("h_pt1pt2_reco_half_test");
  
  TH1D* h_pt_lead_truth_matched_30 = (TH1D*)infile02->Get("h_pt_lead_truth_matched");
  TH1D* h_pt_lead_reco_matched_30 = (TH1D*)infile02->Get("h_pt_lead_reco_matched");

  //combine histos
  response_pt1pt2_full->Add(*response_pt1pt2_full_30);
  response_pt1pt2_half_fill->Add(*response_pt1pt2_half_fill_30);

  h_pt1pt2_truth->Add(h_pt1pt2_truth_30);   
  h_pt1pt2_reco->Add(h_pt1pt2_reco_30);
  h_pt1pt2_truth_half->Add(h_pt1pt2_truth_half_30);
  h_pt1pt2_reco_half_test->Add(h_pt1pt2_reco_half_test_30);
  
  h_pt_lead_truth_matched->Add(h_pt_lead_truth_matched_30);
  h_pt_lead_reco_matched->Add(h_pt_lead_reco_matched_30); 


  //Draw pt

  
 TCanvas* can_ptqa = new TCanvas("can_ptqa", "p_{T1} vs p_{T2}", 800, 800);
  can_ptqa->cd(1);
gPad->SetRightMargin(0.16);
 gPad->SetLogy();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt_lead_truth_matched, "p_{T1, truth} [GeV]", "Counts", 1, 20, 1, 1);
 h_pt_lead_truth_matched->Draw("P");
 h_pt_lead_truth_matched->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_pt_lead_truth_matched->SetAxisRange(pow(10,-9), 1, "Y"); // Adjust Y range if needed
 //h_pt_lead_truth_matched->SetAxisRange(pow(10,-9), 1, "Z"); // Adjust color scale range
 xj_functions::SetHist(h_pt_lead_reco_matched, "p_{T1, truth} [GeV]", "Counts", 2, 20, 1, 1);
 h_pt_lead_reco_matched->Draw("same");
 h_pt_lead_truth_matched->GetZaxis()->SetLabelSize(0.035);
 

  //xj histograms
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

  TH1D *h_xj_reco_full = new TH1D("h_xj_reco_full", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth_full = new TH1D("h_xj_truth_full", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold_full = new TH1D("h_xj_unfold_full", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_reco_half = new TH1D("h_xj_reco_half", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth_half = new TH1D("h_xj_truth_half", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold_half = new TH1D("h_xj_unfold_half", ";x_{J};",nbins, ixj_bins);


  TH1D* h_xj_classical_truth = (TH1D*)infile01->Get("h_xj_classical_truth");
  TH1D* h_xj_classical_truth_30 = (TH1D*)infile02->Get("h_xj_classical_truth");
  h_xj_classical_truth->Add(h_xj_classical_truth_30);
  TH1D* h_xj_classical_reco = (TH1D*)infile01->Get("h_xj_classical_reco");
  TH1D* h_xj_classical_reco_30 = (TH1D*)infile02->Get("h_xj_classical_reco");
  h_xj_classical_reco->Add(h_xj_classical_reco_30);


  //Read in Data File and Set up Histos
  TFile* infile03 = TFile::Open("Hists_xJ_ProjectionTest_TREE_DIJET_v6_1_ana462_2024p010_v001_gl10-00047352-00047733.root");
  TH2D *h_pt1pt2_data = (TH2D*)infile03->Get("h_pt1pt2");
  TH1D *h_xj_data = new TH1D("h_xj_data", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold_data = new TH1D("h_xj_unfold_data", ";x_{J};",nbins, ixj_bins);

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
 

  //Unfolding Process
  TH2D* h_unfolded_pt1pt2_full; TH2D* h_unfolded_pt1pt2_half; TH2D* h_unfolded_pt1pt2_data;
  RooUnfoldBayes *unfolded_pt1pt2_full = new RooUnfoldBayes(response_pt1pt2_full,h_pt1pt2_reco,niterations); 
  RooUnfoldBayes *unfolded_pt1pt2_half = new RooUnfoldBayes(response_pt1pt2_half_fill,h_pt1pt2_reco_half_test,niterations); 
  h_unfolded_pt1pt2_full = (TH2D*)unfolded_pt1pt2_full->Hunfold(RooUnfold::ErrorTreatment::kCovariance);
    h_unfolded_pt1pt2_half = (TH2D*)unfolded_pt1pt2_half->Hunfold(RooUnfold::ErrorTreatment::kCovariance);

  RooUnfoldBayes *unfolded_pt1pt2_data = new RooUnfoldBayes(response_pt1pt2_full,h_pt1pt2_data,niterations);
  h_unfolded_pt1pt2_data = (TH2D*)unfolded_pt1pt2_data->Hunfold(RooUnfold::ErrorTreatment::kCovariance);

  xj_functions::clean_pt1pt2(h_unfolded_pt1pt2_data,nbins,1);
  //project xj
  xj_functions::project_xj(h_pt1pt2_reco, h_xj_reco_full, nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
  xj_functions::project_xj(h_pt1pt2_truth, h_xj_truth_full, nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
  xj_functions::project_xj(h_unfolded_pt1pt2_full, h_xj_unfold_full, nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
   xj_functions::project_xj(h_pt1pt2_reco_half_test, h_xj_reco_half, nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
  xj_functions::project_xj(h_pt1pt2_truth_half, h_xj_truth_half, nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
  xj_functions::project_xj(h_unfolded_pt1pt2_half, h_xj_unfold_half, nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);

  xj_functions::project_xj(h_unfolded_pt1pt2_data, h_xj_unfold_data, nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
   xj_functions::project_xj(h_pt1pt2_data, h_xj_data, nbins, measure_leading_bin, nbins - 7, measure_subleading_bin, nbins - 7);
  
   xj_functions::normalize_histo(h_xj_data, nbins);
   xj_functions::normalize_histo(h_xj_unfold_data,nbins);

  //normalize xj histos
  xj_functions::normalize_histo(h_xj_truth_full, nbins);
  xj_functions::normalize_histo(h_xj_reco_full, nbins);
  xj_functions::normalize_histo(h_xj_truth_half, nbins);
  xj_functions::normalize_histo(h_xj_reco_half, nbins);
  xj_functions::normalize_histo(h_xj_unfold_full, nbins);
  xj_functions::normalize_histo(h_xj_unfold_half, nbins);
  xj_functions::normalize_histo(h_xj_classical_truth, nbins);
  xj_functions::normalize_histo(h_xj_classical_reco, nbins);

 


  //draw histograms
SetsPhenixStyle();
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
  //2D histograms
  TCanvas* can_pt1pt2 = new TCanvas("can_pt1pt2", "p_{T1} vs p_{T2}", 2400, 1600);
  can_pt1pt2->Divide(3,2);
  can_pt1pt2->cd(1);
gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt1pt2_truth, "p_{T1, truth} [GeV]", "p_{T2, truth} [GeV]", 2, 20, 1, 1);
 h_pt1pt2_truth->Draw("Colz");
 h_pt1pt2_truth->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_pt1pt2_truth->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
 h_pt1pt2_truth->SetAxisRange(pow(10,-9), 1, "Z"); // Adjust color scale range
 h_pt1pt2_truth->GetZaxis()->SetLabelSize(0.035);
 

 can_pt1pt2->cd(2);
gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt1pt2_reco, "p_{T1, reco} [GeV]", "p_{T2, reco} [GeV]", 2, 20, 1, 1);
 h_pt1pt2_reco->Draw("Colz");
 h_pt1pt2_reco->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_pt1pt2_reco->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
 h_pt1pt2_reco->SetAxisRange(pow(10,-9), 1, "Z");
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
 h_pt1pt2_reco->GetZaxis()->SetLabelSize(0.035);

 cout << "pointer is: " << h_pt1pt2_reco->GetName() << endl;



 can_pt1pt2->cd(3);
 gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_unfolded_pt1pt2_full, "p_{T1, unfold} [GeV]", "p_{T2, unfold} [GeV]", 2, 20, 1, 1);
 h_unfolded_pt1pt2_full->Draw("Colz");
 h_unfolded_pt1pt2_full->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_unfolded_pt1pt2_full->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
 h_unfolded_pt1pt2_full->SetAxisRange(pow(10,-9), 1, "Z");
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
 h_unfolded_pt1pt2_full->GetZaxis()->SetLabelSize(0.035);

  cout << "made it here 2" << endl;

 can_pt1pt2->cd(4);
gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt1pt2_truth_half, "p_{T1, truth_half} [GeV]", "p_{T2, truth_half} [GeV]", 2, 20, 1, 1);
 h_pt1pt2_truth_half->Draw("Colz");
 h_pt1pt2_truth_half->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_pt1pt2_truth_half->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
  h_pt1pt2_truth_half->SetAxisRange(pow(10,-9), 1, "Z");
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
 h_pt1pt2_truth_half->GetZaxis()->SetLabelSize(0.035);

 can_pt1pt2->cd(5);
gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt1pt2_reco_half_test, "p_{T1, reco_half_test} [GeV]", "p_{T2, reco_half_test} [GeV]", 2, 20, 1, 1);
 h_pt1pt2_reco_half_test->Draw("Colz");
 h_pt1pt2_reco_half_test->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_pt1pt2_reco_half_test->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
  h_pt1pt2_reco_half_test->SetAxisRange(pow(10,-9), 1, "Z");
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
 h_pt1pt2_reco_half_test->GetZaxis()->SetLabelSize(0.035);

 can_pt1pt2->cd(6);
 gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_unfolded_pt1pt2_half, "p_{T1, unfold half} [GeV]", "p_{T2, unfold half} [GeV]", 2, 20, 1, 1);
 h_unfolded_pt1pt2_half->Draw("Colz");
 h_unfolded_pt1pt2_half->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_unfolded_pt1pt2_half->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
 h_unfolded_pt1pt2_half->SetAxisRange(pow(10,-9), 1, "Z");
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
 h_unfolded_pt1pt2_half->GetZaxis()->SetLabelSize(0.035);


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


 TCanvas* can_xj_projected = new TCanvas("can_xj_projected", "p_{T1} vs p_{T2}", 700, 800);
 xj_functions::plotxj_and_ratio(can_xj_projected,h_xj_classical_reco, h_xj_classical_truth,h_xj_truth_full, h_xj_reco_full,true, "x_{J}","Classical/Truth",0, 1.0, 4, 53, 2, 20,  kGreen+1,  20,  8, 53);

 TCanvas* can_xj_full = new TCanvas("can_xj_full", "p_{T1} vs p_{T2}", 700, 800);
 xj_functions::plotxj_and_ratio(can_xj_full,h_xj_reco_full, h_xj_unfold_full,h_xj_truth_full, h_xj_reco_full,false, "x_{J}","Unfold/Truth",0, 1.0, 4, 20, 4, 53,  kGreen+1,  20,  8, 53);
 TLegend* xj_full_unfold;
 xj_functions::draw_xj_legend(xj_full_unfold,"xj_full_unfold", 0.17, 0.45, 0.5, 0.59, h_xj_truth_full,h_xj_reco_full, h_xj_unfold_full, h_xj_unfold_full,"PYTHIA 8","PTYHIA 8 Reco", "Reco Unfolded", " ");
 TLegend* xj_full_2;
 xj_functions::DrawsPHENIXLegend(xj_full_2, "xj_full_2", 0.10, 0.60, 0.5, 0.93,"p+p #sqrt{s}=200 GeV", true,true,false, true);xj_full_2->AddEntry("","Full Closure Test","");

TCanvas* can_xj_half = new TCanvas("can_xj_half", "p_{T1} vs p_{T2}", 700, 800);
 xj_functions::plotxj_and_ratio(can_xj_half,h_xj_reco_half, h_xj_unfold_half,h_xj_truth_half, h_xj_reco_full,false, "x_{J}","Unfold/Truth",0, 1.0, 4, 20, 4, 53,  kGreen+1,  20,  8, 53);
 TLegend* xj_half_unfold;
  xj_functions::draw_xj_legend(xj_half_unfold,"xj_half_unfold", 0.17, 0.45, 0.50, 0.59, h_xj_truth_half, h_xj_reco_half, h_xj_unfold_half, h_xj_unfold_half,"PYTHIA 8","PTYHIA 8 Reco", "Reco Unfolded", " ");
TLegend* xj_half_2;
 xj_functions::DrawsPHENIXLegend(xj_half_2, "xj_half_2", 0.10, 0.60, 0.5, 0.93,"p+p #sqrt{s}=200 GeV", true,true,false,true);xj_half_2->AddEntry("","Half Closure Test","");


 
   //Data Histogram
  TCanvas* can_xj_data = new TCanvas("can_xj_data", "p_{T1} vs p_{T2}", 700, 800);

   // Reduce bottom margin to minimize space between plots
  xj_functions::plotxj_and_ratio(can_xj_data, h_xj_data, h_xj_unfold_data, h_xj_truth_full,h_xj_data,false, "x_{J, Data}","Unfold/Truth",0, 1.0,1 , 20, 1, 53,  kGreen+1,  20,  1, 20);

  can_xj_data->cd(1);
   TLegend* xj_data_unfold;
  xj_functions::draw_xj_legend(xj_data_unfold,"xj_data_unfold", 0.17, 0.45, 0.50, 0.59, h_xj_truth_full, h_xj_data, h_xj_unfold_data, h_xj_unfold_half,"PYTHIA 8","Data", "Data Unfolded", " ");
TLegend* xj_data_2;
 xj_functions::DrawsPHENIXLegend(xj_data_2, "xj_data_2", 0.10, 0.60, 0.5, 0.93,"p+p #sqrt{s}=200 GeV  Runs 47352-47733", true,false,true,true);

   //2D histogram Data

  TCanvas* can_pt1pt2_data = new TCanvas("can_pt1pt2_data", "p_{T1} vs p_{T2}", 2400, 600);
  can_pt1pt2_data->Divide(3,1);
  can_pt1pt2_data->cd(1);
  gPad->SetRightMargin(0.16);
 gPad->SetLogz();
// Draw the 2D histogram
 xj_functions::SetHist(h_pt1pt2_truth, "p_{T1, truth} [GeV]", "p_{T2, truth} [GeV]", 2, 20, 1, 1);
 h_pt1pt2_truth->Draw("Colz");
 h_pt1pt2_truth->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_pt1pt2_truth->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
 h_pt1pt2_truth->SetAxisRange(pow(10,-9), 1, "Z");
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
 h_pt1pt2_truth->GetZaxis()->SetLabelSize(0.035);

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
 xj_functions::SetHist(h_unfolded_pt1pt2_data, "p_{T1, unfold} [GeV]", "p_{T2, unfold} [GeV]", 2, 20, 1, 1);
 h_unfolded_pt1pt2_data->Draw("Colz");
 h_unfolded_pt1pt2_data->SetAxisRange(0, 89, "X"); // Adjust X range if needed
 h_unfolded_pt1pt2_data->SetAxisRange(0, 89, "Y"); // Adjust Y range if needed
 h_unfolded_pt1pt2_data->SetAxisRange(1, pow(10,7), "Z");
//h_pt1pt2->SetAxisRange(1, 375, "Z"); // Adjust color scale range
 h_unfolded_pt1pt2_data->GetZaxis()->SetLabelSize(0.035);

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
   xj_functions::draw_xj_legend(xj_data_unfold_2,"xj_data_unfold_2", 0.17, 0.45, 0.50, 0.59, h_xj_reco_full, h_xj_unfold_half, h_xj_data, h_xj_unfold_data,"PYTHIA 8 Reco","PYTHIA 8 Reco Unfold", "Data", "Data Unfolded", true);
TLegend* xj_data_21;
 xj_functions::DrawsPHENIXLegend(xj_data_21, "xj_data_21", 0.10, 0.60, 0.5, 0.93,"p+p #sqrt{s}=200 GeV  Runs 47352-47733", false,true,true,true);




}//end function
