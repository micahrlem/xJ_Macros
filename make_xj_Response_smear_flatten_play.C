
#include "RooUnfoldResponse.h"
#include "RooUnfold.h"
#include "RooUnfoldBayes.h"

#include "xj_functions.h"
#include "read_binning.h"


#define _USE_MATH_DEFINES

#include <math.h> 
#include <cmath>
#include <iostream>

#include<vector>
#include<array>

void make_xj_Response_smear_flatten_play(const std::string configfile = "binning_original.config", float RValue = 0.4, string infile1 = "pythia8-Jet20-Run21-multiR.root"){

   float Rvaluenew = RValue*10;
  std::ostringstream oss;
  oss << Rvaluenew;
  std::string Rvaluestring = oss.str();
  cout << "R Value is: " << Rvaluestring << endl;
  
  // 
    const double_t pt_bins[]={0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,80};
    const int pt_N = sizeof(pt_bins)/sizeof(pt_bins[0]) - 1;


//define pt and xj bins
 read_binning rb(configfile.c_str());
 Int_t minentries = rb.get_minentries();
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


  

  //random number in each event
  TRandom3 Random;

//Define Histograms

  //Truth Matched QA histograms
TH1D *h_pt_lead_truth_matched = new TH1D("h_pt_lead_truth_matched"," ; Leading Jet p_{T} [GeV]; counts",pt_N,pt_bins);
TH1D *h_eta_lead_truth_matched = new TH1D("h_eta_lead_truth_matched"," ; Leading Jet #eta  ; counts",48,-1.2,1.2);
TH1D *h_phi_lead_truth_matched = new TH1D("h_phi_lead_truth_matched"," ; Leading Jet #Phi  ; counts",nbins,ipt_bins);

TH1D *h_pt_sublead_truth_matched = new TH1D("h_pt_sublead_truth_matched"," ; Subleading Jet p_{T} [GeV]; counts",pt_N,pt_bins);
TH1D *h_eta_sublead_truth_matched = new TH1D("h_eta_sublead_truth_matched"," ; Subleading Jet #eta  ; counts",48,-1.2,1.2);
TH1D *h_phi_sublead_truth_matched = new TH1D("h_phi_sublead_truth_matched"," ; Subleading Jet #Phi  ; counts",32,-TMath::Pi(),TMath::Pi());

//Truth 2D histos and xj base
TH2D *h_pt1pt2_truth =  new TH2D("h_pt1pt2_truth",";p_{T,1};p_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
TH2D *h_pt1pt2_truth_half =  new TH2D("h_pt1pt2_truth_half",";p_{T,1};p_{T,2}", nbins, ipt_bins, nbins, ipt_bins);


TH2D *h_eta_phi_lead_truth_matched = new TH2D("h_eta_phi_lead_truth_matched", "Leading Jet #Phi  ;Leading Jet #eta  ",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);
TH2D *h_eta_phi_sublead_truth_matched = new TH2D("h_eta_phi_sublead_truth_matched", "Subleading Jet #Phi  ;Subleading Jet #eta  ",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);

//Reco Matched QA histograms
TH1D *h_pt_lead_reco_matched = new TH1D("h_pt_lead_reco_matched"," ; Leading Jet p_{T} [GeV]; counts",nbins,ipt_bins);
TH1D *h_eta_lead_reco_matched = new TH1D("h_eta_lead_reco_matched"," ; Leading Jet #eta  ; counts",48,-1.2,1.2);
TH1D *h_phi_lead_reco_matched = new TH1D("h_phi_lead_reco_matched"," ; Leading Jet #Phi  ; counts",nbins,ipt_bins);

TH1D *h_pt_sublead_reco_matched = new TH1D("h_pt_sublead_reco_matched"," ; Subleading Jet p_{T} [GeV]; counts",nbins,ipt_bins);
TH1D *h_eta_sublead_reco_matched = new TH1D("h_eta_sublead_reco_matched"," ; Subleading Jet #eta  ; counts",48,-1.2,1.2);
TH1D *h_phi_sublead_reco_matched = new TH1D("h_phi_sublead_reco_matched"," ; Subleading Jet #Phi  ; counts",32,-TMath::Pi(),TMath::Pi());

//Reco 2D  histos and xj base
TH2D *h_pt1pt2_reco =  new TH2D("h_pt1pt2_reco",";p_{T,1};p_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
TH2D *h_pt1pt2_reco_half_fill =  new TH2D("h_pt1pt2_reco_half_fill",";p_{T,1};p_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
TH2D *h_pt1pt2_reco_half_test =  new TH2D("h_pt1pt2_reco_half_test",";p_{T,1};p_{T,2}", nbins, ipt_bins, nbins, ipt_bins);


TH2D *h_eta_phi_lead_reco_matched = new TH2D("h_eta_phi_lead_reco_matched", "Leading Jet #Phi  ;Leading Jet #eta  ",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);
TH2D *h_eta_phi_sublead_reco_matched = new TH2D("h_eta_phi_sublead_reco_matched", "Subleading Jet #Phi  ;Subleading Jet #eta  ",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);

//Truth vs Reco Histos Matched only
TH2D *h_eta_lead_truth_reco_matched = new TH2D("h_eta_lead_truth_reco_matched", "Leading #eta; Subleading#eta", 48,-1.2,1.2, 48,-1.2,1.2);
TH2D *h_phi_lead_truth_reco_matched = new TH2D("h_phi_lead_truth_reco_matched", "Leading #phi; Subleading#phi", 32,-TMath::Pi(),TMath::Pi(),32,-TMath::Pi(),TMath::Pi());


   //Response Matrix
    RooUnfoldResponse *response_pt1pt2_full = new RooUnfoldResponse("response_pt1pt2_full","");
    RooUnfoldResponse *response_pt1pt2_half_fill = new RooUnfoldResponse("response_pt1pt2_half_fill","");
    RooUnfoldResponse *response_pt1pt2_half_test = new RooUnfoldResponse("response_pt1pt2_half_test","");
    response_pt1pt2_full->Setup(h_pt1pt2_reco, h_pt1pt2_truth);
    response_pt1pt2_half_fill->Setup(h_pt1pt2_reco, h_pt1pt2_truth);
    response_pt1pt2_half_test->Setup(h_pt1pt2_reco, h_pt1pt2_truth);

    //flat and linear Histos
  TH1D *h_linear_truth_xj = new TH1D("h_linear_truth_xj",";x_{J};1/N", 20, 0, 1.0);
  TH1D *h_linear_reco_xj = new TH1D("h_linear_reco_xj",";x_{J};1/N", 20, 0, 1.0);
  TH1D *h_linear_truth_Aj = new TH1D("h_linear_truth_Aj",";A_{J};1/N", 20, 0, 1.0);
  TH1D *h_linear_reco_Aj = new TH1D("h_linear_reco_Aj",";A_{J};1/N", 20, 0, 1.0);

  TH1D *h_xj_classical_truth = new TH1D("h_xj_classical_truth",";x_{J};1/N", nbins, ixj_bins);
  TH1D *h_xj_classical_reco = new TH1D("h_xj_classical_reco",";x_{J};1/N", nbins, ixj_bins);

  TH1D *h_linear_truth_Aj_justmatch = new TH1D("h_linear_truth_Aj_justmatch",";A_{J};1/N", 20, 0, 1.0);
  TH1D *h_linear_truth_xj_justmatch = new TH1D("h_linear_truth_xj_justmatch",";x_{J};1/N", 20, 0, 1.0);


  //30 - 40 cut
  TH1D *h_linear_truth_xj_30_40 = new TH1D("h_linear_truth_xj_30_40",";x_{J};1/N", 20, 0, 1.0);
  TH1D *h_linear_reco_xj_30_40 = new TH1D("h_linear_reco_xj_30_40",";x_{J};1/N", 20, 0, 1.0);
  TH1D *h_linear_truth_Aj_30_40 = new TH1D("h_linear_truth_Aj_30_40",";A_{J};1/N", 20, 0, 1.0);
  TH1D *h_linear_reco_Aj_30_40 = new TH1D("h_linear_reco_Aj_30_40",";A_{J};1/N", 20, 0, 1.0);

  TH1D *h_xj_classical_truth_30_40 = new TH1D("h_xj_classical_truth_30_40",";x_{J};1/N", nbins, ixj_bins);
  TH1D *h_xj_classical_reco_30_40 = new TH1D("h_xj_classical_reco_30_40",";x_{J};1/N", nbins, ixj_bins);


  //40 - 60  cut
   TH1D *h_linear_truth_xj_40_60 = new TH1D("h_linear_truth_xj_40_60",";x_{J};1/N", 20, 0, 1.0);
  TH1D *h_linear_reco_xj_40_60 = new TH1D("h_linear_reco_xj_40_60",";x_{J};1/N", 20, 0, 1.0);
  TH1D *h_linear_truth_Aj_40_60 = new TH1D("h_linear_truth_Aj_40_60",";A_{J};1/N", 20, 0, 1.0);
  TH1D *h_linear_reco_Aj_40_60 = new TH1D("h_linear_reco_Aj_40_60",";A_{J};1/N", 20, 0, 1.0);

  TH1D *h_xj_classical_truth_40_60 = new TH1D("h_xj_classical_truth_40_60",";x_{J};1/N", nbins, ixj_bins);
  TH1D *h_xj_classical_reco_40_60 = new TH1D("h_xj_classical_reco_40_60",";x_{J};1/N", nbins, ixj_bins);
  
  
  // FULL
  TH1D *h_flat_truth_pt1pt2_full = new TH1D("h_truth_flat_pt1pt2_full",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_count_flat_truth_pt1pt2_full = new TH1D("h_truth_count_flat_pt1pt2_full",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_truth_to_response_pt1pt2_full = new TH1D("h_truth_flat_to_response_pt1pt2_full",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  TH1D *h_flat_reco_pt1pt2_full = new TH1D("h_reco_flat_pt1pt2_full",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_count_flat_reco_pt1pt2_full = new TH1D("h_count_reco_flat_pt1pt2_full",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_reco_to_response_pt1pt2_full = new TH1D("h_reco_flat_to_response_pt1pt2_full",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  TH2D *h_flat_response_pt1pt2_full = new TH2D("h_flat_response_pt1pt2_full",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins*nbins, 0, nbins*nbins, nbins*nbins, 0, nbins*nbins);
    int nbin_response = nbins*nbins;
    RooUnfoldResponse rooResponse_dense_pt1pt2_full(nbin_response, 0, nbin_response);

  // HALF_FILL
  TH1D *h_flat_truth_pt1pt2_half_fill = new TH1D("h_truth_flat_pt1pt2_half_fill",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_count_flat_truth_pt1pt2_half_fill = new TH1D("h_truth_count_flat_pt1pt2_half_fill",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_truth_to_response_pt1pt2_half_fill = new TH1D("h_truth_flat_to_response_pt1pt2_half_fill",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  TH1D *h_flat_reco_pt1pt2_half_fill = new TH1D("h_reco_flat_pt1pt2_half_fill",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_count_flat_reco_pt1pt2_half_fill = new TH1D("h_count_reco_flat_pt1pt2_half_fill",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_reco_to_response_pt1pt2_half_fill = new TH1D("h_reco_flat_to_response_pt1pt2_half_fill",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  TH2D *h_flat_response_pt1pt2_half_fill = new TH2D("h_flat_response_pt1pt2_half_fill",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins*nbins, 0, nbins*nbins, nbins*nbins, 0, nbins*nbins);

  RooUnfoldResponse rooResponse_dense_pt1pt2_half_fill(nbin_response, 0, nbin_response);

  // HALF_TEST
  TH1D *h_flat_truth_pt1pt2_half_test = new TH1D("h_truth_flat_pt1pt2_half_test",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_count_flat_truth_pt1pt2_half_test = new TH1D("h_truth_count_flat_pt1pt2_half_test",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_truth_to_response_pt1pt2_half_test = new TH1D("h_truth_flat_to_response_pt1pt2_half_test",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  TH1D *h_flat_reco_pt1pt2_half_test = new TH1D("h_reco_flat_pt1pt2_half_test",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_count_flat_reco_pt1pt2_half_test = new TH1D("h_count_reco_flat_pt1pt2_half_test",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_reco_to_response_pt1pt2_half_test = new TH1D("h_reco_flat_to_response_pt1pt2_half_test",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  TH2D *h_flat_response_pt1pt2_half_test = new TH2D("h_flat_response_pt1pt2_half_test",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins*nbins, 0, nbins*nbins, nbins*nbins, 0, nbins*nbins);




    //connect variables in the code to variables in the tree
   TFile* _file0 = TFile::Open(infile1.c_str());
   TTree* tree = (TTree*)_file0->Get("ttree");

  

   std::vector<double> *tjet_pt ={0}; std::vector<double> *tjet_phi ={0}; std::vector<double> *tjet_eta ={0}; float t_z_vtx = 0;
   std::vector<double> *rjet_pt ={0}; std::vector<double> *rjet_phi ={0}; std::vector<double> *rjet_eta ={0}; float r_z_vtx = 0;
   

        if (Rvaluestring.find("4")  != std::string::npos){
     cout << "Connecting R = 0.4 jets" << endl;
    tree->SetBranchAddress("truth_jet_pt_4",&tjet_pt); tree->SetBranchAddress("truth_jet_eta_4",&tjet_eta); tree->SetBranchAddress("truth_jet_phi_4",&tjet_phi);
    tree->SetBranchAddress("jet_pt_4",&rjet_pt); tree->SetBranchAddress("jet_eta_4",&rjet_eta); tree->SetBranchAddress("jet_phi_4",&rjet_phi);
   }

    if (Rvaluestring.find("2")  != std::string::npos){
       cout << "Connecting R = 0.2 jets" << endl;
    tree->SetBranchAddress("truth_jet_pt_2",&tjet_pt); tree->SetBranchAddress("truth_jet_eta_2",&tjet_eta); tree->SetBranchAddress("truth_jet_phi_2",&tjet_phi);
    tree->SetBranchAddress("jet_pt_2",&rjet_pt); tree->SetBranchAddress("jet_eta_2",&rjet_eta); tree->SetBranchAddress("jet_phi_2",&rjet_phi);
   }
    
tree->SetBranchAddress("truth_vertex_z",&t_z_vtx);tree->SetBranchAddress("mbd_vertex_z",&r_z_vtx);
      //define major variables
      float tleadpt = 0; float tleadphi = 0; float tleadeta = 0; float tzvtx;   
      float rzvtx; float tsubleadpt = 0; float tsubleadphi = 0; float tsubleadeta = 0;
      
      //define cuts
      Double_t dphicut = rb.get_dphicut(); float etacut = 1.1 - RValue; 
      float drminlead = 0.75*RValue; float drminsublead = 0.75*RValue;float zvtxcut = 60;
      float truth_leading_cut = rb.get_truth_leading_cut();
      float truth_subleading_cut = rb.get_truth_subleading_cut();

      float reco_leading_cut = rb.get_reco_leading_cut();
      float reco_subleading_cut = rb.get_reco_subleading_cut();

      float measure_leading_cut = rb.get_measure_leading_cut();
      float measure_subleading_cut = rb.get_measure_subleading_cut();

        std::cout << "Truth1: " << truth_leading_cut << std::endl;
	std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
	std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
	std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
	std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
	std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;
	std::cout << "Dphi Cut: " << dphicut << std::endl;
	std::cout << "Eta Cut: " << etacut << std::endl;
	std::cout << "Zvtx Cut: " << zvtxcut << std::endl;
	std::cout << "R Value is: " << RValue << std::endl;
	std::cout << "drmin is: " << drminlead << std::endl;

	 TF1 *fgaus = new TF1("fgaus", "gaus");
	 fgaus->SetRange(-0.5, 0.5);
	 Double_t JES_sys = rb.get_jes_sys();
	 Double_t JER_sys = rb.get_jer_sys();
	
	 float width = 0.1 + JER_sys;
        
	 fgaus->SetParameters(1, 0, width);
	  std::cout << "JES = " << JES_sys << std::endl;
	 std::cout << "JER = " << JER_sys  << " current width of gaussian: " << width << std::endl;

  if (JER_sys != 0)
    {
      std::cout << "Calculating JER extra = " << JER_sys   << std::endl;
    }

  if (JES_sys != 0)
    {
      std::cout << "Calculating JES extra = " << JES_sys  << std::endl;
    }

	
 // identify file to determine setting
	int setting = 0; float maxpt;float weight;float minpt;float eventscale;

	if (infile1.find("30") != std::string::npos){	  
	  setting = 30; weight = 2.505*pow(10,-9);eventscale = 1;

	   if (Rvaluestring.find("4")  != std::string::npos){
	  minpt = 35;  maxpt = 200;
	   }
	    if (Rvaluestring.find("2")  != std::string::npos){
	  minpt = 31;  maxpt = 200;
	   }
	   
     if(infile1.find("Herwig") != std::string::npos){
        minpt = 30;  maxpt = 200;
       weight = 1.473*pow(10,-12); weight = weight*2.300; cout << "Herwig detected" << endl;
    }
    cout << "30 GEV  JET INFILE" << endl;
   }

  else if (infile1.find("10") != std::string::npos){


    setting = 10;  weight = 3.646*pow(10,-6);eventscale = weight/(2.505*pow(10,-9));

     if (Rvaluestring.find("4")  != std::string::npos){
	  minpt = 14;  maxpt = 22;
	   }
    if (Rvaluestring.find("2")  != std::string::npos){
	  minpt = 12;  maxpt = 20;
	   }
	    
    if(infile1.find("Herwig") != std::string::npos){
      weight = 1.57028*pow(10,-9);  minpt = 10;  maxpt = 30;eventscale = weight/(1.473*pow(10,-12));
      cout << "Herwig detected" << endl;
    }
    cout << "10 GEV JET INFILE" << endl;
   }
     else if (infile1.find("20") != std::string::npos){

       
       setting = 20;weight = 6.218*pow(10,-8);eventscale = weight/(2.505*pow(10,-9));
    cout << "20 GEV MC  JET INFILE" << endl;

      if (Rvaluestring.find("4")  != std::string::npos){
	  minpt = 22;  maxpt = 35;
	   }
      if (Rvaluestring.find("2")  != std::string::npos){
	  minpt = 20;  maxpt = 31;
	   }

   }
  else if (infile1.find("_ana") != std::string::npos){
    setting = 1;maxpt = 200;eventscale = 1;weight=1;
    cout << "DATA JET INFILE" << endl;
   }

	std::cout << "minpt is: " << minpt << " max pt is: " << maxpt << " event scale is: " << eventscale << endl;


	int miss = 0; int match = 0;
      
//matching

   int Nentries = tree->GetEntries();
   //int Nentries = 100000;
    cout << "Nentries = : " << Nentries << endl;
   
    for (int ientry = 0;ientry<Nentries;ientry++){//entry loop start

        //loop over entries
        tree->GetEntry(ientry);
	int ntjets = tjet_pt->size();  int nrjets = rjet_pt->size();
	
	//remove entries without 2 jets]
       

	
	if (ntjets < 2){
	  continue;
	}

	else if (ntjets >= 2){ //can be matched loop
	  tleadpt = 0; tsubleadpt = 0; drminlead = 0.75*RValue; drminsublead = 0.75*RValue;

	  for (int j = 0; j < ntjets; j++){//loop over all truth jets
	      float  tjetpt = tjet_pt->at(j); float tjetphi =  tjet_phi->at(j);  float tjeteta = tjet_eta->at(j); 
	      float rleadpt = 0; float rleadphi = 0; float rleadeta = 0;  float rsubleadpt = 0; float rsubleadphi = 0; float rsubleadeta = 0;
	      float tdeltaphi = 0; float rdeltaphi = 0; float tzvtx = t_z_vtx;float rzvtx = r_z_vtx;

	       //find leading jets
	        if (tjetpt > tleadpt){
			tsubleadpt = tleadpt;	tleadpt = tjetpt;
			tsubleadphi = tleadphi;	tleadphi = tjetphi;
	        	tsubleadeta = tleadeta; tleadeta = tjeteta;
		      }
		  
		 else if (tjetpt > tsubleadpt){
			tsubleadpt = tjetpt;	tsubleadphi = tjetphi;	tsubleadeta = tjeteta;
		      }
	         if (j == (ntjets - 1)){//last entry of event need to match to reco jets
		   
		  

		   for (int k=0; k < nrjets;k++){ //loop over reco jets
		     	       
		      float rjetpt = rjet_pt->at(k);float rjetphi = rjet_phi->at(k); float rjeteta = rjet_eta->at(k); 

		      float detalead = fabs(rjeteta - tleadeta); float dphilead = fabs(rjetphi - tleadphi);
		      float detasublead = fabs(rjeteta - tsubleadeta); float dphisublead = fabs(rjetphi - tsubleadphi); 

		       if (dphilead > TMath::Pi())
			{
			  dphilead = 2*TMath::Pi()-dphilead;
			}

		       if (dphisublead > TMath::Pi())
			{
			  dphisublead = 2*TMath::Pi()-dphisublead;
			}
		      float drlead = TMath::Sqrt(dphilead*dphilead + detalead*detalead);
		      float drsublead = TMath::Sqrt(dphisublead*dphisublead + detasublead*detasublead);

		      if (ientry < 10){

			cout << "Entry: " << ientry <<  " drlead is: " << drlead << " drsublead is: " << drsublead << " drminlead is " << drminlead << endl;


		      }
		     
		      if (drlead < drminlead){//match leading jets
			rleadpt = rjetpt; rleadphi = rjetphi; rleadeta = rjeteta;
			drminlead = drlead;
		      }//end match leading jets
		      if (drsublead < drminsublead){//match subleading jets
			rsubleadpt = rjetpt;rsubleadphi = rjetphi; rsubleadeta = rjeteta;

			drminsublead = drsublead;

		      }//end match subleading jets


		   }// end loop over reco jets

		   // evaluate dphi
		    tdeltaphi = fabs(tleadphi-tsubleadphi); rdeltaphi = fabs(rleadphi - rsubleadphi); 
		    
		    if ( tdeltaphi > TMath::Pi())
			{
			   tdeltaphi = 2*TMath::Pi()-  tdeltaphi;
			}

		       if ( rdeltaphi > TMath::Pi())
			{
			  rdeltaphi = 2*TMath::Pi()- rdeltaphi;
			}
		       // define conditions for fill and for miss



		       //smearing takes place

	  double smear1 = fgaus->GetRandom();
	  double smear2 = fgaus->GetRandom();

	      rleadpt = rleadpt + (JES_sys + smear1)*tleadpt;
	      rsubleadpt = rsubleadpt + (JES_sys + smear2)*tsubleadpt; 
	 
	 
	  //set bins for flat distributions
	    float pt1_truth_bin = nbins;
	  float pt2_truth_bin = nbins;
	  float pt1_reco_bin = nbins;
	  float pt2_reco_bin = nbins;

	  	  for (int ib = 0; ib < nbins; ib++)
	    {
	      if ( tleadpt < ipt_bins[ib+1] && tleadpt >= ipt_bins[ib])
		{
		  pt1_truth_bin = ib;
		}
	      if ( tsubleadpt < ipt_bins[ib+1] && tsubleadpt >= ipt_bins[ib])
		{
		  pt2_truth_bin = ib;
		}
	      if ( rleadpt < ipt_bins[ib+1] && rleadpt >= ipt_bins[ib])
		{
		  pt1_reco_bin = ib;
		}
	      if ( rsubleadpt < ipt_bins[ib+1] && rsubleadpt >= ipt_bins[ib])
		{
		  pt2_reco_bin = ib;
		}
	    }

		       

		       bool truth_good = (tleadpt >= truth_leading_cut && tsubleadpt >= truth_subleading_cut && tdeltaphi > dphicut && fabs(tzvtx) < zvtxcut && fabs(tleadeta) < etacut && fabs(tsubleadeta) < etacut && tleadpt >= minpt && tleadpt <= maxpt);

		       bool reco_good = (rleadpt >= reco_leading_cut && rsubleadpt >= reco_subleading_cut && rdeltaphi > dphicut && fabs(rzvtx) < zvtxcut && fabs(rleadeta) < etacut && fabs(rsubleadeta) < etacut);   
		         double  choice = Random.Rndm();

			 float xjtruth = tsubleadpt/tleadpt;float ajtruth = (tleadpt - tsubleadpt)/(tleadpt + tsubleadpt);

			 if (!truth_good){
			   if (ientry < 10){

			     cout << "Truth Bad in Entry " << ientry << " leadpt is " << tleadpt << " subleadpt is " << tsubleadpt << " deltaphi is " << tdeltaphi << " zvtx is " << rzvtx << " lead eta is " << tleadeta << " subleadeta is " << tsubleadeta << endl;

                           }
			   continue;
			 }
			 //Miss Procedure

			 if (truth_good && !reco_good) {
			     miss++;
			   if (ientry < 100) {
			   
			     cout << "Miss Found: " << ientry
				  << " rleadpt is " << rleadpt
				  << " rsubleadpt is " << rsubleadpt << " rleadeta is: " << rleadeta << " rleadphi is: " << rleadphi << " rsubleadeta is: " << rsubleadeta << " rleadphi is: " << rsubleadphi << " rzvtx is: " << tzvtx   << " choice is: " << choice << endl;
			   }

			   // Fill full set
			   response_pt1pt2_full->Miss(tleadpt, tsubleadpt, eventscale);
			   response_pt1pt2_full->Miss(tsubleadpt, tleadpt, eventscale);
			   h_pt1pt2_truth->Fill(tleadpt, tsubleadpt, eventscale);
			   h_pt1pt2_truth->Fill(tsubleadpt, tleadpt, eventscale);
			   h_xj_classical_truth->Fill(xjtruth, eventscale);
			   if (tleadpt > 20 & tsubleadpt > 5.5){
			   h_linear_truth_xj->Fill(xjtruth,eventscale);  h_linear_truth_Aj->Fill(ajtruth,eventscale);
			   }
			   h_pt_lead_truth_matched->Fill(tleadpt, eventscale);

			   if (tleadpt > 31.18 & tleadpt < 40.7 & tsubleadpt > 9.38){//30_40 hists
			     h_linear_truth_xj_30_40->Fill(xjtruth,eventscale);  h_linear_truth_Aj_30_40->Fill(ajtruth,eventscale);
			     h_xj_classical_truth_30_40->Fill(xjtruth, eventscale);

			   }

			    if (tleadpt > 40.72 & tleadpt < 60.8 & tsubleadpt > 9.38){//40_60 hists
			     h_linear_truth_xj_40_60->Fill(xjtruth,eventscale);  h_linear_truth_Aj_40_60->Fill(ajtruth,eventscale);
			     h_xj_classical_truth_40_60->Fill(xjtruth, eventscale);

			   }
			   
			   // Fill flat + count histograms (FULL)
			   h_flat_truth_to_response_pt1pt2_full->Fill(pt1_truth_bin + nbins * pt2_truth_bin, eventscale);
			   h_flat_truth_to_response_pt1pt2_full->Fill(pt2_truth_bin + nbins * pt1_truth_bin, eventscale);
			   h_flat_truth_pt1pt2_full->Fill(pt1_truth_bin + nbins * pt2_truth_bin, eventscale);
			   h_flat_truth_pt1pt2_full->Fill(pt2_truth_bin + nbins * pt1_truth_bin, eventscale);
			   h_count_flat_truth_pt1pt2_full->Fill(pt1_truth_bin + nbins * pt2_truth_bin);
			   h_count_flat_truth_pt1pt2_full->Fill(pt2_truth_bin + nbins * pt1_truth_bin);
			   rooResponse_dense_pt1pt2_full.Miss(pt1_truth_bin + nbins * pt2_truth_bin, eventscale);
			   rooResponse_dense_pt1pt2_full.Miss(pt2_truth_bin + nbins * pt1_truth_bin, eventscale);

			   // Half closure logic
			   if (choice > 0.5) {
			     response_pt1pt2_half_fill->Miss(tleadpt, tsubleadpt, eventscale);
			     response_pt1pt2_half_fill->Miss(tsubleadpt, tleadpt, eventscale);

			     h_flat_truth_to_response_pt1pt2_half_fill->Fill(pt1_truth_bin + nbins * pt2_truth_bin, eventscale);
			     h_flat_truth_to_response_pt1pt2_half_fill->Fill(pt2_truth_bin + nbins * pt1_truth_bin, eventscale);
			     h_flat_truth_pt1pt2_half_fill->Fill(pt1_truth_bin + nbins * pt2_truth_bin, eventscale);
			     h_flat_truth_pt1pt2_half_fill->Fill(pt2_truth_bin + nbins * pt1_truth_bin, eventscale);
			     h_count_flat_truth_pt1pt2_half_fill->Fill(pt1_truth_bin + nbins * pt2_truth_bin);
			     h_count_flat_truth_pt1pt2_half_fill->Fill(pt2_truth_bin + nbins * pt1_truth_bin);
			     rooResponse_dense_pt1pt2_half_fill.Miss(pt1_truth_bin + nbins * pt2_truth_bin, eventscale);
			     rooResponse_dense_pt1pt2_half_fill.Miss(pt2_truth_bin + nbins * pt1_truth_bin, eventscale);
			   }

			   if (choice < 0.5) {
			     h_pt1pt2_truth_half->Fill(tleadpt, tsubleadpt, eventscale);
			     h_pt1pt2_truth_half->Fill(tsubleadpt, tleadpt, eventscale);

			     h_flat_truth_to_response_pt1pt2_half_test->Fill(pt1_truth_bin + nbins * pt2_truth_bin, eventscale);
			     h_flat_truth_to_response_pt1pt2_half_test->Fill(pt2_truth_bin + nbins * pt1_truth_bin, eventscale);
			     h_flat_truth_pt1pt2_half_test->Fill(pt1_truth_bin + nbins * pt2_truth_bin, eventscale);
			     h_flat_truth_pt1pt2_half_test->Fill(pt2_truth_bin + nbins * pt1_truth_bin, eventscale);
			     h_count_flat_truth_pt1pt2_half_test->Fill(pt1_truth_bin + nbins * pt2_truth_bin);
			     h_count_flat_truth_pt1pt2_half_test->Fill(pt2_truth_bin + nbins * pt1_truth_bin);
			     // No RooUnfoldResponse fill here for half_test
			   }
			 }
//end Miss procedure

			 //Fill Rest of Histograms
			if (truth_good && reco_good) {

			  match++;
			  float xjtruth = tsubleadpt / tleadpt;float ajtruth = (tleadpt - tsubleadpt)/(tleadpt + tsubleadpt);
			  float xjreco  = rsubleadpt / rleadpt;float ajreco = (rleadpt - rsubleadpt)/(rleadpt + rsubleadpt);

			  // Truth info
			  h_pt_lead_truth_matched->Fill(tleadpt, eventscale);
			  h_eta_lead_truth_matched->Fill(tleadeta, eventscale); 
			  h_phi_lead_truth_matched->Fill(tleadphi, eventscale);
			  h_pt_sublead_truth_matched->Fill(tsubleadpt, eventscale); 
			  h_eta_sublead_truth_matched->Fill(tsubleadeta, eventscale);
			  h_phi_sublead_truth_matched->Fill(tsubleadphi, eventscale);
			  h_pt1pt2_truth->Fill(tleadpt, tsubleadpt, eventscale); 
			  h_pt1pt2_truth->Fill(tsubleadpt, tleadpt, eventscale);
			  h_xj_classical_truth->Fill(xjtruth, eventscale);
			    if (tleadpt > 20 & tsubleadpt > 5.5){
			   h_linear_truth_xj->Fill(xjtruth,eventscale);  h_linear_truth_Aj->Fill(ajtruth,eventscale);
			   }
			   h_linear_truth_xj_justmatch->Fill(xjtruth,eventscale); h_linear_truth_Aj_justmatch->Fill(ajtruth,eventscale);

			    if (tleadpt > 31.18 & tleadpt < 40.7 & tsubleadpt > 9.38){//30_40 hists
			     h_linear_truth_xj_30_40->Fill(xjtruth,eventscale);  h_linear_truth_Aj_30_40->Fill(ajtruth,eventscale);
			     h_xj_classical_truth_30_40->Fill(xjtruth, eventscale);

			   }

			   if (tleadpt > 40.72 & tleadpt < 60.8 & tsubleadpt > 9.38){//40_60 hists
			     h_linear_truth_xj_40_60->Fill(xjtruth,eventscale);  h_linear_truth_Aj_40_60->Fill(ajtruth,eventscale);
			     h_xj_classical_truth_40_60->Fill(xjtruth, eventscale);

			   }

			    if (rleadpt > 31.18 & rleadpt < 40.7 & rsubleadpt > 9.38){//30_40 hists
			     h_linear_reco_xj_30_40->Fill(xjreco,eventscale);  h_linear_reco_Aj_30_40->Fill(ajreco,eventscale);
			     h_xj_classical_reco_30_40->Fill(xjreco, eventscale);

			   }

			   if (rleadpt > 40.72 & rleadpt < 60.8 & rsubleadpt > 9.38){//40_60 hists
			     h_linear_reco_xj_40_60->Fill(xjreco,eventscale);  h_linear_reco_Aj_40_60->Fill(ajreco,eventscale);
			     h_xj_classical_reco_40_60->Fill(xjreco, eventscale);

			   }

			  // Reco info
			  h_pt_lead_reco_matched->Fill(rleadpt, eventscale);
			  h_eta_lead_reco_matched->Fill(rleadeta, eventscale); 
			  h_phi_lead_reco_matched->Fill(rleadphi, eventscale);
			  h_pt_sublead_reco_matched->Fill(rsubleadpt, eventscale); 
			  h_eta_sublead_reco_matched->Fill(rsubleadeta, eventscale);
			  h_phi_sublead_reco_matched->Fill(rsubleadphi, eventscale); 
			  h_pt1pt2_reco->Fill(rleadpt, rsubleadpt, eventscale); 
			  h_pt1pt2_reco->Fill(rsubleadpt, rleadpt, eventscale);  
			  h_xj_classical_reco->Fill(xjreco, eventscale);
			   h_linear_reco_xj->Fill(xjreco,eventscale); h_linear_reco_Aj->Fill(ajreco,eventscale);


			  // Full response
			  response_pt1pt2_full->Fill(rleadpt, rsubleadpt, tleadpt, tsubleadpt, eventscale);
			  response_pt1pt2_full->Fill(rsubleadpt, rleadpt, tsubleadpt, tleadpt, eventscale);

			  h_flat_response_pt1pt2_full->Fill(pt1_reco_bin + nbins * pt2_reco_bin, pt1_truth_bin + nbins * pt2_truth_bin, eventscale);
			  h_flat_response_pt1pt2_full->Fill(pt2_reco_bin + nbins * pt1_reco_bin, pt2_truth_bin + nbins * pt1_truth_bin, eventscale);
			  h_flat_reco_to_response_pt1pt2_full->Fill(pt1_reco_bin + nbins * pt2_reco_bin, eventscale);
			  h_flat_reco_to_response_pt1pt2_full->Fill(pt2_reco_bin + nbins * pt1_reco_bin, eventscale);
			  h_flat_reco_pt1pt2_full->Fill(pt1_reco_bin + nbins * pt2_reco_bin, eventscale);
			  h_flat_reco_pt1pt2_full->Fill(pt2_reco_bin + nbins * pt1_reco_bin, eventscale);
			  h_count_flat_reco_pt1pt2_full->Fill(pt1_reco_bin + nbins * pt2_reco_bin);
			  h_count_flat_reco_pt1pt2_full->Fill(pt2_reco_bin + nbins * pt1_reco_bin);
			  rooResponse_dense_pt1pt2_full.Fill(pt1_reco_bin + nbins * pt2_reco_bin, pt1_truth_bin + nbins * pt2_truth_bin, eventscale);
			  rooResponse_dense_pt1pt2_full.Fill(pt2_reco_bin + nbins * pt1_reco_bin, pt2_truth_bin + nbins * pt1_truth_bin, eventscale);

			  if (choice > 0.5) {
			    h_pt1pt2_reco_half_fill->Fill(rleadpt, rsubleadpt, eventscale);
			    h_pt1pt2_reco_half_fill->Fill(rsubleadpt, rleadpt, eventscale);

			    response_pt1pt2_half_fill->Fill(rleadpt, rsubleadpt, tleadpt, tsubleadpt, eventscale);
			    response_pt1pt2_half_fill->Fill(rsubleadpt, rleadpt, tsubleadpt, tleadpt, eventscale);

			    h_flat_response_pt1pt2_half_fill->Fill(pt1_reco_bin + nbins * pt2_reco_bin, pt1_truth_bin + nbins * pt2_truth_bin, eventscale);
			    h_flat_response_pt1pt2_half_fill->Fill(pt2_reco_bin + nbins * pt1_reco_bin, pt2_truth_bin + nbins * pt1_truth_bin, eventscale);
			    h_flat_reco_to_response_pt1pt2_half_fill->Fill(pt1_reco_bin + nbins * pt2_reco_bin, eventscale);
			    h_flat_reco_to_response_pt1pt2_half_fill->Fill(pt2_reco_bin + nbins * pt1_reco_bin, eventscale);
			    h_flat_reco_pt1pt2_half_fill->Fill(pt1_reco_bin + nbins * pt2_reco_bin, eventscale);
			    h_flat_reco_pt1pt2_half_fill->Fill(pt2_reco_bin + nbins * pt1_reco_bin, eventscale);
			    h_count_flat_reco_pt1pt2_half_fill->Fill(pt1_reco_bin + nbins * pt2_reco_bin);
			    h_count_flat_reco_pt1pt2_half_fill->Fill(pt2_reco_bin + nbins * pt1_reco_bin);
			    rooResponse_dense_pt1pt2_half_fill.Fill(pt1_reco_bin + nbins * pt2_reco_bin, pt1_truth_bin + nbins * pt2_truth_bin, eventscale);
			    rooResponse_dense_pt1pt2_half_fill.Fill(pt2_reco_bin + nbins * pt1_reco_bin, pt2_truth_bin + nbins * pt1_truth_bin, eventscale);
			  } 
			  else if (choice < 0.5) {
			    h_pt1pt2_reco_half_test->Fill(rleadpt, rsubleadpt, eventscale);
			    h_pt1pt2_reco_half_test->Fill(rsubleadpt, rleadpt, eventscale);

			    response_pt1pt2_half_test->Fill(rleadpt, rsubleadpt, tleadpt, tsubleadpt, eventscale);
			    response_pt1pt2_half_test->Fill(rsubleadpt, rleadpt, tsubleadpt, tleadpt, eventscale);

			    h_flat_response_pt1pt2_half_test->Fill(pt1_reco_bin + nbins * pt2_reco_bin, pt1_truth_bin + nbins * pt2_truth_bin, eventscale);
			    h_flat_response_pt1pt2_half_test->Fill(pt2_reco_bin + nbins * pt1_reco_bin, pt2_truth_bin + nbins * pt1_truth_bin, eventscale);
			    h_flat_reco_to_response_pt1pt2_half_test->Fill(pt1_reco_bin + nbins * pt2_reco_bin, eventscale);
			    h_flat_reco_to_response_pt1pt2_half_test->Fill(pt2_reco_bin + nbins * pt1_reco_bin, eventscale);
			    h_flat_reco_pt1pt2_half_test->Fill(pt1_reco_bin + nbins * pt2_reco_bin, eventscale);
			    h_flat_reco_pt1pt2_half_test->Fill(pt2_reco_bin + nbins * pt1_reco_bin, eventscale);
			    h_count_flat_reco_pt1pt2_half_test->Fill(pt1_reco_bin + nbins * pt2_reco_bin);
			    h_count_flat_reco_pt1pt2_half_test->Fill(pt2_reco_bin + nbins * pt1_reco_bin);

			    h_pt1pt2_truth_half->Fill(tleadpt, tsubleadpt, eventscale);
			    h_pt1pt2_truth_half->Fill(tsubleadpt, tleadpt, eventscale);
			  }
			}	//End rest of histograms

		

		 }// end last entry loop

		 

	  }// end loop over all truth jets
	}// end can be matched loop
	
	

	
    }// end loop of entries

    cout << "misses: " << miss <<  " matches: " << match <<  " total: " << miss + match << endl;

    //define file for writing
     TString outfilename = infile1;
     std::string nbinsstring = "_Bins";
     TString responsepath = "_Play_Smear_Flatten_";

     if (dphicut > 2.5){
       responsepath += "highdphi_";
     }
     
      if (JES_sys > 0)
	    {
	      responsepath += "posJES_";
	    }
	  else if (JES_sys < 0)
	    {
	      responsepath += "negJES_";
	    }
	  else if (JER_sys > 0)
	    {
	      responsepath += "posJER_";
	    }
	  else if (JER_sys < 0)
	    {
	      responsepath += "negJER_";
	    }
      nbinsstring += std::to_string(nbins);
      // Convert std::string to TString explicitly
      TString nbins_tstr = TString(nbinsstring);         // from std::string nbinsstring
      TString Rvalue_tstr = TString(Rvaluestring);       // from std::string Rvaluestring

      // Proper concatenation — now all components are TString
      outfilename.Prepend("Xj_2D_Response" + nbins_tstr + "_R" + Rvalue_tstr + responsepath);

      TFile *outfile = TFile::Open(outfilename,"RECREATE");


//write histograms
                       ///Truth pt, eta, phi
			   h_pt_lead_truth_matched->Write();
                            h_eta_lead_truth_matched->Write(); 
                            h_phi_lead_truth_matched->Write();

                           h_pt_sublead_truth_matched->Write(); 
                           h_eta_sublead_truth_matched->Write();
                           h_phi_sublead_truth_matched->Write();

                        //Truth 2D histos and xj base
                            h_pt1pt2_truth->Write();                         
                            h_xj_classical_truth->Write();

			       ///Reco pt, eta, phi
			    h_pt_lead_reco_matched->Write();
                            h_eta_lead_reco_matched->Write(); 
                            h_phi_lead_reco_matched->Write();

                           h_pt_sublead_reco_matched->Write(); 
                           h_eta_sublead_reco_matched->Write();
                           h_phi_sublead_reco_matched->Write(); 

			 //Reco 2D histos and xj base
                            h_pt1pt2_reco->Write(); 
                            h_xj_classical_reco->Write(); 

			      //Full Response Matrix 
			       response_pt1pt2_full->Write();
			      
			     
			      //Half Matrices Fill
			      h_pt1pt2_reco_half_fill->Write();
			      h_pt1pt2_truth_half->Write(); 
			       response_pt1pt2_half_fill->Write();

			      //Half Matrices Test
			       h_pt1pt2_reco_half_test->Write();
			       response_pt1pt2_half_test->Write();

			      // ==== WRITE FLAT, AND RESPONSE HISTOGRAMS ====

        
			       // FULL
			       h_flat_truth_pt1pt2_full->Write();
			       h_count_flat_truth_pt1pt2_full->Write();
			       h_flat_reco_pt1pt2_full->Write();
			       h_count_flat_reco_pt1pt2_full->Write();
			       h_flat_response_pt1pt2_full->Write();
			       rooResponse_dense_pt1pt2_full.Write("rooResponse_dense_pt1pt2_full");

			       // HALF_FILL
			       h_flat_truth_pt1pt2_half_fill->Write();
			       h_count_flat_truth_pt1pt2_half_fill->Write();
			       h_flat_reco_pt1pt2_half_fill->Write();
			       h_count_flat_reco_pt1pt2_half_fill->Write();
			       h_flat_response_pt1pt2_half_fill->Write();
			       rooResponse_dense_pt1pt2_half_fill.Write("rooResponse_dense_pt1pt2_half_fill");

			       // HALF_TEST
			       h_flat_truth_pt1pt2_half_test->Write();
			       h_count_flat_truth_pt1pt2_half_test->Write();
			       h_flat_reco_pt1pt2_half_test->Write();
			       h_count_flat_reco_pt1pt2_half_test->Write();
			       h_flat_response_pt1pt2_half_test->Write();
			       // No RooUnfoldResponse for HALF_TEST

			       h_linear_truth_xj->Write();  h_linear_reco_xj->Write();
			       h_linear_truth_xj_justmatch->Write();
			        h_linear_truth_Aj->Write();  h_linear_reco_Aj->Write();
			       h_linear_truth_Aj_justmatch->Write();


			       // 30–40% linear
			       h_linear_truth_xj_30_40->Write();
			       h_linear_reco_xj_30_40->Write();
			       h_linear_truth_Aj_30_40->Write();
			       h_linear_reco_Aj_30_40->Write();

			       // 30–40% classical
			       h_xj_classical_truth_30_40->Write();
			       h_xj_classical_reco_30_40->Write();

			       // 40–60% linear
			       h_linear_truth_xj_40_60->Write();
			       h_linear_reco_xj_40_60->Write();
			       h_linear_truth_Aj_40_60->Write();
			       h_linear_reco_Aj_40_60->Write();

			       // 40–60% classical
			       h_xj_classical_truth_40_60->Write();
			        h_xj_classical_reco_40_60->Write();
			      



}//end function
