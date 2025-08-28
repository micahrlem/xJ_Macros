
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

void make_Data_xJ_Flat(const std::string configfile = "binning_original.config", float RValue = 0.4, string infile1 = "outputestdata20segmentcalofitting.root"){

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


//Reco Matched QA histograms
TH1D *h_pt_lead_data = new TH1D("h_pt_lead_data"," ; Leading Jet p_{T} [GeV]; counts",nbins,ipt_bins);
TH1D *h_eta_lead_data = new TH1D("h_eta_lead_data"," ; Leading Jet #eta  ; counts",48,-1.2,1.2);
TH1D *h_phi_lead_data = new TH1D("h_phi_lead_data"," ; Leading Jet #Phi  ; counts",nbins,ipt_bins);

TH1D *h_pt_sublead_data = new TH1D("h_pt_sublead_data"," ; Subleading Jet p_{T} [GeV]; counts",nbins,ipt_bins);
TH1D *h_eta_sublead_data = new TH1D("h_eta_sublead_data"," ; Subleading Jet #eta  ; counts",48,-1.2,1.2);
TH1D *h_phi_sublead_data = new TH1D("h_phi_sublead_data"," ; Subleading Jet #Phi  ; counts",32,-TMath::Pi(),TMath::Pi());

//Reco 2D  histos and xj base
TH2D *h_pt1pt2_data =  new TH2D("h_pt1pt2_data",";p_{T,1};p_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
TH1D *h_xj_classical_data = new TH1D("h_xj_classical_data",";x_{J};1/N", nbins, ixj_bins);

TH2D *h_eta_phi_lead_data = new TH2D("h_eta_phi_lead_data", "Leading Jet #Phi  ;Leading Jet #eta  ",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);
TH2D *h_eta_phi_sublead_data = new TH2D("h_eta_phi_sublead_data", "Subleading Jet #Phi  ;Subleading Jet #eta  ",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);

    //flat and linear Histos
  TH1D *h_linear_truth_xj = new TH1D("h_lineartruth_xj",";A_{J};1/N", 20, 0, 1.0);
  TH1D *h_linear_data_xj = new TH1D("h_linearreco_xj",";A_{J};1/N", 20, 0, 1.0);

  // FULL
  TH1D *h_flat_data_pt1pt2 = new TH1D("h_data_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_count_flat_data_pt1pt2 = new TH1D("h_count_data_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_data_to_response_pt1pt2 = new TH1D("h_data_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);


    int nbin_response = nbins*nbins;
  

    //connect variables in the code to variables in the tree
   TFile* _file0 = TFile::Open(infile1.c_str());
   TTree* tree = (TTree*)_file0->Get("ttree");

  

  
   std::vector<double> *djet_pt ={0}; std::vector<double> *djet_phi ={0}; std::vector<double> *djet_eta ={0}; float d_z_vtx = 0;
     unsigned long long triggers;tree->SetBranchAddress("gl1_scaled",&triggers);

        if (Rvaluestring.find("4")  != std::string::npos){
     cout << "Connecting R = 0.4 jets" << endl;
   
    tree->SetBranchAddress("jet_pt_4",&djet_pt); tree->SetBranchAddress("jet_eta_4",&djet_eta); tree->SetBranchAddress("jet_phi_4",&djet_phi);
   }

    if (Rvaluestring.find("2")  != std::string::npos){
       cout << "Connecting R = 0.2 jets" << endl;
   
    tree->SetBranchAddress("jet_pt_2",&djet_pt); tree->SetBranchAddress("jet_eta_2",&djet_eta); tree->SetBranchAddress("jet_phi_2",&djet_phi);
   }
    
   tree->SetBranchAddress("mbd_vertex_z",&d_z_vtx);
      //define major variables
      float dleadpt = 0; float dleadphi = 0; float dleadeta = 0; float dzvtx;   
     float dsubleadpt = 0; float dsubleadphi = 0; float dsubleadeta = 0;
      
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
	 std::cout << "JES = " << JES_sys << std::endl;
	 std::cout << "JER = " << JER_sys << std::endl;
	 float width = 0.1 + JER_sys;
	 fgaus->SetParameters(1, 0, width);

  if (JER_sys != 0)
    {
      std::cout << "Calculating JER extra = " << JER_sys  << std::endl;
    }

  if (JES_sys != 0)
    {
      std::cout << "Calculating JES extra = " << JES_sys  << std::endl;
    }

	
 // identify file to determine setting
	int setting = 0; float maxpt;float weight;float minpt;float eventscale;
	if (infile1.find("data") != std::string::npos){
	  setting = 1;maxpt = 200;eventscale = 1; minpt = 5.5;
    cout << "DATA JET INFILE" << endl;
   }

	else if (infile1.find("30") != std::string::npos){
	  setting = 30; minpt = 30;  maxpt = 200;weight = 2.505*pow(10,-9);eventscale = 1;
     if(infile1.find("Herwig") != std::string::npos){
       weight = 1.473*pow(10,-12); cout << "Herwig detected" << endl;
    }
    cout << "30 GEV  JET INFILE" << endl;
   }

  else if (infile1.find("10") != std::string::npos){
    setting = 10; minpt = 10;  maxpt = 20; weight = 3.646*pow(10,-6);eventscale = weight/(2.505*pow(10,-9));
    if(infile1.find("Herwig") != std::string::npos){
      weight = 1.57028*pow(10,-10); maxpt = 30;eventscale = weight/(1.473*pow(10,-12));
      cout << "Herwig detected" << endl;
    }
    cout << "10 GEV JET INFILE" << endl;
   }
     else if (infile1.find("20") != std::string::npos){
       setting = 20;weight = 6.218*pow(10,-8);eventscale = weight/(2.505*pow(10,-9));
    cout << "20 GEV MC  JET INFILE" << endl;
    maxpt = 30; minpt = 20;

   }
 
	std::cout << "minpt is: " << minpt << " max pt is: " << maxpt << " event scale is: " << eventscale << endl;


  
      
//matching

   int Nentries = tree->GetEntries();
   //int Nentries = 100000;
    cout << "Nentries = : " << Nentries << endl;
   
   for (int ientry = 0; ientry < Nentries; ientry++) { // entry loop start

     if ((ientry + 1) % 100000 == 0 || ientry == Nentries - 1) {
        std::cout << "Processed " << (ientry + 1) << " / " << Nentries << " entries..." << std::endl;
    }

     tree->GetEntry(ientry);
     bool flag = false;

     // Check if bit 17 is set
     if (triggers & (1LL << 17)) {
       flag = true; //cout << "trigger found" << endl;
     }

     int ndjets = djet_pt->size();
     if (ndjets < 2) continue;

     if (ndjets >= 2) { // can be matched loop
       dleadpt = 0; dsubleadpt = 0; drminlead = 0.75 * RValue; drminsublead = 0.75 * RValue;

       for (int j = 0; j < ndjets; j++) { // loop over all detector jets
	 float djetpt  = djet_pt->at(j);
	 float djetphi = djet_phi->at(j);
	 float djeteta = djet_eta->at(j);
	 float ddeltaphi = 0;
	 float dzvtx = d_z_vtx;

	 // find leading and subleading detector jets
	 if (djetpt > dleadpt) {
	   dsubleadpt = dleadpt;     dleadpt = djetpt;
	   dsubleadphi = dleadphi;   dleadphi = djetphi;
	   dsubleadeta = dleadeta;   dleadeta = djeteta;
	 } else if (djetpt > dsubleadpt) {
	   dsubleadpt = djetpt; dsubleadphi = djetphi; dsubleadeta = djeteta;
	 }

	 if (j == ndjets - 1) { // last jet in event
	   ddeltaphi = fabs(dleadphi - dsubleadphi);
	   if (ddeltaphi > TMath::Pi()) ddeltaphi = 2 * TMath::Pi() - ddeltaphi;

	   float pt1_data_bin = nbins;
	   float pt2_data_bin = nbins;

	   for (int ib = 0; ib < nbins; ib++) {
	     if (dleadpt < ipt_bins[ib + 1] && dleadpt >= ipt_bins[ib]) pt1_data_bin = ib;
	     if (dsubleadpt < ipt_bins[ib + 1] && dsubleadpt >= ipt_bins[ib]) pt2_data_bin = ib;
	   }

	   bool data_good = ( dleadpt >= reco_leading_cut && dsubleadpt >= reco_subleading_cut && ddeltaphi > dphicut && fabs(dzvtx) < zvtxcut &&
			     fabs(dleadeta) < etacut && fabs(dsubleadeta) < etacut && flag  );

	   double choice = Random.Rndm();

	   if (data_good) {
	     float xjreco = dsubleadpt / dleadpt;

	     h_pt_lead_data->Fill(dleadpt, eventscale);
	     h_eta_lead_data->Fill(dleadeta, eventscale);
	     h_phi_lead_data->Fill(dleadphi, eventscale);
	     h_pt_sublead_data->Fill(dsubleadpt, eventscale);
	     h_eta_sublead_data->Fill(dsubleadeta, eventscale);
	     h_phi_sublead_data->Fill(dsubleadphi, eventscale);
	     h_pt1pt2_data->Fill(dleadpt, dsubleadpt, eventscale);
	     h_pt1pt2_data->Fill(dsubleadpt, dleadpt, eventscale);
	     h_xj_classical_data->Fill(xjreco, eventscale);

	   
	     h_flat_data_pt1pt2->Fill(pt1_data_bin + nbins * pt2_data_bin, eventscale);
	     h_flat_data_pt1pt2->Fill(pt2_data_bin + nbins * pt1_data_bin, eventscale);
	     h_count_flat_data_pt1pt2->Fill(pt1_data_bin + nbins * pt2_data_bin);
	     h_count_flat_data_pt1pt2->Fill(pt2_data_bin + nbins * pt1_data_bin);
	   }
	 } // end last jet
       }   // end jet loop
     }     // end ndjets >= 2
   }       // end event loop


    //define file for writing
     TString outfilename = infile1;
     std::string nbinsstring = "_Bins";
     TString responsepath = "_Data_Flatten_";
      if (JES_sys > 0)
	    {
	      responsepath += "_posJES_";
	    }
	  else if (JES_sys < 0)
	    {
	      responsepath += "_negJES_";
	    }
	  else if (JER_sys > 0)
	    {
	      responsepath += "_posJER_";
	    }
	  else if (JER_sys < 0)
	    {
	      responsepath += "_negJER_";
	    }
      nbinsstring += std::to_string(nbins);
      // Convert std::string to TString explicitly
      TString nbins_tstr = TString(nbinsstring);         // from std::string nbinsstring
      TString Rvalue_tstr = TString(Rvaluestring);       // from std::string Rvaluestring

      // Proper concatenation — now all components are TString
      outfilename.Prepend("Xj_2D_Response" + nbins_tstr + "_R" + Rvalue_tstr + responsepath);

      TFile *outfile = TFile::Open(outfilename,"RECREATE");


//write histograms
                     

			       ///Reco pt, eta, phi
			    h_pt_lead_data->Write();
                            h_eta_lead_data->Write(); 
                            h_phi_lead_data->Write();

                           h_pt_sublead_data->Write(); 
                           h_eta_sublead_data->Write();
                           h_phi_sublead_data->Write(); 

			 //Reco 2D histos and xj base
                            h_pt1pt2_data->Write(); 
                            h_xj_classical_data->Write(); 

	  

			      // ==== WRITE FLAT, AND RESPONSE HISTOGRAMS ====

			       h_flat_data_pt1pt2->Write();
			       h_count_flat_data_pt1pt2->Write();
			       
			      

			   

			    


}//end function
