void Make_Iso_Split_Study_UnmatchTruth(string infile1 = "pythia8-Jet10-Run21-multiR.root", const std::string configfile = "binning_original.config",float Rvalue = 0.4){
  #include <sstream>

  float Rvaluenew = Rvalue*10;
  std::ostringstream oss;
  oss << Rvaluenew;
  std::string Rvaluestring = oss.str();
  cout << "R Value is: " << Rvaluestring << endl;
 
  //#include "xj_functions.h"
#include "read_binning.h"

#define _USE_MATH_DEFINES

#include <math.h> 
#include <cmath>
#include <iostream>

#include<vector>
#include<array>

    const double_t pt_bins[]={0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,80};
    const int pt_N = sizeof(pt_bins)/sizeof(pt_bins[0]) - 1;

    //define pt and xj bins
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


    // Define Truth MC Histos
  // === TH1D Histograms: pt ===
  TH1D* hTruthPtLead_all       = new TH1D("hTruthPtLead_all", "", pt_N, pt_bins);
  TH1D* hTruthPtLead_iso       = new TH1D("hTruthPtLead_iso", "", pt_N, pt_bins);
  TH1D* hTruthPtSubLead_all    = new TH1D("hTruthPtSubLead_all", "", pt_N, pt_bins);
  TH1D* hTruthPtSubLead_iso    = new TH1D("hTruthPtSubLead_iso", "", pt_N, pt_bins);
  TH1D* hTruthPtInclusive_all  = new TH1D("hTruthPtInclusive_all", "", pt_N, pt_bins);
  TH1D* hTruthPtInclusive_iso  = new TH1D("hTruthPtInclusive_iso", "", pt_N, pt_bins);

  // === TH1D Histograms: eta ===
  TH1D* hTruthEtaLead_all      = new TH1D("hTruthEtaLead_all", "", 48, -1.2, 1.2);
  TH1D* hTruthEtaLead_iso      = new TH1D("hTruthEtaLead_iso", "", 48, -1.2, 1.2);
  TH1D* hTruthEtaSubLead_all   = new TH1D("hTruthEtaSubLead_all", "", 48, -1.2, 1.2);
  TH1D* hTruthEtaSubLead_iso   = new TH1D("hTruthEtaSubLead_iso", "", 48, -1.2, 1.2);
  TH1D* hTruthEtaInclusive_all = new TH1D("hTruthEtaInclusive_all", "", 48, -1.2, 1.2);
  TH1D* hTruthEtaInclusive_iso = new TH1D("hTruthEtaInclusive_iso", "", 48, -1.2, 1.2);

  // === TH1D Histograms: phi ===
  TH1D* hTruthPhiLead_all      = new TH1D("hTruthPhiLead_all", "", 32, -TMath::Pi(), TMath::Pi());
  TH1D* hTruthPhiLead_iso      = new TH1D("hTruthPhiLead_iso", "", 32, -TMath::Pi(), TMath::Pi());
  TH1D* hTruthPhiSubLead_all   = new TH1D("hTruthPhiSubLead_all", "", 32, -TMath::Pi(), TMath::Pi());
  TH1D* hTruthPhiSubLead_iso   = new TH1D("hTruthPhiSubLead_iso", "", 32, -TMath::Pi(), TMath::Pi());
  TH1D* hTruthPhiInclusive_all = new TH1D("hTruthPhiInclusive_all", "", 32, -TMath::Pi(), TMath::Pi());
  TH1D* hTruthPhiInclusive_iso = new TH1D("hTruthPhiInclusive_iso", "", 32, -TMath::Pi(), TMath::Pi());

  // === TH1D Histograms: dphi/deta ===
  TH1D* hTruthDeltaPhi_all     = new TH1D("hTruthDeltaPhi_all", "", 50, TMath::Pi() / 2, 2 * TMath::Pi());
  TH1D* hTruthDeltaPhi_iso     = new TH1D("hTruthDeltaPhi_iso", "", 50, TMath::Pi() / 2, 2 * TMath::Pi());
  TH1D* hTruthDeltaEta_all     = new TH1D("hTruthDeltaEta_all", "", 25, 0, 2.4);
  TH1D* hTruthDeltaEta_iso     = new TH1D("hTruthDeltaEta_iso", "", 25, 0, 2.4);

  // === TH2D Histograms: pt over phi ===
  TH2D* hTruthPt_PhiLead_all       = new TH2D("hTruthPt_PhiLead_all", "", 32, -TMath::Pi(), TMath::Pi(), pt_N, pt_bins);
  TH2D* hTruthPt_PhiLead_iso       = new TH2D("hTruthPt_PhiLead_iso", "", 32, -TMath::Pi(), TMath::Pi(), pt_N, pt_bins);
  TH2D* hTruthPt_PhiSubLead_all    = new TH2D("hTruthPt_PhiSubLead_all", "", 32, -TMath::Pi(), TMath::Pi(), pt_N, pt_bins);
  TH2D* hTruthPt_PhiSubLead_iso    = new TH2D("hTruthPt_PhiSubLead_iso", "", 32, -TMath::Pi(), TMath::Pi(), pt_N, pt_bins);
  TH2D* hTruthPt_PhiInclusive_all  = new TH2D("hTruthPt_PhiInclusive_all", "", 32, -TMath::Pi(), TMath::Pi(), pt_N, pt_bins);
  TH2D* hTruthPt_PhiInclusive_iso  = new TH2D("hTruthPt_PhiInclusive_iso", "", 32, -TMath::Pi(), TMath::Pi(), pt_N, pt_bins);

  // === TH2D Histograms: pt over eta ===
  TH2D* hTruthPt_EtaLead_all       = new TH2D("hTruthPt_EtaLead_all", "", 48, -1.2, 1.2, pt_N, pt_bins);
  TH2D* hTruthPt_EtaLead_iso       = new TH2D("hTruthPt_EtaLead_iso", "", 48, -1.2, 1.2, pt_N, pt_bins);
  TH2D* hTruthPt_EtaSubLead_all    = new TH2D("hTruthPt_EtaSubLead_all", "", 48, -1.2, 1.2, pt_N, pt_bins);
  TH2D* hTruthPt_EtaSubLead_iso    = new TH2D("hTruthPt_EtaSubLead_iso", "", 48, -1.2, 1.2, pt_N, pt_bins);
  TH2D* hTruthPt_EtaInclusive_all  = new TH2D("hTruthPt_EtaInclusive_all", "", 48, -1.2, 1.2, pt_N, pt_bins);
  TH2D* hTruthPt_EtaInclusive_iso  = new TH2D("hTruthPt_EtaInclusive_iso", "", 48, -1.2, 1.2, pt_N, pt_bins);

  // === TH2D Histogram: pt1 vs pt2 ===
  TH2D* h_pt1pt2_truth_all = new TH2D("h_pt1pt2_truth_all", ";p_{T,1};p_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D* h_pt1pt2_truth_iso = new TH2D("h_pt1pt2_truth_iso", ";p_{T,1};p_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH1D *h_xj_classical_truth_all = new TH1D("h_xj_classical_truth_all",";x_{J};1/N", nbins, ixj_bins);
    TH1D *h_xj_classical_truth_iso = new TH1D("h_xj_classical_truth_iso",";x_{J};1/N", nbins, ixj_bins);

   // Reco MC Histograms — ALL
    TH1D* hRecoPtLead_all         = new TH1D("hRecoPtLead_all",        "", pt_N, pt_bins);
    TH1D* hRecoPtSubLead_all      = new TH1D("hRecoPtSubLead_all",     "", pt_N, pt_bins);
    TH1D* hRecoPtInclusive_all    = new TH1D("hRecoPtInclusive_all",   "", pt_N, pt_bins);

    TH1D* hRecoEtaLead_all        = new TH1D("hRecoEtaLead_all",       "", 48, -1.2, 1.2);
    TH1D* hRecoEtaSubLead_all     = new TH1D("hRecoEtaSubLead_all",    "", 48, -1.2, 1.2);
    TH1D* hRecoEtaInclusive_all   = new TH1D("hRecoEtaInclusive_all",  "", 48, -1.2, 1.2);

    TH1D* hRecoPhiLead_all        = new TH1D("hRecoPhiLead_all",       "", 32, -TMath::Pi(), TMath::Pi());
    TH1D* hRecoPhiSubLead_all     = new TH1D("hRecoPhiSubLead_all",    "", 32, -TMath::Pi(), TMath::Pi());
    TH1D* hRecoPhiInclusive_all   = new TH1D("hRecoPhiInclusive_all",  "", 32, -TMath::Pi(), TMath::Pi());

    TH1D* hRecoDeltaPhi_all       = new TH1D("hRecoDeltaPhi_all",      "", 50, TMath::Pi()/2, 2*TMath::Pi());
    TH1D* hRecoDeltaEta_all       = new TH1D("hRecoDeltaEta_all",      "", 25, 0, 2.4);

    TH2D* hRecoPt_PhiLead_all     = new TH2D("hRecoPt_PhiLead_all",    "", 32, -TMath::Pi(), TMath::Pi(), pt_N, pt_bins);
    TH2D* hRecoPt_PhiSubLead_all  = new TH2D("hRecoPt_PhiSubLead_all", "", 32, -TMath::Pi(), TMath::Pi(), pt_N, pt_bins);
    TH2D* hRecoPt_PhiInclusive_all= new TH2D("hRecoPt_PhiInclusive_all","", 32, -TMath::Pi(), TMath::Pi(), pt_N, pt_bins);

    TH2D* hRecoPt_EtaLead_all     = new TH2D("hRecoPt_EtaLead_all",    "", 48, -1.2, 1.2, pt_N, pt_bins);
    TH2D* hRecoPt_EtaSubLead_all  = new TH2D("hRecoPt_EtaSubLead_all", "", 48, -1.2, 1.2, pt_N, pt_bins);
    TH2D* hRecoPt_EtaInclusive_all= new TH2D("hRecoPt_EtaInclusive_all","", 48, -1.2, 1.2, pt_N, pt_bins);

    TH2D* h_pt1pt2_reco_all       = new TH2D("h_pt1pt2_reco_all",      ";p_{T,1};p_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
    TH1D* h_xj_classical_reco_all = new TH1D("h_xj_classical_reco_all","", 50, 0, 1.2);


    // Reco MC Histograms — ISO
    TH1D* hRecoPtLead_iso         = new TH1D("hRecoPtLead_iso",        "", pt_N, pt_bins);
    TH1D* hRecoPtSubLead_iso      = new TH1D("hRecoPtSubLead_iso",     "", pt_N, pt_bins);
    TH1D* hRecoPtInclusive_iso    = new TH1D("hRecoPtInclusive_iso",   "", pt_N, pt_bins);

    TH1D* hRecoEtaLead_iso        = new TH1D("hRecoEtaLead_iso",       "", 48, -1.2, 1.2);
    TH1D* hRecoEtaSubLead_iso     = new TH1D("hRecoEtaSubLead_iso",    "", 48, -1.2, 1.2);
    TH1D* hRecoEtaInclusive_iso   = new TH1D("hRecoEtaInclusive_iso",  "", 48, -1.2, 1.2);

    TH1D* hRecoPhiLead_iso        = new TH1D("hRecoPhiLead_iso",       "", 32, -TMath::Pi(), TMath::Pi());
    TH1D* hRecoPhiSubLead_iso     = new TH1D("hRecoPhiSubLead_iso",    "", 32, -TMath::Pi(), TMath::Pi());
    TH1D* hRecoPhiInclusive_iso   = new TH1D("hRecoPhiInclusive_iso",  "", 32, -TMath::Pi(), TMath::Pi());

    TH1D* hRecoDeltaPhi_iso       = new TH1D("hRecoDeltaPhi_iso",      "", 50, TMath::Pi()/2, 2*TMath::Pi());
    TH1D* hRecoDeltaEta_iso       = new TH1D("hRecoDeltaEta_iso",      "", 25, 0, 2.4);

    TH2D* hRecoPt_PhiLead_iso     = new TH2D("hRecoPt_PhiLead_iso",    "", 32, -TMath::Pi(), TMath::Pi(), pt_N, pt_bins);
    TH2D* hRecoPt_PhiSubLead_iso  = new TH2D("hRecoPt_PhiSubLead_iso", "", 32, -TMath::Pi(), TMath::Pi(), pt_N, pt_bins);
    TH2D* hRecoPt_PhiInclusive_iso= new TH2D("hRecoPt_PhiInclusive_iso","", 32, -TMath::Pi(), TMath::Pi(), pt_N, pt_bins);

    TH2D* hRecoPt_EtaLead_iso     = new TH2D("hRecoPt_EtaLead_iso",    "", 48, -1.2, 1.2, pt_N, pt_bins);
    TH2D* hRecoPt_EtaSubLead_iso  = new TH2D("hRecoPt_EtaSubLead_iso", "", 48, -1.2, 1.2, pt_N, pt_bins);
    TH2D* hRecoPt_EtaInclusive_iso= new TH2D("hRecoPt_EtaInclusive_iso","", 48, -1.2, 1.2, pt_N, pt_bins);

    TH2D* h_pt1pt2_reco_iso       = new TH2D("h_pt1pt2_reco_iso",      ";p_{T,1};p_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
    TH1D* h_xj_classical_reco_iso = new TH1D("h_xj_classical_reco_iso","", 50, 0, 1.2);


    // identify file to determine setting
    int setting = 0;      float weight = 0; float tupper = 0; float tlower = 0;
  if (infile1.find("30") != std::string::npos){
    setting = 30;weight =  2.505*pow(10,-9);
    cout << "30 GEV MC JET INFILE" << endl;
    tupper = 200; tlower = 30;
      if(infile1.find("Herwig") != std::string::npos){
       weight =  weight = 1.473*pow(10,-12); cout << "Herwig detected" << endl;
    }

   }
   else if (infile1.find("20") != std::string::npos){
    setting = 20;weight = 6.218*pow(10,-8);
    cout << "20 GEV MC  JET INFILE" << endl;
    tupper = 30; tlower = 20;

   }
  else if (infile1.find("10") != std::string::npos){
    setting = 10;weight = 3.646*pow(10,-6);
    cout << "10 GEV MC  JET INFILE" << endl;
    tupper = 20; tlower = 10;
     if(infile1.find("Herwig") != std::string::npos){
       weight = 1.57028*pow(10,-8);tupper = 30;
      cout << "Herwig detected" << endl;
    }
    

   }
  else if (infile1.find("_ana") != std::string::npos){
    setting = 1;weight = 1;
    cout << "DATA JET INFILE" << endl;
    tupper = 200; tlower = 0;
   }
    //connect variables in the code to variables in the tree
   TFile* _file0 = TFile::Open(infile1.c_str());
   TTree* tree = (TTree*)_file0->Get("ttree");

   //define vectors and variables
   std::vector<double> *tjet_pt ={0}; std::vector<double> *tjet_phi ={0}; std::vector<double> *tjet_eta ={0};
   std::vector<double> *rjet_pt ={0}; std::vector<double> *rjet_phi ={0}; std::vector<double> *rjet_eta ={0};  float z_vtx = 0;float rz_vtx = {0};

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

    
       tree->SetBranchAddress("truth_vertex_z",&z_vtx); tree->SetBranchAddress("mbd_vertex_z",&rz_vtx);
       //  unsigned long long triggers;
      //tree->SetBranchAddress("gl1_scaled",&triggers);
    

   //define variables
   Double_t tleadpt = 0; float tsubleadpt = 0; float tleadphi = 0; float tsubleadphi = 0; 
   float tleadeta = 0; float tsubleadeta = 0;	float tdeltaphi = 0; float tdeltaeta = 0;

   Double_t rleadpt = 0; float rsubleadpt = 0; float rleadphi = 0; float rsubleadphi = 0; 
   float rleadeta = 0; float rsubleadeta = 0;

   float etacut = 1.1 - Rvalue; float entry = 0; float countmiss = 0; float countmatch = 0;

   float drminlead = 0; float drminsublead = 0;

   float isolated = 0; float total = 0;
   int Nentries = tree->GetEntries();
   //int Nentries = 100000;
    cout << "Nentries = : " << Nentries << endl;
    //start loop over entries
    for (int i = 0;i<Nentries;i++){

        //loop over entries
        tree->GetEntry(i);
	entry++; 


	int ntjets = tjet_pt->size();  int nrjets = rjet_pt->size(); float zvtx = z_vtx;float rzvtx = rz_vtx;

	
	
	if (ntjets < 2 || nrjets < 2 ||  fabs(zvtx) > 60 || fabs(rzvtx) > 60 ){
	   continue;
	 }

	 else if (ntjets >= 2 &&  fabs(zvtx) < 60 ){ //open conditional loop
	   
	    tleadpt = 0; tsubleadpt = 0; drminlead = 0.75*Rvalue; drminsublead = 0.75*Rvalue;
	    tleadeta = 0; tleadphi = 0; tsubleadeta = 0; tsubleadphi = 0;
	   for (int j = 0; j < ntjets;j++){ //find leading jet

	     	
              float  tjetpt = tjet_pt->at(j); float tjetphi =  tjet_phi->at(j);  float tjeteta = tjet_eta->at(j);

	      if (i < 20){
	      cout << "entry is: " << i << " leading loop is: " << j << " tjept is: " << tjetpt << " tjetphi is: " << tjetphi << " tjeteta is: " << tjeteta << endl;
	      }
	      //find leading jets
	        if (tjetpt > tleadpt){
			tleadpt = tjetpt;
			tleadphi = tjetphi;
			tleadeta = tjeteta;
		      }


		if (j == ntjets -1){//last entry

		   for (int l=0; l < nrjets;l++){ //loop over reco jets

		  
		     	       
		      float rjetpt = rjet_pt->at(l);float rjetphi = rjet_phi->at(l); float rjeteta = rjet_eta->at(l); 
		      float detalead = fabs(rjeteta - tleadeta); float dphilead = fabs(rjetphi - tleadphi);

		         if (i < 20){
			   cout << "entry is: " << i << " rleading loop is: " << l << " rjept is: " << rjetpt << " tleadpt is: " << tleadpt <<  " rjetphi is: " << rjetphi << " dphi is: " << dphilead <<  " tjeteta is: " << rjeteta << " detalead is: " << detalead << endl;
	      }
		  

		       if (dphilead > TMath::Pi())
			{
			  dphilead = 2*TMath::Pi()-dphilead;
			}

		     
		      float drlead = TMath::Sqrt(dphilead*dphilead + detalead*detalead);
		   
		      
		    
		     
		      if (drlead < drminlead){//match leading jets
			rleadpt = rjetpt; rleadphi = rjetphi; rleadeta = rjeteta;
			drminlead = drlead;
		      }//end match leading jets
		      /* if (i < 20){
	       
			cout << "Entry: " << i  << " drlead is: " << drlead  << " drminlead is " << drminlead<< endl;
			  
			}*/


		   }// end loop over reco jets


		}//end last entry

	    }// close leading jet loop

	   int candidate = 0;
	   for (int k = 0; k < ntjets; k++){//find subleading jets

	    
	      float  tjetpt = tjet_pt->at(k); float tjetphi =  tjet_phi->at(k);  float tjeteta = tjet_eta->at(k);
	      tdeltaphi = fabs(tleadphi - tjetphi); tdeltaeta = fabs(tleadeta - tjeteta);

	        if (i < 20){
	      cout << "entry is: " << i << " subleading loop is: " << k << " tjept is: " << tjetpt << " tjetphi is: " << tjetphi << " tjeteta is: " << tjeteta << endl;
	      }
	     

	       if ( tdeltaphi > TMath::Pi()){
			   tdeltaphi = 2*TMath::Pi()-  tdeltaphi;
			}

	       if(tjetpt > 5.5 && tdeltaphi > 3*TMath::Pi()/4 && fabs(tjeteta) < etacut){//candidate found
		 candidate++;

		 if (tjetpt > tsubleadpt){
		   tsubleadpt = tjetpt; tsubleadeta = tjeteta; tsubleadphi = tjetphi;

		 }
	       }//end candidate

	       	if (k == ntjets -1){//last entry

		   for (int l=0; l < nrjets;l++){ //loop over reco jets
		     	       
		      float rjetpt = rjet_pt->at(l);float rjetphi = rjet_phi->at(l); float rjeteta = rjet_eta->at(l); 
		      float detasublead = fabs(rjeteta - tsubleadeta); float dphisublead = fabs(rjetphi - tsubleadphi);
		      
		         if (i < 20){
			   cout << "entry is: " << i << " rsubleading loop is: " << l << " rjept is: " << rjetpt << " tsubleadpt is: " << tsubleadpt <<  " rjetphi is: " << rjetphi << " dphi is: " << dphisublead <<  " tjeteta is: " << rjeteta << " detasublead is: " << detasublead << endl;
	      }
		  

		       if (dphisublead > TMath::Pi())
			{
			  dphisublead = 2*TMath::Pi()-dphisublead;
			}

		     
		      float drsublead = TMath::Sqrt(dphisublead*dphisublead + detasublead*detasublead);
		   
		      
		    
		     
		      if (drsublead < drminsublead){//match leading jets
			rsubleadpt = rjetpt; rsubleadphi = rjetphi; rsubleadeta = rjeteta;
			drminsublead = drsublead;
		      }//end match leading jets

		      /*   if (i < 20){
		
			cout << "Entry: " << i  << " drsublead is: " << drsublead  << " drminsublead is " << drminsublead<< endl;
		      	  
			}*/
		     


		   }// end loop over reco jets


		}//end last entry
	      
	  
	   }////end subleading jets


	      bool truth_good = (tleadpt > 14  && tsubleadpt > 5.5  && tdeltaphi > 3*TMath::Pi()/4  && fabs(tleadeta) < etacut && fabs(tsubleadeta) < etacut && fabs(zvtx) < 60  && tleadpt >= tlower && tleadpt <= tupper);

	    if(!truth_good){
	      if (i < 40){
		cout << "Entry is: " << i <<   " tleadpt : " << tleadpt << " t subleadpt: " << tsubleadpt << " tdeltaphi: " << tdeltaphi << " tleadeta: " << tleadeta << " tsubleadeta " << tsubleadeta << endl;  

		}
	       tleadpt = 0; tsubleadpt = 0; drminlead = 0.75*Rvalue; drminsublead = 0.75*Rvalue;
	      continue;
	      }

	   

	       float rdeltaphi = fabs(rleadphi - rsubleadphi); float rdeltaeta = fabs(rleadeta - rsubleadeta);
	     

	       if (rdeltaphi > TMath::Pi()){
			   rdeltaphi = 2*TMath::Pi()-  rdeltaphi;
			}

	    
	      bool reco_good = (rleadpt > 18.28  && rsubleadpt > 8.2  && rdeltaphi > 3*TMath::Pi()/4  && fabs(rleadeta) < etacut && fabs(rleadeta) < etacut && fabs(rzvtx) < 60);
	    
	    if(truth_good && reco_good && candidate >= 1) { // fill all hists
	      float xjtruth = tsubleadpt / tleadpt;
	      float xjreco  = rsubleadpt / rleadpt;
	      total++;

	      h_xj_classical_truth_all->Fill(xjtruth, weight); h_xj_classical_reco_all->Fill(xjreco, weight);

	      // === 1D pt ===
	      hTruthPtLead_all->Fill(tleadpt, weight);      hTruthPtSubLead_all->Fill(tsubleadpt, weight);
	      hTruthPtInclusive_all->Fill(tleadpt, weight); hTruthPtInclusive_all->Fill(tsubleadpt, weight);
	      hRecoPtLead_all->Fill(rleadpt, weight);       hRecoPtSubLead_all->Fill(rsubleadpt, weight);
	      hRecoPtInclusive_all->Fill(rleadpt, weight);  hRecoPtInclusive_all->Fill(rsubleadpt, weight);

	      // === 1D eta ===
	      hTruthEtaLead_all->Fill(tleadeta, weight);    hTruthEtaSubLead_all->Fill(tsubleadeta, weight);
	      hTruthEtaInclusive_all->Fill(tleadeta, weight); hTruthEtaInclusive_all->Fill(tsubleadeta, weight);
	      hRecoEtaLead_all->Fill(rleadeta, weight);     hRecoEtaSubLead_all->Fill(rsubleadeta, weight);
	      hRecoEtaInclusive_all->Fill(rleadeta, weight); hRecoEtaInclusive_all->Fill(rsubleadeta, weight);

	      // === 1D phi ===
	      hTruthPhiLead_all->Fill(tleadphi, weight);    hTruthPhiSubLead_all->Fill(tsubleadphi, weight);
	      hTruthPhiInclusive_all->Fill(tleadphi, weight); hTruthPhiInclusive_all->Fill(tsubleadphi, weight);
	      hRecoPhiLead_all->Fill(rleadphi, weight);     hRecoPhiSubLead_all->Fill(rsubleadphi, weight);
	      hRecoPhiInclusive_all->Fill(rleadphi, weight); hRecoPhiInclusive_all->Fill(rsubleadphi, weight);

	      // === 1D dEta/dPhi ===
	      hTruthDeltaPhi_all->Fill(tdeltaphi, weight);  hTruthDeltaEta_all->Fill(tdeltaeta, weight);
	      hRecoDeltaPhi_all->Fill(rdeltaphi, weight);   hRecoDeltaEta_all->Fill(rdeltaeta, weight);

	      // === 2D pt vs phi ===
	      hTruthPt_PhiLead_all->Fill(tleadphi, tleadpt, weight);        hTruthPt_PhiSubLead_all->Fill(tsubleadphi, tsubleadpt, weight);
	      hTruthPt_PhiInclusive_all->Fill(tleadphi, tleadpt, weight);   hTruthPt_PhiInclusive_all->Fill(tsubleadphi, tsubleadpt, weight);
	      hRecoPt_PhiLead_all->Fill(rleadphi, rleadpt, weight);         hRecoPt_PhiSubLead_all->Fill(rsubleadphi, rsubleadpt, weight);
	      hRecoPt_PhiInclusive_all->Fill(rleadphi, rleadpt, weight);    hRecoPt_PhiInclusive_all->Fill(rsubleadphi, rsubleadpt, weight);

	      // === 2D pt vs eta ===
	      hTruthPt_EtaLead_all->Fill(tleadeta, tleadpt, weight);        hTruthPt_EtaSubLead_all->Fill(tsubleadeta, tsubleadpt, weight);
	      hTruthPt_EtaInclusive_all->Fill(tleadeta, tleadpt, weight);   hTruthPt_EtaInclusive_all->Fill(tsubleadeta, tsubleadpt, weight);
	      hRecoPt_EtaLead_all->Fill(rleadeta, rleadpt, weight);         hRecoPt_EtaSubLead_all->Fill(rsubleadeta, rsubleadpt, weight);
	      hRecoPt_EtaInclusive_all->Fill(rleadeta, rleadpt, weight);    hRecoPt_EtaInclusive_all->Fill(rsubleadeta, rsubleadpt, weight);

	      // === pt1 vs pt2 ===
	      h_pt1pt2_truth_all->Fill(tleadpt, tsubleadpt, weight);        h_pt1pt2_truth_all->Fill(tsubleadpt, tleadpt, weight);
	      h_pt1pt2_reco_all->Fill(rleadpt, rsubleadpt, weight);         h_pt1pt2_reco_all->Fill(rsubleadpt, rleadpt, weight);

	      if(candidate == 1) { // fill iso hists
		isolated++;

		h_xj_classical_truth_iso->Fill(xjtruth, weight); h_xj_classical_reco_iso->Fill(xjreco, weight);

		// === 1D pt ===
		hTruthPtLead_iso->Fill(tleadpt, weight);      hTruthPtSubLead_iso->Fill(tsubleadpt, weight);
		hTruthPtInclusive_iso->Fill(tleadpt, weight); hTruthPtInclusive_iso->Fill(tsubleadpt, weight);
		hRecoPtLead_iso->Fill(rleadpt, weight);       hRecoPtSubLead_iso->Fill(rsubleadpt, weight);
		hRecoPtInclusive_iso->Fill(rleadpt, weight);  hRecoPtInclusive_iso->Fill(rsubleadpt, weight);

		// === 1D eta ===
		hTruthEtaLead_iso->Fill(tleadeta, weight);    hTruthEtaSubLead_iso->Fill(tsubleadeta, weight);
		hTruthEtaInclusive_iso->Fill(tleadeta, weight); hTruthEtaInclusive_iso->Fill(tsubleadeta, weight);
		hRecoEtaLead_iso->Fill(rleadeta, weight);     hRecoEtaSubLead_iso->Fill(rsubleadeta, weight);
		hRecoEtaInclusive_iso->Fill(rleadeta, weight); hRecoEtaInclusive_iso->Fill(rsubleadeta, weight);

		// === 1D phi ===
		hTruthPhiLead_iso->Fill(tleadphi, weight);    hTruthPhiSubLead_iso->Fill(tsubleadphi, weight);
		hTruthPhiInclusive_iso->Fill(tleadphi, weight); hTruthPhiInclusive_iso->Fill(tsubleadphi, weight);
		hRecoPhiLead_iso->Fill(rleadphi, weight);     hRecoPhiSubLead_iso->Fill(rsubleadphi, weight);
		hRecoPhiInclusive_iso->Fill(rleadphi, weight); hRecoPhiInclusive_iso->Fill(rsubleadphi, weight);

		// === 1D dEta/dPhi ===
		hTruthDeltaPhi_iso->Fill(tdeltaphi, weight);  hTruthDeltaEta_iso->Fill(tdeltaeta, weight);
		hRecoDeltaPhi_iso->Fill(rdeltaphi, weight);   hRecoDeltaEta_iso->Fill(rdeltaeta, weight);

		// === 2D pt vs phi ===
		hTruthPt_PhiLead_iso->Fill(tleadphi, tleadpt, weight);        hTruthPt_PhiSubLead_iso->Fill(tsubleadphi, tsubleadpt, weight);
		hTruthPt_PhiInclusive_iso->Fill(tleadphi, tleadpt, weight);   hTruthPt_PhiInclusive_iso->Fill(tsubleadphi, tsubleadpt, weight);
		hRecoPt_PhiLead_iso->Fill(rleadphi, rleadpt, weight);         hRecoPt_PhiSubLead_iso->Fill(rsubleadphi, rsubleadpt, weight);
		hRecoPt_PhiInclusive_iso->Fill(rleadphi, rleadpt, weight);    hRecoPt_PhiInclusive_iso->Fill(rsubleadphi, rsubleadpt, weight);

		// === 2D pt vs eta ===
		hTruthPt_EtaLead_iso->Fill(tleadeta, tleadpt, weight);        hTruthPt_EtaSubLead_iso->Fill(tsubleadeta, tsubleadpt, weight);
		hTruthPt_EtaInclusive_iso->Fill(tleadeta, tleadpt, weight);   hTruthPt_EtaInclusive_iso->Fill(tsubleadeta, tsubleadpt, weight);
		hRecoPt_EtaLead_iso->Fill(rleadeta, rleadpt, weight);         hRecoPt_EtaSubLead_iso->Fill(rsubleadeta, rsubleadpt, weight);
		hRecoPt_EtaInclusive_iso->Fill(rleadeta, rleadpt, weight);    hRecoPt_EtaInclusive_iso->Fill(rsubleadeta, rsubleadpt, weight);

		// === pt1 vs pt2 ===
		h_pt1pt2_truth_iso->Fill(tleadpt, tsubleadpt, weight);        h_pt1pt2_truth_iso->Fill(tsubleadpt, tleadpt, weight);
		h_pt1pt2_reco_iso->Fill(rleadpt, rsubleadpt, weight);         h_pt1pt2_reco_iso->Fill(rsubleadpt, rleadpt, weight);
	      }//iso hists
	    }//all hists

	     if(truth_good && !reco_good && candidate >= 1) { // fill all truth  hists
	      float xjtruth = tsubleadpt / tleadpt;
	     
	      h_xj_classical_truth_all->Fill(xjtruth, weight); 

	      // === 1D pt ===
	      hTruthPtLead_all->Fill(tleadpt, weight);      hTruthPtSubLead_all->Fill(tsubleadpt, weight);
	      hTruthPtInclusive_all->Fill(tleadpt, weight); hTruthPtInclusive_all->Fill(tsubleadpt, weight);
	    
	      // === 1D eta ===
	      hTruthEtaLead_all->Fill(tleadeta, weight);    hTruthEtaSubLead_all->Fill(tsubleadeta, weight);
	      hTruthEtaInclusive_all->Fill(tleadeta, weight); hTruthEtaInclusive_all->Fill(tsubleadeta, weight);
	 
	      // === 1D phi ===
	      hTruthPhiLead_all->Fill(tleadphi, weight);    hTruthPhiSubLead_all->Fill(tsubleadphi, weight);
	      hTruthPhiInclusive_all->Fill(tleadphi, weight); hTruthPhiInclusive_all->Fill(tsubleadphi, weight);
	    

	      // === 1D dEta/dPhi ===
	      hTruthDeltaPhi_all->Fill(tdeltaphi, weight);  hTruthDeltaEta_all->Fill(tdeltaeta, weight);
	  
	      // === 2D pt vs phi ===
	      hTruthPt_PhiLead_all->Fill(tleadphi, tleadpt, weight);        hTruthPt_PhiSubLead_all->Fill(tsubleadphi, tsubleadpt, weight);
	      hTruthPt_PhiInclusive_all->Fill(tleadphi, tleadpt, weight);   hTruthPt_PhiInclusive_all->Fill(tsubleadphi, tsubleadpt, weight);
	  

	      // === 2D pt vs eta ===
	      hTruthPt_EtaLead_all->Fill(tleadeta, tleadpt, weight);        hTruthPt_EtaSubLead_all->Fill(tsubleadeta, tsubleadpt, weight);
	      hTruthPt_EtaInclusive_all->Fill(tleadeta, tleadpt, weight);   hTruthPt_EtaInclusive_all->Fill(tsubleadeta, tsubleadpt, weight);
	  

	      // === pt1 vs pt2 ===
	      h_pt1pt2_truth_all->Fill(tleadpt, tsubleadpt, weight);        h_pt1pt2_truth_all->Fill(tsubleadpt, tleadpt, weight);
	   

	      if(candidate == 1) { // fill iso hists
		isolated++;

		h_xj_classical_truth_iso->Fill(xjtruth, weight);

		// === 1D pt ===
		hTruthPtLead_iso->Fill(tleadpt, weight);      hTruthPtSubLead_iso->Fill(tsubleadpt, weight);
		hTruthPtInclusive_iso->Fill(tleadpt, weight); hTruthPtInclusive_iso->Fill(tsubleadpt, weight);
	

		// === 1D eta ===
		hTruthEtaLead_iso->Fill(tleadeta, weight);    hTruthEtaSubLead_iso->Fill(tsubleadeta, weight);
		hTruthEtaInclusive_iso->Fill(tleadeta, weight); hTruthEtaInclusive_iso->Fill(tsubleadeta, weight);
	

		// === 1D phi ===
		hTruthPhiLead_iso->Fill(tleadphi, weight);    hTruthPhiSubLead_iso->Fill(tsubleadphi, weight);
		hTruthPhiInclusive_iso->Fill(tleadphi, weight); hTruthPhiInclusive_iso->Fill(tsubleadphi, weight);
	

		// === 1D dEta/dPhi ===
		hTruthDeltaPhi_iso->Fill(tdeltaphi, weight);  hTruthDeltaEta_iso->Fill(tdeltaeta, weight);
	

		// === 2D pt vs phi ===
		hTruthPt_PhiLead_iso->Fill(tleadphi, tleadpt, weight);        hTruthPt_PhiSubLead_iso->Fill(tsubleadphi, tsubleadpt, weight);
		hTruthPt_PhiInclusive_iso->Fill(tleadphi, tleadpt, weight);   hTruthPt_PhiInclusive_iso->Fill(tsubleadphi, tsubleadpt, weight);	

		// === 2D pt vs eta ===
		hTruthPt_EtaLead_iso->Fill(tleadeta, tleadpt, weight);        hTruthPt_EtaSubLead_iso->Fill(tsubleadeta, tsubleadpt, weight);
		hTruthPt_EtaInclusive_iso->Fill(tleadeta, tleadpt, weight);   hTruthPt_EtaInclusive_iso->Fill(tsubleadeta, tsubleadpt, weight);

		// === pt1 vs pt2 ===
		h_pt1pt2_truth_iso->Fill(tleadpt, tsubleadpt, weight);        h_pt1pt2_truth_iso->Fill(tsubleadpt, tleadpt, weight);
       
	      }//iso hists
	    }//all hists



	  
	 }// close conditional loop


    }// end loop over entries

    cout << "number of isolated jets: " << isolated << " number of split jets: " << total - isolated << endl;
    cout << "proportion of isolated jets: " << isolated*100/(total) << " % proportion of split jets: " << (total - isolated)*100/(total) << endl;

    // cout << "number of misses: " << countmiss << " number of matches: " << countmatch << " percent missed: " << countmiss*100/(countmiss + countmatch) << "%" <<  endl;
    //Save hists to file for later
     TString outfilename = infile1;
     outfilename.Prepend("Hists_Iso_Split_unmatchtruth_R" + Rvaluestring +"_");
     TFile *outfile = TFile::Open(outfilename,"RECREATE");

    

     // === Write TH1D Histograms ===
     hTruthPtLead_all->Write();       hTruthPtLead_iso->Write();
     hTruthPtSubLead_all->Write();    hTruthPtSubLead_iso->Write();
     hTruthPtInclusive_all->Write();  hTruthPtInclusive_iso->Write();

     hTruthEtaLead_all->Write();      hTruthEtaLead_iso->Write();
     hTruthEtaSubLead_all->Write();   hTruthEtaSubLead_iso->Write();
     hTruthEtaInclusive_all->Write(); hTruthEtaInclusive_iso->Write();

     hTruthPhiLead_all->Write();      hTruthPhiLead_iso->Write();
     hTruthPhiSubLead_all->Write();   hTruthPhiSubLead_iso->Write();
     hTruthPhiInclusive_all->Write(); hTruthPhiInclusive_iso->Write();

     hTruthDeltaPhi_all->Write();     hTruthDeltaPhi_iso->Write();
     hTruthDeltaEta_all->Write();     hTruthDeltaEta_iso->Write();

     // === Write TH2D Histograms (pt over phi) ===
     hTruthPt_PhiLead_all->Write();       hTruthPt_PhiLead_iso->Write();
     hTruthPt_PhiSubLead_all->Write();    hTruthPt_PhiSubLead_iso->Write();
     hTruthPt_PhiInclusive_all->Write();  hTruthPt_PhiInclusive_iso->Write();

     // === Write TH2D Histograms (pt over eta) ===
     hTruthPt_EtaLead_all->Write();       hTruthPt_EtaLead_iso->Write();
     hTruthPt_EtaSubLead_all->Write();    hTruthPt_EtaSubLead_iso->Write();
     hTruthPt_EtaInclusive_all->Write();  hTruthPt_EtaInclusive_iso->Write();

     // === Write pt1 vs pt2 Histograms ===
     h_pt1pt2_truth_all->Write();
     h_pt1pt2_truth_iso->Write();


    h_xj_classical_truth_all->Write();
    h_xj_classical_truth_iso->Write();


     // === Write TH1D Histograms ===
     hRecoPtLead_all->Write();       hRecoPtLead_iso->Write();
     hRecoPtSubLead_all->Write();    hRecoPtSubLead_iso->Write();
     hRecoPtInclusive_all->Write();  hRecoPtInclusive_iso->Write();

     hRecoEtaLead_all->Write();      hRecoEtaLead_iso->Write();
     hRecoEtaSubLead_all->Write();   hRecoEtaSubLead_iso->Write();
     hRecoEtaInclusive_all->Write(); hRecoEtaInclusive_iso->Write();

     hRecoPhiLead_all->Write();      hRecoPhiLead_iso->Write();
     hRecoPhiSubLead_all->Write();   hRecoPhiSubLead_iso->Write();
     hRecoPhiInclusive_all->Write(); hRecoPhiInclusive_iso->Write();

     hRecoDeltaPhi_all->Write();     hRecoDeltaPhi_iso->Write();
     hRecoDeltaEta_all->Write();     hRecoDeltaEta_iso->Write();

     // === Write TH2D Histograms (pt over phi) ===
     hRecoPt_PhiLead_all->Write();       hRecoPt_PhiLead_iso->Write();
     hRecoPt_PhiSubLead_all->Write();    hRecoPt_PhiSubLead_iso->Write();
     hRecoPt_PhiInclusive_all->Write();  hRecoPt_PhiInclusive_iso->Write();

     // === Write TH2D Histograms (pt over eta) ===
     hRecoPt_EtaLead_all->Write();       hRecoPt_EtaLead_iso->Write();
     hRecoPt_EtaSubLead_all->Write();    hRecoPt_EtaSubLead_iso->Write();
     hRecoPt_EtaInclusive_all->Write();  hRecoPt_EtaInclusive_iso->Write();

     // === Write pt1 vs pt2 Histograms ===
     h_pt1pt2_reco_all->Write();
     h_pt1pt2_reco_iso->Write();


    h_xj_classical_reco_all->Write();
    h_xj_classical_reco_iso->Write();
}



