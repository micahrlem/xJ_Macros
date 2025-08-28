void Make_QA_Herwig_MultiR(string infile1 = "pythia8-Jet30-Run21-multiR.root", const std::string configfile = "binning_original.config",float Rvalue = 0.4){
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
    // pt
    TH1D* hTruthPtLead  = new TH1D("hTruthPtLead ","",pt_N,pt_bins); 
    TH1D* hTruthPtSubLead  = new TH1D("hTruthPtSubLead ","",pt_N,pt_bins); 
    TH1D* hTruthPtInclusive  = new TH1D("hTruthPtInclusive ","",pt_N,pt_bins);

    // eta
    TH1D* hTruthEtaLead  = new TH1D("hTruthEtaLead ","",48,-1.2,1.2); 
    TH1D* hTruthEtaSubLead  = new TH1D("hTruthEtaSubLead ","",48,-1.2,1.2);  
    TH1D* hTruthEtaInclusive  = new TH1D("hTruthEtaInclusive ","",48,-1.2,1.2);

    // phi
    TH1D* hTruthPhiLead  = new TH1D("hTruthPhiLead ","",32,-TMath::Pi(),TMath::Pi()); 
    TH1D* hTruthPhiSubLead  = new TH1D("hTruthPhiSubLead ","",32,-TMath::Pi(),TMath::Pi());  
    TH1D* hTruthPhiInclusive  = new TH1D("hTruthPhiInclusive ","",32,-TMath::Pi(),TMath::Pi());

    // dphi and deta
    TH1D* hTruthDeltaPhi  = new TH1D("hTruthDeltaPhi ","",50,TMath::Pi()/2,2*TMath::Pi());
    TH1D* hTruthDeltaEta  = new TH1D("hTruthDeltaEta ","",25,0,2.4);

    // pt over eta and phi
    TH2D* hTruthPt_PhiLead  = new TH2D("hTruthPt_PhiLead ","",32,-TMath::Pi(),TMath::Pi(),pt_N,pt_bins); 
    TH2D* hTruthPt_PhiSubLead  = new TH2D("hTruthPt_PhiSubLead ","",32,-TMath::Pi(),TMath::Pi(),pt_N,pt_bins);  
    TH2D* hTruthPt_PhiInclusive  = new TH2D("hTruthPt_PhiInclusive ","",32,-TMath::Pi(),TMath::Pi(),pt_N,pt_bins);

    TH2D* hTruthPt_EtaLead  = new TH2D("hTruthPt_EtaLead ","",48,-1.2,1.2,pt_N,pt_bins); 
    TH2D* hTruthPt_EtaSubLead  = new TH2D("hTruthPt_EtaSubLead ","",48,-1.2,1.2,pt_N,pt_bins);  
    TH2D* hTruthPt_EtaInclusive  = new TH2D("hTruthPt_EtaInclusive ","",48,-1.2,1.2,pt_N,pt_bins);

    //xj classical and pt1pt2
    TH2D *h_pt1pt2_truth =  new TH2D("h_pt1pt2_truth",";p_{T,1};p_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
    TH1D *h_xj_classical_truth = new TH1D("h_xj_classical_truth",";x_{J};1/N", nbins, ixj_bins);

    // Define Reco MC Histos
    TH1D* hRecoPtLead  = new TH1D("hRecoPtLead ","",pt_N,pt_bins); 
    TH1D* hRecoPtSubLead  = new TH1D("hRecoPtSubLead ","",pt_N,pt_bins);  
    TH1D* hRecoPtInclusive  = new TH1D("hRecoPtInclusive ","",pt_N,pt_bins);

    // eta
    TH1D* hRecoEtaLead  = new TH1D("hRecoEtaLead ","",48,-1.2,1.2); 
    TH1D* hRecoEtaSubLead  = new TH1D("hRecoEtaSubLead ","",48,-1.2,1.2);  
    TH1D* hRecoEtaInclusive  = new TH1D("hRecoEtaInclusive ","",48,-1.2,1.2);

    // phi
    TH1D* hRecoPhiLead  = new TH1D("hRecoPhiLead ","",32,-TMath::Pi(),TMath::Pi()); 
    TH1D* hRecoPhiSubLead  = new TH1D("hRecoPhiSubLead ","",32,-TMath::Pi(),TMath::Pi());  
    TH1D* hRecoPhiInclusive  = new TH1D("hRecoPhiInclusive ","",32,-TMath::Pi(),TMath::Pi());

    // dphi and deta
    TH1D* hRecoDeltaPhi  = new TH1D("hRecoDeltaPhi ","",50,TMath::Pi()/2,2*TMath::Pi());
    TH1D* hRecoDeltaEta  = new TH1D("hRecoDeltaEta ","",25,0,2.4);

    // pt over eta and phi
    TH2D* hRecoPt_PhiLead  = new TH2D("hRecoPt_PhiLead ","",32,-TMath::Pi(),TMath::Pi(),pt_N,pt_bins); 
    TH2D* hRecoPt_PhiSubLead  = new TH2D("hRecoPt_PhiSubLead ","",32,-TMath::Pi(),TMath::Pi(),pt_N,pt_bins);  
    TH2D* hRecoPt_PhiInclusive  = new TH2D("hRecoPt_PhiInclusive ","",32,-TMath::Pi(),TMath::Pi(),pt_N,pt_bins);

    TH2D* hRecoPt_EtaLead  = new TH2D("hRecoPt_EtaLead ","",48,-1.2,1.2,pt_N,pt_bins); 
    TH2D* hRecoPt_EtaSubLead  = new TH2D("hRecoPt_EtaSubLead ","",48,-1.2,1.2,pt_N,pt_bins);  
    TH2D* hRecoPt_EtaInclusive  = new TH2D("hRecoPt_EtaInclusive ","",48,-1.2,1.2,pt_N,pt_bins);

   
    //xj classical and pt1pt2
    TH2D *h_pt1pt2_reco =  new TH2D("h_pt1pt2_reco",";p_{T,1};p_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
    TH1D *h_xj_classical_reco = new TH1D("h_xj_classical_reco",";x_{J};1/N", nbins, ixj_bins);


    // identify file to determine setting
    int setting = 0;      float weight = 0; float tupper = 0; float tlower = 0;
  if (infile1.find("30") != std::string::npos){
    setting = 30;weight =  2.505*pow(10,-9);
    cout << "30 GEV MC JET INFILE" << endl;
    tupper = 200; tlower = 30;

   }
  else if (infile1.find("10") != std::string::npos){
    setting = 10;weight = 3.646*pow(10,-6);
    cout << "10 GEV MC  JET INFILE" << endl;
    tupper = 30; tlower = 10;

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
   float tleadeta = 0; float tsubleadeta = 0;

   Double_t rleadpt = 0; float rsubleadpt = 0; float rleadphi = 0; float rsubleadphi = 0; 
   float rleadeta = 0; float rsubleadeta = 0;

   float etacut = 1.1 - Rvalue; float entry = 0; float countmiss = 0; float countmatch = 0;
    float drmin_lead =  Rvalue*0.75; float drmin_sublead =  Rvalue*0.75;
   int Nentries = tree->GetEntries();
   //int Nentries = 100000;
    cout << "Nentries = : " << Nentries << endl;
    //start loop over entries
    for (int i = 0;i<Nentries;i++){

        //loop over entries
        tree->GetEntry(i);
	entry++; 


	int ntjets = tjet_pt->size();  int nrjets = rjet_pt->size();
	
	 if (ntjets < 2 || nrjets < 2){
	   continue;
	 }

	 else if (ntjets >= 2 && nrjets >= 2){ //open conditional loop
	           tleadpt = 0; tsubleadpt = 0;  drmin_lead =  Rvalue*0.75; drmin_sublead =  Rvalue*0.75;
	    for (int j = 0; j < ntjets;j++){ //open ntjets loop

	   
              float  tjetpt = tjet_pt->at(j); float tjetphi =  tjet_phi->at(j);  float tjeteta = tjet_eta->at(j);
	      
	  
	      Double_t reco_matched_leadpt = 0;float reco_matched_leadphi = 0;float reco_matched_leadeta = 0;
	      float reco_matched_subleadpt = 0; float reco_matched_subleadphi = 0; float reco_matched_subleadeta = 0;
	      float tdeltaphi = 0; float tdeltaeta = 0;  float rdeltaphi = 0; float rdeltaeta = 0;
	      float zvtx = z_vtx; float rzvtx = rz_vtx;
	    
	     
	      //find leading jets
	        if (tjetpt > tleadpt){
			tsubleadpt = tleadpt;	tleadpt = tjetpt;
			tsubleadphi = tleadphi;	tleadphi = tjetphi;
	        	tsubleadeta = tleadeta;tleadeta = tjeteta;
		      }
		  
		 else if (tjetpt > tsubleadpt){
			tsubleadpt = tjetpt;	tsubleadphi = tjetphi;	tsubleadeta = tjeteta;
		      }


		  if (j == (ntjets - 1)){//last entry of event

		     for (int k=0; k < nrjets;k++){ //open nrjets loop
		       
	       		       
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
		      
		      //match jets
		      if (drlead < drmin_lead){
			reco_matched_leadpt = rjetpt; reco_matched_leadphi = rjetphi; reco_matched_leadeta = rjeteta;
		

			drmin_lead = drlead;
		      }
		      if (drsublead < drmin_sublead){
			reco_matched_subleadpt = rjetpt;reco_matched_subleadphi = rjetphi;reco_matched_subleadeta = rjeteta;
		

			drmin_sublead = drsublead;

		      }
		     
		     }//close nrjets loop

		   
		     tdeltaphi = fabs(tleadphi-tsubleadphi);tdeltaeta = fabs(tleadeta-tsubleadeta);
		     rdeltaphi = fabs(reco_matched_leadphi - reco_matched_subleadphi); 
		     if ( tdeltaphi > TMath::Pi())
			{
			   tdeltaphi = 2*TMath::Pi()-  tdeltaphi;
			}

		       if ( rdeltaphi > TMath::Pi())
			{
			  rdeltaphi = 2*TMath::Pi()- rdeltaphi;
			}
		     rdeltaeta = fabs(reco_matched_leadeta - reco_matched_subleadeta);
		      
		     bool truth_good = (tleadpt > 14  && tsubleadpt > 7  && tdeltaphi > 3*TMath::Pi()/4  && fabs(tleadeta) < etacut && fabs(tsubleadeta) < etacut && fabs(zvtx) < 60  && tleadpt >= tlower && tleadpt <= tupper);

		       bool reco_good = (reco_matched_leadpt > 18.28  && reco_matched_subleadpt > 8.2  && rdeltaphi > 3*TMath::Pi()/4  && fabs(reco_matched_leadeta) < etacut && fabs(reco_matched_subleadeta) < etacut && fabs(rzvtx) < 60);
		
		     //all cuts on truth loop
		     if(truth_good && reco_good){

		       //Fill  Histograms
		       float xjtruth = tsubleadpt/tleadpt; float xjreco = reco_matched_subleadpt/reco_matched_leadpt;
		
			 //zvtx 60 Histos
			     //truth pt
			   hTruthPtLead ->Fill(tleadpt,weight); hTruthPtSubLead ->Fill(tsubleadpt,weight);
			   hTruthPtInclusive ->Fill(tleadpt,weight); hTruthPtInclusive ->Fill(tsubleadpt,weight);
			     // truth eta
			   hTruthEtaLead ->Fill(tleadeta,weight);hTruthEtaSubLead ->Fill(tsubleadeta,weight);
			   hTruthEtaInclusive ->Fill(tleadeta,weight);hTruthEtaInclusive ->Fill(tsubleadeta,weight);
			     //truth phi
			   hTruthPhiLead ->Fill(tleadphi,weight); hTruthPhiSubLead ->Fill(tsubleadphi,weight); 
			   hTruthPhiInclusive ->Fill(tleadphi,weight);hTruthPhiInclusive ->Fill(tsubleadphi,weight);
			     //truth dphi and deta
			   hTruthDeltaPhi ->Fill(tdeltaphi,weight); hTruthDeltaEta ->Fill(tdeltaeta,weight);

			     //Reco pt
			   hRecoPtLead ->Fill(reco_matched_leadpt,weight); hRecoPtSubLead ->Fill(reco_matched_subleadpt,weight);
			   hRecoPtInclusive ->Fill(reco_matched_leadpt,weight); hRecoPtInclusive ->Fill(reco_matched_subleadpt,weight);
			     // Reco eta
			   hRecoEtaLead ->Fill(reco_matched_leadeta,weight);hRecoEtaSubLead ->Fill(reco_matched_subleadeta,weight);
			   hRecoEtaInclusive ->Fill(reco_matched_leadeta,weight);hRecoEtaInclusive ->Fill(reco_matched_subleadeta,weight);
			     //Reco phi
			   hRecoPhiLead ->Fill(reco_matched_leadphi,weight); hRecoPhiSubLead ->Fill(reco_matched_subleadphi,weight); 
			   hRecoPhiInclusive ->Fill(reco_matched_leadphi,weight);hRecoPhiInclusive ->Fill(reco_matched_subleadphi,weight);
			     //Reco dphi and deta
			   hRecoDeltaPhi ->Fill(rdeltaphi,weight); hRecoDeltaEta ->Fill(rdeltaeta,weight);
	        

			   // Filling 2D histograms for 60vtx
			   hTruthPt_PhiLead ->Fill(tleadphi, tleadpt, weight);
			   hTruthPt_EtaLead ->Fill(tleadeta, tleadpt, weight);
			   hTruthPt_PhiSubLead ->Fill(tsubleadphi, tsubleadpt, weight);
			   hTruthPt_EtaSubLead ->Fill(tsubleadeta, tsubleadpt, weight);
			   hTruthPt_PhiInclusive ->Fill(tleadphi, tleadpt, weight);
			   hTruthPt_PhiInclusive ->Fill(tsubleadphi, tsubleadpt, weight);
			   hTruthPt_EtaInclusive ->Fill(tleadeta, tleadpt, weight);
			   hTruthPt_EtaInclusive ->Fill(tsubleadeta, tsubleadpt, weight);

			   hRecoPt_PhiLead ->Fill(reco_matched_leadphi, reco_matched_leadpt, weight);
			   hRecoPt_EtaLead ->Fill(reco_matched_leadeta, reco_matched_leadpt, weight);
			   hRecoPt_PhiSubLead ->Fill(reco_matched_subleadphi, reco_matched_subleadpt, weight);
			   hRecoPt_EtaSubLead ->Fill(reco_matched_subleadeta, reco_matched_subleadpt, weight);
			   hRecoPt_PhiInclusive ->Fill(reco_matched_leadphi, reco_matched_leadpt, weight);
			   hRecoPt_PhiInclusive ->Fill(reco_matched_subleadphi, reco_matched_subleadpt, weight);
			   hRecoPt_EtaInclusive ->Fill(reco_matched_leadeta, reco_matched_leadpt, weight);
			   hRecoPt_EtaInclusive ->Fill(reco_matched_subleadeta, reco_matched_subleadpt, weight);

			   //xj and pt1pt2
			   h_pt1pt2_truth->Fill(tleadpt,tsubleadpt,weight);h_pt1pt2_truth->Fill(tsubleadpt,tleadpt,weight);
			   h_pt1pt2_reco->Fill(reco_matched_leadpt,reco_matched_subleadpt,weight);h_pt1pt2_reco->Fill(reco_matched_subleadpt,reco_matched_leadpt,weight);
			   h_xj_classical_truth->Fill(xjtruth,weight); h_xj_classical_reco->Fill(xjreco,weight);



		     

		     }// end all cuts loop

		   
		  }// end last entry loop
		     

	    }// close ntjets loop
	   
	 }// close conditional loop

    }// end loop over entries


    // cout << "number of misses: " << countmiss << " number of matches: " << countmatch << " percent missed: " << countmiss*100/(countmiss + countmatch) << "%" <<  endl;
    //Save hists to file for later
     TString outfilename = infile1;
     outfilename.Prepend("Hists_QA_R" + Rvaluestring +"_");
     TFile *outfile = TFile::Open(outfilename,"RECREATE");

     // Write Truth Histos
     hTruthPtLead ->Write();
     hTruthPtSubLead ->Write();
     hTruthPtInclusive ->Write();

     hTruthEtaLead ->Write();
     hTruthEtaSubLead ->Write();
     hTruthEtaInclusive ->Write();

     hTruthPhiLead ->Write();
     hTruthPhiSubLead ->Write();
     hTruthPhiInclusive ->Write();

     hTruthDeltaPhi ->Write();
     hTruthDeltaEta ->Write();

     // Write Reco Histos
     hRecoPtLead ->Write();
     hRecoPtSubLead ->Write();
     hRecoPtInclusive ->Write();

     hRecoEtaLead ->Write();
     hRecoEtaSubLead ->Write();
     hRecoEtaInclusive ->Write();

     hRecoPhiLead ->Write();
     hRecoPhiSubLead ->Write();
     hRecoPhiInclusive ->Write();

     hRecoDeltaPhi ->Write();
     hRecoDeltaEta ->Write();

     // Writing 2D histograms to file (only 60vtx)
     hTruthPt_PhiLead ->Write();
     hTruthPt_PhiSubLead ->Write();
     hTruthPt_PhiInclusive ->Write();
     hTruthPt_EtaLead ->Write();
     hTruthPt_EtaSubLead ->Write();
     hTruthPt_EtaInclusive ->Write();

     hRecoPt_PhiLead ->Write();
     hRecoPt_PhiSubLead ->Write();
     hRecoPt_PhiInclusive ->Write();
     hRecoPt_EtaLead ->Write();
     hRecoPt_EtaSubLead ->Write();
     hRecoPt_EtaInclusive ->Write();

     h_pt1pt2_truth->Write();
     h_pt1pt2_reco->Write();
     h_xj_classical_truth->Write(); h_xj_classical_reco->Write();

}



