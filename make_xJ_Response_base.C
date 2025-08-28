R__LOAD_LIBRARY(/sphenix/user/mmeskowit/Dijet_Analysis/MDC2_Analysis/Xj_Analysis/roounfold/timroounfold/libRooUnfold.so)

#include "xj_functions.h"
#include "read_binning.h"




#define _USE_MATH_DEFINES

#include <math.h> 
#include <cmath>
#include <iostream>

#include<vector>
#include<array>

void make_xJ_Response_base(string infile1 = "Herwig-Jet30-Run21-multiR.root", const std::string configfile = "binning_original.config", float RValue = 0.4){

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

     // Truth No Cut QA histograms
     TH1D *h_pt_lead_truth_nocut = new TH1D("h_pt_lead_truth_nocut"," ; Leading Jet p_{T} [GeV]; counts",pt_N,pt_bins);
     TH1D *h_eta_lead_truth_nocut = new TH1D("h_eta_lead_truth_nocut"," ; Leading Jet #eta  ; counts",48,-1.2,1.2);
     TH1D *h_phi_lead_truth_nocut = new TH1D("h_phi_lead_truth_nocut"," ; Leading Jet #Phi  ; counts",32,-TMath::Pi(),TMath::Pi());

     TH1D *h_pt_sublead_truth_nocut = new TH1D("h_pt_sublead_truth_nocut"," ; Subleading Jet p_{T} [GeV]; counts",pt_N,pt_bins);
     TH1D *h_eta_sublead_truth_nocut = new TH1D("h_eta_sublead_truth_nocut"," ; Subleading Jet #eta  ; counts",48,-1.2,1.2);
     TH1D *h_phi_sublead_truth_nocut = new TH1D("h_phi_sublead_truth_nocut"," ; Subleading Jet #Phi  ; counts",32,-TMath::Pi(),TMath::Pi());

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
    TH1D *h_xj_classical_truth = new TH1D("h_xj_classical_truth",";x_{J};1/N", nbins, ixj_bins);

    TH2D *h_eta_phi_lead_truth_nocut = new TH2D("h_eta_phi_lead_truth_nocut", "Leading Jet #Phi  ;Leading Jet #eta  ",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);TH2D *h_eta_phi_sublead_truth_nocut = new TH2D("h_eta_phi_sublead_truth_nocut", "Subleading Jet #Phi  ;Subleading Jet #eta  ",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);
    TH2D *h_eta_phi_lead_truth_matched = new TH2D("h_eta_phi_lead_truth_matched", "Leading Jet #Phi  ;Leading Jet #eta  ",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);TH2D *h_eta_phi_sublead_truth_matched = new TH2D("h_eta_phi_sublead_truth_matched", "Subleading Jet #Phi  ;Subleading Jet #eta  ",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);


     // Reco No Cut QA histograms
     TH1D *h_pt_lead_reco_nocut = new TH1D("h_pt_lead_reco_nocut"," ; Leading Jet p_{T} [GeV]; counts",nbins,ipt_bins);
     TH1D *h_eta_lead_reco_nocut = new TH1D("h_eta_lead_reco_nocut"," ; Leading Jet #eta  ; counts",48,-1.2,1.2);
     TH1D *h_phi_lead_reco_nocut = new TH1D("h_phi_lead_reco_nocut"," ; Leading Jet #Phi  ; counts",nbins,ipt_bins);

     TH1D *h_pt_sublead_reco_nocut = new TH1D("h_pt_sublead_reco_nocut"," ; Subleading Jet p_{T} [GeV]; counts",nbins,ipt_bins);
     TH1D *h_eta_sublead_reco_nocut = new TH1D("h_eta_sublead_reco_nocut"," ; Subleading Jet #eta  ; counts",48,-1.2,1.2);
     TH1D *h_phi_sublead_reco_nocut = new TH1D("h_phi_sublead_reco_nocut"," ; Subleading Jet #Phi  ; counts",32,-TMath::Pi(),TMath::Pi());

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
    TH1D *h_xj_classical_reco = new TH1D("h_xj_classical_reco",";x_{J};1/N", nbins, ixj_bins);

    TH2D *h_eta_phi_lead_reco_nocut = new TH2D("h_eta_phi_lead_reco_nocut", "Leading Jet #Phi  ;Leading Jet #eta  ",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);TH2D *h_eta_phi_sublead_reco_nocut = new TH2D("h_eta_phi_sublead_reco_nocut", "Subleading Jet #Phi  ;Subleading Jet #eta  ",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);
    TH2D *h_eta_phi_lead_reco_matched = new TH2D("h_eta_phi_lead_reco_matched", "Leading Jet #Phi  ;Leading Jet #eta  ",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);TH2D *h_eta_phi_sublead_reco_matched = new TH2D("h_eta_phi_sublead_reco_matched", "Subleading Jet #Phi  ;Subleading Jet #eta  ",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);
   
   

   //Truth vs Reco Histos No Cut and Matched
    TH2D *h_eta_lead_truth_reco_nocut = new TH2D("h_eta_lead_truth_reco_nocut", "Leading #eta; Subleading#eta", 48,-1.2,1.2, 48,-1.2,1.2);
    TH2D *h_eta_lead_truth_reco_matched = new TH2D("h_eta_lead_truth_reco_matched", "Leading #eta; Subleading#eta", 48,-1.2,1.2, 48,-1.2,1.2);
    TH2D *h_phi_lead_truth_reco_nocut = new TH2D("h_phi_lead_truth_reco_nocut", "Leading #phi; Subleading#phi", 32,-TMath::Pi(),TMath::Pi(),32,-TMath::Pi(),TMath::Pi());
    TH2D *h_phi_lead_truth_reco_matched = new TH2D("h_phi_lead_truth_reco_matched", "Leading #phi; Subleading#phi", 32,-TMath::Pi(),TMath::Pi(),32,-TMath::Pi(),TMath::Pi());


   //Response Matrix
    RooUnfoldResponse *response_pt1pt2_full = new RooUnfoldResponse("response_pt1pt2_full","");
    RooUnfoldResponse *response_pt1pt2_half_fill = new RooUnfoldResponse("response_pt1pt2_half_fill","");
    RooUnfoldResponse *response_pt1pt2_half_test = new RooUnfoldResponse("response_pt1pt2_half_test","");
    response_pt1pt2_full->Setup(h_pt1pt2_reco, h_pt1pt2_truth);
    response_pt1pt2_half_fill->Setup(h_pt1pt2_reco, h_pt1pt2_truth);
    response_pt1pt2_half_test->Setup(h_pt1pt2_reco, h_pt1pt2_truth);



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

	
 // identify file to determine setting
	int setting = 0; float maxpt;float weight;float minpt;
  if (infile1.find("30") != std::string::npos){
    setting = 30; minpt = 30;  maxpt = 200;weight = 2.505*pow(10,-9);
     if(infile1.find("Herwig") != std::string::npos){
       weight = 1.473*pow(10,-12);; cout << "Herwig detected" << endl;
    }
    cout << "30 GEV  JET INFILE" << endl;
   }

  else if (infile1.find("10") != std::string::npos){
    setting = 10; minpt = 10;  maxpt = 20; weight = 3.646*pow(10,-6);
    if(infile1.find("Herwig") != std::string::npos){
       weight = 1.57028*pow(10,-10); maxpt = 30;
      cout << "Herwig detected" << endl;
    }
    cout << "10 GEV JET INFILE" << endl;
   }
     else if (infile1.find("20") != std::string::npos){
    setting = 20;weight = 6.218*pow(10,-8);
    cout << "20 GEV MC  JET INFILE" << endl;
    maxpt = 30; minpt = 20;

   }
  else if (infile1.find("_ana") != std::string::npos){
    setting = 1;maxpt = 200;
    cout << "DATA JET INFILE" << endl;
   }

  std::cout << "minpt is: " << minpt << " max pt is: " << maxpt << endl;
      
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
	      float tdeltaphi = 0; float rdeltaphi = 0; float tzvtx = t_z_vtx;

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
		      float detasublead = fabs(rjeteta - tsubleadeta); float dphisublead = fabs(rjetphi - tsubleadphi); float rzvtx = r_z_vtx;

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

		       bool truth_good = (tleadpt >= truth_leading_cut && tsubleadpt >= truth_subleading_cut && tdeltaphi > dphicut && fabs(tzvtx) < zvtxcut && fabs(tleadeta) < etacut && fabs(tsubleadeta) < etacut && tleadpt >= minpt && tleadpt <= maxpt);

		       bool reco_good = (rleadpt >= reco_leading_cut && rsubleadpt >= reco_subleading_cut && rdeltaphi > dphicut && fabs(rzvtx) < zvtxcut && fabs(rleadeta) < etacut && fabs(rsubleadeta) < etacut);   
		         double  choice = Random.Rndm();

			 float xjtruth = tsubleadpt/tleadpt;

			 if (!truth_good){
			   if (ientry < 10){

			     cout << "Truth Bad in Entry " << ientry << " leadpt is " << tleadpt << " subleadpt is " << tsubleadpt << " deltaphi is " << tdeltaphi << " zvtx is " << tzvtx << " lead eta is " << tleadeta << " subleadeta is " << tsubleadeta << endl;

                           }
			   continue;
			 }
			 //Miss Procedure
			 if(truth_good && !reco_good){
			   if (ientry < 100){
			     cout << "Miss Found: " << ientry << " rleadpt is " << rleadpt << " rsubleadpt is " << rsubleadpt << "choice is: " << choice <<  endl;

			   }
			  

			   response_pt1pt2_full->Miss(tleadpt, tsubleadpt, weight); response_pt1pt2_full->Miss(tsubleadpt, tleadpt, weight);
			   h_pt1pt2_truth->Fill(tleadpt,tsubleadpt,weight);h_pt1pt2_truth->Fill(tsubleadpt,tleadpt,weight);
			   h_xj_classical_truth->Fill(xjtruth,weight);
			   //Half closure test													                    
			   h_pt_lead_truth_matched->Fill(tleadpt,weight);
			   if (choice > 0.5){
			   response_pt1pt2_half_fill->Miss(tleadpt, tsubleadpt, weight); response_pt1pt2_half_fill->Miss(tsubleadpt, tleadpt, weight); 
			
			   }

			   if (choice < 0.5){
			        h_pt1pt2_truth_half->Fill(tleadpt,tsubleadpt,weight); h_pt1pt2_truth_half->Fill(tsubleadpt,tleadpt,weight);

			   }
			 }//end Miss procedure

			 //Fill Rest of Histograms
			 if (truth_good && reco_good){

        

			    xjtruth = tsubleadpt/tleadpt;
			   float xjreco = rsubleadpt/rleadpt;

			   ///Truth pt, eta, phi
			   h_pt_lead_truth_matched->Fill(tleadpt,weight);
			   h_eta_lead_truth_matched->Fill(tleadeta,weight); 
			   h_phi_lead_truth_matched->Fill(tleadphi,weight);

                           h_pt_sublead_truth_matched->Fill(tsubleadpt,weight); 
                           h_eta_sublead_truth_matched->Fill(tsubleadeta,weight);
                           h_phi_sublead_truth_matched->Fill(tsubleadphi,weight);

                        //Truth 2D histos and xj base
			   h_pt1pt2_truth->Fill(tleadpt,tsubleadpt,weight); h_pt1pt2_truth->Fill(tsubleadpt,tleadpt,weight);                         
			   h_xj_classical_truth->Fill(xjtruth,weight);

			       ///Reco pt, eta, phi
			   h_pt_lead_reco_matched->Fill(rleadpt,weight);
			   h_eta_lead_reco_matched->Fill(rleadeta,weight); 
			   h_phi_lead_reco_matched->Fill(rleadphi,weight);

                           h_pt_sublead_reco_matched->Fill(rsubleadpt,weight); 
                           h_eta_sublead_reco_matched->Fill(rsubleadeta,weight);
                           h_phi_sublead_reco_matched->Fill(rsubleadphi,weight); 

			 //Reco 2D histos and xj base
			   h_pt1pt2_reco->Fill(rleadpt,rsubleadpt,weight); h_pt1pt2_reco->Fill(rsubleadpt,rleadpt,weight);  
			   h_xj_classical_reco->Fill(xjreco,weight); 

			      //Full Response Matrix 
			   response_pt1pt2_full->Fill(rleadpt,rsubleadpt,tleadpt,tsubleadpt,weight);response_pt1pt2_full->Fill(rsubleadpt,rleadpt,tsubleadpt,tleadpt,weight);
			      
			     
			      //Half Matrices Fill
			   if (choice > 0.5){
			      h_pt1pt2_reco_half_fill->Fill(rleadpt,rsubleadpt,weight); h_pt1pt2_reco_half_fill->Fill(rsubleadpt,rleadpt,weight);
			    
			      response_pt1pt2_half_fill->Fill(rleadpt,rsubleadpt,tleadpt,tsubleadpt,weight);response_pt1pt2_half_fill->Fill(rsubleadpt,rleadpt,tsubleadpt,tleadpt,weight);
			   } 

			      //Half Matrices Test
			   else if(choice < 0.5){
			     h_pt1pt2_reco_half_test->Fill(rleadpt,rsubleadpt,weight); h_pt1pt2_reco_half_test->Fill(rsubleadpt,rleadpt,weight);
			       response_pt1pt2_half_test->Fill(rleadpt,rsubleadpt,tleadpt,tsubleadpt,weight); response_pt1pt2_half_test->Fill(rsubleadpt,rleadpt,tsubleadpt,tleadpt,weight);

  h_pt1pt2_truth_half->Fill(tleadpt,tsubleadpt,weight); h_pt1pt2_truth_half->Fill(tsubleadpt,tleadpt,weight); 
			   } 
			     

			 }

		

		 }// end last entry loop

		 

	  }// end loop over all truth jets
	}// end can be matched loop
	
	

	
    }// end loop of entries
	
 xj_functions::clean_pt1pt2(h_pt1pt2_reco_half_test,nbins,weight);

    //define file for writing
     TString outfilename = infile1;
     std::string nbinsstring = "_Bins";
      nbinsstring += std::to_string(nbins);
      outfilename.Prepend("Xj_2D_Response_Base"+ nbinsstring + "_R" + Rvaluestring +"_");
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



}//end function
