void Make_xJ_Projection_Test(string infile1 = "TREE_DIJET_v6_1_ana462_2024p010_v001_gl10-00047352-00047733.root"){

#define _USE_MATH_DEFINES

#include <math.h> 
#include <cmath>
#include <iostream>

#include<vector>
#include<array>



  int nbins = 19;
  float minimum =  5.5;
  float  fixed = 14;
  int beforebin = 7;

  float ixj_bins[nbins+1];float ipt_bins[nbins+1];

    float alpha = TMath::Power(fixed/minimum, 1/(float)beforebin);
    float maximum = minimum*TMath::Power(alpha, nbins );
    float bin_xj10 = 1.0;
    float bin_xj1 = 1.0*(minimum/maximum);


     for (int i = 0; i < nbins+1; i++)
      {
	float ipt = minimum*TMath::Power(alpha, (float)i);
	ipt_bins[i] = ipt;
	  }


    for (int i = 0; i < nbins+1; i++)
      {
	float ixj = bin_xj1*TMath::Power(alpha, (float)i);
	ixj_bins[i] = ixj;
      }

    for (int i = 0 ; i < nbins + 1; i++)
    {
      std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
    }


   const double_t pt_bins[]={0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,80};
  const int pt_N = sizeof(pt_bins)/sizeof(pt_bins[0]) - 1;
   
  //Define Data Histos


    //2D Pt and xJ Histos
    TH2D *h_pt1pt2 =  new TH2D("h_pt1pt2",";p_{T,1};p_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
    TH1D *h_xj_classical = new TH1D("h_xj_classical",";x_{J};1/N", nbins, ixj_bins);
    TH1D *h_xj_projected = new TH1D("h_xj_projected",";x_{J};1/N", nbins, ixj_bins);

    //1D Pt histos

    TH1D* hDataPtLead_30vtx = new TH1D("hDataPtLead_30vtx","",nbins,ipt_bins); 
    TH1D* hDataPtLead_60vtx = new TH1D("hDataPtLead_60vtx","",nbins,ipt_bins); 
    TH1D* hDataPtSubLead_30vtx = new TH1D("hDataPtSubLead_30vtx","",nbins,ipt_bins); 
    TH1D* hDataPtSubLead_60vtx = new TH1D("hDataPtSubLead_60vtx","",nbins,ipt_bins);  
    TH1D* hDataPtInclusive_30vtx = new TH1D("hDataPtInclusive_30vtx","",nbins,ipt_bins); 
    TH1D* hDataPtInclusive_60vtx = new TH1D("hDataPtInclusive_60vtx","",nbins,ipt_bins);

    //eta 1D
    TH1D* hDataEtaLead_30vtx = new TH1D("hDataEtaLead_30vtx","",48,-1.2,1.2); 
    TH1D* hDataEtaLead_60vtx = new TH1D("hDataEtaLead_60vtx","",48,-1.2,1.2); 
    TH1D* hDataEtaSubLead_30vtx = new TH1D("hDataEtaSubLead_30vtx","",48,-1.2,1.2); 
    TH1D* hDataEtaSubLead_60vtx = new TH1D("hDataEtaSubLead_60vtx","",48,-1.2,1.2);  
    TH1D* hDataEtaInclusive_30vtx = new TH1D("hDataEtaInclusive_30vtx","",48,-1.2,1.2); 
    TH1D* hDataEtaInclusive_60vtx = new TH1D("hDataEtaInclusive_60vtx","",48,-1.2,1.2);

    //phi 1D
    TH1D* hDataPhiLead_30vtx = new TH1D("hDataPhiLead_30vtx","",32,-TMath::Pi(),TMath::Pi()); 
    TH1D* hDataPhiLead_60vtx = new TH1D("hDataPhiLead_60vtx","",32,-TMath::Pi(),TMath::Pi()); 
    TH1D* hDataPhiSubLead_30vtx = new TH1D("hDataPhiSubLead_30vtx","",32,-TMath::Pi(),TMath::Pi()); 
    TH1D* hDataPhiSubLead_60vtx = new TH1D("hDataPhiSubLead_60vtx","",32,-TMath::Pi(),TMath::Pi());  
    TH1D* hDataPhiInclusive_30vtx = new TH1D("hDataPhiInclusive_30vtx","",32,-TMath::Pi(),TMath::Pi()); 
    TH1D* hDataPhiInclusive_60vtx = new TH1D("hDataPhiInclusive_60vtx","",32,-TMath::Pi(),TMath::Pi());

    //Eta Phi 2D
    TH2D* hDataEtaPhiLead_30vtx = new TH2D("hDataEtaPhiLead_30vtx","",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2); 
    TH2D* hDataEtaPhiLead_60vtx = new TH2D("hDataEtaPhiLead_60vtx","",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2); 
    TH2D* hDataEtaPhiSubLead_30vtx = new TH2D("hDataEtaPhiSubLead_30vtx","",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2); 
    TH2D* hDataEtaPhiSubLead_60vtx = new TH2D("hDataEtaPhiSubLead_60vtx","",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);  
    TH2D* hDataEtaPhiInclusive_30vtx = new TH2D("hDataEtaPhiInclusive_30vtx","",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2); 
    TH2D* hDataEtaPhiInclusive_60vtx = new TH2D("hDataEtaPhiInclusive_60vtx","",32,-TMath::Pi(),TMath::Pi(),48,-1.2,1.2);

    //dphi and deta

    TH1D* hDataDeltaPhi_30vtx = new TH1D("hDataDeltaPhi_30vtx","",50,TMath::Pi()/2,2*TMath::Pi()); 
    TH1D* hDataDeltaPhi_60vtx = new TH1D("hDataDeltaPhi_60vtx","",50,TMath::Pi()/2,2*TMath::Pi());
    TH1D* hDataDeltaEta_30vtx = new TH1D("hDataDeltaEta_30vtx","",25,0,2.4); 
    TH1D* hDataDeltaEta_60vtx = new TH1D("hDataDeltaEta_60vtx","",25,0,2.4);

   
    // identify file to determine setting
  int setting = 0; 
  if (infile1.find("_30_") != std::string::npos){
    setting = 30;
    cout << "30 GEV PYTHIA JET INFILE" << endl;
   }
  else if (infile1.find("_10_") != std::string::npos){
    setting = 10;
    cout << "10 GEV PYTHIA JET INFILE" << endl;
   }
  else if (infile1.find("_ana") != std::string::npos){
    setting = 1;
    cout << "DATA JET INFILE" << endl;
   }

    //connect variables in the code to variables in the tree
   TFile* _file0 = TFile::Open(infile1.c_str());
   TTree* tree = (TTree*)_file0->Get("ttree");

   //define vectors and variables

    std::vector<float> *rjet_pt ={0}; std::vector<float> *rjet_phi ={0}; std::vector<float> *rjet_eta ={0};   std::vector<float> *rjet_et ={0}; std::vector<float> *rjet_e ={0}; float z_vtx = 0;
    unsigned long long triggers;
    tree->SetBranchAddress("jet_pt_4",&rjet_pt); tree->SetBranchAddress("jet_eta_4",&rjet_eta); tree->SetBranchAddress("jet_phi_4",&rjet_phi);tree->SetBranchAddress("jet_et_4",&rjet_et);tree->SetBranchAddress("jet_e_4",&rjet_e); tree->SetBranchAddress("mbd_vertex_z",&z_vtx);tree->SetBranchAddress("gl1_scaled",&triggers);


   //define variables

    float rleadphi = 0; float rsubleadphi = 0; 
   float rleadeta = 0; float rsubleadeta = 0; float rleadE = 0; float rleadET = 0; 
   float rdeltaphi = 0; float rdeltaeta = 0;

   float Rvalue = 0.4; float etacut = 1.1 - Rvalue; float entry = 0;

   	float nrjets = 0;
   //QA Variables

	float leadptpass = 0; float leadetapass = 0; float subleadetapass = 0; float dphipass = 0; float subleadptpass = 0; float endgame = 0; float sixtyfilled = 0; float thirtyfilled = 0;float lastevent = 0; float smallsize = 0; float zvtxpass = 0; 

	float inclusivesixty = 0; float inclusivethirty = 0;

   int Nentries = tree->GetEntries(); 
   //int Nentries = 100000;
    cout << "Nentries = : " << Nentries << endl;
    //start loop over entries
    for (int i = 0;i<Nentries;i++){

        //loop over entries
        tree->GetEntry(i);
	entry++; 
	bool flag = false;  

	// Check if bit 17 is set
	if (triggers & (1LL << 17)) {
	  flag = true; 
	}


	float rleadpt = 0; float rsubleadpt = 0;
	  nrjets = rjet_pt->size();
	  /*
	  if (nrjets < 2){
	   smallsize++;
	   continue;
	   }*/

	  //  else if (nrjets >= 2){ //open conditional loop

	    for (int j = 0; j < nrjets;j++){ //open nrjets loop

	       float  rjetpt = rjet_pt->at(j); float rjetphi =  rjet_phi->at(j);  float rjeteta = rjet_eta->at(j);
	       float zvtx = z_vtx;	     
	       float weight = 0;


	       //inclusive jets

	       if (rjetpt > 10 && fabs(rjeteta) < etacut && fabs(zvtx) < 60 && flag){

		 weight = 1; inclusivesixty++;
		 hDataPtInclusive_60vtx->Fill(rjetpt,weight); 
		 hDataEtaInclusive_60vtx->Fill(rjeteta,weight);
		 hDataPhiInclusive_60vtx->Fill(rjetphi,weight);
		 hDataEtaPhiInclusive_60vtx->Fill(rjetphi,rjeteta);
		      
		 if (fabs(zvtx) < 30){
		   inclusivethirty++;
		 hDataPtInclusive_30vtx->Fill(rjetpt,weight); 
		 hDataEtaInclusive_30vtx->Fill(rjeteta,weight);
		 hDataPhiInclusive_30vtx->Fill(rjetphi,weight);
		 hDataEtaPhiInclusive_30vtx->Fill(rjetphi,rjeteta);
		 }

	       }
	     
	        if (rjetpt > rleadpt){
			rsubleadpt = rleadpt;	rleadpt = rjetpt;
			rsubleadphi = rleadphi;	rleadphi = rjetphi;
	        	rsubleadeta = rleadeta; rleadeta = rjeteta;
			
		      }
		  
		 else if (rjetpt > rsubleadpt){
			rsubleadpt = rjetpt;	rsubleadphi = rjetphi;	rsubleadeta = rjeteta;
		      }
		 if (j == (nrjets - 1)){//last entry of event

		     rdeltaphi = fabs(rleadphi-rsubleadphi);rdeltaeta = fabs(rleadeta-rsubleadeta);
		        if ( rdeltaphi > TMath::Pi())
			{
			  rdeltaphi = 2*TMath::Pi()- rdeltaphi;
			}
	
		     lastevent++;

		     if(rleadpt > 15){
		       leadptpass++;
		     }
		     if (rsubleadpt > 8){
		       subleadptpass++;
		     }
		     if (rdeltaphi > 3*TMath::Pi()/4){
		       dphipass++;
		     }
		     if (fabs(rleadeta) < etacut ){
		       leadetapass++;

		     }
		       if (fabs(rsubleadeta) < etacut ){
		       subleadetapass++;

		     }

		       if (fabs(zvtx) < 60){
			 zvtxpass++;
		       }
		     
		     //reco cuts loop
		       if(rleadpt > 18.2836  && rsubleadpt > 8.2087 && rdeltaphi > 3*TMath::Pi()/4 && fabs(rleadeta) < etacut && fabs(rsubleadeta) < etacut && fabs(zvtx) < 60 && flag){
			endgame++;
		       
		        //Fill DATA Histograms
			 weight = 1;
			  //zvtx 60 Histos	   
			    sixtyfilled++;
			    //Data pt
			   hDataPtLead_60vtx->Fill(rleadpt,weight); hDataPtSubLead_60vtx->Fill(rsubleadpt,weight);
			     // Data eta
			   hDataEtaLead_60vtx->Fill(rleadeta,weight);hDataEtaSubLead_60vtx->Fill(rsubleadeta,weight);	       
			     //Data phi
			   hDataPhiLead_60vtx->Fill(rleadphi,weight); hDataPhiSubLead_60vtx->Fill(rsubleadphi,weight); 	
			    //Data eta phi
			    hDataEtaPhiLead_60vtx->Fill(rleadphi,rleadeta);hDataEtaPhiSubLead_60vtx->Fill(rsubleadphi,rsubleadeta);
			     //Data dphi and deta
			   hDataDeltaPhi_60vtx->Fill(rdeltaphi,weight); hDataDeltaEta_60vtx->Fill(rdeltaeta,weight);
			   
			   float xj = rsubleadpt/rleadpt;
			   h_pt1pt2->Fill(rleadpt,rsubleadpt);h_pt1pt2->Fill(rsubleadpt,rleadpt);
			   h_xj_classical->Fill(xj);

			  
			 //zvtx 30
			 if(fabs(zvtx) < 30){ 
			    thirtyfilled++;
			    //Data pt
			   hDataPtLead_30vtx->Fill(rleadpt,weight); hDataPtSubLead_30vtx->Fill(rsubleadpt,weight);
			   // Data eta
			   hDataEtaLead_30vtx->Fill(rleadeta,weight);hDataEtaSubLead_30vtx->Fill(rsubleadeta,weight);
			     //Data phi
			   hDataPhiLead_30vtx->Fill(rleadphi,weight); hDataPhiSubLead_30vtx->Fill(rsubleadphi,weight); 
			     //Data eta phi
			    hDataEtaPhiLead_30vtx->Fill(rleadphi,rleadeta);hDataEtaPhiSubLead_30vtx->Fill(rsubleadphi,rsubleadeta);
			     //Data dphi and deta
			   hDataDeltaPhi_30vtx->Fill(rdeltaphi,weight); hDataDeltaEta_30vtx->Fill(rdeltaeta,weight);
		

			 }//end zvtx 30


		      }//end reco cuts loop

		      // rleadpt = 0; rsubleadpt = 0;
		 }//close last entry loop

	    }//close nrjets loops
	    //  }// close conditional loop

    }// end loop over entries

 //xj projection making loop

			    TH1D *h_unc = (TH1D*) h_xj_projected->Clone();
			    h_unc->Reset();

			    TH2D *h_asym_pt1pt2 = (TH2D*) h_pt1pt2->Clone();

			      for (int ix = 0; ix < nbins; ix++)
				{
				  for (int iy = 0; iy < nbins; iy++)
				    {
				      int bin = h_asym_pt1pt2->GetBin(ix+1, iy+1);

				      if (ix > iy)
					{
					  h_asym_pt1pt2->SetBinContent(bin, h_asym_pt1pt2->GetBinContent(bin)*2.);
					  h_asym_pt1pt2->SetBinError(bin, h_asym_pt1pt2->GetBinError(bin)*2);
					}
				      else if (ix < iy)
					{
					  h_asym_pt1pt2->SetBinContent(bin, 0);
					  h_asym_pt1pt2->SetBinError(bin, 0);
					}
	
	      
				    }
				}

			        for (int ix = 0; ix < nbins; ix++)
				  {
				    for (int iy = 0; iy <= ix; iy++)
				      {
					int low =  iy - ix - 1;
					int high = iy - ix + 1;

					int xjbin_low = nbins + low + 1;
					int xjbin_high = nbins + low + 2;
					//std::cout << ix << " / " << iy << " -- > " << xjbin_low << "--"<<xjbin_high<<std::endl;
					int bin = h_asym_pt1pt2->GetBin(ix+1, iy+1);

			       

					if (ix == iy)
					  {

					    h_xj_projected->Fill(h_xj_projected->GetBinCenter(xjbin_low), h_asym_pt1pt2->GetBinContent(bin));
					    h_unc->Fill(h_xj_projected->GetBinCenter(xjbin_low),TMath::Power( h_asym_pt1pt2->GetBinError(bin), 2));
					  }
					else
					  {
					    h_xj_projected->Fill(h_xj_projected->GetBinCenter(xjbin_high), h_asym_pt1pt2->GetBinContent(bin)/2.);
					    h_xj_projected->Fill(h_xj_projected->GetBinCenter(xjbin_low), h_asym_pt1pt2->GetBinContent(bin)/2.);
					    h_unc->Fill(h_xj_projected->GetBinCenter(xjbin_high),TMath::Power( h_asym_pt1pt2->GetBinError(bin)/sqrt(2), 2));
					    h_unc->Fill(h_xj_projected->GetBinCenter(xjbin_low),TMath::Power( h_asym_pt1pt2->GetBinError(bin)/sqrt(2), 2));
					  }
				      }
				  }

				for (int ix = 0; ix < nbins; ix++)
				  {
				    h_xj_projected->SetBinError(ix+1, sqrt(h_unc->GetBinContent(ix+1)));
				  }  for (int ix = 0; ix < nbins; ix++)
				       {
					 for (int iy = 0; iy <= ix; iy++)
					   {
					     int low =  iy - ix - 1;
					     int high = iy - ix + 1;

					     int xjbin_low = nbins + low + 1;
					     int xjbin_high = nbins + low + 2;
					     //std::cout << ix << " / " << iy << " -- > " << xjbin_low << "--"<<xjbin_high<<std::endl;
					     int bin = h_asym_pt1pt2->GetBin(ix+1, iy+1);

					  
					     if (ix == iy)
					       {

						 h_xj_projected->Fill(h_xj_projected->GetBinCenter(xjbin_low), h_asym_pt1pt2->GetBinContent(bin));
						 h_unc->Fill(h_xj_projected->GetBinCenter(xjbin_low),TMath::Power( h_asym_pt1pt2->GetBinError(bin), 2));
					       }
					     else
					       {
						 h_xj_projected->Fill(h_xj_projected->GetBinCenter(xjbin_high), h_asym_pt1pt2->GetBinContent(bin)/2.);
						 h_xj_projected->Fill(h_xj_projected->GetBinCenter(xjbin_low), h_asym_pt1pt2->GetBinContent(bin)/2.);
						 h_unc->Fill(h_xj_projected->GetBinCenter(xjbin_high),TMath::Power( h_asym_pt1pt2->GetBinError(bin)/sqrt(2), 2));
						 h_unc->Fill(h_xj_projected->GetBinCenter(xjbin_low),TMath::Power( h_asym_pt1pt2->GetBinError(bin)/sqrt(2), 2));
					       }
					   }
				       }

				for (int ix = 0; ix < nbins; ix++)
				  {
				    h_xj_projected->SetBinError(ix+1, sqrt(h_unc->GetBinContent(ix+1)));
				  }



    cout << " lead pt passed: " << leadptpass << " lead eta pass: " << leadetapass << " sublead pt passed: " << subleadptpass << " sublead eta pass: " << subleadetapass << " zvtx pass: " << zvtxpass << endl;
 cout << " dijet cut passed: " << dphipass << " endgame entered: " << endgame << " 60 vtx filled: " << sixtyfilled << " 30 vtx filled: " << thirtyfilled << " last event reached: " << lastevent << " jet too small: " << smallsize << "total entries gone through: " << lastevent + smallsize <<  endl;


    //Save hists to file for later
      TString outfilename = infile1;
     outfilename.Prepend("Hists_xJ_ProjectionTest_");
     TFile *outfile = TFile::Open(outfilename,"RECREATE");

     


        //Write DATA Histos
     h_pt1pt2->Write(); h_xj_classical->Write(); h_xj_projected->Write();

      hDataPtLead_60vtx->Write(); hDataPtSubLead_60vtx->Write();
      hDataPtInclusive_60vtx->Write(); hDataPtInclusive_60vtx->Write();
      hDataEtaLead_60vtx->Write();hDataEtaSubLead_60vtx->Write();
      hDataEtaInclusive_60vtx->Write();hDataEtaInclusive_60vtx->Write();
      hDataPhiLead_60vtx->Write(); hDataPhiSubLead_60vtx->Write(); 
      hDataPhiInclusive_60vtx->Write();hDataPhiInclusive_60vtx->Write();
      hDataDeltaPhi_60vtx->Write(); hDataDeltaEta_60vtx->Write();
       hDataEtaPhiLead_60vtx->Write();hDataEtaPhiSubLead_60vtx->Write();
      hDataEtaPhiInclusive_60vtx->Write();


      hDataPtLead_30vtx->Write(); hDataPtSubLead_30vtx->Write();
      hDataPtInclusive_30vtx->Write(); hDataPtInclusive_30vtx->Write();
      hDataEtaLead_30vtx->Write();hDataEtaSubLead_30vtx->Write();
      hDataEtaInclusive_30vtx->Write();hDataEtaInclusive_30vtx->Write();
      hDataPhiLead_30vtx->Write(); hDataPhiSubLead_30vtx->Write(); 
      hDataPhiInclusive_30vtx->Write();hDataPhiInclusive_30vtx->Write();
      hDataDeltaPhi_30vtx->Write(); hDataDeltaEta_30vtx->Write();
      hDataEtaPhiLead_30vtx->Write();hDataEtaPhiSubLead_30vtx->Write();
      hDataEtaPhiInclusive_30vtx->Write();



       

}


  
