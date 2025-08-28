#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>

#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllUtils.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <g4centrality/PHG4CentralityReco.h>

#include <HIJetReco.C>
#include <JetValidationData.h>
#include <jetbase/JetReco.h>

#include <Calo_Calib.C>

#include <mbd/MbdReco.h>
#include <globalvertex/GlobalVertexReco.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libJetValidation.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libmbd.so)


 
#endif
void Fun4All_JetVal_Data(const char *filelistcalofit = "calofittest.list", const char *outname = "outputestdata2segmentcalofittingattempt2.root")
{

 
  
  Fun4AllServer *se = Fun4AllServer::instance();
  int verbosity = 0;

  se->Verbosity(verbosity);
  recoConsts *rc = recoConsts::instance();
  
  /*const std::string &fname = filelistcalo;
  
  pair<int, int> runseg = Fun4AllUtils::GetRunSegment(fname);
  int runnumber = runseg.first;*/

   int runnumber;
  std::ifstream file(filelistcalofit);
  std::string line;
  std::getline(file, line);
  file.close();
  std::regex pattern(R"((\d+)-\d+\.root)");
  std::smatch matches;    
  if (std::regex_search(line, matches, pattern) && matches.size() > 1) {
    runnumber = std::stoi(matches[1].str());
    std::cout << "Run number: " << runnumber << std::endl;
  } else {
    std::cerr << "Run number is not clear." << std::endl;
    return;
  }
  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc-> set_uint64Flag("TIMESTAMP",runnumber);
  CDBInterface::instance()-> Verbosity(1);
  //rc -> set_uint64Flag("TIMESTAMP",runnumber);

    gSystem->Load("libg4dst");

// MBD/BBC Reconstruction
  MbdReco *mbdreco = new MbdReco();
  se->registerSubsystem(mbdreco);

  // Official vertex storage
  GlobalVertexReco *gvertex = new GlobalVertexReco();
  se->registerSubsystem(gvertex);

  

    // Geometry 
  std::cout << "Adding Geometry file" << std::endl;
  Fun4AllInputManager *intrue2 = new Fun4AllRunNodeInputManager("DST_GEO");
  std::string geoLocation = CDBInterface::instance()->getUrl("calo_geo");
  intrue2->AddFile(geoLocation);
  se->registerInputManager(intrue2);

  /*
  Fun4AllInputManager *in0 = new Fun4AllDstInputManager("DSTJet");
  in0->AddListFile(filelistjet,1);
  se->registerInputManager(in0);*/

  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTCaloFitting");
  in1->AddListFile(filelistcalofit,1);
  se->registerInputManager(in1);

  Process_Calo_Calib();
  HIJetReco();



  PHG4CentralityReco *cent = new PHG4CentralityReco();
  cent->Verbosity(verbosity);
  cent->GetCalibrationParameters().ReadFromFile("centrality", "xml", 0, 0, string(getenv("CALIBRATIONROOT")) + string("/Centrality/"));
  se->registerSubsystem(cent);
  
  Enable::VERBOSITY = verbosity;
 

   
  JetValidation *myJetVal = new JetValidation("AntiKt_Tower_r04", "AntiKt_Truth_r04", outname);

  myJetVal->setPtRange(5, 100);
  myJetVal->setEtaRange(-1.1, 1.1);
  myJetVal->doUnsub(0);
  myJetVal->doTruth(0);
  myJetVal->doSeeds(0);
  se->registerSubsystem(myJetVal);


  /*
  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("G4Hits");
  in1->AddFile(filelistG4);
  se->registerInputManager(in1);

   Fun4AllInputManager *intrue = new Fun4AllDstInputManager("DSTtruth");
  intrue->AddFile(filelisttruth);
  se->registerInputManager(intrue);

  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DSTcalo");
  in2->AddFile(filelistcalo);
  se->registerInputManager(in2);

   Fun4AllInputManager *in3 = new Fun4AllDstInputManager("DSTglobal");
  in3->AddFile(filelistglobal);
  se->registerInputManager(in3);*/

  
  
  se->run(-1);
  se->End();

  gSystem->Exit(0);
  //  return 0;

}
