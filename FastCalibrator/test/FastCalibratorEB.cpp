#include "FastCalibratorEB.h"
#include <iostream>
#include <fstream>
#include <stdlib.h> 

#include "ConfigParser.h"
#include "ntpleUtils.h"
#include "CalibrationUtils.h"


int main (int argc, char ** argv) {
    
  ///Check if all nedeed arguments to parse are there
  if(argc != 2){
    std::cerr << ">>>>> FastCalibrator::usage: " << argv[0] << " configFileName" << std::endl ;
    return 1;
  }
    
  /// Parse the config file
  parseConfigFile (argv[1]) ;
 
  // txt file with the list of input root files
  std::string inputList = gConfigParser -> readStringOption("Input::inputList");

  // input tree name 
  std::string inputTree = "NULL";
  try{ inputTree = gConfigParser -> readStringOption("Input::inputTree");}
  catch(char const* exceptionString ){ inputTree = "simpleNtupleEoverP/SimpleNtupleEoverP";}

  // input dead xtal name --> switch off by hand 
  std::string inputFileDeadXtal = "NULL";
  try{ inputFileDeadXtal = gConfigParser -> readStringOption("Input::inputFileDeadXtal");}
  catch( char const* exceptionString ){ inputFileDeadXtal = "NULL"; }

  // input dead xtal name --> switch off by hand 
  bool isDeadTriggerTower =false ;
  try{ if(gConfigParser -> readIntOption("Input::isDeadTriggerTower") == 1) isDeadTriggerTower = true;}
  catch( char const* exceptionString ){ isDeadTriggerTower = false; }
  std::cout<<" isDeadTriggerTower "<<isDeadTriggerTower<<std::endl;  
  // jsonFileName  
  std::string jsonFileName ="NULL";
  try{  jsonFileName  = gConfigParser -> readStringOption("Input::jsonFileName");}
  catch( char const* exceptionString ){ jsonFileName = "json/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt";}

  std::map<int, std::vector<std::pair<int, int> > > jsonMap;
  jsonMap = readJSONFile(jsonFileName);

  // map with scalibration
  std::string miscalibMap = "NULL";
  try{ miscalibMap = gConfigParser -> readStringOption("Input::miscalibMap");}
  catch(char const* exceptionString ){ miscalibMap = "/gwteray/users/brianza/scalibMap2.txt";}
  
  // Miscalibration 
  bool isMiscalib ;
  try{isMiscalib = gConfigParser -> readBoolOption("Input::isMiscalib");}
  catch( char const* exceptionString ){ isMiscalib = false;}

  // Save EoverP distribution
  bool isSaveEPDistribution ;
  try{ isSaveEPDistribution = gConfigParser -> readBoolOption("Input::isSaveEPDistribution");}
  catch( char const* exceptionString ){ isSaveEPDistribution = false; }

  // Do the E/P selection
  bool isEPselection ;
  try{ isEPselection = gConfigParser -> readBoolOption("Input::isEPselection");}
  catch( char const* exceptionString ){ isEPselection = false; }
  
  // Pt treshold bool and cut
  bool isPtCut ;
  try{ isPtCut = gConfigParser -> readBoolOption("Input::isPtCut"); }
  catch( char const* exceptionString ){ isPtCut = false;}

  float PtMin ;
  try{ PtMin = gConfigParser -> readFloatOption("Input::PtMin");}
  catch( char const* exceptionString ){ PtMin = 0.;}
 
  // fbrem treshold bool and cut
  bool isfbrem ;
  try { isfbrem = gConfigParser -> readBoolOption("Input::isfbrem"); }
  catch( char const* exceptionString ){ isfbrem = false;}
 
  float fbremMax ;
  try { fbremMax = gConfigParser -> readFloatOption("Input::fbremMax"); }
  catch( char const* exceptionString ){ fbremMax = 100.;}

  // R9 treshold bool and cut
  bool isR9selection ;
  try{ isR9selection = gConfigParser -> readBoolOption("Input::isR9selection");}
  catch( char const* exceptionString ){ isR9selection = false; }

  float R9Min ;
  try{ R9Min = gConfigParser -> readFloatOption("Input::R9Min");}
  catch( char const* exceptionString ){ R9Min = 0.; }

  // Run Calibration on E/Etrue instead of E/P --> MC only
  bool isMCTruth ;
  try { isMCTruth = gConfigParser -> readBoolOption("Input::isMCTruth"); }
  catch( char const* exceptionString ){ isMCTruth = false; }

  // method for the miscalibration
  float miscalibMethod ;
  try { miscalibMethod = gConfigParser -> readFloatOption("Input::miscalibMethod"); }
  catch( char const* exceptionString ){ miscalibMethod = false; }

  // Momentum scale file
  std::string inputMomentumScale ;
  try{ inputMomentumScale =  gConfigParser -> readStringOption("Input::inputMomentumScale"); }
  catch( char const* exceptionString ) { inputMomentumScale = "output/MomentumCalibrationCombined_2011AB-2012ABC.root";}

  std::string typeEB ;
  try{ typeEB = gConfigParser -> readStringOption("Input::typeEB"); }
  catch( char const* exceptionString ) { typeEB = "none" ; }

  std::string typeEE ;
  try{ typeEE = gConfigParser -> readStringOption("Input::typeEE"); }
  catch( char const* exceptionString ) { typeEE = "none" ; }
  
  int nRegionsEB = GetNRegionsEB(typeEB);
  
  // Name of the output calib file
  std::string outputPath;
  try{ outputPath = gConfigParser -> readStringOption("Output::outputPath");
       system(("mkdir -p "+outputPath).c_str());}
  catch( char const* exceptionString ) { outputPath = "output/Oct22_Run2012ABC_Cal_Dic2012/";
                                         system(("mkdir -p "+outputPath).c_str());
  }
  
  std::string outputFile ;
  try{ outputFile = gConfigParser -> readStringOption("Output::outputFile"); }
  catch( char const* exceptionString ) { outputFile = "FastCalibrator_Oct22_Run2012ABC_Cal_Dic2012" ;
  }
  
  // Other options for the L3 algo
  int numberOfEvents ;
  try { numberOfEvents = gConfigParser -> readIntOption("Options::numberOfEvents"); }
  catch( char const* exceptionString ) { numberOfEvents = -1 ; }

  int useZ ;
  try { useZ = gConfigParser -> readIntOption("Options::useZ"); }
  catch( char const* exceptionString ) { useZ = 1 ; }

  int useW ;
  try{ useW = gConfigParser -> readIntOption("Options::useW"); }
  catch( char const* exceptionString ) { useW = 1 ; }

  int splitStat ;
  try{ splitStat = gConfigParser -> readIntOption("Options::splitStat"); }
  catch( char const* exceptionString ) { splitStat = 0 ; }

  int nLoops ;
  try{ nLoops = gConfigParser -> readIntOption("Options::nLoops"); }
  catch( char const* exceptionString ) { nLoops = 20 ; }
  
  
  /// open ntupla of data or MC
  TChain * tree = new TChain (inputTree.c_str());
  FillChain(*tree,inputList); 
  
  /// open calibration momentum graph
  TFile* momentumscale = new TFile((inputMomentumScale+"_"+typeEB+"_"+typeEE+".root").c_str());
  std::vector<TGraphErrors*> g_EoC_EB;
  
  for(int i = 0; i < nRegionsEB; ++i){
    TString Name = Form("g_EoC_EB_%d",i);
    g_EoC_EB.push_back( (TGraphErrors*)(momentumscale->Get(Name)) );
  }
  
  ///Use the whole sample statistics if numberOfEvents < 0
  if ( numberOfEvents < 0 ) numberOfEvents = tree->GetEntries(); 
  
  /// run in normal mode: full statistics
  if ( splitStat == 0 ) {
   
    TString name ;
    TString name_tmp;

    if(isMiscalib == true && useZ == 1 && isR9selection ==true && isEPselection == false && isfbrem == false && isPtCut ==false ) 
     name_tmp = Form ("%s_WZ_R9_miscalib_EB",outputFile.c_str());
    else if(isMiscalib == true && useZ == 1 && isR9selection ==false && isEPselection == false && isfbrem == true && isPtCut ==false ) 
     name_tmp = Form ("%s_WZ_Fbrem_miscalib_EB",outputFile.c_str());
    else if(isMiscalib == true && useZ == 1 && isR9selection ==true && isEPselection == false && isfbrem == false && isPtCut ==true ) 
     name_tmp = Form ("%s_WZ_PT_miscalib_EB",outputFile.c_str());
    else if(isMiscalib == true && useZ == 1 && isR9selection ==true && isEPselection == true && isfbrem == false && isPtCut ==false ) 
     name_tmp = Form ("%s_WZ_EP_miscalib_EB",outputFile.c_str());
    else if(isMiscalib == true && useZ == 1 && isEPselection ==false && isR9selection==false && isPtCut ==false && isfbrem ==false  ) 
     name_tmp =Form ("%s_WZ_noEP_miscalib_EB",outputFile.c_str());

    else if(isMiscalib == false && useZ == 1 && isR9selection ==true && isEPselection == false && isfbrem == false && isPtCut ==false ) 
     name_tmp = Form ("%s_WZ_R9_EB",outputFile.c_str());
    else if(isMiscalib == false && useZ == 1 && isR9selection ==false && isEPselection == false && isfbrem == true && isPtCut ==false ) 
     name_tmp = Form ("%s_WZ_Fbrem_EB",outputFile.c_str());
    else if(isMiscalib == false && useZ == 1 && isR9selection ==true && isEPselection == false && isfbrem == false && isPtCut ==true ) 
     name_tmp = Form ("%s_WZ_PT_EB",outputFile.c_str());
    else if(isMiscalib == false && useZ == 1 && isR9selection ==true && isEPselection == true && isfbrem == false && isPtCut ==false ) 
     name_tmp = Form ("%s_WZ_EP_EB",outputFile.c_str());
    else if(isMiscalib == false && useZ == 1 && isEPselection ==false && isR9selection==false && isPtCut ==false && isfbrem ==false  ) 
     name_tmp =Form ("%s_WZ_noEP_EB",outputFile.c_str());

    else if(isMiscalib == true && useZ == 0 && isR9selection ==true && isEPselection == false && isfbrem == false && isPtCut ==false ) 
     name_tmp = Form ("%s_W_R9_miscalib_EB",outputFile.c_str());
    else if(isMiscalib == true && useZ == 0 && isR9selection ==false && isEPselection == false && isfbrem == true && isPtCut ==false ) 
     name_tmp = Form ("%s_W_Fbrem_miscalib_EB",outputFile.c_str());
    else if(isMiscalib == true && useZ == 0 && isR9selection ==true && isEPselection == false && isfbrem == false && isPtCut ==true ) 
     name_tmp = Form ("%s_W_PT_miscalib_EB",outputFile.c_str());
    else if(isMiscalib == true && useZ == 0 && isR9selection ==true && isEPselection == true && isfbrem == false && isPtCut ==false ) 
     name_tmp = Form ("%s_W_EP_miscalib_EB",outputFile.c_str());
    else if(isMiscalib == true && useZ == 0 && isEPselection ==false && isR9selection==false && isPtCut ==false && isfbrem ==false  ) 
     name_tmp =Form ("%s_W_noEP_miscalib_EB",outputFile.c_str());

    else if(isMiscalib == false && useZ == 0 && isR9selection ==true && isEPselection == false && isfbrem == false && isPtCut ==false ) 
     name_tmp = Form ("%s_WZ_R9_EB",outputFile.c_str());
    else if(isMiscalib == false && useZ == 0 && isR9selection ==false && isEPselection == false && isfbrem == true && isPtCut ==false ) 
     name_tmp = Form ("%s_W_Fbrem_EB",outputFile.c_str());
    else if(isMiscalib == false && useZ == 0 && isR9selection ==true && isEPselection == false && isfbrem == false && isPtCut ==true ) 
     name_tmp = Form ("%s_W_PT_EB",outputFile.c_str());
    else if(isMiscalib == false && useZ == 0 && isR9selection ==true && isEPselection == true && isfbrem == false && isPtCut ==false ) 
     name_tmp = Form ("%s_W_EP_EB",outputFile.c_str());
    else if(isMiscalib == false && useZ == 0 && isEPselection ==false && isR9selection==false && isPtCut ==false && isfbrem ==false  ) 
     name_tmp =Form ("%s_W_noEP_EB",outputFile.c_str());
    else { std::cout<<" Option not considered --> exit "<<std::endl; return -1 ;}

    name = Form("%s%s.root",outputPath.c_str(),name_tmp.Data());
    TFile *outputName = new TFile(name,"RECREATE");
    
    TString outEPDistribution = "Weight_"+name;
    
    TString DeadXtal = Form("%s",inputFileDeadXtal.c_str());    

    if(isSaveEPDistribution == true){
      FastCalibratorEB analyzer(tree, g_EoC_EB, typeEB, outEPDistribution);
      analyzer.bookHistos(nLoops);
      analyzer.AcquireDeadXtal(DeadXtal,isDeadTriggerTower);
      analyzer.Loop(numberOfEvents, useZ, useW, splitStat, nLoops, isMiscalib,isSaveEPDistribution,isEPselection,isR9selection,R9Min,isfbrem,fbremMax,isPtCut,PtMin,isMCTruth,jsonMap, miscalibMethod, miscalibMap);
      analyzer.saveHistos(outputName);
    }
    else{
      FastCalibratorEB analyzer(tree, g_EoC_EB, typeEB);
      analyzer.bookHistos(nLoops);
      analyzer.AcquireDeadXtal(DeadXtal,isDeadTriggerTower);
      analyzer.Loop(numberOfEvents, useZ, useW, splitStat, nLoops, isMiscalib,isSaveEPDistribution,isEPselection,isR9selection,R9Min,isfbrem,fbremMax,isPtCut,PtMin,isMCTruth,jsonMap, miscalibMethod, miscalibMap);
      analyzer.saveHistos(outputName);
    }
    
  }
  
  /// run in even-odd mode: half statistics
  else if ( splitStat == 1 ) {
    
    /// Prepare the outputs
    TString name;
    TString name2;
    
    if(isMiscalib == true && useZ == 1 && isR9selection==true && isEPselection==false && isfbrem ==false && isPtCut==false){
      name  = Form ("%s_WZ_R9_miscalib_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_WZ_R9_miscalib_EB_odd.root", outputFile.c_str());
    }
    else if(isMiscalib == true && useZ == 1 && isR9selection==false && isEPselection==true && isfbrem ==false && isPtCut==false){
      name  = Form ("%s_WZ_EP_miscalib_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_WZ_EP_miscalib_EB_odd.root", outputFile.c_str());
    }
    else if(isMiscalib == true && useZ == 1 && isR9selection==true && isEPselection==false && isfbrem ==true && isPtCut==false){ 
      name  = Form ("%s_WZ_Fbrem_miscalib_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_WZ_Fbrem_miscalib_EB_odd.root", outputFile.c_str());
    }
    else if(isMiscalib == true && useZ == 1 && isR9selection==true && isEPselection==false && isfbrem ==false && isPtCut==true){
      name  = Form ("%s_WZ_PT_miscalib_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_WZ_PT_miscalib_EB_odd.root", outputFile.c_str());
    }
    else if(isMiscalib == true && useZ == 1 && isR9selection==false && isEPselection==false && isfbrem ==false && isPtCut==false){
      name  = Form ("%s_WZ_noEP_miscalib_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_WZ_noEP_miscalib_EB_odd.root", outputFile.c_str());
    }



    else if(isMiscalib == false && useZ == 1 && isR9selection==true && isEPselection==false && isfbrem ==false && isPtCut==false){
      name  = Form ("%s_WZ_R9_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_WZ_R9_EB_odd.root", outputFile.c_str());
    }
    else if(isMiscalib == false && useZ == 1 && isR9selection==false && isEPselection==true && isfbrem ==false && isPtCut==false){
      name  = Form ("%s_WZ_EP_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_WZ_EP_EB_odd.root", outputFile.c_str());
    }
    else if(isMiscalib == false && useZ == 1 && isR9selection==false && isEPselection==false && isfbrem ==true && isPtCut==false){
      name  = Form ("%s_WZ_Fbrem_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_WZ_Fbrem_EB_odd.root", outputFile.c_str());
    }
    else if(isMiscalib == false && useZ == 1 && isR9selection==false && isEPselection==false && isfbrem ==false && isPtCut==true){
      name  = Form ("%s_WZ_PT_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_WZ_PT_EB_odd.root", outputFile.c_str());
    }
    else if(isMiscalib == false && useZ == 1 && isR9selection==false && isEPselection==false && isfbrem ==false && isPtCut==false){
      name  = Form ("%s_WZ_noEP_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_WZ_noEP_EB_odd.root", outputFile.c_str());
    }
    

    else if(isMiscalib == true && useZ == 0 && isR9selection==true && isEPselection==false && isfbrem ==false && isPtCut==false){
      name  = Form ("%s_W_R9_miscalib_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_W_R9_miscalib_EB_odd.root", outputFile.c_str());
    }
    
    else if(isMiscalib == true && useZ == 0 && isR9selection==false && isEPselection==true && isfbrem ==false && isPtCut==false){
      name  = Form ("%s_W_EP_miscalib_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_W_EP_miscalib_EB_odd.root", outputFile.c_str());
    }
    else if(isMiscalib == true && useZ == 0 && isR9selection==false && isEPselection==false && isfbrem ==true && isPtCut==false){
      name  = Form ("%s_W_Fbrem_miscalib_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_W_Fbrem_miscalib_EB_odd.root", outputFile.c_str());
    }
    else if(isMiscalib == true && useZ == 0 && isR9selection==false && isEPselection==false && isfbrem ==false && isPtCut==true){
      name  = Form ("%s_W_PT_miscalib_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_W_PT_miscalib_EB_odd.root", outputFile.c_str());
    }
    else if(isMiscalib == true && useZ == 0 && isR9selection==false && isEPselection==false && isfbrem ==false && isPtCut==false){
      name  = Form ("%s_W_noEP_miscalib_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_W_noEP_miscalib_EB_odd.root", outputFile.c_str());
    }

    else if(isMiscalib == false && useZ == 0 && isR9selection==true && isEPselection==false && isfbrem ==false && isPtCut==false){
      name  = Form ("%s_W_R9_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_W_R9_EB_odd.root", outputFile.c_str());
    }

    
    else if(isMiscalib == false && useZ == 0 && isR9selection==true && isEPselection==false && isfbrem ==false && isPtCut==false){
      name  = Form ("%s_W_EP_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_W_EP_EB_odd.root", outputFile.c_str());
    }
    else if(isMiscalib == false && useZ == 0 && isR9selection==false && isEPselection==true && isfbrem ==false && isPtCut==false){
       name  = Form ("%s_W_EP_EB_even.root",outputFile.c_str());
       name2 = Form ("%s_W_EP_EB_odd.root", outputFile.c_str());
    }
    else if(isMiscalib == false && useZ == 0 && isR9selection==false && isEPselection==false && isfbrem ==true && isPtCut==false){
      name  = Form ("%s_W_Fbrem_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_W_Fbrem_EB_odd.root", outputFile.c_str());
    }
    else if(isMiscalib == false && useZ == 0 && isR9selection==false && isEPselection==false && isfbrem ==false && isPtCut==true){
      name  = Form ("%s_W_PT_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_W_PT_EB_odd.root", outputFile.c_str());
    }
    else if(isMiscalib == false && useZ == 0 && isR9selection==true && isEPselection==false && isfbrem ==false && isPtCut==false){
      name  = Form ("%s_W_noEP_EB_even.root",outputFile.c_str());
      name2 = Form ("%s_W_noEP_EB_odd.root", outputFile.c_str());
    }
    else { std::cout<<" Option not considered --> exit "<<std::endl; return -1 ;}

    TFile *outputName1 = new TFile(outputPath+name,"RECREATE");
    TFile *outputName2 = new TFile(outputPath+name2,"RECREATE");

    TString DeadXtal = Form("%s",inputFileDeadXtal.c_str());
     
    /// Run on odd
    FastCalibratorEB analyzer_even(tree, g_EoC_EB, typeEB);
    analyzer_even.bookHistos(nLoops);
    analyzer_even.AcquireDeadXtal(DeadXtal,isDeadTriggerTower);
    analyzer_even.Loop(numberOfEvents, useZ, useW, splitStat, nLoops,isMiscalib,isSaveEPDistribution,isEPselection,isR9selection,R9Min,isfbrem,fbremMax,isPtCut,PtMin,isMCTruth,jsonMap, miscalibMethod, miscalibMap);
    analyzer_even.saveHistos(outputName1);
  
    /// Run on even
    FastCalibratorEB analyzer_odd(tree, g_EoC_EB, typeEB);
    analyzer_odd.bookHistos(nLoops);
    analyzer_odd.AcquireDeadXtal(DeadXtal,isDeadTriggerTower);
    analyzer_odd.Loop(numberOfEvents, useZ, useW, splitStat*(-1), nLoops,isMiscalib,isSaveEPDistribution,isEPselection,isR9selection,R9Min,isfbrem,fbremMax,isPtCut,PtMin,isMCTruth,jsonMap, miscalibMethod, miscalibMap);
    analyzer_odd.saveHistos(outputName2);
    
  }
  
  delete tree;
  return 0;
}
