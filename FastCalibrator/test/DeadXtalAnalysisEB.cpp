#include <vector>
#include <utility>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>

#include "TFile.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TApplication.h"
#include "TEndcapRings.h"

#include "ConfigParser.h"
#include "ntpleUtils.h"

int main (int argc, char** argv){

  gStyle->SetOptFit(111);

  if(argc<2){ std::cout<<" Not correct number of input parameter --> Need Just one cfg file exit "<<std::endl; return -1; }

  parseConfigFile(argv[1]);

  // txt file with the list of input root files                                                                                                                                                 
  std::string inputList = gConfigParser -> readStringOption("Input::inputList");

  std::string DeadChannelMapName ;
  try{ DeadChannelMapName = gConfigParser -> readStringOption("Input::DeadChannelMapName"); }
  catch(char const* exceptionString ){ DeadChannelMapName = "h_map_Dead_Channels";}

  std::string ICMapName ;
  try{ ICMapName = gConfigParser -> readStringOption("Input::ICMapName"); }
  catch(char const* exceptionString ){ ICMapName = "h_scale_EB";}

  int IPhiWindow ;
  try{ IPhiWindow = gConfigParser -> readIntOption("Input::IPhiWindow"); }
  catch(char const* exceptionString ){ IPhiWindow = 3;}

  int IEtaWindow ;
  try{ IEtaWindow = gConfigParser -> readIntOption("Input::IEtaWindow"); }
  catch(char const* exceptionString ){ IEtaWindow = 3;}

  std::string outputCanvasPlot ;
  try{ outputCanvasPlot = gConfigParser -> readStringOption("Output::outputCanvasPlot"); }
  catch(char const* exceptionString ){ outputCanvasPlot = "output/outDeadXtalPlots/";}

  std::cout<<" inputList          : "<<inputList<<std::endl;
  std::cout<<" DeadChannelMapName : "<<DeadChannelMapName<<std::endl;
  std::cout<<" ICMapName          : "<<ICMapName<<std::endl;
  std::cout<<" IPhiWindow         : "<<IPhiWindow<<std::endl;
  std::cout<<" IEtaWindow         : "<<IEtaWindow<<std::endl;
  std::cout<<" outputCanvasPlot   : "<<outputCanvasPlot<<std::endl;

  system(("mkdir -p "+outputCanvasPlot).c_str());
  system(("rm "+outputCanvasPlot+"*").c_str());
  
  std::ifstream inFile(inputList.c_str());
  std::string buffer;

  std::vector<TFile*> inputFileList ; 

  if(!inFile.is_open()){ std::cout << "** ERROR: Can't open '" << inputList << "' for input" << std::endl; return -1;}

  std::cout<<" Input Files : "<<std::endl;  

  while(!inFile.eof()){

    inFile >> buffer;
    if( buffer.at(0) == '#' ) continue;
    std::cout<<buffer<<std::endl;    
    inputFileList.push_back(new TFile(buffer.c_str(),"READ"));
  }

  std::vector<std::vector<TH1F*> > ICCrystalEB(IPhiWindow, std::vector<TH1F*> (IEtaWindow) );

  std::vector<TH2F*> DeadCrystalEB(inputFileList.size());
  std::vector<TH2F*> ICMapEB(inputFileList.size());  

  for( int iPhi = 0 ; iPhi < IPhiWindow ; iPhi ++){
    for( int iEta = 0 ; iEta < IEtaWindow ; iEta ++) {
     (ICCrystalEB.at(iPhi)).at(iEta) = new TH1F(std::string(Form("ICCrystalEB_%d_%d",iPhi,iEta)).c_str(),"",50,0.95,1.4);
    }  
 }
     
  for( unsigned int iFile = 0 ; iFile < inputFileList.size() ; iFile ++){

    DeadCrystalEB.push_back((TH2F*) inputFileList.at(iFile)->Get(DeadChannelMapName.c_str()));
    ICMapEB.push_back((TH2F*) inputFileList.at(iFile)->Get(ICMapName.c_str()));

    for( int iPhi = 0 ; iPhi < 360 ; iPhi ++){
      for( int iEta = 0 ; iEta < 170 ; iEta ++){
 
	if(DeadCrystalEB.back()->GetBinContent(iPhi+1,iEta+1)!=0 && iEta!=0){
          ICCrystalEB[int(IPhiWindow/2)][int(IEtaWindow/2)]->Fill( ICMapEB.back()->GetBinContent(iPhi+1,iEta+1));

	  for( int IPHI = iPhi - int((IPhiWindow-1)/2) ; IPHI <= iPhi + int((IPhiWindow-1)/2) ; IPHI ++){
	    for( int IETA = iEta - int((IEtaWindow-1)/2) ; IETA <= iEta + int((IEtaWindow-1)/2) ; IETA ++){
	      (ICCrystalEB.at(int(IPhiWindow/2)+(IPHI-iPhi))).at(int((IEtaWindow)/2)+(IETA-iEta))->Fill(ICMapEB.back()->GetBinContent(IPHI+1,IETA+1));
          }
	 }	
       }
      }
    }
  }

  std::vector<std::vector<TCanvas*> >Can(IPhiWindow);
  for( int iPhi = 0 ; iPhi < IPhiWindow ; iPhi ++){
    for(int iEta = 0 ; iEta < IPhiWindow ; iEta ++)
      (Can.at(iPhi)).push_back( new TCanvas(std::string(Form("Can_%d_%d",iPhi,iEta)).c_str(),"",500,500));
  }

  std::vector<std::vector<TF1*> > GaussianFits (IPhiWindow, std::vector<TF1*> (IEtaWindow));

  
  for( int iPhi = 0 ; iPhi < IPhiWindow ; iPhi ++){
   for( int iEta = 0 ; iEta < IEtaWindow ; iEta ++){
     GaussianFits[iPhi][iEta] = new TF1(std::string(Form("Gaus_%d_%d",iPhi,iEta)).c_str(),"gaus",0.9,1.5);
     (Can.at(iPhi)).at(iEta)->cd();
     (Can.at(iPhi)).at(iEta)->SetGridx();
     (Can.at(iPhi)).at(iEta)->SetGridy();
     ICCrystalEB[iPhi][iEta]->SetLineColor(kBlack);
     ICCrystalEB[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB[iPhi][iEta]->SetMarkerStyle(20);
     ICCrystalEB[iPhi][iEta]->SetMarkerSize(1.);
     ICCrystalEB[iPhi][iEta]->SetMarkerColor(kRed);
     GaussianFits[iPhi][iEta]->SetLineColor(kBlue);
     GaussianFits[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB[iPhi][iEta]->Fit(GaussianFits[iPhi][iEta],"RMEQ");
     ICCrystalEB[iPhi][iEta]->Draw("E");
     (Can.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can.at(iPhi)).at(iEta)->GetName())+".pdf").c_str(),"pdf");
     (Can.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can.at(iPhi)).at(iEta)->GetName())+".png").c_str(),"png");
   }
  }
  return 0; 

}
