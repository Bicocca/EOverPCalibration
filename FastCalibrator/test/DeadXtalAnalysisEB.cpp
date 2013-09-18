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
  gStyle->SetOptStat(1100);

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

  std::vector<std::vector<TH1F*> > ICCrystalEB_Mod1(IPhiWindow, std::vector<TH1F*> (IEtaWindow) );
  std::vector<std::vector<TH1F*> > ICCrystalEB_Mod2(IPhiWindow, std::vector<TH1F*> (IEtaWindow) );
  std::vector<std::vector<TH1F*> > ICCrystalEB_Mod3(IPhiWindow, std::vector<TH1F*> (IEtaWindow) );
  std::vector<std::vector<TH1F*> > ICCrystalEB_Mod4(IPhiWindow, std::vector<TH1F*> (IEtaWindow) );

  std::vector<std::vector<TH1F*> > ICCrystalEB_Phi1(IPhiWindow, std::vector<TH1F*> (IEtaWindow) );
  std::vector<std::vector<TH1F*> > ICCrystalEB_Phi2(IPhiWindow, std::vector<TH1F*> (IEtaWindow) );
  std::vector<std::vector<TH1F*> > ICCrystalEB_Phi3(IPhiWindow, std::vector<TH1F*> (IEtaWindow) );
  std::vector<std::vector<TH1F*> > ICCrystalEB_Phi4(IPhiWindow, std::vector<TH1F*> (IEtaWindow) );

  std::vector<TH2F*> DeadCrystalEB(inputFileList.size());
  std::vector<TH2F*> ICMapEB(inputFileList.size());  

  for( int iPhi = 0 ; iPhi < IPhiWindow ; iPhi ++){
    for( int iEta = 0 ; iEta < IEtaWindow ; iEta ++) {
     (ICCrystalEB.at(iPhi)).at(iEta)      = new TH1F(std::string(Form("ICCrystalEB_%d_%d",iPhi,iEta)).c_str(),"",50,0.95,1.4);
     (ICCrystalEB_Mod1.at(iPhi)).at(iEta) = new TH1F(std::string(Form("ICCrystalEB_Mod1_%d_%d",iPhi,iEta)).c_str(),"",50,0.95,1.4);
     (ICCrystalEB_Mod2.at(iPhi)).at(iEta) = new TH1F(std::string(Form("ICCrystalEB_Mod2_%d_%d",iPhi,iEta)).c_str(),"",50,0.95,1.4);
     (ICCrystalEB_Mod3.at(iPhi)).at(iEta) = new TH1F(std::string(Form("ICCrystalEB_Mod3_%d_%d",iPhi,iEta)).c_str(),"",50,0.95,1.4);
     (ICCrystalEB_Mod4.at(iPhi)).at(iEta) = new TH1F(std::string(Form("ICCrystalEB_Mod4_%d_%d",iPhi,iEta)).c_str(),"",50,0.95,1.4);
     (ICCrystalEB_Phi1.at(iPhi)).at(iEta) = new TH1F(std::string(Form("ICCrystalEB_Phi1_%d_%d",iPhi,iEta)).c_str(),"",50,0.95,1.4);
     (ICCrystalEB_Phi2.at(iPhi)).at(iEta) = new TH1F(std::string(Form("ICCrystalEB_Phi2_%d_%d",iPhi,iEta)).c_str(),"",50,0.95,1.4);
     (ICCrystalEB_Phi3.at(iPhi)).at(iEta) = new TH1F(std::string(Form("ICCrystalEB_Phi3_%d_%d",iPhi,iEta)).c_str(),"",50,0.95,1.4);
     (ICCrystalEB_Phi4.at(iPhi)).at(iEta) = new TH1F(std::string(Form("ICCrystalEB_Phi4_%d_%d",iPhi,iEta)).c_str(),"",50,0.95,1.4);
    }  
 }
     
  for( unsigned int iFile = 0 ; iFile < inputFileList.size() ; iFile ++){

    DeadCrystalEB.push_back((TH2F*) inputFileList.at(iFile)->Get(DeadChannelMapName.c_str()));
    ICMapEB.push_back((TH2F*) inputFileList.at(iFile)->Get(ICMapName.c_str()));

    for( int iPhi = 0 ; iPhi < 360 ; iPhi ++){
      for( int iEta = 0 ; iEta < 170 ; iEta ++){
 
	if(DeadCrystalEB.back()->GetBinContent(iPhi+1,iEta+1)!=0 && iEta!=0){
          ICCrystalEB[int(IPhiWindow/2)][int(IEtaWindow/2)]->Fill( ICMapEB.back()->GetBinContent(iPhi+1,iEta+1));

          if(fabs(iEta-85)<=20)                          ICCrystalEB_Mod1[int(IPhiWindow/2)][int(IEtaWindow/2)]->Fill( ICMapEB.back()->GetBinContent(iPhi+1,iEta+1));
          else if(fabs(iEta-85)>20 && fabs(iEta-85)<=40) ICCrystalEB_Mod2[int(IPhiWindow/2)][int(IEtaWindow/2)]->Fill( ICMapEB.back()->GetBinContent(iPhi+1,iEta+1));
          else if(fabs(iEta-85)>40 && fabs(iEta-85)<=60) ICCrystalEB_Mod3[int(IPhiWindow/2)][int(IEtaWindow/2)]->Fill( ICMapEB.back()->GetBinContent(iPhi+1,iEta+1));
          else                                           ICCrystalEB_Mod4[int(IPhiWindow/2)][int(IEtaWindow/2)]->Fill( ICMapEB.back()->GetBinContent(iPhi+1,iEta+1));

          if( int(iPhi%20)<=5 )                          ICCrystalEB_Phi1[int(IPhiWindow/2)][int(IEtaWindow/2)]->Fill( ICMapEB.back()->GetBinContent(iPhi+1,iEta+1));
          else if( int(iPhi%20)>5  && int(iPhi%20)<=10 ) ICCrystalEB_Phi2[int(IPhiWindow/2)][int(IEtaWindow/2)]->Fill( ICMapEB.back()->GetBinContent(iPhi+1,iEta+1)); 
	  else if( int(iPhi%20)>10 && int(iPhi%20)<=15 ) ICCrystalEB_Phi3[int(IPhiWindow/2)][int(IEtaWindow/2)]->Fill( ICMapEB.back()->GetBinContent(iPhi+1,iEta+1)); 
	  else if( int(iPhi%20)>15 && int(iPhi%20)<=20 ) ICCrystalEB_Phi4[int(IPhiWindow/2)][int(IEtaWindow/2)]->Fill( ICMapEB.back()->GetBinContent(iPhi+1,iEta+1));
          else continue ;

	  for( int IPHI = iPhi - int((IPhiWindow-1)/2) ; IPHI <= iPhi + int((IPhiWindow-1)/2) ; IPHI ++){
	    for( int IETA = iEta - int((IEtaWindow-1)/2) ; IETA <= iEta + int((IEtaWindow-1)/2) ; IETA ++){
	      (ICCrystalEB.at(int(IPhiWindow/2)+(IPHI-iPhi))).at(int((IEtaWindow)/2)+(IETA-iEta))->Fill(ICMapEB.back()->GetBinContent(IPHI+1,IETA+1));

              if(fabs(iEta-85)<=20)                          
                (ICCrystalEB_Mod1.at(int(IPhiWindow/2)+(IPHI-iPhi))).at(int((IEtaWindow)/2)+(IETA-iEta))->Fill(ICMapEB.back()->GetBinContent(IPHI+1,IETA+1));
              else if(fabs(iEta-85)>20 && fabs(iEta-85)<=40) 
                (ICCrystalEB_Mod2.at(int(IPhiWindow/2)+(IPHI-iPhi))).at(int((IEtaWindow)/2)+(IETA-iEta))->Fill(ICMapEB.back()->GetBinContent(IPHI+1,IETA+1));
              else if(fabs(iEta-85)>40 && fabs(iEta-85)<=60) 
                (ICCrystalEB_Mod3.at(int(IPhiWindow/2)+(IPHI-iPhi))).at(int((IEtaWindow)/2)+(IETA-iEta))->Fill(ICMapEB.back()->GetBinContent(IPHI+1,IETA+1));
              else                                           
                (ICCrystalEB_Mod4.at(int(IPhiWindow/2)+(IPHI-iPhi))).at(int((IEtaWindow)/2)+(IETA-iEta))->Fill(ICMapEB.back()->GetBinContent(IPHI+1,IETA+1));


              if( int(iPhi%20)<=5 )  
		(ICCrystalEB_Phi1.at(int(IPhiWindow/2)+(IPHI-iPhi))).at(int((IEtaWindow)/2)+(IETA-iEta))->Fill(ICMapEB.back()->GetBinContent(IPHI+1,IETA+1));                                   
              else if( int(iPhi%20)>5  && int(iPhi%20)<=10) 
                (ICCrystalEB_Phi2.at(int(IPhiWindow/2)+(IPHI-iPhi))).at(int((IEtaWindow)/2)+(IETA-iEta))->Fill(ICMapEB.back()->GetBinContent(IPHI+1,IETA+1));
              else if( int(iPhi%20)>10 && int(iPhi%20)<=15) 
                (ICCrystalEB_Phi3.at(int(IPhiWindow/2)+(IPHI-iPhi))).at(int((IEtaWindow)/2)+(IETA-iEta))->Fill(ICMapEB.back()->GetBinContent(IPHI+1,IETA+1));
              else if( int(iPhi%20)>15 && int(iPhi%20)<=20) 
                (ICCrystalEB_Phi4.at(int(IPhiWindow/2)+(IPHI-iPhi))).at(int((IEtaWindow)/2)+(IETA-iEta))->Fill(ICMapEB.back()->GetBinContent(IPHI+1,IETA+1));
              else continue ;

          }
	 }	
       }
      }
    }
  }

  std::vector<std::vector<TCanvas*> >Can(IPhiWindow, std::vector<TCanvas*> (IEtaWindow));
  std::vector<std::vector<TCanvas*> >Can_Mod1(IPhiWindow,std::vector<TCanvas*> (IEtaWindow));
  std::vector<std::vector<TCanvas*> >Can_Mod2(IPhiWindow,std::vector<TCanvas*> (IEtaWindow));
  std::vector<std::vector<TCanvas*> >Can_Mod3(IPhiWindow,std::vector<TCanvas*> (IEtaWindow));
  std::vector<std::vector<TCanvas*> >Can_Mod4(IPhiWindow,std::vector<TCanvas*> (IEtaWindow));
  std::vector<std::vector<TCanvas*> >Can_Phi1(IPhiWindow,std::vector<TCanvas*> (IEtaWindow));
  std::vector<std::vector<TCanvas*> >Can_Phi2(IPhiWindow,std::vector<TCanvas*> (IEtaWindow));
  std::vector<std::vector<TCanvas*> >Can_Phi3(IPhiWindow,std::vector<TCanvas*> (IEtaWindow));
  std::vector<std::vector<TCanvas*> >Can_Phi4(IPhiWindow,std::vector<TCanvas*> (IEtaWindow));

  for( int iPhi = 0 ; iPhi < IPhiWindow ; iPhi ++){
    for(int iEta = 0 ; iEta < IPhiWindow ; iEta ++){
      (Can.at(iPhi)).at(iEta) = new TCanvas(std::string(Form("Can_%d_%d",iPhi,iEta)).c_str(),"",500,500);
      (Can_Mod1.at(iPhi)).at(iEta) = new TCanvas(std::string(Form("Can_Mod1_%d_%d",iPhi,iEta)).c_str(),"",500,500);
      (Can_Mod2.at(iPhi)).at(iEta) = new TCanvas(std::string(Form("Can_Mod2_%d_%d",iPhi,iEta)).c_str(),"",500,500);
      (Can_Mod3.at(iPhi)).at(iEta) = new TCanvas(std::string(Form("Can_Mod3_%d_%d",iPhi,iEta)).c_str(),"",500,500);
      (Can_Mod4.at(iPhi)).at(iEta) = new TCanvas(std::string(Form("Can_Mod4_%d_%d",iPhi,iEta)).c_str(),"",500,500);
      (Can_Phi1.at(iPhi)).at(iEta) = new TCanvas(std::string(Form("Can_Phi1_%d_%d",iPhi,iEta)).c_str(),"",500,500);
      (Can_Phi2.at(iPhi)).at(iEta) = new TCanvas(std::string(Form("Can_Phi2_%d_%d",iPhi,iEta)).c_str(),"",500,500);
      (Can_Phi3.at(iPhi)).at(iEta) = new TCanvas(std::string(Form("Can_Phi3_%d_%d",iPhi,iEta)).c_str(),"",500,500);
      (Can_Phi4.at(iPhi)).at(iEta) = new TCanvas(std::string(Form("Can_Phi4_%d_%d",iPhi,iEta)).c_str(),"",500,500);
    }
  }
  std::vector<std::vector<TF1*> > GaussianFits (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > GaussianFits_Mod1 (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > GaussianFits_Mod2 (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > GaussianFits_Mod3 (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > GaussianFits_Mod4 (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > GaussianFits_Phi1 (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > GaussianFits_Phi2 (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > GaussianFits_Phi3 (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > GaussianFits_Phi4 (IPhiWindow, std::vector<TF1*> (IEtaWindow));

  std::vector<std::vector<TH1F*> > MeanIC_vs_Eta (IPhiWindow, std::vector<TH1F*> (IEtaWindow));
  std::vector<std::vector<TH1F*> > Gaus_MeanIC_vs_Eta (IPhiWindow, std::vector<TH1F*> (IEtaWindow));
  std::vector<std::vector<TH1F*> > MeanIC_vs_Phi (IPhiWindow, std::vector<TH1F*> (IEtaWindow));
  std::vector<std::vector<TH1F*> > Gaus_MeanIC_vs_Phi (IPhiWindow, std::vector<TH1F*> (IEtaWindow));

  std::vector<std::vector<TF1*> > Pol0_MeanIC_vs_Eta (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > Pol1_MeanIC_vs_Eta (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > Pol2_MeanIC_vs_Eta (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > Pol0_Gaus_MeanIC_vs_Eta (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > Pol1_Gaus_MeanIC_vs_Eta (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > Pol2_Gaus_MeanIC_vs_Eta (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > Pol0_MeanIC_vs_Phi (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > Pol1_MeanIC_vs_Phi (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > Pol2_MeanIC_vs_Phi (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > Pol0_Gaus_MeanIC_vs_Phi (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > Pol1_Gaus_MeanIC_vs_Phi (IPhiWindow, std::vector<TF1*> (IEtaWindow));
  std::vector<std::vector<TF1*> > Pol2_Gaus_MeanIC_vs_Phi (IPhiWindow, std::vector<TF1*> (IEtaWindow));

  std::vector<std::vector<TCanvas*> > Can_MeanIC_vs_Eta(IPhiWindow,std::vector<TCanvas*> (IEtaWindow));
  std::vector<std::vector<TCanvas*> > Can_Gaus_MeanIC_vs_Eta(IPhiWindow,std::vector<TCanvas*> (IEtaWindow));
  std::vector<std::vector<TCanvas*> > Can_MeanIC_vs_Phi(IPhiWindow,std::vector<TCanvas*> (IEtaWindow));
  std::vector<std::vector<TCanvas*> > Can_Gaus_MeanIC_vs_Phi(IPhiWindow,std::vector<TCanvas*> (IEtaWindow));

  TH1F* htemp1; 
  TH1F* htemp2; 
  TPaveStats* pave1, *pave2, *pave3;
  
  for( int iPhi = 0 ; iPhi < IPhiWindow ; iPhi ++){
   for( int iEta = 0 ; iEta < IEtaWindow ; iEta ++){

     if(iPhi==1 && iEta ==1) continue ;

     GaussianFits[iPhi][iEta] = new TF1(std::string(Form("Gaus_%d_%d",iPhi,iEta)).c_str(),"gaus",0.9,1.5);
     GaussianFits_Mod1[iPhi][iEta] = new TF1(std::string(Form("Gaus_Mod1_%d_%d",iPhi,iEta)).c_str(),"gaus",0.9,1.5);
     GaussianFits_Mod2[iPhi][iEta] = new TF1(std::string(Form("Gaus_Mod2_%d_%d",iPhi,iEta)).c_str(),"gaus",0.9,1.5);
     GaussianFits_Mod3[iPhi][iEta] = new TF1(std::string(Form("Gaus_Mod3_%d_%d",iPhi,iEta)).c_str(),"gaus",0.9,1.5);
     GaussianFits_Mod4[iPhi][iEta] = new TF1(std::string(Form("Gaus_Mod4_%d_%d",iPhi,iEta)).c_str(),"gaus",0.9,1.5);

     GaussianFits_Phi1[iPhi][iEta] = new TF1(std::string(Form("Gaus_Phi1_%d_%d",iPhi,iEta)).c_str(),"gaus",0.9,1.5);
     GaussianFits_Phi2[iPhi][iEta] = new TF1(std::string(Form("Gaus_Phi2_%d_%d",iPhi,iEta)).c_str(),"gaus",0.9,1.5);
     GaussianFits_Phi3[iPhi][iEta] = new TF1(std::string(Form("Gaus_Phi3_%d_%d",iPhi,iEta)).c_str(),"gaus",0.9,1.5);
     GaussianFits_Phi4[iPhi][iEta] = new TF1(std::string(Form("Gaus_Phi4_%d_%d",iPhi,iEta)).c_str(),"gaus",0.9,1.5);

     MeanIC_vs_Eta[iPhi][iEta]     = new TH1F(std::string(Form("MeanIC_vs_Eta_%d_%d",iPhi,iEta)).c_str(),"",4,0,4); 
     Can_MeanIC_vs_Eta[iPhi][iEta] = new TCanvas(std::string(Form("Can_MeanIC_vs_Eta_%d_%d",iPhi,iEta)).c_str(),"",500,500);

     MeanIC_vs_Phi[iPhi][iEta]     = new TH1F(std::string(Form("MeanIC_vs_Phi_%d_%d",iPhi,iEta)).c_str(),"",4,0,4); 
     Can_MeanIC_vs_Phi[iPhi][iEta] = new TCanvas(std::string(Form("Can_MeanIC_vs_Phi_%d_%d",iPhi,iEta)).c_str(),"",500,500);

     Gaus_MeanIC_vs_Eta[iPhi][iEta]     = new TH1F(std::string(Form("Gaus_MeanIC_vs_Eta_%d_%d",iPhi,iEta)).c_str(),"",4,0,4); 
     Can_Gaus_MeanIC_vs_Eta[iPhi][iEta] = new TCanvas(std::string(Form("Can_Gaus_MeanIC_vs_Eta_%d_%d",iPhi,iEta)).c_str(),"",500,500);

     Gaus_MeanIC_vs_Phi[iPhi][iEta]     = new TH1F(std::string(Form("Gaus_MeanIC_vs_Phi_%d_%d",iPhi,iEta)).c_str(),"",4,0,4); 
     Can_Gaus_MeanIC_vs_Phi[iPhi][iEta] = new TCanvas(std::string(Form("Can_Gaus_MeanIC_vs_Phi_%d_%d",iPhi,iEta)).c_str(),"",500,500);

     Pol0_MeanIC_vs_Eta[iPhi][iEta]  = new TF1(std::string(Form("Pol0_MeanIC_Eta_%d_%d",iPhi,iEta)).c_str(),"pol0",0,4);
     Pol1_MeanIC_vs_Eta[iPhi][iEta]  = new TF1(std::string(Form("Pol1_MeanIC_Eta_%d_%d",iPhi,iEta)).c_str(),"pol1",0,4);
     Pol2_MeanIC_vs_Eta[iPhi][iEta]  = new TF1(std::string(Form("Pol2_MeanIC_Eta_%d_%d",iPhi,iEta)).c_str(),"pol2",0,4);
     Pol0_Gaus_MeanIC_vs_Eta[iPhi][iEta]  = new TF1(std::string(Form("Pol0_Gaus_MeanIC_Eta_%d_%d",iPhi,iEta)).c_str(),"pol0",0,4);
     Pol1_Gaus_MeanIC_vs_Eta[iPhi][iEta]  = new TF1(std::string(Form("Pol1_Gaus_MeanIC_Eta_%d_%d",iPhi,iEta)).c_str(),"pol1",0,4);
     Pol2_Gaus_MeanIC_vs_Eta[iPhi][iEta]  = new TF1(std::string(Form("Pol2_Gaus_MeanIC_Eta_%d_%d",iPhi,iEta)).c_str(),"pol2",0,4);

     Pol0_MeanIC_vs_Phi[iPhi][iEta]  = new TF1(std::string(Form("Pol0_MeanIC_Phi_%d_%d",iPhi,iEta)).c_str(),"pol0",0,4);
     Pol1_MeanIC_vs_Phi[iPhi][iEta]  = new TF1(std::string(Form("Pol1_MeanIC_Phi_%d_%d",iPhi,iEta)).c_str(),"pol1",0,4);
     Pol2_MeanIC_vs_Phi[iPhi][iEta]  = new TF1(std::string(Form("Pol2_MeanIC_Phi_%d_%d",iPhi,iEta)).c_str(),"pol2",0,4);
     Pol0_Gaus_MeanIC_vs_Phi[iPhi][iEta]  = new TF1(std::string(Form("Pol0_Gaus_MeanIC_Phi_%d_%d",iPhi,iEta)).c_str(),"pol0",0,4);
     Pol1_Gaus_MeanIC_vs_Phi[iPhi][iEta]  = new TF1(std::string(Form("Pol1_Gaus_MeanIC_Phi_%d_%d",iPhi,iEta)).c_str(),"pol1",0,4);
     Pol2_Gaus_MeanIC_vs_Phi[iPhi][iEta]  = new TF1(std::string(Form("Pol2_Gaus_MeanIC_Phi_%d_%d",iPhi,iEta)).c_str(),"pol2",0,4);


     /// Original IC distribution

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
     ICCrystalEB[iPhi][iEta]->GetXaxis()->SetTitle("IC");
     ICCrystalEB[iPhi][iEta]->GetYaxis()->SetTitle("Entries");
     (Can.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can.at(iPhi)).at(iEta)->GetName())+".pdf").c_str(),"pdf");
     (Can.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can.at(iPhi)).at(iEta)->GetName())+".png").c_str(),"png");

     /// Original IC distribution Mod 1

     (Can_Mod1.at(iPhi)).at(iEta)->cd();
     (Can_Mod1.at(iPhi)).at(iEta)->SetGridx();
     (Can_Mod1.at(iPhi)).at(iEta)->SetGridy();
     ICCrystalEB_Mod1[iPhi][iEta]->SetLineColor(kBlack);
     ICCrystalEB_Mod1[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Mod1[iPhi][iEta]->SetMarkerStyle(20);
     ICCrystalEB_Mod1[iPhi][iEta]->SetMarkerSize(1.);
     ICCrystalEB_Mod1[iPhi][iEta]->SetMarkerColor(kRed);
     GaussianFits_Mod1[iPhi][iEta]->SetLineColor(kBlue);
     GaussianFits_Mod1[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Mod1[iPhi][iEta]->Fit(GaussianFits_Mod1[iPhi][iEta],"RMEQ");
     ICCrystalEB_Mod1[iPhi][iEta]->Draw("E");
     ICCrystalEB_Mod1[iPhi][iEta]->GetXaxis()->SetTitle("IC");
     ICCrystalEB_Mod1[iPhi][iEta]->GetYaxis()->SetTitle("Entries");
     (Can_Mod1.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Mod1.at(iPhi)).at(iEta)->GetName())+".pdf").c_str(),"pdf");
     (Can_Mod1.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Mod1.at(iPhi)).at(iEta)->GetName())+".png").c_str(),"png");

     /// Original IC distribution Mod 2

     (Can_Mod2.at(iPhi)).at(iEta)->cd();
     (Can_Mod2.at(iPhi)).at(iEta)->SetGridx();
     (Can_Mod2.at(iPhi)).at(iEta)->SetGridy();
     ICCrystalEB_Mod2[iPhi][iEta]->SetLineColor(kBlack);
     ICCrystalEB_Mod2[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Mod2[iPhi][iEta]->SetMarkerStyle(20);
     ICCrystalEB_Mod2[iPhi][iEta]->SetMarkerSize(1.);
     ICCrystalEB_Mod2[iPhi][iEta]->SetMarkerColor(kRed);
     GaussianFits_Mod2[iPhi][iEta]->SetLineColor(kBlue);
     GaussianFits_Mod2[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Mod2[iPhi][iEta]->Fit(GaussianFits_Mod2[iPhi][iEta],"RMEQ");
     ICCrystalEB_Mod2[iPhi][iEta]->Draw("E");
     ICCrystalEB_Mod2[iPhi][iEta]->GetXaxis()->SetTitle("IC");
     ICCrystalEB_Mod2[iPhi][iEta]->GetYaxis()->SetTitle("Entries");
     (Can_Mod2.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Mod2.at(iPhi)).at(iEta)->GetName())+".pdf").c_str(),"pdf");
     (Can_Mod2.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Mod2.at(iPhi)).at(iEta)->GetName())+".png").c_str(),"png");

     /// Original IC distribution Mod 3

     (Can_Mod3.at(iPhi)).at(iEta)->cd();
     (Can_Mod3.at(iPhi)).at(iEta)->SetGridx();
     (Can_Mod3.at(iPhi)).at(iEta)->SetGridy();
     ICCrystalEB_Mod3[iPhi][iEta]->SetLineColor(kBlack);
     ICCrystalEB_Mod3[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Mod3[iPhi][iEta]->SetMarkerStyle(20);
     ICCrystalEB_Mod3[iPhi][iEta]->SetMarkerSize(1.);
     ICCrystalEB_Mod3[iPhi][iEta]->SetMarkerColor(kRed);
     GaussianFits_Mod3[iPhi][iEta]->SetLineColor(kBlue);
     GaussianFits_Mod3[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Mod3[iPhi][iEta]->Fit(GaussianFits_Mod3[iPhi][iEta],"RMEQ");
     ICCrystalEB_Mod3[iPhi][iEta]->Draw("E");
     ICCrystalEB_Mod3[iPhi][iEta]->GetXaxis()->SetTitle("IC");
     ICCrystalEB_Mod3[iPhi][iEta]->GetYaxis()->SetTitle("Entries");
     (Can_Mod3.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Mod3.at(iPhi)).at(iEta)->GetName())+".pdf").c_str(),"pdf");
     (Can_Mod3.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Mod3.at(iPhi)).at(iEta)->GetName())+".png").c_str(),"png");

     /// Original IC distribution Mod 4

     (Can_Mod4.at(iPhi)).at(iEta)->cd();
     (Can_Mod4.at(iPhi)).at(iEta)->SetGridx();
     (Can_Mod4.at(iPhi)).at(iEta)->SetGridy();
     ICCrystalEB_Mod4[iPhi][iEta]->SetLineColor(kBlack);
     ICCrystalEB_Mod4[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Mod4[iPhi][iEta]->SetMarkerStyle(20);
     ICCrystalEB_Mod4[iPhi][iEta]->SetMarkerSize(1.);
     ICCrystalEB_Mod4[iPhi][iEta]->SetMarkerColor(kRed);
     GaussianFits_Mod4[iPhi][iEta]->SetLineColor(kBlue);
     GaussianFits_Mod4[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Mod4[iPhi][iEta]->Fit(GaussianFits_Mod4[iPhi][iEta],"RMEQ");
     ICCrystalEB_Mod4[iPhi][iEta]->Draw("E");
     ICCrystalEB_Mod4[iPhi][iEta]->GetXaxis()->SetTitle("IC");
     ICCrystalEB_Mod4[iPhi][iEta]->GetYaxis()->SetTitle("Entries");
     (Can_Mod4.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Mod4.at(iPhi)).at(iEta)->GetName())+".pdf").c_str(),"pdf");
     (Can_Mod4.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Mod4.at(iPhi)).at(iEta)->GetName())+".png").c_str(),"png");

     /// Mean IC vs Eta
    
     (Can_MeanIC_vs_Eta.at(iPhi)).at(iEta)->cd();
     (Can_MeanIC_vs_Eta.at(iPhi)).at(iEta)->SetGridx();
     (Can_MeanIC_vs_Eta.at(iPhi)).at(iEta)->SetGridy();

     MeanIC_vs_Eta[iPhi][iEta]->SetLineColor(kBlack);
     MeanIC_vs_Eta[iPhi][iEta]->SetLineWidth(2);
     MeanIC_vs_Eta[iPhi][iEta]->SetMarkerStyle(20);
     MeanIC_vs_Eta[iPhi][iEta]->SetMarkerSize(1.);
     //     MeanIC_vs_Eta[iPhi][iEta]->SetMarkerColor(kRed);
     MeanIC_vs_Eta[iPhi][iEta]->SetBinContent(1,ICCrystalEB_Mod1[iPhi][iEta]->GetMean());
     MeanIC_vs_Eta[iPhi][iEta]->SetBinContent(2,ICCrystalEB_Mod2[iPhi][iEta]->GetMean());
     MeanIC_vs_Eta[iPhi][iEta]->SetBinContent(3,ICCrystalEB_Mod3[iPhi][iEta]->GetMean());
     MeanIC_vs_Eta[iPhi][iEta]->SetBinContent(4,ICCrystalEB_Mod4[iPhi][iEta]->GetMean());
     MeanIC_vs_Eta[iPhi][iEta]->SetBinError(1,ICCrystalEB_Mod1[iPhi][iEta]->GetMeanError());
     MeanIC_vs_Eta[iPhi][iEta]->SetBinError(2,ICCrystalEB_Mod2[iPhi][iEta]->GetMeanError());
     MeanIC_vs_Eta[iPhi][iEta]->SetBinError(3,ICCrystalEB_Mod3[iPhi][iEta]->GetMeanError());
     MeanIC_vs_Eta[iPhi][iEta]->SetBinError(4,ICCrystalEB_Mod4[iPhi][iEta]->GetMeanError());

     Pol0_MeanIC_vs_Eta[iPhi][iEta] -> SetLineColor(kBlue);
     Pol1_MeanIC_vs_Eta[iPhi][iEta] -> SetLineColor(kGreen+1);
     Pol2_MeanIC_vs_Eta[iPhi][iEta] -> SetLineColor(kRed);
     Pol0_MeanIC_vs_Eta[iPhi][iEta] -> SetLineWidth(2);
     Pol1_MeanIC_vs_Eta[iPhi][iEta] -> SetLineWidth(2);
     Pol2_MeanIC_vs_Eta[iPhi][iEta] -> SetLineWidth(2);

     MeanIC_vs_Eta[iPhi][iEta]->Fit(Pol2_MeanIC_vs_Eta[iPhi][iEta],"RMEQ0");
     htemp1 = (TH1F*) MeanIC_vs_Eta[iPhi][iEta]->Clone("htemp1");
     htemp1->Fit(Pol1_MeanIC_vs_Eta[iPhi][iEta],"RMEQ0");
     htemp2 = (TH1F*) MeanIC_vs_Eta[iPhi][iEta]->Clone("htemp2");
     htemp2->Fit(Pol0_MeanIC_vs_Eta[iPhi][iEta],"RMEQ0");

     MeanIC_vs_Eta[iPhi][iEta]->Draw("E");
     MeanIC_vs_Eta[iPhi][iEta]->GetXaxis()->SetTitle("#eta module");
     MeanIC_vs_Eta[iPhi][iEta]->GetYaxis()->SetTitle("IC Mean");
     MeanIC_vs_Eta[iPhi][iEta]->GetYaxis()->SetRangeUser(MeanIC_vs_Eta[iPhi][iEta]->GetMinimum()*0.95,MeanIC_vs_Eta[iPhi][iEta]->GetMaximum()*1.15);

     (Can_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Update();
     pave1=(TPaveStats*)(MeanIC_vs_Eta[iPhi][iEta]->GetListOfFunctions()->FindObject("stats"));
     pave1->SetTextColor(kRed); 
     pave1->SetBorderSize(2); 
     pave1->SetLineColor(kRed); 
     pave1->SetFillStyle(0); 
     pave1->SetY1NDC(0.6);
     pave1->SetY2NDC(0.9);
     (Can_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Modified();

     htemp1->Draw("Esames");
     (Can_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Update();
     pave2=(TPaveStats*)(htemp1->GetListOfFunctions()->FindObject("stats"));
     pave2->SetOptStat(000);
     pave2->SetBorderSize(2); 
     pave2->SetLineColor(kGreen+1); 
     pave2->SetFillStyle(0); 
     pave2->SetTextSize(pave2->GetTextSize()*1.7); 
     pave2->SetTextColor(kGreen+1); 
     pave2->SetX1NDC(0.35);
     pave2->SetX2NDC(0.6);
     pave2->SetY1NDC(0.6);
     pave2->SetY2NDC(0.9);
     (Can_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Modified();

     htemp2->Draw("Esames");
     (Can_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Update();
     pave3=(TPaveStats*)(htemp2->GetListOfFunctions()->FindObject("stats"));
     pave3->SetOptStat(000);
     pave3->SetBorderSize(2); 
     pave3->SetLineColor(kBlue); 
     pave3->SetFillStyle(0); 
     pave3->SetTextSize(pave3->GetTextSize()*1.7); 
     pave3->SetTextColor(kBlue); 
     pave3->SetX1NDC(0.11);
     pave3->SetX2NDC(0.33);
     pave3->SetY1NDC(0.6);
     pave3->SetY2NDC(0.9);
     (Can_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Modified();
     
     Pol0_MeanIC_vs_Eta[iPhi][iEta] -> Draw("Lsame");
     Pol1_MeanIC_vs_Eta[iPhi][iEta] -> Draw("Lsame");
     Pol2_MeanIC_vs_Eta[iPhi][iEta] -> Draw("Lsame");

     (Can_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_MeanIC_vs_Eta.at(iPhi)).at(iEta)->GetName())+".pdf").c_str(),"pdf");
     (Can_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_MeanIC_vs_Eta.at(iPhi)).at(iEta)->GetName())+".png").c_str(),"png");
     
     if(htemp1!=0) delete htemp1 ;
     if(htemp2!=0) delete htemp2 ; 
         
     /// Gauss Mean  vs Eta
    
     (Can_Gaus_MeanIC_vs_Eta.at(iPhi)).at(iEta)->cd();
     (Can_Gaus_MeanIC_vs_Eta.at(iPhi)).at(iEta)->SetGridx();
     (Can_Gaus_MeanIC_vs_Eta.at(iPhi)).at(iEta)->SetGridy();

     Gaus_MeanIC_vs_Eta[iPhi][iEta]->SetLineColor(kBlack);
     Gaus_MeanIC_vs_Eta[iPhi][iEta]->SetLineWidth(2);
     Gaus_MeanIC_vs_Eta[iPhi][iEta]->SetMarkerStyle(20);
     Gaus_MeanIC_vs_Eta[iPhi][iEta]->SetMarkerSize(1.);
     //     Gaus_MeanIC_vs_Eta[iPhi][iEta]->SetMarkerColor(kRed);
     Gaus_MeanIC_vs_Eta[iPhi][iEta]->SetBinContent(1,GaussianFits_Mod1[iPhi][iEta]->GetParameter(1));
     Gaus_MeanIC_vs_Eta[iPhi][iEta]->SetBinContent(2,GaussianFits_Mod2[iPhi][iEta]->GetParameter(1));
     Gaus_MeanIC_vs_Eta[iPhi][iEta]->SetBinContent(3,GaussianFits_Mod3[iPhi][iEta]->GetParameter(1));
     Gaus_MeanIC_vs_Eta[iPhi][iEta]->SetBinContent(4,GaussianFits_Mod4[iPhi][iEta]->GetParameter(1));
     Gaus_MeanIC_vs_Eta[iPhi][iEta]->SetBinError(1,GaussianFits_Mod1[iPhi][iEta]->GetParError(1));
     Gaus_MeanIC_vs_Eta[iPhi][iEta]->SetBinError(2,GaussianFits_Mod2[iPhi][iEta]->GetParError(1));
     Gaus_MeanIC_vs_Eta[iPhi][iEta]->SetBinError(3,GaussianFits_Mod3[iPhi][iEta]->GetParError(1));
     Gaus_MeanIC_vs_Eta[iPhi][iEta]->SetBinError(4,GaussianFits_Mod4[iPhi][iEta]->GetParError(1));

     Pol0_Gaus_MeanIC_vs_Eta[iPhi][iEta] -> SetLineColor(kBlue);
     Pol1_Gaus_MeanIC_vs_Eta[iPhi][iEta] -> SetLineColor(kGreen+1);
     Pol2_Gaus_MeanIC_vs_Eta[iPhi][iEta] -> SetLineColor(kRed);
     Pol0_Gaus_MeanIC_vs_Eta[iPhi][iEta] -> SetLineWidth(2);
     Pol1_Gaus_MeanIC_vs_Eta[iPhi][iEta] -> SetLineWidth(2);
     Pol2_Gaus_MeanIC_vs_Eta[iPhi][iEta] -> SetLineWidth(2);

     Gaus_MeanIC_vs_Eta[iPhi][iEta]->Fit(Pol2_Gaus_MeanIC_vs_Eta[iPhi][iEta],"RMEQ0");
     htemp1 = (TH1F*) Gaus_MeanIC_vs_Eta[iPhi][iEta]->Clone("htemp1");
     htemp1->Fit(Pol1_Gaus_MeanIC_vs_Eta[iPhi][iEta],"RMEQ0");
     htemp2 = (TH1F*) Gaus_MeanIC_vs_Eta[iPhi][iEta]->Clone("htemp2");
     htemp2->Fit(Pol0_Gaus_MeanIC_vs_Eta[iPhi][iEta],"RMEQ0");

     Gaus_MeanIC_vs_Eta[iPhi][iEta]->Draw("E");
     Gaus_MeanIC_vs_Eta[iPhi][iEta]->GetXaxis()->SetTitle("#eta module");
     Gaus_MeanIC_vs_Eta[iPhi][iEta]->GetYaxis()->SetTitle("IC Gaussian Mean");
     Gaus_MeanIC_vs_Eta[iPhi][iEta]->GetYaxis()->SetRangeUser(Gaus_MeanIC_vs_Eta[iPhi][iEta]->GetMinimum()*0.95,Gaus_MeanIC_vs_Eta[iPhi][iEta]->GetMaximum()*1.15);

     (Can_Gaus_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Update();
     pave1=(TPaveStats*)(Gaus_MeanIC_vs_Eta[iPhi][iEta]->GetListOfFunctions()->FindObject("stats"));
     pave1->SetTextColor(kRed); 
     pave1->SetBorderSize(2); 
     pave1->SetLineColor(kRed); 
     pave1->SetFillStyle(0); 
     pave1->SetY1NDC(0.6);
     pave1->SetY2NDC(0.9);
     (Can_Gaus_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Modified();
        
     htemp1->Draw("Esames");
     (Can_Gaus_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Update();
     pave2=(TPaveStats*)(htemp1->GetListOfFunctions()->FindObject("stats"));
     pave2->SetOptStat(000);
     pave2->SetBorderSize(2); 
     pave2->SetLineColor(kGreen+1); 
     pave2->SetFillStyle(0); 
     pave2->SetTextSize(pave2->GetTextSize()*1.7); 
     pave2->SetTextColor(kGreen+1); 
     pave2->SetX1NDC(0.35);
     pave2->SetX2NDC(0.6);
     pave2->SetY1NDC(0.6);
     pave2->SetY2NDC(0.9);
     (Can_Gaus_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Modified();

     htemp2->Draw("Esames");
     (Can_Gaus_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Update();
     pave3=(TPaveStats*)(htemp2->GetListOfFunctions()->FindObject("stats"));
     pave3->SetOptStat(000);
     pave3->SetBorderSize(2); 
     pave3->SetLineColor(kBlue); 
     pave3->SetFillStyle(0); 
     pave3->SetTextSize(pave3->GetTextSize()*1.7); 
     pave3->SetTextColor(kBlue); 
     pave3->SetX1NDC(0.11);
     pave3->SetX2NDC(0.33);
     pave3->SetY1NDC(0.6);
     pave3->SetY2NDC(0.9);
     (Can_Gaus_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Modified();

     Pol0_Gaus_MeanIC_vs_Eta[iPhi][iEta] -> Draw("Lsame");
     Pol1_Gaus_MeanIC_vs_Eta[iPhi][iEta] -> Draw("Lsame");
     Pol2_Gaus_MeanIC_vs_Eta[iPhi][iEta] -> Draw("Lsame");

     (Can_Gaus_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Gaus_MeanIC_vs_Eta.at(iPhi)).at(iEta)->GetName())+".pdf").c_str(),"pdf");
     (Can_Gaus_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Gaus_MeanIC_vs_Eta.at(iPhi)).at(iEta)->GetName())+".png").c_str(),"png");
     (Can_Gaus_MeanIC_vs_Eta.at(iPhi)).at(iEta)->Close();

     delete htemp1 ;
     delete htemp2 ;
     
     /// Original IC %20 < 5

     (Can_Phi1.at(iPhi)).at(iEta)->cd();
     (Can_Phi1.at(iPhi)).at(iEta)->SetGridx();
     (Can_Phi1.at(iPhi)).at(iEta)->SetGridy();
     ICCrystalEB_Phi1[iPhi][iEta]->SetLineColor(kBlack);
     ICCrystalEB_Phi1[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Phi1[iPhi][iEta]->SetMarkerStyle(20);
     ICCrystalEB_Phi1[iPhi][iEta]->SetMarkerSize(1.);
     ICCrystalEB_Phi1[iPhi][iEta]->SetMarkerColor(kRed);
     GaussianFits_Phi1[iPhi][iEta]->SetLineColor(kBlue);
     GaussianFits_Phi1[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Phi1[iPhi][iEta]->Fit(GaussianFits_Phi1[iPhi][iEta],"RMEQ");
     ICCrystalEB_Phi1[iPhi][iEta]->Draw("E");
     ICCrystalEB_Phi1[iPhi][iEta]->GetXaxis()->SetTitle("IC");
     ICCrystalEB_Phi1[iPhi][iEta]->GetYaxis()->SetTitle("Entries");
     (Can_Phi1.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Phi1.at(iPhi)).at(iEta)->GetName())+".pdf").c_str(),"pdf");
     (Can_Phi1.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Phi1.at(iPhi)).at(iEta)->GetName())+".png").c_str(),"png");

     /// Original IC %20 in [5,10]

     (Can_Phi2.at(iPhi)).at(iEta)->cd();
     (Can_Phi2.at(iPhi)).at(iEta)->SetGridx();
     (Can_Phi2.at(iPhi)).at(iEta)->SetGridy();
     ICCrystalEB_Phi2[iPhi][iEta]->SetLineColor(kBlack);
     ICCrystalEB_Phi2[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Phi2[iPhi][iEta]->SetMarkerStyle(20);
     ICCrystalEB_Phi2[iPhi][iEta]->SetMarkerSize(1.);
     ICCrystalEB_Phi2[iPhi][iEta]->SetMarkerColor(kRed);
     GaussianFits_Phi2[iPhi][iEta]->SetLineColor(kBlue);
     GaussianFits_Phi2[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Phi2[iPhi][iEta]->Fit(GaussianFits_Phi2[iPhi][iEta],"RMEQ");
     ICCrystalEB_Phi2[iPhi][iEta]->Draw("E");
     ICCrystalEB_Phi2[iPhi][iEta]->GetXaxis()->SetTitle("IC");
     ICCrystalEB_Phi2[iPhi][iEta]->GetYaxis()->SetTitle("Entries");
     (Can_Phi2.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Phi2.at(iPhi)).at(iEta)->GetName())+".pdf").c_str(),"pdf");
     (Can_Phi2.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Phi2.at(iPhi)).at(iEta)->GetName())+".png").c_str(),"png");

     /// Original IC %20 in [10,15]

     (Can_Phi3.at(iPhi)).at(iEta)->cd();
     (Can_Phi3.at(iPhi)).at(iEta)->SetGridx();
     (Can_Phi3.at(iPhi)).at(iEta)->SetGridy();
     ICCrystalEB_Phi3[iPhi][iEta]->SetLineColor(kBlack);
     ICCrystalEB_Phi3[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Phi3[iPhi][iEta]->SetMarkerStyle(20);
     ICCrystalEB_Phi3[iPhi][iEta]->SetMarkerSize(1.);
     ICCrystalEB_Phi3[iPhi][iEta]->SetMarkerColor(kRed);
     GaussianFits_Phi3[iPhi][iEta]->SetLineColor(kBlue);
     GaussianFits_Phi3[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Phi3[iPhi][iEta]->Fit(GaussianFits_Phi3[iPhi][iEta],"RMEQ");
     ICCrystalEB_Phi3[iPhi][iEta]->Draw("E");
     ICCrystalEB_Phi3[iPhi][iEta]->GetXaxis()->SetTitle("IC");
     ICCrystalEB_Phi3[iPhi][iEta]->GetYaxis()->SetTitle("Entries");
     (Can_Phi3.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Phi3.at(iPhi)).at(iEta)->GetName())+".pdf").c_str(),"pdf");
     (Can_Phi3.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Phi3.at(iPhi)).at(iEta)->GetName())+".png").c_str(),"png");

     /// Original IC %20 in [15,20]

     (Can_Phi4.at(iPhi)).at(iEta)->cd();
     (Can_Phi4.at(iPhi)).at(iEta)->SetGridx();
     (Can_Phi4.at(iPhi)).at(iEta)->SetGridy();
     ICCrystalEB_Phi4[iPhi][iEta]->SetLineColor(kBlack);
     ICCrystalEB_Phi4[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Phi4[iPhi][iEta]->SetMarkerStyle(20);
     ICCrystalEB_Phi4[iPhi][iEta]->SetMarkerSize(1.);
     ICCrystalEB_Phi4[iPhi][iEta]->SetMarkerColor(kRed);
     GaussianFits_Phi4[iPhi][iEta]->SetLineColor(kBlue);
     GaussianFits_Phi4[iPhi][iEta]->SetLineWidth(2);
     ICCrystalEB_Phi4[iPhi][iEta]->Fit(GaussianFits_Phi4[iPhi][iEta],"RMEQ");
     ICCrystalEB_Phi4[iPhi][iEta]->Draw("E");
     ICCrystalEB_Phi4[iPhi][iEta]->GetXaxis()->SetTitle("IC");
     ICCrystalEB_Phi4[iPhi][iEta]->GetYaxis()->SetTitle("Entries");
     (Can_Phi4.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Phi4.at(iPhi)).at(iEta)->GetName())+".pdf").c_str(),"pdf");
     (Can_Phi4.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Phi4.at(iPhi)).at(iEta)->GetName())+".png").c_str(),"png");
     

     /// Mean IC vs Phi

     (Can_MeanIC_vs_Phi.at(iPhi)).at(iEta)->cd();
     (Can_MeanIC_vs_Phi.at(iPhi)).at(iEta)->SetGridx();
     (Can_MeanIC_vs_Phi.at(iPhi)).at(iEta)->SetGridy();

     MeanIC_vs_Phi[iPhi][iEta]->SetLineColor(kBlack);
     MeanIC_vs_Phi[iPhi][iEta]->SetLineWidth(2);
     MeanIC_vs_Phi[iPhi][iEta]->SetMarkerStyle(20);
     MeanIC_vs_Phi[iPhi][iEta]->SetMarkerSize(1.);
     //     MeanIC_vs_Phi[iPhi][iEta]->SetMarkerColor(kRed);
     MeanIC_vs_Phi[iPhi][iEta]->SetBinContent(1,ICCrystalEB_Phi1[iPhi][iEta]->GetMean());
     MeanIC_vs_Phi[iPhi][iEta]->SetBinContent(2,ICCrystalEB_Phi2[iPhi][iEta]->GetMean());
     MeanIC_vs_Phi[iPhi][iEta]->SetBinContent(3,ICCrystalEB_Phi3[iPhi][iEta]->GetMean());
     MeanIC_vs_Phi[iPhi][iEta]->SetBinContent(4,ICCrystalEB_Phi4[iPhi][iEta]->GetMean());
     MeanIC_vs_Phi[iPhi][iEta]->SetBinError(1,ICCrystalEB_Phi1[iPhi][iEta]->GetMeanError());
     MeanIC_vs_Phi[iPhi][iEta]->SetBinError(2,ICCrystalEB_Phi2[iPhi][iEta]->GetMeanError());
     MeanIC_vs_Phi[iPhi][iEta]->SetBinError(3,ICCrystalEB_Phi3[iPhi][iEta]->GetMeanError());
     MeanIC_vs_Phi[iPhi][iEta]->SetBinError(4,ICCrystalEB_Phi4[iPhi][iEta]->GetMeanError());

     Pol0_MeanIC_vs_Phi[iPhi][iEta] -> SetLineColor(kBlue);
     Pol1_MeanIC_vs_Phi[iPhi][iEta] -> SetLineColor(kGreen+1);
     Pol2_MeanIC_vs_Phi[iPhi][iEta] -> SetLineColor(kRed);
     Pol0_MeanIC_vs_Phi[iPhi][iEta] -> SetLineWidth(2);
     Pol1_MeanIC_vs_Phi[iPhi][iEta] -> SetLineWidth(2);
     Pol2_MeanIC_vs_Phi[iPhi][iEta] -> SetLineWidth(2);

     MeanIC_vs_Phi[iPhi][iEta]->Fit(Pol2_MeanIC_vs_Phi[iPhi][iEta],"RMEQ0");
     htemp1 = (TH1F*) MeanIC_vs_Phi[iPhi][iEta]->Clone("htemp1");
     htemp1->Fit(Pol1_MeanIC_vs_Phi[iPhi][iEta],"RMEQ0");
     htemp2 = (TH1F*) MeanIC_vs_Phi[iPhi][iEta]->Clone("htemp2");
     htemp2->Fit(Pol0_MeanIC_vs_Phi[iPhi][iEta],"RMEQ0");

     MeanIC_vs_Phi[iPhi][iEta]->Draw("E");
     MeanIC_vs_Phi[iPhi][iEta]->GetXaxis()->SetTitle("#phi bin SM");
     MeanIC_vs_Phi[iPhi][iEta]->GetYaxis()->SetTitle("IC Mean");
     MeanIC_vs_Phi[iPhi][iEta]->GetYaxis()->SetRangeUser(MeanIC_vs_Phi[iPhi][iEta]->GetMinimum()*0.95,MeanIC_vs_Phi[iPhi][iEta]->GetMaximum()*1.15);

     (Can_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Update();
     pave1=(TPaveStats*)(MeanIC_vs_Phi[iPhi][iEta]->GetListOfFunctions()->FindObject("stats"));
     pave1->SetTextColor(kRed); 
     pave1->SetBorderSize(2); 
     pave1->SetLineColor(kRed); 
     pave1->SetFillStyle(0); 
     pave1->SetY1NDC(0.6);
     pave1->SetY2NDC(0.9);
     (Can_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Modified();

     htemp1->Draw("Esames");
     (Can_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Update();
     pave2=(TPaveStats*)(htemp1->GetListOfFunctions()->FindObject("stats"));
     pave2->SetOptStat(000);
     pave2->SetBorderSize(2); 
     pave2->SetLineColor(kGreen+1); 
     pave2->SetFillStyle(0); 
     pave2->SetTextSize(pave2->GetTextSize()*1.7); 
     pave2->SetTextColor(kGreen+1); 
     pave2->SetX1NDC(0.35);
     pave2->SetX2NDC(0.6);
     pave2->SetY1NDC(0.6);
     pave2->SetY2NDC(0.9);
     (Can_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Modified();

     htemp2->Draw("Esames");
     (Can_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Update();
     pave3=(TPaveStats*)(htemp2->GetListOfFunctions()->FindObject("stats"));
     pave3->SetOptStat(000);
     pave3->SetBorderSize(2); 
     pave3->SetLineColor(kBlue); 
     pave3->SetFillStyle(0); 
     pave3->SetTextSize(pave3->GetTextSize()*1.7); 
     pave3->SetTextColor(kBlue); 
     pave3->SetX1NDC(0.11);
     pave3->SetX2NDC(0.33);
     pave3->SetY1NDC(0.6);
     pave3->SetY2NDC(0.9);
     (Can_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Modified();
     
     Pol0_MeanIC_vs_Phi[iPhi][iEta] -> Draw("Lsame");
     Pol1_MeanIC_vs_Phi[iPhi][iEta] -> Draw("Lsame");
     Pol2_MeanIC_vs_Phi[iPhi][iEta] -> Draw("Lsame");

     (Can_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_MeanIC_vs_Phi.at(iPhi)).at(iEta)->GetName())+".pdf").c_str(),"pdf");
     (Can_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_MeanIC_vs_Phi.at(iPhi)).at(iEta)->GetName())+".png").c_str(),"png");
     
     if(htemp1!=0) delete htemp1 ;
     if(htemp2!=0) delete htemp2 ; 

     /// Gauss Mean  vs Phi
    
     (Can_Gaus_MeanIC_vs_Phi.at(iPhi)).at(iEta)->cd();
     (Can_Gaus_MeanIC_vs_Phi.at(iPhi)).at(iEta)->SetGridx();
     (Can_Gaus_MeanIC_vs_Phi.at(iPhi)).at(iEta)->SetGridy();

     Gaus_MeanIC_vs_Phi[iPhi][iEta]->SetLineColor(kBlack);
     Gaus_MeanIC_vs_Phi[iPhi][iEta]->SetLineWidth(2);
     Gaus_MeanIC_vs_Phi[iPhi][iEta]->SetMarkerStyle(20);
     Gaus_MeanIC_vs_Phi[iPhi][iEta]->SetMarkerSize(1.);
     //     Gaus_MeanIC_vs_Phi[iPhi][iEta]->SetMarkerColor(kRed);
     Gaus_MeanIC_vs_Phi[iPhi][iEta]->SetBinContent(1,GaussianFits_Phi1[iPhi][iEta]->GetParameter(1));
     Gaus_MeanIC_vs_Phi[iPhi][iEta]->SetBinContent(2,GaussianFits_Phi2[iPhi][iEta]->GetParameter(1));
     Gaus_MeanIC_vs_Phi[iPhi][iEta]->SetBinContent(3,GaussianFits_Phi3[iPhi][iEta]->GetParameter(1));
     Gaus_MeanIC_vs_Phi[iPhi][iEta]->SetBinContent(4,GaussianFits_Phi4[iPhi][iEta]->GetParameter(1));
     Gaus_MeanIC_vs_Phi[iPhi][iEta]->SetBinError(1,GaussianFits_Phi1[iPhi][iEta]->GetParError(1));
     Gaus_MeanIC_vs_Phi[iPhi][iEta]->SetBinError(2,GaussianFits_Phi2[iPhi][iEta]->GetParError(1));
     Gaus_MeanIC_vs_Phi[iPhi][iEta]->SetBinError(3,GaussianFits_Phi3[iPhi][iEta]->GetParError(1));
     Gaus_MeanIC_vs_Phi[iPhi][iEta]->SetBinError(4,GaussianFits_Phi4[iPhi][iEta]->GetParError(1));

     Pol0_Gaus_MeanIC_vs_Phi[iPhi][iEta] -> SetLineColor(kBlue);
     Pol1_Gaus_MeanIC_vs_Phi[iPhi][iEta] -> SetLineColor(kGreen+1);
     Pol2_Gaus_MeanIC_vs_Phi[iPhi][iEta] -> SetLineColor(kRed);
     Pol0_Gaus_MeanIC_vs_Phi[iPhi][iEta] -> SetLineWidth(2);
     Pol1_Gaus_MeanIC_vs_Phi[iPhi][iEta] -> SetLineWidth(2);
     Pol2_Gaus_MeanIC_vs_Phi[iPhi][iEta] -> SetLineWidth(2);

     Gaus_MeanIC_vs_Phi[iPhi][iEta]->Fit(Pol2_Gaus_MeanIC_vs_Phi[iPhi][iEta],"RMEQ0");
     htemp1 = (TH1F*) Gaus_MeanIC_vs_Phi[iPhi][iEta]->Clone("htemp1");
     htemp1->Fit(Pol1_Gaus_MeanIC_vs_Phi[iPhi][iEta],"RMEQ0");
     htemp2 = (TH1F*) Gaus_MeanIC_vs_Phi[iPhi][iEta]->Clone("htemp2");
     htemp2->Fit(Pol0_Gaus_MeanIC_vs_Phi[iPhi][iEta],"RMEQ0");

     Gaus_MeanIC_vs_Phi[iPhi][iEta]->Draw("E");
     Gaus_MeanIC_vs_Phi[iPhi][iEta]->GetXaxis()->SetTitle("#eta module");
     Gaus_MeanIC_vs_Phi[iPhi][iEta]->GetYaxis()->SetTitle("IC Gaussian Mean");
     Gaus_MeanIC_vs_Phi[iPhi][iEta]->GetYaxis()->SetRangeUser(Gaus_MeanIC_vs_Phi[iPhi][iEta]->GetMinimum()*0.95,Gaus_MeanIC_vs_Phi[iPhi][iEta]->GetMaximum()*1.15);

     (Can_Gaus_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Update();
     pave1=(TPaveStats*)(Gaus_MeanIC_vs_Phi[iPhi][iEta]->GetListOfFunctions()->FindObject("stats"));
     pave1->SetTextColor(kRed); 
     pave1->SetBorderSize(2); 
     pave1->SetLineColor(kRed); 
     pave1->SetFillStyle(0); 
     pave1->SetY1NDC(0.6);
     pave1->SetY2NDC(0.9);
     (Can_Gaus_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Modified();
        
     htemp1->Draw("Esames");
     (Can_Gaus_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Update();
     pave2=(TPaveStats*)(htemp1->GetListOfFunctions()->FindObject("stats"));
     pave2->SetOptStat(000);
     pave2->SetBorderSize(2); 
     pave2->SetLineColor(kGreen+1); 
     pave2->SetFillStyle(0); 
     pave2->SetTextSize(pave2->GetTextSize()*1.7); 
     pave2->SetTextColor(kGreen+1); 
     pave2->SetX1NDC(0.35);
     pave2->SetX2NDC(0.6);
     pave2->SetY1NDC(0.6);
     pave2->SetY2NDC(0.9);
     (Can_Gaus_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Modified();

     htemp2->Draw("Esames");
     (Can_Gaus_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Update();
     pave3=(TPaveStats*)(htemp2->GetListOfFunctions()->FindObject("stats"));
     pave3->SetOptStat(000);
     pave3->SetBorderSize(2); 
     pave3->SetLineColor(kBlue); 
     pave3->SetFillStyle(0); 
     pave3->SetTextSize(pave3->GetTextSize()*1.7); 
     pave3->SetTextColor(kBlue); 
     pave3->SetX1NDC(0.11);
     pave3->SetX2NDC(0.33);
     pave3->SetY1NDC(0.6);
     pave3->SetY2NDC(0.9);
     (Can_Gaus_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Modified();

     Pol0_Gaus_MeanIC_vs_Phi[iPhi][iEta] -> Draw("Lsame");
     Pol1_Gaus_MeanIC_vs_Phi[iPhi][iEta] -> Draw("Lsame");
     Pol2_Gaus_MeanIC_vs_Phi[iPhi][iEta] -> Draw("Lsame");

     (Can_Gaus_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Gaus_MeanIC_vs_Phi.at(iPhi)).at(iEta)->GetName())+".pdf").c_str(),"pdf");
     (Can_Gaus_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Print((outputCanvasPlot+std::string((Can_Gaus_MeanIC_vs_Phi.at(iPhi)).at(iEta)->GetName())+".png").c_str(),"png");
     (Can_Gaus_MeanIC_vs_Phi.at(iPhi)).at(iEta)->Close();

     delete htemp1 ;
     delete htemp2 ;

     
   }

 }


 return 0; 

}
