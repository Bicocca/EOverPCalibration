#include "TEndcapRings.h"
#include "ntpleUtils.h"
#include "treeReader.h"
#include "CalibrationUtils.h"
#include "../CommonTools/histoFunc.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "ConfigParser.h"
#include "TMath.h"

using namespace std;

 
bool IsEtaGap(float eta)
{
  float feta = fabs(eta);
  if( fabs(feta - 0 ) < 3 ) return true;
  if( fabs(feta - 25) < 3 ) return true;
  if( fabs(feta - 45) < 3 ) return true;
  if( fabs(feta - 65) < 3 ) return true;
  if( fabs(feta - 85) < 3 ) return true;
  return false;
}






//**************  MAIN PROGRAM **************************************************************
int main(int argc, char** argv)
{
  // Acquisition from cfg file
  if(argc != 2)
  {
    std::cerr << ">>>>> CalibrationMomentum.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return 1;
  }
  
  parseConfigFile(argv[1]);
  
  std::string TreeName    = gConfigParser -> readStringOption("Input::TreeName");
  std::string infileDATA  = gConfigParser -> readStringOption("Input::infileDATA");
  std::string infileMC    = gConfigParser -> readStringOption("Input::infileMC");
  std::string WeightforMC = gConfigParser -> readStringOption("Input::WeightforMC");
  
  std::string typeEB = gConfigParser -> readStringOption("Input::typeEB");
  std::string typeEE = gConfigParser -> readStringOption("Input::typeEE");
  int  nPhiBinsEB = gConfigParser -> readIntOption("Input::nPhiBinsEB");
  int  nPhiBinsEE = gConfigParser -> readIntOption("Input::nPhiBinsEE");
  int  nPhiBinsTempEB = gConfigParser -> readIntOption("Input::nPhiBinsTempEB");
  int  nPhiBinsTempEE = gConfigParser -> readIntOption("Input::nPhiBinsTempEE");
  int  rebinEB = gConfigParser -> readIntOption("Input::rebinEB");
  int  rebinEE = gConfigParser -> readIntOption("Input::rebinEE");
  bool usePUweights = gConfigParser -> readBoolOption("Input::usePUweights");
  std::string outputFile = gConfigParser -> readStringOption("Output::outputFile");
  
  int nRegionsEB = GetNRegionsEB(typeEB);
  int nRegionsEE = GetNRegionsEE(typeEE);

  std::cout<<"REGIONI: "<<nRegionsEE<<std::endl;
  getchar();

  cout <<" Basic Configuration " <<endl;
  cout <<" Tree Name = "<<TreeName<<endl;
  cout <<" infileDATA = "<<infileDATA<<endl;
  cout <<" infileMC = "<<infileMC<<endl;
  cout <<" WeightforMC = "<<WeightforMC<<endl;
  cout <<" nRegionsEB = "<<nRegionsEB<<endl;
  cout <<" nRegionsEE = "<<nRegionsEE<<endl;
  cout <<" nPhiBinsEB = "<<nPhiBinsEB<<endl;
  cout <<" nPhiBinsEE = "<<nPhiBinsEE<<endl;
  cout <<" nPhiBinsTempEB = "<<nPhiBinsTempEB<<endl;
  cout <<" nPhiBinsTempEE = "<<nPhiBinsTempEE<<endl;
  cout <<" rebinEB = "<<rebinEB<<endl;
  cout <<" rebinEE = "<<rebinEE<<endl;
  cout <<" usePUweights = "<<usePUweights<<endl;
  
  cout << "Momentum correction "<< endl;
  
  
  
  //---- variables for selection
  float etaMax  = 2.5;
  float eta2Max = 2.5;

  //histos to get the bin in phi given the electron phi
  TH1F* hPhiBinEE = new TH1F("hphiEE","hphiEE",nPhiBinsEE, -1.*TMath::Pi(),1.*TMath::Pi());
  
  //----- NTUPLES--------------------
  TChain *ntu_DA = new TChain(TreeName.c_str());
  
  if(!FillChain(*ntu_DA, infileDATA.c_str())) return 1;
  
  std::cout << "     DATA: " << ntu_DA->GetEntries() << " entries in Data sample" << std::endl;
  
  // observables  
  int isW;
  float mZ;
  float scEta,  scPhi;
  float scEta2, scPhi2;
  float eleEta,  elePhi;
  float eleEta2, elePhi2;
  float scEne,  scEneReg,  scEt,  scERaw,  e3x3,  R9;
  float scEne2, scEneReg2, scEt2, scERaw2, e3x32, R92;
  float charge, charge2;
  float pTK,pTK2; 
  int iphiSeed,  ele1_ix, ele1_iy, ele1_iz; 
  int iphiSeed2, ele2_ix, ele2_iy, ele2_iz;
  int npu;
  
  // Set branch addresses for Data  
  ntu_DA->SetBranchStatus("*",0);
  ntu_DA->SetBranchStatus("isW", 1);                 ntu_DA->SetBranchAddress("isW", &isW);
  ntu_DA->SetBranchStatus("ele1_eta", 1);            ntu_DA->SetBranchAddress("ele1_eta", &eleEta);
  ntu_DA->SetBranchStatus("ele2_eta", 1);            ntu_DA->SetBranchAddress("ele2_eta", &eleEta2);
  ntu_DA->SetBranchStatus("ele1_phi", 1);            ntu_DA->SetBranchAddress("ele1_phi", &elePhi);
  ntu_DA->SetBranchStatus("ele2_phi", 1);            ntu_DA->SetBranchAddress("ele2_phi", &elePhi2);
  ntu_DA->SetBranchStatus("ele1_scEta", 1);          ntu_DA->SetBranchAddress("ele1_scEta", &scEta);
  ntu_DA->SetBranchStatus("ele2_scEta", 1);          ntu_DA->SetBranchAddress("ele2_scEta", &scEta2);
  ntu_DA->SetBranchStatus("ele1_scPhi", 1);          ntu_DA->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_DA->SetBranchStatus("ele2_scPhi", 1);          ntu_DA->SetBranchAddress("ele2_scPhi", &scPhi2);
  ntu_DA->SetBranchStatus("ele1_scE", 1);            ntu_DA->SetBranchAddress("ele1_scE", &scEne);
  ntu_DA->SetBranchStatus("ele2_scE", 1);            ntu_DA->SetBranchAddress("ele2_scE", &scEne2);
  ntu_DA->SetBranchStatus("ele1_scE_regression", 1); ntu_DA->SetBranchAddress("ele1_scE_regression", &scEneReg);
  ntu_DA->SetBranchStatus("ele2_scE_regression", 1); ntu_DA->SetBranchAddress("ele2_scE_regression", &scEneReg2);
  ntu_DA->SetBranchStatus("ele1_scEt", 1);           ntu_DA->SetBranchAddress("ele1_scEt",&scEt);
  ntu_DA->SetBranchStatus("ele2_scEt", 1);           ntu_DA->SetBranchAddress("ele2_scEt",&scEt2);
  ntu_DA->SetBranchStatus("ele1_scERaw", 1);         ntu_DA->SetBranchAddress("ele1_scERaw", &scERaw);
  ntu_DA->SetBranchStatus("ele2_scERaw", 1);         ntu_DA->SetBranchAddress("ele2_scERaw", &scERaw2);
  ntu_DA->SetBranchStatus("ele1_e3x3", 1);           ntu_DA->SetBranchAddress("ele1_e3x3", &e3x3);
  ntu_DA->SetBranchStatus("ele2_e3x3", 1);           ntu_DA->SetBranchAddress("ele2_e3x3", &e3x32);
  ntu_DA->SetBranchStatus("ele1ele2_scM", 1);        ntu_DA->SetBranchAddress("ele1ele2_scM", &mZ);
  ntu_DA->SetBranchStatus("ele1_charge", 1);         ntu_DA->SetBranchAddress("ele1_charge", &charge);
  ntu_DA->SetBranchStatus("ele2_charge", 1);         ntu_DA->SetBranchAddress("ele2_charge", &charge2);
  ntu_DA->SetBranchStatus("ele1_tkP", 1);            ntu_DA->SetBranchAddress("ele1_tkP", &pTK);
  ntu_DA->SetBranchStatus("ele2_tkP", 1);            ntu_DA->SetBranchAddress("ele2_tkP", &pTK2);
  ntu_DA->SetBranchStatus("ele1_seedIphi", 1);       ntu_DA->SetBranchAddress("ele1_seedIphi", &iphiSeed);
  ntu_DA->SetBranchStatus("ele2_seedIphi", 1);       ntu_DA->SetBranchAddress("ele2_seedIphi", &iphiSeed2);
  ntu_DA->SetBranchStatus("ele1_seedIx", 1);         ntu_DA->SetBranchAddress("ele1_seedIx", &ele1_ix);
  ntu_DA->SetBranchStatus("ele2_seedIx", 1);         ntu_DA->SetBranchAddress("ele2_seedIx", &ele2_ix);
  ntu_DA->SetBranchStatus("ele1_seedIy", 1);         ntu_DA->SetBranchAddress("ele1_seedIy", &ele1_iy);
  ntu_DA->SetBranchStatus("ele2_seedIy", 1);         ntu_DA->SetBranchAddress("ele2_seedIy", &ele2_iy);
  ntu_DA->SetBranchStatus("ele1_seedZside", 1);      ntu_DA->SetBranchAddress("ele1_seedZside", &ele1_iz);
  ntu_DA->SetBranchStatus("ele2_seedZside", 1);      ntu_DA->SetBranchAddress("ele2_seedZside", &ele2_iz);
  
    
  // histogram definition in EE and fit functions
  std::vector<std::vector<TH1F*> > h_Phi_EE(nPhiBinsEE); // used to map iEta (as defined for Barrel and Endcap geom) into eta 
  std::vector<std::vector<TH1F*> > h_EoC_EE(nPhiBinsEE);  
  std::vector<std::vector<TF1*> > f_EoC_EE(nPhiBinsEE);
  
    
  nRegionsEE=2;
  
  // Initializate histos in EE
  std::cout << ">>> Initialize EE histos" << std::endl;
  for(int i = 0; i < nPhiBinsEE; ++i)
  {
    for(int j = 0; j < nRegionsEE; ++j)
    {
      TString histoName;
      histoName= Form("EE_EoP_%d_%d", i,j);
      TH1F* temp = new TH1F (histoName, histoName, 2200, 0., 100.);
      temp->Sumw2();
      temp->SetFillColor(kGreen+2);
      temp->SetLineColor(kGreen+2);
      temp->SetFillStyle(3004);
      (h_EoC_EE.at(i)).push_back(temp);
      
      histoName=Form("EE_Phi_%d_%d", i,j);
      temp = new TH1F(histoName, histoName, 360, 0., 360.); 
      (h_Phi_EE.at(i)).push_back(temp); 
    }
  }
  
  
  
  
  TH1F** h_phi_data_EE = new TH1F*[nRegionsEE];

  for(int index = 0; index < nRegionsEE; ++index)
  {
    TString name;
    name=Form("EE_h_phi_data_%d",index);
    h_phi_data_EE[index] = new TH1F(name,"h_phi_data",100,-TMath::Pi(),TMath::Pi());
  }
  
  TH1F* h_et_data = new TH1F("h_et_data","h_et_data",100,0.,100.);
    
  
  // Initialize endcap geometry
  TEndcapRings *eRings = new TEndcapRings(); 
  
  // Map for conversion (ix,iy) into Eta for EE
  TH2F * mapConversionEEp = new TH2F ("mapConversionEEp","mapConversionEEp",101,1,101,101,1,101);
  TH2F * mapConversionEEm = new TH2F ("mapConversionEEm","mapConversionEEm",101,1,101,101,1,101);
  
  for(int ix =0; ix<mapConversionEEp->GetNbinsX(); ix++)
    for(int iy =0; iy<mapConversionEEp->GetNbinsY(); iy++)
    {
      mapConversionEEp->SetBinContent(ix+1,iy+1,0);
      mapConversionEEm->SetBinContent(ix+1,iy+1,0);
    }
  
  
  
  // fill MC templates
  
  std::vector<int> refIdEE;
  refIdEE.assign(nPhiBinsEE,0);
  
  for(int iphi = 0; iphi < nPhiBinsEE; ++iphi)
  {
    float phi = hPhiBinEE->GetBinCenter(iphi+1);
    
    phi = 2.*TMath::Pi() + phi + TMath::Pi()*10./180.;
    phi -= int(phi/2./TMath::Pi()) * 2.*TMath::Pi();
    
    int modPhi = int(phi/(2.*TMath::Pi()/nPhiBinsTempEE));
    if( modPhi == nPhiBinsTempEE ) modPhi = 0;
    refIdEE.at(iphi) = modPhi;
  }
  
    
  
  //**************************** loop on DATA
  
  std::cout << "Loop in Data events " << endl; 
  
  for(int entry = 0; entry < ntu_DA->GetEntries(); ++entry)
  {
    if( entry%10000 == 0 ) std::cout << "reading saved entry " << entry << "\r" << std::flush;
    
    ntu_DA->GetEntry(entry);
    
    //    if( isW == 1 )               continue;
    if( fabs(scEta)  > etaMax )  continue;
    if( fabs(scEta2) > eta2Max ) continue;
    if( scEt  < 20. ) continue;
    if( scEt2 < 20. ) continue;    
    
    R9  = e3x3  / scERaw;
    R92 = e3x32 / scERaw2;

    if (R9<0.8 && R92<0.8) continue;
        
    float ww = 1.;
    int index=0;

    // DATA - ENDCAP - ele1
    if (ele1_iz!=0 && R9>=0.8)
    {
      if( ele1_iz ==  1 ) mapConversionEEp -> SetBinContent(ele1_ix,ele1_iy,scEta);
      if( ele1_iz == -1 ) mapConversionEEm -> SetBinContent(ele1_ix,ele1_iy,scEta);
      
      int PhibinEE = hPhiBinEE->FindBin(elePhi) - 1;
      if( PhibinEE == nPhiBinsEE ) PhibinEE = 0;
      
      //      int regionId = templIndexEE(typeEE,eleEta,charge,R9);
      //      if( regionId == -1 ) continue;

      if (ele1_iz==1)  index = 0;
      if (ele1_iz==-1) index = 1;
      
      (h_EoC_EE.at(PhibinEE)).at(index) -> Fill(pTK,ww);  // This is DATA
      (h_Phi_EE.at(PhibinEE)).at(index) -> Fill(elePhi); 
      h_phi_data_EE[index] -> Fill(elePhi);
    }
    
    
        
    // DATA  - ele2
    if(ele2_iz!=0 && R92>=0.8)
    {     
      if( ele2_iz ==  1 ) mapConversionEEp -> SetBinContent(ele2_ix,ele2_iy,scEta2);
      if( ele2_iz == -1 ) mapConversionEEm -> SetBinContent(ele2_ix,ele2_iy,scEta2);
      
      int PhibinEE = hPhiBinEE->FindBin(elePhi2) - 1;
      if( PhibinEE == nPhiBinsEE ) PhibinEE = 0;
      
      //      int regionId = templIndexEE(typeEE,eleEta2,charge2,R92);
      //      if( regionId == -1 ) continue;

      if (ele2_iz==1)  index = 0;
      if (ele2_iz==-1) index = 1;
      
      (h_EoC_EE.at(PhibinEE)).at(index) -> Fill(pTK2,ww);  // This is DATA
      (h_Phi_EE.at(PhibinEE)).at(index) -> Fill(elePhi2); 
      h_phi_data_EE[index] ->Fill(elePhi2);
    }
    
    h_et_data ->Fill(scEt);
    h_et_data ->Fill(scEt2);
  }
  
  std::cout << "End loop: Analyze events " << endl; 
  
  
  
  
  
  
  //----------------
  // Initializations
  
  // initialize TGraphs
  TFile* o = new TFile((outputFile+"_"+typeEB+"_"+typeEE+".root").c_str(),"RECREATE");
    
  TGraphErrors** g_EoC_EE = new TGraphErrors*[nRegionsEE];
 
  for(int j = 0; j < nRegionsEE; ++j)
  {
    g_EoC_EE[j]= new TGraphErrors();
  }
  
  
  //-------------------
  // Template Fit in EE
  
  if( typeEE != "none" )
  {
    float pVector[nPhiBinsEE][2];
    float pVectorErr[nPhiBinsEE][2];

    for(int i = 0; i < nPhiBinsEE; ++i)
    {
      for(int j = 0; j < nRegionsEE; ++j)
      {
        float flPhi = hPhiBinEE->GetXaxis()->GetBinCenter(i);
        
        (h_EoC_EE.at(i)).at(j) -> Rebin(rebinEE);    
        
        
        // define the fitting function
        // N.B. [0] * ( [1] * f( [1]*(x-[2]) ) )
        
        char funcName[50];
        sprintf(funcName,"f_EoP_%d_%d_Ref_%d_%d_EE",i,j,refIdEE.at(i),j);
        
        sprintf(funcName,"f_EoC_%d_%d_Ref_%d_%d_EE",i,j,refIdEE.at(i),j);
	//        (f_EoC_EE.at(i)).push_back( new TF1(funcName, (templateHistoFuncEE.at(refIdEE.at(i))).at(j), 0.7, 1.1, 3, "histoFunc") );
        (f_EoC_EE.at(i)).push_back( new TF1(funcName, "gaus", 0., 100 ) );
        
        (f_EoC_EE.at(i)).at(j) -> SetParName(0,"Norm"); 
        (f_EoC_EE.at(i)).at(j) -> SetParName(1,"Scale factor"); 
        
        (f_EoC_EE.at(i)).at(j) -> SetLineWidth(1); 
        (f_EoC_EE.at(i)).at(j) -> SetLineColor(kGreen+2); 
        (f_EoC_EE.at(i)).at(j) -> SetNpx(10000);
        
	//        (f_EoC_EE.at(i)).at(j) -> FixParameter(0, xNorm);
	//        (f_EoC_EE.at(i)).at(j) -> FixParameter(2, 0.);
        
        
        std::cout << "***** Fitting DATA EE " << flPhi << " (" << i << "," << j << "):   ";
        TFitResultPtr rp;
        int fStatus;

        for(int trial = 0; trial < 10; ++trial)
        {
	  //          (f_EoC_EE.at(i)).at(j) -> SetParameter(1, 0.99);
          rp = (h_EoC_EE.at(i)).at(j) -> Fit(funcName, "QR+");
          if( fStatus !=4 && (f_EoC_EE.at(i)).at(j)->GetParError(1) != 0. )
          {
            std::cout << "fit OK    ";
            
            double k = (f_EoC_EE.at(i)).at(j)->GetParameter(1);
            double eee = (f_EoC_EE.at(i)).at(j)->GetParError(1);
	    pVector[i][j] = k;

	    //            g_EoC_EE[j] -> SetPoint(i, flPhi, k);
	    //            g_EoC_EE[j] -> SetPointError(i, 0., eee);
            
            break;
          }
          else if( trial == 9 )
          {
	    pVector[i][j]=-1;
          }
        }
	pVector[i][j] = (h_EoC_EE.at(i)).at(j) -> GetMean();
	std::cout<<"mean p of ring "<<i<<" region: "<<j<<" is "<<pVector[i][j]<<" entries: "<<(h_EoC_EE.at(i)).at(j) -> GetEntries()<<std::endl;
	pVectorErr[i][j] = ((h_EoC_EE.at(i)).at(j) -> GetRMS())/sqrt((h_EoC_EE.at(i)).at(j) -> GetEntries());
      }
      
      std::cout << std::endl;
    }
    getchar();

	///////
    float pMean[2];
    float pMeanErr[2];

    for(int jc = 0; jc < nRegionsEE; ++jc)
      {	
	float sum=0.;
	float sumErr=0.;
	int n=0;

	for (int c=0; c<nPhiBinsEE; c++)
	  {
	    if (pVector[c][jc]==-1) continue;
	    sum+=pVector[c][jc];
	    sumErr+=(1/(pVectorErr[c][jc]*pVectorErr[c][jc]));
	    n++;  
	  }
	pMean[jc] = sum/(float)n;
	pMeanErr[jc] = sqrt(1/sumErr);

	std::cout<<"pMEan: "<<pMean[jc]<<std::endl;
	std::cout<<"pMeanErr: "<<pMeanErr[jc]<<std::endl;
	getchar();
      }


    for(int jc = 0; jc < nRegionsEE; ++jc)
      {	
	for (int c=0; c<nPhiBinsEE; c++)
	  {
	    float flPhi = hPhiBinEE->GetXaxis()->GetBinCenter(c);
	    g_EoC_EE[jc] -> SetPoint(c,c,pVector[c][jc]/pMean[jc]);
	    float err=(pVectorErr[c][jc]/pMean[jc])*(pVectorErr[c][jc]/pMean[jc])+(pVector[c][jc]/(pMean[jc]*pMean[jc])*(pMeanErr[jc]*pMeanErr[jc]))*(pVector[c][jc]/(pMean[jc]*pMean[jc])*(pMeanErr[jc]*pMeanErr[jc]));
	    g_EoC_EE[jc] -> SetPointError(c,0,err);
	    std::cout<<flPhi<<" "<<pVector[c][jc]/pMean[jc]<<" "<<err<<std::endl;
	  }
      }  
  ////////        


  }
  else
  {
    for(int i = 0; i < nPhiBinsEE; ++i)
    {  
      for(int j = 0; j < nRegionsEE; ++j)
      {
        float flPhi = hPhiBinEE->GetXaxis()->GetBinCenter(i+1);
        g_EoC_EE[j] -> SetPoint(i, flPhi, 1.);
      }
    }
  }
  
  
  
  
  
  
  //-------
  // Output
   
  o -> cd();
  
  
  for(int j = 0; j < nRegionsEE; ++j)
  {
    TString Name;
    //Name = Form("g_EoP_EE_%d",j);
    //if(g_EoP_EE[j]->GetN()!=0) g_EoP_EE[j] -> Write(Name);
    Name = Form("g_EoC_EE_%d",j);
    if(g_EoC_EE[j]->GetN()!=0) g_EoC_EE[j] -> Write(Name);
    //Name = Form("g_Rat_EE_%d",j);
    //if(g_Rat_EE[j]->GetN()!=0) g_Rat_EE[j] -> Write(Name);
  }
    
  for(int j =0; j< nRegionsEE; ++j)
  {
    if( h_phi_data_EE[j] -> GetEntries() !=0 ) h_phi_data_EE[j] -> Write();
  }
  
  h_et_data->Write();
  
  mapConversionEEp -> Write();
  mapConversionEEm -> Write();
  
  o -> Close();
  
  
  
  return 0;
}
