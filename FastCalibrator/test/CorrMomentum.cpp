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


//**************  MAIN PROGRAM **************************************************************
int main(int argc, char** argv)
{
  // Acquisition from cfg file
  if(argc != 2)
  {
    std::cerr << ">>>>> CorrMomentum.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return 1;
  }
  
  parseConfigFile(argv[1]);
  
  std::string TreeName    = gConfigParser -> readStringOption("Input::TreeName");
  std::string infileDATA  = gConfigParser -> readStringOption("Input::infileDATA");
  std::string infileMC    = gConfigParser -> readStringOption("Input::infileMC");
  
  std::string typeEE = gConfigParser -> readStringOption("Input::typeEE");
  int  nRegionsEE = gConfigParser -> readIntOption("Input::nRegionsEE");
  int  nPhiBinsEE = gConfigParser -> readIntOption("Input::nPhiBinsEE");
  int  nPhiBinsTempEE = gConfigParser -> readIntOption("Input::nPhiBinsTempEE");
  int  nEtaBinsEE = gConfigParser -> readIntOption("Input::nEtaBinsEE");
  int  nEtaBinsTempEE = gConfigParser -> readIntOption("Input::nEtaBinsTempEE");
  int  rebinEE = gConfigParser -> readIntOption("Input::rebinEE");
  std::string outputFile = gConfigParser -> readStringOption("Output::outputFile");
  int nEntriesMC = gConfigParser -> readIntOption("Input::nEntriesMC");
  int nEntriesData = gConfigParser -> readIntOption("Input::nEntriesData");
  
  //  int nRegionsEE = GetNRegionsEE(typeEE);

  std::cout<<"REGIONI: "<<nRegionsEE<<std::endl;

  cout <<" Basic Configuration " <<endl;
  cout <<" Tree Name = "<<TreeName<<endl;
  cout <<" infileDATA = "<<infileDATA<<endl;
  cout <<" infileMC = "<<infileMC<<endl;
  cout <<" nRegionsEE = "<<nRegionsEE<<endl;
  cout <<" nPhiBinsEE = "<<nPhiBinsEE<<endl;
  cout <<" nPhiBinsTempEE = "<<nPhiBinsTempEE<<endl;
  cout <<" nEtaBinsEE = "<<nEtaBinsEE<<endl;
  cout <<" nEtaBinsTempEE = "<<nEtaBinsTempEE<<endl;
  cout <<" rebinEE = "<<rebinEE<<endl;
  
  cout << "Momentum correction "<< endl;
  
  
  
  //---- variables for selection
  float etaMax  = 2.5;
  float eta2Max = 2.5;

  float etaMin = 1.479;

  //  nEtaBinsEE=5;

  //histos to get the bin in phi given the electron phi
  TH1F* hPhiBinEE = new TH1F("hphiEE","hphiEE",nPhiBinsEE, -1.*TMath::Pi(),1.*TMath::Pi());
  TH1F* hEtaBinEE = new TH1F("hetaEE","hetaEE",nEtaBinsEE, 1.4,2.5);
  
  //----- NTUPLES--------------------
  TChain *ntu_DA = new TChain(TreeName.c_str());  
  if(!FillChain(*ntu_DA, infileDATA.c_str())) return 1;

  TChain *ntu_MC = new TChain(TreeName.c_str());  
  if(!FillChain(*ntu_MC, infileMC.c_str())) return 1;
  
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
  float pEle,pEle2;
  
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
  ntu_DA->SetBranchStatus("ele1_p", 1);            ntu_DA->SetBranchAddress("ele1_p", &pEle);
  ntu_DA->SetBranchStatus("ele2_p", 1);            ntu_DA->SetBranchAddress("ele2_p", &pEle2);
  ntu_DA->SetBranchStatus("ele1_seedIphi", 1);       ntu_DA->SetBranchAddress("ele1_seedIphi", &iphiSeed);
  ntu_DA->SetBranchStatus("ele2_seedIphi", 1);       ntu_DA->SetBranchAddress("ele2_seedIphi", &iphiSeed2);
  ntu_DA->SetBranchStatus("ele1_seedIx", 1);         ntu_DA->SetBranchAddress("ele1_seedIx", &ele1_ix);
  ntu_DA->SetBranchStatus("ele2_seedIx", 1);         ntu_DA->SetBranchAddress("ele2_seedIx", &ele2_ix);
  ntu_DA->SetBranchStatus("ele1_seedIy", 1);         ntu_DA->SetBranchAddress("ele1_seedIy", &ele1_iy);
  ntu_DA->SetBranchStatus("ele2_seedIy", 1);         ntu_DA->SetBranchAddress("ele2_seedIy", &ele2_iy);
  ntu_DA->SetBranchStatus("ele1_seedZside", 1);      ntu_DA->SetBranchAddress("ele1_seedZside", &ele1_iz);
  ntu_DA->SetBranchStatus("ele2_seedZside", 1);      ntu_DA->SetBranchAddress("ele2_seedZside", &ele2_iz);

  // Set branch addresses for MC
  ntu_MC->SetBranchStatus("*",0);
  ntu_MC->SetBranchStatus("isW", 1);                 ntu_MC->SetBranchAddress("isW", &isW);
  ntu_MC->SetBranchStatus("ele1_eta", 1);            ntu_MC->SetBranchAddress("ele1_eta", &eleEta);
  ntu_MC->SetBranchStatus("ele2_eta", 1);            ntu_MC->SetBranchAddress("ele2_eta", &eleEta2);
  ntu_MC->SetBranchStatus("ele1_phi", 1);            ntu_MC->SetBranchAddress("ele1_phi", &elePhi);
  ntu_MC->SetBranchStatus("ele2_phi", 1);            ntu_MC->SetBranchAddress("ele2_phi", &elePhi2);
  ntu_MC->SetBranchStatus("ele1_scEta", 1);          ntu_MC->SetBranchAddress("ele1_scEta", &scEta);
  ntu_MC->SetBranchStatus("ele2_scEta", 1);          ntu_MC->SetBranchAddress("ele2_scEta", &scEta2);
  ntu_MC->SetBranchStatus("ele1_scPhi", 1);          ntu_MC->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_MC->SetBranchStatus("ele2_scPhi", 1);          ntu_MC->SetBranchAddress("ele2_scPhi", &scPhi2);
  ntu_MC->SetBranchStatus("ele1_scE", 1);            ntu_MC->SetBranchAddress("ele1_scE", &scEne);
  ntu_MC->SetBranchStatus("ele2_scE", 1);            ntu_MC->SetBranchAddress("ele2_scE", &scEne2);
  ntu_MC->SetBranchStatus("ele1_scE_regression", 1); ntu_MC->SetBranchAddress("ele1_scE_regression", &scEneReg);
  ntu_MC->SetBranchStatus("ele2_scE_regression", 1); ntu_MC->SetBranchAddress("ele2_scE_regression", &scEneReg2);
  ntu_MC->SetBranchStatus("ele1_scEt", 1);           ntu_MC->SetBranchAddress("ele1_scEt",&scEt);
  ntu_MC->SetBranchStatus("ele2_scEt", 1);           ntu_MC->SetBranchAddress("ele2_scEt",&scEt2);
  ntu_MC->SetBranchStatus("ele1_scERaw", 1);         ntu_MC->SetBranchAddress("ele1_scERaw", &scERaw);
  ntu_MC->SetBranchStatus("ele2_scERaw", 1);         ntu_MC->SetBranchAddress("ele2_scERaw", &scERaw2);
  ntu_MC->SetBranchStatus("ele1_e3x3", 1);           ntu_MC->SetBranchAddress("ele1_e3x3", &e3x3);
  ntu_MC->SetBranchStatus("ele2_e3x3", 1);           ntu_MC->SetBranchAddress("ele2_e3x3", &e3x32);
  ntu_MC->SetBranchStatus("ele1ele2_scM", 1);        ntu_MC->SetBranchAddress("ele1ele2_scM", &mZ);
  ntu_MC->SetBranchStatus("ele1_charge", 1);         ntu_MC->SetBranchAddress("ele1_charge", &charge);
  ntu_MC->SetBranchStatus("ele2_charge", 1);         ntu_MC->SetBranchAddress("ele2_charge", &charge2);
  ntu_MC->SetBranchStatus("ele1_tkP", 1);            ntu_MC->SetBranchAddress("ele1_tkP", &pTK);
  ntu_MC->SetBranchStatus("ele2_tkP", 1);            ntu_MC->SetBranchAddress("ele2_tkP", &pTK2);
  ntu_MC->SetBranchStatus("ele1_p", 1);            ntu_MC->SetBranchAddress("ele1_p", &pEle);
  ntu_MC->SetBranchStatus("ele2_p", 1);            ntu_MC->SetBranchAddress("ele2_p", &pEle2);
  ntu_MC->SetBranchStatus("ele1_seedIphi", 1);       ntu_MC->SetBranchAddress("ele1_seedIphi", &iphiSeed);
  ntu_MC->SetBranchStatus("ele2_seedIphi", 1);       ntu_MC->SetBranchAddress("ele2_seedIphi", &iphiSeed2);
  ntu_MC->SetBranchStatus("ele1_seedIx", 1);         ntu_MC->SetBranchAddress("ele1_seedIx", &ele1_ix);
  ntu_MC->SetBranchStatus("ele2_seedIx", 1);         ntu_MC->SetBranchAddress("ele2_seedIx", &ele2_ix);
  ntu_MC->SetBranchStatus("ele1_seedIy", 1);         ntu_MC->SetBranchAddress("ele1_seedIy", &ele1_iy);
  ntu_MC->SetBranchStatus("ele2_seedIy", 1);         ntu_MC->SetBranchAddress("ele2_seedIy", &ele2_iy);
  ntu_MC->SetBranchStatus("ele1_seedZside", 1);      ntu_MC->SetBranchAddress("ele1_seedZside", &ele1_iz);
  ntu_MC->SetBranchStatus("ele2_seedZside", 1);      ntu_MC->SetBranchAddress("ele2_seedZside", &ele2_iz);

  std::cout<<"QUI"<<std::endl;  

  // histogram definition in EE and fit functions                                                                      
  //  std::vector<std::vector<TH1F*> > h_Phi_EE(nPhiBinsEE); // used to map iEta (as defined for Barrel and Endcap geom) into eta          
  //  std::vector<std::vector<TH1F*> > h_Eta_EE(nEtaBinsEE); // used to map iEta (as defined for Barrel and Endcap geom) into eta          
  TH1F* h_pData_EE[nPhiBinsEE][nEtaBinsEE][nRegionsEE];
  //  std::vector<std::vector<std::vector<TH1F*> > > h_pData_EE(nPhiBinsEE);
  TF1* f_pData_EE[nPhiBinsEE][nEtaBinsEE][nRegionsEE];
  
  TH1F* histoPull_EE[nEtaBinsEE][nRegionsEE];
  
  //  nRegionsEE=2; //EE- and EE+
  //  std::vector<TH1F* > vect1(nEtaBinsEE);

  // Initializate histos in EE
  std::cout << ">>> Initialize EE histos" << std::endl;
    //    std::vector<std::vector<TH1F*> >tempVect(nEtaBinsEE);
  for (int k=0; k<nEtaBinsEE; ++k)
    {
      TString histoName;
      TH1F* temp;
      //      std::cout<<i<<" "<<k<<" "<<j<<std::endl;
      for(int j = 0; j < nRegionsEE; ++j)
	{
	  
	  for(int i = 0; i < nPhiBinsEE; ++i)
	    {
	      
	      histoName= Form("EE_pData_%d_%d_%d", i,k,j);
	      temp = new TH1F (histoName, histoName, 50, 0., 500.);
	      temp->Sumw2();
	      temp->SetFillColor(kGreen+2);
	      temp->SetLineColor(kGreen+2);
	      temp->SetFillStyle(3004);
	      h_pData_EE[i][k][j] = temp;
	      //      (tempVect.at(k)).push_back(temp);
	      
	      //      histoName=Form("EE_Phi_%d_%d_%d", i,k,j);
	      //      temp = new TH1F(histoName, histoName, 360, 0., 360.); 
	      //      (h_Phi_EE.at(i)).push_back(temp); 
	    }
	  
	  //    std::cout<<"qui?"<<std::endl;
	  histoName=Form("histoPull_%d_%d", k,j);
	  temp = new TH1F(histoName, histoName, 100, -10, 10); 
	  histoPull_EE[k][j]=temp;
	  //      (h_Eta_EE.at(k)).push_back(temp); 
	  //    std::cout<<"qui2?"<<std::endl;
	}
      //    (h_pData_EE).push_back(tempVect);
      
    }


 // Template in EE
  //  std::vector<std::vector<TH1F*> > h_template_EE(nPhiBinsTempEE);
  TH1F* h_template_EE[nPhiBinsTempEE][nEtaBinsTempEE][nRegionsEE];
    
  std::cout << ">>> Initialize EE template" << std::endl;
  for(int mod = 0; mod < nPhiBinsTempEE; ++mod)
  {
  for(int k = 0; k < nEtaBinsEE; ++k)
  {
    for(int j = 0; j < nRegionsEE; ++j)
    {
      TString histoName;
      histoName=Form("EE_template_%d_%d_%d",mod,k,j);
      TH1F* temp = new TH1F(histoName,"",50,0.,500.);
      h_template_EE[mod][k][j] = temp;
      //      std::cout<<"mah: "<<mod<<" "<<j<<std::endl;
    }
  }  
  }  
 
  TH1F** h_phi_data_EE = new TH1F*[nRegionsEE];
  TH1F** h_eta_data_EE = new TH1F*[nRegionsEE];
  TH1F** h_phi_mc_EE = new TH1F*[nRegionsEE];

  for(int index = 0; index < nRegionsEE; ++index)
  {
    TString name;
    name=Form("EE_h_phi_data_%d",index);
    h_phi_data_EE[index] = new TH1F(name,"h_phi_data",100,-TMath::Pi(),TMath::Pi());
    name=Form("EE_h_phi_mc_%d",index);
    h_phi_mc_EE[index] =  new TH1F(name,"h_phi_mc",100,-TMath::Pi(),TMath::Pi());
    name=Form("EE_h_eta_data_%d",index);
    h_eta_data_EE[index] = new TH1F(name,"h_eta_data",100,1.479,2.5);
  }
  
  TH1F* h_et_data = new TH1F("h_et_data","h_et_data",100,0.,100.);
  TH1F* h_et_mc   = new TH1F("h_et_mc",  "h_et_mc",  100,0.,100.);
    
  
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
    //    std::cout<<iphi<<" "<<modPhi<<std::endl;
  }
    
 //**************************** loop on MC
  
  std::cout << "first loop: fill template histo" << endl; 
  
  //  for(int entry = 0; entry < ntu_MC->GetEntries(); ++entry)
  if (nEntriesMC<0) nEntriesMC = ntu_MC->GetEntries();
  for(int entry = 0; entry < nEntriesMC; ++entry)
  {
    if( entry%10000 == 0 ) 
      std::cout << "reading saved entry " << entry << "\r" << std::flush;
    
    ntu_MC->GetEntry(entry);
    //    std::cout<<fabs(scEta)<<" "<<fabs(scEta2)<<" "<<scEt<<" "<<scEt2<<std::endl;
    
    //    if( isW == 1 )               continue;
    if( fabs(scEta)  > etaMax )  continue;
    if( fabs(scEta2) > eta2Max ) continue;

    if(fabs(eleEta)<etaMin) continue;
    if(fabs(eleEta2)<etaMin) continue;

    if( scEt  < 20. ) continue;
    if( scEt2 < 20. ) continue;    
    
    R9  = e3x3  / scERaw;
    R92 = e3x32 / scERaw2;

    float ww = 1.;
    int index=0;

   //   std::cout<<ele1_iz<<std::endl;
    // MC - ENDCAP - ele1
    if (ele1_iz!=0)
    {
      if( ele1_iz ==  1 ) mapConversionEEp -> SetBinContent(ele1_ix,ele1_iy,scEta);
      if( ele1_iz == -1 ) mapConversionEEm -> SetBinContent(ele1_ix,ele1_iy,scEta);
      
      int iphi = eRings->GetEndcapIphi(ele1_ix,ele1_iy,ele1_iz);
      
      // fill MC templates
      int modPhi = int (iphi/(360./nPhiBinsTempEE));
      if( modPhi == nPhiBinsTempEE ) modPhi = 0;

      int EtabinEE = hEtaBinEE->FindBin(fabs(eleEta)) - 1;
      if( EtabinEE == nEtaBinsEE ) EtabinEE = 0;
      
//      int regionId =  templIndexEE(typeEE,eleEta1,charge2,R92);
//      if(regionId == -1) continue;
 
      if (ele1_iz==1)  index = 0;
      if (ele1_iz==-1) index = 1;
     
      h_template_EE[modPhi][EtabinEE][index] ->  Fill(pEle,ww);
      
      // fill MC histos in eta bins
      int PhibinEE = hPhiBinEE->FindBin(scPhi) - 1;
      if(PhibinEE==nPhiBinsEE) PhibinEE = 0;
      
      //      std::cout<<"MC: fill with "<<pTK<<" "<<ww<<std::endl;
      //      (h_pMC_EE.at(PhibinEE)).at(index) -> Fill(pTK,ww);  // This is MC
      h_phi_mc_EE[index]->Fill(scPhi,ww);
    }

    // MC - ENDCAP - ele2
    if (ele2_iz!=0)
    {
      if( ele2_iz ==  1 ) mapConversionEEp -> SetBinContent(ele2_ix,ele2_iy,scEta2);
      if( ele2_iz == -1 ) mapConversionEEm -> SetBinContent(ele2_ix,ele2_iy,scEta2);
      
      int iphi = eRings->GetEndcapIphi(ele2_ix,ele2_iy,ele2_iz);
      
      // fill MC templates
      int modPhi = int (iphi/(360./nPhiBinsTempEE));
      if( modPhi == nPhiBinsTempEE ) modPhi = 0;

      int EtabinEE = hEtaBinEE->FindBin(fabs(eleEta2)) - 1;
      if( EtabinEE == nEtaBinsEE ) EtabinEE = 0;
      
//      int regionId =  templIndexEE(typeEE,eleEta1,charge2,R92);
//      if(regionId == -1) continue;
 
      if (ele2_iz==1)  index = 0;
      if (ele2_iz==-1) index = 1;

      h_template_EE[modPhi][EtabinEE][index] ->  Fill(pEle2,ww);
      
      // fill MC histos in eta bins
      int PhibinEE = hPhiBinEE->FindBin(scPhi2) - 1;
      if(PhibinEE==nPhiBinsEE) PhibinEE = 0;
      
      //      (h_pMC_EE.at(PhibinEE)).at(index) -> Fill(pTK2,ww);  // This is MC
      h_phi_mc_EE[index]->Fill(scPhi2,ww);
    }
    
    h_et_mc ->Fill(scEt, ww);
    h_et_mc ->Fill(scEt2,ww);
  }
  


  
  //**************************** loop on DATA
  
  std::cout << "Loop in Data events " << endl; 

  /*  for (int i=0; i<nPhiBinsEE; ++i) {
    for (int k=0; k<nEtaBinsEE; ++k) {
      for (int j=0; j<nRegionsEE; ++j) {
	std::cout<<i<<" "<<j<<" "<<k<<std::endl;
	((h_pData_EE.at(i)).at(k)).at(j) -> Fill(0);  // This is DATA
      }
    }
  }
  */
  //  for(int entry = 0; entry < ntu_DA->GetEntries(); ++entry)
  if (nEntriesData<0) nEntriesData = ntu_DA->GetEntries();
  for(int entry = 0; entry < nEntriesData; ++entry)
  {
    if( entry%10000 == 0 ) 
      std::cout << "reading saved entry " << entry << "\r" << std::flush;
      
    ntu_DA->GetEntry(entry);
    //    if( isW == 1 )               continue;
    if( fabs(scEta)  > etaMax )  continue;
    if( fabs(scEta2) > eta2Max ) continue;

    if(fabs(eleEta)<etaMin) continue;
    if(fabs(eleEta2)<etaMin) continue;

    if( scEt  < 20. ) continue;
    if( scEt2 < 20. ) continue;    
    
    R9  = e3x3  / scERaw;
    R92 = e3x32 / scERaw2;

    float ww = 1.;
    int index=0;

    // DATA - ENDCAP - ele1
    if (ele1_iz!=0)
    {
      if( ele1_iz ==  1 ) mapConversionEEp -> SetBinContent(ele1_ix,ele1_iy,scEta);
      if( ele1_iz == -1 ) mapConversionEEm -> SetBinContent(ele1_ix,ele1_iy,scEta);
      
      int PhibinEE = hPhiBinEE->FindBin(elePhi) - 1;
      if( PhibinEE == nPhiBinsEE ) PhibinEE = 0;

      int EtabinEE = hEtaBinEE->FindBin(fabs(eleEta)) - 1;
      if( EtabinEE == nEtaBinsEE ) EtabinEE = 0;
      
      //      int regionId = templIndexEE(typeEE,eleEta,charge,R9);
      //      if( regionId == -1 ) continue;

      if (ele1_iz==1)  index = 0;
      if (ele1_iz==-1) index = 1;
      
      h_pData_EE[PhibinEE][EtabinEE][index] -> Fill(pTK,ww);  // This is DATA
      //      (h_Phi_EE.at(PhibinEE)).at(index) -> Fill(elePhi); 

      h_phi_data_EE[index] -> Fill(elePhi);
      h_eta_data_EE[index] -> Fill(fabs(eleEta));
    }
    
        
    // DATA  - ele2
    if(ele2_iz!=0)
    {     
      if( ele2_iz ==  1 ) mapConversionEEp -> SetBinContent(ele2_ix,ele2_iy,scEta2);
      if( ele2_iz == -1 ) mapConversionEEm -> SetBinContent(ele2_ix,ele2_iy,scEta2);
      
      int PhibinEE = hPhiBinEE->FindBin(elePhi2) - 1;
      if( PhibinEE == nPhiBinsEE ) PhibinEE = 0;

      int EtabinEE = hEtaBinEE->FindBin(fabs(eleEta2)) - 1;
      if( EtabinEE == nEtaBinsEE ) EtabinEE = 0;
      
      //      int regionId = templIndexEE(typeEE,eleEta2,charge2,R92);
      //      if( regionId == -1 ) continue;

      if (ele2_iz==1)  index = 0;
      if (ele2_iz==-1) index = 1;
      //      std::cout<<"qui2 "<<EtabinEE<<" "<<fabs(eleEta2)<<std::endl;
      
      h_pData_EE[PhibinEE][EtabinEE][index] -> Fill(pTK2,ww);  // This is DATA
      //      (h_Phi_EE.at(PhibinEE)).at(index) -> Fill(elePhi2); 
      h_phi_data_EE[index] ->Fill(elePhi2);
      h_eta_data_EE[index] -> Fill(fabs(eleEta2));
    }
    
    h_et_data ->Fill(scEt);
    h_et_data ->Fill(scEt2);
  }
  
  std::cout << "End loop: Analyze events " << endl; 
  
  
  
  
  
  
  //----------------
  // Initializations
  
  // initialize TGraphs
  TFile* o = new TFile((outputFile+"_"+typeEE+".root").c_str(),"RECREATE");
    
  TGraphErrors* g_pData_EE[nEtaBinsEE][nRegionsEE];// = new TGraphErrors**[nEtaBinsEE][nRegionsEE];
  TGraphErrors* g_pAbs_EE[nEtaBinsEE][nRegionsEE];// = new TGraphErrors**[nEtaBinsEE][nRegionsEE];

  for (int a=0; a<nEtaBinsEE; ++a)
    {
  for(int j = 0; j < nRegionsEE; ++j)
  {
    g_pData_EE[a][j]= new TGraphErrors();
    g_pAbs_EE[a][j]= new TGraphErrors();
  }
    } 
  
  // initialize template functions  
  //  std::vector<std::vector<histoFunc*> > templateHistoFuncEE(nPhiBinsTempEE);
  histoFunc* templateHistoFuncEE[nPhiBinsTempEE][nEtaBinsEE][nRegionsEE];

  for(int mod = 0; mod < nPhiBinsTempEE; ++mod)
  {
  for(int k = 0; k < nEtaBinsEE; ++k)
  {
    for(int j = 0; j < nRegionsEE; ++j)
    {
      //      h_template_EE[mod][k][j] -> Rebin(rebinEE);
      templateHistoFuncEE[mod][k][j] = new histoFunc(h_template_EE[mod][k][j]);
    }
  }
  }

  //-------------------
  // Template Fit in EE
  
  if( typeEE != "none" )
  {
    float pVector[nPhiBinsEE][nEtaBinsEE][2];
    float pVectorErr[nPhiBinsEE][nEtaBinsEE][2];

    for(int i = 0; i < nPhiBinsEE; ++i)
    {
    for(int k = 0; k < nEtaBinsEE; ++k)
    {
      for(int j = 0; j < nRegionsEE; ++j)
      {
        float flPhi = hPhiBinEE->GetXaxis()->GetBinCenter(i);
        
	//        (h_pMC_EE.at(i)).at(j) -> Rebin(rebinEE);
	//        h_pData_EE[i][k][j] -> Rebin(rebinEE);    
        
        
        // define the fitting function
        // N.B. [0] * ( [1] * f( [1]*(x-[2]) ) )
        
        char funcName[50];        
        
        sprintf(funcName,"f_pData_%d_%d_%d_Ref_%d_%d_%d_EE",i,k,j,refIdEE.at(i),k,j);
	TF1* temp;
	temp = new TF1(funcName, templateHistoFuncEE[refIdEE.at(i)][k][j], 0., 500, 3, "histoFunc");
        f_pData_EE[i][k][j] =  temp;
        
        f_pData_EE[i][k][j] -> SetParName(0,"Norm"); 
        f_pData_EE[i][k][j] -> SetParName(1,"Scale factor"); 
        
        f_pData_EE[i][k][j] -> SetLineWidth(1); 
        f_pData_EE[i][k][j] -> SetLineColor(kGreen+2); 
        f_pData_EE[i][k][j] -> SetNpx(10000);
        
	//	f_pData_EE[i][k][j] -> SetParameter(0, xNorm);
	f_pData_EE[i][k][j] -> SetParameter(0, 1.);
	f_pData_EE[i][k][j] -> SetParameter(1, 1);

	float shift=8.;
	f_pData_EE[i][k][j] -> SetParameter(2, shift);
                               
        std::cout << "***** Fitting DATA EE " << flPhi << " (" << i << "," << j << "):   ";

	TFitResultPtr rp;
	int fStatus; 

        for(int trial = 0; trial < 100; ++trial)
        {
          rp = h_pData_EE[i][k][j] -> Fit(funcName, "QR+");
          fStatus = rp;

          if( fStatus !=4 && f_pData_EE[i][k][j]->GetParError(1) != 0. && f_pData_EE[i][k][j] -> GetMaximumX(0.,500.)>30. )
          {
            std::cout << "fit OK    ";
            
            double coeff = f_pData_EE[i][k][j]->GetParameter(1);
            double eee = f_pData_EE[i][k][j]->GetParError(1);
	    pVector[i][k][j] = coeff;

            break;
          }
          else if( trial %10 == 0 )
          {
	    pVector[i][k][j]=-1;
	    //	    std::cout<<" BAD FIT. Make another try with different parameters.. "<<std::endl;
	    shift-=5;
	    f_pData_EE[i][k][j] -> SetParameter(2, shift);
	    if (shift==-25)  shift+=50.;
	    //	    trial = 0;
	    //	    getchar();
          }
        }

	std::cout<<"media del bin "<<i<<" : "<<f_pData_EE[i][k][j] -> GetMaximumX(0.,500.)<<std::endl;
	//(f_pData_EE.at(0)).at(0)->GetParameter(2)*(f_pData_EE.at(0)).at(0)->GetParameter(1)+(h_template_EE.at(0)).at(j)->GetMean()<<std::endl;

	if( fStatus !=4 && f_pData_EE[i][k][j]->GetParError(1) != 0. && f_pData_EE[i][k][j] -> GetMaximumX(0.,500.)>30. ) {
	  pVector[i][k][j] = f_pData_EE[i][k][j] -> GetMaximumX(0.,500.);
	  pVectorErr[i][k][j] = (h_pData_EE[i][k][j] -> GetRMS())/sqrt(h_pData_EE[i][k][j] -> GetEntries());
	}
	else {
	  std::cout<<"BAD FIT!!!"<<std::endl;
	  pVector[i][k][j] = -1;  //if fit fails
	  pVectorErr[i][k][j] = 0;
	}
      }
    
  }
      
      std::cout << std::endl;
    }

	///////
    float pMean[nEtaBinsEE][nRegionsEE];
    float pMeanErr2[nEtaBinsEE][nRegionsEE];

    for(int jc = 0; jc < nRegionsEE; ++jc)
      {	
	float sum=0.;
	float sumErr2=0.;
	int n=0;
	for (int a=0; a<nEtaBinsEE; a++)
	  {
	sum=0.;
	sumErr2=0.;
	n=0;

	for (int c=0; c<nPhiBinsEE; c++)
	  {
	    if (pVector[c][a][jc]==-1) continue;
	    sum+=pVector[c][a][jc];
	    sumErr2+=(1/(pVectorErr[c][a][jc]*pVectorErr[c][a][jc]));
	    n++;  
	  }
	pMean[a][jc] = sum/(float)n;
	pMeanErr2[a][jc] = sqrt(1/sumErr2);
	std::cout<<"pMEan: "<<pMean[a][jc]<<std::endl;
	std::cout<<"pMeanErr2: "<<pMeanErr2[a][jc]<<std::endl;

	  }
      }


    for(int jc = 0; jc < nRegionsEE; ++jc)
      {	
	for (int a=0; a<nEtaBinsEE; a++)
	  {
	for (int c=0; c<nPhiBinsEE; c++)
	  {
	    float flPhi = hPhiBinEE->GetXaxis()->GetBinCenter(c);
	    if (pVector[c][a][jc]==-1) {
	      pVector[c][a][jc]=pMean[a][jc]; //if fit has failed, fill with the mean value
	      std::cout<<"be careful!! ("<<c<<","<<a<<","<<jc<<") has a bad value! Fill with default value (1)."<<std::endl;
	    }
	    if ( (pVector[c][a][jc]/pMean[a][jc])>1.1) 
	      g_pData_EE[a][jc] -> SetPoint(c,c*(int(360/nPhiBinsEE)),1.1);
	    else if ( (pVector[c][a][jc]/pMean[a][jc])<0.9 )
	      g_pData_EE[a][jc] -> SetPoint(c,c*(int(360/nPhiBinsEE)),0.9);
	    else	    
	      g_pData_EE[a][jc] -> SetPoint(c,c*(int(360/nPhiBinsEE)),pVector[c][a][jc]/pMean[a][jc]);

	    g_pAbs_EE[a][jc] -> SetPoint(c,c*(int(360/nPhiBinsEE)),pVector[c][a][jc]);
	    histoPull_EE[a][jc] -> Fill((pVector[c][a][jc]-pMean[a][jc])/pVectorErr[c][a][jc]);
            
	    //	    float err=(pVectorErr[c][a][jc]/pMean[a][jc])*(pVectorErr[c][a][jc]/pMean[a][jc])+(pVector[c][a][jc]/(pMean[a][jc]*pMean[a][jc])*(pMeanErr2[a][jc]*pMeanErr2[a][jc]))*(pVector[c][a][jc]/(pMean[a][jc]*pMean[a][jc])*(pMeanErr2[a][jc]*pMeanErr2[a][jc]));
	    float err=(pVectorErr[c][a][jc]/pMean[a][jc]);
	    g_pData_EE[a][jc] -> SetPointError(c,0,err);
	    g_pAbs_EE[a][jc] -> SetPointError(c,0,pVectorErr[c][a][jc]);
	    //	    std::cout<<flPhi<<" "<<pVector[c][a][jc]/pMean[a][jc]<<" "<<err<<std::endl;
	  }
	  }
      }  
  ////////        


  }
  else
  {
    for(int i = 0; i < nPhiBinsEE; ++i)
    {  
      for (int a=0; a<nEtaBinsEE; ++a)
	{
      for(int j = 0; j < nRegionsEE; ++j)
      {
        float flPhi = hPhiBinEE->GetXaxis()->GetBinCenter(i+1);
        g_pData_EE[a][j] -> SetPoint(i, flPhi, 1.);
      }
	}
    }
  }
  
  
  
  
  
  
  //-------
  // Output
   
  o -> cd();
  
  for (int a=0; a<nEtaBinsEE; ++a)
    {  
  for(int j = 0; j < nRegionsEE; ++j)
  {
    TString Name;
    //Name = Form("g_pMC_EE_%d",j);
    //if(g_pMC_EE[j]->GetN()!=0) g_pMC_EE[j] -> Write(Name);
    Name = Form("g_pData_EE_%d_%d",a,j);
    if(g_pData_EE[a][j]->GetN()!=0) g_pData_EE[a][j] -> Write(Name);
    Name = Form("g_pAbs_EE_%d_%d",a,j);
    if(g_pAbs_EE[a][j]->GetN()!=0) g_pAbs_EE[a][j] -> Write(Name);
    //Name = Form("g_Rat_EE_%d",j);
    //if(g_Rat_EE[j]->GetN()!=0) g_Rat_EE[j] -> Write(Name);
  }
    }
    
  for(int j =0; j< nRegionsEE; ++j)
  {
    if( h_phi_data_EE[j] -> GetEntries() !=0 ) h_phi_data_EE[j] -> Write();
  }
  
  h_et_data->Write();
  
  mapConversionEEp -> Write();
  mapConversionEEm -> Write();

  h_template_EE[0][0][0] -> Write();
  h_template_EE[0][0][1] -> Write();
  //  h_template_EE[0][1][0] -> Write();
  //  h_template_EE[0][1][1] -> Write();
  //  h_template_EE[0][2][0] -> Write();
  //  h_template_EE[0][2][1] -> Write();
  //  h_template_EE[0][4][0] -> Write();
  //  h_template_EE[0][4][1] -> Write();

  for (int k=0; k<nEtaBinsEE; ++k)
    {
      for(int j = 0; j < nRegionsEE; ++j)
	{
	  for(int i = 0; i < nPhiBinsEE; ++i)
	    {  
	      h_pData_EE[i][k][j] -> Write();
	    }	  
	}
    }

  for (int k=0; k<nEtaBinsEE; ++k)
    {
      for(int j = 0; j < nRegionsEE; ++j)
	{
	  histoPull_EE[k][j]->Write();
	}
    }

  o -> Close();
  
  
  
  return 0;
}
