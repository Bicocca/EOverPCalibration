#include "TEndcapRings.h"
#include <iostream>
#include <fstream>
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
#include "TMath.h"

using namespace std;


//**************  MAIN PROGRAM **************************************************************                                                                     
int main(int argc, char** argv)
{
  // Set style options
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(1110); 
  gStyle->SetOptFit(0); 
  gStyle->SetFitFormat("6.3g"); 
  gStyle->SetPalette(1); 
 
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.05);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetLabelSize(0.05);
  gROOT->ForceStyle();

  TFile* o = new TFile("oldMap.root","RECREATE");

  double x,y,z,IC;
  double x2,y2,z2,IC2;
  TH2F* map[2];
  TH2F* mapNorm[2];

  map[0] = new TH2F("mapEEM","mapEEM",100,0,100,100,0,100);
  map[1] = new TH2F("mapEEP","mapEEP",100,0,100,100,0,100);
  mapNorm[0] = new TH2F("mapNormEEM","mapNormEEM",100,0,100,100,0,100);
  mapNorm[1] = new TH2F("mapNormEEP","mapNormEEP",100,0,100,100,0,100);

  ifstream inputCfg ("IC_digi.txt");;
  ifstream inputCfg2 ("IC_mc.txt");;

  ofstream inputCfg3 ("IC_digi_su_reco.txt");;
  //  inputCfg.open("ICRunD_EE2.txt");

  TEndcapRings* eRings = new TEndcapRings();

  while(!inputCfg.eof())
    {
      inputCfg >> x >> y >> z >> IC;
      inputCfg2 >> x2 >> y2 >> z2 >> IC2;
	  int ring = eRings->GetEndcapRing(x,y,z);
	  //	  if (ring<2) std::cout<<x<<" "<<y<<" "<<IC<<" "<<IC2<<" "<<IC/IC2<<std::endl;

      if (z==-1)
	map[0]->Fill(x,y,IC2);      
      else if (z==1)
	map[1]->Fill(x,y,IC2);      
      inputCfg3 << x <<"\t"<< y <<"\t"<< z <<"\t"<< IC/IC2<<"\t"<<-1<<std::endl;
    }
  //  map->Draw("COLZ");
  //  getchar();

  std::map<int,std::vector<float> > sumIC;
  std::map<int,std::vector<int> > numIC;

  (sumIC[0]).assign(40,0.);
  (sumIC[1]).assign(40,0.);

  (numIC[0]).assign(40,0);
  (numIC[1]).assign(40,0);

  // mean over phi corrected skipping dead channel                                                                                                                
  for(int k = 0; k < 2; ++k) {    
    for(int ix = 1; ix <= map[k] -> GetNbinsX(); ++ix) {
      for(int iy = 1; iy <= map[k] -> GetNbinsY(); ++iy)
	{
	  //	  std::cout<<" "<<ix<<" "<<iy<<" "<<k<<std::endl;

	  int zside;
	  if (k==0) zside=-1;
	  else if (k==1) zside=1;

	  int ring = eRings->GetEndcapRing(ix,iy,k);
	  //	  if (k==0) std::cout<<ring<<std::endl;
	  if( ring == -1 ) continue;

	  int iRing = 85 + eRings -> GetEndcapRing(ix,iy,k);
	  float eta = eRings -> GetEtaFromIRing(iRing);

	  if( map[k]->GetBinContent(ix,iy) == -1 || map[k]->GetBinContent(ix,iy) == 0) continue;

	  (sumIC[k]).at(ring) += map[k]->GetBinContent(ix,iy); 
	  (numIC[k]).at(ring) += 1;
	  //	  if (ring==1) std::cout<<ring<<" "<<(numIC[k]).at(ring)<<" "<<(sumIC[k]).at(ring)<<" "<<map[k]->GetBinContent(ix,iy)<<std::endl;
	    //	    std::cout<<"boh"<<std::endl;
	}
    }
  }
  //  getchar();

  int n=0;
  for(int k = 0; k < 2; ++k) {    
    for(int ix = 1; ix <= map[k] -> GetNbinsX(); ++ix) {
      for(int iy = 1; iy <= map[k] -> GetNbinsY(); ++iy)
	{
	  int zside;
	  if (k==0)  { zside=-1; }

	  else if (k==1) zside=1;

	  int ring = eRings->GetEndcapRing(ix,iy,k);
	  //	  if (k==1) std::cout<<"ring: "<<ring<<std::endl;
	  if( ring == -1 ) continue;

	  if( ring > 33)
	    {
	      mapNorm[k] -> Fill(ix,iy,0.);
	      continue;
	    }
	  else
	    {
	      //	      if (k==0) std::cout<<(numIC[k]).at(ring)<<" "<<(sumIC[k]).at(ring)<<std::endl;

	      if( (numIC[k]).at(ring) != 0 && (sumIC[k]).at(ring) != 0 ) {
		mapNorm[k] -> Fill(ix,iy,map[k]->GetBinContent(ix,iy)/((sumIC[k]).at(ring)/(numIC[k]).at(ring)));
		//		if (ring==1) std::cout<<ix<<" "<<iy<<" "<<map[k]->GetBinContent(ix,iy)<<" "<<map[k]->GetBinContent(ix,iy)/((sumIC[k]).at(ring)/(numIC[k]).at(ring))<<std::endl;
		//if (ring==1) std::cout<<ring<<" "<<(numIC[k]).at(ring)<<" "<<(sumIC[k]).at(ring)<<std::endl;

		//		if (k==0) { std::cout<<"riempio"<<std::endl; n++; }
	      }
	    }
	}
    }
  }

  //  std::cout<<n<<std::endl;

  o -> cd();

  map[0]->Draw("COLZ");
  map[1]->Draw("COLZ");

  map[0]->Write();
  map[1]->Write();
  mapNorm[0]->Write();
  mapNorm[1]->Write();

  o -> Close();


  return(0);

}
