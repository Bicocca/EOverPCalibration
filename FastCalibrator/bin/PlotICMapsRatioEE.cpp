/// Standalone program to mormalize IC EB by the mean on a eta ring + skipping xtal near dead channels and TT
/// in the normalization procedure
/// Folded Plots for Spread IC, Statistical Precision and spread
/// Correct IC near cracks and for momentum scale and produce txt IC values

#include <vector>
#include <utility>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include "TLatex.h"
#include "TFile.h"
#include "TGraph.h"

#include "TEndcapRings.h"

int main(int argc, char **argv)
{

  std::ifstream io1, io2, rms1, rms2;

  //  std::ifstream io3, io4, io6;
  io1.open ("output_EE_runD_SISCALIB_GAUSS_NOETA_STRAWEAK/IC_Run2012ABC_22JanuaryRereco_WZ_R9_EE_relative.txt");
  io2.open ("output_EE_runD_NOSCALIB/IC_Run2012ABC_22JanuaryRereco_WZ_R9_EE_relative.txt");

  TEndcapRings *eRings = new TEndcapRings();

  //  io3.open ("output_runD_10ITER_SISCALIB_ETABIN3/IC_Run2012ABC_22JanuaryRereco_WZ_Fbrem_EB_SISCALIB_relative.txt");
  //  io4.open ("output_runD_10ITER_SISCALIB_ETABIN4/IC_Run2012ABC_22JanuaryRereco_WZ_Fbrem_EB_SISCALIB_relative.txt");
  //  io6.open ("output_runD_10ITER_SISCALIB_ETABIN6/IC_Run2012ABC_22JanuaryRereco_WZ_Fbrem_EB_SISCALIB_relative.txt");
  //  rms1.open ("output_runD_10ITER_SISCALIB_ETALINEAR/RMSFile.txt");
  //  rms2.open ("output_runD_10ITER_NOSCALIB/RMSFile.txt");
  
  float status, IC, err;
  float status2, IC2, err2;

  //  float IC3,IC4,IC6;

  int x, y, x2, y2;

  TH2F *mapRatioEEp = new TH2F ("mapEEp", "mapEEp", 100, 0.5, 100.5, 100, 0.5, 100.5);
  TH2F *mapRatioEEm = new TH2F ("mapEEm", "mapEEm", 100, 0.5, 100.5, 100, 0.5, 100.5);

  TH1F *histoEEp = new TH1F ("histoEEp", "histoEEp", 100, 0.97, 1.03);
  TH1F *histoEEm = new TH1F ("histoEEm", "histoEEm", 100, 0.97, 1.03);

    TH2F *ICmap1EEp = new TH2F ("map1EEp", "map1EEp", 100, 0.5, 100.5, 100, 0.5, 100.5);
    TH2F *ICmap2EEp = new TH2F ("map2EEp", "map2EEp", 100, 0.5, 100.5, 100, 0.5, 100.5);
    TH2F *ICmap1EEm = new TH2F ("map1EEm", "map1EEm", 100, 0.5, 100.5, 100, 0.5, 100.5);
    TH2F *ICmap2EEm = new TH2F ("map2EEm", "map2EEm", 100, 0.5, 100.5, 100, 0.5, 100.5);

  TH1F *hIC1EEp = new TH1F ("hIC1EEp", "hIC1EEp", 100, 0.7, 1.3);
  TH1F *hIC2EEp = new TH1F ("hIC2EEp", "hIC2EEp", 100, 0.7, 1.3);
  TH1F *hIC1EEm = new TH1F ("hIC1EEm", "hIC1EEm", 100, 0.7, 1.3);
  TH1F *hIC2EEm = new TH1F ("hIC2EEm", "hIC2EEm", 100, 0.7, 1.3);

  TGraph *g_RMSEEp = new TGraph();
  TGraph *g_RMSEEm = new TGraph();

  //  TGraph *g_RMS3 = new TGraph();
  //  TGraph *g_RMS4 = new TGraph();
  //  TGraph *g_RMS6 = new TGraph();

  TH1F *histoEtaRingEEp[40];
  char histoNameEEp[100];
  char funcNameEEp[100];
  TH1F *histoEtaRingEEm[40];
  char histoNameEEm[100];
  char funcNameEEm[100];

  /*  TH1F *histoEtaRing4[86];
  char histoName4[100];
  TH1F *histoEtaRing3[86];
  char histoName3[100];
  TH1F *histoEtaRing6[86];
  char histoName6[100];
  */
  float mapEEp[100][100];
  float map2EEp[100][100];
  float mapEEm[100][100];
  float map2EEm[100][100];

  for (int e=0; e<40; e++) {
    sprintf(histoNameEEp,"h_ratio_EEp_%d",e);
    histoEtaRingEEp[e] = new TH1F(histoNameEEp,"",150,0.0,2.0);
    sprintf(histoNameEEm,"h_ratio_EEm_%d",e);
    histoEtaRingEEm[e] = new TH1F(histoNameEEm,"",150,0.0,2.0);

    /*    sprintf(histoName4,"h_ratio4_%d",e);
    histoEtaRing4[e] = new TH1F(histoName4,"",50,0.99,1.01);
    sprintf(histoName3,"h_ratio3_%d",e);
    histoEtaRing3[e] = new TH1F(histoName3,"",50,0.99,1.01);
    sprintf(histoName6,"h_ratio6_%d",e);
    histoEtaRing6[e] = new TH1F(histoName6,"",50,0.99,1.01);
    */  }

 for (int x=0; x<100; x++) {
  for (int y=0; y<100; y++) {
      mapEEp[x][y]=0.;
      map2EEp[x][y]=0.;
      mapEEm[x][y]=0.;
      map2EEm[x][y]=0.;
    }
  }


  //  mapRatio->SetDrawOption ("colz");
  
  x=0;
  float r;

  while (!io1.eof())
    {
      io1>>x>>y>>status>>IC>>err;
      io2>>x2>>y2>>status2>>IC2>>err2;
      //      std::cout<<x<<" "<<y<<" "<<IC<<" "<<status<<std::endl;
      //   std::cout<<x2<<" "<<y2<<" "<<IC2<<" "<<status2<<std::endl;


      //      io3>>x>>phi>>status>>IC3>>err;
      //      io4>>x>>phi>>status>>IC4>>err;
      //      io6>>x>>phi>>status>>IC6>>err;
	    //            r = sqrt ((x-50)*(x-50) + (y-50)*(y-50));

      if (status==-1) {
	  mapEEm[x-1][y-1]=IC;
	  map2EEm[x2-1][y2-1]=IC2;
      }
      else if (status==1) {
	  mapEEp[x-1][y-1]=IC;
	  map2EEp[x2-1][y2-1]=IC2;
      }

      //                  if (eRings->GetEndcapRing(x,y,status)==32)
      //	    std::cout<<eRings->GetEndcapRing(x2,y2,status2)<<" "<<IC<<" "<<IC2<<" "<<status2<<std::endl;

      if ( (status==-1) && (IC!=-1) && (status2==-1) && (IC2!=-1)) {
	if ((x==x2) && (y==y2)) {

	    histoEtaRingEEm[int(eRings->GetEndcapRing(x,y,0))]->Fill(IC/IC2);
	    //	    if (x==34 && y==47)  std::cout<<"anello: "<<eRings->GetEndcapRing(x,y,status)<<std::endl;
	      if (eRings->GetEndcapRing(x,y,0)==16)
	      std::cout<<x<<" "<<y<<" "<<IC<<" "<<IC2<<std::endl;
	  //	  std::cout<<x<<" "<<y<<" "<<etaRing<<std::endl;
	  //	  histoEtaRing3[int(fabs(eta))]->Fill(IC3/IC2);
	  //	  histoEtaRing4[int(fabs(eta))]->Fill(IC4/IC2);
	  //	  histoEtaRing6[int(fabs(eta))]->Fill(IC6/IC2);

	  mapRatioEEm->SetBinContent (x-1,y-1,IC/IC2);
	  histoEEm->Fill(IC/IC2);
	  //	  std::cout<<IC/IC2<<std::endl;

	  ICmap1EEm->SetBinContent (x-1,y-1,IC);
	  ICmap2EEm->SetBinContent (x2-1,y2-1,IC2);
	  //	  std::cout<<IC-IC2<<std::endl;
	  hIC1EEm->Fill(IC);
	  hIC2EEm->Fill(IC2);
	
	}
	else
	  std::cout<<"Problem: incoherent x or y "<<std::endl;
	
       }

      else if ((status==1) && (IC!=-1) && (status2==1) && (IC2!=-1)) {
	if ((x==x2) && (y==y2)) {

	  histoEtaRingEEp[int(eRings->GetEndcapRing(x,y,status))]->Fill(IC/IC2);

	  //	  histoEtaRing3[int(fabs(eta))]->Fill(IC3/IC2);
	  //	  histoEtaRing4[int(fabs(eta))]->Fill(IC4/IC2);
	  //	  histoEtaRing6[int(fabs(eta))]->Fill(IC6/IC2);

	  mapRatioEEp->SetBinContent (x-1,y-1,IC/IC2);
	  histoEEp->Fill(IC/IC2);
	  //	  std::cout<<IC/IC2<<std::endl;

	  ICmap1EEp->SetBinContent (x-1,y-1,IC);
	  ICmap2EEp->SetBinContent (x2-1,y2-1,IC2);
	  //	  std::cout<<IC-IC2<<std::endl;
	  hIC1EEp->Fill(IC);
	  hIC2EEp->Fill(IC2);
	
	}
	else
	  std::cout<<"Problem: incoherent x or y "<<std::endl;
	
       }

    }



  TFile f1 ("confronti.root", "RECREATE");
  f1.cd();

    for (int e=0; e<40; e++) {

      sprintf(funcNameEEp,"f_EEp_%d",e);
      TF1* fgausEEp = new TF1(funcNameEEp,"gaus",0.1,1.9);
      fgausEEp -> SetParameter(1,histoEtaRingEEp[e]->GetMean());
      fgausEEp -> SetParameter(2,histoEtaRingEEp[e]->GetRMS());
      //      histoEtaRingEEp[e] -> Fit(funcNameEEp,"QS+","",1-histoEtaRingEEp[e]->GetRMS(),1+histoEtaRingEEp[e]->GetRMS());
      //      g_RMSEEp->SetPoint (e, float(e), fgausEEp->GetParameter(2));

      sprintf(funcNameEEm,"f_EEm_%d",e);
      TF1* fgausEEm = new TF1(funcNameEEm,"gaus",0.1,1.9);
      fgausEEm -> SetParameter(1,histoEtaRingEEm[e]->GetMean());
      fgausEEm -> SetParameter(2,histoEtaRingEEm[e]->GetRMS());
      //      histoEtaRingEEm[e] -> Fit(funcNameEEm,"QS+","",1-histoEtaRingEEm[e]->GetRMS(),1+histoEtaRingEEm[e]->GetRMS());
      //      g_RMSEEm->SetPoint (e, float(e), fgausEEm->GetParameter(2));

          g_RMSEEp->SetPoint (e, float(e), histoEtaRingEEp[e]->GetRMS());
          g_RMSEEm->SetPoint (e, float(e), histoEtaRingEEm[e]->GetRMS());    

    //    g_RMS3->SetPoint (e-1, float(e), histoEtaRing3[e]->GetRMS());
    //    g_RMS4->SetPoint (e-1, float(e), histoEtaRing4[e]->GetRMS());
    //    g_RMS6->SetPoint (e-1, float(e), histoEtaRing6[e]->GetRMS());
    //    g_RMS->SetPointError (e-1, 0, 0);
    if (e==1 || e==2 || e==4 || e==19 || e==26 || e==30 || e==31 || e==32 || e==33 ) {
      histoEtaRingEEp[e]->Draw();
      //      fgausEEp->Draw("same");

      histoEtaRingEEm[e]->Draw();
      //      fgausEEm->Draw("same");

      histoEtaRingEEp[e]->Write();
      histoEtaRingEEm[e]->Write();
    }
  }
    
  TCanvas *c1 = new TCanvas("c1");
  c1->cd();
  g_RMSEEp -> GetXaxis() -> SetTitle("i|#eta|");
  g_RMSEEp -> GetYaxis() -> SetTitle("RMS");
  g_RMSEEp -> SetMinimum(0.00000);
  //      g_RMSEEp -> SetMaximum(0.015);
  //    g_RMSEEm -> SetMaximum(0.015);
  g_RMSEEp -> SetMarkerStyle(20);
  g_RMSEEp -> SetMarkerSize(1.0);
  g_RMSEEp -> SetMarkerColor(kBlue+1);
  c1->SetGrid();
  g_RMSEEp -> Draw("AP");

    g_RMSEEp -> GetXaxis() -> SetRangeUser (0,32.5);
  
  /*  g_RMS3 -> SetMarkerStyle(20);
  g_RMS3 -> SetMarkerSize(1.0);
  g_RMS3 -> SetMarkerColor(51+1);
  g_RMS4 -> SetMarkerStyle(20);
  g_RMS4 -> SetMarkerSize(1.0);
  g_RMS4 -> SetMarkerColor(kRed+1);
  g_RMS6 -> SetMarkerStyle(20);
  g_RMS6 -> SetMarkerSize(1.0);
  g_RMS6 -> SetMarkerColor(kOrange+1);

      g_RMS3 -> Draw("PLsame");
     g_RMS4 -> Draw("PLsame");
    g_RMS6 -> Draw("PLsame");
  
    TLegend* leg = new TLegend(0.15,0.72,0.43,0.89);
    leg -> SetFillColor(0);
    leg -> SetTextFont(42);
    leg -> SetTextSize(0.05);
    leg -> AddEntry(g_RMS,"Miscalib. 2%","P");
    leg -> AddEntry(g_RMS3,"Miscalib. 3%","P");
    leg -> AddEntry(g_RMS4,"Miscalib. 4%","P");
    leg -> AddEntry(g_RMS6,"Miscalib. 6%","P");
    leg -> Draw("same");
  */
    c1->Print("g_RMS_EE+.png","png");


  TCanvas *c2 = new TCanvas("c2");
  c2->cd();
  g_RMSEEm -> GetXaxis() -> SetTitle("i|#eta|");
  g_RMSEEm -> GetYaxis() -> SetTitle("RMS");
  g_RMSEEm -> SetMinimum(0.00000);
  //    g_RMS -> SetMaximum(0.0025);
  g_RMSEEm -> SetMarkerStyle(20);
  g_RMSEEm -> SetMarkerSize(1.0);
  g_RMSEEm -> SetMarkerColor(kBlue+1);
  c2->SetGrid();
  g_RMSEEm -> Draw("AP");
  g_RMSEEm -> GetXaxis() -> SetRangeUser (0,33.5);
  c2->Print("g_RMS_EE-.png","png");

  mapRatioEEp -> GetXaxis() -> SetTitle("ix");
  mapRatioEEp -> GetYaxis() -> SetRangeUser(-100.5,100.5);
  mapRatioEEp -> GetYaxis() -> SetTitle("iy");
  histoEEp -> GetXaxis() -> SetTitle("IC1/IC2");
  histoEEp -> GetYaxis() -> SetTitle("N");
  histoEEp->SetStats(1);
  mapRatioEEp->GetZaxis()->SetRangeUser(0.98, 1.02);

  mapRatioEEm -> GetXaxis() -> SetTitle("ix");
  mapRatioEEm -> GetYaxis() -> SetRangeUser(-100.5,100.5);
  mapRatioEEm -> GetYaxis() -> SetTitle("iy");
  histoEEm -> GetXaxis() -> SetTitle("IC1/IC2");
  histoEEm -> GetYaxis() -> SetTitle("N");
  histoEEm->SetStats(1);
  mapRatioEEm->GetZaxis()->SetRangeUser(0.98, 1.02);

  mapRatioEEp->Write();
  histoEEp->Write();
  hIC1EEp->Write();
  hIC2EEp->Write();
  ICmap1EEp->Write();
  ICmap2EEp->Write();
  g_RMSEEp->Write("g_RMSEEp");
  g_RMSEEm->Write("g_RMSEEm");
  mapRatioEEm->Write();
  histoEEm->Write();
  hIC1EEm->Write();
  hIC2EEm->Write();
  ICmap1EEm->Write();
  ICmap2EEm->Write();

  f1.Close();

      return 0;
}
