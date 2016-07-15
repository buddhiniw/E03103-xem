
/*
A script to compare the data and simulation histograms from dummy_ana.f
This compares the simulated and calculated acceptance parameters
Buddhini - 03-19-2016
Buddhini - 07-14-2016 Add 2-D plots for delta vs ytar and xptar vs yptar

To run:
in root do-
.L compare_dat_sim.C
compare_dat_sim("18deg_3740_dummy_us")
 */


#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TVectorD.h>
#include <dirent.h> 
#include <stdio.h> 
#include <TStyle.h> 
#include <TDirectory.h> 
#include <TCanvas.h> 
#include <TGraph.h> 
#include <TAxis.h> 
#include <TF1.h> 
#include <TH1.h> 
#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TString.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TFile.h>
#include <THStack.h>
#include <TCanvas.h>

void compare_dat_sim(TString rfile){

  // histo parameters
  gStyle->SetTitleYOffset(1.6);
  gStyle->SetTitleXOffset(1.6);
  gStyle->SetLabelSize(0.07,"x");
  gStyle->SetLabelSize(0.07,"y");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetTitleX(0.3);
  gStyle->SetTitleW(0.6);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleFontSize(0.09);

 
  // Fit and stat parameters
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetStatY(0.99);
  gStyle->SetStatX(0.99);
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.15);
  
  Char_t filename[200]; 

  //Get the root file
  //sprintf(filename,"../hbook/replayed/%s.root",rfile.Data());
  //sprintf(filename,"../hbook/%s.root",rfile.Data());
  //sprintf(filename,"../hbook/carbon/%s.root",rfile.Data());
  sprintf(filename,"../hbook/aluminum/pass0/%s.root",rfile.Data());

  TFile * f = new TFile(filename);
  if(!f->IsOpen())
    exit(1);
  if (f->IsZombie())
    exit(1);

  std::cout<<"Opened "<<filename<<"\n";

  /******************
   Draw all the 1D histograms
  *******************/
 TCanvas * c = new TCanvas("c", "",0,0,1800,1000);  
  c->Divide(4,4);

  
  // e data histos idbase = 1
  // e sim histos idbase = 12
  /*          call hfill(3000+idbase,sngl(hsdelta),zero,sngl(weight))
              call hfill(3100+idbase,sngl(hsxptar),zero,sngl(weight))
              call hfill(3200+idbase,sngl(hsyptar),zero,sngl(weight))
              call hfill(3300+idbase,sngl(hsytar),zero,sngl(weight))
              call hfill(3400+idbase,sngl(hsxfp),zero,sngl(weight))
              call hfill(3500+idbase,sngl(hsxpfp),zero,sngl(weight))
              call hfill(3600+idbase,sngl(hsyfp),zero,sngl(weight))
              call hfill(3700+idbase,sngl(hsypfp),zero,sngl(weight))
              call hfill(3800+idbase,sngl(W),zero,sngl(weight))
              call hfill(3900+idbase,sngl(Q2),zero,sngl(weight))
              call hfill(4200+idbase,sngl(hstheta-th_rad),zero,sngl(weight))
              call hfill(4300+idbase,sngl(eprime),zero,sngl(weight))
              call hfill(4400+idbase,sngl(hsdelta),sngl(hstheta-th_rad),sngl(weight))         
              call hfill(4100+idbase,sngl(xi),zero,sngl(weight))
              call hfill(4000+idbase,sngl(x_bj),zero,sngl(weight))
	      call hfill(4500+idbase,sngl(hsdelta),sngl(hsytar),sngl(weight))	      
	      call hfill(4600+idbase,sngl(hsxptar),sngl(hsyptar),sngl(weight))	      


C     HISTOGRAM KEY
C     1X01     DATA 
C     1X02     DATA (ELCLEAN)
C     1X03     DATA(ELREAL)/DATA(ELCLEAN)
C     1X04     POSITRONS
C     1X05     DATA*(1-POSITRONS/DATA(ELCLEAN)) = DATA_POSCOR
C     1X06     DUMMY
C     1X07     DUMMY (ELCLEAN)
C     1X08     DUMMY(ELREAL)/DUMMY(ELCLEAN)
C     1X09     DUMMY POSITRONS
C     1X10     DUMMY*( 1-(DUMMY POSITRONS)/DUMMY(ELCLEAN) ) = DUMMY_POSCOR
C     1X11     DATA_POSCOR - DUMMY_POSCOR/7.2 * (DATA_CHARGE/DUMMY_CHARGE)
C     1X12     SIMC
C     1X13     DATA/SIMC
C     1X14     (DATA/SIMC)*MODEL
C     1X15     POSITRONS/DATA(ELCLEAN)


What we should compare to check if the data matches simulation is SIMC (12) to DATA*(1-POSITRONS/DATA(ELCLEAN)) = DATA_POSCOR (5), the positron corrected dat or DATA_POSCOR - DUMMY_POSCOR/7.2 * (DATA_CHARGE/DUMMY_CHARGE) (11) same as positron corrected data since for the dummy target theres is DUMMY_POSCOR = 0.

So compare idbase = 5/11 to sim idbase = 12
  */

  THStack *hs[15];
  Double_t didbase = 0;
  Double_t sidbase = 0;

  Double_t histid[15]={3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4200,4300,4400,4100,4000};
  TString histnam[15]={"hsdelta","hsxptar","hsyptar","hsytar","hsxfp","hsxpfp","hsyfp","hsypfp","W","Q2","hstheta","eprime","th_vs_delta","xi","x_bj"};

  // loop over the 15 types of histos from data and sim and stack them
  for (Int_t i=0;i<15; i++){
    // std::cout<<i<<std::endl;
    c->cd(i+1);
    // gPad->SetLogy();
    hs[i] = new THStack("fs",Form("%s",histnam[i].Data()));
    // first get positron corrected electron data. This could be idbase = 5 or idbase = 11 (for dummy targets only)
    didbase =5;
    TH1F *hdat = (TH1F*)f->Get(Form("h%2.0f",histid[i]+didbase));
    //std::cout<<Form("h%2.0f",histid[i]+didbase)<<std::endl;
    // next get sim histos
    sidbase = 12;
    TH1F *hsim = (TH1F*)f->Get(Form("h%2.0f",histid[i]+sidbase));
    hsim->SetLineColor(kRed);
    // std::cout<<Form("h%2.0f",histid[i]+sidbase)<<std::endl;
  
    hs[i]->Add(hsim);
    hs[i]->Add(hdat);
    hs[i]->Draw("nostack,hist");
    gPad->Update();   

    TLegend * leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(hdat,Form("data(1x%2.0f)",didbase),"l");
    leg->AddEntry(hsim,"simc(1x12)","l");
    leg->Draw();

    // c1->cd(i+1);
    // idbase = 1;
    // TH1F *hdat1 = (TH1F*)f->Get(Form("h%2.0f",histid[i]+idbase));
    // //hdat1->SetLineColor(kGreen+2);
    // //hdat1->Add(hsim,0.0);
    // hdat1->Draw();
    // hdat1->GetXaxis()->SetLabelSize(0.07);
    // hdat1->GetYaxis()->SetLabelSize(0.07);
    // gPad->Update();
    
  }

  c->cd(16);
  TPaveText *pt = new TPaveText(.05,.1,.95,.8);
  
  pt->AddText("ROOT file used:");
  pt->AddText(rfile);
  pt->AddLine(.0,.5,1.,.5);
  pt->Draw();


  /******************
   Draw the 2D histograms
  *******************/



  //c->Print(Form("~/MyWork/HallC/SRC/xem/scripts/compare_dat_sim/stack_dat_%1.0f_sim_12_%s.png",didbase,rfile.Data()));
  //c->Print(Form("~/MyWork/HallC/SRC/xem/scripts/compare_dat_sim/stack_dat_%1.0f_sim_12_%s.C",didbase,rfile.Data()));

  //c1->Print(Form("diff_dat_sim_%s.png",rfile.Data()));
  //c1->Print(Form("diff_dat_sim_%s.C",rfile.Data()));

}
