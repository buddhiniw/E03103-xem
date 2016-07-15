
/*
A script to extract the sig_raw, sig_qe, sig_dis from the xsec model from the radiative corrections 
files and compare them to data for all the dummy angles


Buddhini - 05-14-2015 

To run:
in root do-
.L compare_mod_data_xsec_for_all_angles.c
compare_mod_data_xsec_for_all_angles("us")
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
#include <TROOT.h>
#include "compare_mod_data_xsec.C"

using namespace std;

void compare_mod_data_xsec_for_all_angles(TString tgt_pos)
{
  // histo parameters
  gStyle->SetTitleYOffset(1.4);
  gStyle->SetTitleXOffset(1.6);
  gStyle->SetLabelSize(0.07,"x");
  gStyle->SetLabelSize(0.07,"y");
  gStyle->SetTitleSize(0.07,"x");
  gStyle->SetTitleSize(0.07,"y");
  gStyle->SetTitleX(0.3);
  gStyle->SetTitleW(0.6);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleFontSize(0.09);
  gStyle->SetStripDecimals(0);

 
  // Fit and stat parameters
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  gStyle->SetStatY(0.99);
  gStyle->SetStatX(0.99);
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.15);

  //Pad parameters
  gStyle->SetPadColor(0); 
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetPadTopMargin(0.18);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  
 //Delete all the objects stored in the current directory memmory
  gDirectory->Delete("*");
  
  // Open a canvas
  TCanvas * c = new TCanvas("c", "",0,0,1600,1000);  
  c->Draw();
  c->Divide(2,3);

  Double_t angle[6] = {18,22,26,32,40,50};

  for(int i=0;i<6;i++)
    { 
      c->cd(i+1);
      compare_mod_data_xsec(angle[i],tgt_pos);
      gPad->Update();
    }

  //  c->Print(Form("xsec_rad_exp_mod_ratio_%s_dummy_all_angles_x_geq_%2.1f_and_x_leq_%2.1f_fit_%s.png",tgt_pos.Data(),xlow,xup,fit.Data()));
}


