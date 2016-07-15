/*
A script to the scaling function from several nuclei using model params provided in Nadia's thesis page 143
From the model:
f(y) = (f0-B)*alpha**2*exp(-(a*y)**2)/(alpha**2+y**2)+bigB*exp((-b*abs(y))^2) 


Buddhini - 06-03-2016 
To run:
in root do-
.L fy_vs_y_for_several_nuclei.C+g
draw_fy_vs_y()

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
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TString.h>
#include <TMath.h>

using namespace std;

Double_t H2_fy_from_model(Double_t x);
Double_t He3_fy_from_model(Double_t x);
Double_t He4_fy_from_model(Double_t x);
Double_t Be9_fy_from_model(Double_t x);
Double_t C12_fy_from_model(Double_t x);
Double_t Cu63_fy_from_model(Double_t x);
Double_t Au197_fy_from_model(Double_t x);

Double_t al_fy_from_model(Double_t x);

void draw_fy_vs_y(){

  // histo parameters
  gStyle->SetTitleYOffset(1.4);
  gStyle->SetTitleXOffset(1.0);
  gStyle->SetLabelSize(0.06,"x");
  gStyle->SetLabelSize(0.06,"y");
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetTitleX(0.4);
  gStyle->SetTitleW(0.6);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleFontSize(0.09);
  gStyle->SetStripDecimals(0);

 
  // Fit and stat parameters
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  gStyle->SetStatFont(132);
  gStyle->SetStatY(0.99);
  gStyle->SetStatX(0.99);
  gStyle->SetStatW(0.20);
  gStyle->SetStatH(0.20);

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
  

  // Open canvas to draw the comparision between Fy vs Y using sig_exp and sig_exp corrected for sig_dis
  TCanvas * c = new TCanvas("c", "",0,0,1000,600);  
  c->Draw();

  // 2H:
  TF1 *fh2 = new TF1("fh2"," H2_fy_from_model(x)",-0.6,0.2); 
  fh2->SetLineColor(kBlack);
  fh2->SetLineStyle(1);
  fh2->Draw();

  //3He
  TF1 *fhe3 = new TF1("fhe3"," He3_fy_from_model(x)",-0.6,0.2); 
  fhe3->SetLineColor(kBlue);
  fhe3->SetLineStyle(1);
  fhe3->Draw("same");
  fhe3->GetXaxis()->SetTitle("Y(GeV^-1)");
  //4He
  TF1 *fhe4 = new TF1("fhe4"," He4_fy_from_model(x)",-0.6,0.2); 
  fhe4->SetLineColor(kRed);
  fhe4->SetLineStyle(1);
  fhe4->Draw("same");

  //9Be
  TF1 *fbe9 = new TF1("fbe9"," Be9_fy_from_model(x)",-0.6,0.2); 
  fbe9->SetLineColor(kGreen);
  fbe9->SetLineStyle(1);
  fbe9->Draw("same");

  //12C
  TF1 *fc12 = new TF1("fc12"," C12_fy_from_model(x)",-0.6,0.2); 
  fc12->SetLineColor(kOrange);
  fc12->SetLineStyle(1);
  fc12->Draw("same");

  //63Cu
  TF1 *fcu63 = new TF1("fcu63"," Cu63_fy_from_model(x)",-0.6,0.2); 
  fcu63->SetLineColor(kRed-8);
  fcu63->SetLineStyle(1);
  fcu63->Draw("same");

  //197Au
  TF1 *fau197 = new TF1("fau197"," Au197_fy_from_model(x)",-0.6,0.2); 
  fau197->SetLineColor(kBlue+9);
  fau197->SetLineStyle(1);
  fau197->Draw("same");

  TLegend * leg= new TLegend(0.2,0.4,0.35,0.8);
  leg->AddEntry(fh2,"2H","l");
  leg->AddEntry(fhe3,"3He","l");
  leg->AddEntry(fhe4,"4He","l");
  leg->AddEntry(fbe9,"9Be","l");
  leg->AddEntry(fc12,"12C","l");
  leg->AddEntry(fcu63,"63Cu","l");
  leg->AddEntry(fau197,"197Au","l");

  leg->Draw();

  // c->Print(Form("xsec_rad_exp_mod_ratio_x_geq_%2.1f_and_x_leq_%2.1f_%s_dummy_%s_%s.png",xlow,xup,tgt_pos.Data(),angle.Data(),fit.Data()));

  
}


Double_t H2_fy_from_model(Double_t x){

  Double_t f_0= 8.742;// GeV^-1
  Double_t B = 0.8239;// GeV^-1
  Double_t a = 7.727;// GeV
  Double_t b = 9.394;// GeV^-1
  Double_t alpha =  45.3E-03; // GeV

  Double_t val = (f_0-B)*(TMath::Power(alpha,2)*TMath::Exp(-1*TMath::Power((a*x),2)))/(TMath::Power(alpha,2)+TMath::Power(x,2)) + B*TMath::Exp(-1*TMath::Power((b*TMath::Abs(x)),2));

  return val;
}


Double_t He3_fy_from_model(Double_t x){

  // From Nadia's thesis page 143. Changed units from MeV^-1 in Nadia's table to GeV^-1
  Double_t f_0 = 5.309; //GeV^-1
  Double_t B = 2.184;// GeV^-1
  Double_t a = 2.886;// GeV^-1
  Double_t b =1.035;// GeV^-1
  Double_t alpha = 64.2e-3;// GeV

  Double_t val = (f_0-B)*(TMath::Power(alpha,2)*TMath::Exp(-1*TMath::Power((a*x),2)))/(TMath::Power(alpha,2)+TMath::Power(x,2)) + B*TMath::Exp(-1*TMath::Power((b*TMath::Abs(x)),2));
 return val;
} 


Double_t He4_fy_from_model(Double_t x){

  // From Nadia's thesis page 143. Changed units from MeV^-1 in Nadia's table to GeV^-1
  Double_t f_0 = 4.020; //GeV^-1
  Double_t B = 1.345;// GeV^-1
  Double_t a = 2.699;// GeV^-1
  Double_t b =7.494;// GeV^-1
  Double_t alpha = 100.2e-3;// GeV

  Double_t val = (f_0-B)*(TMath::Power(alpha,2)*TMath::Exp(-1*TMath::Power((a*x),2)))/(TMath::Power(alpha,2)+TMath::Power(x,2)) + B*TMath::Exp(-1*TMath::Power((b*TMath::Abs(x)),2));
 return val;
} 


Double_t Be9_fy_from_model(Double_t x){

  // From Nadia's thesis page 143. Changed units from MeV^-1 in Nadia's table to GeV^-1
  Double_t f_0 = 3.481; //GeV^-1
  Double_t B = 1.161;// GeV^-1
  Double_t a = 3.120;// GeV^-1
  Double_t b =7.840;// GeV^-1
  Double_t alpha = 110.0e-3;// GeV

  Double_t val = (f_0-B)*(TMath::Power(alpha,2)*TMath::Exp(-1*TMath::Power((a*x),2)))/(TMath::Power(alpha,2)+TMath::Power(x,2)) + B*TMath::Exp(-1*TMath::Power((b*TMath::Abs(x)),2));
 return val;
} 

Double_t C12_fy_from_model(Double_t x){
  
  Double_t f_0 = 3.182;
  Double_t B = 1.359;
  Double_t a = 3.027;
  Double_t b =7.050;
  Double_t alpha = 137.2e-3;

  Double_t val = (f_0-B)*(TMath::Power(alpha,2)*TMath::Exp(-1*TMath::Power((a*x),2)))/(TMath::Power(alpha,2)+TMath::Power(x,2)) + B*TMath::Exp(-1*TMath::Power((b*TMath::Abs(x)),2));
 return val;
} 


Double_t Cu63_fy_from_model(Double_t x){
  
  Double_t f_0 = 2.874;
  Double_t B = 0.8866;
  Double_t a = 3.096;
  Double_t b =7.094;
  Double_t alpha = 132.4e-3;

  Double_t val = (f_0-B)*(TMath::Power(alpha,2)*TMath::Exp(-1*TMath::Power((a*x),2)))/(TMath::Power(alpha,2)+TMath::Power(x,2)) + B*TMath::Exp(-1*TMath::Power((b*TMath::Abs(x)),2));
 return val;
} 


Double_t Au197_fy_from_model(Double_t x){
  
  Double_t f_0 = 2.642;
  Double_t B = 0.7632;
  Double_t a = 3.065;
  Double_t b =6.768;
  Double_t alpha = 132.4e-3;

  Double_t val = (f_0-B)*(TMath::Power(alpha,2)*TMath::Exp(-1*TMath::Power((a*x),2)))/(TMath::Power(alpha,2)+TMath::Power(x,2)) + B*TMath::Exp(-1*TMath::Power((b*TMath::Abs(x)),2));
 return val;
} 



Double_t al_fy_from_model(Double_t x){

  // From Aji's radiative correction files
  Double_t f_0= 3.610000000000000;// GeV^-1
  Double_t B = 5.999999999999999;// GeV^-1
  Double_t a = 4.930000000000001;// GeV
  Double_t b = 6.000000000000000;// GeV^-1
  Double_t alpha =  156.000000000000E-003; // GeV


  Double_t val = (f_0-B)*(TMath::Power(alpha,2)*TMath::Exp(-1*TMath::Power((a*x),2)))/(TMath::Power(alpha,2)+TMath::Power(x,2)) + B*TMath::Exp(-1*TMath::Power((b*TMath::Abs(x)),2));

  return val;
}




