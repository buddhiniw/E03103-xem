/*
A script to draw get the scaling function parameters by fitting Fy vs y determined by dummy_ana.f for A>2 targets.

Buddhini - 05-21-2015 
Buddhini - 05-17-2015 Fixing the fit
Buddhini - 06-15-2016 Add 5% fixed errors to each Fy point
To run:
in root do-
.L fit_fy_vs_y.C+g
fit_fy_vs_y(18,"c12")

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
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TString.h>
#include <TMath.h>

using namespace std;

void myfunc();

Bool_t upstream = false;


Double_t myfunction(Double_t *y, Double_t *par);

Double_t calculate_fy_including_dis(Double_t sig, Double_t fact);
Double_t c12_fy_from_model(Double_t x);
Double_t al_fy_from_model(Double_t x);
Double_t al_fy_from_pass0(Double_t x);
void fit_fy_vs_y(Double_t angle, TString tgt_pos, TCanvas *c);


void fit_fy_vs_y_for_al(Double_t angle){
  TCanvas * c1 = new TCanvas("c1", "",0,0,800,1000);  
  TCanvas * c2 = new TCanvas("c2", "",0,0,800,1000);  
  std::cout<<"US target\n";
  upstream = true;
  fit_fy_vs_y(angle, "us",c1);

  std::cout<<"DS target\n";
  upstream = false;
  fit_fy_vs_y(angle, "ds",c2);

}


void fit_fy_vs_y(Double_t angle, TString tgt_pos,TCanvas *c){

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
  
  // Open a file to output the fy values
  FILE *fp;
  fp = fopen(Form("%2.0f_al_fy_%s.txt",angle,tgt_pos.Data()), "w+");
  // fprintf(fp, This is testing for fprintf...\n");
  // fputs("This is testing for fputs...\n", fp);
  
  //Delete all the objects stored in the current directory memmory
  gDirectory->Delete("*");
  
  // select files for a target with the given angle
  DIR *dir;
  struct dirent *ent;
  std::vector<string> filenames;

  // TString dirtry = "/home/buddhini/MyWork/HallC/SRC/xem/OUT/with_fy/";
  //TString dirtry = "/home/buddhini/MyWork/HallC/SRC/xem/OUT/with_dis_kine_fact/";
  //TString dirtry = "/home/buddhini/MyWork/HallC/SRC/xem/OUT/from_singles_ana/";
  //TString dirtry = "/home/buddhini/MyWork/HallC/SRC/xem/OUT/using_ajis_hist_files/";
  ///TString dirtry = "/home/buddhini/MyWork/HallC/SRC/xem/OUT/form_fac_fix/";
  //TString dirtry = "/home/buddhini/MyWork/HallC/SRC/xem/OUT/using_my_simhisto_for_c12/";
  // TString dirtry = "/home/buddhini/MyWork/HallC/SRC/xem/OUT/using_my_simhisto_for_al/pass0/";
  // TString dirtry = "/home/buddhini/MyWork/HallC/SRC/xem/OUT/using_my_simhisto_for_al/using_pass0_seperate_fits/";
  TString dirtry = "/home/buddhini/MyWork/HallC/SRC/xem/OUT/using_my_simhisto_for_al/pass0/"; //Aji's model
  //TString dirtry = "/home/buddhini/MyWork/HallC/SRC/xem/OUT/using_my_simhisto_for_al/pass1_using_full_dataset/"; //pass0 Fy fit. US and DS are treated seperately


  if ((dir = opendir (dirtry)) != NULL) 
    {
      std::cout<<"Opened dir = "<<dirtry<<std::endl;
      
      /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL) 
	{  
	  if((TString (ent->d_name)).Contains(Form("%2.0fdeg",angle)) && (TString (ent->d_name)).Contains(tgt_pos))
	    {
	      filenames.push_back(Form("%s/%s",dirtry.Data(),ent->d_name));
	      std::cout<<"Found : "<<ent->d_name<<std::endl;
	    }	
	}
      
      closedir (dir);
    } 
  else {
    
    /* could not open directory */
    std::cout<<"Unable to open directory!\n";
    perror ("");
  }
  


 
  // Get fy and y from *.deltabin.xsec files
  std::string line;
  std::string ys;
  std::string fys;
  std::string efys;
  std::string q2s;
  std::string sigrads;
  std::string sigdiss;
  std::string facts;


  TVectorD fyv(0),yv(0),q2v(0);
  TVectorD fyev(0), yev(0);

  TVectorD sigdisv(0),sigradv(0),factv(0);
  TVectorD fyuncorrected(0); //fy with dis


  //get the size of the file array
  int j = filenames.size();
  int k;
  int kv = 0;
  for(int i=0;i<j;i++)
    { 
      ifstream myfile ((char*)(filenames.at(i)).c_str());
      k=0;

      if (myfile.is_open())
	{

	  while ( getline (myfile,line) )
	    {
	      if(k<3) 
		{
		  //skip the first two lines
		  k++;
		}
	      else 
		{
		  k++;

		  // first read y (GeV/c). If its in the range -1<y<-0.1
		  yv.ResizeTo(kv+1);
		  yev.ResizeTo(kv+1);
		  fyv.ResizeTo(kv+1);
		  fyev.ResizeTo(kv+1);
		  q2v.ResizeTo(kv+1);
		  sigradv.ResizeTo(kv+1);
		  sigdisv.ResizeTo(kv+1);
		  factv.ResizeTo(kv+1);
		  fyuncorrected.ResizeTo(kv+1);

		  //std::cout << line << '\n';
		  line.erase (1,35);
		  q2s =line.substr (1,7);
		  q2v.operator()(kv) = atof((char*)(q2s.c_str()));
		  //std::cout << "Q2 = "<<q2s << '\n';

		  //std::cout << line << '\n';
		  line.erase (1,7);
		  ys =line.substr (1,10);
		  yv.operator()(kv) = atof((char*)(ys.c_str())); // y is in units of GeV/c.
		  yev.operator()(kv) = 0.0; //no errors on y for now
		  // std::cout << "y = "<<ys << '\n';

		  line.erase (1,10);
		  //std::cout << line << '\n';
		  fys = line.substr (1,11);
		  fyv.operator()(kv) = (atof((char*)(fys.c_str())))*1e-3; //fy is read in as TeV^-1. So convert to GeV^-1
		  //std::cout << "fy = "<<fyv(kv) << '\n';

		  line.erase (1,12);
		  // std::cout << line << '\n';
		  efys = line.substr (1,11);
		  // stat error
		  fyev.operator()(kv) = (atof((char*)(efys.c_str())))*1e-3 + fyv(kv)*0.05; 
		  // efy is read in as TeV^-1. So convert to GeV^-1
		  // added a fy*5% systematic to each point too.
		  

		  //std::cout << "(y,fys) = "<<yv(kv)<<","<<fyv(kv)<< '\n';
		  //fyev.operator()(kv) = fyv(kv)*0.05; //add 5% errors
		  std::cout << "(y,fy,fye) = "<<yv(kv)<<","<<fyv(kv)<<","<<fyev(kv)<< '\n';

		  // write to file
		  fputs(Form("%f, %f,  %f\n",yv(kv),fyv(kv),fyev(kv)), fp);

		  // factor including the kinematics and the n and p cross-sections:
		  /*
		    dwdy = q3(i)/sqrt(mp**2+q3(i)**2+y(i)**2+2.0*q3(i)*y(i))
		    
		    if(targ_a.gt.1.0) then
		    fact(i) = (targ_Z*sig_p(i)+(targ_a-targ_Z)*sig_n(i))/dwdy
		    else
		    fact(i) = targ_Z*sig_p(i)/dwdy

		   */
		  line.erase (1,12);
		  //std::cout << line << '\n';
		  facts = line.substr (1,12);
		  factv.operator()(kv) = atof((char*)(facts.c_str())); //units of ster/mB
		  //std::cout << " factor = "<<factv(kv)*0.0001<< '\n';

		  // Now to get the information needed to calculate fy without the 
		  // dis subtraction.
		  // radiated experimental cross section
		  line.erase (1,12);
		  //std::cout << line << '\n';
		  sigrads = line.substr (1,12);
		  sigradv.operator()(kv) = atof((char*)(sigrads.c_str())); //units of nB(ster MeV)
		  // std::cout << "sig rad = "<<sigradv(kv)<< '\n';

		  // dis cross section from the XEM model
		  //line.erase (1,12);
		  //std::cout << line << '\n';
		  sigdiss = line.substr (1,12);
		  sigdisv.operator()(kv) = atof((char*)(sigdiss.c_str()));//units of nB(ster MeV)
		  //std::cout << "sig dis = "<<sigdisv(kv)<< '\n';

		  //calculated fy with dis using the 
		  fyuncorrected.operator()(kv) = calculate_fy_including_dis(sigradv(kv), factv(kv));// get into units of MeV
		  //std::cout << "fy using dis uncorrected exp xsection = "<<fyuncorrected(kv)<< '\n';
		  kv++;

		}
	    }
	  myfile.close();
	}
      else std::cout << "Unable to open file :"<<filenames.at(i) <<std::endl; 


    }

  myfunc();

  // Open canvas to draw the comparision between Fy vs Y using sig_exp and sig_exp corrected for sig_dis
  //TCanvas * c = new TCanvas("c", "",0,0,1000,1000);  
  c->Draw();
  c->Divide(1,2);

  // From model:
  TF1 *ft = new TF1("ft"," al_fy_from_model(x)",-0.2,0.1); 
  TF1 *ft0 = new TF1("ft0"," al_fy_from_pass0(x)",-0.2,0.1); 

  ft->SetLineColor(kBlack);
  ft->SetLineStyle(1);
  ft0->SetLineColor(kMagenta);
  ft0->SetLineStyle(9);


  TMultiGraph * mg = new TMultiGraph("mg",Form("%s Al, Y scaling at %2.0f deg with and without DIS;y(GeV/c);F(y,q)(GeV/c)^{-1}",tgt_pos.Data(),angle));

  // draw Fy with and without dis subtraction
  c->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();

  // with DIS subtraction
  TGraphErrors* gr = new TGraphErrors(yv,fyv,yev,fyev); 
  gr->SetMarkerColor(kBlue);
  gr->SetLineStyle(1);
  gr->SetMarkerStyle(25);
  gr->SetMarkerSize(1);
  gr->Fit("myfunc","EBMR");
  TF1 *tfit = gr->GetFunction("myfunc");
  mg->Add(gr);

  // without dis subtraction
  TGraph* gr1 = new TGraph(yv,fyuncorrected);
  gr1->SetMarkerColor(kRed);
  gr1->SetLineStyle(1);
  gr1->SetMarkerStyle(24);
  gr1->SetMarkerSize(1);
  mg->Add(gr1);

  mg->Draw("AEP");
  mg->GetYaxis()->SetRangeUser(0.0,4);
  mg->GetXaxis()->SetRangeUser(-1.0,0.2);

  // draw the model
  ft->Draw("same");

  //draw the pass0 fit
  ///ft0->Draw("same");

  TLegend * leg1= new TLegend(0.2,0.6,0.6,0.8);
  leg1->AddEntry(gr1,"without DIS subtraction","p");
  leg1->AddEntry(gr,"with DIS subtraction","p");
  leg1->AddEntry(ft,"model","l");
  leg1->AddEntry(tfit,"data","l");

  //leg1->AddEntry(tfit,"pass1","l");
  leg1->Draw();

  
  // Draw FY with DIS subtraction but in Log scale
  c->cd(2);
  gPad->SetLogy();
  gPad->SetGridy();
  gPad->SetGridx();

  // with DIS subtraction
  TGraphErrors * gr2 = new TGraphErrors(yv,fyv,yev,fyev);
  gr2->SetMarkerColor(kBlue);
  gr2->SetLineStyle(2);
  gr2->SetMarkerStyle(25);
  gr2->SetMarkerSize(1);
  gr2->Draw("aep");
  gr2->GetXaxis()->SetRangeUser(-1.0,0.2);
  gr2->GetYaxis()->SetRangeUser(0.0,4.0);

  gr2->Fit("myfunc","QEBMR");
  TF1 *tfit1 = gr2->GetFunction("myfunc");


  gr2->GetXaxis()->SetTitle("y(GeV/c)");
  gr2->GetYaxis()->SetTitle("F(y,q)(GeV/c)^{-1}");
  gr2->SetTitle("");
  // draw model
  ft->Draw("same");
  //draw pass0 results
  //ft0->Draw("same");
 

  TLegend * leg2= new TLegend(0.2,0.65,0.5,0.8);
  leg2->AddEntry(gr2,"with DIS subtraction","p");
  leg2->AddEntry(ft,"model","l");
  //leg2->AddEntry(ft0,"pass0","l");
  leg2->AddEntry(tfit1,"data","l");

  leg2->Draw();

  //close the text file
  fclose(fp);

  // c->Print(Form("xsec_rad_exp_mod_ratio_x_geq_%2.1f_and_x_leq_%2.1f_%s_dummy_%s_%s.png",xlow,xup,tgt_pos.Data(),angle.Data(),fit.Data()));

  
}

Double_t myfunction(Double_t *y, Double_t *par)
{
  // From model:
  // f(y) = (f0-B)*alpha**2*exp(-(a*y)**2)/(alpha**2+y**2)+bigB*exp((-b*abs(y))^2) 
  //par0 = f0
  //par1 = B
  //par2 = alpha
  //par3 = a
  //par4 = b

  Float_t x =y[0];
  Double_t f = (par[0]-par[1])*(TMath::Power(par[2],2)*TMath::Exp(-1*TMath::Power((par[3]*x),2)))/(TMath::Power(par[2],2)+TMath::Power(x,2)) + par[1]*TMath::Exp(-1*TMath::Power((par[4]*TMath::Abs(x)),2));
 return f;
}


void myfunc()
{
  // The limits on the fit parameters are set 
  // according to the trend seen in the Fy of other nuclei

  TF1 *f1 = new TF1("myfunc",myfunction,-0.2,0.1,5);
  //f1->SetParameters(3.112,1.465,111E-3,3.944,5.784);
  f1->SetParameters(3.61,5.999,156E-3,4.93,6.0);


  f1->SetParName(0,"f0");
  f1->SetParLimits(0,0,8); 

  f1->SetParName(1,"B");
  f1->SetParLimits(1,0,6);

  f1->SetParName(2,"#alpha");
  f1->SetParLimits(2,0.0,2.0);
  
  f1->SetParName(3,"a");
  f1->SetParLimits(3,0,10);
 
  f1->SetParName(4,"b");
  f1->SetParLimits(4,5,10);

  f1->SetLineColor(kBlue);
  //f1->SetLineStyle(7);

}



Double_t calculate_fy_including_dis(Double_t sig, Double_t fact){
  /*
    For A>2 targets,
    fy(i) = (sigexp(i)-sigdis(i))/fact(i) in units of TeV^-1 (as outputed from dummy_ana.f)
  */
  return (sig/fact)*1e-3; //multiply by 1e-3 to get into units of GeV^-1 
}

Double_t c12_fy_from_model(Double_t x){

  // From the model:
  // f(y) = (f0-B)*alpha**2*exp(-(a*y)**2)/(alpha**2+y**2)+bigB*exp((-b*abs(y))^2) 
  // These values were obtained from the In the radcor_xem/target_info.f file 
  // According to the fy.f file the units of these fit parameters should be in MeV. I converted to GeV
  // from radiative correction file and Ajis thesis.
  Double_t f_0 = 3.20668; // GeV^-1
  Double_t B = 1.14392;// GeV^-1
  Double_t a = 3.41618;//GeV^-1 
  Double_t b = 6.63266;// GeV^-1 
  Double_t alpha = 134.148E-3;//  GeV 


  // // From J. Arrington's thesis page 160
  // Double_t f_0 = 3.3; // GeV^-1
  // Double_t B = 0.4; // GeV^-1
  // Double_t a = 3.88; // GeV^-1
  // Double_t b = 10; //GeV^-1
  // Double_t alpha = 0.140; // GeV^-1

  // From Nadia's thesis page 143
  // Double_t f_0 = 3.182e-03;
  // Double_t B = 1.359e-03;
  // Double_t a = 3.027e-03;
  // Double_t b =7.050e-03;
  // Double_t alpha = 137.2e-3;

  Double_t val = (f_0-B)*(TMath::Power(alpha,2)*TMath::Exp(-1*TMath::Power((a*x),2)))/(TMath::Power(alpha,2)+TMath::Power(x,2)) + B*TMath::Exp(-1*TMath::Power((b*TMath::Abs(x)),2));
 return val;
} 


Double_t al_fy_from_model(Double_t x){



  // From the model:
  // f(y) = (f0-B)*alpha**2*exp(-(a*y)**2)/(alpha**2+y**2)+bigB*exp((-b*abs(y))^2) 
  //These values were obtained from the radcor_xem/target_info.f file 
  Double_t f_0= 3.610000000000000;// GeV^-1
  Double_t B = 5.999999999999999;// GeV^-1
  Double_t a = 4.930000000000001;// GeV
  Double_t b = 6.000000000000000;// GeV^-1
  Double_t alpha =  156.000000000000E-003; // GeV determines the width of the peak

  // //from C12
  // Double_t f_0 = 3.20668; // GeV^-1
  // Double_t B = 1.14392;// GeV^-1
  // Double_t a = 3.41618;//GeV^-1 
  // Double_t b = 6.63266;// GeV^-1 
  // Double_t alpha = 134.148E-3;//  GeV 


  Double_t val = (f_0-B)*(TMath::Power(alpha,2)*TMath::Exp(-1*TMath::Power((a*x),2)))/(TMath::Power(alpha,2)+TMath::Power(x,2)) + B*TMath::Exp(-1*TMath::Power((b*TMath::Abs(x)),2));

  return val;
}


Double_t al_fy_from_pass0(Double_t x){

  Double_t f_0 = 0  ; // GeV^-1
  Double_t B = 0;// GeV^-1
  Double_t a = 0;//GeV^-1 
  Double_t b = 0;// GeV^-1 
  Double_t alpha = 0;//  GeV 


  // From the model:
  // f(y) = (f0-B)*alpha**2*exp(-(a*y)**2)/(alpha**2+y**2)+bigB*exp((-b*abs(y))^2) 


  // from the fit on pass0 18 deg US target in the y = -1.0 to 0.1 region
  if(upstream){
    f_0 = 3.101  ; // GeV^-1
    B = 1.113;// GeV^-1
    a = 3.637;//GeV^-1 
    b = 6.319;// GeV^-1 
    alpha = 119.0E-3;//  GeV 
  }
  else{
    //from the fit on the pass0 18 deg  DS target in the y = -1.0 to 0.1 region
    f_0 =3.038;
    B = 1.751;// GeV^-1
    a = 3.46;//GeV^-1 
    b = 6.421;// GeV^-1 
    alpha = 148.1E-3;//  GeV 
    
  }


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


