/*
A script to draw the radiated experimental cross section

To run:
in root do-
.L draw_cross_sections_from_data.C+g
draw_cross_sections(18,"c12")

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


Bool_t upstream = false;



void draw_cross_sections(Double_t angle, TString tgt_pos, TCanvas *c);


void draw_xsec_for_all(Double_t angle){
  TCanvas * c1 = new TCanvas("c1", "",0,0,800,500);  
  TCanvas * c2 = new TCanvas("c2", "",0,0,800,500);  
  std::cout<<"US target\n";
  upstream = true;
  draw_cross_sections(angle, "us",c1);

  std::cout<<"DS target\n";
  upstream = false;
  draw_cross_sections(angle, "ds",c2);

}


void draw_cross_sections(Double_t angle, TString tgt_pos,TCanvas *c){

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
  


 
  // Get sig_rad and esig from *.deltabin.xsec files
  std::string line;
  std::string ys;
  std::string fys;
  std::string efys;
  std::string q2s;
  std::string sigrads;
  std::string esigrads;
  std::string sigdiss;
  std::string facts;
  std::string xbjs;


  TVectorD fyv(0),yv(0),q2v(0);
  TVectorD fyev(0), yev(0);

  TVectorD sigdisv(0),sigradv(0),esigradv(0),factv(0),x_bj(0);
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
		  esigradv.ResizeTo(kv+1);
		  sigdisv.ResizeTo(kv+1);
		  factv.ResizeTo(kv+1);
		  fyuncorrected.ResizeTo(kv+1);
		  x_bj.ResizeTo(kv+1);

		  //std::cout << line << '\n';

		  line.erase (1,7);
		  xbjs =line.substr (1,7);
		  x_bj.operator()(kv) = atof((char*)(xbjs.c_str()));
		  //std::cout << "xbj = "<<xbjs << '\n';


		  line.erase (1,28);
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
		  
		  line.erase (1,12);
		  //std::cout << line << '\n';
		  facts = line.substr (1,12);
		  factv.operator()(kv) = atof((char*)(facts.c_str())); //units of ster/mB
		  //std::cout << " factor = "<<factv(kv)*0.0001<< '\n';

		  // dis subtraction.
		  // radiated experimental cross section
		  line.erase (1,12);
		  //std::cout << line << '\n';
		  sigrads = line.substr (1,12);
		  sigradv.operator()(kv) = atof((char*)(sigrads.c_str())); //units of nB(ster MeV)
		  //std::cout << "sig rad = "<<sigradv(kv)<< '\n';

		  
		  // dis cross section from the XEM model
		  line.erase (1,12);
		  //std::cout << line << '\n';
		  sigdiss = line.substr (1,12);
		  sigdisv.operator()(kv) = atof((char*)(sigdiss.c_str()));//units of nB(ster MeV)
		  //std::cout << "sig dis = "<<sigdisv(kv)<< '\n';

		  //error on the radiated cross-section
		  line.erase (1,24);
		  //std::cout << line << '\n';
		  esigrads = line.substr (1,12);
		  esigradv.operator()(kv) = atof((char*)(esigrads.c_str())); //units of nB(ster MeV)
		  //std::cout << "esig = "<<esigradv(kv)<< '\n';

		  kv++;

		}
	    }
	  myfile.close();
	}
      else std::cout << "Unable to open file :"<<filenames.at(i) <<std::endl; 


    }


  // Open canvas 
  //c = new TCanvas("c", "",0,0,600,500);  
  c->Draw();
  TMultiGraph * mg = new TMultiGraph("mg",Form("%s Al, cross-sections;x_{bj};#frac{d#sigma}{d#Omega dE}",tgt_pos.Data(),angle));

  // total cross sections
  TGraphErrors* gr = new TGraphErrors(x_bj,sigradv,yev,esigradv); 
  gr->SetMarkerColor(kBlue);
  gr->SetLineStyle(1);
  gr->SetMarkerStyle(25);
  gr->SetMarkerSize(1);
  // gr->Draw("AEP");
  // gr->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega dE} mB");
  // gr->GetXaxis()->SetTitle("x_{bj}");

  // TF1 *tfit = gr->GetFunction("myfunc");
  mg->Add(gr);

  // dis cross-section
  TGraph* gr1 = new TGraph(x_bj,sigdisv);
  gr1->SetMarkerColor(kRed);
  gr1->SetLineStyle(1);
  gr1->SetMarkerStyle(24);
  gr1->SetMarkerSize(1);
  mg->Add(gr1);

  mg->Draw("AEP");
  //mg->GetYaxis()->SetRangeUser(0.0,4);
  mg->GetXaxis()->SetRangeUser(-1.0,2.0);


  TLegend * leg= new TLegend(0.2,0.6,0.6,0.8);
  leg->AddEntry(gr1,"DIS","p");
  leg->AddEntry(gr,"total","p");
  leg->Draw();


  // //close the text file
  // fclose(fp);

  // c->Print(Form("xsec_rad_exp_mod_ratio_x_geq_%2.1f_and_x_leq_%2.1f_%s_dummy_%s_%s.png",xlow,xup,tgt_pos.Data(),angle.Data(),fit.Data()));

  
}

