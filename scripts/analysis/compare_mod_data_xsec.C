
/*
A script to draw the model cross-section (tot, QE and DIS peices) and the data xsec for the 
Al dummy target from the RAD files.
This will allow us to see how good the comparision between the model and the data is.

Buddhini - 05-13-2015 

06-30-3016 add radiative correction plot

To run:
in root do-
.L compare_mod_data_xsec.C
compare_mod_data_xsec("40","up")
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

using namespace std;
//void myfunc();

void compare_mod_data_xsec(Double_t angle, TString tgt_pos){

  // histo parameters
  gStyle->SetTitleYOffset(1.0);
  gStyle->SetTitleXOffset(1.0);
  gStyle->SetLabelSize(0.07,"x");
  gStyle->SetLabelSize(0.07,"y");
  gStyle->SetTitleSize(0.07,"x");
  gStyle->SetTitleSize(0.07,"y");
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
  
  // select rad files of dummy tgt_pos with the given angle
  DIR *dir;
  struct dirent *ent;
  std::vector<string> filenames;
  if ((dir = opendir ("/home/buddhini/MyWork/HallC/SRC/E03103/radcor_xem/OUT/pass0")) != NULL) 
    {
      std::cout<<"Found the directory with rad files\n";
      
      /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL) 
	{
	  if((TString (ent->d_name)).Contains(Form("%2.0f",angle)) && (TString (ent->d_name)).Contains(tgt_pos))
	    {
	      // pass0 contains rad files generated with the QE model based on an Fy fit over the
	      // pass0 18 deg US and DS Fy seperately. 
	      filenames.push_back(Form("/home/buddhini/MyWork/HallC/SRC/E03103/radcor_xem/OUT/pass1_using_full_dataset/%s",(ent->d_name)));
	      //filenames.push_back(Form("/home/buddhini/MyWork/HallC/SRC/E03103/radcor_xem/OUT/pass0/%s",(ent->d_name)));

	      printf ("%s\n", ent->d_name);
	      std::cout<<ent->d_name<<std::endl;
	    }
	  else{
	    std::cout<<"unable to find rad file for "<<tgt_pos<<" target at angle "<<angle<<std::endl;
	    //  exit(1);
	  }
	}
      closedir (dir);
    } else 
    {

      /* could not open directory */
      perror ("");
    }


  // // Open a canvas
  TCanvas * c = new TCanvas("c", "",0,0,1600,800);  
  c->Draw();
 
  // Get the sig_rad and sig_dis, sig_qe and sigraw from the text file
  std::string line;
  std::string x_bjs;
  std::string thet;
  std::string sig_raw;
  std::string sig_qe;
  std::string sig_dis;
  std::string rcs;

  double theta = 0;
  // double xsec_dis = 0;
  // double xsec_qe = 0;
  // double xsec_raw = 0;

  TVectorD xsec_dis(0),xsec_qe(0),xsec_raw(0),x_bj(0);
  TVectorD rc(0);



  //get the size of the file array
  //int j = filenames.size();
  int k;
  int kv = 0;
  for(int i=0;i<1;i++)
    { 

      ifstream myfile ((char*)(filenames.at(i)).c_str());
      k=1;
      if (myfile.is_open())
	{

	  while ( getline (myfile,line) )
	    {
	      if(k<32) 
		{
		  //skip the first 31 lines
		  k++;
		}
	      else 
		{
		  k++;
		  //if(k==35) break;

	
		  //std::cout << line << '\n';
		  line.erase (1,15);
		  thet =line.substr (1,11);
		  theta = atof((char*)(thet.c_str()));
		  //std::cout << "theta = "<<theta << '\n';

		  if(theta == angle){
		    //std::cout << "theta = "<<theta << '\n';

		    x_bj.ResizeTo(kv+1);
		    xsec_dis.ResizeTo(kv+1);
		    xsec_qe.ResizeTo(kv+1);
		    xsec_raw.ResizeTo(kv+1);
		    rc.ResizeTo(kv+1);


		    line.erase (1,17);
		    //std::cout << line << '\n';
		    x_bjs = line.substr (1,11);
		    x_bj.operator()(kv) = atof((char*)(x_bjs.c_str()));
		    //std::cout << "xbj = "<<x_bj(kv) << '\n';
		    
		    line.erase (1,12);
		    //std::cout << line << '\n';
		    sig_dis = line.substr (1,12);
		    xsec_dis.operator()(kv) = atof((char*)(sig_dis.c_str()));
		    //std::cout << "sig_dis = "<<xsec_dis(k) << '\n';
		    
		    line.erase (1,12);
		    // std::cout << line << '\n';
		    sig_qe = line.substr (1,12);
		    xsec_qe.operator()(kv)  = atof((char*)(sig_qe.c_str()));
		    //std::cout << "sig_qe = "<<xsec_qe(k) << '\n';
		    
		    line.erase (1,12);
		    // std::cout << line << '\n';
		    sig_raw = line.substr (1,12);
		    xsec_raw.operator()(kv)  = atof((char*)(sig_raw.c_str()));
		    //std::cout << "sig_raw = "<<xsec_raw(k) << '\n';

		    line.erase (1,27);
		    //std::cout << line << '\n';
		    rcs = line.substr (1,12);
		    rc.operator()(kv)  = atof((char*)(rcs.c_str()));
		    //std::cout << "rc = "<<rc(kv) << '\n';

		    // if((xsec_dis(k)==0) || (xsec_qe==0)||(xsec_raw==0))
		    //   std::cout<<"zero!!\n";

		    kv++;
		  }
		}
	    }
	  myfile.close();
	  std::cout<<" read "<<kv<<" lines"<<std::endl;
	}
      else std::cout << "Unable to open file!\n"; 


    }
  // myfunc();


  TMultiGraph * mg = new TMultiGraph("mg",Form("%s Al dummy Model cross-section at %2.0f deg;x_{bj};#frac{d#sigma}{d#Omega dE}",tgt_pos.Data(),angle));

  TGraph* gr = new TGraph(x_bj,xsec_raw);
  gr->SetMarkerColor(kGreen);
  gr->SetLineStyle(2);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1);
  mg->Add(gr);

  TGraph* gq = new TGraph(x_bj,xsec_qe);
  gq->SetMarkerColor(kBlue);
  gq->SetLineStyle(2);
  gq->SetMarkerStyle(4);
  gq->SetMarkerSize(1);
  mg->Add(gq);

  TGraph* gd = new TGraph(x_bj,xsec_dis);
  gd->SetMarkerColor(kRed);
  gd->SetLineStyle(2);
  gd->SetMarkerStyle(4);
  gd->SetMarkerSize(1);
  mg->Add(gd);

  mg->Draw("AP");
  mg->GetXaxis()->SetRangeUser(0.2,2);
  //mg->GetYaxis()->SetRangeUser(0.2,2);

  TLegend * leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(gr,"Total","p");
  leg->AddEntry(gd,"DIS","p");
  leg->AddEntry(gq,"QE","p");
  leg->Draw();


  // Draw radiative corrections vs x
  TCanvas * c1 = new TCanvas("c1", "",0,0,1600,800);  
  c1->Draw();

  TGraph* grc = new TGraph(x_bj,rc);
  grc->SetMarkerColor(kGreen);
  grc->SetLineStyle(2);
  grc->SetMarkerStyle(20);
  grc->SetMarkerSize(1);
  grc->Draw("AEP");


  // // c->Print(Form("xsec_rad_exp_mod_ratio_x_geq_%2.1f_and_x_leq_%2.1f_%s_dummy_%s_%s.png",xlow,xup,tgt_pos.Data(),angle.Data(),fit.Data()));


}

