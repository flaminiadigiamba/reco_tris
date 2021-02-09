//This in particular compile using  g++ Analyzer.cxx Spotsize.cxx -o spotsize.exe `root-config --libs --cflags`
//Then use as ./spotsize.exe path_to_rootfile

#include <iostream>
#include <string>
#include <vector>
#include "Analyzer.h"
#include <TTree.h>
#include <TFile.h>
#include <TArrayF.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TLine.h>
#include <TMath.h>

using namespace std;



void ScIndicesElem(int nSc,int* ScNelements, vector<int>* B, vector<int> *E){
  B->clear();
  E->clear();
  
  int parcount=0;
 
  for(int i=0;i<nSc;i++){
    B->push_back(parcount);
    E->push_back(parcount+ScNelements[i]);

    parcount+=ScNelements[i];
  }
}



int main(int argc, char** argv){

  ////////////////////////////////////Get File //////////////////////////////////////
  TFile* f = TFile::Open(Form("%s",argv[1]));
  TTree* tree = (TTree*)f->Get("Events");
  
  ///////////////////////////////////Set Branches and Define Variables////////////////////////////////////
  int nmax=250000;
  int nscmax=50;
  int npixel=2304;
  int npixelsmall=250;
  float slimnesslimit=0.6;
  
  unsigned int nSc;
  int run;
  int event;
             //Pixels
  
  int scID[nmax];
  int scIDall[nmax];
  int ScNpixels[nscmax];
  float XPix[nmax];
  float YPix[nmax];
  float ZPix[nmax];
  int ScNpixelsall[nscmax];
  float XPixall[nmax];
  float YPixall[nmax];
  float ZPixall[nmax];
  float width[nscmax];
  float length[nscmax];
  
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("event",&event);           
  tree->SetBranchAddress("nSc",&nSc);
  tree->SetBranchAddress("sc_ID",scID);
  tree->SetBranchAddress("sc_nintpixels",ScNpixels);
  tree->SetBranchAddress("sc_ypixelcoord",XPix);
  tree->SetBranchAddress("sc_xpixelcoord",YPix);
  tree->SetBranchAddress("sc_zpixel",ZPix);
  tree->SetBranchAddress("sc_IDall",scIDall);
  tree->SetBranchAddress("sc_nallintpixels",ScNpixelsall);
  tree->SetBranchAddress("sc_yallpixelcoord",XPixall);
  tree->SetBranchAddress("sc_xallpixelcoord",YPixall);
  tree->SetBranchAddress("sc_zallpixel",ZPixall);
  tree->SetBranchAddress("sc_width",width);
  tree->SetBranchAddress("sc_length",length);


  /////////////////////////////////Analysis Variables ////////////////////////////////////////////////
  vector<int> BeginScPix;
  vector<int> EndScPix;
  vector<int> BeginScallPix;
  vector<int> EndScallPix;

  
  /////////////////////////////////Analysis //////////////////////////////////////////////////////////
  int counter=0;
  int counterall=0;
  TH2F* Track=new TH2F("Track","Track",npixel,0,npixel,npixel,0,npixel);
  TH2F* Trackall=new TH2F("Trackall","Trackall",npixel,0,npixel,npixel,0,npixel);
  TFile* fout = new TFile(Form("Tracks.root"),"recreate");
  //fout->cd();
  for(int k=0;k<tree->GetEntries();k++)
  {
    tree->GetEntry(k);
    cout << "Nev: "<< k << endl;
    
    ScIndicesElem(nSc,ScNpixels,&BeginScPix,&EndScPix);
    ScIndicesElem(nSc,ScNpixelsall,&BeginScallPix,&EndScallPix);
	
	//Start the cycle on the supercluster of the event
    for(int i=0;i<nSc;i++)
    {   
	   if(scID[BeginScPix[i]]>=0)
	   {
		   if(width[i]/length[i]>0.6 && width[i]/length[i]<1)
		   {
			for(int j=BeginScPix[i];j<EndScPix[i];j++)
			{
				Track->SetBinContent(XPix[j],YPix[j],ZPix[j]);
			}  //chiudo for j (fill histos)
			
			Track->SetName(Form("Track%i_run%i_evt%i",counter,run,event));
			Track->SetTitle(Form("Track%i_run%i_evt%i",counter,run,event));
			//Track->Write();
			
			Analyzer Traccia("Track",npixelsmall,Track,npixel);
			Traccia.SavetoFile(Form("Track_small%i_run%i_evt%i",counter,run,event));
			
			//test for the centre and radius
			TCanvas *c1=new TCanvas("c1","c1",1000,1000);
			TH2F *test=Traccia.GetHisto();
			double x,y;
			Traccia.GetCentre(x,y);
			TLine *ltest=new TLine(x,y,x+Traccia.GetRadius()/sqrt(2),y+Traccia.GetRadius()/sqrt(2));
			cout<<x<<"   "<<y<<"     "<<Traccia.GetRadius()<<endl;
			ltest->SetLineWidth(2);
			ltest->SetLineColor(kRed);
			test->Draw("colz");
			ltest->Draw("same");
			c1->Write();
			
			counter++;
			Track->Reset();
			//cout<<counter<<endl;
		   }	
		}		//end if on overthreshold sc
	   
	   if(scIDall[BeginScallPix[i]]>=0)
	   {
		   if(width[i]/length[i]>0.6 && width[i]/length[i]<1)
		   {
			for(int j=BeginScallPix[i];j<EndScallPix[i];j++)
			{
				Trackall->SetBinContent(XPixall[j],YPixall[j],ZPixall[j]);
			}//chiudo for j (fill histos)
			Trackall->SetName(Form("Trackall%i_run%i_evt%i",counter,run,event));
			Trackall->SetTitle(Form("Trackall%i_run%i_evt%i",counter,run,event));
			//Trackall->Write();
			Analyzer Traccia("Trackall",npixelsmall,Trackall,npixel);
			Traccia.SavetoFile(Form("Trackall_small%i_run%i_evt%i",counterall,run,event));
			
			
			counterall++;
			Trackall->Reset();
			//cout<<counter<<endl;
		   }	
		}		//end if on allhit sc
	
	}	//end if on supercluster number
	
  }//chiudo for k (loop on the events)
  

////////////////////////Save Plot On File////////////////////////////////////////

      fout->Save();
      fout->Close();
    
    delete Track;
  
  return 0;
}

