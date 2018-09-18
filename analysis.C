#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TVector2.h>
#include <TMath.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <TLorentzVector.h>

#include "HttStyles.h"

#include "helperFxns.hh"
#include "drawFxns.hh"

using namespace std;

//not great practice but ok
const Double_t MUON_MASS = 0.105658369;
const Double_t ELE_MASS =  0.000511;
const Double_t Z_MASS =   91.1876;
const Float_t LUMI=3000;

//Float_t selection(const TString inputfile, HistHolder hists, float &l3, float &l2, float &l1, float &l0, float &errorl3, float &errorl2, float &errorl1, float &errorl0, float &evtWeight);
Float_t selection(const TString inputfile, HistHolder hists, int &l3, int &l2, int &l1, int &l0, float &evtWeight);
//Float_t selection(const TString inputfile, HistHolder hists, int &lll, int &eej, int &emj, int &mmj, int &ejj, int &mjj, int &jjj, float &evtWeight, TDirectory *cdir[7]);
//Float_t selection(const TString inputfile, HistHolder hists, int &lll, int &eej, int &emj, int &mmj, int &ejj, int &mjj, int &jjj, float &evtWeight);

void analysis() {
  
  HistHolder signal; initHistHolder(signal, "sig");
  HistHolder ttbar;  initHistHolder(ttbar, "tt");
  HistHolder otherb; initHistHolder(otherb, "bgd");
  
  int l3[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int l2[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int l1[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int l0[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  float scale[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  
  /*
  int lll[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int eej[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int emj[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int mmj[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int ejj[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int mjj[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int jjj[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  TString category[7] = {"lll", "eej", "emj", "mmj", "ejj", "mjj", "jjj"}; //category name, like l3 or l2 etc
  TFile *f = new TFile("test_input.root","recreate");//new file
  TDirectory *cdir[7];
  for (int i=0; i<7; i++)
    {
      cdir[i] = f->mkdir(category[i]);// directory by the name of the category
      cdir[i]->cd();//switch to that directory
      TH1D *data_obs = new TH1D("data_obs","data_obs",30,230,1500); //"data", has to be called "data_obs"
      TH1D *ZZZ = new TH1D("ZZZ","ZZZ",30,230,1500);//signal
      TH1D *bkg = new TH1D("bkg","bkg",30,230,1500);
      f->cd();
    }
  f->Write();
  
  
  float zzz=selection("events2/ZZZ.root", signal, lll[0], eej[0], emj[0], mmj[0], ejj[0], mjj[0], jjj[0], scale[0], cdir);
  float ttz=selection("events2/TTZ.root", otherb, lll[1], eej[1], emj[1], mmj[1], ejj[1], mjj[1], jjj[1], scale[1], cdir);
  float tt =selection("events2/TT_350.root", otherb, lll[2], eej[2], emj[2], mmj[2], ejj[2], mjj[2], jjj[2], scale[2], cdir);
  float www=selection("events2/WWW.root", otherb, lll[3], eej[3], emj[3], mmj[3], ejj[3], mjj[3], jjj[3], scale[3], cdir);
  float wzz=selection("events2/WZZ.root", otherb, lll[4], eej[4], emj[4], mmj[4], ejj[4], mjj[4], jjj[4], scale[4], cdir);
  float wwz=selection("events2/WWZ.root", otherb, lll[5], eej[5], emj[5], mmj[5], ejj[5], mjj[5], jjj[5], scale[5], cdir);
  float ww =selection("events2/WW.root", otherb, lll[6], eej[6], emj[6], mmj[6], ejj[6], mjj[6], jjj[6], scale[6], cdir);
  float zz2l =selection("events2/ZZTo2L2Q.root", otherb, lll[7], eej[7], emj[7], mmj[7], ejj[7], mjj[7], jjj[7], scale[7], cdir);
  float hzz=selection("events2/HZZ.root", otherb, lll[8], eej[8], emj[8], mmj[8], ejj[8], mjj[8], jjj[8], scale[8], cdir);
  float zhzz=selection("events2/ZHToZZ.root", otherb, lll[9], eej[9], emj[9], mmj[9], ejj[9], mjj[9], jjj[9], scale[9], cdir);
  float zz4l=selection("events2/ZZTo4L.root", ttbar, lll[10], eej[10], emj[10], mmj[10], ejj[10], mjj[10], jjj[10], scale[10], cdir);
  // calculate the effective evtWeight for bkg in each channel
  float sample[7] = {0, 0, 0, 0, 0, 0, 0};
  float yield[7] = {0, 0, 0, 0, 0, 0, 0};
  for (int i = 1; i < 11; i++)
    {
      sample[0] += lll[i];
      yield[0] += lll[i] * scale[i];
      sample[1] += eej[i];
      yield[1] += eej[i] * scale[i];
      sample[2] += emj[i];
      yield[2] += emj[i] * scale[i];
      sample[3] += mmj[i];
      yield[3] += mmj[i] * scale[i];
      sample[4] += ejj[i];
      yield[4] += ejj[i] * scale[i];
      sample[5] += mjj[i];
      yield[5] += mjj[i] * scale[i];
      sample[6] += jjj[i];
      yield[6] += jjj[i] * scale[i];
    }
  float weight[7] = {0, 0, 0, 0, 0, 0, 0};
  for (int i=0; i<7; i++)
    weight[i] = yield[i] / sample[i];
  int s[7] = {0, 0, 0, 0, 0, 0, 0};
  for (int i=1; i<11;i++)
    {
      s[0] += lll[i];
      s[1] += eej[i];
      s[2] += emj[i];
      s[3] += mmj[i];
      s[4] += ejj[i];
      s[5] += mjj[i];
      s[6] += jjj[i];
    }
  */
  /*
  float zzz=selection("events2/ZZZ.root", signal, lll[0], eej[0], emj[0], mmj[0], ejj[0], mjj[0], jjj[0], scale[0]);
  float ttz=selection("events2/TTZ.root", otherb, lll[1], eej[1], emj[1], mmj[1], ejj[1], mjj[1], jjj[1], scale[1]);
  float tt =selection("events2/TT_350.root", otherb, lll[2], eej[2], emj[2], mmj[2], ejj[2], mjj[2], jjj[2], scale[2]);
  float www=selection("events2/WWW.root", otherb, lll[3], eej[3], emj[3], mmj[3], ejj[3], mjj[3], jjj[3], scale[3]);
  float wzz=selection("events2/WZZ.root", otherb, lll[4], eej[4], emj[4], mmj[4], ejj[4], mjj[4], jjj[4], scale[4]);
  float wwz=selection("events2/WWZ.root", otherb, lll[5], eej[5], emj[5], mmj[5], ejj[5], mjj[5], jjj[5], scale[5]);
  float ww =selection("events2/WW.root", otherb, lll[6], eej[6], emj[6], mmj[6], ejj[6], mjj[6], jjj[6], scale[6]);
  float zz2l =selection("events2/ZZTo2L2Q.root", otherb, lll[7], eej[7], emj[7], mmj[7], ejj[7], mjj[7], jjj[7], scale[7]);
  float hzz=selection("events2/HZZ.root", otherb, lll[8], eej[8], emj[8], mmj[8], ejj[8], mjj[8], jjj[8], scale[8]);
  float zhzz=selection("events2/ZHToZZ.root", otherb, lll[9], eej[9], emj[9], mmj[9], ejj[9], mjj[9], jjj[9], scale[9]);
  float zz4l=selection("events2/ZZTo4L.root", ttbar, lll[10], eej[10], emj[10], mmj[10], ejj[10], mjj[10], jjj[10], scale[10]);
  */
  /*
  float zzz=selection("events2/ZZZ.root", signal, l3[0], l2[0], l1[0], l0[0], scale[0]);
  float ttz=selection("events2/TTZ.root", otherb, l3[1], l2[1], l1[1], l0[1], scale[1]);
  float tt =selection("events2/TT_350.root", otherb, l3[2], l2[2], l1[2], l0[2], scale[2]);
  float www=selection("events2/WWW.root", otherb, l3[3], l2[3], l1[3], l0[3], scale[3]);
  float wzz=selection("events2/WZZ.root", otherb, l3[4], l2[4], l1[4], l0[4], scale[4]);
  float wwz=selection("events2/WWZ.root", otherb, l3[5], l2[5], l1[5], l0[5], scale[5]);
  float ww =selection("events2/WW.root", ttbar, l3[6], l2[6], l1[6], l0[6], scale[6]);
  float zz2l =selection("events2/ZZTo2L2Q.root", otherb, l3[7], l2[7], l1[7], l0[7], scale[7]);
  float hzz=selection("events2/HZZ.root", otherb, l3[8], l2[8], l1[8], l0[8], scale[8]);
  //float zhzz=selection("events2/ZHToZZ.root", otherb, l3[9], l2[9], l1[9], l0[9], scale[9]);
  float zz4l=selection("events2/ZZTo4L.root", otherb, l3[10], l2[10], l1[10], l0[10], scale[10]);
  */
  
  float zzz=selection("events4/ZZZ.root", signal, l3[0], l2[0], l1[0], l0[0], scale[0]);
  float ttx=selection("events4/TTXX_comb.root", otherb, l3[1], l2[1], l1[1], l0[1], scale[1]);
  float tt =selection("events4/TT.root", ttbar, l3[2], l2[2], l1[2], l0[2], scale[2]);
  float www=selection("events4/WWW.root", otherb, l3[3], l2[3], l1[3], l0[3], scale[3]);
  float wzz=selection("events4/WZZ.root", otherb, l3[4], l2[4], l1[4], l0[4], scale[4]);
  float wwz=selection("events4/WWZ.root", otherb, l3[5], l2[5], l1[5], l0[5], scale[5]);
  float ww =selection("events4/WW.root", otherb, l3[6], l2[6], l1[6], l0[6], scale[6]);
  float ssww =selection("events4/ssWW_comb.root", otherb, l3[7], l2[7], l1[7], l0[7], scale[7]);
  float higgs =selection("events4/Higgs_comb.root", otherb, l3[8], l2[8], l1[8], l0[8], scale[8]);
  float zz2l =selection("events4/ZZTo2L2Q.root", otherb, l3[9], l2[9], l1[9], l0[9], scale[9]);
  float zz4l=selection("events4/ZZTo4L.root", otherb, l3[10], l2[10], l1[10], l0[10], scale[10]);
  float wz=selection("events4/WZ3LNu_0j.root", otherb, l3[11], l2[11], l1[11], l0[11], scale[11]);
  float dy=selection("events4/DYJets_comb.root", otherb, l3[12], l2[12], l1[12], l0[12], scale[12]);

  cout << "expected yields --------------------------" << endl;
  for (int i=0; i<13; i++)
    {
      cout << scale[i]*l3[i] << " " << scale[i]*l2[i] << " " << scale[i]*l1[i] << " " << scale[i]*l0[i] << endl;
    }
  
  
  // output the uncertainties of the expected yields
  cout << "systematic uncertainties ------------------------------" << endl;
  for (int i=0; i<13; i++)
    {
      //systematic uncertainties
      float u3 = scale[i]*l3[i]*TMath::Sqrt(0.01*0.01+6*0.005*0.005);
      float u2 = scale[i]*l2[i]*TMath::Sqrt(0.01*0.01+4*0.005*0.005+2*0.022*0.022);
      float u1 = scale[i]*l1[i]*TMath::Sqrt(0.01*0.01+2*0.005*0.005+4*0.022*0.022);
      float u0 = scale[i]*l0[i]*TMath::Sqrt(0.01*0.01+6*0.022*0.022);
      cout << u3 << " " << u2 << " " << u1 << " " << u0 << " " << TMath::Sqrt(u3*u3+u2*u2+u1*u1+u0*u0) << endl;    
    }
  cout << "statistical uncertainties ------------------------------" << endl;
  for (int i=0; i<13; i++)
    {   
      //MC uncertainties
      int temp = l3[i] + l2[i] + l1[i] + l0[i];
      cout << scale[i]*l3[i]/TMath::Sqrt(l3[i]+1) << " " << scale[i]*l2[i]/TMath::Sqrt(l2[i]+1) << " " << scale[i]*l1[i]/TMath::Sqrt(l1[i]+1) << " " << scale[i]*l0[i]/TMath::Sqrt(l0[i]+1) << " " << scale[i]*temp/TMath::Sqrt(temp+1) << endl;     
    }
  
  
  /*
  cdir[0]->cd();
  TH1D *ZZZ = (TH1D*) cdir[0]->Get("ZZZ");
  TH1D *data_obs = (TH1D*) cdir[0]->Get("data_obs");
  TH1D *bkg = (TH1D*) cdir[0]->Get("bkg");
  data_obs->Add(ZZZ); 
  data_obs->Add(bkg);
  f->Write();
  ofstream cout("test_datacard_lll.txt");
  cout << "max 1" << endl;
  cout << "jmax 1" << endl;
  cout << "kmax *" << endl;
  cout << "--------------------------" << endl;
  cout << "shapes * * test_input.root $CHANNEL/$PROCESS" << endl;
  cout << "--------------------------" << endl;
  cout << "bin         " << category[0] << endl;
  cout << "observation " << data_obs->Integral() << endl;
  cout << "-------------------------" << endl;
  cout << "bin " << category[0] << " " << category[0] << endl;
  cout << "process         0               1" << endl;
  cout << "process         ZZZ             bkg" << endl;
  cout << "rate " << ZZZ->Integral() << " " << bkg->Integral() << endl; //get the rate from each histogram
  cout << "--------------------------" << endl;
  cout << "MC00 gmN " << lll[0] << " " << scale[0] << " -" << endl;
  cout << "MC01 gmN " << s[0] << " - " << weight[0] << endl;
  cout << "MC01 gmN " << lll[1] << " - " << scale[1] << endl;
  cout << "MC02 gmN " << lll[2] << " - " << scale[2] << endl;
  cout << "MC03 gmN " << lll[3] << " - " << scale[3] << endl;
  cout << "MC04 gmN " << lll[4] << " - " << scale[4] << endl;
  cout << "MC05 gmN " << lll[5] << " - " << scale[5] << endl;
  cout << "MC06 gmN " << lll[6] << " - " << scale[6] << endl;
  cout << "MC07 gmN " << lll[7] << " - " << scale[7] << endl;
  cout << "MC08 gmN " << lll[8] << " - " << scale[8] << endl;
  cout << "MC09 gmN " << lll[9] << " - " << scale[9] << endl;
  cout << "MC010 gmN " << lll[10] << " - " << scale[10] << endl;
  cout << "lumi lnN 1.01 1.01" << endl;
  cout << "e lnN " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << endl;
  cout << "--------------------------" << endl;
  f->cd();
  cdir[1]->cd();
  ZZZ = (TH1D*) cdir[1]->Get("ZZZ");
  data_obs = (TH1D*) cdir[1]->Get("data_obs");
  bkg = (TH1D*) cdir[1]->Get("bkg");
  data_obs->Add(ZZZ); 
  data_obs->Add(bkg);
  f->Write();
  ofstream cout1("test_datacard_eej.txt");
  cout1 << "max 1" << endl;
  cout1 << "jmax 1" << endl;
  cout1 << "kmax *" << endl;
  cout1 << "--------------------------" << endl;
  cout1 << "shapes * * test_input.root $CHANNEL/$PROCESS" << endl;
  cout1 << "--------------------------" << endl;
  cout1 << "bin         " << category[1] << endl;
  cout1 << "observation " << data_obs->Integral() << endl;
  cout1 << "-------------------------" << endl;
  cout1 << "bin " << category[1] << " " << category[1] << endl;
  cout1 << "process         0               1" << endl;
  cout1 << "process         ZZZ             bkg" << endl;
  cout1 << "rate " << ZZZ->Integral() << " " << bkg->Integral() << endl; //get the rate from each histogram
  cout1 << "--------------------------" << endl;
  cout1 << "MC10 gmN " << eej[0] << " " << scale[0] << " -" << endl;
  cout1 << "MC11 gmN " << s[1] << " - " << weight[1] << endl;
  cout1 << "MC11 gmN " << eej[1] << " - " << scale[1] << endl;
  cout1 << "MC12 gmN " << eej[2] << " - " << scale[2] << endl;
  cout1 << "MC13 gmN " << eej[3] << " - " << scale[3] << endl;
  cout1 << "MC14 gmN " << eej[4] << " - " << scale[4] << endl;
  cout1 << "MC15 gmN " << eej[5] << " - " << scale[5] << endl;
  cout1 << "MC16 gmN " << eej[6] << " - " << scale[6] << endl;
  cout1 << "MC17 gmN " << eej[7] << " - " << scale[7] << endl;
  cout1 << "MC18 gmN " << eej[8] << " - " << scale[8] << endl;
  cout1 << "MC19 gmN " << eej[9] << " - " << scale[9] << endl;
  cout1 << "MC110 gmN " << eej[10] << " - " << scale[10] << endl;
  cout1 << "lumi lnN 1.01 1.01" << endl;
  cout1 << "e lnN " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << endl;
  cout1 << "j lnN " << 1+sqrt(2*(0.02)*(0.02)) << " " << 1+sqrt(2*(0.02)*(0.02)) << endl;
  cout1 << "--------------------------" << endl;
  f->cd();
  cdir[2]->cd();
  ZZZ = (TH1D*) cdir[2]->Get("ZZZ");
  data_obs = (TH1D*) cdir[2]->Get("data_obs");
  bkg = (TH1D*) cdir[2]->Get("bkg");
  data_obs->Add(ZZZ); 
  data_obs->Add(bkg);
  f->Write();
  ofstream cout2("test_datacard_emj.txt");
  cout2 << "max 1" << endl;
  cout2 << "jmax 1" << endl;
  cout2 << "kmax *" << endl;
  cout2 << "--------------------------" << endl;
  cout2 << "shapes * * test_input.root $CHANNEL/$PROCESS" << endl;
  cout2 << "--------------------------" << endl;
  cout2 << "bin         " << category[2] << endl;
  cout2 << "observation " << data_obs->Integral() << endl;
  cout2 << "-------------------------" << endl;
  cout2 << "bin " << category[2] << " " << category[2] << endl;
  cout2 << "process         0               1" << endl;
  cout2 << "process         ZZZ             bkg" << endl;
  cout2 << "rate " << ZZZ->Integral() << " " << bkg->Integral() << endl; //get the rate from each histogram
  cout2 << "--------------------------" << endl;
  cout2 << "MC20 gmN " << emj[0] << " " << scale[0] << " -" << endl;
  cout2 << "MC21 gmN " << s[2] << " - " << weight[2] << endl;
  cout2 << "MC21 gmN " << emj[1] << " - " << scale[1] << endl;
  cout2 << "MC22 gmN " << emj[2] << " - " << scale[2] << endl;
  cout2 << "MC23 gmN " << emj[3] << " - " << scale[3] << endl;
  cout2 << "MC24 gmN " << emj[4] << " - " << scale[4] << endl;
  cout2 << "MC25 gmN " << emj[5] << " - " << scale[5] << endl;
  cout2 << "MC26 gmN " << emj[6] << " - " << scale[6] << endl;
  cout2 << "MC27 gmN " << emj[7] << " - " << scale[7] << endl;
  cout2 << "MC28 gmN " << emj[8] << " - " << scale[8] << endl;
  cout2 << "MC29 gmN " << emj[9] << " - " << scale[9] << endl;
  cout2 << "MC210 gmN " << emj[10] << " - " << scale[10] << endl;
  cout2 << "lumi lnN 1.01 1.01" << endl;
  cout2 << "e lnN " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << endl;
  cout2 << "m lnN " << 1+sqrt(2*(0.005)*(0.005)) << " " << 1+sqrt(2*(0.005)*(0.005)) << endl;
  cout2 << "j lnN " << 1+sqrt(2*(0.02)*(0.02)) << " " << 1+sqrt(2*(0.02)*(0.02)) << endl;
  cout2 << "--------------------------" << endl;
  f->cd();
  cdir[3]->cd();
  ZZZ = (TH1D*) cdir[3]->Get("ZZZ");
  data_obs = (TH1D*) cdir[3]->Get("data_obs");
  bkg = (TH1D*) cdir[3]->Get("bkg");
  data_obs->Add(ZZZ); 
  data_obs->Add(bkg);
  f->Write();
  ofstream cout3("test_datacard_mmj.txt");
  cout3 << "max 1" << endl;
  cout3 << "jmax 1" << endl;
  cout3 << "kmax *" << endl;
  cout3 << "--------------------------" << endl;
  cout3 << "shapes * * test_input.root $CHANNEL/$PROCESS" << endl;
  cout3 << "--------------------------" << endl;
  cout3 << "bin         " << category[3] << endl;
  cout3 << "observation " << data_obs->Integral() << endl;
  cout3 << "-------------------------" << endl;
  cout3 << "bin " << category[3] << " " << category[3] << endl;
  cout3 << "process         0               1" << endl;
  cout3 << "process         ZZZ             bkg" << endl;
  cout3 << "rate " << ZZZ->Integral() << " " << bkg->Integral() << endl; //get the rate from each histogram
  cout3 << "--------------------------" << endl;
  cout3 << "MC30 gmN " << mmj[0] << " " << scale[0] << " -" << endl;
  cout3 << "MC31 gmN " << s[3] << " - " << weight[3] << endl;
  cout3 << "MC31 gmN " << mmj[1] << " - " << scale[1] << endl;
  cout3 << "MC32 gmN " << mmj[2] << " - " << scale[2] << endl;
  cout3 << "MC33 gmN " << mmj[3] << " - " << scale[3] << endl;
  cout3 << "MC34 gmN " << mmj[4] << " - " << scale[4] << endl;
  cout3 << "MC35 gmN " << mmj[5] << " - " << scale[5] << endl;
  cout3 << "MC36 gmN " << mmj[6] << " - " << scale[6] << endl;
  cout3 << "MC37 gmN " << mmj[7] << " - " << scale[7] << endl;
  cout3 << "MC38 gmN " << mmj[8] << " - " << scale[8] << endl;
  cout3 << "MC39 gmN " << mmj[9] << " - " << scale[9] << endl;
  cout3 << "MC310 gmN " << mmj[10] << " - " << scale[10] << endl;
  cout3 << "lumi lnN 1.01 1.01" << endl;
  cout3 << "m lnN " << 1+sqrt(4*(0.005)*(0.005)) << " " << 1+sqrt(4*(0.005)*(0.005)) << endl;
  cout3 << "j lnN " << 1+sqrt(2*(0.02)*(0.02)) << " " << 1+sqrt(2*(0.02)*(0.02)) << endl;
  cout3 << "--------------------------" << endl;
  f->cd();
  cdir[4]->cd();
  ZZZ = (TH1D*) cdir[4]->Get("ZZZ");
  data_obs = (TH1D*) cdir[4]->Get("data_obs");
  bkg = (TH1D*) cdir[4]->Get("bkg");
  data_obs->Add(ZZZ); 
  data_obs->Add(bkg);
  f->Write();
  ofstream cout4("test_datacard_ejj.txt");
  cout4 << "max 1" << endl;
  cout4 << "jmax 1" << endl;
  cout4 << "kmax *" << endl;
  cout4 << "--------------------------" << endl;
  cout4 << "shapes * * test_input.root $CHANNEL/$PROCESS" << endl;
  cout4 << "--------------------------" << endl;
  cout4 << "bin         " << category[4] << endl;
  cout4 << "observation " << data_obs->Integral() << endl;
  cout4 << "-------------------------" << endl;
  cout4 << "bin " << category[4] << " " << category[4] << endl;
  cout4 << "process         0               1" << endl;
  cout4 << "process         ZZZ             bkg" << endl;
  cout4 << "rate " << ZZZ->Integral() << " " << bkg->Integral() << endl; //get the rate from each histogram
  cout4 << "--------------------------" << endl;
  cout4 << "MC40 gmN " << ejj[0] << " " << scale[0] << " -" << endl;
  cout4 << "MC41 gmN " << s[4] << " - " << weight[4] << endl;
  cout4 << "MC41 gmN " << ejj[1] << " - " << scale[1] << endl;
  cout4 << "MC42 gmN " << ejj[2] << " - " << scale[2] << endl;
  cout4 << "MC43 gmN " << ejj[3] << " - " << scale[3] << endl;
  cout4 << "MC44 gmN " << ejj[4] << " - " << scale[4] << endl;
  cout4 << "MC45 gmN " << ejj[5] << " - " << scale[5] << endl;
  cout4 << "MC46 gmN " << ejj[6] << " - " << scale[6] << endl;
  cout4 << "MC47 gmN " << ejj[7] << " - " << scale[7] << endl;
  cout4 << "MC48 gmN " << ejj[8] << " - " << scale[8] << endl;
  cout4 << "MC49 gmN " << ejj[9] << " - " << scale[9] << endl;
  cout4 << "MC410 gmN " << ejj[10] << " - " << scale[10] << endl;
  cout4 << "lumi lnN 1.01 1.01" << endl;
  cout4 << "e lnN " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << endl;
  cout4 << "j lnN " << 1+sqrt(4*(0.02)*(0.02)) << " " << 1+sqrt(4*(0.02)*(0.02)) << endl;
  cout4 << "--------------------------" << endl;
  f->cd();
  cdir[5]->cd();
  ZZZ = (TH1D*) cdir[5]->Get("ZZZ");
  data_obs = (TH1D*) cdir[5]->Get("data_obs");
  bkg = (TH1D*) cdir[5]->Get("bkg");
  data_obs->Add(ZZZ); 
  data_obs->Add(bkg);
  f->Write();
  ofstream cout5("test_datacard_mjj.txt");
  cout5 << "max 1" << endl;
  cout5 << "jmax 1" << endl;
  cout5 << "kmax *" << endl;
  cout5 << "--------------------------" << endl;
  cout5 << "shapes * * test_input.root $CHANNEL/$PROCESS" << endl;
  cout5 << "--------------------------" << endl;
  cout5 << "bin         " << category[5] << endl;
  cout5 << "observation " << data_obs->Integral() << endl;
  cout5 << "-------------------------" << endl;
  cout5 << "bin " << category[5] << " " << category[5] << endl;
  cout5 << "process         0               1" << endl;
  cout5 << "process         ZZZ             bkg" << endl;
  cout5 << "rate " << ZZZ->Integral() << " " << bkg->Integral() << endl; //get the rate from each histogram
  cout5 << "--------------------------" << endl;
  cout5 << "MC50 gmN " << mjj[0] << " " << scale[0] << " -" << endl;
  cout5 << "MC51 gmN " << s[5] << " - " << weight[5] << endl;
  cout5 << "MC51 gmN " << mjj[1] << " - " << scale[1] << endl;
  cout5 << "MC52 gmN " << mjj[2] << " - " << scale[2] << endl;
  cout5 << "MC53 gmN " << mjj[3] << " - " << scale[3] << endl;
  cout5 << "MC54 gmN " << mjj[4] << " - " << scale[4] << endl;
  cout5 << "MC55 gmN " << mjj[5] << " - " << scale[5] << endl;
  cout5 << "MC56 gmN " << mjj[6] << " - " << scale[6] << endl;
  cout5 << "MC57 gmN " << mjj[7] << " - " << scale[7] << endl;
  cout5 << "MC58 gmN " << mjj[8] << " - " << scale[8] << endl;
  cout5 << "MC59 gmN " << mjj[9] << " - " << scale[9] << endl;
  cout5 << "MC510 gmN " << mjj[10] << " - " << scale[10] << endl;
  cout5 << "lumi lnN 1.01 1.01" << endl;
  cout5 << "m lnN " << 1+sqrt(2*(0.005)*(0.005)) << " " << 1+sqrt(2*(0.005)*(0.005)) << endl;
  cout5 << "j lnN " << 1+sqrt(4*(0.02)*(0.02)) << " " << 1+sqrt(4*(0.02)*(0.02)) << endl;
  cout5 << "--------------------------" << endl;
  f->cd();
  cdir[6]->cd();
  ZZZ = (TH1D*) cdir[6]->Get("ZZZ");
  data_obs = (TH1D*) cdir[6]->Get("data_obs");
  bkg = (TH1D*) cdir[6]->Get("bkg");
  data_obs->Add(ZZZ); 
  data_obs->Add(bkg);
  f->Write();
  ofstream cout6("test_datacard_jjj.txt");
  cout6 << "max 1" << endl;
  cout6 << "jmax 1" << endl;
  cout6 << "kmax *" << endl;
  cout6 << "--------------------------" << endl;
  cout6 << "shapes * * test_input.root $CHANNEL/$PROCESS" << endl;
  cout6 << "--------------------------" << endl;
  cout6 << "bin         " << category[6] << endl;
  cout6 << "observation " << data_obs->Integral() << endl;
  cout6 << "-------------------------" << endl;
  cout6 << "bin " << category[6] << " " << category[6] << endl;
  cout6 << "process         0               1" << endl;
  cout6 << "process         ZZZ             bkg" << endl;
  cout6 << "rate " << ZZZ->Integral() << " " << bkg->Integral() << endl; //get the rate from each histogram
  cout6 << "--------------------------" << endl;
  cout6 << "MC60 gmN " << jjj[0] << " " << scale[0] << " -" << endl;
  cout6 << "MC61 gmN " << s[6] << " - " << weight[6] << endl;
  cout6 << "MC61 gmN " << jjj[1] << " - " << scale[1] << endl;
  cout6 << "MC62 gmN " << jjj[2] << " - " << scale[2] << endl;
  cout6 << "MC63 gmN " << jjj[3] << " - " << scale[3] << endl;
  cout6 << "MC64 gmN " << jjj[4] << " - " << scale[4] << endl;
  cout6 << "MC65 gmN " << jjj[5] << " - " << scale[5] << endl;
  cout6 << "MC66 gmN " << jjj[6] << " - " << scale[6] << endl;
  cout6 << "MC67 gmN " << jjj[7] << " - " << scale[7] << endl;
  cout6 << "MC68 gmN " << jjj[8] << " - " << scale[8] << endl;
  cout6 << "MC69 gmN " << jjj[9] << " - " << scale[9] << endl;
  cout6 << "MC610 gmN " << jjj[10] << " - " << scale[10] << endl;
  cout6 << "lumi lnN 1.01 1.01" << endl;
  cout6 << "j lnN " << 1+sqrt(6*(0.02)*(0.02)) << " " << 1+sqrt(6*(0.02)*(0.02)) << endl;
  cout6 << "--------------------------" << endl;
  f->cd();
  
  f->Close();
  */
  /*
  float zzz=selection("events3/ZZZ.root", signal, l3[0], l2[0], l1[0], l0[0], scale[0]);
  float ttz=selection("events3/TTZToLLNuNu.root", otherb, l3[1], l2[1], l1[1], l0[1], scale[1]);
  ttz+=selection("events3/TTZToQQ.root", otherb, l3[1], l2[1], l1[1], l0[1], scale[1]);
  float tt =selection("events3/TTJets.root", otherb, l3[2], l2[2], l1[2], l0[2], scale[2]);
  float www=selection("events3/WWW.root", otherb, l3[3], l2[3], l1[3], l0[3], scale[3]);
  float wzz=selection("events3/WZZ.root", otherb, l3[4], l2[4], l1[4], l0[4], scale[4]);
  float wwz=selection("events3/WWZ.root", otherb, l3[5], l2[5], l1[5], l0[5], scale[5]);
  float ww =selection("events3/WW.root", ttbar, l3[6], l2[6], l1[6], l0[6], scale[6]);
  float zz =selection("events3/ZZ.root", otherb, l3[7], l2[7], l1[7], l0[7], scale[7]);
  float wz =selection("events3/WZ.root", otherb, l3[8], l2[8], l1[8], l0[8], scale[8]);
  float dy =selection("events3/DYJetsToLL_HT_1200_2500.root", otherb, l3[9], l2[9], l1[9], l0[9], scale[9]);
  dy+=selection("events3/DYJetsToLL_HT_200_400.root", otherb, l3[9], l2[9], l1[9], l0[9], scale[9]);
  dy+=selection("events3/DYJetsToLL_HT_2500_Inf.root", otherb, l3[9], l2[9], l1[9], l0[9], scale[9]);
  float ttw=selection("events3/TTWJetsToLNu.root", otherb, l3[10], l2[10], l1[10], l0[10], scale[10]);
  ttw+=selection("events3/TTWJetsToQQ.root", otherb, l3[10], l2[10], l1[10], l0[10], scale[10]);
  */
  
  TCanvas *c = MakeCanvas("c","c",800, 600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleX(0.1f);
  gStyle->SetTitleW(0.8f);
  gStyle->SetTitleSize(0.03,"t");
  c->SetRightMargin(0.09);
  
  //Drawing plots (look at drawFxns.hh)
  DrawHists(c, "constructed Z mass from electron.png", signal.zMe, ttbar.zMe, otherb.zMe);
  DrawHists(c, "constructed Z mass from muon.png", signal.zMm, ttbar.zMm, otherb.zMm);
  DrawHists(c, "constructed Z mass from jet.png", signal.zMj, ttbar.zMj, otherb.zMj);

  DrawHists(c, "number of lepton pairs for constructed Z.png", signal.nl, ttbar.nl, otherb.nl);
  DrawHists(c, "MET.png", signal.met, ttbar.met, otherb.met);
  DrawHists(c, "max PT of leptons.png", signal.maxPT, ttbar.maxPT, otherb.maxPT);
  DrawHists(c, "constructed Z PT from electron.png", signal.zPTe, ttbar.zPTe, otherb.zPTe);
  DrawHists(c, "constructed Z PT from muon.png", signal.zPTm, ttbar.zPTm, otherb.zPTm);
  DrawHists(c, "constructed Z PT from jet.png", signal.zPTj, ttbar.zPTj, otherb.zPTj);
   
  DrawHists(c, "number of electrons.png", signal.nume, ttbar.nume, otherb.nume);
  DrawHists(c, "number of muons.png", signal.numm, ttbar.numm, otherb.numm);
  DrawHists(c, "number of jets.png", signal.numj, ttbar.numj, otherb.numj);

  DrawHists(c, "mass of ZZZ.png", signal.mZZZ, ttbar.mZZZ, otherb.mZZZ);
  DrawHists(c, "PT of ZZZ.png", signal.ptZZZ, ttbar.ptZZZ, otherb.ptZZZ);

  DrawHists(c, "deltaR between electrons for the same Z.png", signal.dre, ttbar.dre, otherb.dre);
  DrawHists(c, "deltaR between muons for the same Z.png", signal.drm, ttbar.drm, otherb.drm);
  DrawHists(c, "min deltaR among the ZZZ.png", signal.drZZZ, ttbar.drZZZ, otherb.drZZZ);
  /*DrawHists(c, "min deltaR among the ZZZ.png", signal.min, ttbar.min, otherb.min);
  DrawHists(c, "max deltaR among the ZZZ.png", signal.max, ttbar.max, otherb.max);
  DrawHists(c, "deltaR among the ZZZ.png", signal.all, ttbar.all, otherb.all);*/
  
  /*
  ofstream myfile;
  myfile.open ("card_cutcount.txt");
  myfile << "imax 7 number of channels" <<"\n";
  myfile << "jmax 10 number of backgrounds" <<"\n";
  myfile << "kmax * number of nuisance parameters (sources of systematical uncertainties)" <<"\n";
  myfile << "--------------------------------" <<"\n";
  myfile << "bin lll eej emj mmj ejj mjj jjj" <<"\n";
  myfile << "observation 0 0 0 0 0 0 0" <<"\n";
  myfile << "--------------------------------" <<"\n";
  myfile << "bin lll lll lll lll lll lll lll lll lll lll lll eej eej eej eej eej eej eej eej eej eej eej emj emj emj emj emj emj emj emj emj emj emj mmj mmj mmj mmj mmj mmj mmj mmj mmj mmj mmj ejj ejj ejj ejj ejj ejj ejj ejj ejj ejj ejj mjj mjj mjj mjj mjj mjj mjj mjj mjj mjj mjj jjj jjj jjj jjj jjj jjj jjj jjj jjj jjj jjj" <<"\n";
  myfile << "process zzz ttz tt www wzz wwz ww zz2l hzz zhzz zz4l zzz ttz tt www wzz wwz ww zz2l hzz zhzz zz4l zzz ttz tt www wzz wwz ww zz2l hzz zhzz zz4l zzz ttz tt www wzz wwz ww zz2l hzz zhzz zz4l zzz ttz tt www wzz wwz ww zz2l hzz zhzz zz4l zzz ttz tt www wzz wwz ww zz2l hzz zhzz zz4l zzz ttz tt www wzz wwz ww zz2l hzz zhzz zz4l" <<"\n";
  myfile << "process 0 1 2 3 4 5 6 7 8 9 10 0 1 2 3 4 5 6 7 8 9 10 0 1 2 3 4 5 6 7 8 9 10 0 1 2 3 4 5 6 7 8 9 10 0 1 2 3 4 5 6 7 8 9 10 0 1 2 3 4 5 6 7 8 9 10 0 1 2 3 4 5 6 7 8 9 10" <<"\n";

  //myfile << "rate " << l3[0] << " " << l3[1] << " " << l3[2] << " " << l3[3] << " " << l3[4] << " " << l3[5] << " " <<	l3[6] << " " << l3[7] << " " <<	l2[0] << " " << l2[1] << " " << l2[2] << " " << l2[3] << " " <<	l2[4] << " " << l2[5] << " " <<	l2[6] << " " << l2[7] << " " << l1[0] << " " << l1[1] << " " << l1[2] << " " << l1[3] << " " <<	l1[4] << " " << l1[5] << " " <<	l1[6] << " " << l1[7] << " " << l0[0] << " " << l0[1] << " " << l0[2] << " " << l0[3] << " " <<	l0[4] << " " << l0[5] << " " <<	l0[6] << " " << l0[7] << "\n";
  myfile << "rate " << scale[0]*lll[0] << " " << scale[1]*lll[1] << " " << scale[2]*lll[2] << " " << scale[3]*lll[3] << " " <<  scale[4]*lll[4] << " " << scale[5]*lll[5] << " " << scale[6]*lll[6] << " " << scale[7]*lll[7] << " " << scale[8]*lll[8] << " " << scale[9]*lll[9] << " " << scale[10]*lll[10] << " " << scale[0]*eej[0] << " " << scale[1]*eej[1] << " " << scale[2]*eej[2] << " " << scale[3]*eej[3] << " " << scale[4]*eej[4] << " " << scale[5]*eej[5] << " " << scale[6]*eej[6] << " " << scale[7]*eej[7] << " " << scale[8]*eej[8] << " " << scale[9]*eej[9] << " " << scale[10]*eej[10] << " " << scale[0]*emj[0] << " " << scale[1]*emj[1] << " " << scale[2]*emj[2] << " " << scale[3]*emj[3] << " " << scale[4]*emj[4] << " " << scale[5]*emj[5] << " " << scale[6]*emj[6] << " " << scale[7]*emj[7] << " " << scale[8]*emj[8] << " " << scale[9]*emj[9] << " " << scale[10]*emj[10] << " " << scale[0]*mmj[0] << " " << scale[1]*mmj[1] << " " << scale[2]*mmj[2] << " " << scale[3]*mmj[3] << " " << scale[4]*mmj[4] << " " << scale[5]*mmj[5] << " " << scale[6]*mmj[6] << " " << scale[7]*mmj[7] << " " << scale[8]*mmj[8] << " " << scale[9]*mmj[9] << " " << scale[10]*mmj[10] << " " << scale[0]*ejj[0] << " " << scale[1]*ejj[1] << " " << scale[2]*ejj[2] << " " << scale[3]*ejj[3] << " " << scale[4]*ejj[4] << " " << scale[5]*ejj[5] << " " << scale[6]*ejj[6] << " " << scale[7]*ejj[7] << " " << scale[8]*ejj[8] << " " << scale[9]*ejj[9] << " " << scale[10]*ejj[10] << " " << scale[0]*mjj[0] << " " << scale[1]*mjj[1] << " " << scale[2]*mjj[2] << " " << scale[3]*mjj[3] << " " << scale[4]*mjj[4] << " " << scale[5]*mjj[5] << " " << scale[6]*mjj[6] << " " << scale[7]*mjj[7] << " " << scale[8]*mjj[8] << " " << scale[9]*mjj[9] << " " << scale[10]*mjj[10] << " " << scale[0]*jjj[0] << " " << scale[1]*jjj[1] << " " << scale[2]*jjj[2] << " " << scale[3]*jjj[3] << " " << scale[4]*jjj[4] << " " << scale[5]*jjj[5] << " " << scale[6]*jjj[6] << " " << scale[7]*jjj[7] << " " << scale[8]*jjj[8] << " " << scale[9]*jjj[9] << " " << scale[10]*jjj[10] << "\n";
  
  myfile << "--------------------------------" <<"\n";
  */
  /*
  myfile << "Ezzzl3 gmN " << l3[0] << " " << scale[0] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezzzl2 gmN " << l2[0] << " - - - - - - - - " << scale[0] << " - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezzzl1 gmN " << l1[0] << " - - - - - - - - - - - - - - - - " << scale[0] << " - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezzzl0 gmN " << l0[0] << " - - - - - - - - - - - - - - - - - - - - - - - - " << scale[0] << " - - - - - - -" <<"\n";
  myfile << "Ettzl3 gmN " << l3[1] << " - " << scale[1] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ettzl2 gmN " << l2[1] << " - - - - - - - - - " << scale[1] << " - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ettzl1 gmN " << l1[1] << " - - - - - - - - - - - - - - - - - " << scale[1] << " - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ettzl0 gmN " << l0[1] << " - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[1] << " - - - - - -" <<"\n";
  myfile << "Ettl3 gmN " << l3[2] << " - - " << scale[2] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ettl2 gmN " << l2[2] << " - - - - - - - - - - " << scale[2] << " - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ettl1 gmN " << l1[2] << " - - - - - - - - - - - - - - - - - - " << scale[2] << " - - - - - - - - - - - - -" <<"\n";
  myfile << "Ettl0 gmN " << l0[2] << " - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[2] << " - - - - -" <<"\n";
  myfile << "Ewwwl3 gmN " << l3[3] << " - - - " << scale[3] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewwwl2 gmN " << l2[3] << " - - - - - - - - - - - " << scale[3] << " - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewwwl1 gmN " << l1[3] << " - - - - - - - - - - - - - - - - - - - " << scale[3] << " - - - - - - - - - - - -" <<"\n";
  myfile << "Ewwwl0 gmN " << l0[3] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[3] << " - - - -" <<"\n";
  myfile << "Ewzzl3 gmN " << l3[4] << " - - - - " << scale[4] << " - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewzzl2 gmN " << l2[4] << " - - - - - - - - - - - - " << scale[4] << " - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewzzl1 gmN " << l1[4] << " - - - - - - - - - - - - - - - - - - - - " << scale[4] << " - - - - - - - - - - -" <<"\n";
  myfile << "Ewzzl0 gmN " << l0[4] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[4] << " - - -" <<"\n";
  myfile << "Ewwzl3 gmN " << l3[5] << " - - - - - " << scale[5] << " - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewwzl2 gmN " << l2[5] << " - - - - - - - - - - - - - " << scale[5] << " - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewwzl1 gmN " << l1[5] << " - - - - - - - - - - - - - - - - - - - - - " << scale[5] << " - - - - - - - - - -" <<"\n";
  myfile << "Ewwzl0 gmN " << l0[5] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[5] << " - -" <<"\n";
  myfile << "Ewwl3 gmN " << l3[6] << " - - - - - - " << scale[6] << " - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewwl2 gmN " << l2[6] << " - - - - - - - - - - - - - - " << scale[6] << " - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewwl1 gmN " << l1[6] << " - - - - - - - - - - - - - - - - - - - - - - " << scale[6] << " - - - - - - - - -" <<"\n";
  myfile << "Ewwl0 gmN " << l0[6] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[6] << " -" <<"\n";
  myfile << "Ezzl3 gmN " << l3[7] << " - - - - - - - " << scale[7] << " - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezzl2 gmN " << l2[7] << " - - - - - - - - - - - - - - - " << scale[7] << " - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezzl1 gmN " << l1[7] << " - - - - - - - - - - - - - - - - - - - - - - - " << scale[7] << " - - - - - - - -" <<"\n";
  myfile << "Ezzl0 gmN " << l0[7] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[7] <<"\n";
  */

  /*
  myfile << "Ezzz0 gmN " << lll[0] << " " << scale[0] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezzz1 gmN " << eej[0] << " - - - - - - - - - - - " << scale[0] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezzz2 gmN " << emj[0] << " - - - - - - - - - - - - - - - - - - - - - - " << scale[0] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezzz3 gmN " << mmj[0] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[0] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezzz4 gmN " << ejj[0] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[0] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezzz5 gmN " << mjj[0] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[0] << " - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezzz6 gmN " << jjj[0] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[0] << " - - - - - - - - - -" <<"\n";  
  myfile << "Ettz0 gmN " << lll[1] << " - " << scale[1] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ettz1 gmN " << eej[1] << " - - - - - - - - - - - - " << scale[1] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ettz2 gmN " << emj[1] << " - - - - - - - - - - - - - - - - - - - - - - - " << scale[1] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ettz3 gmN " << mmj[1] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[1] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ettz4 gmN " << ejj[1] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[1] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ettz5 gmN " << mjj[1] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[1] << " - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ettz6 gmN " << jjj[1] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[1] << " - - - - - - - - -" <<"\n";							   
  myfile << "Ett0 gmN " << lll[2] << " - - " << scale[2] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ett1 gmN " << eej[2] << " - - - - - - - - - - - - - " << scale[2] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ett2 gmN " << emj[2] << " - - - - - - - - - - - - - - - - - - - - - - - - " << scale[2] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ett3 gmN " << mmj[2] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[2] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ett4 gmN " << ejj[2] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[2] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ett5 gmN " << mjj[2] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[2] << " - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ett6 gmN " << jjj[2] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[2] << " - - - - - - - -" <<"\n";																		   
  myfile << "Ewww0 gmN " << lll[3] << " - - - " << scale[3] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewww1 gmN " << eej[3] << " - - - - - - - - - - - - - - " << scale[3] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewww2 gmN " << emj[3] << " - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[3] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewww3 gmN " << mmj[3] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[3] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewww4 gmN " << ejj[3] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[3] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewww5 gmN " << mjj[3] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[3] << " - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewww6 gmN " << jjj[3] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[3] << " - - - - - - -" <<"\n";																					   
  myfile << "Ewzz0 gmN " << lll[4] << " - - - - " << scale[4] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewzz1 gmN " << eej[4] << " - - - - - - - - - - - - - - - " << scale[4] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewzz2 gmN " << emj[4] << " - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[4] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewzz3 gmN " << mmj[4] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[4] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewzz4 gmN " << ejj[4] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[4] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewzz5 gmN " << mjj[4] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[4] << " - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewzz6 gmN " << jjj[4] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[4] << " - - - - - -" <<"\n";																							   
  myfile << "Ewwz0 gmN " << lll[5] << " - - - - - " << scale[5] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewwz1 gmN " << eej[5] << " - - - - - - - - - - - - - - - - " << scale[5] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewwz2 gmN " << emj[5] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[5] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewwz3 gmN " << mmj[5] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[5] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewwz4 gmN " << ejj[5] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[5] << " - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewwz5 gmN " << mjj[5] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[5] << " - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ewwz6 gmN " << jjj[5] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[5] << " - - - - -" <<"\n";																						   
  myfile << "Eww0 gmN " << lll[6] << " - - - - - - " << scale[6] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Eww1 gmN " << eej[6] << " - - - - - - - - - - - - - - - - - " << scale[6] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Eww2 gmN " << emj[6] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[6] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Eww3 gmN " << mmj[6] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[6] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Eww4 gmN " << ejj[6] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[6] << " - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Eww5 gmN " << mjj[6] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[6] << " - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Eww6 gmN " << jjj[6] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[6] << " - - - -" <<"\n";																							   
  myfile << "Ezz2l0 gmN " << lll[7] << " - - - - - - - " << scale[7] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezz2l1 gmN " << eej[7] << " - - - - - - - - - - - - - - - - - - " << scale[7] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezz2l2 gmN " << emj[7] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[7] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezz2l3 gmN " << mmj[7] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[7] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezz2l4 gmN " << ejj[7] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[7] << " - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezz2l5 gmN " << mjj[7] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[7] << " - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezz2l6 gmN " << jjj[7] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[7] << " - - -" <<"\n";																						   
  myfile << "Ehzz0 gmN " << lll[8] << " - - - - - - - - " << scale[8] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ehzz1 gmN " << eej[8] << " - - - - - - - - - - - - - - - - - - - " << scale[8] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ehzz2 gmN " << emj[8] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[8] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ehzz3 gmN " << mmj[8] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[8] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ehzz4 gmN " << ejj[8] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[8] << " - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ehzz5 gmN " << mjj[8] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[8] << " - - - - - - - - - - - - -" <<"\n";
  myfile << "Ehzz6 gmN " << jjj[8] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[8] << " - -" <<"\n";																							   
  myfile << "Ezhzz0 gmN " << lll[9] << " - - - - - - - - - " << scale[9] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezhzz1 gmN " << eej[9] << " - - - - - - - - - - - - - - - - - - - - " << scale[9] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezhzz2 gmN " << emj[9] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[9] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezhzz3 gmN " << mmj[9] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[9] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezhzz4 gmN " << ejj[9] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[9] << " - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezhzz5 gmN " << mjj[9] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[9] << " - - - - - - - - - - - -" <<"\n";
  myfile << "Ezhzz6 gmN " << jjj[9] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[9] << " -" <<"\n";																						   
  myfile << "Ezz4l0 gmN " << lll[10] << " - - - - - - - - - - " << scale[10] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezz4l1 gmN " << eej[10] << " - - - - - - - - - - - - - - - - - - - - - " << scale[10] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezz4l2 gmN " << emj[10] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[10] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezz4l3 gmN " << mmj[10] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[10] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezz4l4 gmN " << ejj[10] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[10] << " - - - - - - - - - - - - - - - - - - - - - -" <<"\n";
  myfile << "Ezz4l5 gmN " << mjj[10] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[10] << " - - - - - - - - - - -" <<"\n";
  myfile << "Ezz4l6 gmN " << jjj[10] << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << scale[10] <<"\n";
  myfile << "lumi lnN 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01" << endl;
  myfile << "e lnN " << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " - - - - - - - - - - - " << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " - - - - - - - - - - - - - - - - - - - - - -" << "\n";
  myfile << "m lnN - - - - - - - - - - - - - - - - - - - - - - " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " - - - - - - - - - - - " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " - - - - - - - - - - - " << "\n";
  myfile << "j lnN - - - - - - - - - - - " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(2*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " "  << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(4*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << " " << 1+sqrt(6*(0.01)*(0.01)) << "\n";
																							   
  myfile.close();
  */
}

Float_t selection(const TString inputfile, HistHolder hists, int &l3, int &l2, int &l1, int &l0, float &evtWeight) {
//Float_t selection(const TString inputfile, HistHolder hists, int &lll, int &eej, int &emj, int &mmj, int &ejj, int &mjj, int &jjj, float &evtWeight, TDirectory *cdir[7]) {
  //Float_t selection(const TString inputfile, HistHolder hists, int &lll, int &eej, int &emj, int &mmj, int &ejj, int &mjj, int &jjj, float &evtWeight) {
  
  Float_t eventWeight;
  uint nEl, nMu, nJet, nBJet;
  Float_t metPt;
  Float_t metPhi;

  TClonesArray *vEl  = new TClonesArray("TLepton");
  TClonesArray *vMu  = new TClonesArray("TLepton");
  TClonesArray *vJet = new TClonesArray("THadJet");
  
  TFile *inFile = new TFile(inputfile, "READ");
  TTree *inTree = (TTree*) inFile->Get("Events");
  inTree->SetBranchAddress("eventWeightN", &eventWeight);
  inTree->SetBranchAddress("nEl", &nEl);
  inTree->SetBranchAddress("nMu", &nMu);
  inTree->SetBranchAddress("nJet", &nJet);
  inTree->SetBranchAddress("nBJet", &nBJet);
  inTree->SetBranchAddress("metPt", &metPt);
  inTree->SetBranchAddress("metPhi", &metPhi);

  inTree->SetBranchAddress("vEl", &vEl);
  inTree->SetBranchAddress("vMu", &vMu);
  inTree->SetBranchAddress("vJet", &vJet);

  Long64_t allEntries = inTree->GetEntries();
  TH1D *hTot = (TH1D*) inFile->Get("hTot");
  float nTot = hTot->GetBinContent(1);

  cout << "** File contains " << allEntries << " events" << endl;
  
  float nCount = 0;

  // Loop over all events
  for(Long64_t entry = 0; entry < allEntries; ++entry)
  {
    inTree->GetEntry(entry);
    Int_t ide[6] = {-1, -1, -1, -1, -1, -1};
    Int_t idm[6] = {-1, -1, -1, -1, -1, -1};
    Float_t ideeta[6] = {-1, -1, -1, -1, -1, -1};
    Float_t idmeta[6] = {-1, -1, -1, -1, -1, -1};
    Float_t idephi[6] = {-1, -1, -1, -1, -1, -1};
    Float_t idmphi[6] = {-1, -1, -1, -1, -1, -1};

    Int_t idj[6] = {-1, -1, -1, -1, -1, -1};
    Double_t min[9] = {15, 15, 15, 15, 15, 15, 15, 15, 15};
    Double_t pt[9] = {-15, -15, -15, -15, -15, -15, -15, -15, -15};
    Double_t eta[9] = {15, 15, 15, 15, 15, 15, 15, 15, 15};
    Double_t phi[9] = {15, 15, 15, 15, 15, 15, 15, 15, 15};
    Double_t m[9] = {15, 15, 15, 15, 15, 15, 15, 15, 15};
  
    if (nBJet>0) continue;
    //if (nEl==1) continue;
    //if (nMu==1) continue;

    //evtWeight = (eventWeight * LUMI) / nTot;
    evtWeight = eventWeight * LUMI;
    
    float maxpt = 0;

    Int_t num = nEl / 2; // number of electron pairs we can make for each event
   
    for(Int_t k = 0; k < num and k < 3; ++k)
      {
	// Loop over all electrons in event
	for(Int_t i = 0; i < nEl; ++i)
	  {
	    for(Int_t j = i+1; j < nEl; ++j)
	      {
		if ((k == 1 or k == 2) and i == ide[0] and j == ide[1])
		  continue;
		if (k == 2 and i == ide[2] and j == ide[3])
                  continue;

		TLepton *e1 = (TLepton*)vEl->At(i);                                                
		TLepton *e2 = (TLepton*)vEl->At(j); 

		if (e1->PT > maxpt) maxpt = e1->PT;
		TLorentzVector v1;                                                                                                                                                  
		v1.SetPtEtaPhiM(e1->PT, e1->Eta, e1->Phi, ELE_MASS);                                                                                                                         
		TLorentzVector v2;                                                                                                                                                           
		v2.SetPtEtaPhiM(e2->PT, e2->Eta, e2->Phi, ELE_MASS);                                                                                                                                       
		TLorentzVector v = v1 + v2;   
		//if (v.Pt() < 40) continue;
		float temp = deltaR(e1->Eta, e2->Eta, e1->Phi, e2->Phi);
		//if (temp > 2) continue;
		if (TMath::Abs(v.M()-Z_MASS) < min[k])
                  {
                    min[k] = TMath::Abs(v.M()-Z_MASS);
                    ide[2*k] = i;
                    ide[2*k+1] = j;
		    ideeta[2*k] = e1->Eta;
                    ideeta[2*k+1] = e2->Eta;
		    idephi[2*k] = e1->Phi;
                    idephi[2*k+1] = e2->Phi;
		    pt[k] = v.Pt();
		    eta[k] = v.Eta();
		    phi[k] = v.Phi();
		    m[k] = v.M();
		  }
	      }
	  }
      }
    	
    num = nMu / 2; // number of muon pairs we can make for each event                                                                                                    

    for(Int_t k = 0; k < num and k < 3; ++k)
      {
      
        for(Int_t i = 0; i < nMu; ++i)
          {
            for(Int_t j = i+1; j < nMu; ++j)
              {
                if ((k == 1 or k == 2) and i == idm[0] and j == idm[1])
                  continue;
                if (k == 2 and i == idm[2] and j == idm[3])
                  continue;

		TLepton *e1 = (TLepton*)vMu->At(i);
                TLepton *e2 = (TLepton*)vMu->At(j);

		if (e1->PT > maxpt) maxpt = e1->PT;
		TLorentzVector v1;
		v1.SetPtEtaPhiM(e1->PT, e1->Eta, e1->Phi, MUON_MASS);
		TLorentzVector v2;
		v2.SetPtEtaPhiM(e2->PT, e2->Eta, e2->Phi, MUON_MASS);
		TLorentzVector v = v1 + v2;
		//if (v.Pt() < 40) continue;
		float temp = deltaR(e1->Eta, e2->Eta, e1->Phi, e2->Phi);
		//if (temp > 2) continue;
		if (TMath::Abs(v.M()-Z_MASS) < min[k+3])
                  {
                    min[k+3] = TMath::Abs(v.M()-Z_MASS);
                    idm[2*k] = i;
                    idm[2*k+1] = j;
		    idmeta[2*k] = e1->Eta;
                    idmeta[2*k+1] = e2->Eta;
                    idmphi[2*k] = e1->Phi;
                    idmphi[2*k+1] = e2->Phi;
		    pt[k+3] = v.Pt();
                    eta[k+3] = v.Eta();
		    phi[k+3] = v.Phi();
                    m[k+3] = v.M();
		  }

              }
          }
      }

    //if (maxpt < 50) continue;
    
    num = nJet / 2; // number of jet pairs we can make for each event                                                                                                                  

    for(Int_t k = 0; k < num and k < 3; ++k)
      {

        for(Int_t i = 0; i < nJet; ++i)
          {
            for(Int_t j = i+1; j < nJet; ++j)
              {
		if ((k == 1 or k == 2) and i == idj[0] and j == idj[1])
                  continue;
                if (k == 2 and i == idj[2] and j == idj[3])
                  continue;
		
		THadJet *e1 = (THadJet*)vJet->At(i);
                THadJet *e2 = (THadJet*)vJet->At(j);
		
		TLorentzVector v1;
		v1.SetPtEtaPhiM(e1->PT, e1->Eta, e1->Phi, e1->Mass);
                TLorentzVector v2;
		v2.SetPtEtaPhiM(e2->PT, e2->Eta, e2->Phi, e2->Mass);
		TLorentzVector v = v1 + v2;
		
                if (TMath::Abs(v.M()-Z_MASS) < min[k+6])
                  {
                    min[k+6] = TMath::Abs(v.M()-Z_MASS);
                    idj[2*k] = i;
                    idj[2*k+1] = j;
		    pt[k+6] = v.Pt();
                    eta[k+6] = v.Eta();
		    phi[k+6] = v.Phi();
                    m[k+6] = v.M();
		    }
              }
          }
      }

    Int_t ctre = 0;
    Int_t ctrm = 0;
    Int_t ctrh = 0;
    Int_t index[3] = {-1, -1, -1};
    for(Int_t i = 0; i < 3; ++i)
      {
	Double_t temp = 15;
	//Double_t temp;
	//if (i < 2)
	  //temp = 15;
	//else
	//temp = 100;
	
	for(Int_t k = 0; k < 9; ++k)
	  {
	    if (pt[k] >= 0 and min[k] < temp)
	      {
		temp = min[k];
		index[i] = k;
	      }
	  }
	
	if (index[i] < 3 and index[i] >= 0)
	  ctre++;
	else if (index[i] < 6 and index[i] >= 3)
	  ctrm++;
	else if (index[i] >= 6)
	  ctrh++;
	if (index[i] >= 0) min[index[i]] = 100;	
      }

    if (ctre+ctrm+ctrh == 3)
      {

	TLorentzVector v1;
	v1.SetPtEtaPhiM(pt[index[0]], eta[index[0]], phi[index[0]], m[index[0]]);
        TLorentzVector v2;
        v2.SetPtEtaPhiM(pt[index[1]], eta[index[1]], phi[index[1]], m[index[1]]);
        TLorentzVector v3;
        v3.SetPtEtaPhiM(pt[index[2]], eta[index[2]], phi[index[2]], m[index[2]]);
        TLorentzVector v = v1 + v2 + v3;
	/*
	if (ctre+ctrm == 3) lll++;
	else if (ctre == 2 and ctrh == 1) eej++;
	else if (ctre == 1 and ctrm == 1 and ctrh == 1) emj++;
	else if (ctrm == 2 and ctrh == 1) mmj++;
	else if (ctre == 1 and ctrh == 2) ejj++;
	else if (ctrm == 1 and ctrh == 2) mjj++;
	else jjj++;
	*/
	/*
	if (ctre+ctrm == 3)
	  {
	    lll++;
	    cdir[0]->cd();
	    if (inputfile == "events2/ZZZ.root")
	      {
		TH1D *ZZZ = (TH1D*) cdir[0]->Get("ZZZ");
		ZZZ->Fill(v.M(), evtWeight);
	      }
	    else
	      {
		TH1D *bkg = (TH1D*) cdir[0]->Get("bkg");
		bkg->Fill(v.M(), evtWeight);
	      }
	  }
	else if (ctre == 2 and ctrh == 1)
	  {
	    eej++;
	    cdir[1]->cd();
	    if (inputfile == "events2/ZZZ.root")
	      {
		TH1D *ZZZ = (TH1D*) cdir[1]->Get("ZZZ");
		ZZZ->Fill(v.M(), evtWeight);
	      }
	    else
	      {
		TH1D *bkg = (TH1D*) cdir[1]->Get("bkg");
		bkg->Fill(v.M(), evtWeight);
	      }
	  }
	else if (ctre == 1 and ctrm == 1 and ctrh == 1)
	  {
	    emj++;
	    cdir[2]->cd();
	    if (inputfile == "events2/ZZZ.root")
	      {
		TH1D *ZZZ = (TH1D*) cdir[2]->Get("ZZZ");
		ZZZ->Fill(v.M(), evtWeight);
	      }
	    else
	      {
		TH1D *bkg = (TH1D*) cdir[2]->Get("bkg");
		bkg->Fill(v.M(), evtWeight);
	      }
	  }
	else if (ctrm == 2 and ctrh == 1)
	  {
	    mmj++;
	    cdir[3]->cd();
	    if (inputfile == "events2/ZZZ.root")
	      {
		TH1D *ZZZ = (TH1D*) cdir[3]->Get("ZZZ");
		ZZZ->Fill(v.M(), evtWeight);
	      }
	    else
	      {
		TH1D *bkg = (TH1D*) cdir[3]->Get("bkg");
		bkg->Fill(v.M(), evtWeight);
	      }
	  }
	else if (ctre == 1 and ctrh == 2)
	  {
	    ejj++;
	     cdir[4]->cd();
	    if (inputfile == "events2/ZZZ.root")
	      {
		TH1D *ZZZ = (TH1D*) cdir[4]->Get("ZZZ");
		ZZZ->Fill(v.M(), evtWeight);
	      }
	    else
	      {
		TH1D *bkg = (TH1D*) cdir[4]->Get("bkg");
		bkg->Fill(v.M(), evtWeight);
	      }
	  }
	else if (ctrm == 1 and ctrh == 2)
	  {
	    mjj++;
	    cdir[5]->cd();
	    if (inputfile == "events2/ZZZ.root")
	      {
		TH1D *ZZZ = (TH1D*) cdir[5]->Get("ZZZ");
		ZZZ->Fill(v.M(), evtWeight);
	      }
	    else
	      {
		TH1D *bkg = (TH1D*) cdir[5]->Get("bkg");
		bkg->Fill(v.M(), evtWeight);
	      }
	  }
	else
	  {
	    jjj++;
	    cdir[6]->cd();
	    if (inputfile == "events2/ZZZ.root")
	      {
		TH1D *ZZZ = (TH1D*) cdir[6]->Get("ZZZ");
		ZZZ->Fill(v.M(), evtWeight);
	      }
	    else
	      {
		TH1D *bkg = (TH1D*) cdir[6]->Get("bkg");
		bkg->Fill(v.M(), evtWeight);
	      }
	  }
	*/
	
	float temp1 = deltaR(v1.Eta(), v2.Eta(), v1.Phi(), v2.Phi());
	float temp2 = deltaR(v1.Eta(), v3.Eta(), v1.Phi(), v3.Phi());
	float temp3 = deltaR(v3.Eta(), v2.Eta(), v3.Phi(), v2.Phi());
	float mindR;
	float maxdR;
	if (temp1<=temp2 and temp1<=temp3) mindR = temp1;
	else if (temp2<=temp1 and temp2<=temp3) mindR = temp2;
	else mindR = temp3;
	if (temp1>=temp2 and temp1>=temp3) maxdR = temp1;
	else if (temp2>=temp1 and temp2>=temp3) maxdR = temp2;
	else maxdR = temp3;
	
	//if (v.M() > 350 and mindR > 1.5) {
	  
	nCount+=evtWeight;
	hists.nume->Fill(nEl, evtWeight);
	hists.numm->Fill(nMu, evtWeight);
	hists.numj->Fill(nJet, evtWeight);
	hists.maxPT->Fill(maxpt, evtWeight);
	if (ctrh == 0)
	  {
	    hists.nl->Fill(3.0, evtWeight);
	    l3++;
	  }
	else if (ctrh == 1)
	  {
	    hists.nl->Fill(2.0, evtWeight);
	    l2++;
	    /*hists.min->Fill(mindR, evtWeight);
	    hists.max->Fill(maxdR, evtWeight);
	    hists.all->Fill(temp1, evtWeight);
	    hists.all->Fill(temp2, evtWeight);
	    hists.all->Fill(temp3, evtWeight);*/
	    
	  }
	else if (ctrh == 2)
	  {
	    hists.nl->Fill(1.0, evtWeight);
	    l1++;
	  }
	else
	  {
	    hists.nl->Fill(0.0, evtWeight);
	    l0++;
	  }
	
	for(Int_t i = 0; i < 3; ++i)
	  {
	    if (index[i]<3)
	      {
		hists.zMe->Fill(m[index[i]], evtWeight);
		hists.zPTe->Fill(pt[index[i]], evtWeight);
		float temp = deltaR(ideeta[2*index[i]], ideeta[2*index[i]+1], idephi[2*index[i]], idephi[2*index[i]+1]);
		hists.dre->Fill(temp, evtWeight);
	      }
	    else if (index[i]>=3 and index[i]<6)
	      {
		hists.zMm->Fill(m[index[i]], evtWeight);
		hists.zPTm->Fill(pt[index[i]], evtWeight);
		float temp = deltaR(idmeta[2*(index[i]-3)], idmeta[2*(index[i]-3)+1], idmphi[2*(index[i]-3)], idmphi[2*(index[i]-3)+1]);
                hists.drm->Fill(temp, evtWeight);
	      }
	    else
	      {
		hists.zMj->Fill(m[index[i]], evtWeight);
		hists.zPTj->Fill(pt[index[i]], evtWeight);
	      }
	  }
	
	hists.mZZZ->Fill(v.M(), evtWeight);
	hists.ptZZZ->Fill(v.Pt(), evtWeight);
	hists.drZZZ->Fill(mindR, evtWeight);
	}
    //}
    hists.met->Fill(metPt, evtWeight);

  }

  return nCount;
}
