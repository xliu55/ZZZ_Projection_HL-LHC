#ifndef DRAW_FXNS_HH
#define DRAW_FXNS_HH

struct HistHolder {
  TH1D *zMe; 
  TH1D *zMm;
  TH1D *zMj;

  TH1D *zPTe; 
  TH1D *zPTm;
  TH1D *zPTj;
  
  TH1D *nl;
  TH1D *met;
  TH1D *maxPT;

  TH1D *nume;
  TH1D *numm;
  TH1D *numj;

  TH1D *ptZZZ; 
  TH1D *mZZZ; 

  TH1D *dre;
  TH1D *drm;
  TH1D *min;
  TH1D *max;
  TH1D *all;
  TH1D *drZZZ;
  
};

void initHistHolder(HistHolder &h, TString id) {

  const int nBins_zMass=20, xMin_zMass=75, xMax_zMass=107;
  const int nBins_lPT  =25, xMin_lPT  = 0, xMax_lPT  =200;
  const int nBins_nObj =4, xMin_nObj = 0, xMax_nObj = 4;
  const int nBins_3zpt=25, xMin_3zpt=0, xMax_3zpt=500;
  const int nBins_3zMass=25, xMin_3zMass=220, xMax_3zMass=1200;
  const int nBins_dr=20, xMin_dr=0, xMax_dr=5;
  const int nBins_met=20, xMin_met=0, xMax_met=300;
  
  h.zMe = new TH1D(Form("%s_zMe", id.Data()), Form("zMe"), 
		   nBins_zMass, xMin_zMass, xMax_zMass);
  h.zMe->GetXaxis()->SetTitle("m_{ee} [GeV]");
  h.zMe->GetYaxis()->SetTitle("Events");
  h.zMe->Sumw2();
  //h.zMe->SetTitle("");
  
  h.zMm = new TH1D(Form("%s_zMm", id.Data()), Form("zMm"),
                   nBins_zMass, xMin_zMass, xMax_zMass);
  h.zMm->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
  h.zMm->GetYaxis()->SetTitle("Events");
  h.zMm->Sumw2();

  h.zMj = new TH1D(Form("%s_zMj", id.Data()), Form("zMj"),
                   nBins_zMass, xMin_zMass, xMax_zMass);
  h.zMj->GetXaxis()->SetTitle("m_{jj} [GeV]");
  h.zMj->GetYaxis()->SetTitle("Events");
  h.zMj->Sumw2();

  h.zPTe = new TH1D(Form("%s_zPTe", id.Data()), Form("zPTe"),
		    nBins_lPT, xMin_lPT, 500);

  h.zPTe->GetXaxis()->SetTitle("PT of Z (Z to ee) [GeV]");
  h.zPTe->GetYaxis()->SetTitle("Events");
  h.zPTe->Sumw2();

  h.zPTm = new TH1D(Form("%s_zPTm", id.Data()), Form("zPTm"),
                    nBins_lPT, xMin_lPT, 500);

  h.zPTm->GetXaxis()->SetTitle("PT of Z (Z to #mu#mu) [GeV]");
  h.zPTm->GetYaxis()->SetTitle("Events");
  h.zPTm->Sumw2();

  h.zPTj = new TH1D(Form("%s_zPTj", id.Data()), Form("zPTj"),
                    nBins_lPT, xMin_lPT, xMax_lPT);

  h.zPTj->GetXaxis()->SetTitle("PT of Z (Z to jj) [GeV]");
  h.zPTj->GetYaxis()->SetTitle("Events");
  h.zPTj->Sumw2();

  h.maxPT = new TH1D(Form("%s_maxPT", id.Data()), Form("maxPT"),
                    25, 0, 300);

  h.maxPT->GetXaxis()->SetTitle("max PT of leptons [GeV]");
  h.maxPT->GetYaxis()->SetTitle("Events");
  h.maxPT->Sumw2();

  h.nl = new TH1D(Form("%s_nl", id.Data()), Form("nl"), 
		  nBins_nObj, xMin_nObj, xMax_nObj);

  h.nl->GetXaxis()->SetTitle("number of lepton pairs per event");
  h.nl->GetYaxis()->SetTitle("Events");
  h.nl->Sumw2();

  h.nume = new TH1D(Form("%s_nume", id.Data()), Form("nume"),
		    6, 0, 6);

  h.nume->GetXaxis()->SetTitle("number of electrons per event");
  h.nume->GetYaxis()->SetTitle("Events");
  h.nume->Sumw2();

  h.numm = new TH1D(Form("%s_numm", id.Data()), Form("numm"),
		    6, 0, 6);

  h.numm->GetXaxis()->SetTitle("number of muons per event");
  h.numm->GetYaxis()->SetTitle("Events");
  h.numm->Sumw2();

  h.numj = new TH1D(Form("%s_numj", id.Data()), Form("numj"),
		    19, 2, 20);

  h.numj->GetXaxis()->SetTitle("number of jets  per event");
  h.numj->GetYaxis()->SetTitle("Events");
  h.numj->Sumw2();


  h.met = new TH1D(Form("%s_met", id.Data()), Form("met"),
                  nBins_met, xMin_met, xMax_met);

  h.met->GetXaxis()->SetTitle("MET [GeV]");
  h.met->GetYaxis()->SetTitle("Events");
  h.met->Sumw2();

  h.mZZZ = new TH1D(Form("%s_mZZZ", id.Data()), Form("mZZZ"),
                   nBins_3zMass, xMin_3zMass, xMax_3zMass);
  h.mZZZ->GetXaxis()->SetTitle("m of ZZZ [GeV]");
  h.mZZZ->GetYaxis()->SetTitle("Events");
  h.mZZZ->Sumw2();

  h.ptZZZ = new TH1D(Form("%s_ptZZZ", id.Data()), Form("ptZZZ"),
                    nBins_3zpt, xMin_3zpt, xMax_3zpt);

  h.ptZZZ->GetXaxis()->SetTitle("PT of ZZZ [GeV]");
  h.ptZZZ->GetYaxis()->SetTitle("Events");
  h.ptZZZ->Sumw2();

  h.dre = new TH1D(Form("%s_dre", id.Data()), Form("dre"),
                  nBins_dr, xMin_dr, xMax_dr);

  h.dre->GetXaxis()->SetTitle("deltaR_{ee}");
  h.dre->GetYaxis()->SetTitle("Events");
  h.dre->Sumw2();

  h.drm = new TH1D(Form("%s_drm", id.Data()), Form("drm"),
                  nBins_dr, xMin_dr, xMax_dr);

  h.drm->GetXaxis()->SetTitle("deltaR_{#mu#mu}");
  h.drm->GetYaxis()->SetTitle("Events");
  h.drm->Sumw2();

  h.drZZZ = new TH1D(Form("%s_drZZZ", id.Data()), Form("drZZZ"),
                  nBins_dr, xMin_dr, xMax_dr);

  h.drZZZ->GetXaxis()->SetTitle("min_deltaR_{ZZZ}");
  h.drZZZ->GetYaxis()->SetTitle("Events");
  h.drZZZ->Sumw2();
  
  h.min = new TH1D(Form("%s_mindR", id.Data()), Form("mindR"),
                  nBins_dr, 1.5, xMax_dr);

  h.min->GetXaxis()->SetTitle("min_deltaR_{ZZZ}");
  h.min->GetYaxis()->SetTitle("Events");
  h.min->Sumw2();

  h.max = new TH1D(Form("%s_maxdR", id.Data()), Form("maxdR"),
		  nBins_dr, 1.5, xMax_dr);

  h.max->GetXaxis()->SetTitle("max_deltaR_{ZZZ}");
  h.max->GetYaxis()->SetTitle("Events");
  h.max->Sumw2();

  h.all = new TH1D(Form("%s_alldR", id.Data()), Form("alldR"),
                  nBins_dr, 1.5, xMax_dr);

  h.all->GetXaxis()->SetTitle("all_deltaR_{ZZZ}");
  h.all->GetYaxis()->SetTitle("Events");
  h.all->Sumw2();

}

void DrawHists(TCanvas *c, TString name, TH1D* signal, TH1D* ttbar, TH1D* otherb) {

  signal->SetLineColor(kBlue); 
  ttbar->SetLineColor(kRed); ttbar->SetFillColor(kRed); //ttbar->SetFillStyle(1001);
  otherb->SetLineColor(kGreen); otherb->SetFillColor(kGreen); //otherb->SetFillStyle(1001);

  otherb->Add(ttbar);
  //  signal->Add(otherb);

  TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetShadowColor(0);

  leg->AddEntry(signal, "ZZZ", "lf");
  leg->AddEntry(ttbar, "TT", "lf");
  leg->AddEntry(otherb, "Other", "lf");

  otherb->Draw("hist E");
  ttbar->Draw("hist E same");
  c->Update();

  Float_t rightmax = 1.1*signal->GetMaximum();
  Float_t scale = gPad->GetUymax()/rightmax;
  signal->Scale(scale);
  signal->SetTitle("");
  signal->Draw("hist E same");

  //draw an axis on the right side
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510, "+L");
  axis->SetLineColor(kBlue);
  axis->SetLabelColor(kBlue);
  axis->Draw();

  leg->Draw();

  c->SaveAs(name);
  
}

#endif
