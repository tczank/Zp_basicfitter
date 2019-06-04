//C++
#include <iostream>
#include <fstream>
using namespace std;

//ROOT
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TText.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TGaxis.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLorentzRotation.h"
#include <TGenPhaseSpace.h>

//Dark Style
#include "darklib/NiceGr.c"
#include "darklib/Nice1D.c"
#include "darklib/Nice2D.c"

//External checkxs reader

//## 2018/03/22
// currently used program to evaluate the gp coupling between muons and zp, based on the cross section 90% upper limit from the simulated signal in madgraph with the biggest background sample aafh mumumumu
// last update was use of the "theoretical" cross section, output cross section from madgraph 5, this time without specifying the decay channel of zp, previously set as zp > mu+ mu-, which gave us the incomplete "theoretical" cross section


void gplim_db_realdat_v0() {


  TFile *babargp = new TFile("babargprime.root");
  TGraph *babares = new TGraph();
  TFile *babarg2 = new TFile("babarg2.root");
  TGraph *babarg2p = new TGraph();
  TFile *trident = new TFile("trident.root");
  TGraph *tridentp = new TGraph();


  babares = (TGraph*)babargp->Get("BaBar");
  babares->SetLineColor(2);
  babarg2p = (TGraph*)babarg2->Get("BaBar");
  babarg2p->SetLineColor(3);
  tridentp = (TGraph*)trident->Get("BaBar");
  tridentp->SetLineColor(4);





  cout << endl << "Plot 90% CL upper limit on cross section" << endl;
  cout <<         "========================================" << endl;

  gROOT->Reset();

  gROOT->SetStyle("Bold");
  gStyle->SetCanvasColor(0);
  gStyle->SetLabelColor(1);
  gStyle->SetLabelColor(1,"Y");
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineWidth(1);
  gStyle->SetNdivisions(505);
  gStyle->SetNdivisions(505,"Y");
  gStyle->SetHistFillColor(999);
  gROOT->SetStyle("Plain");  // white as bg
  gStyle->SetOptStat("111111");
  gStyle->SetFrameLineWidth(1);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleBorderSize(0);

  //Definitions
  Double_t smallBetween1 = .1;
  Double_t smallBetween2 = .1;
  Double_t smallBetween3 = .1;
  Double_t smallBetween4 = .1;

  Double_t small = .00001;

  TLine TLine;
  TLatex *t = new TLatex();
  t->SetTextSize(0.08);
  t->SetTextFont(42);
  t->SetTextAlign(12);
  t->SetTextColor(1);
  t->SetTextFont(12);

  TCanvas *C1;

  TPad *SmallC;
  TGaxis *XAxis,*YAxis;
  TLatex XTitle,YTitle,XtopTitle;

  TLegend *legend;

  ifstream in;
  ofstream out;

  //// READ BR
  TFile * ifile = new TFile("../../MEGA/MEGAsync/part-phys/rootfiles/newdarkz/Zp_BR.root");
  TGraph * gr_mu = (TGraph*) ifile->Get("gr_mu");
  TGraph * gr_tau = (TGraph*) ifile->Get("gr_tau");
  TGraph * gr_nu = (TGraph*) ifile->Get("gr_nu");

  //// READ CROSS SECTION X BR UL
  TFile * xsth = new TFile("madgraphxs_nodecaymode.root"); // no zp decay mode specified
  //TFile * newif = new TFile("deteff_all.root");
  TFile * newif = new TFile("real_onres_xs_all.root");


  TGraphErrors *newxsbr = (TGraphErrors*) newif->Get("gr_xs");

  TGraphErrors *madZpxs = (TGraphErrors*) xsth->Get("madgraphxsZpm_gppone");
 newxsbr->Sort();



  //// Create TF1 of the theoretical cross section
  double gp_th = 0.1;

  //// Build xs and g UL
  //// Loop on mass

  TGraph * gr_xs = new TGraph();
  TGraph * gr_gp = new TGraph();
  TGraph * gr_modelind = new TGraph();


  //double masscand[107] = {0.212, 0.22, 0.23, 0.24, 0.25, 0.3, 0.41, 0.42, 0.43, 0.44, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1};
  //double xsval[107];

cout << " the number of entries in newxsbr " << newxsbr->GetN() << endl;
 double xsval[90013];
 // double masscand[10013];

 double masscand = 0.212;
 int i = 0;

 while (masscand < 9.3){
   xsval[i] = newxsbr->Eval(masscand);
    double Zpmass = masscand;
    double xsbr_val = xsval[i];
    double br = gr_mu->Eval(Zpmass);
    //double xs_val = xsbr_val * br;
    double xs_val = xsbr_val;
    double xs_th = (madZpxs->Eval(Zpmass)) * 1e6; //pb -> ab
    double gp_val = gp_th * sqrt(xs_val  / xs_th);
    gr_xs->SetPoint(i, Zpmass, newxsbr->Eval(masscand));
    gr_gp->SetPoint(i, Zpmass, gp_val);
    // gr_modelind->SetPoint(i, Zpmass, gp_val_ind);
    cout << "The mass for Zp is = " << Zpmass << " The corresponding BR is = " << br << " the old xs is = " << xsbr_val << " the corrected xs is = " << xs_val << " the theoretical xs = " << xs_th  << " gp is  = " << gp_val << endl;
    i++;
    if(masscand < 1.){masscand = masscand + 0.0005;}
    else {masscand = masscand + 0.001;}
 }

  //  for (int i = 0; i < newxsbr->GetN(); i ++) {
 /* for (int i = 0;i < 9461; i ++) {
    xsval[i] = newxsbr->Eval(masscand[i]);
    double Zpmass = masscand[i];
    double xsbr_val = xsval[i];
    double br = gr_mu->Eval(Zpmass);
        double xs_val = xsbr_val * br;
	//    double xs_val_ind = xsbr_val;
	// "CORRECTED" CROSS SECTION
	//	double xs_th = (madZpxs->Eval(Zpmass)/10.) * 1e6; //pb -> ab
	double xs_th = (madZpxs->Eval(Zpmass)) * 1e6; //pb -> ab
    double gp_val = gp_th * sqrt(xs_val  / xs_th);
    //  double gp_val_ind = gp_th * sqrt(xs_val_ind / xs_th);
    // gr_xs->SetPoint(i, Zpmass, xsbr_val);
    gr_xs->SetPoint(i, Zpmass, newxsbr->Eval(masscand[i]));
    gr_gp->SetPoint(i, Zpmass, gp_val);
    // gr_modelind->SetPoint(i, Zpmass, gp_val_ind);
    cout << "The mass for Zp is = " << Zpmass << " The corresponding BR is = " << br << " the old xs is = " << xsbr_val << " the corrected xs is = " << xs_val << " the theoretical xs = " << xs_th  << " gp is  = " << gp_val << endl;
    }*/

 C1 = new TCanvas("teste","teste",10,10,1000,1000);
 C1->SetLogy();
 C1->SetLogx();
 gr_gp->GetYaxis()->SetRangeUser(1e-4,1e-1);
 gr_gp->GetXaxis()->SetRangeUser(1e-1,10);
 gr_gp->SetLineColor(5);
 gr_gp->Draw();
 babares->SetLineColor(1);
 babares->Draw("same");
 babarg2p->Draw("same");
 tridentp->Draw("same");
 //   gr_modelind->SetLineColor(3);
 // gr_modelind->Draw("same");
 C1->Update();

  TH2F * Draw_br = new TH2F("Draw_br","",10,0.,11.99,10,0.,0.99);
  Nice2D(Draw_br,0.05,0.05,42,505,1.1,1.2,"","#font[42]{m_{Z'} [GeV/#it{c}^{2}]}","#font[42]{BR}");

  smallBetween1 = .15;
  smallBetween2 = .05;
  smallBetween3 = .05;
  smallBetween4 = .15;

  C1 = new TCanvas("Z_BR","Z_BR",10,10,800,800);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);

  Draw_br->Draw();

  gr_nu->SetMarkerStyle(20);
  gr_nu->SetMarkerSize(1.2);
  gr_nu->SetMarkerColor(1);
  gr_nu->SetLineColor(1);
  gr_nu->SetLineWidth(2);

  gr_mu->SetMarkerStyle(21);
  gr_mu->SetMarkerSize(1.2);
  gr_mu->SetMarkerColor(2);
  gr_mu->SetLineColor(2);
  gr_mu->SetLineWidth(2);

  gr_tau->SetMarkerStyle(22);
  gr_tau->SetMarkerSize(1.2);
  gr_tau->SetMarkerColor(4);
  gr_tau->SetLineColor(4);
  gr_tau->SetLineWidth(2);

  gr_nu->Draw("csame");
  gr_mu->Draw("csame");
  gr_tau->Draw("csame");

  legend=new TLegend(0.7,0.75,0.93,0.93);
  legend->AddEntry(gr_nu,"#font[42]{Z' #rightarrow #nu_{l}#bar{#nu}_{l}}","l");
  legend->AddEntry(gr_mu,"#font[42]{Z' #rightarrow #mu^{-}#mu^{+}}","l");
  legend->AddEntry(gr_tau,"#font[42]{Z' #rightarrow #tau^{-}#tau^{+}}","l");
  legend->SetFillColor(0);
  legend->SetTextFont(22);
  legend->SetTextSize(.05);
  legend->SetLineColor(0);
  legend->Draw("same");

  C1->Print("Z_BR.eps");

   TH2F * DrawXS_BR = new TH2F("DrawXS_BR", "", 20, 0., 10.19, 20, 0., 1200.00);
   Nice2D(DrawXS_BR, 0.05, 0.05, 42, 505, 1.2, 1.5, "", "#font[42]{#it{m}_{Z'} [GeV/#it{c}^{2}]}","#font[42]{#sigma x BR [ab]}");

  C1 = new TCanvas("XS_times_BR_Plots","XS_times_BR_Plots",10,10,900,900);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);

  DrawXS_BR->Draw();

 newxsbr->SetLineColor(1);
  newxsbr->Draw("CS");

   C1->Print("Z_XS_times_BR.eps");

  TH2F * DrawXS = new TH2F("DrawXS", "", 20, 0., 10.11, 20, -100., 2000.00);
  Nice2D(DrawXS, 0.05, 0.05, 42, 505, 1.2, 1.5, "", "#font[42]{#it{m}_{Z'} [GeV/#it{c}^{2}]}","#font[42]{#sigma [ab]}");

  C1 = new TCanvas("XS_Plots","XS_Plots",10,10,900,900);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);

  DrawXS->Draw();

  gr_xs->SetLineColor(1);
  gr_xs->Draw("CS");

  C1->Print("Z_XS.eps");

    TH2F * Draw_GP = new TH2F("Draw_GP", "", 20, 0., 10.19, 20, 1e-5, 5);
   Nice2D(Draw_GP, 0.05, 0.05, 42, 505, 1.2, 1.5, "", "#font[42]{#it{m}_{Z'} [GeV/#it{c}^{2}]}","#font[42]{g'}");

  C1 = new TCanvas("gp_Plots","gp_Plots",10,10,900,900);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);

  Draw_GP->Draw();

  gr_gp->SetLineColor(1);
  gr_gp->Draw("CS");
  babares->Draw("same");
  //  gr_modelind->Draw("CS same");

  gPad->SetLogy();

  C1->Print("Z_GP.eps");

}
