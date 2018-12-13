#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1F.h"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFile.h"
#include <cmath>



void realdata_test(TString signalfilename, TString signalfilename_other, TString backfilename) {
  //  void realdat_test_v1(TString signalfilename, TString signalfilename_other, TString backfilename, TString backpsimu, TString backpsipi) {
  // TFile * br_fil = new TFile("Zp_BR.root");

  //  TGraph *gr_mu = new TGraph();

  //gr_mu = (TGraph*)br_fil->Get("gr_mu");

  gStyle->SetOptFit();
  TFile *signal = new TFile(signalfilename); // tauskim
  TFile *signal_other = new TFile(signalfilename_other); // hadronbj
  TFile *bg = new TFile(backfilename); // main back 4mu
  //  TFile *bgpsimu = new TFile(backpsimu);
  // TFile *bgpsipi = new TFile(backpsipi);
  TH1F *dpinvmasslm[4];
  TH1F *dpinvmasslm_oth[4];
  TH1F *bginvmasslm[4];
  //TH1F *bgpsimunvmasslm[4];
  // TH1F *bgpsipinvmasslm[4];
  TH1F * reserva;
  

  // TF1 *parawidth = new TF1("parawidth","cheb9",0.5,10);

  for(int i = 3; i < 4; i ++) {
    dpinvmasslm[i] = (TH1F*)signal->Get(TString::Format("h_mycombitrigeff_3"));
   reserva = (TH1F*)bg->Get(TString::Format("h_mycombitrigeff_3"));
    dpinvmasslm_oth[i] = (TH1F*)signal_other->Get(TString::Format("h_mycombitrigeff_3"));
    bginvmasslm[i] = (TH1F*)bg->Get(TString::Format("h_babarpjpsicut_5"));
    // bgpsimunvmasslm[i] = (TH1F*)bgpsimu->Get(TString::Format("h_mycomtrigeff_3"));
    //  bgpsipinvmasslm[i] = (TH1F*)bgpsipi->Get(TString::Format("h_mycomtrigeff_3"));
  }

  TCanvas *C1 = new TCanvas("C1", "", 10, 10, 800, 800);
  for(int i = 3; i < 4; i ++) {
     dpinvmasslm[i]->SetLineColor(1);
     dpinvmasslm_oth[i]->SetLineColor(2);
     bginvmasslm[i]->SetLineColor(3);
     reserva->SetLineColor(4);

     //considering the luminosity of 42 fb{-1}
     //  bginvmasslm[i]->Scale(0.0140928);
     // considering the luminosity of 30.35 fb{-1} (averaged between hadron bj and tauab)
     // bginvmasslm[i]->Scale(0.0101837);
      // considering the luminosity of 58.70 fb{-1} ( sum of hadrob bj and tauab)
     //     bginvmasslm[i]->Scale(0.0196996);

       //disconsidering jpsi background
      // bgpsimunvmasslm[i]->Scale(0.0);
      // bgpsipinvmasslm[i]->Scale(0.0);

      //bginvmasslm[i]->Add(bgpsimunvmasslm[i]);
      //bginvmasslm[i]->Add(bgpsipinvmasslm[i]);

      C1->cd(i + 1);

      
      //dpinvmasslm[i]->Rebin(105);
      // dpinvmasslm_oth[i]->Rebin(105);

      //  dpinvmasslm[i]->Add(dpinvmasslm_oth[i]);
      //   bginvmasslm[i]->Rebin(105);
      //  reserva->Rebin(105);

        //scaling the 4 muon main background for the luminosity of the hadronBJ and tauskims separetly //
        bginvmasslm[i]->Scale(0.0112748);
        reserva->Scale(0.00909591);



          dpinvmasslm[i]->Draw("hist");
          bginvmasslm[i]->Draw("same, hist");
          //    dpinvmasslm_oth[i]->Draw("same, hist");
          //  reserva->Draw("same, hist");
  }


}
