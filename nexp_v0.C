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

//2019/06/07 Merging different CMS E Belle samples procedure
// Calculation of the total number of expected events by "guessing" different gZ'
// scanning from 10^-4 to 1

// Luminosities list fb-1
// Up1S 4.77836
// Up2S 3.5135
// Up3S 16.89427
// Up4S 690.555
// Up5S 123.81655
// continuum 85.73205


void nexp_v0() {

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  // gROOT->SetBatch(1);

  //## Loading the different CMS E "Theoretical" Cross section Br and det eff
  TFile * br_fil = new TFile("../merging_energies/Zp_BR.root");
  TFile * up1s_xs = new TFile("/home/tczank/MEGA/MEGAsync/part-phys/rootfiles/newdarkz/gplimproc/xslist_1s.root");
  TFile * up3s_xs = new TFile("/home/tczank/MEGA/MEGAsync/part-phys/rootfiles/newdarkz/gplimproc/xslist_3s.root");
  TFile * up4s_xs = new TFile("/home/tczank/MEGA/MEGAsync/part-phys/rootfiles/newdarkz/gplimproc/madgraphxs_nodecaymode.root");
  TFile * up5s_xs = new TFile("/home/tczank/MEGA/MEGAsync/part-phys/rootfiles/newdarkz/gplimproc/xslist_5s.root");
  TFile * deteff_plot = new TFile("../fit_reborn/db_parfits_new/detefffit.root");
  //#############################################################//

  //## Setting madgraph cross sections ##//
  TGraph *gr_mu = (TGraph*)br_fil->Get("gr_mu");
  TGraphErrors *up1sxs = (TGraphErrors*) up1s_xs->Get("Cross_section_Y1S");
  TGraphErrors *up3sxs = (TGraphErrors*) up3s_xs->Get("Cross_section_Y3S");
  TGraphErrors *up4sxs = (TGraphErrors*) up4s_xs->Get("madgraphxsZpm_gppone");
  TGraphErrors *up5sxs = (TGraphErrors*) up5s_xs->Get("Cross_section_Y5S");
  TF1 * deteff_fit = (TF1*) deteff_plot->Get("PrevFitTMP");
  // ###############################################################//

  //Output File //
  TFile *f = new TFile("nexp_merge.root","RECREATE");
  TH3F *nexp = new TH3F("h_nexp_gp_m", "number of expected events by Z' coupling strength and mass;m_{Z'}[GeV/c^{2}];g';numer of expected events;", 10000,0.0,10.0,10000,0.0001,1.,1000,1e-4,10);
  // ########################################################## //

  // Number of Expected Events Calculation //
  double_t x = 0.018500; // di muon threshold mass
  double_t mupdg = 4.*pow(0.1056583745,2);
  double_t i = sqrt(pow(x,2) + mupdg); //mz
  double_t j = 0.0001; // gp

  while (j < 0.1){
    while( i < 9.21){
      cout << " the value of i and j is " << i << " " << j << endl;
      double brlumdet = gr_mu->Eval(i) * deteff_fit->Eval(i) * ((up1sxs->Eval(i)*1e3*4.77836) +(up3sxs->Eval(i)*1e3*16.89427) + (up4sxs->Eval(i)*1e3*690.555) + (up5sxs->Eval(i)*1e3*123.81655)  );
      cout << " the deteff mubr and brlumdet are " << gr_mu->Eval(i) << " " << deteff_fit->Eval(i) << " " << brlumdet << endl;
      double nexp_n = pow(j,2)*brlumdet;
      cout << " the number of expected events is " << nexp_n << endl;
      nexp->Fill(i,j,nexp_n);
      i = i + 0.001;
    }
    i = sqrt(pow(x,2) + mupdg);
    j = j + 0.001;
  }


  // Saving Output File //
  nexp->SetName("Number_exp_dist");
  nexp->Write();
  f->Write();
  f->Close();
}
