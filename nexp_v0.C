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

void nexp_v0() {

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  // gROOT->SetBatch(1);

  //## Loading the different CMS E "Theoretical" Cross section Br and det eff
  TFile * br_fil = new TFile("../merging_energies/Zp_BR.root");
  TFile * up1s_xs = new TFile("../fit_reborn/real_Up1s_xs_all.root");
  TFile * up3s_xs = new TFile("../fit_reborn/real_Up3s_xs_all.root");
  TFile * up4s_xs = new TFile("../fit_reborn/madgraphxs_nodecaymode.root");
  TFile * up5s_xs = new TFile("../fit_reborn/real_Up5s_xs_all.root");
  TFile * deteff_plot = new TFile("../fit_reborn/db_parfits_new/detefffit.root")
  //#############################################################//

  //## Setting madgraph cross sections ##//
  TGraph *gr_mu = (TGraph*)br_fil->Get("gr_mu");
  TGraphErrors *up1sxs = (TGraphErrors*) up1s_xs->Get("Cross_section_Y1S");
  TGraphErrors *up3sxs = (TGraphErrors*) up3s_xs->Get("Cross_section_Y3S");
  TGraphErrors *up4sxs = (TGraphErrors*) up1s_xs->Get("madgraphxsZpm_gppone");
  TGraphErrors *up5sxs = (TGraphErrors*) up1s_xs->Get("Cross_section_Y5S");
  TGraph * deteff_fit = (TGraph*) deteff_plot->Get("PrevFitTMP")
  // ###############################################################//

  //Output File //
  TFile *f = new TFile("nexp_merge.root","RECREATE");
  TH2F *nexp = new TH2F("h_nexp_gp_m", "number of expected events by Z' coupling strength and mass;m_{Z'}[GeV/c^{2}];g';", 10000,0.0,10.0,10000,0.0001,0.1);
  // ########################################################## //

  // Number of Expected Events Calculation //
  double_t x = 0.018500; // di muon threshold mass
  double_t mupdg = 4.*pow(0.1056583745,2);
  double_t i = sqrt(pow(x,2) + mupdg);
  double_t j = 0.0001;

  while (j < 0.1){
    while( i < 9.31){
      double nexp_n = pow(j,2)
        }
  }


  // Saving Output File //
  nexp->SetName("Number_exp_dist");
  nexp->Write();
  f->Write();
  f->Close();
}
