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


  //## Loading the different CMS E "Theoretical" Cross section Br and det eff
  TFile * br_fil = new TFile("../merging_energies/Zp_BR.root");
  TFile * up1s_xs = new TFile("/home/tczank/MEGA/MEGAsync/part-phys/rootfiles/newdarkz/gplimproc/xslist_1s.root");
  TFile * up2s_xs = new TFile("/home/tczank/MEGA/MEGAsync/part-phys/rootfiles/newdarkz/gplimproc/xslist_2s.root");
  TFile * up3s_xs = new TFile("/home/tczank/MEGA/MEGAsync/part-phys/rootfiles/newdarkz/gplimproc/xslist_3s.root");
  TFile * up4s_xs = new TFile("/home/tczank/MEGA/MEGAsync/part-phys/rootfiles/newdarkz/gplimproc/madgraphxs_nodecaymode.root");
  TFile * up5s_xs = new TFile("/home/tczank/MEGA/MEGAsync/part-phys/rootfiles/newdarkz/gplimproc/xslist_5s.root");
  TFile * deteff_plot = new TFile("../fit_reborn/db_parfits_new/detefffit.root");
  TFile * sqrs_scale = new TFile("./sqrts_scalefit.root");
  //#############################################################//

  //## Setting madgraph cross sections ##//
  TGraph *gr_mu = (TGraph*)br_fil->Get("gr_mu");
  TGraphErrors *up1sxs = (TGraphErrors*) up1s_xs->Get("Cross_section_Y1S");
  TGraphErrors *up2sxs = (TGraphErrors*) up2s_xs->Get("Cross_section_Y2S");
  TGraphErrors *up3sxs = (TGraphErrors*) up3s_xs->Get("Cross_section_Y3S");
  TGraphErrors *up4sxs = (TGraphErrors*) up4s_xs->Get("madgraphxsZpm_gppone");
  TGraphErrors *up5sxs = (TGraphErrors*) up5s_xs->Get("Cross_section_Y5S");
  TF1 * deteff_fit = (TF1*) deteff_plot->Get("PrevFitTMP");
  TF1 * sqrs_scale_fit = (TF1*) sqrs_scale->Get("PrevFitTMP");
  // ###############################################################//

  //Output File //
  TFile * f = new TFile("nexp_merge.root","RECREATE");
  TH2F * nexp = new TH2F("h_nexp_gp_m", "number of expected events by Z' coupling strength and mass;m_{Z'}[GeV/c^{2}];g';number of expected events;", 10000,0.0,10.0,10000,0.0001,1.0);
  TH1F * nexp_x = new TH1F("h_nexp_gp_m_x", "number of expected events by Z' mass Xproject;m_{Z'}[GeV/c^{2}];number of expected events;", 10000,0.0,9.22);
  TH1F * nexp_y = new TH1F("h_nexp_gp_m_y", "number of expected events by Z' coupling strength Yproject;g';number of expected events;", 10000,0.0,0.2);
  // ########################################################## //

  // Continuum sample theoretical cross section scaling//
  double_t continuum_sqrts[]={10.5177,10.5176,10.5162,10.5183,10.5193,10.5195,10.5198,10.5204,10.5227,10.5243,10.5763,10.5805,10.5159,10.5153,10.5143};
  double_t continuum_entries[]={2105.48,1293.5,1397.98,1027.81,374.049,27.7626,48.6592,174.039,21.7921,45.674,1,1,218.817,615.852,329.271};
  double_t continuum_norm[14];
  double_t continuum_th_lum;
  double_t continuum_th_lum_all;
  for(int l = 0; l < 15; l++){
    continuum_norm[l] = continuum_entries[l]/continuum_entries[0];
    continuum_th_lum = 1e3*85.73205*(continuum_norm[l]*sqrs_scale_fit->Eval(continuum_sqrts[l]));
    continuum_th_lum_all = continuum_th_lum_all + continuum_th_lum;
    // cout << " the weighted and scaled luminosity is " << continuum_th_lum_all << endl;
  }
  //###############################################

  // Number of Expected Events Calculation //
  double_t x = 0.018500; // reduced di muon threshold mass
  double_t mupdg = 4.*pow(0.1056583745,2); // reduced mass correction to invariant
  double_t mass = sqrt(pow(x,2) + mupdg); //mz
  double_t gz = 0.0001; // gp
  int i = 0;
  int j = 0;

  int xbin = nexp_x->FindBin(0.212);
  int ybin = nexp_y->FindBin(0.0001);
  //  cout << " the bin corresponding to the muon threshold is " << xbin << endl;

  i = xbin;
  j = ybin;

    while (gz < 1.){
    while( mass < 9.21){
      //  cout << " the value of i and j is " << i << " " << j << endl;
      double brlumdet = gr_mu->Eval(mass) * deteff_fit->Eval(mass) * ((up1sxs->Eval(mass)*1e3*4.77836) + (up2sxs->Eval(mass)*1e3*3.5135) + (up3sxs->Eval(mass)*1e3*16.89427) + (up4sxs->Eval(mass)*1e3*690.555) + (up5sxs->Eval(mass)*1e3*123.81655) + continuum_th_lum_all );
      //       cout << " the deteff mubr and brlumdet are " << deteff_fit->Eval(mass) << " " << gr_mu->Eval(mass) << " " << brlumdet << endl;
      double nexp_n = pow(gz,2)*brlumdet;
      //  cout << " the number of expected events is " << nexp_n << endl;
      nexp->SetBinContent(i,j,nexp_n);
      mass = mass + 0.001;
      i = i + 1;
    }
    mass = sqrt(pow(x,2) + mupdg);
    i = xbin;
    gz = gz + 0.001;
    j = j + 1;
    }

  //  nexp->Draw("contz4");
  //   gPad->SetLogz();

  // Saving Output File //
  nexp->SetName("Number_exp_dist");
  nexp->Write();
  f->Write();
  f->Close();
}
