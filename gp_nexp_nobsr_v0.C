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
// Up4S 690.555 (without exp65)
// Up5S 123.81655
// continuum 85.73205


void gp_nexp_nobsr_v0() {


  //## Loading the different CMS E "Theoretical" Cross section Br and det eff
  TFile * nobs_file = new TFile("../fit_reborn/real_merged_xs_all.root");
  TFile * nexp_file = new TFile("./nexp_merge.root");
  //#############################################################//

  //## Setting madgraph cross sections ##//
  TGraph *nobs = (TGraph*)nobs_file->Get("gr_obs");
  TH2F * nexp = (TH2F*)nexp_file->Get("Number_exp_dist");
  // ###############################################################//

  //Output File //
  TFile * f = new TFile("gp.root","RECREATE");
  TH1F * gp = new TH1F("h_gp_m_gz", "Z' coupling strength by mass and g'z;m_{Z'}[GeV/c^{2}];g';", 10000,0.0,10.0);
  //  TH1F * gp_x = new TH1F("h_nexp_gp_m_x", "number of expected events by Z' mass Xproject;m_{Z'}[GeV/c^{2}];number of expected events;", 10000,0.0,9.22);
  TH1F * nexp_y = new TH1F("h_nexp_gp_m_y", "number of expected events by Z' coupling strength Yproject;g';number of expected events;", 10000,0.0,0.2);
  // ########################################################## //

  // Number of Expected Events Number of Observed Events Ratio Calculation //
  double_t x = 0.018500; // reduced di muon threshold mass
  double_t mupdg = 4.*pow(0.1056583745,2); // reduced mass correction to invariant
  double_t mass = sqrt(pow(x,2) + mupdg); //mz
  double_t gz = 0.0001; // gp
  int i = 0;
  int j = 0;

  int xbin = gp->FindBin(0.212125);
   int ybin = nexp_y->FindBin(0.0001);
  //  cout << " the bin corresponding to the muon threshold is " << xbin << endl;

  i = xbin;
  j = ybin;

    while (gz < 1.){
    while( mass < 9.21){
      double gp_val = nobs->Eval(mass)/nexp->GetBinContent(i,j);
      if(gp_val < 1.){gp->SetBinContent(i,gp_val);}
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
  gp->SetName("Nexp_Nobs_ratio_g_dist");
  gp->Write();
  f->Write();
  f->Close();
}
