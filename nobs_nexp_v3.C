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


void nobs_nexp_v3() {


  //## Loading the different CMS E "Theoretical" Cross section Br and det eff
  TFile * br_fil = new TFile("../merging_energies/Zp_BR.root");
  TFile * nobs_file = new TFile("./real_pioff_all.root");
  TFile * up1s_xs = new TFile("~tczank/MEGA/MEGAsync/part-phys/rootfiles/newdarkz/gplimproc/xslist_1s.root");
  TFile * up2s_xs = new TFile("~tczank/MEGA/MEGAsync/part-phys/rootfiles/newdarkz/gplimproc/xslist_2s.root");
  TFile * up3s_xs = new TFile("~tczank/MEGA/MEGAsync/part-phys/rootfiles/newdarkz/gplimproc/xslist_3s.root");
  TFile * up4s_xs = new TFile("~tczank/MEGA/MEGAsync/part-phys/rootfiles/newdarkz/gplimproc/madgraphxs_nodecaymode.root");
  TFile * up5s_xs = new TFile("~tczank/MEGA/MEGAsync/part-phys/rootfiles/newdarkz/gplimproc/xslist_5s.root");
  // thczank or tczank get your cd pwd first

  TFile * deteff_plot = new TFile("./pionveto_on_off/pioff_pars/pioff_deteff.root");
  // path for detefffit might include a folder or not

  TFile * sqrs_scale = new TFile("./sqrts_scalefit.root");
  //#############################################################//

  //## Setting madgraph cross sections ##//
  TGraph *gr_mu = (TGraph*)br_fil->Get("gr_mu");
  TGraphErrors *up1sxs = (TGraphErrors*) up1s_xs->Get("Cross_section_Y1S");
  TGraphErrors *up2sxs = (TGraphErrors*) up2s_xs->Get("Cross_section_Y2S");
  TGraphErrors *up3sxs = (TGraphErrors*) up3s_xs->Get("Cross_section_Y3S");
  TGraphErrors *up4sxs = (TGraphErrors*) up4s_xs->Get("madgraphxsZpm_gppone");
  TGraphErrors *up5sxs = (TGraphErrors*) up5s_xs->Get("Cross_section_Y5S");
  TGraph *nobs = (TGraph*)nobs_file->Get("gr_obs");
  TF1 * deteff_fit = (TF1*) deteff_plot->Get("PrevFitTMP");
  TF1 * sqrs_scale_fit = (TF1*) sqrs_scale->Get("PrevFitTMP");
  // ###############################################################//

  //Output File //
  TFile * f = new TFile("real_pioff_nexp_nobs.root","RECREATE");
  TH2D * nexp = new TH2D("h_nexp_gp_m", "number of expected events by Z' coupling strength and mass;m_{Z'}[GeV/c^{2}];g';number of expected events;", 11064,0.0,10.0,10000,-5.,0.0);
  TH1F * nexp_x = new TH1F("h_nexp_gp_m_x", "number of expected events by Z' mass Xproject;m_{Z'}[GeV/c^{2}];number of expected events;", 11064,0.0,10.0);
  TH1F * nexp_y = new TH1F("h_nexp_gp_m_y", "number of expected events by Z' coupling strength Yproject;g';number of expected events;", 10000,-5.,0.);
  TH2D * gp = new TH2D("h_gp_m_gz", "g' coupling strength by mass and g'z;m_{Z'}[GeV/c^{2}];g';", 11064,0.0,10.0,10000,-5.,0.0);
  // ########################################################## //

  // Continuum sample theoretical cross section scaling//
  double_t continuum_sqrts[]={10.5177,10.5176,10.5162,10.5183,10.5193,10.5195,10.5198,10.5204,10.5227,10.5243,10.5763,10.5805,10.5159,10.5153,10.5143};
  double_t continuum_entries[]={2105.48,1293.5,1397.98,1027.81,374.049,27.7626,48.6592,174.039,21.7921,45.674,1,1,218.817,615.852,329.271};
  double_t continuum_norm[14];
  double_t continuum_th_lum;
  double_t continuum_th_lum_all;
   //###############################################

  // Number of Expected Events Calculation //
  double_t x = 0.018500; // reduced di muon threshold mass
  double_t mupdg = 4.*pow(0.1056583745,2); // reduced mass correction to invariant
  // double_t mass = sqrt(pow(x,2) + mupdg); //mz
  // double_t gz = 0.0001; // gp
  double_t vXout, vYout, deteff;

  int xbin = nexp_x->FindBin(0.212125);
  int ybin = nexp_y->FindBin(0.0000);
  //  cout << " the bin corresponding to the muon threshold is " << xbin << endl;


  for(int j = 0;j<10000;j++){
    for(int i = 0;i<11064;i++){
     if(nexp_x->GetBinCenter(i+1) >= 0.212125){
       double_t mass = nexp_x->GetBinCenter(i+1);
      //      double_t gz = exp(nexp_y->GetBinCenter(j+1)*log(10));
      double_t gz = pow(10.,nexp_y->GetBinCenter(j+1));
      // cout << " the value of i and j is " << i << " " << j << endl;
           // cout << " mass is " << mass << " and the gz " << gz << endl;
      for(int l = 0; l < 15; l++){
        continuum_norm[l] = continuum_entries[l]/(continuum_entries[0]+continuum_entries[1]+continuum_entries[2]+continuum_entries[3]+continuum_entries[4]+continuum_entries[5]+continuum_entries[6]+continuum_entries[7]+continuum_entries[8]+continuum_entries[9]+continuum_entries[10]+continuum_entries[11]+continuum_entries[12]+continuum_entries[13]+continuum_entries[14]) ;
        continuum_th_lum = 1e3*85.73205*(continuum_norm[l]*up4sxs->Eval(mass));
        continuum_th_lum_all = continuum_th_lum_all + continuum_th_lum;
        //  cout << " the weighted and scaled luminosity is " << continuum_th_lum_all << endl;
      }
      deteff = deteff_fit->Eval(mass);
      if(deteff < 0){ deteff = 0;}
      double brlumdet = gr_mu->Eval(mass) * deteff * ((up1sxs->Eval(mass)*1e3*4.77836) + (up2sxs->Eval(mass)*1e3*3.5135) + (up3sxs->Eval(mass)*1e3*16.89427) + (up4sxs->Eval(mass)*1e3*690.555) + (up5sxs->Eval(mass)*1e3*123.81655) + continuum_th_lum_all );
      //     cout << " the deteff mubr and brlumdet are " << deteff_fit->Eval(mass) << " " << gr_mu->Eval(mass) << " " << brlumdet << endl;
      if(i-235>=0){nobs->GetPoint(i-235,vXout,vYout);}
      else{nobs->GetPoint(0,vXout,vYout);}
      double nexp_n = (brlumdet*pow(gz,2))/pow(0.1,2);
      if(mass > 8.42 && mass < 8.43){nexp_n = 0;}
      if(mass > 9.64 && mass < 9.66){nexp_n = 0;}
      if(mass > 9.68 && mass < 9.74){nexp_n = 0;}
      double gp_val = nexp_n/vYout;
      if(gp_val >= 1.){
        gp->SetBinContent(i+1,j+1,gp_val);
      }
      // cout << " the number of expected events is " << nexp_n << endl;
      //  cout << "  g' " << gp_val << endl;
      nexp->SetBinContent(i+1,j+1,nexp_n);
      continuum_th_lum = 0;
      continuum_th_lum_all = 0;
         }
    }
    }

  //  nexp->Draw("contz4");
  //   gPad->SetLogz();

  // Saving Output File //
  gp->SetName("Nexp_Nobs_ratio_g_dist");
  gp->Write();

  nexp->SetName("Number_exp_dist");
  nexp->Write();
  f->Write();
  f->Close();
}
