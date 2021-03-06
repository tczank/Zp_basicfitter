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


void babar_gp_db() {


  //## Loading the different CMS E "Theoretical" Cross section Br and det eff
  TFile * br_fil = new TFile("../merging_energies/Zp_BR.root");
  TFile * nobs_file = new TFile("./babarrebin_newpion.root");
  TFile * upws_xs = new TFile("./Official-weighted-theoretical-cross-section.root");
  // thczank or tczank get your cd pwd first

  // TFile * deteff_plot = new TFile("./weighted_pionveto_on/pion_par/pion_now_deteff.root");
  TFile * deteff_plot = new TFile("../signal_mc_dist_reader/newpionvetodeteff.root");
  // path for detefffit might include a folder or not

  TFile * sqrs_scale = new TFile("./sqrts_scalefit.root");
  //#############################################################//

  //## Setting madgraph cross sections ##//
  TGraph *gr_mu = (TGraph*)br_fil->Get("gr_mu");
  TGraphErrors *upwsxs = (TGraphErrors*) upws_xs->Get("gr_xs_isr");
  TGraph *nobs = (TGraph*)nobs_file->Get("gr_obs");
  TGraph *deteff_fit = (TGraph*)deteff_plot->Get("gr_det_isr");
  // TF1 * deteff_fit = (TF1*) deteff_plot->Get("PrevFitTMP");
  TF1 * sqrs_scale_fit = (TF1*) sqrs_scale->Get("PrevFitTMP");
  // ###############################################################//

  //Output File //
  TFile * f = new TFile("babar_gp_dp_rebin_newpion_jpsiveto.root","RECREATE");
  TH2F * nexp = new TH2F("h_nexp_gp_m", "number of expected events by Z' coupling strength and mass;m_{Z'}[GeV/c^{2}];g';number of expected events;", 2070,0.212,9.99923,10000,-5.,0.0);
  TH1F * nexp_x = new TH1F("h_nexp_gp_m_x", "number of expected events by Z' mass Xproject;m_{Z'}[GeV/c^{2}];number of expected events;", 2070,0.212,9.99923);
  TH1F * nexp_y = new TH1F("h_nexp_gp_m_y", "number of expected events by Z' coupling strength Yproject;g';number of expected events;", 10000,-5.,0.);
  TH2I * gp = new TH2I("h_gp_m_gz", "g' coupling strength by mass and g'z;m_{Z'}[GeV/c^{2}];g';", 2070,0.212,9.99923,10000,-5.,0.0);
  //  TGraph * gr_gp[11064];
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
  // double_t mass_ar[11064]; //mz
  //  double_t gz_ar[10000]; // gp
  double_t vXout, vYout, deteff;

  int xbin = nexp_x->FindBin(0.212125);
  int ybin = nexp_y->FindBin(0.0000);
  //  cout << " the bin corresponding to the muon threshold is " << xbin << endl;

  // THOMAS FROM THE PAST MESSAGE//
  // AFTER SO MANY CHANGES ON THE BG FIT THE NUMBER OF BINS HAS DECREASED, REMEMBER TO CHECK IT BEFORE USING 11064//
  for(int j = 0;j<10000;j++){
    for(int i = 0;i<2070;i++){
      // if(nexp_x->GetBinCenter(i+1) >= 0.212125){
       double_t mass = nexp_x->GetBinCenter(i+1);
      //      double_t gz = exp(nexp_y->GetBinCenter(j+1)*log(10));
      double_t gz = pow(10.,nexp_y->GetBinCenter(j+1));
      //    cout << " the value of i and j is " << i << " " << j << endl;
      //      cout << " mass is " << mass << " and the gz " << gz << endl;
      for(int l = 0; l < 15; l++){
        continuum_norm[l] = continuum_entries[l]/(continuum_entries[0]+continuum_entries[1]+continuum_entries[2]+continuum_entries[3]+continuum_entries[4]+continuum_entries[5]+continuum_entries[6]+continuum_entries[7]+continuum_entries[8]+continuum_entries[9]+continuum_entries[10]+continuum_entries[11]+continuum_entries[12]+continuum_entries[13]+continuum_entries[14]) ;
        continuum_th_lum = 1e3*85.73205*(continuum_norm[l]*upwsxs->Eval(mass));
        continuum_th_lum_all = continuum_th_lum_all + continuum_th_lum;
        //  cout << " the weighted and scaled luminosity is " << continuum_th_lum_all << endl;
      }
      deteff = deteff_fit->Eval(mass);
      if(deteff < 0){ deteff = 0;}
      //double brlumdet = gr_mu->Eval(mass) * deteff * ((up1sxs->Eval(mass)*1e3*4.77836) + (up2sxs->Eval(mass)*1e3*3.5135) + (up3sxs->Eval(mass)*1e3*16.89427) + (up4sxs->Eval(mass)*1e3*690.555) + (up5sxs->Eval(mass)*1e3*123.81655) + continuum_th_lum_all );
      // cout << " the deteff mubr and brlumdet are " << deteff_fit->Eval(mass) << " " << gr_mu->Eval(mass) << " " << brlumdet << endl;
      double brlumdet = gr_mu->Eval(mass) * deteff * 0.92528973*(upwsxs->Eval(mass));
        continuum_th_lum = 0;
        continuum_th_lum_all = 0;
      //double brlumdet =  925.28973*((up1sxs->Eval(mass)*1e3*4.77836/925.28973) + (up2sxs->Eval(mass)*1e3*3.5135/925.28973) + (up3sxs->Eval(mass)*1e3*16.89427/925.28973) + (up4sxs->Eval(mass)*1e3*690.555/925.28973) + (up5sxs->Eval(mass)*1e3*123.81655/925.28973) + continuum_th_lum_all/925.28973 );
        // cout << " the deteff mubr and brlumdet are " << deteff_fit->Eval(mass) << " " << gr_mu->Eval(mass) << " " << brlumdet << endl;
        //# for the previous number of bins on x 11064 the dimuon mass threshold is on the bin # 235
        //  if(i-222>=0){nobs->GetPoint(i-222,vXout,vYout);}
      nobs->GetPoint(i,vXout,vYout);
        // else{nobs->GetPoint(0,vXout,vYout);}
      //  gz = 0.001;
           // cout << " normalized brlumdet by gz " << (brlumdet*pow(gz,2))/pow(0.1,2) << endl;
      double nexp_n = (brlumdet*pow(gz,2))/pow(0.1,2);
      nexp->SetBinContent(i+1,j+1,nexp_n);
      // if(mass > 8.42 && mass < 8.43){nexp_n = 0;}
      //if(mass > 9.64 && mass < 9.66){nexp_n = 0;}
      //if(mass > 9.68 && mass < 9.74){nexp_n = 0;}
      double gp_val = nexp_n/vYout;
      //    cout << " gp_val " << gp_val << " and vYout " << vYout << " and nexp_n " << nexp_n << " and vXout " << vXout << endl;
      if(gp_val >= 1.0 && (mass < 3.05 || mass > 3.13)){
        gp->SetBinContent(i+1,j+1,gp_val);
      }
      else{
        nexp->SetBinContent(i+1,j+1,0);
        gp->SetBinContent(i+1,j+1,0);
      }
      //  if(nexp_n < vYout){gp->SetBinContent(i+1,j+1,0.);}

      //  else gp->Fill(i+1,j+1,gp_val);
      // cout << " the number of expected events is " << nexp_n << endl;
      //  cout << "  g' " << gp_val << endl;
     // mass_ar[i] = vXout;
     //   gz_ar[j] = gz;
     // gr_gp[j]= new TGraph()
     //}
     }
    }

  //  nexp->Draw("contz4");
  //   gPad->SetLogz();

  // Saving Output File //
  //  gp->SetMinimum(1.0);
  gp->SetName("gp");
  gp->Write();

  nexp->SetName("Number_exp_dist");
  nexp->Write();
  f->Write();
  f->Close();
}
