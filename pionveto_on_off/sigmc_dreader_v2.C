#include "TFile.h"
#include "TString.h"



void sigmc_dreader_v2 () {

  TString isr_name[] = {"../cor_signal_mc_isr/darkz_isr_212_all.root", "../cor_signal_mc_isr/darkz_isr_712_all.root","../cor_signal_mc_isr/darkz_isr_1212_all.root","../cor_signal_mc_isr/darkz_isr_1712_all.root","../cor_signal_mc_isr/darkz_isr_2212_all.root","../cor_signal_mc_isr/darkz_isr_2712_all.root","../cor_signal_mc_isr/darkz_isr_3212_all.root","../cor_signal_mc_isr/darkz_isr_3712_all.root","../cor_signal_mc_isr/darkz_isr_4212_all.root","../cor_signal_mc_isr/darkz_isr_4712_all.root","../cor_signal_mc_isr/darkz_isr_5212_all.root", "../cor_signal_mc_isr/darkz_isr_5712_all.root","../cor_signal_mc_isr/darkz_isr_6212_all.root","../cor_signal_mc_isr/darkz_isr_6712_all.root","../cor_signal_mc_isr/darkz_isr_7212_all.root","../cor_signal_mc_isr/darkz_isr_7712_all.root","../cor_signal_mc_isr/darkz_isr_8212_all.root","../cor_signal_mc_isr/darkz_isr_8712_all.root","../cor_signal_mc_isr/darkz_isr_9212_all.root","../cor_signal_mc_isr/darkz_isr_9712_all.root","../cor_signal_mc_isr/darkz_isr_10000_all.root"};


  TH1F * redmu_isr[21];
  TH1F * redmu_isr_4mu[21];

  TFile * sig_isrfile[21];

  double  redmass_isr[21];
  double redmass_isr_4mu[21];
  double redmass_isr_width[21];
  double redmass_isr_width_4mu[21];

  for (int i = 0; i < 21; i++){
    sig_isrfile[i] = new TFile(isr_name[i]);
    redmu_isr[i] = (TH1F*)sig_isrfile[i]->Get("h_genidredmu_0");
    redmu_isr_4mu[i] = (TH1F*)sig_isrfile[i]->Get("h_mycombitrigeff_3");
    //cout << " for the file " << isr_name[i] << " number " << i << " the zp resonance is in the reduced mass of " << redmu_isr[i]->GetBinCenter(redmu_isr[i]->GetMaximumBin()) << endl;
    redmass_isr[i] = redmu_isr[i]->GetBinCenter(redmu_isr[i]->GetMaximumBin());
    redmass_isr_width[i] = redmu_isr[i]->GetStdDev();

    redmass_isr_4mu[i] = redmu_isr_4mu[i]->GetBinCenter(redmu_isr_4mu[i]->GetMaximumBin());
    redmass_isr_width_4mu[i] = redmu_isr_4mu[i]->GetStdDev();

  }

  TGraph * gr_w_isr = new TGraph(21,redmass_isr, redmass_isr_width);
  gr_w_isr->SetTitle("Evolution of tru tag Z' width by its reduced mass");
  gr_w_isr->SetName("gr_w_isr");

  //  TGraph * gr_w_isr_4mu = new TGraph(21,redmass_isr_4mu, redmass_isr_width_4mu);
  // gr_w_isr_4mu->SetTitle("Evolution of tru tag Z' width by its reduced mass");
  // gr_w_isr_4mu->SetName("gr_w_isr_4mu");

  TFile * dists = new TFile("all_4mu_redmu_cor_isr.root", "RECREATE");
  gr_w_isr->Write();
  //  gr_w_isr_4mu->Write();

   for(int i = 0; i < 21; i++){
     // redmu_isr[i]->SetName(TString::Format("genid_redmu_isr_%d",i));
     // redmu_isr[i]->Write();
     redmu_isr_4mu[i]->SetName(TString::Format("redmu_4mu_isr_%d",i));
     redmu_isr_4mu[i]->Write();
   }


  dists->Write();
  dists->Close();

}
