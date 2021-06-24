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
//#include "CrystalBall.C"

void smear_oldnewfit_paropt(TString signalfilename) {

  TFile * br_fil = new TFile("../../merging_energies/Zp_BR.root");
  TFile * genid_w = new TFile("../../signal_mc_dist_reader/nonisr_smeared_deteff.root");

  TGraph *gr_mu = new TGraph();
  TGraph *gr_w = new TGraph();

  gr_mu = (TGraph*)br_fil->Get("gr_mu");
  gr_w = (TGraph*)genid_w->Get("gr_w");

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);

  //TFile * par = new TFile("./smear_trippar_all.root");
     TFile * par = new TFile("./smear_optpar_all_10.root");

  TGraphErrors* frac1;
  TGraphErrors* frac2;
  TGraphErrors* tgww;
  TGraphErrors* tg1w;
  TGraphErrors* tg2w;
  TGraphErrors* tg3w;

  frac1 = (TGraphErrors*)par->Get("gr_frac1");
  frac2 = (TGraphErrors*)par->Get("gr_frac2");
  tg1w = (TGraphErrors*)par->Get("gr_fw");
  tg2w = (TGraphErrors*)par->Get("gr_sw");
  tg3w = (TGraphErrors*)par->Get("gr_tw");
  tgww = (TGraphErrors*)par->Get("gr_ww");

  // gROOT->SetBatch(1);
  TFile *signal = new TFile(signalfilename);

  TH1F *dpinvmasslm;

  TH1F * genid_redmass;

  TF1 *parawidth = new TF1("parawidth", "pol9",0.,10.);
  TF1 *parawei = new TF1("parawei", "pol5",0., 10.);

  genid_redmass = (TH1F*)signal->Get(TString::Format("h_genidredmu_0"));
  dpinvmasslm = (TH1F*)signal->Get(TString::Format("h_babarpjpsicut_5"));

  TCanvas *C1 = new TCanvas("C1", "", 10, 10, 800, 800);
     dpinvmasslm->SetLineColor(1);
     C1->cd(1);

     dpinvmasslm->GetYaxis()->SetTitleOffset(1.6);
    dpinvmasslm->Draw();

  TF1 * singpartripgaus;
  TF1 * trigtripo;

  double hist_mean;
  double peakwidth;
  double entriesatmean;
  double rms;
  double std_dev;
  double significance;
  double thirdpolchi;
  double gausnobschi;

  double zpgenideff = genid_redmass->GetEntries();

  double tppar[7];
  hist_mean = dpinvmasslm->GetBinCenter(dpinvmasslm->GetMaximumBin());

  tppar[0] = frac1->Eval(hist_mean);
  tppar[1] = hist_mean;
  tppar[2] = tg1w->Eval(hist_mean);
  tppar[3] = frac2->Eval(hist_mean);
  tppar[4] = tg2w->Eval(hist_mean);
  tppar[5] = tg3w->Eval(hist_mean);
  tppar[6] = 1. - tppar[0] - tppar[3];

  peakwidth = tg1w->Eval(hist_mean);
  entriesatmean = dpinvmasslm->GetBinContent(dpinvmasslm->GetMaximumBin());
  rms = dpinvmasslm->GetBinWidth(1);
  std_dev = dpinvmasslm->GetStdDev(1);

  double lowerfit = hist_mean - 100*peakwidth;
  if(lowerfit < 0){lowerfit =0;}
  double higherfit = hist_mean + 300*peakwidth;
  if(higherfit > 10.0){higherfit = 10.;}

      singpartripgaus = new TF1("single_par_triple_gaus","(([0]/((sqrt(2*TMath::Pi()*[2]*[2]))))*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) +(([3])/((sqrt(2*TMath::Pi()*[4]*[4]))))*exp(-((x-[1])*(x-[1]))/(2*[4]*[4])) + (([6])/((sqrt(2*TMath::Pi()*[5]*[5]))))*exp(-((x-[1])*(x-[1]))/(2*[5]*[5])) )", lowerfit, higherfit);
    singpartripgaus->SetLineColor(1);
    // defining parameters names //

 for(int i = 0; i < 7; i++){
      singpartripgaus->SetParameter(i,tppar[i]);
      //  double_crystalball->SetParLimits(i,dbpar[i]/2, dbpar[i]*2);
    }


    singpartripgaus->SetParName(0,"gaus_height1");
    singpartripgaus->SetParName(1,"mean");
    singpartripgaus->SetParName(2,"first_width");
    singpartripgaus->SetParName(3,"gaus_height2");
    singpartripgaus->SetParName(4,"second_width");
     singpartripgaus->SetParName(5,"third_width");
     singpartripgaus->SetParName(6,"gaus_height3");
     singpartripgaus->SetNpx(1000);

     singpartripgaus->SetParLimits(0,0.0,entriesatmean*tppar[0]*2.5);
    singpartripgaus->FixParameter(1,hist_mean);
    singpartripgaus->SetParLimits(2,rms,2.5*tppar[2]);
    singpartripgaus->SetParLimits(3,0.0,entriesatmean*tppar[3]*2.5);
    singpartripgaus->SetParLimits(4,rms,2.5*tppar[4]);
     singpartripgaus->SetParLimits(5,rms,2.5*tppar[5]);
     singpartripgaus->SetParLimits(6,0.0,entriesatmean*tppar[6]*2.5);

     float rangelow = 25;
     float rangehigh = 25;

        if(hist_mean > 7.515 && hist_mean < 7.715){
       rangelow = 80;
       rangehigh = 80;
       singpartripgaus->SetRange(hist_mean-rangelow*peakwidth,hist_mean+rangehigh*peakwidth);
     }
        /*  else{
     rangelow = 5.;
     rangehigh = 5.;
     singpartripgaus->SetRange(hist_mean-rangelow*peakwidth,hist_mean+rangehigh*peakwidth);
     }*/

     singpartripgaus->SetRange(hist_mean-rangelow*peakwidth,hist_mean+rangehigh*peakwidth);


  Int_t nbinsdp;
  Int_t nbinsbg;

  Double_t x;
  Double_t y;

  Double_t z;
  Double_t w;

  Double_t er;

  double par_sigfit[3];
  TH2F * sigres;
  TH2F * sigres_alt;
  TH2F * pull;

  TRandom * r1 = new TRandom();

  TH1D * h_pull[10000];
  TH1D * h_pull_res[2];
  h_pull_res[0]= new TH1D("pull distribution_0","pull;pull;entries;", 100,-5,5);
  h_pull_res[1]= new TH1D("pull distribution_1","pull;pull;entries;", 1000,-200,500);

  //2018/02/23

  //##After fitting the single normalized gaussian, we build a new function with the defined parameters for the single normalized gaussian + a 3rd order polynomial, THEN, we perform a fit of this function over the Background and we take out the new height which is allowed to float, the error of the height and then we divide it by the bin width and multiply it by 1.4*height error

  double normtripg[3];
  double normer[3];
  double heightlist[3];
  double herr[3];

  //   hist_mean = dpinvmasslm->GetBinCenter(dpinvmasslm->GetMaximumBin());
  // std_dev = dpinvmasslm->GetStdDev();
  // peakwidth = gr_w->Eval(hist_mean);


      singpartripgaus->SetNpx(1000);
      TFitResultPtr finalfit = dpinvmasslm->Fit(singpartripgaus,"SOQ+");
      singpartripgaus->ReleaseParameter(1);
      finalfit = dpinvmasslm->Fit(singpartripgaus, "ORMBQ");
      finalfit = dpinvmasslm->Fit(singpartripgaus, "ORMBQ"); 
      finalfit = dpinvmasslm->Fit(singpartripgaus, "ORMBQ"); 
      finalfit = dpinvmasslm->Fit(singpartripgaus, "ORMBQ"); 
      finalfit = dpinvmasslm->Fit(singpartripgaus, "ORMBQ"); 
      finalfit = dpinvmasslm->Fit(singpartripgaus, "ORMBQ"); 
      finalfit = dpinvmasslm->Fit(singpartripgaus, "RMBQS+"); 

      // finalfit->Print();

      heightlist[0] = singpartripgaus->GetParameter(0);
      heightlist[1] = singpartripgaus->GetParameter(3);
      heightlist[2] = singpartripgaus->GetParameter(6);

      herr[0] = singpartripgaus->GetParError(0);
      herr[1] = singpartripgaus->GetParError(3);
      herr[2] = singpartripgaus->GetParError(6);

      normer[0] = sqrt(pow(herr[0]/(heightlist[0] + heightlist[1] + heightlist[2]),2) + pow(heightlist[0]*herr[1]/(pow(heightlist[0] + heightlist[1] + heightlist[2],2)),2) + pow(heightlist[0]*herr[2]/(pow(heightlist[0] + heightlist[1] + heightlist[2],2)),2));

      normer[1] = sqrt(pow(herr[1]/(heightlist[0] + heightlist[1] + heightlist[2]),2) + pow(heightlist[1]*herr[0]/(pow(heightlist[0] + heightlist[1] + heightlist[2],2)),2) + pow(heightlist[1]*herr[2]/(pow(heightlist[0] + heightlist[1] + heightlist[2],2)),2));

      normer[2] = sqrt(pow(herr[2]/(heightlist[0] + heightlist[1] + heightlist[2]),2) + pow(heightlist[2]*herr[1]/(pow(heightlist[0] + heightlist[1] + heightlist[2],2)),2) + pow(heightlist[2]*herr[0]/(pow(heightlist[0] + heightlist[1] + heightlist[2],2)),2));

      for (int l = 0; l < 3 ; l++){
        normtripg[l] = (heightlist[l]/(heightlist[0] + heightlist[1] +heightlist[2]));
        //   cout << " the normtripg is " << normtripg[l] << " and the heightlist is " << heightlist[l] << endl;
      }

      TF1 *only_one_gaus[3];

      for(int l = 0; l < 3; l++){
        only_one_gaus[l] = new TF1(TString::Format("only_one_gaus_%d",l),"gausn(0)",hist_mean-20*peakwidth, hist_mean+20*peakwidth);
      }

      //### John's correction following BN 1229 ################################ //

      double f_one = singpartripgaus->GetParameter(0)/(singpartripgaus->GetParameter(3) + singpartripgaus->GetParameter(0)+singpartripgaus->GetParameter(6));
      double f_two = singpartripgaus->GetParameter(3)/(singpartripgaus->GetParameter(3) + singpartripgaus->GetParameter(0)+singpartripgaus->GetParameter(6));

      double tripww = sqrt(f_one * singpartripgaus->GetParameter(2)*singpartripgaus->GetParameter(2) + f_two * singpartripgaus->GetParameter(4)*singpartripgaus->GetParameter(4) + (1.0 - f_one - f_two)*singpartripgaus->GetParameter(5)*singpartripgaus->GetParameter(5));

       double triperww = tripww * sqrt(TMath::Power(singpartripgaus->GetParError(2),2) + TMath::Power(singpartripgaus->GetParError(4),2) + TMath::Power(singpartripgaus->GetParError(5),2));
      double tripwh = sqrt((TMath::Power(singpartripgaus->GetParameter(0)*singpartripgaus->GetParameter(2),2)+TMath::Power(singpartripgaus->GetParameter(3)*singpartripgaus->GetParameter(4),2) + TMath::Power(singpartripgaus->GetParameter(6)*singpartripgaus->GetParameter(5),2))/(singpartripgaus->GetParameter(2)*singpartripgaus->GetParameter(2) + singpartripgaus->GetParameter(4)*singpartripgaus->GetParameter(4) +  singpartripgaus->GetParameter(5)*singpartripgaus->GetParameter(5)));
      double triperwh = tripww * sqrt(TMath::Power(singpartripgaus->GetParError(0),2) + TMath::Power(singpartripgaus->GetParError(4),2) + TMath::Power(singpartripgaus->GetParError(5),2));



      double tripger[] = {triperwh,peakwidth,triperww};


    //naming triple gaussian parameters

    double firstwidth = singpartripgaus->GetParameter(2);
    double firstwidther = singpartripgaus->GetParError(2);

    double secondwidth = singpartripgaus->GetParameter(4);
    double secondwidther = singpartripgaus->GetParError(4);

    double thirdwidth = singpartripgaus->GetParameter(5);
    double thirdwidther = singpartripgaus->GetParError(5);

    TAxis * xaxis = dpinvmasslm->GetXaxis();
    TAxis * yaxis = dpinvmasslm->GetYaxis();
    Int_t binxl = xaxis->FindBin(singpartripgaus->GetParameter(1)-rangelow*peakwidth);
    Double_t binxlerror = dpinvmasslm->GetBinError(binxl);
    Int_t binxh = xaxis->FindBin( singpartripgaus->GetParameter(1)+rangehigh*peakwidth);
     
    dpinvmasslm->GetXaxis()->SetRange(binxl,binxh);


    TString signalplotname = signalfilename + string(".eps");
    C1->Print(signalplotname);



       //to parametrize the triple gaussian completely
    cout << hist_mean << " " << firstwidth << " " << firstwidther << " " << secondwidth << " " << secondwidther << " " << thirdwidth << " " << thirdwidther << " " << normtripg[0] << " " << normer[0] << " " << normtripg[1] << " " << normer[1] << " " << 1. - normtripg[0] - normtripg[1] << " " << normer[2] << " " << tripww << " " << triperww << endl;

}
