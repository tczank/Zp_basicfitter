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

//2019/03/25
//## Working fit of a double crystal ball over the isr Z' signal mc samples

void dcball_paropt(TString signalfilename) {

  TFile * br_fil = new TFile("Zp_BR.root");
  TFile * genid_isr_w = new TFile("smear_skim_deteff.root");

  TGraph *gr_mu = new TGraph();
  TGraph *gr_isr_w = new TGraph();

  // TFile * dbin = new TFile("dcball_par.root");
   TFile * dbin = new TFile("dcball_optpar_8.root");

  TGraphErrors *dbfrac1;
  TGraphErrors *dbfrac2;
  TGraphErrors *dbww;
  TGraphErrors *db1w;
  TGraphErrors *db2w;
  TGraphErrors *db1al;
  TGraphErrors *db2al;
  TGraphErrors *db1n;
  TGraphErrors *db2n;
  TGraph *deteff;

  double dbpar[10];

  dbfrac1 = (TGraphErrors *)dbin->Get("gr_dbfrac1");
  dbfrac2 = (TGraphErrors *)dbin->Get("gr_dbfrac2");
  db1w = (TGraphErrors *)dbin->Get("gr_fw");
  db2w = (TGraphErrors *)dbin->Get("gr_sw");
  db1al = (TGraphErrors *)dbin->Get("gr_fal");
  db2al = (TGraphErrors *)dbin->Get("gr_sal");
  db1n = (TGraphErrors *)dbin->Get("gr_fn");
  db2n = (TGraphErrors *)dbin->Get("gr_sn");
  dbww = (TGraphErrors *)dbin->Get("gr_ww");
  deteff = (TGraph *)genid_isr_w->Get("gr_det_isr");

  gr_mu = (TGraph*)br_fil->Get("gr_mu");
  gr_isr_w = (TGraph*)genid_isr_w->Get("gr_w");

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  // gROOT->SetBatch(1);
  TFile *signal = new TFile(signalfilename);

  TH1F *dpinvmasslm;

  TH1F * genid_invmass;


    genid_invmass = (TH1F*)signal->Get("h_babarpjpsicut_5");

    dpinvmasslm = genid_invmass;

  TCanvas *C1 = new TCanvas("C1", "", 100, 100, 1400, 1400);
     dpinvmasslm->SetLineColor(1);

      C1->cd(1);

      dpinvmasslm->GetYaxis()->SetTitleOffset(1.6);
    dpinvmasslm->Draw();

    TF1 * double_crystalball;

  double hist_mean;
  double peakwidth;
  double entriesatmean;
  double rms;
  double std_dev;
  double significance;

  double zpgenideff = genid_invmass->GetEntries();

  hist_mean = dpinvmasslm->GetBinCenter(dpinvmasslm->GetMaximumBin());
  // peakwidth = gr_isr_w->Eval(hist_mean);
  entriesatmean = dpinvmasslm->GetBinContent(dpinvmasslm->GetMaximumBin());
  rms = dpinvmasslm->GetBinWidth(1);
  std_dev = dpinvmasslm->GetStdDev(1);

    dbpar[0] = dbfrac1->Eval(hist_mean);
  dbpar[1] = hist_mean;
  dbpar[2] = db1w->Eval(hist_mean);
  dbpar[3] = db1al->Eval(hist_mean);
  dbpar[4] = db1n->Eval(hist_mean);
  dbpar[5] =  1. - dbpar[0];
  dbpar[6] =  hist_mean;
  dbpar[7] = db2w->Eval(hist_mean);
  dbpar[8] = db2al->Eval(hist_mean);
  dbpar[9] = db2n->Eval(hist_mean);

  peakwidth = dbww->Eval(hist_mean);
  double peakwidth_2 = dbpar[7];
  double al_1 = dbpar[3];
  double n_1 = dbpar[4];
  double al_2 = dbpar[8];
  double n_2 = dbpar[9];


  double lowerfit;
  double higherfit;

  if(hist_mean < 0.212){
    lowerfit = hist_mean -1*peakwidth;
    if(lowerfit < 0){lowerfit =0;}
    higherfit = hist_mean + 3*peakwidth;
    if(higherfit > 10.0){higherfit = 10.;}
  }

  else{
    lowerfit = hist_mean -10*peakwidth;
    if(lowerfit < 0){lowerfit =0;}
    higherfit = hist_mean + 20*peakwidth;
    if(higherfit > 10.0){higherfit = 10.;}
  }


    double_crystalball = new TF1("double_crystalball", "crystalball(0) + crystalball(5) ", hist_mean - 500*peakwidth , hist_mean + 500*peakwidth);

     for(int i = 0; i < 10; i++){
      double_crystalball->SetParameter(i,dbpar[i]);
      double_crystalball->SetParLimits(i,dbpar[i]/10, dbpar[i]*10);
    }

    double_crystalball->SetParName(0,"Constant_1");
    double_crystalball->SetParName(1,"Mean_1");
    double_crystalball->SetParName(2,"Sigma_1");
    double_crystalball->SetParName(3,"Alpha_1");
    double_crystalball->SetParName(4,"N_1");
    double_crystalball->SetParName(5,"Constant_2");
    double_crystalball->SetParName(6,"Mean_2");
    double_crystalball->SetParName(7,"Sigma_2");
    double_crystalball->SetParName(8,"Alpha_2");
    double_crystalball->SetParName(9,"N_2");

    double_crystalball->SetNpx(1000);

    double_crystalball->SetParLimits(0,entriesatmean*dbpar[0],entriesatmean*dbpar[0]+500 );
    double_crystalball->FixParameter(1,hist_mean);
    double_crystalball->FixParameter(6,hist_mean);
    double_crystalball->SetParLimits(2,peakwidth/2,peakwidth*6);
    double_crystalball->SetParLimits(3,-5.,0);
    double_crystalball->SetParLimits(4,0.,8);
    double_crystalball->SetParLimits(5,entriesatmean*dbpar[5],entriesatmean*dbpar[5]+500 );
    double_crystalball->SetParLimits(7,peakwidth/2,peakwidth*6);
    double_crystalball->SetParLimits(8,0.0,5);
    double_crystalball->SetParLimits(9,-20.,20);


    double_crystalball->SetLineColor(4);
    double_crystalball->SetRange(hist_mean - 25*peakwidth, hist_mean + 25*peakwidth);


  Int_t nbinsdp;
  Int_t nbinsbg;

  Double_t x;
  Double_t y;

  Double_t z;
  Double_t w;

  Double_t er;

TF1 * opt_dbcball = new TF1("optmized double_crystalball", "crystalball(0) + crystalball(5) ", 0. , 10.5);
 opt_dbcball->SetRange(lowerfit,  higherfit);


  // here it is ok
      // TFitResultPtr rebornfit = dpinvmasslm->Fit(triplegexp,"RBQS+");
      TFitResultPtr cballfit = dpinvmasslm->Fit(double_crystalball,"S0Q+");
      double chi2 = cballfit->Chi2();

        for (int j = 0; j < 10 ; j++){
        opt_dbcball->SetParameter(j,double_crystalball->GetParameter(j));
        opt_dbcball->SetParName(j,double_crystalball->GetParName(j));
        dbpar[j] = double_crystalball->GetParameter(j);

      }

      opt_dbcball->FixParameter(1, hist_mean);
      opt_dbcball->FixParameter(6, hist_mean);

      double_crystalball->ReleaseParameter(1);
      double_crystalball->ReleaseParameter(6);

       cballfit = dpinvmasslm->Fit(double_crystalball,"0RMBQ");
       cballfit = dpinvmasslm->Fit(double_crystalball,"0RMBQ");
       cballfit = dpinvmasslm->Fit(double_crystalball,"0RMBQ");
       cballfit = dpinvmasslm->Fit(double_crystalball,"0RMBQ");
       cballfit = dpinvmasslm->Fit(double_crystalball,"0RMBQ");
       cballfit = dpinvmasslm->Fit(double_crystalball,"0RMBQ");
       cballfit = dpinvmasslm->Fit(double_crystalball,"0RMBQ");
       cballfit = dpinvmasslm->Fit(double_crystalball,"RMBQS+");



       //       cballfit->Print();

      // cout << "the limits of the signal window after checkin are " << binxl_val << " and " << binxh_val << endl;
      //   cout << " number of backgroun entries in the signal windows is " << entriesinint << endl;

      double dbfrac_1 = double_crystalball->GetParameter(0)/(double_crystalball->GetParameter(0) + double_crystalball->GetParameter(5));

      //cout << "frac1 " << dbfrac_1 << endl;

      double dbfrac_1_er = sqrt(pow(double_crystalball->GetParError(0)/(double_crystalball->GetParameter(0) + double_crystalball->GetParameter(5)) - (double_crystalball->GetParameter(0)*double_crystalball->GetParError(0))/pow(double_crystalball->GetParameter(0) + double_crystalball->GetParameter(5),2),2) +pow((-1*double_crystalball->GetParameter(0)*double_crystalball->GetParError(5))/(pow(double_crystalball->GetParameter(0) + double_crystalball->GetParameter(5),2)),2));

      //cout << "frac1 error " << dbfrac_1_er << endl;

      double dbfrac_2 = double_crystalball->GetParameter(5) /
                        (double_crystalball->GetParameter(0) +
                         double_crystalball->GetParameter(5));

      double dbfrac_2_er =
          sqrt(pow(double_crystalball->GetParError(5) /
                   (double_crystalball->GetParameter(0) +
                                   double_crystalball->GetParameter(5)
                               ) -
                       (double_crystalball->GetParameter(5) *
                        double_crystalball->GetParError(5)) /
                           pow(double_crystalball->GetParameter(0) +
                                   double_crystalball->GetParameter(5),
                               2),
                   2) +
               pow((-1 * double_crystalball->GetParameter(5) *
                    double_crystalball->GetParError(0)) /
                       (pow(double_crystalball->GetParameter(0) +
                                double_crystalball->GetParameter(5),
                            2)),
                   2));


      //      cout << "frac2 " << dbfrac_2 << endl;

      //    double dbfrac_2_er = double_crystalball->GetParError(5)/(double_crystalball->GetParameter(0) + double_crystalball->GetParameter(5)) - double_crystalball->GetParameter(5)*double_crystalball->GetParError(0)/(pow(double_crystalball->GetParameter(0) + double_crystalball->GetParameter(5),2));

      //  cout << "frac2 er " << dbfrac_2_er << endl;

      double dbw = sqrt(dbfrac_1*pow(double_crystalball->GetParameter(2),2) + dbfrac_2*pow(double_crystalball->GetParameter(7),2) );

      // cout << "weighted width " << dbw << endl;

      double dbw_er = (dbfrac_1_er*pow(double_crystalball->GetParameter(2),2) + dbfrac_2_er*pow(double_crystalball->GetParameter(7),2))/(2*dbw) + (double_crystalball->GetParameter(2)*double_crystalball->GetParError(2)*dbfrac_1 + double_crystalball->GetParameter(7)*double_crystalball->GetParError(7)*dbfrac_2)/dbw;

      // cout << "weighted width er " << dbw_er << endl;

      double tripSeff = zpgenideff / 100000;
      //      cout << " the efficiency from the gen id dist is " << tripSeff << endl;
      double tripSeffer = sqrt((tripSeff*(1-tripSeff))/100000);

      double fitfeff = (double_crystalball->Integral(hist_mean - 3*dbw, hist_mean + 3*dbw)/dpinvmasslm->GetBinWidth(0))/100000;
      double fitfeffer = ((double_crystalball->IntegralError(hist_mean-3*dbw,hist_mean+3*dbw,cballfit->GetParams(), cballfit->GetCovarianceMatrix().GetMatrixArray()))/dpinvmasslm->GetBinWidth(0))/100000;



    //Assuming that all possible combinations of final state muons are equivalent, the background count will be divided by 4, this will be corrected later by properly identifying the most energetic muons, this pair is the one coming out of the Z', this sould be checked in detail before filling the Z' invariant mass histograms//****************************************************************************************************************************

      TAxis *xaxis = dpinvmasslm->GetXaxis();
      TAxis *yaxis = dpinvmasslm->GetYaxis();
      Int_t binxl =  xaxis->FindBin(double_crystalball->GetParameter(1) - 30 * peakwidth);
      Double_t binxlerror = dpinvmasslm->GetBinError(binxl);
      Int_t binxh =  xaxis->FindBin(double_crystalball->GetParameter(1) + 30 * peakwidth);

      dpinvmasslm->GetXaxis()->SetRange(binxl, binxh);

      // cout << tripS << " " << tripSer << " " << B << " " << Ber << " " <<
      // intestep << endl;

            cout << hist_mean << " " << dbw << " " << dbw_er << " " << dbfrac_1 << " "
           << dbfrac_1_er << " " << dbfrac_2 << " " << dbfrac_2_er << " "
           << double_crystalball->GetParameter(2) << " "
           << double_crystalball->GetParError(2) << " "
           << double_crystalball->GetParameter(3) << " "
           << double_crystalball->GetParError(3) << " "
           << double_crystalball->GetParameter(4) << " "
           << double_crystalball->GetParError(4) << " "
           << double_crystalball->GetParameter(7) << " "
           << double_crystalball->GetParError(7) << " "
           << double_crystalball->GetParameter(8) << " "
           << double_crystalball->GetParError(8) << " "
           << double_crystalball->GetParameter(9) << " "
           << double_crystalball->GetParError(9) << endl;


      TString signalplotname = signalfilename + string(".eps");
        C1->Print(signalplotname);


}
