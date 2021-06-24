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

void opt_fit(double fitstep, TString fitstepname) {

  // Branching Ratio //
  TFile * br_fil = new TFile("../../merging_energies/Zp_BR.root");

  TGraph *gr_mu = new TGraph();
  gr_mu = (TGraph*)br_fil->Get("gr_mu");
  // **** //



  /// Fitted parameters functions ////
  TFile * db1fracfit = new TFile("../weighted_pionveto_on/redcball_paropt_3/dbfrac1.root");
  TFile * dbwwfit = new TFile("../weighted_pionveto_on/redcball_paropt_3/ww.root");
  TFile * db1wfit = new TFile("../weighted_pionveto_on/redcball_paropt_3/fw.root");
  TFile * db1alfit = new TFile("../weighted_pionveto_on/redcball_paropt_3/fal.root");
  TFile * db2alfit = new TFile("../weighted_pionveto_on/redcball_paropt_3/sal.root");
  TFile * db1nfit = new TFile("../weighted_pionveto_on/redcball_paropt_3/fn.root");
  TFile * db2nfit = new TFile("../weighted_pionveto_on/redcball_paropt_3/sn.root");
  TFile * deteffit = new TFile("../../dcball_optmize/kekcc_pion_par/newpionvetodeteff.root");

   TF1 * dbfrac1;
   TF1 * dbfrac2;
   TF1 * dbww;
   TF1 * db1w;
   TF1 * db1al;
   TF1 * db2al;
   TF1 * db1n;
   TF1 * db2n;
   TGraph * deteff;

   dbfrac1 = (TF1*)db1fracfit->Get("PrevFitTMP");
   db1w = (TF1*)db1wfit->Get("PrevFitTMP");
   db1al = (TF1*)db1alfit->Get("PrevFitTMP");
   db2al = (TF1*)db2alfit->Get("PrevFitTMP");
   db1n = (TF1*)db1nfit->Get("PrevFitTMP");
   db2n = (TF1*)db2nfit->Get("PrevFitTMP");
   dbww = (TF1*)dbwwfit->Get("PrevFitTMP");
   deteff = (TGraph*)deteffit->Get("gr_det_isr");
   /// ******************************** ////


  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  // gROOT->SetBatch(1);
  TFile *signal = new TFile("./merged_newpion.root");

  TH1F *dpinvmasslm;

   dpinvmasslm = (TH1F*)signal->Get(TString::Format("h_mycombitrigeff_3"));
   dpinvmasslm->Rebin(2);

  TCanvas *C1 = new TCanvas("C1", "", 100, 100, 1400, 1400);
     dpinvmasslm->SetLineColor(1);
     //   bginvmasslm->SetLineColor(2);
     //  bginvmasslm->Scale(-1.23171);

      C1->cd(1);
      dpinvmasslm->GetYaxis()->SetTitleOffset(1.6);
      dpinvmasslm->Draw();
    //    bginvmasslm->Draw("same");

    TF1 * pol3n;
    TF1 * crystalball;
    TF1 * double_crystalball;

  double hist_mean = fitstep;
  double peakwidth;
  double entriesatmean;
  double rms;
  double std_dev;
  double significance;
  double thirdpolchi;
  double dballnobschi;

    peakwidth = dbww->Eval(hist_mean);
   double dbw = dbww->Eval(hist_mean);
    rms = dpinvmasslm->GetBinWidth(1);
    std_dev = dpinvmasslm->GetStdDev(1);

    double db2w = sqrt((pow(db1w->Eval(hist_mean),2)*(abs(dbfrac1->Eval(hist_mean) - 1)))/ ( 1 - dbfrac1->Eval(hist_mean)));

    double lowerfit;
    double higherfit;

    if(hist_mean <= 0.212){
        lowerfit = hist_mean - 30*peakwidth;
        higherfit = hist_mean + 30*peakwidth;}
    else{
        lowerfit = hist_mean - 50*peakwidth;
        if(lowerfit <= 0.212){lowerfit = 0.212;}
        higherfit = hist_mean + 50*peakwidth;
        if(higherfit >= 10.){higherfit = 10.;}
        }

    // real pol3//
      pol3n = new TF1("norm pol3", "[0]*([1]+ [2]*x +[3]*x*x +[4]*x*x*x)", 0., 10.5);
    // pol6 pol3n = new TF1("norm pol3", "[0]*([1]+ [2]*x +[3]*x*x +[4]*x*x*x + [5]*x*x*x*x + [6]*x*x*x*x*x + [7]*x*x*x*x*x*x)", 0., 10.5);
  //  pol3n = new TF1("norm pol3", "[0]*([1]+ [2]*x +[3]*x*x +[4]*x*x*x + gaus(5) + gaus(8) + [11]*x*x*x*x + [12]*x*x*x*x*x)", 0., 10.5);

//	pol3n->SetParameter(6,hist_mean);
//	pol3n->SetParameter(9,hist_mean);
//	pol3n->SetParameter(12,0);

    pol3n->SetRange(lowerfit,higherfit);
      pol3n->SetParameter(0,1);

      TCanvas *C2 = new TCanvas("C2", "", 10, 10, 800, 800);
      C2->cd(1);

    pol3n->SetLineColor(3);
    pol3n->SetLineStyle(2);

  Int_t nbinsdp;
  Int_t nbinsbg;

  Double_t x;
  Double_t y;

  Double_t z;
  Double_t w;

  Double_t er;


  TH2F *sigres;
  TH2F *sigres_alt;
  TH2F *pull;

  TRandom *r1 = new TRandom();

  TH1D * h_pull[10000];
  TH1D * h_pull_res[2];
  h_pull_res[0]= new TH1D("pull distribution_0","pull;pull;entries;", 100,-5,5);
  h_pull_res[1]= new TH1D("pull distribution_1","pull;pull;entries;", 1000,-200,500);

  // here it is ok
      TFitResultPtr normpol3 = dpinvmasslm->Fit(pol3n,"0Q");
      normpol3 = dpinvmasslm->Fit(pol3n,"0Q");
      normpol3 = dpinvmasslm->Fit(pol3n,"0Q");
      normpol3 = dpinvmasslm->Fit(pol3n,"0Q");
      normpol3 = dpinvmasslm->Fit(pol3n,"0RMBQSW+");
      normpol3 = dpinvmasslm->Fit(pol3n,"RMBQSW+");


      thirdpolchi = normpol3->Chi2();
      double_t thirdpoldof = normpol3->Ndf();

      TAxis * xaxis = dpinvmasslm->GetXaxis();
      TAxis * yaxis = dpinvmasslm->GetYaxis();
      Int_t binxl = xaxis->FindBin(lowerfit);
      Double_t binxlerror = dpinvmasslm->GetBinError(binxl);
      Int_t binxh = xaxis->FindBin(higherfit);

    Int_t meanbin = xaxis->FindBin(hist_mean);

      Double_t binxl_val = dpinvmasslm->GetBinCenter(binxl);
      Double_t binxh_val = dpinvmasslm->GetBinCenter(binxh);
      if(binxl_val < 0.00000){
        binxl = xaxis->FindBin(0.0);
        binxl_val = dpinvmasslm->GetBinCenter(binxl);
      }
      if(binxh_val > 10.500){
        binxh = xaxis->FindBin(10.5);
        binxh_val = dpinvmasslm->GetBinCenter(binxh);
    }

    Int_t binint = binxh - binxl;
     Double_t sigwinbin = (binxh_val - binxl_val)/dpinvmasslm->GetBinWidth(1);
      double_t entriesinint = dpinvmasslm->Integral(binxl,binxh);

       TF1 * dbcrysnpol3 = new TF1("double crystal ball with a 3rd order poly", " [0]*(crystalball(1)  + crystalball(6)) + [11]*([12]+[13]*x +[14]*x*x + [15]*x*x*x)", 0., 10.5);



      dbcrysnpol3->SetParName(0,"Constant_1");
      dbcrysnpol3->SetParName(1,"Mean_1");
      dbcrysnpol3->SetParName(2,"Sigma_1");
      dbcrysnpol3->SetParName(3,"Alpha_1");
      dbcrysnpol3->SetParName(4,"N_1");
      dbcrysnpol3->SetParName(5,"Constant_2");
      dbcrysnpol3->SetParName(6,"Mean_2");
      dbcrysnpol3->SetParName(7,"Sigma_2");
      dbcrysnpol3->SetParName(8,"Alpha_2");
      dbcrysnpol3->SetParName(9,"N_2");

      dbcrysnpol3->SetParName(10,"pol3_norm");
      dbcrysnpol3->SetParName(11,"x_0");
      dbcrysnpol3->SetParName(12,"a");
      dbcrysnpol3->SetParName(13,"b");
      dbcrysnpol3->SetParName(14,"c");
      dbcrysnpol3->SetParName(15,"dbcrys_norm");

      dbcrysnpol3->SetNpx(1000);

      dbcrysnpol3->FixParameter(1,dbfrac1->Eval(hist_mean));
      dbcrysnpol3->FixParameter(2,hist_mean);
      dbcrysnpol3->FixParameter(3,db1w->Eval(hist_mean));
      dbcrysnpol3->FixParameter(4,db1al->Eval(hist_mean));
      dbcrysnpol3->FixParameter(5,db1n->Eval(hist_mean));
      dbcrysnpol3->FixParameter(6,(1. - dbfrac1->Eval(hist_mean)));
      dbcrysnpol3->FixParameter(7,hist_mean);
      dbcrysnpol3->FixParameter(8,db2w);
      dbcrysnpol3->FixParameter(9,db2al->Eval(hist_mean));
      dbcrysnpol3->FixParameter(10,db2n->Eval(hist_mean));
      dbcrysnpol3->SetParameter(11,pol3n->GetParameter(0));
      dbcrysnpol3->SetParameter(12,pol3n->GetParameter(1));
      dbcrysnpol3->SetParameter(13,pol3n->GetParameter(2));
      dbcrysnpol3->SetParameter(14,pol3n->GetParameter(3));
      dbcrysnpol3->SetParameter(15,pol3n->GetParameter(4));
      dbcrysnpol3->SetParameter(0,entriesinint);
      //  dbcrysnpol3->SetParameter(16,pol3n->GetParameter(5));
    //   dbcrysnpol3->FixParameter(17,pol3n->GetParameter(6));
    //   dbcrysnpol3->SetParameter(18,pol3n->GetParameter(7));
    //   dbcrysnpol3->SetParameter(19,pol3n->GetParameter(8));
    //   dbcrysnpol3->FixParameter(20,pol3n->GetParameter(9));
    //   dbcrysnpol3->SetParameter(21,pol3n->GetParameter(10));
    //   dbcrysnpol3->SetParameter(22,pol3n->GetParameter(11));
    //   dbcrysnpol3->SetParameter(23,pol3n->GetParameter(12));

     dbcrysnpol3->SetRange(lowerfit,higherfit);
     cout << "I came just before the double crystal ball fit " << endl;
     TFitResultPtr dballpolfit = dpinvmasslm->Fit(dbcrysnpol3,"Q0");
        dballpolfit = dpinvmasslm->Fit(dbcrysnpol3,"Q0");
        dballpolfit = dpinvmasslm->Fit(dbcrysnpol3,"Q0");
        dballpolfit = dpinvmasslm->Fit(dbcrysnpol3,"Q0");
        dballpolfit = dpinvmasslm->Fit(dbcrysnpol3,"Q0");
        dballpolfit = dpinvmasslm->Fit(dbcrysnpol3,"RMBQ0");
        dballpolfit = dpinvmasslm->Fit(dbcrysnpol3,"RMBQS+");

        dballnobschi = dballpolfit->Chi2();
        double_t dballnobsdof = dballpolfit->Ndf();

      ///####Toy Montecarlo##################/////


        // db +pol3//
        //  TF1 * dbcrysnpol3_forpull = new TF1("double crystal ball with a 3rd order poly", " [15]*(crystalball(0)  + crystalball(5)) + [10]*([11]+[12]*x +[13]*x*x + [14]*x*x*x)", 0., 10.5);
        // db + pol6 TF1 * dbcrysnpol3_forpull = new TF1("double crystal ball with a 3rd order poly", " [15]*(crystalball(0)  + crystalball(5)) + [10]*([11]+[12]*x +[13]*x*x + [14]*x*x*x + [16]*x*x*x*x + [17]*x*x*x*x*x + [18]*x*x*x*x*x*x)", 0., 10.5);
     // TF1 * dbcrysnpol3_forpull = new TF1("double crystal ball with a 3rd order poly", " [15]*(crystalball(0)  + crystalball(5)) + [10]*([11]+[12]*x +[13]*x*x + [14]*x*x*x + gaus(16) + gaus(19) + [22]*x*x*x*x + [23]*x*x*x*x*x )", 0., 10.5);

        TF1 * dbcrysnpol3_forpull = new TF1("3rd order poly for pull", "[0]*([1] + [2]*x + [3]*x*x + [4]*x*x*x)", 0., 10.5);


/*      dbcrysnpol3_forpull->SetParName(0,"Constant_1");
      dbcrysnpol3_forpull->SetParName(1,"Mean_1");
      dbcrysnpol3_forpull->SetParName(2,"Sigma_1");
      dbcrysnpol3_forpull->SetParName(3,"Alpha_1");
      dbcrysnpol3_forpull->SetParName(4,"N_1");
      dbcrysnpol3_forpull->SetParName(5,"Constant_2");
      dbcrysnpol3_forpull->SetParName(6,"Mean_2");
      dbcrysnpol3_forpull->SetParName(7,"Sigma_2");
      dbcrysnpol3_forpull->SetParName(8,"Alpha_2");
      dbcrysnpol3_forpull->SetParName(9,"N_2");

      dbcrysnpol3_forpull->SetParName(10,"pol3_norm");
      dbcrysnpol3_forpull->SetParName(11,"x_0");
      dbcrysnpol3_forpull->SetParName(12,"a");
      dbcrysnpol3_forpull->SetParName(13,"b");
      dbcrysnpol3_forpull->SetParName(14,"c");
      dbcrysnpol3_forpull->SetParName(15,"dbcrys_norm");
*/
      dbcrysnpol3_forpull->SetNpx(1000);

      dbcrysnpol3_forpull->SetParameter(0,dbcrysnpol3->GetParameter(10));
      dbcrysnpol3_forpull->SetParameter(1,dbcrysnpol3->GetParameter(12));
      dbcrysnpol3_forpull->SetParameter(2,dbcrysnpol3->GetParameter(13));
      dbcrysnpol3_forpull->SetParameter(3,dbcrysnpol3->GetParameter(14));
      dbcrysnpol3_forpull->SetParameter(4,dbcrysnpol3->GetParameter(15));
 //     dbcrysnpol3_forpull->SetParameter(5,dbcrysnpol3->GetParameter(5));
  /*
      dbcrysnpol3_forpull->FixParameter(6,dbcrysnpol3->GetParameter(6));
      dbcrysnpol3_forpull->FixParameter(7,dbcrysnpol3->GetParameter(7));
      dbcrysnpol3_forpull->FixParameter(8,dbcrysnpol3->GetParameter(8));
      dbcrysnpol3_forpull->FixParameter(9,dbcrysnpol3->GetParameter(9));
      dbcrysnpol3_forpull->FixParameter(10,dbcrysnpol3->GetParameter(10));
      dbcrysnpol3_forpull->FixParameter(11,dbcrysnpol3->GetParameter(11));
      dbcrysnpol3_forpull->FixParameter(12,dbcrysnpol3->GetParameter(12));
      dbcrysnpol3_forpull->FixParameter(13,dbcrysnpol3->GetParameter(13));
      dbcrysnpol3_forpull->FixParameter(14,dbcrysnpol3->GetParameter(14));
      dbcrysnpol3_forpull->SetParameter(16,dbcrysnpol3->GetParameter(16));
      dbcrysnpol3_forpull->SetParameter(17,dbcrysnpol3->GetParameter(17));
      dbcrysnpol3_forpull->SetParameter(18,dbcrysnpol3->GetParameter(18));
      dbcrysnpol3_forpull->SetParameter(19,dbcrysnpol3->GetParameter(19));
      dbcrysnpol3_forpull->SetParameter(20,dbcrysnpol3->GetParameter(20));
      dbcrysnpol3_forpull->SetParameter(21,dbcrysnpol3->GetParameter(21));
      dbcrysnpol3_forpull->SetParameter(22,dbcrysnpol3->GetParameter(22));
      dbcrysnpol3_forpull->SetParameter(23,dbcrysnpol3->GetParameter(23));
*/

      double low_range = lowerfit;
      double high_range = higherfit;

      dbcrysnpol3_forpull->SetRange(low_range,high_range);

      if(low_range < 0.0000){
        low_range = 0.00;
      }

      if(high_range > 10.500){
        high_range = 10.5;
      }

           TTimeStamp * c = new TTimeStamp();
           r1->SetSeed(c->GetNanoSec());
        
        


          for(int l = 0; l < 1000; l++){
           h_pull[l] = new TH1D("Pull distribution", "Toy MC reduced dimuon mass [GeV/c^{1}];m_{R};entries;", sigwinbin, low_range, high_range);
           h_pull[l]->Sumw2();
           c = new TTimeStamp();
           double_t timeseed = c->GetNanoSec();
           r1->SetSeed(timeseed);
           h_pull[l]->FillRandom("norm pol3",r1->Poisson(entriesinint));
      //     h_pull[l]->Rebin(1);
            TFitResultPtr dbnpol_forpull = h_pull[l]->Fit(dbcrysnpol3_forpull,"Q0");
            dbnpol_forpull = h_pull[l]->Fit(dbcrysnpol3_forpull,"RBQS+");
                    double signyield_alt = dbcrysnpol3_forpull->GetParameter(0);
           double signyield_alter = dbcrysnpol3_forpull->GetParError(0);
       //      cout << " the alternate signal yield is " << signyield_alt << " +/- " << signyield_alter << endl;
           h_pull_res[0]->Fill(signyield_alt/signyield_alter);
           }
         h_pull_res[0]->Fit("gaus", "BQS+");
        
///################################################################################################################################################################################////

        TF1 * dbball = new TF1("number of events with double crystal ball", "[10]*(crystalball(0) + crystalball(5))", lowerfit, higherfit);

        dbball->SetParameter(0,dbcrysnpol3->GetParameter(0));
        dbball->SetParameter(1,dbcrysnpol3->GetParameter(1));
        dbball->SetParameter(1,dbcrysnpol3->GetParameter(2));
        dbball->SetParameter(3,dbcrysnpol3->GetParameter(3));
        dbball->SetParameter(4,dbcrysnpol3->GetParameter(4));
        dbball->SetParameter(5,dbcrysnpol3->GetParameter(5));
        dbball->SetParameter(6,dbcrysnpol3->GetParameter(6));
        dbball->SetParameter(7,dbcrysnpol3->GetParameter(7));
        dbball->SetParameter(8,dbcrysnpol3->GetParameter(8));
        dbball->SetParameter(9,dbcrysnpol3->GetParameter(9));
        dbball->SetParameter(10, dbcrysnpol3->GetParameter(15));

        dbball->SetParError(0,dbcrysnpol3->GetParError(0));
        dbball->SetParError(1,dbcrysnpol3->GetParError(1));
        dbball->SetParError(2,dbcrysnpol3->GetParError(2));
        dbball->SetParError(3,dbcrysnpol3->GetParError(3));
        dbball->SetParError(4,dbcrysnpol3->GetParError(4));
        dbball->SetParError(5,dbcrysnpol3->GetParError(5));
        dbball->SetParError(6,dbcrysnpol3->GetParError(6));
        dbball->SetParError(7,dbcrysnpol3->GetParError(7));
        dbball->SetParError(8,dbcrysnpol3->GetParError(8));
        dbball->SetParError(9,dbcrysnpol3->GetParError(9));
        dbball->SetParError(10, dbcrysnpol3->GetParError(15));

        double sigtrigS = (dbball->Integral(hist_mean-3*dbw,hist_mean+3*dbw))/dpinvmasslm->GetBinWidth(1);
        double sigtrigSer = (dbball->GetParError(10)/dbball->GetParameter(10))*sigtrigS;

         double Nobs = sigtrigS;
         double Nobser = sigtrigSer;

          significance = ROOT::Math::Sign(Nobs)*sqrt(abs(dballnobschi - thirdpolchi));

          //          significance = ROOT::Math::Sign(Nobs)*sqrt(pow(dballnobschi/dballnobsdof - thirdpolchi/thirdpoldof,2));
    // cout << " the dball pol3 fit chi2 " << dballnobschi << " the pol3 by itself fit chi2 is << " << thirdpolchi << " the ratio is " << significance << endl;
   //   cout << " the number of observed events is " << ROOT::Math::Sign(Nobs) << endl;



   TCanvas * C8 = new TCanvas("nobspdf"," ",10,10,800,800);
   C8->cd(1);
   dbball->Draw("sames");

     TF1 * Nobspdff = new TF1("PDF for the number of observed events", " gaus ", Nobs - 20*Nobser, Nobs + 20*Nobser);
    Nobspdff->SetParameters(1.,Nobs,Nobser);
    Nobspdff->SetTitle("PDF for the number of observed events; Number of Observed Events ; Probability");
    Nobspdff->SetNpx(1000);
    Nobspdff->Draw();

    //   TString pdfplotname = signalfilename + string("_pdf.eps");
    TString pdfplotname = string("./dbonres") + fitstepname + string("_pdf.eps");
 C8->Print(pdfplotname);

   double ninetylim;
   double ninetylimer;
   double corninetylim = 0;
   double corninetylimer;
   double infinitylim;
   double intestep = 0.0001*Nobser;
   double intesteper = sqrt(intestep);
   double testep = 0.001*Nobser;
   double wholerange = 100*Nobser;
   double ninetyrange = 0.9*wholerange;

    ninetylim = Nobspdff->Integral(0.,0.00000001*Nobser);
   infinitylim = Nobspdff->Integral(0.,ninetyrange);
   corninetylim = ninetylim/(Nobspdff->Integral(0,Nobser));
   ninetylimer = sqrt(ninetylim);
   corninetylimer = ninetylimer/(Nobspdff->Integral(0,Nobser));


   while(corninetylim <= 0.90 ){
     ninetylim = Nobspdff->Integral(0,intestep);
     corninetylim = ninetylim/(Nobspdff->Integral(0,100*Nobser));
     intestep = intestep + 0.00001*Nobser ;
     testep = testep + 0.0001*Nobser;
     ninetylimer = sqrt(ninetylim);
     intesteper = sqrt(intestep);
     //     cout << " the 90% limit point is " << intestep <<  " ninety is " << ninetylim << " probability = " << corninetylim  << endl;
   }


   //cout << " the number of observed events from 0 to infinity is " << Nobspdff->Integral(0,100*Nobser) << endl;

        dpinvmasslm->GetXaxis()->SetRangeUser(lowerfit, higherfit);

    double entries = 0;
    double intnorm = 1/(sqrt(2*TMath::Pi())*dbw);
    sigres  = new TH2F("signal_fit_residual","res;X_mass[GeV/c^2];residuals;", 100,binxl_val,binxh_val, 100, -300, 300);
    sigres_alt = new TH2F("background_fit","res;X_mass[GeV/c^2];residuals;", 100,binxl_val-0.1,binxh_val+0.1, 100, -10, 800);
    pull = new TH2F("pull_fit","pull;X_mass[GeV/c^2];pull;",100,binxl_val,binxh_val,100,-300,300);

    for(int k = binxl-50; k < binxh+50; k++){
      x = dpinvmasslm->GetBinCenter(k);
      y = (dpinvmasslm->GetBinContent(k) - dbcrysnpol3->Eval(x));//realwidth;
     sigres->Fill(x,y);
     z = dpinvmasslm->GetBinCenter(k);
     //w = (bginvmasslm->GetBinContent(k) - gausnpol3->Eval(z));
     w = ((dpinvmasslm->GetBinContent(k) - dbcrysnpol3->Eval(z))*(dpinvmasslm->GetBinContent(k) - dbcrysnpol3->Eval(z)));
     sigres_alt->Fill(z,w);
     er = (dpinvmasslm->GetBinError(k)*dpinvmasslm->GetBinError(k));
     pull->Fill(z,w/er);
    }

     double sighistint = dpinvmasslm->Integral(binxl,binxh);
    Double_t histonorm = sighistint*dpinvmasslm->GetBinWidth(binxh);
    Double_t binxherror = dpinvmasslm->GetBinError(binxh);

    double bcount = (dpinvmasslm->IntegralAndError(binxl,binxh,binxlerror));
    double bcounter = sqrt(bcount);
    double bcountfit = pol3n->Integral(binxl_val,binxh_val);
    double bcountfiter = pol3n->IntegralError(binxl_val,binxh_val,normpol3->GetParams(),normpol3->GetCovarianceMatrix().GetMatrixArray());
    double scount = (dpinvmasslm->Integral(binxl,binxh));

    //double crystalball trial
    double_t tripSeff = deteff->Eval(hist_mean);

    double correctedmass = sqrt(pow(hist_mean,2)+4.*pow(0.1056583745,2) );

    //690.555 fb-1 is the Up4s lum
    //4.77836 fb-1 is the Up1s lum
    //3.5135 fb-1 is the Up2s lum
    //123.81655 fb-1 is the Up5s lum
    //85.73205 fb-1 is the continuum lum
    // all merged together is 925.28973 fb-1

   cout << hist_mean << " " << tripSeff  << " " << intestep << " " << ROOT::Math::Sign(Nobs)*intestep/(0.92528973*(gr_mu->Eval(hist_mean)*tripSeff)) << " " << significance << endl;


 //    cout << hist_mean << " " << correctedmass << " " << tripSeff  << " " << intestep << " " << intestep/(0.92528973*(gr_mu->Eval(hist_mean)*tripSeff)) << " " << significance << endl;

       TString signalfilename = "./dbonres" + fitstepname;

       TString signalplotname = signalfilename + string(".eps");
        C1->Print(signalplotname);

 TCanvas * C7 = new TCanvas("residuals"," ",10,10,800,800);
 //C7->Divide(2,1);

   C7->cd(1);
   sigres_alt->SetMarkerStyle(3);
   sigres_alt->Draw("*");
   //   C7->cd(2);
   // pull->SetMarkerStyle(3);
   // pull->Draw("*");

   //TString signalresname = signalfilename + string("res.eps");
   TString signalresname = "./dbonres" + fitstepname + string("res.eps");
 C7->Print(signalresname);

    TCanvas * C89 = new TCanvas("randomtest","",10,10,800,800);
 C89->Divide(2,1);
 C89->cd(1);
 h_pull[0]->Draw();

  C89->cd(2);
  h_pull_res[0]->Draw();

  // TString pullplotname = signalfilename + string("pull.eps");
  TString pullplotname = "./dbonres" + fitstepname + string("pull.eps");
 C89->Print(pullplotname);

}
