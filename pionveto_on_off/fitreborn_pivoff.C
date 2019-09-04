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

void fitreborn_pivoff(TString signalfilename) {

  TFile * br_fil = new TFile("../../merging_energies/Zp_BR.root");
  TFile * genid_isr_w = new TFile("./all_4mu_redmu_cor_isr_pvoff.root");

  TGraph *gr_mu = new TGraph();
  TGraph *gr_isr_w = new TGraph();

  gr_mu = (TGraph*)br_fil->Get("gr_mu");
  gr_isr_w = (TGraph*)genid_isr_w->Get("gr_w_isr");

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(1);
  // gROOT->SetBatch(1);
  TFile *signal = new TFile(signalfilename);
  TFile *bg = new TFile("../../pion_veto_off/aafhmumumumu_pioff.root");

  TH1F *dpinvmasslm;
  TH1F *bginvmasslm;

  TH1F * genid_invmass;

  // TF1 *parawidth = new TF1("parawidth","cheb9",0.5,10);
  TF1 *parawidth = new TF1("parawidth", "pol9",0.,10.);
  TF1 *parawei = new TF1("parawei", "pol5",0., 10.);

  TF1 * dbwidth = new TF1("dbwidth_evol", "pol3", 0,10);
  dbwidth->SetParameters(0.003697,0.001426,-0.0001015,-7.062e-06);

  dpinvmasslm = (TH1F*)signal->Get(TString::Format("h_mycombitrigeff_3"));
    bginvmasslm = (TH1F*)bg->Get(TString::Format("h_mycombitrigeff_3"));
    genid_invmass = (TH1F*)signal->Get("h_genidredmu_0");

    //Considering Ishikawa-san's correction on getting signal shape parameters based on truth tagged
    dpinvmasslm = genid_invmass;

  TCanvas *C1 = new TCanvas("C1", "", 100, 100, 1400, 1400);
     dpinvmasslm->SetLineColor(1);
     bginvmasslm->SetLineColor(2);
     bginvmasslm->Scale(0.23171);

      C1->cd(1);

      dpinvmasslm->GetYaxis()->SetTitleOffset(1.6);
    dpinvmasslm->Draw();
    bginvmasslm->Draw("same");

    TF1 * pol3n;
    TF1 * crystalball;
    TF1 * double_crystalball;

  double hist_mean;
  double peakwidth;
  double entriesatmean;
  double rms;
  double std_dev;
  double significance;
  double thirdpolchi;
  double dballnobschi;

  double zpgenideff = genid_invmass->GetEntries();

  hist_mean = dpinvmasslm->GetBinCenter(dpinvmasslm->GetMaximumBin());
  //  peakwidth = gr_isr_w->Eval(hist_mean);
  peakwidth = dbwidth->Eval(hist_mean);
  entriesatmean = dpinvmasslm->GetBinContent(dpinvmasslm->GetMaximumBin());
  rms = dpinvmasslm->GetBinWidth(1);
  std_dev = dpinvmasslm->GetStdDev(1);

    double lowerfit = hist_mean - 70*peakwidth;
    if(lowerfit < 0){lowerfit =0;}
    double higherfit = hist_mean + 30*peakwidth;
    if(higherfit > 10.5){higherfit = 10.5;}



    double_crystalball = new TF1("double_crystalball", "crystalball(0) + crystalball(5) ", 0. , 10.5);

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

    double_crystalball->SetParLimits(0,10.,entriesatmean);
    double_crystalball->SetParameter(1,hist_mean);
    double_crystalball->SetParameter(6,hist_mean);
    double_crystalball->SetParLimits(2,rms,3*peakwidth);
    double_crystalball->SetParLimits(3,-6.,0.);
    double_crystalball->SetParLimits(4,0.,4.);
    double_crystalball->SetParLimits(5,10.,entriesatmean);
    double_crystalball->SetParLimits(7,rms,3*peakwidth);
    double_crystalball->SetParLimits(8,0.0,6);
    // double_crystalball->SetParLimits(9,50.0,entriesatmean);
    double_crystalball->SetParLimits(9,0,4);

    double_crystalball->SetLineColor(4);
    double_crystalball->SetRange(lowerfit,  higherfit);
    //if(hist_mean >= 8.212){double_crystalball->SetRange(hist_mean - 1, hist_mean +1);}

      pol3n = new TF1("norm pol3", "[0]*([1]+ [2]*x +[3]*x*x +[4]*x*x*x)", 0., 10.5);
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
      TFitResultPtr normpol3 = bginvmasslm->Fit(pol3n,"RMEBQGIS+");
      thirdpolchi = normpol3->Chi2();

      // TFitResultPtr rebornfit = dpinvmasslm->Fit(triplegexp,"RBQS+");
       TFitResultPtr cballfit = dpinvmasslm->Fit(double_crystalball,"WWRMEGINBQS+");
       cballfit = dpinvmasslm->Fit(double_crystalball,"RMNQBS+");
       cballfit = dpinvmasslm->Fit(double_crystalball,"BSQ+");
       //       cballfit->Print();

      TAxis * xaxis = bginvmasslm->GetXaxis();
      TAxis * yaxis = bginvmasslm->GetYaxis();
      Int_t binxl = xaxis->FindBin(lowerfit);
      Double_t binxlerror = bginvmasslm->GetBinError(binxl);
      Int_t binxh = xaxis->FindBin(higherfit);
      Double_t binxl_val = dpinvmasslm->GetBinCenter(binxl);
      Double_t binxh_val = dpinvmasslm->GetBinCenter(binxh);
      if(binxl_val < 0.00000){
        binxl = xaxis->FindBin(0.);
        binxl_val = dpinvmasslm->GetBinCenter(binxl);
      }
      if(binxh_val > 10.500){
        binxh = xaxis->FindBin(10.5);
        binxh_val = dpinvmasslm->GetBinCenter(binxh);
    }
      Int_t binint = binxh - binxl;
     Double_t sigwinbin = (binxh_val - binxl_val)/bginvmasslm->GetBinWidth(0);
      //  cout << "the limits of the signal window are " << binxl_val << " and " << binxh_val << endl;
      // cout << "the limits of the signal window after checkin are " << binxl_val << " and " << binxh_val << endl;
      double_t entriesinint = bginvmasslm->Integral(binxl,binxh);
      //   cout << " number of backgroun entries in the signal windows is " << entriesinint << endl;

      double dbfrac_1 = double_crystalball->GetParameter(0)/(double_crystalball->GetParameter(0) + double_crystalball->GetParameter(5));

      double dbfrac_1_er = double_crystalball->GetParError(0)/(double_crystalball->GetParameter(0) + double_crystalball->GetParameter(5)) - double_crystalball->GetParameter(0)*double_crystalball->GetParError(5)/(pow(double_crystalball->GetParameter(0) + double_crystalball->GetParameter(5),2));

      double dbfrac_2 = double_crystalball->GetParameter(5)/(double_crystalball->GetParameter(0) + double_crystalball->GetParameter(5));

      double dbfrac_2_er = double_crystalball->GetParError(5)/(double_crystalball->GetParameter(0) + double_crystalball->GetParameter(5)) - double_crystalball->GetParameter(5)*double_crystalball->GetParError(0)/(pow(double_crystalball->GetParameter(0) + double_crystalball->GetParameter(5),2));

      double dbw = sqrt(dbfrac_1*pow(double_crystalball->GetParameter(2),2) + dbfrac_2*pow(double_crystalball->GetParameter(7),2) );

      double dbw_er = (dbfrac_1_er*pow(double_crystalball->GetParameter(2),2) + dbfrac_2_er*pow(double_crystalball->GetParameter(7),2))/(2*dbw) + (double_crystalball->GetParameter(2)*double_crystalball->GetParError(2)*dbfrac_1 + double_crystalball->GetParameter(7)*double_crystalball->GetParError(7)*dbfrac_2)/dbw;

       double tripS = double_crystalball->Integral(hist_mean-3*dbw, hist_mean+3*dbw);
      double tripSer = (double_crystalball->IntegralError(hist_mean-3*dbw, hist_mean+3*dbw, cballfit->GetParams(), cballfit->GetCovarianceMatrix().GetMatrixArray()));

      double tripSeff = zpgenideff/100000;
      //      cout << " the efficiency from the gen id dist is " << tripSeff << endl;
      double tripSeffer = sqrt((tripSeff*(1-tripSeff))/100000);

      double fitfeff = (double_crystalball->Integral(hist_mean - 3*dbw, hist_mean + 3*dbw)/dpinvmasslm->GetBinWidth(0))/100000;
      double fitfeffer = ((double_crystalball->IntegralError(hist_mean-3*dbw,hist_mean+3*dbw,cballfit->GetParams(), cballfit->GetCovarianceMatrix().GetMatrixArray()))/dpinvmasslm->GetBinWidth(0))/100000;


      TF1 * dbcrysnpol3 = new TF1("double crystal ball with a 3rd order poly", " [15]*(crystalball(0)  + crystalball(5)) + [10]*([11]+[12]*x +[13]*x*x + [14]*x*x*x)", 0., 10.5);

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

      dbcrysnpol3->FixParameter(0,double_crystalball->GetParameter(0));
      dbcrysnpol3->FixParameter(1,double_crystalball->GetParameter(1));
      dbcrysnpol3->FixParameter(2,double_crystalball->GetParameter(2));
      dbcrysnpol3->FixParameter(3,double_crystalball->GetParameter(3));
      dbcrysnpol3->FixParameter(4,double_crystalball->GetParameter(4));
      dbcrysnpol3->FixParameter(5,double_crystalball->GetParameter(5));
      dbcrysnpol3->FixParameter(6,double_crystalball->GetParameter(6));
      dbcrysnpol3->FixParameter(7,double_crystalball->GetParameter(7));
      dbcrysnpol3->FixParameter(8,double_crystalball->GetParameter(8));
      dbcrysnpol3->FixParameter(9,double_crystalball->GetParameter(9));
      dbcrysnpol3->SetParameter(10,pol3n->GetParameter(0));
      dbcrysnpol3->SetParameter(11,pol3n->GetParameter(1));
      dbcrysnpol3->SetParameter(12,pol3n->GetParameter(2));
      dbcrysnpol3->SetParameter(13,pol3n->GetParameter(3));
      dbcrysnpol3->SetParameter(14,pol3n->GetParameter(4));

      dbcrysnpol3->SetRange(lowerfit,higherfit);
        TFitResultPtr dballpolfit = bginvmasslm->Fit(dbcrysnpol3,"RMBQS+");

        dballnobschi = dballpolfit->Chi2();

      ///####Toy Montecarlo##################/////


      TF1 * dbcrysnpol3_forpull = new TF1("double crystal ball with a 3rd order poly", " [15]*(crystalball(0)  + crystalball(5)) + [10]*([11]+[12]*x +[13]*x*x + [14]*x*x*x)", 0., 10.5);


      dbcrysnpol3_forpull->SetParName(0,"Constant_1");
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

      dbcrysnpol3_forpull->SetNpx(1000);

      dbcrysnpol3_forpull->FixParameter(0,double_crystalball->GetParameter(0));
      dbcrysnpol3_forpull->FixParameter(1,double_crystalball->GetParameter(1));
      dbcrysnpol3_forpull->FixParameter(2,double_crystalball->GetParameter(2));
      dbcrysnpol3_forpull->FixParameter(3,double_crystalball->GetParameter(3));
      dbcrysnpol3_forpull->FixParameter(4,double_crystalball->GetParameter(4));
      dbcrysnpol3_forpull->FixParameter(5,double_crystalball->GetParameter(5));
      dbcrysnpol3_forpull->FixParameter(6,double_crystalball->GetParameter(6));
      dbcrysnpol3_forpull->FixParameter(7,double_crystalball->GetParameter(7));
      dbcrysnpol3_forpull->FixParameter(8,double_crystalball->GetParameter(8));
      dbcrysnpol3_forpull->FixParameter(9,double_crystalball->GetParameter(9));
      dbcrysnpol3_forpull->SetParameter(10,pol3n->GetParameter(0));
      dbcrysnpol3_forpull->SetParameter(11,pol3n->GetParameter(1));
      dbcrysnpol3_forpull->SetParameter(12,pol3n->GetParameter(2));
      dbcrysnpol3_forpull->SetParameter(13,pol3n->GetParameter(3));
      dbcrysnpol3_forpull->SetParameter(14,pol3n->GetParameter(4));

      double low_range = lowerfit;
      double high_range = higherfit;

      dbcrysnpol3_forpull->SetRange(low_range,high_range);

      if(low_range < 0.0000){
        low_range = 0.00;
      }

      if(high_range > 10.500){
        high_range = 10.5;
      }

            for(int l = 0; l < 1000; l++){
           h_pull[l] = new TH1D("Pull distribution", "Toy MC reduced dimuon mass [GeV/c^{2}];m_{R};entries;", sigwinbin, low_range, high_range);
           h_pull[l]->Sumw2();
           TTimeStamp * c = new TTimeStamp();
           TTimeStamp * d = new TTimeStamp();
           double_t timeseed = c->GetNanoSec();
           double_t timeseed2 = d->GetNanoSec();
           r1->SetSeed(timeseed);
           h_pull[l]->FillRandom("norm pol3",r1->Poisson(entriesinint));
            TFitResultPtr dbnpol_forpull = h_pull[l]->Fit(dbcrysnpol3_forpull,"RNMBGQS+");
            dbnpol_forpull = h_pull[l]->Fit(dbcrysnpol3_forpull,"RMBGQS+");
           double signyield_alt = dbcrysnpol3_forpull->GetParameter(15);
           double signyield_alter = dbcrysnpol3_forpull->GetParError(15);
           //  cout << " the alternate signal yield is " << signyield_alt << " +/- " << signyield_alter << endl;
           h_pull_res[0]->Fill(signyield_alt/signyield_alter);
           }
         ////############################################////
            TFitResultPtr pull_result =  h_pull_res[0]->Fit("gaus","EWWMNISGQ+");
            pull_result =  h_pull_res[0]->Fit("gaus", "");

        TF1 * dbball = new TF1("number of events with double crystal ball", "[10]*(crystalball(0) + crystalball(5))", hist_mean-50*peakwidth, hist_mean+70*peakwidth);

        dbball->SetParameter(0,double_crystalball->GetParameter(0));
        dbball->SetParameter(1,double_crystalball->GetParameter(1));
        dbball->SetParameter(2,double_crystalball->GetParameter(2));
        dbball->SetParameter(3,double_crystalball->GetParameter(3));
        dbball->SetParameter(4,double_crystalball->GetParameter(4));
        dbball->SetParameter(5,double_crystalball->GetParameter(5));
        dbball->SetParameter(6,double_crystalball->GetParameter(6));
        dbball->SetParameter(7,double_crystalball->GetParameter(7));
        dbball->SetParameter(8,double_crystalball->GetParameter(8));
        dbball->SetParameter(9,double_crystalball->GetParameter(9));
        dbball->SetParameter(10, dbcrysnpol3->GetParameter(15));

        dbball->SetParError(0,double_crystalball->GetParError(0));
        dbball->SetParError(1,double_crystalball->GetParError(1));
        dbball->SetParError(2,double_crystalball->GetParError(2));
        dbball->SetParError(3,double_crystalball->GetParError(3));
        dbball->SetParError(4,double_crystalball->GetParError(4));
        dbball->SetParError(5,double_crystalball->GetParError(5));
        dbball->SetParError(6,double_crystalball->GetParError(6));
        dbball->SetParError(7,double_crystalball->GetParError(7));
        dbball->SetParError(8,double_crystalball->GetParError(8));
        dbball->SetParError(9,double_crystalball->GetParError(9));
        dbball->SetParError(10, dbcrysnpol3->GetParError(15));

        double sigtrigS = (dbball->Integral(hist_mean-3*dbw,hist_mean+3*dbw))/bginvmasslm->GetBinWidth(0);
        double sigtrigSer = (dbball->GetParError(10)/dbball->GetParameter(10))*sigtrigS;

         double Nobs = sigtrigS;
          double Nobser = sigtrigSer;

    significance = ROOT::Math::Sign(Nobs)*sqrt(dballnobschi/thirdpolchi);
   //   cout << " the number of observed events is " << ROOT::Math::Sign(Nobs) << endl;

   dbball->Draw("sames");


   TCanvas * C8 = new TCanvas("nobspdf"," ",10,10,800,800);
   C8->cd(1);


     TF1 * Nobspdff = new TF1("PDF for the number of observed events", " gaus ", Nobs - 20*Nobser, Nobs + 20*Nobser);
    Nobspdff->SetParameters(1.,Nobs,Nobser);
    Nobspdff->SetTitle("PDF for the number of observed events; Number of Observed Events ; Probability");
    Nobspdff->SetNpx(1000);
    Nobspdff->Draw();

   TString pdfplotname = signalfilename + string("_pdf.eps");
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
   double B = (pol3n->Integral(hist_mean-3*dbw,hist_mean+3*dbw));
   double Ber= (pol3n->IntegralError(double_crystalball->GetParameter(1)-3*dbw,double_crystalball->GetParameter(1)+3*dbw,normpol3->GetParams(), normpol3->GetCovarianceMatrix().GetMatrixArray()));

        dpinvmasslm->GetXaxis()->SetRangeUser(lowerfit, higherfit);

    double entries = 0;
    double intnorm = 1/(sqrt(2*TMath::Pi())*dbw);
    sigres  = new TH2F("signal_fit_residual","res;X_mass[GeV/c^2];residuals;", 100,binxl_val,binxh_val, 100, -100, 100);
    sigres_alt = new TH2F("background_fit","res;X_mass[GeV/c^2];residuals;", 100,binxl_val-0.1,binxh_val+0.1, 100, -10, 100);
    pull = new TH2F("pull_fit","pull;X_mass[GeV/c^2];pull;",100,binxl_val,binxh_val,100,-300,300);

    for(int k = binxl-50; k < binxh+50; k++){
      x = dpinvmasslm->GetBinCenter(k);
      y = (dpinvmasslm->GetBinContent(k) - double_crystalball->Eval(x));//realwidth;
     sigres->Fill(x,y);
     z = bginvmasslm->GetBinCenter(k);
     //w = (bginvmasslm->GetBinContent(k) - gausnpol3->Eval(z));
     w = ((bginvmasslm->GetBinContent(k) - dbcrysnpol3->Eval(z))*(bginvmasslm->GetBinContent(k) - dbcrysnpol3->Eval(z)));
     sigres_alt->Fill(z,w);
     er = (bginvmasslm->GetBinError(k)*bginvmasslm->GetBinError(k));
     pull->Fill(z,w/er);
    }

     double sighistint = dpinvmasslm->Integral(binxl,binxh);
    Double_t histonorm = sighistint*dpinvmasslm->GetBinWidth(binxh);
    Double_t binxherror = bginvmasslm->GetBinError(binxh);
    //Assuming that all possible combinations of final state muons are equivalent, the background count will be divided by 4, this will be corrected later by properly identifying the most energetic muons, this pair is the one coming out of the Z', this sould be checked in detail before filling the Z' invariant mass histograms//****************************************************************************************************************************
    double bcount = (bginvmasslm->IntegralAndError(binxl,binxh,binxlerror));
    double bcounter = sqrt(bcount);
    double bcountfit = pol3n->Integral(binxl_val,binxh_val);
    double bcountfiter = pol3n->IntegralError(binxl_val,binxh_val,normpol3->GetParams(),normpol3->GetCovarianceMatrix().GetMatrixArray());
    double scount = (dpinvmasslm->Integral(binxl,binxh));


    //double crystalball trial

    //        cout << hist_mean << " " << dbw << " " << dbw_er << " " << dbfrac_1 << " " << dbfrac_1_er << " " << dbfrac_2 << " " << dbfrac_2_er << " " << fitfeff << " " << fitfeffer << " " << intestep << " " << intestep/(0.690555*(gr_mu->Eval(hist_mean)*tripSeff)) << " " << significance << endl;

     cout << hist_mean << " " << dbw << " " << dbw_er << " " << dbfrac_1 << " " << dbfrac_1_er << " " << dbfrac_2 << " " << dbfrac_2_er << " " << tripSeff << " " << tripSeffer << " " << intestep << " " << intestep/(0.690555*(gr_mu->Eval(hist_mean)*tripSeff)) << " " << significance << " " << double_crystalball->GetParameter(2) << " " << double_crystalball->GetParError(2) << " " << double_crystalball->GetParameter(3) << " " << double_crystalball->GetParError(3) << " " << double_crystalball->GetParameter(4) << " " << double_crystalball->GetParError(4) << " " << double_crystalball->GetParameter(7) << " " << double_crystalball->GetParError(7) << " " << double_crystalball->GetParameter(8) << " " << double_crystalball->GetParError(8) << " " << double_crystalball->GetParameter(9) << " " << double_crystalball->GetParError(9) << endl;


    //    cout << hist_mean << " " << tripSeff << " " << tripSeffer << " " << fitfeff << " " << fitfeffer << endl;

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

TString signalresname = signalfilename + string("res.eps");
 C7->Print(signalresname);

    TCanvas * C89 = new TCanvas("randomtest","",10,10,800,800);
 C89->Divide(2,1);
 C89->cd(1);
 h_pull[0]->Draw();

  C89->cd(2);
  h_pull_res[0]->Draw();

  TString pullplotname = signalfilename + string("pull.eps");
 C89->Print(pullplotname);

}
