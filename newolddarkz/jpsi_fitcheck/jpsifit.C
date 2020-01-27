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

void jpsifit(TString signalfilename) {

  TFile * br_fil = new TFile("../../../merging_energies/Zp_BR.root");
  TFile * genid_w = new TFile("../../../signal_mc_dist_reader/newolddarkz_detandw.root");

  TGraph *gr_mu = new TGraph();
  TGraph *gr_w = new TGraph();

  gr_mu = (TGraph*)br_fil->Get("gr_mu");
  gr_w = (TGraph*)genid_w->Get("gr_w");

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  // gROOT->SetBatch(1);
  TFile *signal = new TFile(signalfilename);
  TFile *bg = new TFile("../../../new_pion_veto/smearbackvphojpsimumu.root");

  TH1F *dpinvmasslm;
  TH1F *bginvmasslm;

  TH1F * genid_redmass;

  TF1 *parawidth = new TF1("parawidth", "pol9",0.,10.);
  TF1 *parawei = new TF1("parawei", "pol5",0., 10.);

  genid_redmass = (TH1F*)signal->Get(TString::Format("h_mycombitrigeffw_0"));
    dpinvmasslm = (TH1F*)signal->Get(TString::Format("h_mycombitrigeffw_0"));
    bginvmasslm = (TH1F*)bg->Get(TString::Format("h_mycombitrigeffw_3"));

  TCanvas *C1 = new TCanvas("C1", "", 10, 10, 800, 800);
     dpinvmasslm->SetLineColor(1);
     bginvmasslm->SetLineColor(2);
     bginvmasslm->Scale(0.23171);
     C1->cd(1);

     dpinvmasslm->GetYaxis()->SetTitleOffset(1.6);
    dpinvmasslm->Draw();
    bginvmasslm->Draw("same");

  TF1 * pol3n;
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

  hist_mean = dpinvmasslm->GetBinCenter(dpinvmasslm->GetMaximumBin());
  peakwidth = gr_w->Eval(hist_mean);
  entriesatmean = dpinvmasslm->GetBinContent(dpinvmasslm->GetMaximumBin());
  rms = dpinvmasslm->GetBinWidth(1);
  std_dev = dpinvmasslm->GetStdDev(1);

  double lowerfit = hist_mean -60*peakwidth;
  if(lowerfit < 0){lowerfit =0;}
  double higherfit = hist_mean + 70*peakwidth;
  if(higherfit > 10.0){higherfit = 10.;}

      singpartripgaus = new TF1("single_par_triple_gaus","(([0]/((sqrt(2*TMath::Pi()*[2]*[2]))))*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) +(([3])/((sqrt(2*TMath::Pi()*[4]*[4]))))*exp(-((x-[1])*(x-[1]))/(2*[4]*[4])) + (([6])/((sqrt(2*TMath::Pi()*[5]*[5]))))*exp(-((x-[1])*(x-[1]))/(2*[5]*[5])) )", 0., 10.5);
    singpartripgaus->SetLineColor(1);
    // defining parameters names //

    singpartripgaus->SetParName(0,"gaus_height1");
    singpartripgaus->SetParName(1,"mean");
    singpartripgaus->SetParName(2,"first_width");
    singpartripgaus->SetParName(3,"gaus_height2");
    singpartripgaus->SetParName(4,"second_width");
     singpartripgaus->SetParName(5,"third_width");
     singpartripgaus->SetParName(6,"gaus_height3");
     singpartripgaus->SetNpx(1000);

    singpartripgaus->SetParLimits(0,0.0,entriesatmean);
    singpartripgaus->SetParLimits(1,hist_mean-3*peakwidth,hist_mean+3*peakwidth);
    singpartripgaus->SetParLimits(2,rms,20*peakwidth);
    singpartripgaus->SetParLimits(3,0.0,entriesatmean);
    singpartripgaus->SetParLimits(4,rms,20*peakwidth);
     singpartripgaus->SetParLimits(5,rms,20*peakwidth);
     singpartripgaus->SetParLimits(6,0.0,entriesatmean);

    singpartripgaus->SetRange(hist_mean-20*peakwidth,hist_mean+20*peakwidth);

      pol3n = new TF1("norm pol3", "[0]*([1]+ [2]*x +[3]*x*x +[4]*x*x*x)", 0., 10.5);
      pol3n->SetRange(hist_mean-20*peakwidth,hist_mean+20*peakwidth);
      pol3n->FixParameter(0,1);

       pol3n->SetLineColor(3);
    pol3n->SetLineStyle(2);

  TH1F *fom = new TH1F("figure of merit with B only","fom;# of loose muons;",4,0,3);
  TH1F *foma = new TH1F( "figure of merit with S+B","foma;# of loose muons;",4,0,3);
  TH1F *ponzifom = new TH1F("figure of merit following ponzi model","ponzi fom; # of loose muons;",4,0,3);

  Int_t nbinsdp;
  Int_t nbinsbg;

  Double_t x;
  Double_t y;

  Double_t z;
  Double_t w;

  Double_t er;

  double gaus_pol3_ParEr[7];
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

   hist_mean = dpinvmasslm->GetBinCenter(dpinvmasslm->GetMaximumBin());
   std_dev = dpinvmasslm->GetStdDev();
   peakwidth = gr_w->Eval(hist_mean);
   TFitResultPtr normpol3 = bginvmasslm->Fit(pol3n,"RBQS+");
      thirdpolchi = normpol3->Chi2();


      singpartripgaus->SetNpx(1000);
      TFitResultPtr finalfit = dpinvmasslm->Fit(singpartripgaus,"RBQS+");
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


      TAxis * xaxis = bginvmasslm->GetXaxis();
      TAxis * yaxis = bginvmasslm->GetYaxis();
      Int_t binxl = xaxis->FindBin(singpartripgaus->GetParameter(1)-20*peakwidth);
      Double_t binxlerror = bginvmasslm->GetBinError(binxl);
      Int_t binxh = xaxis->FindBin( singpartripgaus->GetParameter(1)+20*peakwidth);
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
     Double_t sigwinbin = (binxh_val - binxl_val)/bginvmasslm->GetBinWidth(0);
      //  cout << "the limits of the signal window are " << binxl_val << " and " << binxh_val << endl;
      // cout << "the limits of the signal window after checkin are " << binxl_val << " and " << binxh_val << endl;
      double_t entriesinint = bginvmasslm->Integral(binxl,binxh);
      //   cout << " number of backgroun entries in the signal windows is " << entriesinint << endl;

      double tripger[] = {triperwh,peakwidth,triperww};

       double tripS = singpartripgaus->Integral(hist_mean-3*peakwidth, hist_mean+3*peakwidth);
      double tripSer = (singpartripgaus->IntegralError(hist_mean-3*peakwidth, hist_mean+3*peakwidth, finalfit->GetParams(), finalfit->GetCovarianceMatrix().GetMatrixArray()));
      //      double tripSeff = (tripS/(dpinvmasslm[i]->GetBinWidth(0)))/100000;
      //   double tripSeffer = (tripSer/(dpinvmasslm[i]->GetBinWidth(0)))/100000;

      // CORRECTED EFFICIENCY WITH ISR FACTOR 0.848148 +/- 0.094711635
      // double tripSeff = (dpinvmasslm[i]->Integral(binxl ,binxh))/100000;
      double tripSeff = 0.848148*(dpinvmasslm->GetEntries())/100000;
      //cout << " the new efficiency is " << tripSeff << " between " << hist_mean - 3*tripww << " and " << hist_mean + 3*tripww << endl;
      double tripSeffer = sqrt((tripSeff*(1-tripSeff))/100000);

           TF1 * tripgtripo = new TF1("tripgtripo", " [12]*((([0]))*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) +(([3]))*exp(-((x-[1])*(x-[1]))/(2*[4]*[4])) + (([6]))*exp(-((x-[1])*(x-[1]))/(2*[5]*[5]))) + [7]*([8]+[9]*x +[10]*x*x + [11]*x*x*x)", 0., 10.5);

           tripgtripo->SetNpx(1000);

        tripgtripo->FixParameter(0,normtripg[0]);
        tripgtripo->FixParameter(1,singpartripgaus->GetParameter(1));
        tripgtripo->SetParError(1,singpartripgaus->GetParError(1));
         tripgtripo->FixParameter(2,singpartripgaus->GetParameter(2));
         //  tripgtripo->SetParameter(2,singpartripgaus->GetParameter(2));
        tripgtripo->SetParError(2,triperww);
      tripgtripo->FixParameter(3,normtripg[1]);
      tripgtripo->SetParError(3,singpartripgaus->GetParError(3));
      tripgtripo->FixParameter(4,singpartripgaus->GetParameter(4));
      //tripgtripo->SetParameter(4,singpartripgaus->GetParameter(4));
          tripgtripo->FixParameter(5,singpartripgaus->GetParameter(5));
      //  tripgtripo->SetParameter(5,singpartripgaus->GetParameter(5));
          tripgtripo->FixParameter(6,1. -normtripg[0] -normtripg[1]);
      tripgtripo->SetParameter(7,pol3n->GetParameter(0));

                tripgtripo->SetParameter(8,pol3n->GetParameter(1));
                tripgtripo->SetParameter(9,pol3n->GetParameter(2));
                tripgtripo->SetParError(9,pol3n->GetParError(2));
                tripgtripo->SetParameter(10,pol3n->GetParameter(3));
                tripgtripo->SetParError(10,pol3n->GetParError(3));
                tripgtripo->SetParameter(11,pol3n->GetParameter(4));
                tripgtripo->SetParError(11,pol3n->GetParError(4));

                tripgtripo->SetRange(hist_mean-20*peakwidth,hist_mean+20*peakwidth);
                TFitResultPtr trigtrip = bginvmasslm->Fit(tripgtripo,"RBQS+");
                //trigtrip->Print();

                double tripochi = trigtrip->Chi2();

                //    cout << " the chi square for the 3rd order poly fit over the background is " << thirdpolchi << " the chi square for the gaus 3rd poly over background is " << gausnobschi << " and the significance is = " << significance << endl;
      ///####Toy Montecarlo##################/////

                /*TF1 * gausnpol3_forpull = new TF1("norm double gaus with pol3 forpull", " gaus(0)  + [3]*([4]+[5]*x +[6]*x*x + [7]*x*x*x)", binxl_val, binxh_val);

      gausnpol3_forpull->SetLineColor(6);
      gausnpol3_forpull->SetParName(0, "height1");
      gausnpol3_forpull->SetParName(1, "mean1");
      gausnpol3_forpull->SetParName(2, "width1");
      gausnpol3_forpull->SetParName(3,"norm");
      gausnpol3_forpull->SetParName(4, "p[0]");
      gausnpol3_forpull->SetParName(5, "p[1]");
      gausnpol3_forpull->SetParName(6, "p[2]");
      gausnpol3_forpull->SetParName(7, "p[3]");

      gausnpol3_forpull->SetParameter(0,tripwh);
      gausnpol3_forpull->FixParameter(1,hist_mean);
      gausnpol3_forpull->SetParError(1,triple_gaus->GetParError(1));
      gausnpol3_forpull->FixParameter(2,tripww);
      gausnpol3_forpull->SetParError(2,triperww);
      gausnpol3_forpull->FixParameter(4,pol3n->GetParameter(1));
      gausnpol3_forpull->SetParError(4,pol3n->GetParError(1));
      gausnpol3_forpull->FixParameter(5,pol3n->GetParameter(2));
      gausnpol3_forpull->SetParError(5,pol3n->GetParError(2));
      gausnpol3_forpull->FixParameter(6,pol3n->GetParameter(3));
      gausnpol3_forpull->SetParError(6,pol3n->GetParError(3));
      gausnpol3_forpull->FixParameter(7,pol3n->GetParameter(4));
      gausnpol3_forpull->SetParError(7,pol3n->GetParError(4));

              gausnpol3_forpull->SetParameter(4,pol3n->GetParameter(1));
      gausnpol3_forpull->SetParError(4,pol3n->GetParError(1));
      gausnpol3_forpull->SetParameter(5,pol3n->GetParameter(2));
      gausnpol3_forpull->SetParError(5,pol3n->GetParError(2));
      gausnpol3_forpull->SetParameter(6,pol3n->GetParameter(3));
      gausnpol3_forpull->SetParError(6,pol3n->GetParError(3));
      gausnpol3_forpull->SetParameter(7,pol3n->GetParameter(4));
      gausnpol3_forpull->SetParError(7,pol3n->GetParError(4));


      double low_range = hist_mean - 20*peakwidth;
      double high_range = hist_mean + 20*peakwidth;

      if(low_range < 0.0000){
        low_range = 0.00;
      }

      if(high_range > 10.500){
        high_range = 10.5;
      }

            for(int l = 0; l < 100; l++){
           h_pull[l] = new TH1D("Pull distribution", "Toy MC reduced dimuon mass [GeV/c^{2}];m_{R};entries;", sigwinbin, low_range, high_range);
           h_pull[l]->Sumw2();
           TTimeStamp * c = new TTimeStamp();
           TTimeStamp * d = new TTimeStamp();
           double_t timeseed = c->GetNanoSec();
           double_t timeseed2 = d->GetNanoSec();
           r1->SetSeed(timeseed);
           h_pull[l]->FillRandom("norm pol3",r1->Poisson(entriesinint));
            TFitResultPtr gausnpol_forpull = h_pull[l]->Fit(gausnpol3_forpull,"RBQS+");
           double signyield_alt = gausnpol3_forpull->GetParameter(0);
           double signyield_alter = gausnpol3_forpull->GetParError(0);
           //  cout << " the alternate signal yield is " << signyield_alt << " +/- " << signyield_alter << endl;
           h_pull_res[0]->Fill(signyield_alt/signyield_alter);
           }
         ////############################################////
         h_pull_res[0]->Fit("gaus", "BQS+");

                */
                TF1 * tripsiggaus = new TF1("number of events with trip gaus", "[9]*(gaus(0) + gaus(3) + gaus(6))", hist_mean-40*peakwidth, hist_mean+40*peakwidth );

                //area
                tripsiggaus->SetParameter(9,tripgtripo->GetParameter(12));
                tripsiggaus->SetParError(9, tripgtripo->GetParError(12));

                // single mean
                tripsiggaus->SetParameter(1,tripgtripo->GetParameter(1));
                tripsiggaus->SetParameter(4,tripgtripo->GetParameter(1));
                tripsiggaus->SetParameter(7,tripgtripo->GetParameter(1));

                //different heights
                tripsiggaus->SetParameter(0,tripgtripo->GetParameter(0));
                tripsiggaus->SetParError(0,tripgtripo->GetParError(0));
                tripsiggaus->SetParameter(3,tripgtripo->GetParameter(3));
                tripsiggaus->SetParError(3,tripgtripo->GetParError(3));
                tripsiggaus->SetParameter(6,tripgtripo->GetParameter(6));
                tripsiggaus->SetParError(6,tripgtripo->GetParError(6));

                //different widths
                tripsiggaus->SetParameter(2,tripgtripo->GetParameter(2));
                tripsiggaus->SetParError(2,tripgtripo->GetParError(2));
                tripsiggaus->SetParameter(4,tripgtripo->GetParameter(4));
                tripsiggaus->SetParError(4,tripgtripo->GetParError(4));
                tripsiggaus->SetParameter(8,tripgtripo->GetParameter(5));
                tripsiggaus->SetParError(8,tripgtripo->GetParError(5));


                double sigtrigS = tripsiggaus->Integral(hist_mean-3*tripsiggaus->GetParameter(4),hist_mean+3*tripsiggaus->GetParameter(4))/bginvmasslm->GetBinWidth(0);
                 double sigtrigSer = (tripsiggaus->GetParError(9)/tripsiggaus->GetParameter(9))*sigtrigS;

                 // cout << " the number of observed events for the trip = " << sigtrigS << " +/- " << sigtrigSer << endl;


   double Nobs = sigtrigS;
   double Nobser = sigtrigSer;

   significance = ROOT::Math::Sign(Nobs)*sqrt(abs(tripochi - thirdpolchi));
   //   cout << " the number of observed events is " << ROOT::Math::Sign(Nobs) << endl;

   // TCanvas *C2 = new TCanvas("C2", "", 10, 10, 800, 800);
   // C2->cd(1);

   //   pol3n->Draw();
   // siggaus->Draw("same");

   TCanvas * C8 = new TCanvas("nobspdf"," ",10,10,800,800);
   C8->cd(1);

   TF1 * Nobspdff = new TF1("PDF for the number of observed events", " gaus(0) ", Nobs - 20*Nobser, Nobs + 20*Nobser);
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
      // cout << " the 90% limit point is " << intestep <<  " ninety is " << ninetylim << " probability = " << corninetylim  << endl;
   }


   //cout << " the number of observed events from 0 to infinity is " << Nobspdff->Integral(0,100*Nobser) << endl;
    double B = (pol3n->Integral(singpartripgaus->GetParameter(1)-3*tripww,singpartripgaus->GetParameter(1)+3*tripww));
    double Ber= (pol3n->IntegralError(singpartripgaus->GetParameter(1)-3*tripww,singpartripgaus->GetParameter(1)+3*tripww,normpol3->GetParams(), normpol3->GetCovarianceMatrix().GetMatrixArray()));

    dpinvmasslm->GetXaxis()->SetRange(binxl,binxh);
    double entries = 0;
    double intnorm = 1/(sqrt(2*TMath::Pi())*tripww);
    sigres  = new TH2F("signal_fit_residual","res;X_mass[GeV/c^2];residuals;", 100,binxl_val,binxh_val, 100, -300, 300);
    sigres_alt = new TH2F("background_fit","res;X_mass[GeV/c^2];residuals;", 100,binxl_val,binxh_val, 100, -300, 300);
    pull = new TH2F("pull_fit","pull;X_mass[GeV/c^2];pull;",100,binxl_val,binxh_val,100,-300,300);
    for(int k = binxl-50; k < binxh+50; k++){
      x = dpinvmasslm->GetBinCenter(k);
      y = (dpinvmasslm->GetBinContent(k) - singpartripgaus->Eval(x));//realwidth;
     sigres->Fill(x,y);
     z = bginvmasslm->GetBinCenter(k);
     w = (bginvmasslm->GetBinContent(k) - singpartripgaus->Eval(z));
     sigres_alt->Fill(z,w);
     er = bginvmasslm->GetBinError(k);
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


    //model dependent with the cross section
    // cout << tripSeff << " " << tripSeffer << " " << gr_mu->Eval(singpartripgaus->GetParameter(1)) << " " << 0.005 << " " <<  bcountfit/(dpinvmasslm->GetBinWidth(0)) << " " << bcountfiter/(dpinvmasslm->GetBinWidth(0)) << " " << intestep << " " << (intestep)/(0.977*(gr_mu->Eval(singpartripgaus->GetParameter(1)))*tripSeff) << " " << intesteper/(0.977*(gr_mu->Eval(singpartripgaus->GetParameter(1)))*tripSeff) << " " << tripww << " " << triperww << " " << significance << endl;

    //naming triple gaussian parameters

    double firstwidth = singpartripgaus->GetParameter(2);
    double firstwidther = singpartripgaus->GetParError(2);

    double secondwidth = singpartripgaus->GetParameter(4);
    double secondwidther = singpartripgaus->GetParError(4);

    double thirdwidth = singpartripgaus->GetParameter(5);
    double thirdwidther = singpartripgaus->GetParError(5);

       //to parametrize the double gaussian completely
    cout << hist_mean << " " << firstwidth << " " << firstwidther << " " << secondwidth << " " << secondwidther << " " << thirdwidth << " " << thirdwidther << " " << normtripg[0] << " " << normer[0] << " " << normtripg[1] << " " << normer[1] << " " << 1 - normtripg[0] - normtripg[1] << " " << normer[2] << " " << intestep << " " << intestep/(0.92528973*(gr_mu->Eval(hist_mean)*tripSeff)) << " " << tripww << " " << triperww << " " << significance  << endl;



    //model independent
    // cout << tripSeff << " " << tripSeffer << " " << gr_mu->Eval(singpartripgaus[i]->GetParameter(1)) << " " << 0.005 << " " <<  bcountfit/(dpinvmasslm[i]->GetBinWidth(0)) << " " << bcountfiter/(dpinvmasslm[i]->GetBinWidth(0)) << " " << intestep << " " << (intestep)/(0.977*tripSeff) << " " << intesteper/(0.977*tripSeff) << " " << tripww << " " << triperww << " " << significance << endl;



      TString signalplotname = signalfilename + string(".eps");
        C1->Print(signalplotname);

 TCanvas * C7 = new TCanvas("residuals"," ",10,10,800,800);
 //C7->Divide(2,1);

   C7->cd(1);
   sigres_alt->SetMarkerStyle(3);
   sigres_alt->Draw("*");
   //   C7->cd(2);
   // pull[i]->SetMarkerStyle(3);
   // pull[i]->Draw("*");

   TString signalresname = signalfilename + string("res.eps");
 C7->Print(signalresname);

 /*   TCanvas * C89 = new TCanvas("randomtest","",10,10,800,800);
 C89->Divide(2,1);
 C89->cd(1);
 h_pull[0]->Draw();

 for(int l = 1; l < 1; l++){
   h_pull[l]->Draw("same");
 }

 C89->cd(2);
 h_pull_res[0]->Draw();

 TString pullplotname = signalfilename + string("pull.eps");
 //C89->Print(pullplotname);
 */

}
