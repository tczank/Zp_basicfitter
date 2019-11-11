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

//2018/02/09
//## This new fit macro is to be used with histograms of the reduced dimuon mass, using muons from the list which was created using the Four Constraint Fitter, so no conservation cuts can be applied to this list ##//
//## It needs the gaussian width parametrization from newfit_v3_4.C ## //


//2018/02/13
//## New parametrization for double gaussian widths including new Z' masses 0.41 0.42 0.43 0.44 0.45

//2018/02/20
//## Calculation of the visible cross section for the background

//2018/02/22
//## Fit now using weighted width of the double gaussian function, it increased the detection efficiency slightly, but it also increased the number of background events on the same interval

//2018/02/23
//## Using single normalized gaussian to fit Signal MC sample and obtain signal shape, from fitting process we build a new function, normalized gaussian + 3rd order polynomial where the normalized gaussian part has same mean and width as the one fitted to the signal mc sample. Fitting the new function to the background mc sample we obtain a gaussian height and its error, from that we multiply it by the signal shape 1.64*width and by the binwidth

//2018/03/02
//## newest correction on the background histogram, as well as calculation of the number of observed events that was not present earlier

//2018/03/31
//## inclusion of jpsi background properly scaled fit region in auto zoom !!!!!! still missing eps generation for each fit


//2018/04/04
//## fitting of the jpsi peak to attempt yet another background rejection

//2018/04/06
//## getting number of observed event through the integration of  Nobs pdf

//2018/04/18
//##inclusion of pull and res distributions, trial of triple gaussian with only one mean, as well as, signalgaussianfit plus 3rd order poly replaced by triple gaussian

//2018/04/19
//##triple gaussian with single mean used for the detection efficiency

void newerfit_v2_2(TString signalfilename, TString backfilename, TString backpsimu, TString backpsipi) {

  TFile * br_fil = new TFile("Zp_BR.root");

  TGraph *gr_mu = new TGraph();

  gr_mu = (TGraph*)br_fil->Get("gr_mu");

  gStyle->SetOptFit();
  TFile *signal = new TFile(signalfilename);
  TFile *bg = new TFile(backfilename);
  TFile *bgpsimu = new TFile(backpsimu);
  TFile *bgpsipi = new TFile(backpsipi);
  TH1F *dpinvmasslm[4];
  TH1F *bginvmasslm[4];
  TH1F *bgpsimunvmasslm[4];
  TH1F *bgpsipinvmasslm[4];

  // TF1 *parawidth = new TF1("parawidth","cheb9",0.5,10);
  TF1 *parawidth = new TF1("parawidth", "pol9",0.,10.);

  //newest parametrization//
  parawidth->SetParameters(0.002562,-0.0001796,0.001217,-0.0008057,0.000269,-4.043e-05,9.323e-07,4.47e-07,-5.013e-08,1.652e-09);
  double parawidth_er[] ={1.063e-05,6.199e-06,1.157e-06,1.415e-07,1.641e-08,1.861e-09,2.058e-10,2.202e-11,2.245e-12,2.122e-13};
   parawidth->SetParErrors(parawidth_er);

  for(int i = 3; i < 4; i ++) {
    dpinvmasslm[i] = (TH1F*)signal->Get(TString::Format("h_babarpjpsicut_2"));
    bginvmasslm[i] = (TH1F*)bg->Get(TString::Format("h_babarpjpsicut_2"));
    bgpsimunvmasslm[i] = (TH1F*)bgpsimu->Get(TString::Format("h_babarpjpsicut_2"));
    bgpsipinvmasslm[i] = (TH1F*)bgpsipi->Get(TString::Format("h_babarpjpsicut_2"));
    // dpinvmasslm[i] = (TH1F*)signal->Get(TString::Format("h_babarkf_0"));
    //  bginvmasslm[i] = (TH1F*)bg->Get(TString::Format("h_babarkf_0"));
    //  bgpsimunvmasslm[i] = (TH1F*)bgpsimu->Get(TString::Format("h_babarkf_0"));
    //  bgpsipinvmasslm[i] = (TH1F*)bgpsipi->Get(TString::Format("h_babarkf_0"));
    // dpinvmasslm[i] = (TH1F*)signal->Get(TString::Format("h_mZpsing"));
    // bginvmasslm[i] = (TH1F*)bg->Get(TString::Format("h_mZpsing"));
    //dpinvmasslm[i] = (TH1F*)signal->Get(TString::Format("h_mRedsing"));
    //bginvmasslm[i] = (TH1F*)bg->Get(TString::Format("h_mRedsing"));
  }

  TCanvas *C1 = new TCanvas("C1", "", 10, 10, 800, 800);
  //  C1->Divide(2,2);
  for(int i = 3; i < 4; i ++) {
     dpinvmasslm[i]->SetLineColor(1);
     bginvmasslm[i]->SetLineColor(2);
     bginvmasslm[i]->Scale(0.327825);
     bgpsimunvmasslm[i]->Scale(2.98285);
     bgpsipinvmasslm[i]->Scale(2.98285);
     bginvmasslm[i]->Add(bgpsimunvmasslm[i]);
     bginvmasslm[i]->Add(bgpsipinvmasslm[i]);

      C1->cd(i + 1);

    dpinvmasslm[i]->Draw();
    bginvmasslm[i]->Draw("same");
  }

 TF1 * gaus_pol3[4];
  TF1 * gaus3[4];
  TF1 * pol0[4];
  TF1 * double_gaus[4];
  TF1 * triple_gaus[4];
  TF1 * quad_gaus[4];
  TF1 * gaus_pol1[4];
  TF1 * pol3[4];
  TF1 * gausn[4];
  TF1 * pol3n[4];
  TF1 * singpartripgaus[4];


  for(int i = 3; i < 4; i ++){
  double hist_mean = dpinvmasslm[i]->GetBinCenter(dpinvmasslm[i]->GetMaximumBin());
  double peakwidth = parawidth->Eval(hist_mean);
  //cout << "peak location is = " << hist_mean << endl;
   double entriesatmean = dpinvmasslm[i]->GetBinContent(dpinvmasslm[i]->GetMaximumBin());
   //  cout << " the entries at the mean is = " << entriesatmean << endl;
   double rms = dpinvmasslm[i]->GetBinWidth(1);
  //cout << "hist rms is = " << rms << endl;

  /// Scaling of the background due to the proper luminosity of 327.646 fb^-1 for the background and the full Belle luminosity of 977 fb^-1 /////
   //// New scaling factor recalculated by Igal is actually 334.396 fb^-1 (no) / ///////
    gaus_pol3[i]  = new TF1(TString::Format("gaus_pol3_%d",i),"gausn + [3] + [4]*x + [5]*x*x + [6]*x*x*x",hist_mean-0.25,hist_mean+0.25);
    gaus3[i]  = new TF1(TString::Format("gaus3_%d",i),"gausn", hist_mean-0.25,hist_mean+0.25);

    double_gaus[i] = new TF1(TString::Format("double_gaussian_%d",i)," (gausn(0) + gausn(3))/2", hist_mean-15*peakwidth, hist_mean+15*peakwidth);
    triple_gaus[i] = new TF1(TString::Format("triple_gaussian_%d",i)," (gausn(0)+ gausn(3) + gausn(6))",  hist_mean-15*peakwidth, hist_mean+15*peakwidth);
    // singpartripgaus[i] = new TF1(TString::Format("single_par_triple_gaus_%d",i),"(([0]/((sqrt(2*TMath::Pi()*[2]*[2]))))*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) +(((1-[0])*[3])/((sqrt(2*TMath::Pi()*[4]*[4]))))*exp(-((x-[1])*(x-[1]))/(2*[4]*[4])) + (((1-[0])*(1-[3]))/((sqrt(2*TMath::Pi()*[5]*[5]))))*exp(-((x-[1])*(x-[1]))/(2*[5]*[5])) )", hist_mean - 15*peakwidth, hist_mean + 15*peakwidth);
    //    singpartripgaus[i] = new TF1(TString::Format("single_par_triple_gaus_%d",i),"( ([0]/(3*(sqrt(2*TMath::Pi()*[2]*[2]))))*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) +(([0])/(3*(sqrt(2*TMath::Pi()*[3]*[3]))))*exp(-((x-[1])*(x-[1]))/(2*[3]*[3])) + ([0]/(3*(sqrt(2*TMath::Pi()*[4]*[4]))))*exp(-((x-[1])*(x-[1]))/(2*[4]*[4])) )", hist_mean - 15*peakwidth, hist_mean + 15*peakwidth);
      singpartripgaus[i] = new TF1(TString::Format("single_par_triple_gaus_%d",i),"(([0]/((sqrt(2*TMath::Pi()*[2]*[2]))))*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) +(([3])/((sqrt(2*TMath::Pi()*[4]*[4]))))*exp(-((x-[1])*(x-[1]))/(2*[4]*[4])) + (([6])/((sqrt(2*TMath::Pi()*[5]*[5]))))*exp(-((x-[1])*(x-[1]))/(2*[5]*[5])) )", hist_mean - 15*peakwidth, hist_mean + 15*peakwidth);

    quad_gaus[i] = new TF1(TString::Format("quad_gaussian_%d",i)," 0.25*(gausn(0)+ gausn(3) + gausn(6) + gausn(9) )",  hist_mean-15*peakwidth, hist_mean+15*peakwidth);
    gausn[i] = new TF1(TString::Format("norm_gaussian_%d",i)," gausn(0) ", hist_mean-10*peakwidth, hist_mean+10*peakwidth);

    gaus_pol1[i] = new TF1(TString::Format("gaussian_pol1_%d",i),"gausn + pol1", hist_mean-0.25, hist_mean+0.25);
    gaus_pol1[i]->SetLineColor(5);
    double_gaus[i]->SetLineColor(4);
    triple_gaus[i]->SetLineColor(8);
    quad_gaus[i]->SetLineColor(7);
    singpartripgaus[i]->SetLineColor(1);

    // defining parameters names //

    gaus_pol3[i]->SetParName(0,"gaus_height");
    gaus_pol3[i]->SetParName(1,"mean");
    gaus_pol3[i]->SetParName(2,"std");
    gaus_pol3[i]->SetParName(3,"x0");
    gaus_pol3[i]->SetParName(4,"x1");
    gaus_pol3[i]->SetParName(5,"x2");
    gaus_pol3[i]->SetParName(6,"x3");

    gausn[i]->SetParName(0,"gaus_height");
    gausn[i]->SetParName(1,"mean");
    gausn[i]->SetParName(2,"width");

    gausn[i]->SetParLimits(0,0.0,entriesatmean);
    gausn[i]->SetParLimits(1,hist_mean-3*peakwidth, hist_mean+3*peakwidth);
    gausn[i]->SetParLimits(2,rms,5*peakwidth);

    double_gaus[i]->SetParName(0,"gaus_height_1");
    double_gaus[i]->SetParName(1,"mean_1");
    double_gaus[i]->SetParName(2,"std_1");
    double_gaus[i]->SetParName(3,"gaus_height_2");
    double_gaus[i]->SetParName(4,"mean_2");
    double_gaus[i]->SetParName(5,"std_2");

    double_gaus[i]->SetParLimits(0,0.0,entriesatmean);//,dpinvmasslm[i]->GetEntries());
    double_gaus[i]->SetParLimits(1,hist_mean-3*peakwidth,hist_mean+3*peakwidth);
    double_gaus[i]->SetParLimits(2,rms,5*peakwidth);
    double_gaus[i]->SetParLimits(3,0.0,entriesatmean);//,dpinvmasslm[i]->GetEntries());
    double_gaus[i]->SetParLimits(4,hist_mean-3*peakwidth,hist_mean+3*peakwidth);
    double_gaus[i]->SetParLimits(5,rms,10*peakwidth);

    double_gaus[i]->SetRange(hist_mean-15*peakwidth,hist_mean+15*peakwidth);

    triple_gaus[i]->SetParName(0,"gaus_height_1");
    triple_gaus[i]->SetParName(1,"mean_1");
    triple_gaus[i]->SetParName(2,"std_1");
    triple_gaus[i]->SetParName(3,"gaus_height_2");
    triple_gaus[i]->SetParName(4,"mean_2");
    triple_gaus[i]->SetParName(5,"std_2");
    triple_gaus[i]->SetParName(6,"gaus_height_3");
    triple_gaus[i]->SetParName(7,"mean_3");
    triple_gaus[i]->SetParName(8,"std_3");

    triple_gaus[i]->SetParLimits(0,0.0,entriesatmean);//,dpinvmasslm[i]->GetEntries());
    triple_gaus[i]->SetParLimits(1,hist_mean-3*peakwidth,hist_mean+3*peakwidth);
    triple_gaus[i]->SetParLimits(2,rms,10*peakwidth);
    triple_gaus[i]->SetParLimits(3,0.0,entriesatmean);//,dpinvmasslm[i]->GetEntries());
    triple_gaus[i]->SetParLimits(4,hist_mean-3*peakwidth,hist_mean+3*peakwidth);
    triple_gaus[i]->SetParLimits(5,rms,10*peakwidth);
    triple_gaus[i]->SetParLimits(6,0.0,entriesatmean);//,dpinvmasslm[i]->GetEntries());
    triple_gaus[i]->SetParLimits(7,hist_mean-3*peakwidth,hist_mean+3*peakwidth);
    triple_gaus[i]->SetParLimits(8,rms,10*peakwidth);

    triple_gaus[i]->SetRange(hist_mean-20*peakwidth,hist_mean+20*peakwidth);

    singpartripgaus[i]->SetParName(0,"gaus_height1");
    singpartripgaus[i]->SetParName(1,"mean");
    singpartripgaus[i]->SetParName(2,"first_width");
    singpartripgaus[i]->SetParName(3,"gaus_height2");
    singpartripgaus[i]->SetParName(4,"second_width");
     singpartripgaus[i]->SetParName(5,"third_width");
     singpartripgaus[i]->SetParName(6,"gaus_height3");
     singpartripgaus[i]->SetNpx(1000);
     


    singpartripgaus[i]->SetParLimits(0,0.0,entriesatmean);
    singpartripgaus[i]->SetParLimits(1,hist_mean-3*peakwidth,hist_mean+3*peakwidth);
    singpartripgaus[i]->SetParLimits(2,rms,30*peakwidth);
    singpartripgaus[i]->SetParLimits(3,0.0,entriesatmean);
    singpartripgaus[i]->SetParLimits(4,rms,30*peakwidth);
     singpartripgaus[i]->SetParLimits(5,rms,30*peakwidth);
     singpartripgaus[i]->SetParLimits(6,0.0,entriesatmean);
     

    singpartripgaus[i]->SetRange(hist_mean-20*peakwidth,hist_mean+20*peakwidth);
    

    quad_gaus[i]->SetParName(0,"gaus_height_1");
    quad_gaus[i]->SetParName(1,"mean_1");
    quad_gaus[i]->SetParName(2,"std_1");
    quad_gaus[i]->SetParName(3,"gaus_height_2");
    quad_gaus[i]->SetParName(4,"mean_2");
    quad_gaus[i]->SetParName(5,"std_2");
    quad_gaus[i]->SetParName(6,"gaus_height_3");
    quad_gaus[i]->SetParName(7,"mean_3");
    quad_gaus[i]->SetParName(8,"std_3");
    quad_gaus[i]->SetParName(9,"gaus_height_4");
    quad_gaus[i]->SetParName(10,"mean_4");
    quad_gaus[i]->SetParName(11,"std_4");


    quad_gaus[i]->SetParLimits(0,0.0,entriesatmean);//,dpinvmasslm[i]->GetEntries());
    quad_gaus[i]->SetParLimits(1,hist_mean-3*peakwidth,hist_mean+3*peakwidth);
    quad_gaus[i]->SetParLimits(2,rms,5*peakwidth);
    quad_gaus[i]->SetParLimits(3,0.0,entriesatmean);//,dpinvmasslm[i]->GetEntries());
    quad_gaus[i]->SetParLimits(4,hist_mean-3*peakwidth,hist_mean+3*peakwidth);
    quad_gaus[i]->SetParLimits(5,rms,10*peakwidth);
    quad_gaus[i]->SetParLimits(6,0.0,entriesatmean);//,dpinvmasslm[i]->GetEntries());
    quad_gaus[i]->SetParLimits(7,hist_mean-3*peakwidth,hist_mean+3*peakwidth);
    quad_gaus[i]->SetParLimits(8,rms,10*peakwidth);
    quad_gaus[i]->SetParLimits(9,0.0,entriesatmean);//,dpinvmasslm[i]->GetEntries());
    quad_gaus[i]->SetParLimits(10,hist_mean-3*peakwidth,hist_mean+3*peakwidth);
    quad_gaus[i]->SetParLimits(11,rms,10*peakwidth);

    quad_gaus[i]->SetRange(hist_mean-15*peakwidth,hist_mean+15*peakwidth);

    gausn[i]->SetRange(hist_mean-15*peakwidth,hist_mean+15*peakwidth);

    gaus_pol1[i]->SetParName(0,"gaus_height");
    gaus_pol1[i]->SetParName(1,"mean");
    gaus_pol1[i]->SetParName(2,"std");
    gaus_pol1[i]->SetParName(3,"a");
    gaus_pol1[i]->SetParName(4,"b");

    gaus_pol1[i]->SetParLimits(0,0.0,dpinvmasslm[i]->GetEntries());
    gaus_pol1[i]->SetParLimits(1,hist_mean-0.05,hist_mean+0.05);
    gaus_pol1[i]->SetParLimits(2,0.001,0.05);

      gaus_pol3[i]->SetParLimits(0,0.0,dpinvmasslm[i]->GetEntries());
      gaus_pol3[i]->SetParLimits(1,hist_mean-10*rms,hist_mean+10*rms);
      gaus_pol3[i]->SetParLimits(2,0.001,1);
      pol0[i] = new TF1("pol0", "pol0",hist_mean-5*peakwidth,hist_mean+5*peakwidth);
      pol0[i]->SetRange(hist_mean-9*peakwidth,hist_mean+9*peakwidth);
      pol3[i] = new TF1("pol3", "pol3", hist_mean-15*peakwidth,hist_mean+15*peakwidth);
      pol3n[i] = new TF1("norm pol3", "[0]*([1]+ [2]*x +[3]*x*x +[4]*x*x*x)", hist_mean-15*peakwidth, hist_mean+15*peakwidth);

      pol3n[i]->FixParameter(0,1);

    }
  for(int i = 3; i<4; i++){
    gaus_pol3[i]->SetLineColor(1);
    pol0[i]->SetLineColor(2);
    pol3[i]->SetLineColor(3);
    pol3[i]->SetLineStyle(2);
    pol3n[i]->SetLineColor(5);
    pol3n[i]->SetLineStyle(2);
    gausn[i]->SetLineColor(2);
    }

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
  TH2F *sigres[4];
  TH2F *sigres_alt[4];
  TH2F *pull[4];
  

  //2018/02/23

  //##After fitting the single normalized gaussian, we build a new function with the defined parameters for the single normalized gaussian + a 3rd order polynomial, THEN, we perform a fit of this function over the Background and we take out the new height which is allowed to float, the error of the height and then we divide it by the bin width and multiply it by 1.4*height error

 for(int i = 3 ;i < 4;i++){

   double hist_mean = dpinvmasslm[i]->GetBinCenter(dpinvmasslm[i]->GetMaximumBin());
   //  cout << " the mean is = " << hist_mean << endl;
   double peakwidth = parawidth->Eval(hist_mean);
   //   TFitResultPtr s = bginvmasslm[i]->Fit(pol3[i],"MIERBQS+");
      nbinsdp = dpinvmasslm[i]->GetSize();
      // TFitResultPtr doubgausfit = dpinvmasslm[i]->Fit(double_gaus[i],"MIERBQS+");
      // doubgausfit->Print();
      TFitResultPtr gauspol1 = bginvmasslm[i]->Fit(pol3[i],"RBQS+");
      //  TFitResultPtr gausny = dpinvmasslm[i]->Fit(gausn[i],"RBQS+");
      TFitResultPtr normpol3 = bginvmasslm[i]->Fit(pol3n[i],"RBQS+");
      triple_gaus[i]->SetNpx(1000);
      TFitResultPtr tripgaus = dpinvmasslm[i]->Fit(triple_gaus[i],"LRBQS+");
      //    tripgaus->Print();
      // cout << " chisquare from tip gaus fit divided by ndf " << tripgaus->Chi2()/tripgaus->Ndf() << endl;
      //  cout << " the number of degrees of freedom is " << tripgaus->Ndf() << endl;
      //TFitResultPtr quadgaus = dpinvmasslm[i]->Fit(quad_gaus[i],"RBQS+");
      // quadgaus->Print();

      TFitResultPtr finalfit = dpinvmasslm[i]->Fit(singpartripgaus[i],"RBQS+");
      //      finalfit->Print();
      

      TF1 *only_one_gaus[3];

      for(int l = 0; l < 3; l++){
        only_one_gaus[l] = new TF1(TString::Format("only_one_gaus_%d",l),"gausn(0)",hist_mean-15*peakwidth, hist_mean+15*peakwidth);
      }

      ////////////################# weighted width with cubic root ########################////
         double tripalww = cbrt((TMath::Power(triple_gaus[i]->GetParameter(0)*triple_gaus[i]->GetParameter(2),3)+TMath::Power(triple_gaus[i]->GetParameter(3)*triple_gaus[i]->GetParameter(5),3) + TMath::Power(triple_gaus[i]->GetParameter(6)*triple_gaus[i]->GetParameter(8),3) )/(triple_gaus[i]->GetParameter(0)*triple_gaus[i]->GetParameter(0)*triple_gaus[i]->GetParameter(0) + triple_gaus[i]->GetParameter(3)*triple_gaus[i]->GetParameter(3)*triple_gaus[i]->GetParameter(3) + triple_gaus[i]->GetParameter(6)*triple_gaus[i]->GetParameter(6)*triple_gaus[i]->GetParameter(6)));
      double tripalerww = tripalww * cbrt(TMath::Power(triple_gaus[i]->GetParError(2),3) + TMath::Power(triple_gaus[i]->GetParError(5),3) + TMath::Power(triple_gaus[i]->GetParError(8),3));
      double tripalwh = cbrt((TMath::Power(triple_gaus[i]->GetParameter(0)*triple_gaus[i]->GetParameter(2),3)+TMath::Power(triple_gaus[i]->GetParameter(3)*triple_gaus[i]->GetParameter(5),3) + TMath::Power(triple_gaus[i]->GetParameter(6)*triple_gaus[i]->GetParameter(8),3))/(triple_gaus[i]->GetParameter(2)*triple_gaus[i]->GetParameter(2)*triple_gaus[i]->GetParameter(2) + triple_gaus[i]->GetParameter(5)*triple_gaus[i]->GetParameter(5)*triple_gaus[i]->GetParameter(5) + triple_gaus[i]->GetParameter(8) * triple_gaus[i]->GetParameter(8)*triple_gaus[i]->GetParameter(8)));
      double tripalerwh = tripalww * cbrt(TMath::Power(triple_gaus[i]->GetParError(0),3) + TMath::Power(triple_gaus[i]->GetParError(3),3) + TMath::Power(triple_gaus[i]->GetParError(6),3));

      //   cout << " the weighted width of the signal shape is " << tripalww << " +/- " << tripalerww << endl;
      // cout << " the weighted height of the signal shape is " << tripalwh << " +/- " << tripalerwh << endl;


      ////##### weighted width with square root ###############///// possibly right
      double tripww = sqrt((TMath::Power(triple_gaus[i]->GetParameter(0)*triple_gaus[i]->GetParameter(2),2)+TMath::Power(triple_gaus[i]->GetParameter(3)*triple_gaus[i]->GetParameter(5),2) + TMath::Power(triple_gaus[i]->GetParameter(6)*triple_gaus[i]->GetParameter(8),2) )/(triple_gaus[i]->GetParameter(0)*triple_gaus[i]->GetParameter(0) + triple_gaus[i]->GetParameter(3)*triple_gaus[i]->GetParameter(3) + triple_gaus[i]->GetParameter(6)*triple_gaus[i]->GetParameter(6)));
      double triperww = tripww * sqrt(TMath::Power(triple_gaus[i]->GetParError(2),2) + TMath::Power(triple_gaus[i]->GetParError(5),2) + TMath::Power(triple_gaus[i]->GetParError(8),2));
      double tripwh = sqrt((TMath::Power(triple_gaus[i]->GetParameter(0)*triple_gaus[i]->GetParameter(2),2)+TMath::Power(triple_gaus[i]->GetParameter(3)*triple_gaus[i]->GetParameter(5),2) + TMath::Power(triple_gaus[i]->GetParameter(6)*triple_gaus[i]->GetParameter(8),2))/(triple_gaus[i]->GetParameter(2)*triple_gaus[i]->GetParameter(2) + triple_gaus[i]->GetParameter(5)*triple_gaus[i]->GetParameter(5) +  triple_gaus[i]->GetParameter(8)*triple_gaus[i]->GetParameter(8)));
      double triperwh = tripww * sqrt(TMath::Power(triple_gaus[i]->GetParError(0),2) + TMath::Power(triple_gaus[i]->GetParError(3),2) + TMath::Power(triple_gaus[i]->GetParError(6),2));



      // cout << " the weighted width of the signal shape is " << tripww << " +/- " << triperww << endl;
      // cout << " the weighted height of the signal shape is " << tripwh << " +/- " << triperwh << endl;
      double tripger[] = {triperwh,peakwidth,triperww};

      TF1 * newgaus = new TF1("triple gaus sig", " gausn", hist_mean-9*tripww,hist_mean+9*tripww);
      newgaus->SetParameters(tripwh,hist_mean,tripww);
      newgaus->SetParErrors(tripger);

      double newgausint = newgaus->Integral(hist_mean-9*tripww,hist_mean+9*tripww);
      // cout << " the new gauss integral = " << newgausint << endl;
      double newgausinteff = (newgausint/(dpinvmasslm[i]->GetBinWidth(0)))/100000;

       double tripS = triple_gaus[i]->Integral(hist_mean-3*tripww, hist_mean+3*tripww);
      double tripSer = (triple_gaus[i]->IntegralError(hist_mean-3*tripww, hist_mean+3*tripww, tripgaus->GetParams(), tripgaus->GetCovarianceMatrix().GetMatrixArray()));
      double tripSeff = (tripS/(dpinvmasslm[i]->GetBinWidth(0)))/100000;
      double tripSeffer = (tripSer/(dpinvmasslm[i]->GetBinWidth(0)))/100000;


      //**** using single mean triple gaussian ***** ///////

      /*   double tripS = singpartripgaus[i]->Integral(hist_mean-3*tripww, hist_mean+3*tripww);
      double tripSer = (singpartripgaus[i]->IntegralError(hist_mean-3*tripww, hist_mean+3*tripww, tripgaus->GetParams(), tripgaus->GetCovarianceMatrix().GetMatrixArray()));
      double tripSeff = (tripS/(dpinvmasslm[i]->GetBinWidth(0)))/100000;
      double tripSeffer = (tripSer/(dpinvmasslm[i]->GetBinWidth(0)))/100000;*/

      TF1 * tripgauspol = new TF1("triple gaus with single mean plus 3rd poly", "(([0]/((sqrt(2*TMath::Pi()*[2]*[2]))))*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) +(([3])/((sqrt(2*TMath::Pi()*[4]*[4]))))*exp(-((x-[1])*(x-[1]))/(2*[4]*[4])) + (([6])/((sqrt(2*TMath::Pi()*[5]*[5]))))*exp(-((x-[1])*(x-[1]))/(2*[5]*[5]))) + [11]*([7] + [8]*x + [9]*x*x + [10]*x*x*x)", hist_mean - 15*peakwidth, hist_mean + 15*peakwidth);

      tripgauspol->SetLineColor(kOrange);
      tripgauspol->SetParName(0,"height1");
      tripgauspol->SetParName(1,"mean");
      tripgauspol->SetParName(2,"width1");
      tripgauspol->SetParName(3,"height2");
      tripgauspol->SetParName(4,"width2");
      tripgauspol->SetParName(5,"width3");
      tripgauspol->SetParName(6,"height3");
      tripgauspol->SetParName(7,"pol0");
      tripgauspol->SetParName(8,"pol1");
      tripgauspol->SetParName(9,"pol2");
      tripgauspol->SetParName(10,"pol3");
      tripgauspol->SetParName(11,"pol_frac");
      

      tripgauspol->SetParameter(0,singpartripgaus[i]->GetParameter(0));
      tripgauspol->FixParameter(1,singpartripgaus[i]->GetParameter(1));
      tripgauspol->SetParError(1,singpartripgaus[i]->GetParError(1));
      tripgauspol->FixParameter(2,singpartripgaus[i]->GetParameter(2));
      tripgauspol->SetParError(2,singpartripgaus[i]->GetParError(2));
      tripgauspol->SetParameter(3,singpartripgaus[i]->GetParameter(3));
      tripgauspol->FixParameter(4,singpartripgaus[i]->GetParameter(4));
      tripgauspol->SetParError(4,singpartripgaus[i]->GetParError(4));
      tripgauspol->FixParameter(5,singpartripgaus[i]->GetParameter(5));
      tripgauspol->SetParError(5,singpartripgaus[i]->GetParError(5));
      tripgauspol->SetParameter(6,singpartripgaus[i]->GetParameter(6));
      tripgauspol->FixParameter(7,pol3n[i]->GetParameter(1));
      tripgauspol->SetParError(7,pol3n[i]->GetParError(1));
      tripgauspol->FixParameter(8,pol3n[i]->GetParameter(2));
      tripgauspol->SetParError(8,pol3n[i]->GetParError(2));
      tripgauspol->FixParameter(9,pol3n[i]->GetParameter(3));
      tripgauspol->SetParError(9,pol3n[i]->GetParError(3));
      tripgauspol->FixParameter(10,pol3n[i]->GetParameter(4));
      tripgauspol->SetParError(10,pol3n[i]->GetParError(4));
      tripgauspol->SetParameter(11,pol3n[i]->GetParameter(0));
      


      TF1 * gausnpol3 = new TF1("norm double gaus with pol3", " gausn(0)  + [3]*([4]+[5]*x +[6]*x*x + [7]*x*x*x)", hist_mean-15*peakwidth, hist_mean+15*peakwidth);

      gausnpol3->SetLineColor(6);
      gausnpol3->SetParName(0, "height1");
      gausnpol3->SetParName(1, "mean1");
      gausnpol3->SetParName(2, "width1");
      gausnpol3->SetParName(3,"norm");
      gausnpol3->SetParName(4, "p[0]");
      gausnpol3->SetParName(5, "p[1]");
      gausnpol3->SetParName(6, "p[2]");
      gausnpol3->SetParName(7, "p[3]");

      gausnpol3->SetParameter(0,tripwh);
      gausnpol3->FixParameter(1,hist_mean);
      gausnpol3->SetParError(1,triple_gaus[i]->GetParError(1));
      gausnpol3->FixParameter(2,tripww);
      gausnpol3->SetParError(2,triperww);
      gausnpol3->FixParameter(4,pol3n[i]->GetParameter(1));
      gausnpol3->SetParError(4,pol3n[i]->GetParError(1));
      gausnpol3->FixParameter(5,pol3n[i]->GetParameter(2));
      gausnpol3->SetParError(5,pol3n[i]->GetParError(2));
      gausnpol3->FixParameter(6,pol3n[i]->GetParameter(3));
      gausnpol3->SetParError(6,pol3n[i]->GetParError(3));
      gausnpol3->FixParameter(7,pol3n[i]->GetParameter(4));
      gausnpol3->SetParError(7,pol3n[i]->GetParError(4));

      // You've stopped here 2018/04/03
      /*   for(int k = 0; k < 100; k++){
           gausnpol3->SetParameter(0,tripwh);*/
        TFitResultPtr gausnpol = bginvmasslm[i]->Fit(gausnpol3,"RBQS+");

         TFitResultPtr finalfitre = bginvmasslm[i]->Fit(tripgauspol,"RBQS+");
         //     finalfitre->Print();

   TF1 * siggaus = new TF1("number of events gauss", "gausn", hist_mean-15*peakwidth, hist_mean+15*peakwidth);

   siggaus->SetParameters(gausnpol3->GetParameter(0),hist_mean,tripww);
   siggaus->SetParError(0,gausnpol3->GetParError(0));
   siggaus->SetLineColor(1);
   siggaus->SetLineStyle(2);

   double sigS = siggaus->Integral(hist_mean-3*tripww,hist_mean+3*tripww);
   double sigSer = tripww*sigS;

   double Nobs = (siggaus->Integral(siggaus->GetParameter(1) - siggaus->GetParameter(2), siggaus->GetParameter(1)+siggaus->GetParameter(2)))/bginvmasslm[i]->GetBinWidth(0);
   double Nobser = Nobs*(siggaus->GetParError(0)/siggaus->GetParameter(0));

   //   cout << " number of observed events = " << Nobs << "+/-" << Nobser << endl;
   pol3[i]->Draw("same");
   siggaus->Draw("same");

   TF1 * sigtripgaus = new TF1("number of obs events trip gaus", " (([0]/((sqrt(2*TMath::Pi()*[2]*[2]))))*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) +(([3])/((sqrt(2*TMath::Pi()*[4]*[4]))))*exp(-((x-[1])*(x-[1]))/(2*[4]*[4])) + (([6])/((sqrt(2*TMath::Pi()*[5]*[5]))))*exp(-((x-[1])*(x-[1]))/(2*[5]*[5]))) ", hist_mean - 15*peakwidth, hist_mean + 15*peakwidth);

   sigtripgaus->SetParameters(tripgauspol->GetParameter(0), tripgauspol->GetParameter(1), tripgauspol->GetParameter(2), tripgauspol->GetParameter(3), tripgauspol->GetParameter(4), tripgauspol->GetParameter(5), tripgauspol->GetParameter(6));
   double sigtripgaus_er[] = {tripgauspol->GetParError(0), tripgauspol->GetParError(1), tripgauspol->GetParError(2), tripgauspol->GetParError(3), tripgauspol->GetParError(4), tripgauspol->GetParError(5), tripgauspol->GetParError(6)};
   sigtripgaus->SetParErrors(sigtripgaus_er);

   sigtripgaus->SetLineColor(kOrange);
   sigtripgaus->SetLineStyle(2);

   //  TFitResultPtr sigtripgaus_fit = bginvmasslm[i]->Fit(sigtripgaus,"RBQS+");
   

   TCanvas * C13 = new TCanvas("nobs_alt", " ",10,10,800,800);

   sigtripgaus->Draw();

   // double Nobs_alt = (sigtripgaus->Integral(sigtripgaus->GetParameter(1) - sigtripgaus->GetParameter(2), sigtripgaus->GetParameter(1) + sigtripgaus->GetParameter(2)))/bginvmasslm[i]->GetBinWidth(0);
   double Nobs_alt = (sigtripgaus->GetParameter(0)*sigtripgaus->GetParameter(2) + sigtripgaus->GetParameter(3)*sigtripgaus->GetParameter(4) + sigtripgaus->GetParameter(6)*sigtripgaus->GetParameter(5))/bginvmasslm[i]->GetBinWidth(0);
   
     // double Nobs_alter = abs(Nobs_alt)*((sigtripgaus->GetParError(0)+sigtripgaus->GetParError(3)+sigtripgaus->GetParError(6))/(sigtripgaus->GetParameter(0)+sigtripgaus->GetParameter(3)+sigtripgaus->GetParameter(6)));
   // double Nobs_alter = abs(Nobs_alt)*(sigtripgaus->GetParError(0)/sigtripgaus->GetParameter(0) + sigtripgaus->GetParError(3)/sigtripgaus->GetParameter(3) + sigtripgaus->GetParError(6)/sigtripgaus->GetParameter(6));
   double Nobs_alter = abs(Nobs_alt*((sigtripgaus->GetParError(0))/(sigtripgaus->GetParameter(0))));

   //  cout << " alternate number of events " << Nobs_alt << "+/-" << Nobs_alter << endl;

   TCanvas * C8 = new TCanvas("nobspdf"," ",10,10,800,800);
   // C8->Divide(2,1);
   C8->cd(1);

   TH1F * Nobspdf = new TH1F("PDF for the Number of Obs", "Number of Observed Events; Observed Events; entries;", 500, -100,100);

   /* TF1 * Nobspdff_alt = new TF1("PDF for the number of obs trip gaus", "gaus", Nobs_alt - 3*Nobs_alter, Nobs_alt + 3*Nobs_alter);
    Nobspdff_alt->SetParameters(1.,Nobs_alt,Nobs_alter);
    Nobspdff_alt->SetNpx(1000);
    Nobspdff_alt->Draw();


 double ninetylim_alt;
   double ninetylimer_alt;
   double corninetylim_alt = 0;
   double corninetylimer_alt;
   double infinitylim_alt;
   double testep_alt = 0.001*Nobs_alter;
   double intestep_alt = 0.00000001*Nobs_alter;
   double wholerange_alt = 100*Nobs_alter;
   double ninetyrange_alt = 0.9*wholerange_alt;
   // cout << " the integration step is = " << intestep_alt << " the wholerange is = " << wholerange_alt << " and the ninety range is = " << ninetyrange_alt << endl;

     //  cout << " the 90% limit point is = " << 1.64*Nobs_alter << endl;
   //ninetylim_alt = Nobspdff_alt->Integral(0.,1.64*Nobs_alter);
   //  cout << " ninetylim_alt " << ninetyrange_alt << endl;
   infinitylim_alt = Nobspdff_alt->Integral(0.,ninetyrange_alt);
     corninetylim_alt = ninetylim_alt/(Nobspdff_alt->Integral(0,1.64*Nobs_alter));
   ninetylimer_alt = sqrt(ninetylim_alt);
   corninetylimer_alt = ninetylimer_alt/(Nobspdff_alt->Integral(0,1.64*Nobs_alter));

   //   cout << " nobs is = " << Nobs_alt << " and the nobs error is = " << Nobs_alter << endl;

     while(corninetylim_alt <= 0.90){
       //  while(intestep_alt < 6*Nobs_alter){
       //      cout << " ninety_alt is " << ninetylim_alt << " and corninetylim_alt  " << Nobspdff_alt->Integral(0,100*Nobs_alter) << endl;
     ninetylim_alt = Nobspdff_alt->Integral(0,intestep_alt);
     intestep_alt = intestep_alt + 0.0001*Nobs_alter ;
     testep_alt = testep_alt + 0.01*Nobs_alter;
     corninetylim_alt = ninetylim_alt/(Nobspdff_alt->Integral(0,100*Nobs_alter));
     ninetylimer_alt = sqrt(ninetylim_alt);
     // cout << " the n90 point is " << intestep_alt << " testep_al " << corninetylim_alt << endl;
   }

     // cout << " ninety_alt is " << ninetylim_alt << " and corninetylim_alt  " << corninetylim_alt << endl;
     C8->cd(2);*/
   TF1 * Nobspdff = new TF1("PDF for the number of obs", "gaus", Nobs - 15*Nobser, Nobs + 15*Nobser);
   Nobspdff->SetParameters(1.,Nobs,Nobser);
   Nobspdff->SetNpx(1000);
   Nobspdff->Draw();

   double ninetylim;
   double ninetylimer;
   double corninetylim = 0;
   double corninetylimer;
   double infinitylim;
   double intestep = 0.0001*Nobser;
   double testep = 0.001*Nobser;
   double wholerange = 100*Nobser;
   double ninetyrange = 0.9*wholerange;
   // cout << " the integration step is = " << intestep << " the wholerange is = " << wholerange << " and the ninety range is = " << ninetyrange << endl;

   //  if(Nobs < 0 ){
     //  cout << " the 90% limit point is = " << 1.64*Nobser << endl;
    ninetylim = Nobspdff->Integral(0.,0.00000001*Nobser);
   infinitylim = Nobspdff->Integral(0.,ninetyrange);
   corninetylim = ninetylim/(Nobspdff->Integral(0,Nobser));
   ninetylimer = sqrt(ninetylim);
   corninetylimer = ninetylimer/(Nobspdff->Integral(0,Nobser));

   //  cout << " nobs is = " << Nobs << " and the nobs error is = " << Nobser << endl;

   while(corninetylim <= 0.90 ){
   //while(intestep < 15*Nobser){
     // cout << " ninety is " << ninetylim << " probability = " << corninetylim  << endl;
     ninetylim = Nobspdff->Integral(0,intestep);
     corninetylim = ninetylim/(Nobspdff->Integral(0,100*Nobser));
     intestep = intestep + 0.00001*Nobser ;
     testep = testep + 0.0001*Nobser;
      ninetylimer = sqrt(ninetylim);
     //    cout << " the n90 point is " << intestep << endl;
   }

   //   cout << " error for ninetylim = " << ninetylimer << " or the full range = " << sqrt(corninetylim) << endl;
   
   // cout << " ninety is " << ninetylim << " and corninetylim  " << corninetylim << endl;
   // cout << " the number of observed events was smaller than 0 = " << ninetylim << " +/- " << ninetylimer << endl;
    // cout << " the corrected number of observed events if smaller than 0 is " << corninetylim << endl;
  // }

    double B = (pol3[i]->Integral(triple_gaus[i]->GetParameter(1)-3*tripww,triple_gaus[i]->GetParameter(1)+3*tripww));
    double Ber= (pol3[i]->IntegralError(triple_gaus[i]->GetParameter(1)-3*tripww,triple_gaus[i]->GetParameter(1)+3*tripww,gauspol1->GetParams(), gauspol1->GetCovarianceMatrix().GetMatrixArray()));

    TAxis * xaxis = bginvmasslm[i]->GetXaxis();
    TAxis * yaxis = bginvmasslm[i]->GetYaxis();
    Int_t binxl = xaxis->FindBin(triple_gaus[i]->GetParameter(1)-15*tripww);
    Double_t binxlerror = bginvmasslm[i]->GetBinError(binxl);
    Int_t binxh = xaxis->FindBin( triple_gaus[i]->GetParameter(1)+15*tripww);
    Int_t binint = binxh - binxl;
    Double_t binxl_val = dpinvmasslm[i]->GetBinCenter(binxl);
    Double_t binxh_val = dpinvmasslm[i]->GetBinCenter(binxh);
    dpinvmasslm[i]->GetXaxis()->SetRange(binxl,binxh);
    double entries = 0;
    double intnorm = 1/(sqrt(2*TMath::Pi())*tripww);
    sigres[i]  = new TH2F(TString("signal_fit_residual_%d",i),"res;X_mass[GeV/c^2];residuals;", 100,binxl_val,binxh_val, 100, -300, 300);
    sigres_alt[i] = new TH2F(TString("background_fit_%d",i),"res;X_mass[GeV/c^2];residuals;", 100,binxl_val,binxh_val, 100, -300, 300);
    pull[i] = new TH2F(TString("pull_fit_%d",i),"pull;X_mass[GeV/c^2];pull;",100,binxl_val,binxh_val,100,-300,300);
    for(int k = binxl-50; k < binxh+50; k++){
      x = dpinvmasslm[i]->GetBinCenter(k);
      y = (dpinvmasslm[i]->GetBinContent(k) - triple_gaus[i]->Eval(x));//realwidth;
     sigres[i]->Fill(x,y);
     z = bginvmasslm[i]->GetBinCenter(k);
     w = (bginvmasslm[i]->GetBinContent(k) - gausnpol3->Eval(z));
     sigres_alt[i]->Fill(z,w);
     er = bginvmasslm[i]->GetBinError(k);
     pull[i]->Fill(z,w/er);
    }
     double sighistint = dpinvmasslm[i]->Integral(binxl,binxh);
    Double_t histonorm = sighistint*dpinvmasslm[i]->GetBinWidth(binxh);
    Double_t binxherror = bginvmasslm[i]->GetBinError(binxh);
    //Assuming that all possible combinations of final state muons are equivalent, the background count will be divided by 4, this will be corrected later by properly identifying the most energetic muons, this pair is the one coming out of the Z', this sould be checked in detail before filling the Z' invariant mass histograms//****************************************************************************************************************************
    double bcount = (bginvmasslm[i]->IntegralAndError(binxl,binxh,binxlerror));
    double bcounter = sqrt(bcount);
    double bcountfit = pol3[i]->Integral(binxl_val,binxh_val);
    double bcountfiter = pol3[i]->IntegralError(binxl_val,binxh_val,gauspol1->GetParams(),gauspol1->GetCovarianceMatrix().GetMatrixArray());
    double scount = (dpinvmasslm[i]->Integral(binxl,binxh));

    // assuming 0 Nobs using full number of background *** //
    //  cout << maybecorrectEff << " " << maybecorrectEffer << " " << gr_mu->Eval(double_gaus[i]->GetParameter(1)) << " " << 0.005 << " " << bcountfit/(dpinvmasslm[i]->GetBinWidth(0)) << " " << bcountfiter/(dpinvmasslm[i]->GetBinWidth(0)) << " " << (bcountfit/(dpinvmasslm[i]->GetBinWidth(0))) << endl; // " " << (bcountfit/(dpinvmasslm[i]->GetBinWidth(0)))/(0.977*(gr_mu->Eval(double_gaus[i]->GetParameter(1)))*maybecorrectEff) << endl;


    //## same from before but with triple gaussian fit
    //       cout << tripSeff << " " << tripSeffer << " " << gr_mu->Eval(double_gaus[i]->GetParameter(1)) << " " << 0.005 << " " << bcountfit/(dpinvmasslm[i]->GetBinWidth(0)) << " " << bcountfiter/(dpinvmasslm[i]->GetBinWidth(0)) << " " << (bcountfit/(dpinvmasslm[i]->GetBinWidth(0))) << endl; // " " << (bcountfit/(dpinvmasslm[i]->GetBinWidth(0)))/(0.977*(gr_mu->Eval(double_gaus[i]->GetParameter(1)))*maybecorrectEff) << endl;


    //  cout << maybecorrectEff << " " << maybecorrectEffer << " " << gr_mu->Eval(double_gaus[i]->GetParameter(1)) << " " << 0.005 << " " << Nobs << " " << Nobser << " " << Nobs << " " << Nobs/(0.977*(gr_mu->Eval(double_gaus[i]->GetParameter(1)))*maybecorrectEff) << endl;

    //**** correct calculated number of observed events with triple gaussian **** ///
    //         cout << tripSeff << " " << tripSeffer << " " << gr_mu->Eval(triple_gaus[i]->GetParameter(1)) << " " << 0.005 << " " <<  bcountfit/(dpinvmasslm[i]->GetBinWidth(0)) << " " << bcountfiter/(dpinvmasslm[i]->GetBinWidth(0)) << " " << Nobs << " " << (Nobs)/(0.977*(gr_mu->Eval(triple_gaus[i]->GetParameter(1)))*tripSeff) << " " << Nobser/(0.977*(gr_mu->Eval(triple_gaus[i]->GetParameter(1)))*tripSeff) << endl;


    //**** same but with 90% confidence limit **** ///
      cout << tripSeff << " " << tripSeffer << " " << gr_mu->Eval(triple_gaus[i]->GetParameter(1)) << " " << 0.005 << " " <<  bcountfit/(dpinvmasslm[i]->GetBinWidth(0)) << " " << bcountfiter/(dpinvmasslm[i]->GetBinWidth(0)) << " " << ninetylim << " " << (ninetylim)/(0.977*(gr_mu->Eval(triple_gaus[i]->GetParameter(1)))*tripSeff) << " " << ninetylimer/(0.977*(gr_mu->Eval(triple_gaus[i]->GetParameter(1)))*tripSeff) << endl;

        //** using single mean triple gaussian for the fit over the background to extract 90% lim :****////
    //  cout << tripSeff << " " << tripSeffer << " " << gr_mu->Eval(triple_gaus[i]->GetParameter(1)) << " " << 0.005 << " " <<  bcountfit/(dpinvmasslm[i]->GetBinWidth(0)) << " " << bcountfiter/(dpinvmasslm[i]->GetBinWidth(0)) << " " << ninetylim_alt << " " << (ninetylim_alt)/(0.977*(gr_mu->Eval(triple_gaus[i]->GetParameter(1)))*tripSeff) << " " << ninetylimer_alt/(0.977*(gr_mu->Eval(triple_gaus[i]->GetParameter(1)))*tripSeff) << endl;

        //**** inclusion of number of events error
    //  cout << tripSeff << " " << tripSeffer << " " << gr_mu->Eval(triple_gaus[i]->GetParameter(1)) << " " << 0.005 << " " <<  bcountfit/(dpinvmasslm[i]->GetBinWidth(0)) << " " << bcountfiter/(dpinvmasslm[i]->GetBinWidth(0)) << " " << ninetylim << " " << ninetylimer << " " << (ninetylim)/(0.977*(gr_mu->Eval(triple_gaus[i]->GetParameter(1)))*tripSeff) << " " << ninetylimer/(0.977*(gr_mu->Eval(triple_gaus[i]->GetParameter(1)))*tripSeff) << endl;


             //**** same but with 90% confidence limit corrected with normalization factor **** ///
    //  cout << tripSeff << " " << tripSeffer << " " << gr_mu->Eval(triple_gaus[i]->GetParameter(1)) << " " << 0.005 << " " <<  bcountfit/(dpinvmasslm[i]->GetBinWidth(0)) << " " << bcountfiter/(dpinvmasslm[i]->GetBinWidth(0)) << " " << corninetylim << " " << (corninetylim)/(0.977*(gr_mu->Eval(triple_gaus[i]->GetParameter(1)))*tripSeff) << " " << corninetylimer/(0.977*(gr_mu->Eval(triple_gaus[i]->GetParameter(1)))*tripSeff) << endl;


    // double_gaus[i]->GetParameter(1) << endl;

        TString signalplotname = signalfilename + string(".eps");
        C1->Print(signalplotname);

        /*    pol3[i]->Draw("same");
              siggaus->Draw("same");*/
    //            cout << " the siggaus height is = " << siggaus->GetParameter(0) << " the mean is " << siggaus->GetParameter(1) << " and the width is " << siggaus->GetParameter(2) << endl;
 }
 TCanvas * C7 = new TCanvas("residuals"," ",10,10,800,800);
 C7->Divide(2,1);

 for(int i = 3; i < 4 ; i++){
   C7->cd(1);
   // sigres[i]->SetMarkerStyle(3);
   //sigres[i]->Draw("*");
   // C7->cd(2);
   sigres_alt[i]->SetMarkerStyle(3);
   sigres_alt[i]->Draw("*");
   C7->cd(2);
   pull[i]->SetMarkerStyle(3);
   pull[i]->Draw("*");
 }
 TString signalresname = signalfilename + string("res.pdf");
 /* TCanvas *C3 = new TCanvas("C3", "", 10, 10, 800, 800);

  C3->Divide(3,1);
  C3->cd(1);
  fom->Draw();
  C3->cd(2);
  foma->Draw();
  C3->cd(3);
  ponzifom->Draw();*/

  TCanvas * C9 = new TCanvas("jpsicut", "", 10, 10, 800, 800);
 C9->cd(1);



 TFitResultPtr tripgausbak = bginvmasslm[3]->Fit(triple_gaus[3],"LRBQS+");

 //triple_gaus[3]->Set

 triple_gaus[3]->SetLineColor(2);

 double tripww = cbrt((TMath::Power(triple_gaus[3]->GetParameter(0)*triple_gaus[3]->GetParameter(2),3)+TMath::Power(triple_gaus[3]->GetParameter(3)*triple_gaus[3]->GetParameter(5),3) + TMath::Power(triple_gaus[3]->GetParameter(6)*triple_gaus[3]->GetParameter(8),3) )/(triple_gaus[3]->GetParameter(0)*triple_gaus[3]->GetParameter(0)*triple_gaus[3]->GetParameter(0) + triple_gaus[3]->GetParameter(3)*triple_gaus[3]->GetParameter(3)*triple_gaus[3]->GetParameter(3) + triple_gaus[3]->GetParameter(6)*triple_gaus[3]->GetParameter(6)*triple_gaus[3]->GetParameter(6)));
 double triperww = tripww * cbrt(TMath::Power(triple_gaus[3]->GetParError(2),3) + TMath::Power(triple_gaus[3]->GetParError(5),3) + TMath::Power(triple_gaus[3]->GetParError(8),3));

 cout << " the weighted width of the background shape is " << tripww << "+/-" << triperww << endl;

 bginvmasslm[3]->Draw();




}
