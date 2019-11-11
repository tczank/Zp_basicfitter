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

//2018/04/24
//##triple gaussian with single mean used to define weighted width and then define pdf to extract 90% upper limit

//2018//04/25
//## inclusion of the toy montecarlo study and eps figures printing command

//2018/06/05
//## gen id width parametrization

//2018/06/06
//## weighted width peak widht parametrization

//#### TO CALCULATE THE DIFFERENT DETECTION EFFICIENCIES FOR HADRONBJ and TAU AB SKIMS ( REMEMBER TO CHANGE THE NAME OF THE HISTOGRAMS!!!!!!!!!!!!!!)

void newerfit_v3(TString signalfilename, TString backfilename, TString backpsimu, TString backpsipi) {

  TFile * br_fil = new TFile("Zp_BR.root");

  TGraph *gr_mu = new TGraph();

  gr_mu = (TGraph*)br_fil->Get("gr_mu");

  gStyle->SetOptFit();
  // gROOT->SetBatch(1);
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
  TF1 *parawei = new TF1("parawei", "pol5",0., 10.);

  parawei->SetParameters(0.002846,0.003433,-0.001757,0.0005654,-7.123e-05,3.021e-06);
  double parawei_er[] = {1.326e-06,4.6e-06,3.722e-06,1.145e-06,1.466e-07,6.59e-09};
  parawei->SetParErrors(parawei_er);

   parawidth->SetParameters(0.001918,0.002347,-0.001771,0.0008198,-0.0001845,2.042e-05,-6.609e-07,-8.292e-08,8.569e-09,-2.27e-10);
   double parawidth_er[] ={6.664e-06,4.493e-06,8.904e-07,1.094e-07,1.275e-08,1.45e-09,1.607e-10,1.722e-11,1.758e-12,1.665e-13};
   parawidth->SetParErrors(parawidth_er);

   /* h_hadronBeff[i] = new TH1F(TString::Format("h_hadronBeff_%d",i),"Reduced Di muon mass [GeV/c^{2}];m_{R}[GeV/c^{2}];entries;",10500,0.,10.5);
  h_hadronBeff[i]->Sumw2();
  h_tauskimeff[i] = new TH1F(TString::Format("h_tauskimeff_%d",i),"Reduced Di muon mass [GeV/c^{2}];m_{R}[GeV/c^{2}];entries;",10500,0.,10.5);
  h_tauskimeff[i]->Sumw2();
  h_mycombitrigeff[i] = new TH1F(TString::Format("h_mycombitrigeff_%d",i),"Reduced Di muon mass [GeV/c^{2}];m_{R}[GeV/c^{2}];entries;",10500,0.,10.5);
  h_mycombitrigeff[i]->Sumw2();
  }*/


  for(int i = 3; i < 4; i ++) {
    //    dpinvmasslm[i] = (TH1F*)signal->Get(TString::Format("h_babarpjpsicut_5"));
        dpinvmasslm[i] = (TH1F*)signal->Get(TString::Format("h_genidredmu_0"));
    bginvmasslm[i] = (TH1F*)bg->Get(TString::Format("h_babarpjpsicut_6"));
    bgpsimunvmasslm[i] = (TH1F*)bgpsimu->Get(TString::Format("h_babarpjpsicut_6"));
    bgpsipinvmasslm[i] = (TH1F*)bgpsipi->Get(TString::Format("h_babarpjpsicut_6"));
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
  TF1 * trigtripo;

  double hist_mean;
  double peakwidth;
  double entriesatmean;
  double rms;
  double std_dev;
  double significance;
  double thirdpolchi;
  double gausnobschi;

  for(int i = 3; i < 4; i ++){
  hist_mean = dpinvmasslm[i]->GetBinCenter(dpinvmasslm[i]->GetMaximumBin());
  //peakwidth = parawidth->Eval(hist_mean);
  peakwidth = parawei->Eval(hist_mean);
  entriesatmean = dpinvmasslm[i]->GetBinContent(dpinvmasslm[i]->GetMaximumBin());
  rms = dpinvmasslm[i]->GetBinWidth(1);
  std_dev = dpinvmasslm[i]->GetStdDev(1);

  // cout << " the std is " << std_dev << " and the rms is " << rms << endl;

  /// Scaling of the background due to the proper luminosity of 327.646 fb^-1 for the background and the full Belle luminosity of 977 fb^-1 /////
   //// New scaling factor recalculated by Igal is actually 334.396 fb^-1 (no) / ///////
    gaus_pol3[i]  = new TF1(TString::Format("gaus_pol3_%d",i),"gausn + [3] + [4]*x + [5]*x*x + [6]*x*x*x",hist_mean-0.25,hist_mean+0.25);
    gaus3[i]  = new TF1(TString::Format("gaus3_%d",i),"gausn", hist_mean-0.25,hist_mean+0.25);

    double_gaus[i] = new TF1(TString::Format("double_gaussian_%d",i)," (gausn(0) + gausn(3))/2", hist_mean-20*peakwidth, hist_mean+20*peakwidth);
    triple_gaus[i] = new TF1(TString::Format("triple_gaussian_%d",i)," (gausn(0)+ gausn(3) + gausn(6))",  hist_mean-20*peakwidth, hist_mean+20*peakwidth);
      singpartripgaus[i] = new TF1(TString::Format("single_par_triple_gaus_%d",i),"(([0]/((sqrt(2*TMath::Pi()*[2]*[2]))))*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) +(([3])/((sqrt(2*TMath::Pi()*[4]*[4]))))*exp(-((x-[1])*(x-[1]))/(2*[4]*[4])) + (([6])/((sqrt(2*TMath::Pi()*[5]*[5]))))*exp(-((x-[1])*(x-[1]))/(2*[5]*[5])) )", hist_mean - 20*peakwidth, hist_mean + 20*peakwidth);

    quad_gaus[i] = new TF1(TString::Format("quad_gaussian_%d",i)," 0.25*(gausn(0)+ gausn(3) + gausn(6) + gausn(9) )",  hist_mean-20*peakwidth, hist_mean+20*peakwidth);
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

    double_gaus[i]->SetRange(hist_mean-20*peakwidth,hist_mean+20*peakwidth);

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
    singpartripgaus[i]->SetParLimits(2,rms,20*peakwidth);
    singpartripgaus[i]->SetParLimits(3,0.0,entriesatmean);
    singpartripgaus[i]->SetParLimits(4,rms,20*peakwidth);
     singpartripgaus[i]->SetParLimits(5,rms,20*peakwidth);
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

    quad_gaus[i]->SetRange(hist_mean-20*peakwidth,hist_mean+20*peakwidth);

    gausn[i]->SetRange(hist_mean-20*peakwidth,hist_mean+20*peakwidth);

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
      pol3[i] = new TF1("pol3", "pol3", hist_mean-20*peakwidth,hist_mean+20*peakwidth);
      pol3n[i] = new TF1("norm pol3", "[0]*([1]+ [2]*x +[3]*x*x +[4]*x*x*x)", hist_mean-20*peakwidth, hist_mean+20*peakwidth);
      pol3n[i]->SetRange(hist_mean-20*peakwidth,hist_mean+20*peakwidth);
      
      pol3n[i]->FixParameter(0,1);

    }

    TCanvas *C2 = new TCanvas("C2", "", 10, 10, 800, 800);

  for(int i = 3; i<4; i++){
    gaus_pol3[i]->SetLineColor(1);
    pol0[i]->SetLineColor(2);
    pol3[i]->SetLineColor(5);
    pol3[i]->SetLineStyle(2);
    pol3n[i]->SetLineColor(3);
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

  TRandom *r1 = new TRandom();

  TH1D * h_pull[10000];
  TH1D * h_pull_res[2];
  h_pull_res[0]= new TH1D("pull distribution_0","pull;pull;entries;", 100,-5,5);
  h_pull_res[1]= new TH1D("pull distribution_1","pull;pull;entries;", 1000,-200,500);

  //2018/02/23

  //##After fitting the single normalized gaussian, we build a new function with the defined parameters for the single normalized gaussian + a 3rd order polynomial, THEN, we perform a fit of this function over the Background and we take out the new height which is allowed to float, the error of the height and then we divide it by the bin width and multiply it by 1.4*height error



  double normtripg[3];
  double heightlist[3];

  for(int i = 3 ;i < 4;i++){

   hist_mean = dpinvmasslm[i]->GetBinCenter(dpinvmasslm[i]->GetMaximumBin());
   std_dev = dpinvmasslm[i]->GetStdDev();
   //cout << " the range is " << hist_mean - 3*std_dev << " " << hist_mean + 3*std_dev << endl;
   // peakwidth = parawidth->Eval(hist_mean);
   peakwidth = parawei->Eval(hist_mean);
   TFitResultPtr normpol3 = bginvmasslm[i]->Fit(pol3n[i],"RBQS+");
      thirdpolchi = normpol3->Chi2();

      triple_gaus[i]->SetNpx(1000);

      singpartripgaus[i]->SetNpx(1000);
      TFitResultPtr finalfit = dpinvmasslm[i]->Fit(singpartripgaus[i],"RBQS+");
      // finalfit->Print();

      heightlist[0] = singpartripgaus[i]->GetParameter(0);
      heightlist[1] = singpartripgaus[i]->GetParameter(3);
      heightlist[2] = singpartripgaus[i]->GetParameter(6);

      for (int l = 0; l < 3 ; l++){
        normtripg[l] = (heightlist[l]/(heightlist[0] + heightlist[1] +heightlist[2]));
        cout << " the normtripg is " << normtripg[l] << " and the heightlist is " << heightlist[l] << endl;
      }

      TF1 *only_one_gaus[3];

      for(int l = 0; l < 3; l++){
        only_one_gaus[l] = new TF1(TString::Format("only_one_gaus_%d",l),"gausn(0)",hist_mean-20*peakwidth, hist_mean+20*peakwidth);
      }

      //### John's correction following BN 1229 ################################ //

      double f_one = singpartripgaus[i]->GetParameter(0)/(singpartripgaus[i]->GetParameter(3) + singpartripgaus[i]->GetParameter(0)+singpartripgaus[i]->GetParameter(6));
      double f_two = singpartripgaus[i]->GetParameter(3)/(singpartripgaus[i]->GetParameter(3) + singpartripgaus[i]->GetParameter(0)+singpartripgaus[i]->GetParameter(6));

      double tripww = sqrt(f_one * singpartripgaus[i]->GetParameter(2)*singpartripgaus[i]->GetParameter(2) + f_two * singpartripgaus[i]->GetParameter(4)*singpartripgaus[i]->GetParameter(4) + (1.0 - f_one - f_two)*singpartripgaus[i]->GetParameter(5)*singpartripgaus[i]->GetParameter(5));

      ////##### weighted width with square root ###############///// possibly right
      //        double tripww = sqrt((TMath::Power(singpartripgaus[i]->GetParameter(0)*singpartripgaus[i]->GetParameter(2),2)+TMath::Power(singpartripgaus[i]->GetParameter(3)*singpartripgaus[i]->GetParameter(4),2) + TMath::Power(singpartripgaus[i]->GetParameter(6)*singpartripgaus[i]->GetParameter(5),2) )/(singpartripgaus[i]->GetParameter(0)*singpartripgaus[i]->GetParameter(0) + singpartripgaus[i]->GetParameter(3)*singpartripgaus[i]->GetParameter(3) + singpartripgaus[i]->GetParameter(6)*singpartripgaus[i]->GetParameter(6)));
       double triperww = tripww * sqrt(TMath::Power(singpartripgaus[i]->GetParError(2),2) + TMath::Power(singpartripgaus[i]->GetParError(4),2) + TMath::Power(singpartripgaus[i]->GetParError(5),2));
      double tripwh = sqrt((TMath::Power(singpartripgaus[i]->GetParameter(0)*singpartripgaus[i]->GetParameter(2),2)+TMath::Power(singpartripgaus[i]->GetParameter(3)*singpartripgaus[i]->GetParameter(4),2) + TMath::Power(singpartripgaus[i]->GetParameter(6)*singpartripgaus[i]->GetParameter(5),2))/(singpartripgaus[i]->GetParameter(2)*singpartripgaus[i]->GetParameter(2) + singpartripgaus[i]->GetParameter(4)*singpartripgaus[i]->GetParameter(4) +  singpartripgaus[i]->GetParameter(5)*singpartripgaus[i]->GetParameter(5)));
      double triperwh = tripww * sqrt(TMath::Power(singpartripgaus[i]->GetParError(0),2) + TMath::Power(singpartripgaus[i]->GetParError(4),2) + TMath::Power(singpartripgaus[i]->GetParError(5),2));


      TAxis * xaxis = bginvmasslm[i]->GetXaxis();
      TAxis * yaxis = bginvmasslm[i]->GetYaxis();
      Int_t binxl = xaxis->FindBin(singpartripgaus[i]->GetParameter(1)-20*peakwidth);
      Double_t binxlerror = bginvmasslm[i]->GetBinError(binxl);
      Int_t binxh = xaxis->FindBin( singpartripgaus[i]->GetParameter(1)+20*peakwidth);
      Double_t binxl_val = dpinvmasslm[i]->GetBinCenter(binxl);
      Double_t binxh_val = dpinvmasslm[i]->GetBinCenter(binxh);
      if(binxl_val < 0.00000){
        binxl = xaxis->FindBin(0.0);
        binxl_val = dpinvmasslm[i]->GetBinCenter(binxl);
      }
      if(binxh_val > 10.500){
        binxh = xaxis->FindBin(10.5);
        binxh_val = dpinvmasslm[i]->GetBinCenter(binxh);
    }
      Int_t binint = binxh - binxl;
     Double_t sigwinbin = (binxh_val - binxl_val)/bginvmasslm[i]->GetBinWidth(0);
      //  cout << "the limits of the signal window are " << binxl_val << " and " << binxh_val << endl;
      // cout << "the limits of the signal window after checkin are " << binxl_val << " and " << binxh_val << endl;
      double_t entriesinint = bginvmasslm[i]->Integral(binxl,binxh);
      //   cout << " number of backgroun entries in the signal windows is " << entriesinint << endl;

      double tripger[] = {triperwh,peakwidth,triperww};

      TF1 * newgaus = new TF1("triple gaus sig", " gausn", hist_mean-20*peakwidth,hist_mean+20*peakwidth);
      newgaus->SetParameters(tripwh,hist_mean,peakwidth);
      newgaus->SetParErrors(tripger);

      double newgausint = newgaus->Integral(hist_mean-9*peakwidth,hist_mean+9*peakwidth);
      double newgausinteff = (newgausint/(dpinvmasslm[i]->GetBinWidth(0)))/100000;

       double tripS = singpartripgaus[i]->Integral(hist_mean-3*peakwidth, hist_mean+3*peakwidth);
      double tripSer = (singpartripgaus[i]->IntegralError(hist_mean-3*peakwidth, hist_mean+3*peakwidth, finalfit->GetParams(), finalfit->GetCovarianceMatrix().GetMatrixArray()));
      //      double tripSeff = (tripS/(dpinvmasslm[i]->GetBinWidth(0)))/100000;
      //   double tripSeffer = (tripSer/(dpinvmasslm[i]->GetBinWidth(0)))/100000;

      // CORRECTED EFFICIENCY WITH ISR FACTOR 0.848148 +/- 0.094711635
      // double tripSeff = (dpinvmasslm[i]->Integral(binxl ,binxh))/100000;
       double tripSeff = 0.848148*(dpinvmasslm[i]->Integral(binxl ,binxh))/100000;
      //cout << " the new efficiency is " << tripSeff << " between " << hist_mean - 3*tripww << " and " << hist_mean + 3*tripww << endl;
      double tripSeffer = sqrt((tripSeff*(1-tripSeff))/100000);


      TF1 * gausnpol3 = new TF1("norm double gaus with pol3", " gaus(0)  + [3]*([4]+[5]*x +[6]*x*x + [7]*x*x*x)", hist_mean-20*peakwidth, hist_mean+20*peakwidth);

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
      gausnpol3->SetParError(1,singpartripgaus[i]->GetParError(1));
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

      gausnpol3->SetRange(hist_mean-20*peakwidth,hist_mean+20*peakwidth);
        TFitResultPtr gausnpol = bginvmasslm[i]->Fit(gausnpol3,"RBQS+");

        gausnobschi = gausnpol->Chi2();


        //TF1 * tripgtripo = new TF1("tripgtripo", " [12]*(([0]/((sqrt(2*TMath::Pi()*[2]*[2]))))*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) +(([3])/((sqrt(2*TMath::Pi()*[4]*[4]))))*exp(-((x-[1])*(x-[1]))/(2*[4]*[4])) + (([6])/((sqrt(2*TMath::Pi()*[5]*[5]))))*exp(-((x-[1])*(x-[1]))/(2*[5]*[5]))) + [7]*([8]+[9]*x +[10]*x*x + [11]*x*x*x)", hist_mean-20*peakwidth, hist_mean+20*peakwidth);


        TF1 * tripgtripo = new TF1("tripgtripo", " [12]*((([0]))*exp(-((x-[1])*(x-[1]))/(2*[2]*[2])) +(([3]))*exp(-((x-[1])*(x-[1]))/(2*[4]*[4])) + (([6]))*exp(-((x-[1])*(x-[1]))/(2*[5]*[5]))) + [7]*([8]+[9]*x +[10]*x*x + [11]*x*x*x)", hist_mean-20*peakwidth, hist_mean+20*peakwidth);
        

        tripgtripo->FixParameter(0,normtripg[0]);
        tripgtripo->FixParameter(1,singpartripgaus[i]->GetParameter(1));
        tripgtripo->SetParError(1,singpartripgaus[i]->GetParError(1));
         tripgtripo->FixParameter(2,singpartripgaus[i]->GetParameter(2));
         //  tripgtripo->SetParameter(2,singpartripgaus[i]->GetParameter(2));
        tripgtripo->SetParError(2,triperww);
      tripgtripo->FixParameter(3,normtripg[1]);
      tripgtripo->SetParError(3,singpartripgaus[i]->GetParError(3));
      tripgtripo->FixParameter(4,singpartripgaus[i]->GetParameter(4));
      //tripgtripo->SetParameter(4,singpartripgaus[i]->GetParameter(4));
          tripgtripo->FixParameter(5,singpartripgaus[i]->GetParameter(5));
      //  tripgtripo->SetParameter(5,singpartripgaus[i]->GetParameter(5));
                                 tripgtripo->FixParameter(6,normtripg[2]);
      tripgtripo->SetParameter(7,pol3n[i]->GetParameter(0));
                tripgtripo->FixParameter(8,pol3n[i]->GetParameter(1));
                tripgtripo->FixParameter(9,pol3n[i]->GetParameter(2));
                tripgtripo->SetParError(9,pol3n[i]->GetParError(2));
                tripgtripo->FixParameter(10,pol3n[i]->GetParameter(3));
                tripgtripo->SetParError(10,pol3n[i]->GetParError(3));
                tripgtripo->FixParameter(11,pol3n[i]->GetParameter(4));
                tripgtripo->SetParError(11,pol3n[i]->GetParError(4));

                tripgtripo->SetRange(hist_mean-20*peakwidth,hist_mean+20*peakwidth);
                TFitResultPtr trigtrip = bginvmasslm[i]->Fit(tripgtripo,"RBQS+");
                trigtrip->Print();


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
      gausnpol3_forpull->SetParError(1,triple_gaus[i]->GetParError(1));
      gausnpol3_forpull->FixParameter(2,tripww);
      gausnpol3_forpull->SetParError(2,triperww);
      gausnpol3_forpull->FixParameter(4,pol3n[i]->GetParameter(1));
      gausnpol3_forpull->SetParError(4,pol3n[i]->GetParError(1));
      gausnpol3_forpull->FixParameter(5,pol3n[i]->GetParameter(2));
      gausnpol3_forpull->SetParError(5,pol3n[i]->GetParError(2));
      gausnpol3_forpull->FixParameter(6,pol3n[i]->GetParameter(3));
      gausnpol3_forpull->SetParError(6,pol3n[i]->GetParError(3));
      gausnpol3_forpull->FixParameter(7,pol3n[i]->GetParameter(4));
      gausnpol3_forpull->SetParError(7,pol3n[i]->GetParError(4));

              gausnpol3_forpull->SetParameter(4,pol3n[i]->GetParameter(1));
      gausnpol3_forpull->SetParError(4,pol3n[i]->GetParError(1));
      gausnpol3_forpull->SetParameter(5,pol3n[i]->GetParameter(2));
      gausnpol3_forpull->SetParError(5,pol3n[i]->GetParError(2));
      gausnpol3_forpull->SetParameter(6,pol3n[i]->GetParameter(3));
      gausnpol3_forpull->SetParError(6,pol3n[i]->GetParError(3));
      gausnpol3_forpull->SetParameter(7,pol3n[i]->GetParameter(4));
      gausnpol3_forpull->SetParError(7,pol3n[i]->GetParError(4));


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


                double sigtrigS = tripsiggaus->Integral(hist_mean-3*tripsiggaus->GetParameter(4),hist_mean+3*tripsiggaus->GetParameter(4))/bginvmasslm[i]->GetBinWidth(0);
                 double sigtrigSer = (tripsiggaus->GetParError(9)/tripsiggaus->GetParameter(9))*sigtrigS;

                  cout << " the number of observed events for the trip = " << sigtrigS << " +/- " << sigtrigSer << endl;

         TF1 * siggaus = new TF1("number of events gauss", "gaus", hist_mean-20*peakwidth, hist_mean+20*peakwidth);

   siggaus->SetParameters(gausnpol3->GetParameter(0),hist_mean,tripww);
   siggaus->SetParError(0,gausnpol3->GetParError(0));
   siggaus->SetLineColor(1);
   siggaus->SetLineStyle(2);

   double sigS = siggaus->Integral(hist_mean-3*tripww,hist_mean+3*tripww)*1000;
   double sigSer = (siggaus->GetParError(0)/siggaus->GetParameter(0))*sigS;

 cout << " the number of observed events for the single = " << sigS << " +/- " << sigSer << endl;

    double Nobs = (siggaus->Integral(siggaus->GetParameter(1) - siggaus->GetParameter(2), siggaus->GetParameter(1)+siggaus->GetParameter(2)))/bginvmasslm[i]->GetBinWidth(0);
    double Nobser = Nobs*(siggaus->GetParError(0)/siggaus->GetParameter(0));

 //    double Nobs = sigtrigS;
 //  double Nobser = sigtrigSer;

   significance = ROOT::Math::Sign(Nobs)*sqrt(thirdpolchi - gausnobschi);
   //   cout << " the number of observed events is " << ROOT::Math::Sign(Nobs) << endl;

   pol3[i]->Draw("same");
   siggaus->Draw("same");


   TCanvas * C8 = new TCanvas("nobspdf"," ",10,10,800,800);
   C8->cd(1);


   TF1 * Nobspdff = new TF1("PDF for the number of observed events", " gaus ", Nobs - 40*Nobser, Nobs + 40*Nobser);
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
    double B = (pol3n[i]->Integral(singpartripgaus[i]->GetParameter(1)-3*tripww,singpartripgaus[i]->GetParameter(1)+3*tripww));
    double Ber= (pol3n[i]->IntegralError(singpartripgaus[i]->GetParameter(1)-3*tripww,singpartripgaus[i]->GetParameter(1)+3*tripww,normpol3->GetParams(), normpol3->GetCovarianceMatrix().GetMatrixArray()));

    dpinvmasslm[i]->GetXaxis()->SetRange(binxl,binxh);
    double entries = 0;
    double intnorm = 1/(sqrt(2*TMath::Pi())*tripww);
    sigres[i]  = new TH2F(TString("signal_fit_residual_%d",i),"res;X_mass[GeV/c^2];residuals;", 100,binxl_val,binxh_val, 100, -300, 300);
    sigres_alt[i] = new TH2F(TString("background_fit_%d",i),"res;X_mass[GeV/c^2];residuals;", 100,binxl_val,binxh_val, 100, -300, 300);
    pull[i] = new TH2F(TString("pull_fit_%d",i),"pull;X_mass[GeV/c^2];pull;",100,binxl_val,binxh_val,100,-300,300);
    for(int k = binxl-50; k < binxh+50; k++){
      x = dpinvmasslm[i]->GetBinCenter(k);
      y = (dpinvmasslm[i]->GetBinContent(k) - singpartripgaus[i]->Eval(x));//realwidth;
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
    double bcountfit = pol3n[i]->Integral(binxl_val,binxh_val);
    double bcountfiter = pol3n[i]->IntegralError(binxl_val,binxh_val,normpol3->GetParams(),normpol3->GetCovarianceMatrix().GetMatrixArray());
    double scount = (dpinvmasslm[i]->Integral(binxl,binxh));


    //model dependent with the cross section
       cout << tripSeff << " " << tripSeffer << " " << gr_mu->Eval(singpartripgaus[i]->GetParameter(1)) << " " << 0.005 << " " <<  bcountfit/(dpinvmasslm[i]->GetBinWidth(0)) << " " << bcountfiter/(dpinvmasslm[i]->GetBinWidth(0)) << " " << intestep << " " << (intestep)/(0.977*(gr_mu->Eval(singpartripgaus[i]->GetParameter(1)))*tripSeff) << " " << intesteper/(0.977*(gr_mu->Eval(singpartripgaus[i]->GetParameter(1)))*tripSeff) << " " << tripww << " " << triperww << " " << significance << endl;

    //model independent
    // cout << tripSeff << " " << tripSeffer << " " << gr_mu->Eval(singpartripgaus[i]->GetParameter(1)) << " " << 0.005 << " " <<  bcountfit/(dpinvmasslm[i]->GetBinWidth(0)) << " " << bcountfiter/(dpinvmasslm[i]->GetBinWidth(0)) << " " << intestep << " " << (intestep)/(0.977*tripSeff) << " " << intesteper/(0.977*tripSeff) << " " << tripww << " " << triperww << " " << significance << endl;



      TString signalplotname = signalfilename + string(".eps");
        C1->Print(signalplotname);

 }
 TCanvas * C7 = new TCanvas("residuals"," ",10,10,800,800);
 //C7->Divide(2,1);

 for(int i = 3; i < 4 ; i++){
   C7->cd(1);
   sigres_alt[i]->SetMarkerStyle(3);
   sigres_alt[i]->Draw("*");
   //   C7->cd(2);
   // pull[i]->SetMarkerStyle(3);
   // pull[i]->Draw("*");
 }
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
