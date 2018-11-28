#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1F.h"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFile.h"
//#include "CrystalBall.C"

//2018/02/09
//## This new fit macro is to be used with histograms of the reduced dimuon mass, using muons from the list which was created using the Four Constraint Fitter, so no conservation cuts can be applied to this list ##//
//## The name of the histogram containing the dimuon combinations from the Four Constraint Fitter is h_babarkf
//## This macro is to parametrize the double gaussian function used to fit

//2018/02/13
//## A new parametrization of the double gaussian was made using the new Z' samples for 0.41 0.42 0.43 0.44 0.45

void init_para_mc(TString signalfilename, TString backfilename) {

  TFile * br_fil = new TFile("Zp_BR.root");

  TGraph *gr_mu = new TGraph();

  gr_mu = (TGraph*)br_fil->Get("gr_mu");


  gStyle->SetOptFit();
  TFile *signal = new TFile(signalfilename);
  TFile *bg = new TFile(backfilename);
  TH1F *dpinvmasslm[4];
  TH1F *bginvmasslm[4];

  //  TF1 *parawidth = new TF1("parawidth","cheb9",0,10);
  TF1 *parawidth = new TF1("parawidth", "pol9",0.,10.);
  /*parawidth->SetParameters(0.00231,0.001966,-0.0004716,8.519e-05,-4.81e-06,-1.143e-07,1.164e-08,4.464e-10,-4.195e-11,7.143e-13);
    double parawidth_er[] ={6.678e-06,4.1e-06,4.06e-07,2.559e-08,1.497e-09,8.514e-11,4.714e-12,2.52e-13,1.283e-14,6.049e-16};*/
  /*  parawidth->SetParameters(0.01438, -0.01662, 0.006814, -0.001145, 7.408e-05, 9.52e-07, -2.175e-07, -4.053e-09, 7.668e-10, -1.686e-11);
  double parawidth_er[] ={8.317e-06, 1.835e-06, 1.105e-07, 6.409e-09, 3.618e-10, 2.005e-11, 1.094e-12, 5.836e-14, 3.048e-15, 1.54e-16};
  parawidth->SetParErrors(parawidth_er);*/

  /* parawidth->SetParameters(0.003874, -0.002188, 0.002796, -0.0002553, -0.00029, 0.0001182, -1.765e-05, 1.002e-06, 2.143e-09, -1.455e-09);
  double parawidth_er[] ={1.782e-05, 1.578e-05, 4.624e-06, 8.121e-07, 1.199e-07, 1.572e-08, 1.934e-09, 2.245e-10, 2.391e-11, 2.177e-12};
  parawidth->SetParErrors(parawidth_er);*/

  /*  parawidth->SetParameters(0.003572,-0.003622,0.005392,-0.003133,0.0009685,-0.0001549,1.022e-05,2.504e-07,-6.811e-08,2.527e-09);
  double parawidth_er[] ={1.211e-05,6.302e-06,1.129e-06,1.423e-07,1.669e-08,1.909e-09,2.132e-10,2.298e-11,2.353e-12,2.219e-13};
  parawidth->SetParErrors(parawidth_er);*/


  //newest parametrization//
  parawidth->SetParameters(0.002562,-0.0001796,0.001217,-0.0008057,0.000269,-4.043e-05,9.323e-07,4.47e-07,-5.013e-08,1.652e-09);
  double parawidth_er[] ={1.063e-05,6.199e-06,1.157e-06,1.415e-07,1.641e-08,1.861e-09,2.058e-10,2.202e-11,2.245e-12,2.122e-13};
  parawidth->SetParErrors(parawidth_er);


  for(int i = 3; i < 4; i ++) {
    dpinvmasslm[i] = (TH1F*)signal->Get(TString::Format("h_babarkf_0"));
    bginvmasslm[i] = (TH1F*)bg->Get(TString::Format("h_babarkf_0"));
  }

  TCanvas *C1 = new TCanvas("C1", "", 10, 10, 800, 800);
  //  C1->Divide(2,2);
  for(int i = 3; i < 4; i ++) {
     dpinvmasslm[i]->SetLineColor(1);
     bginvmasslm[i]->SetLineColor(2);

      C1->cd(i + 1);

    dpinvmasslm[i]->Draw();
    bginvmasslm[i]->Draw("same");
  }

  TF1 * gaus_pol3[4];
  TF1 * gaus3[4];
  TF1 * pol0[4];
  TF1 * pol3[4];
  TF1 * double_gaus[4];
  TF1 * gaus_pol1[4];

  for(int i = 3; i < 4; i ++){
  double hist_mean = dpinvmasslm[i]->GetBinCenter(dpinvmasslm[i]->GetMaximumBin());
  double peakwidth = parawidth->Eval(hist_mean);
  //cout << "peak location is = " << hist_mean << endl;
   double entriesatmean = dpinvmasslm[i]->GetBinContent(dpinvmasslm[i]->GetMaximumBin());
   //  cout << " the entries at the mean is = " << entriesatmean << endl;
   double rms = dpinvmasslm[i]->GetBinWidth(1);
  //cout << "hist rms is = " << rms << endl;
  /// Scaling of the background due to the proper luminosity of 327.646 fb^-1 for the background and the full Belle luminosity of 977 fb^-1 /////
    bginvmasslm[i]->Scale(3.27825);
    gaus_pol3[i]  = new TF1(TString::Format("gaus_pol3_%d",i),"gausn + [3] + [4]*x + [5]*x*x + [6]*x*x*x",hist_mean-0.25,hist_mean+0.25);
    gaus3[i]  = new TF1(TString::Format("gaus3_%d",i),"gausn",hist_mean-0.25,hist_mean+0.25);

    double_gaus[i] = new TF1(TString::Format("double_gaussian_%d",i)," gausn(0) + gausn(3) ", hist_mean-5*peakwidth, hist_mean+5*peakwidth);

    gaus_pol1[i] = new TF1(TString::Format("gaussian_pol1_%d",i),"gausn + pol1", hist_mean-0.25, hist_mean+0.25);
    gaus_pol1[i]->SetLineColor(5);
    double_gaus[i]->SetLineColor(4);
    // defining parameters names //

    gaus_pol3[i]->SetParName(0,"gaus_height");
    gaus_pol3[i]->SetParName(1,"mean");
    gaus_pol3[i]->SetParName(2,"std");
    gaus_pol3[i]->SetParName(3,"x0");
    gaus_pol3[i]->SetParName(4,"x1");
    gaus_pol3[i]->SetParName(5,"x2");
    gaus_pol3[i]->SetParName(6,"x3");

    double_gaus[i]->SetParName(0,"gaus_height_1");
    double_gaus[i]->SetParName(1,"mean_1");
    double_gaus[i]->SetParName(2,"std_1");
    double_gaus[i]->SetParName(3,"gaus_height_2");
    double_gaus[i]->SetParName(4,"mean_2");
    double_gaus[i]->SetParName(5,"std_2");

    /* double_gaus[i]->SetParameter(0,entriesatmean);//,dpinvmasslm[i]->GetEntries());
    double_gaus[i]->SetParameter(1,hist_mean);
    double_gaus[i]->SetParameter(2,peakwidth);
    double_gaus[i]->SetParameter(3,entriesatmean);//,dpinvmasslm[i]->GetEntries());
    double_gaus[i]->SetParameter(4,hist_mean);
    double_gaus[i]->SetParameter(5,peakwidth);*/


    double_gaus[i]->SetParLimits(0,0.0,entriesatmean);//,dpinvmasslm[i]->GetEntries());
    double_gaus[i]->SetParLimits(1,hist_mean-3*peakwidth,hist_mean+3*peakwidth);
    double_gaus[i]->SetParLimits(2,rms,15*peakwidth);
    // double onestwidth = double_gaus[i]->GetParameter(2);
    //cout << " the current 1st width of the double gaussian is = " << onestwidth << endl;
    double_gaus[i]->SetParLimits(3,0.0,entriesatmean);//,dpinvmasslm[i]->GetEntries());
    double_gaus[i]->SetParLimits(4,hist_mean-2*peakwidth,hist_mean+3*peakwidth);
    double_gaus[i]->SetParLimits(5,rms,15*peakwidth);

    //  double_gaus[i]->GetParameters()->Release();
    //cout << " the lower range is = " << hist_mean -10*peakwidth << " the higher range is = " << hist_mean +10*peakwidth << endl;
    double_gaus[i]->SetRange(hist_mean-6*peakwidth,hist_mean+6*peakwidth);

    gaus_pol1[i]->SetParName(0,"gaus_height");
    gaus_pol1[i]->SetParName(1,"mean");
    gaus_pol1[i]->SetParName(2,"std");
    gaus_pol1[i]->SetParName(3,"a");
    gaus_pol1[i]->SetParName(4,"b");

    gaus_pol1[i]->SetParLimits(0,0.0,dpinvmasslm[i]->GetEntries());
    gaus_pol1[i]->SetParLimits(1,hist_mean-15*peakwidth,hist_mean+15*peakwidth);
    gaus_pol1[i]->SetParLimits(2,0.001,0.05);

      gaus_pol3[i]->SetParLimits(0,0.0,dpinvmasslm[i]->GetEntries());
      gaus_pol3[i]->SetParLimits(1,hist_mean-15*peakwidth,hist_mean+15*peakwidth);
      gaus_pol3[i]->SetParLimits(2,0.001,1);
      pol0[i] = new TF1("pol0", "pol0",hist_mean-5*peakwidth,hist_mean+5*peakwidth);
      pol3[i] = new TF1("pol3", "pol3", hist_mean-6*peakwidth,hist_mean+6*peakwidth);}
  for(int i = 3; i<4; i++){
    gaus_pol3[i]->SetLineColor(1);
    pol0[i]->SetLineColor(2);
    pol3[i]->SetLineColor(3);
  }

    //  TFitter * doublegausfit;
    //doublegausfit->SetFunction(double_gaus[3],false);
    //doublegausfit->ParSettings(0)->Release();
    //  TMinuit *gMinuit = new TMinuit(6);
    // gMinuit->SetFCN(double_gaus[3]);




  TH1F *fom = new TH1F("figure of merit with B only","fom;# of loose muons;",4,0,3);
  TH1F *foma = new TH1F( "figure of merit with S+B","foma;# of loose muons;",4,0,3);
  TH1F *ponzifom = new TH1F("figure of merit following ponzi model","ponzi fom; # of loose muons;",4,0,3);

  Int_t nbinsdp;
  Int_t nbinsbg;

  Double_t x;
  Double_t y;

  double gaus_pol3_ParEr[7];
  double par_sigfit[3];
  TH2F *sigres[4];
 for(int i = 3 ;i < 4;i++){
    // dpinvmasslm[i]->Rebin(1);
    // bginvmasslm[i]->Rebin(1);
    //TFitResultPtr r = dpinvmasslm[i]->Fit(gaus_pol3[i],"MIERBQS+");
   //TFitResultPtr s = bginvmasslm[i]->Fit(pol3[i],"RBQS+");
      nbinsdp = dpinvmasslm[i]->GetSize();

      // Here there is the following error//
      /** FUNCTION VALUE DOES NOT SEEM TO DEPEND ON ANY OF THE 1 VARIABLE PARAMETERS.
      VERIFY THAT STEP SIZES ARE BIG ENOUGH AND CHECK FCN LOGIC **/
      TFitResultPtr doubgausfit = dpinvmasslm[i]->Fit(double_gaus[i],"RBQS+");
      TFitResultPtr gauspol1 = bginvmasslm[i]->Fit(pol3[i],"RBQS+");

      //cout << " sigma from gaus pol1 is = " << gaus_pol1[i]->GetParameter(2) << endl;
      //cout << " The number of bins in sig hist = " << nbinsdp << endl;
      for(int n = 0; n < 3; n++){
        double par_transf;
        double par_err_transf;
        par_transf = gaus_pol3[i]->GetParameter(n);
        par_err_transf = gaus_pol3[i]->GetParError(n);
        gaus3[i]->SetParameter(n, par_transf);
        gaus3[i]->SetParError(n,par_err_transf);
      }
      double S = (double_gaus[i]->Integral(double_gaus[i]->GetParameter(1)-5*double_gaus[i]->GetParameter(2),double_gaus[i]->GetParameter(1)+5*double_gaus[i]->GetParameter(2)));
      double Ser = (double_gaus[i]->IntegralError(double_gaus[i]->GetParameter(1)-5*double_gaus[i]->GetParameter(2),double_gaus[i]->GetParameter(1)+5*double_gaus[i]->GetParameter(2),doubgausfit->GetParams(), doubgausfit->GetCovarianceMatrix().GetMatrixArray()));
    double B = (pol3[i]->Integral(double_gaus[i]->GetParameter(1)-5*double_gaus[i]->GetParameter(2),double_gaus[i]->GetParameter(1)+5*double_gaus[i]->GetParameter(2)));
    double Ber= (pol3[i]->IntegralError(double_gaus[i]->GetParameter(1)-5*double_gaus[i]->GetParameter(2),double_gaus[i]->GetParameter(1)+5*double_gaus[i]->GetParameter(2),gauspol1->GetParams(), gauspol1->GetCovarianceMatrix().GetMatrixArray()));

    TAxis * xaxis = bginvmasslm[i]->GetXaxis();
    TAxis * yaxis = bginvmasslm[i]->GetYaxis();
    Int_t binxl = xaxis->FindBin( double_gaus[i]->GetParameter(1)-3*double_gaus[i]->GetParameter(2));
    Double_t binxlerror = bginvmasslm[i]->GetBinError(binxl);
    Int_t binxh = xaxis->FindBin( double_gaus[i]->GetParameter(1)+3*double_gaus[i]->GetParameter(2));
    Int_t binint = binxh - binxl;
    Double_t binxl_val = dpinvmasslm[i]->GetBinCenter(binxl);
    Double_t binxh_val = dpinvmasslm[i]->GetBinCenter(binxh);
    double entries = 0;
    double intnorm = 1/(sqrt(2*TMath::Pi())*double_gaus[i]->GetParameter(2));
    sigres[i]  = new TH2F(TString("signal_fit_residual_%d",i),"res;X_mass[GeV/c^2];residuals;", 100,binxl_val,binxh_val, 100, -300, 300);
  for(int k = binxl; k < binxh+1; k++){
      x = dpinvmasslm[i]->GetBinCenter(k);
      y = dpinvmasslm[i]->GetBinContent(k) - double_gaus[i]->Eval(x);
     sigres[i]->Fill(x,y);}
  //    C7->cd(i+1);
     double sighistint = dpinvmasslm[i]->Integral(binxl,binxh);
    Double_t histonorm = sighistint*dpinvmasslm[i]->GetBinWidth(binxh);
    Double_t binxherror = bginvmasslm[i]->GetBinError(binxh);
    //Assuming that all possible combinations of final state muons are equivalent, the background count will be divided by 4, this will be corrected later by properly identifying the most energetic muons, this pair is the one coming out of the Z', this sould be checked in detail before filling the Z' invariant mass histograms//****************************************************************************************************************************
    double bcount = (bginvmasslm[i]->IntegralAndError(binxl,binxh,binxlerror));
    double bcounter = sqrt(bcount);
    double bcountfit = pol3[i]->Integral(binxl_val,binxh_val);
    double bcountfiter = pol3[i]->IntegralError(binxl_val,binxh_val,gauspol1->GetParams(),gauspol1->GetCovarianceMatrix().GetMatrixArray());
    double scount = (dpinvmasslm[i]->Integral(binxl,binxh));
    double maybecorrectS = S/(dpinvmasslm[i]->GetBinWidth(binxh));
    //cout << " maybecorrectS is  = " << maybecorrectS << endl;
    //cout << " the bin width is = " << dpinvmasslm[i]->GetBinWidth(binxh);
    double maybecorrectSer =Ser/(dpinvmasslm[i]->GetBinWidth(binxh));
    double maybecorrectEff = (maybecorrectS)/100000;
    double maybecorrectEffer = (maybecorrectSer)/100000;
    double FOM = S / sqrt(B);
    double FOMA = S / sqrt(S+B);
    double PONZIFOM = S/(3*3/8 + 9*3*3/13 + 3*sqrt(B)+3*sqrt(3*3 + 4 * 3 * sqrt(B)+4*sqrt(B))/2);
    ponzifom->SetBinContent(i+1,PONZIFOM);
    fom->SetBinContent(i+1,FOM);
    foma->SetBinContent(i+1,FOMA);
    double efficiency = S*100/100000;
    double geoeff = (S*100000)/(binxh-binxl)/100000;
    double geoeffer = (Ser*100000)/(binxh-binxl)/100000;
    //    cout << "the 1st gaussian mean is " << double_gaus[i]->GetParameter(1) << "+/-" << double_gaus[i]->GetParError(1) <<  " the 1st gaussian std is " << double_gaus[i]->GetParameter(2) << "+/-" << double_gaus[i]->GetParError(2) << endl;
    //cout << "the 2nd gaussian mean is " << double_gaus[i]->GetParameter(4) << "+/-" << double_gaus[i]->GetParError(4) << " the 2nd gaussian std is " << double_gaus[i]->GetParameter(5) << "+/-" << double_gaus[i]->GetParError(5) << endl;
    //cout << maybecorrectEff << " " << maybecorrectEffer << " " << gr_mu->Eval(double_gaus[i]->GetParameter(1)) << " " << 0.005 << " " << bcountfit/(dpinvmasslm[i]->GetBinWidth(0)) << " " << bcountfiter/(dpinvmasslm[i]->GetBinWidth(0)) << " " << bcountfit/(dpinvmasslm[i]->GetBinWidth(0)) << endl;

    //1st sigma/[to get double gaussian fitted function First Mean and First width]/////////////////
    // cout << double_gaus[i]->GetParameter(1) << " "<< double_gaus[i]->GetParError(1) <<" " << double_gaus[i]->GetParameter(2) << " "  << double_gaus[i]->GetParError(2) << endl;

      //2nd sigma//[to get double gaussian fitted function Second Mean and Second width]//////////////
          cout  << double_gaus[i]->GetParameter(4) << " " << double_gaus[i]->GetParError(4) << " " << double_gaus[i]->GetParameter(5) << " " << double_gaus[i]->GetParError(5) << endl;

        //1st height/////[to get double gaussian fitted function First Mean and First Height ]
    //      cout << double_gaus[i]->GetParameter(1) << " " << double_gaus[i]->GetParError(1) << " " << double_gaus[i]->GetParameter(0) << " " << double_gaus[i]->GetParError(0) << endl;

        //2nd height/////[to get double gaussian fitted function Second mean and second height]
    //     cout << double_gaus[i]->GetParameter(4) << " " << double_gaus[i]->GetParError(4) << " " << double_gaus[i]->GetParameter(3) << " " << double_gaus[i]->GetParError(3) << endl;
 }
 // TCanvas * C7 = new TCanvas("residuals"," ",10,10,800,800);

 //C7->Divide(2,2);
 /*for(int i = 3; i < 4 ; i++){
   C7->cd(i+1);
   sigres[i]->SetMarkerStyle(3);
   sigres[i]->Draw("*");
   }*/
 TString signalplotname = signalfilename + string(".pdf");
 TString signalresname = signalfilename + string("res.pdf");
 //C1->Print(signalplotname);
 //C7->Print(signalresname);
 // TCanvas *C3 = new TCanvas("C3", "", 10, 10, 800, 800);

  /* C3->Divide(3,1);
  C3->cd(1);
  fom->Draw();
  C3->cd(2);
  foma->Draw();
  C3->cd(3);
  ponzifom->Draw();*/

}