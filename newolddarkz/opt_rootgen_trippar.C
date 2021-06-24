#include "Riostream.h"

void opt_rootgen_trippar() {
  ifstream in;
  in.open(Form("opt_trippar_all.dat"));

  Float_t x,fw,fwer,sw,swer,tw,twer,frac1,frac1er,frac2,frac2er,frac3,frac3er,nobs,ww,wwer;

  double_t m[107];
  Double_t ex[107] = {0};

  double_t firstw[107];
  double_t firstwer[107];

  double_t secondw[107];
  double_t secondwer[107];

  double_t thirdw[107];
  double_t thirdwer[107];

  double_t tf1[107];
  double_t tf1er[107];

  double_t tf2[107];
  double_t tf2er[107];

  Double_t tf3[107];
  Double_t tf3er[107];

  double_t weiwi[107];
 double_t weiwier[107];

  Int_t nlines = 0;
   TFile *f = new TFile("opt_trippar_all_2.root","RECREATE");

   while (nlines < 107) {
     in >> x >> fw >> fwer >> sw >> swer >> tw >> twer >> frac1 >> frac1er >> frac2 >> frac2er >> frac3 >> frac3er >>  ww >> wwer;
     //  if (!in.good()) break;
     if (nlines < 108) printf(" x=%8f, fw=%8f, fwer=%8f, sw=%8f, swer=%8f, tw=%8f, twer=%8f, frac1=%8f, frac1er=%8f, frac2=%8f, frac2er=%8f, frac3=%8f, frac3er=%8f, ww=%8f, wwer=%8f", x,fw, fwer, sw, swer, tw, twer, frac1, frac1er, frac2, frac2er, frac3, frac3er, ww, wwer );
     m[nlines] = x;

     firstw[nlines] = fw;
     firstwer[nlines] = fwer;
     secondw[nlines] = sw;
     secondwer[nlines] = swer;
     thirdw[nlines] = tw;
     thirdwer[nlines] = twer;
     if(fw > sw){
     firstw[nlines] = sw;
     firstwer[nlines] = swer;
     secondw[nlines] = fw;
     secondwer[nlines] = fwer;}
     if(fw > tw){
       firstw[nlines] = tw;
       firstwer[nlines] = twer;
       thirdw[nlines] = fw;
       thirdwer[nlines] = fwer;
     }
     if(sw > tw){
       secondw[nlines] = tw;
       secondwer[nlines] = twer;
       thirdw[nlines] = sw;
       thirdwer[nlines] = swer;
     }
     tf1[nlines] = frac1;
     tf1er[nlines] = frac1er;
     tf2[nlines] = frac2;
     tf2er[nlines] = frac2er;
     tf3[nlines] = frac3;
     tf3er[nlines] = frac3er;
     if(frac1 > frac2){
     tf1[nlines] = frac2;
     tf1er[nlines] = frac2er;
     tf2[nlines] = frac1;
     tf2er[nlines] = frac1er;}
     if(frac1 > frac3){
       tf1[nlines] = frac3;
       tf1er[nlines] = frac3er;
       tf3[nlines] = frac1;
       tf3er[nlines] = frac1er;
     }
     if(frac2 > frac3){
       tf2[nlines] = frac3;
       tf2er[nlines] = frac3er;
       tf3[nlines] = frac2;
       tf3er[nlines] = frac2er;
     }
     weiwi[nlines] = ww;
     weiwier[nlines] = wwer;
     nlines++;}
     printf(" found %d points\n",nlines);


  TGraphErrors * gr_frac1 = new TGraphErrors(nlines,m,tf1,ex,tf1er);
  gr_frac1->SetTitle("triple gaussian fit function normalized fraction 1 parameterized;M[GeV/c^{2}];Fraction1;");

  TGraphErrors * gr_frac2 = new TGraphErrors(nlines,m,tf2,ex,tf2er);
  gr_frac2->SetTitle("triple gaussian fit function normalized fraction 2 parameterized;M[GeV/c^{2}];Fraction2;");

  TGraphErrors * gr_frac3 = new TGraphErrors(nlines,m,tf3,ex,tf3er);
  gr_frac3->SetTitle("triple gaussian fit function normalized fraction 3 parameterized;M[GeV/c^{2}];Fraction3;");

  TGraphErrors * gr_ww = new TGraphErrors(nlines,m,weiwi,ex,weiwier);
  gr_ww->SetTitle("triple gaussian fit function weighted width parametrized;M[GeV/c^{2}];Weighted Width[GeV/c^{2}];");

  TGraphErrors * gr_fw = new TGraphErrors(nlines,m,firstw,ex,firstwer);
  gr_fw->SetTitle("1st Width from triple gaussian ff;M[GeV/c^{2}];#sigma_{1}");

  TGraphErrors * gr_sw = new TGraphErrors(nlines,m,secondw,ex,secondwer);
  gr_sw->SetTitle("2nd Width from triple gaussian ff;M[GeV/c^{2}];#sigma_{2}");

  TGraphErrors * gr_tw = new TGraphErrors(nlines,m,thirdw,ex,thirdwer);
  gr_tw->SetTitle("3rd Width from triple gaussian ff;M[GeV/c^{2}];#sigma_{2}");

  gr_frac1->SetName("gr_frac1");
  gr_frac1->Write();

  gr_frac2->SetName("gr_frac2");
  gr_frac2->Write();

  gr_frac3->SetName("gr_frac3");
  gr_frac3->Write();

  gr_ww->SetName("gr_ww");
  gr_ww->Write();

  gr_fw->SetName("gr_fw");
  gr_fw->Write();

  gr_sw->SetName("gr_sw");
  gr_sw->Write();

  gr_tw->SetName("gr_tw");
  gr_tw->Write();

  f->Write();
   f->Close();
   delete f;
   in.close();

}
