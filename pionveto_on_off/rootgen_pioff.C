#include "Riostream.h"

void rootgen_pioff() {
  ifstream in;
  in.open(Form("./pioff_eff_all.dat"));

  Float_t x,z,ez,f1,ef1,f2,ef2,y,ey,o,xs,s;

  double_t m[21];
  Double_t ex[21] = {0};

  double_t ww[21];
  double_t ww_er[21];

  double_t db_f1[21];
  double_t e_db_f1[21];

  double_t db_f2[21];
  double_t e_db_f2[21];

  Double_t deteff[21];
  Double_t deteffer[21];

  double_t nobs[21];
  double_t cross[21];
  double_t significance[21];

  Int_t nlines = 0;
   TFile *f = new TFile("pioffeff_all.root","RECREATE");
   TNtuple *ntuple = new TNtuple("ntuple","data from parameters_all","x:z:ez:f1:ef1:f2:ef2:y:ey:o:xs:s");

   while (1) {
     in >> x >> z >> ez >> f1 >> ef1 >> f2 >> ef2 >> y >> ey >> o >> xs >> s;
     if (!in.good()) break;
     if (nlines < 21) printf(" x=%8f, z=%8f, ez=%8f, f1=%8f, ef1=%8f, f2=%8f, ef2=%8f, y=%8f, ey=%8f, o=%8f, xs=%8f, s=%8f \n ",x, z, ez, f1,ef1, f2, ef2, y, ey, o, xs, s);
     ntuple->Fill(x,z,ez,f1,ef1,f2,ef2,y,ey,o,xs,s);
     m[nlines] = x;
     deteff[nlines] = y;
     deteffer[nlines] = ey;
     ww[nlines] = z;
     ww_er[nlines] = ez;
     if(f1 < f2){
       db_f1[nlines] = f1;
       e_db_f1[nlines] =ef1;
       db_f2[nlines] =f2;
       e_db_f2[nlines] =ef2;
     }
     if(f1 > f2){
     db_f2[nlines] =f1;
     e_db_f2[nlines] =ef1;
     db_f1[nlines] =f2;
     e_db_f1[nlines] =ef2;}
     nobs[nlines] = o;
     cross[nlines] = xs;
     significance[nlines] = s;
     nlines++;
   }
    printf(" found %d points\n",nlines);

  TGraphErrors * gr_deteff = new TGraphErrors(nlines,m,deteff,ex,deteffer);
  gr_deteff->SetTitle(" detection efficiency as function of the Z' mass;M[GeV/c^{2}];Detection Efficiency;");

  TGraphErrors * gr_dbfrac1 = new TGraphErrors(nlines,m,db_f1,ex,e_db_f1);
  gr_dbfrac1->SetTitle("double crystal ball fit function normalized fraction 1 parameterized;M[GeV/c^{2}];Fraction;");

  TGraphErrors * gr_dbfrac2 = new TGraphErrors(nlines,m,db_f2,ex,e_db_f2);
  gr_dbfrac2->SetTitle("double crystal ball fit function normalized fraction 2 parameterized;M[GeV/c^{2}];Fraction;");

  TGraphErrors * gr_ww = new TGraphErrors(nlines,m,ww,ex,ww_er);
  gr_ww->SetTitle("double crystal ball fit function weighted width parametrized;M[GeV/c^{2}];Weighted Width[GeV/c^{2}];");

  TGraph * gr_obs = new TGraph(nlines,m,nobs);
  gr_obs->SetTitle("number of observed events as function of the Z' mass;M[GeV/c^{2}];Nobs;");

  TGraph * gr_xs = new TGraph(nlines,m,cross);
  gr_xs->SetTitle("90% CL cross section for Z' based on MC background as function of the Z' mass;M[GeV/c^{2}];#sigma[ab]");

  TGraph * gr_s = new TGraph(nlines,m,significance);
  gr_s->SetTitle("Significance of the Fit double crystal ball compared with the background pol3; M[GeV/c^{2}];Significance");

  gr_deteff->SetName("gr_deteff");
  gr_deteff->Write();

  gr_dbfrac1->SetName("gr_dbfrac1");
  gr_dbfrac1->Write();

  gr_dbfrac2->SetName("gr_dbfrac2");
  gr_dbfrac2->Write();

  gr_ww->SetName("gr_ww");
  gr_ww->Write();

  gr_obs->SetName("gr_obs");
  gr_obs->Write();

  gr_xs->SetName("gr_xs");
  gr_xs->Write();

  gr_s->SetName("gr_s");
  gr_s->Write();

  f->Write();
   f->Close();
   delete f;
   in.close();

}
