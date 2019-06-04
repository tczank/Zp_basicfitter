#include "Riostream.h"

void rootgen_isr_realdat_v0_Up1s() {
  ifstream in;
  in.open(Form("./real_onres_xs_all.dat"));

  Float_t x,y,o,xs,s;

  double_t m[20000];
  Double_t ex[20000] = {0};

  Double_t deteff[20000];

  double_t nobs[20000];
  double_t cross[20000];
  double_t significance[20000];

  Int_t nlines = 0;
   TFile *f = new TFile("real_onres_xs_all.root","RECREATE");
   TNtuple *ntuple = new TNtuple("ntuple","data from parameters_all","x:y:ey:o:xs:s");

   while (1) {
     in >> x >> y >> o >> xs >> s;
     if (!in.good()) break;
     // if (nlines < 19) printf(" x=%8f, z=%8f, ez=%8f, f1=%8f, ef1=%8f, f2=%8f, ef2=%8f, y=%8f, ey=%8f, o=%8f, xs=%8f, s=%8f",x, z, ez, f1,ef1, f2, ef2, y, ey, o, xs, s);
     ntuple->Fill(x,y,o,xs,s);
     m[nlines] = x;
     deteff[nlines] = y;
     nobs[nlines] = o;
     cross[nlines] = xs;
     significance[nlines] = s;
     nlines++;
   }
  printf(" found %d points\n",nlines);


  TGraph * gr_obs = new TGraph(nlines,m,nobs);
  gr_obs->SetTitle("number of observed events as function of the Z' mass;M[GeV/c^{2}];Nobs;");

  TGraph * gr_xs = new TGraph(nlines,m,cross);
  gr_xs->SetTitle("90% CL cross section for Z' based on MC background as function of the Z' mass;M[GeV/c^{2}];#sigma[ab]");

  TGraph * gr_s = new TGraph(nlines,m,significance);
  gr_s->SetTitle("Significance of the Fit double crystal ball compared with the background pol3; M[GeV/c^{2}];Significance");

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
