#include "Riostream.h"

void rootgen_biast() {
  ifstream in;
  in.open(Form("./bias_test.dat"));

  Float_t x,z,y;

  double_t m[21];
  Double_t ex[21] = {0};

  double_t ww[21];
  double_t ww_er[21];

  Double_t deteff[21];
  Double_t deteffer[21];

  Int_t nlines = 0;
   TFile *f = new TFile("biasmap.root","RECREATE");
   TNtuple *ntuple = new TNtuple("ntuple","data from parameters_all","x:z:y");

   while (1) {
     in >> x >> z >> y;
     if (!in.good()) break;
     if (nlines < 21) printf(" x=%8f, z=%8f, y=%8f",x, z, y);
     ntuple->Fill(x,z,y);
     m[nlines] = x;
     deteff[nlines] = y;
     ww[nlines] = z;
     nlines++;
   }
     printf(" found %d points\n",nlines);

  TGraph * gr_deteff = new TGraph(nlines,m,deteff);
  gr_deteff->SetTitle(" pull dist fit sigma as function of the Z' mass;M[GeV/c^{2}];pull fit mean;");

  TGraph * gr_ww = new TGraph(nlines,m,ww);
  gr_ww->SetTitle("pull dist fit mean as function of Z' mass;M[GeV/c^{2}];pull fit sigma;");


  gr_deteff->SetName("gr_deteff");
  gr_deteff->Write();

  gr_ww->SetName("gr_ww");
  gr_ww->Write();


  f->Write();
   f->Close();
   delete f;
   in.close();

}
