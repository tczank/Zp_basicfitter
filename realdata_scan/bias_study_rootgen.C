#include "Riostream.h"

void bias_study_rootgen() {
  ifstream in;
  in.open(Form("./bias_study_all.dat"));

  Float_t x,y,o,xs,s;

  double_t m[2219];

  double_t toymc_sigma[2219];
  double_t toymc_mean[2219];

  Int_t nlines = 0;
   TFile *f = new TFile("bias_study_ww.root","RECREATE");
   TNtuple *ntuple = new TNtuple("ntuple","data from parameters_all","x:y:ey:o:xs:s");

   while (1) {
     in >> x >> y >> o ;
     if (!in.good()) break;
     // if (nlines < 19) printf(" x=%8f, z=%8f, ez=%8f, f1=%8f, ef1=%8f, f2=%8f, ef2=%8f, y=%8f, ey=%8f, o=%8f, xs=%8f, s=%8f",x, z, ez, f1,ef1, f2, ef2, y, ey, o, xs, s);
     ntuple->Fill(x,y,o);
//     m[nlines] = sqrt(pow(x,2) -4*pow(0.1056583745,2));
     m[nlines] = x;
     toymc_mean[nlines] = y;
     toymc_sigma[nlines] = o;
     nlines++;
   }
  printf(" found %d points\n",nlines);

 
  TGraph * gr_toymc_gm = new TGraph(nlines,m,toymc_mean);
  gr_toymc_gm->SetTitle("mean of the gaussian from toy mc");

  TGraph * gr_toymc_sigma = new TGraph(nlines,m,toymc_sigma);
  gr_toymc_sigma->SetTitle("sigma of the gaussian from toy mc");

   gr_toymc_gm->SetName("gr_toymc_gm");
  gr_toymc_gm->Write();

  gr_toymc_sigma->SetName("gr_toymc_sigma");
  gr_toymc_sigma->Write();

   f->Write();
   f->Close();
   delete f;
   in.close();

}
