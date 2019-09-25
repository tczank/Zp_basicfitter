#include "Riostream.h"

void rootgen_pioff_pars() {
  ifstream in;
  in.open(Form("pion_eff_all.dat"));

  Float_t x,z,ez,f1,ef1,f2,ef2,y,ey,o,xs,s, fw, fw_er, fal, fal_er, fn, fn_er, sw, sw_er, sal, sal_er, sn, sn_er;

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

  double_t w1[21];
  double_t w1_er[21];

  double_t al1[21];
  double_t al1_er[21];

  double_t n1[21];
  double_t n1_er[21];

  double_t w2[21];
  double_t w2_er[21];

  double_t al2[21];
  double_t al2_er[21];

  double_t n2[21];
 double_t n2_er[21];

  Int_t nlines = 0;
   TFile *f = new TFile("pioneff_all.root","RECREATE");
   TNtuple *ntuple = new TNtuple("ntuple","data from parameters_all","x:z:ez:f1:ef1:f2:ef2:y:ey:o:xs:s");

   while (1) {
     in >> x >> z >> ez >> f1 >> ef1 >> f2 >> ef2 >> y >> ey >> o >> xs >> s >> fw >> fw_er >> fal >> fal_er >> fn >> fn_er >> sw >> sw_er >> sal >> sal_er >> sn >> sn_er;
     if (!in.good()) break;
     if (nlines < 21) printf(" x=%8f, z=%8f, ez=%8f, f1=%8f, ef1=%8f, f2=%8f, ef2=%8f, y=%8f, ey=%8f, o=%8f, xs=%8f, s=%8f",x, z, ez, f1,ef1, f2, ef2, y, ey, o, xs, s);
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
     if (fw < sw){
       w1[nlines] = fw;
       w1_er[nlines] = fw_er;
       w2[nlines] = sw;
       w2_er[nlines] = sw_er;
     }
     if (fw > sw){
       w1[nlines] = sw;
       w1_er[nlines] = sw_er;
       w2[nlines] = fw;
       w2_er[nlines] = fw_er;
     }
     if( fal < sal){
       al1[nlines] = fal;
       al1_er[nlines] = fal_er;
       al2[nlines] = sal;
       al2_er[nlines] = sal_er;
     }
     if ( fal > sal){
       al1[nlines] = sal;
       al1_er[nlines] = sal_er;
       al2[nlines] = fal;
       al2_er[nlines] = fal_er;
     }
     if ( fn < sn){
       n1[nlines] = fn;
       n1_er[nlines] = fn_er;
       n2[nlines] = sn;
       n2_er[nlines] = sn_er;
     }
     if ( fn > sn){
       n1[nlines] = sn;
       n1_er[nlines] = sn_er;
       n2[nlines] = fn;
       n2_er[nlines] = fn_er;
     }
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

  TGraphErrors * gr_fw = new TGraphErrors(nlines,m,w1,ex,w1_er);
  gr_fw->SetTitle("1st Width from double gaussian ff;M[GeV/c^{2}];#sigma_{1}");

  TGraphErrors * gr_sw = new TGraphErrors(nlines,m,w2,ex,w2_er);
  gr_sw->SetTitle("2nd Width from double gaussian ff;M[GeV/c^{2}];#sigma_{2}");

  TGraphErrors * gr_fal = new TGraphErrors(nlines,m,al1,ex,al1_er);
  gr_fal->SetTitle("1st alpha from double gaussian ff;M[GeV/c^{2}];#alpha_{1}");

  TGraphErrors * gr_sal = new TGraphErrors(nlines,m,al2,ex,al2_er);
  gr_sal->SetTitle("2nd alpha from double gaussian ff;M[GeV/c^{2}];#alpha_{2}");

  TGraphErrors * gr_fn = new TGraphErrors(nlines,m,n1,ex,n1_er);
  gr_fn->SetTitle("1st n from double gaussian ff;M[GeV/c^{2]}];N_{1}");

  TGraphErrors * gr_sn = new TGraphErrors(nlines,m,n2,ex,n2_er);
  gr_sn->SetTitle("2nd n from double gaussian ff;M[GeV/c^{2]}];N_{2}");

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

  gr_fw->SetName("gr_fw");
  gr_fw->Write();

  gr_sw->SetName("gr_sw");
  gr_sw->Write();

  gr_fal->SetName("gr_fal");
  gr_fal->Write();

  gr_sal->SetName("gr_sal");
  gr_sal->Write();

  gr_fn->SetName("gr_fn");
  gr_fn->Write();

  gr_sn->SetName("gr_sn");
  gr_sn->Write();

  f->Write();
   f->Close();
   delete f;
   in.close();

}
