#include "Riostream.h"

void rootgen_dcball_par() {
  ifstream in;
  in.open(Form("dcball_par_all.dat"));

  Float_t x,z,ez,f1,ef1,f2,ef2, fw, fw_er, fal, fal_er, fn, fn_er, sw, sw_er, sal, sal_er, sn, sn_er;

  double_t m[54];
  Double_t ex[54] = {0};

  string line;

  double_t ww[54];
  double_t ww_er[54];

  double_t db_f1[54];
  double_t e_db_f1[54];

  double_t db_f2[54];
  double_t e_db_f2[54];

  double_t w1[54];
  double_t w1_er[54];

  double_t al1[54];
  double_t al1_er[54];

  double_t n1[54];
  double_t n1_er[54];

  double_t w2[54];
  double_t w2_er[54];

  double_t al2[54];
  double_t al2_er[54];

  double_t n2[54];
 double_t n2_er[54];

  Int_t nlines = 0;
   TFile *f = new TFile("dcball_par.root","RECREATE");

   while (nlines < 53) {
     in >> x >> z >> ez >> f1 >> ef1 >> f2 >> ef2 >> fw >> fw_er >> fal >> fal_er >> fn >> fn_er >> sw >> sw_er >> sal >> sal_er >> sn >> sn_er;
     // if (!in.good()) break;
     if (nlines < 20) printf(" x=%8f, z=%8f, ez=%8f, f1=%8f, ef1=%8f, f2=%8f, ef2=%8f",x, z, ez, f1, ef1, f2, ef2);
     m[nlines] = x;
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
     cout << " number of lines " << nlines << endl;
     nlines++;
   }
  printf(" found %d points\n",nlines);

  TGraphErrors * gr_dbfrac1 = new TGraphErrors(nlines,m,db_f1,ex,e_db_f1);
  gr_dbfrac1->SetTitle("double crystal ball fit function normalized fraction 1 parameterized;M[GeV/c^{2}];Fraction;");

  TGraphErrors * gr_dbfrac2 = new TGraphErrors(nlines,m,db_f2,ex,e_db_f2);
  gr_dbfrac2->SetTitle("double crystal ball fit function normalized fraction 2 parameterized;M[GeV/c^{2}];Fraction;");

  TGraphErrors * gr_ww = new TGraphErrors(nlines,m,ww,ex,ww_er);
  gr_ww->SetTitle("double crystal ball fit function weighted width parametrized;M[GeV/c^{2}];Weighted Width[GeV/c^{2}];");

  TGraphErrors * gr_fw = new TGraphErrors(nlines,m,w1,ex,w1_er);
  gr_fw->SetTitle("1st Width from double crystalball ff;M[GeV/c^{2}];#sigma_{1}");

  TGraphErrors * gr_sw = new TGraphErrors(nlines,m,w2,ex,w2_er);
  gr_sw->SetTitle("2nd Width from double crystalball ff;M[GeV/c^{2}];#sigma_{2}");

  TGraphErrors * gr_fal = new TGraphErrors(nlines,m,al1,ex,al1_er);
  gr_fal->SetTitle("1st alpha from double crystalball ff;M[GeV/c^{2}];#alpha_{1}");

  TGraphErrors * gr_sal = new TGraphErrors(nlines,m,al2,ex,al2_er);
  gr_sal->SetTitle("2nd alpha from double crystalball ff;M[GeV/c^{2}];#alpha_{2}");

  TGraphErrors * gr_fn = new TGraphErrors(nlines,m,n1,ex,n1_er);
  gr_fn->SetTitle("1st n from double crystalball ff;M[GeV/c^{2]}];N_{1}");

  TGraphErrors * gr_sn = new TGraphErrors(nlines,m,n2,ex,n2_er);
  gr_sn->SetTitle("2nd n from double crystalball ff;M[GeV/c^{2]}];N_{2}");

  gr_dbfrac1->SetName("gr_dbfrac1");
  gr_dbfrac1->Write();

  gr_dbfrac2->SetName("gr_dbfrac2");
  gr_dbfrac2->Write();

  gr_ww->SetName("gr_ww");
  gr_ww->Write();

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
