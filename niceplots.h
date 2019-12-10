#include "Nice2D.c"
#include "Nice3D2D.c"
#include "Nice1D.c"


void niceplots{
gROOT->Reset();

  gROOT->SetStyle("Bold");
  gStyle->SetCanvasColor(0);
  gStyle->SetLabelColor(1);
  gStyle->SetLabelColor(1,"Y");
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineWidth(1);
  gStyle->SetNdivisions(505);
  gStyle->SetNdivisions(505,"Y");
  //gROOT->Macro("setcolor2.c");
  gStyle->SetHistFillColor(999);
  gROOT->SetStyle("Plain");  // white as bg
  gStyle->SetOptStat("111111");
  gStyle->SetFrameLineWidth(1);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleBorderSize(0);

  //Definitions
  Double_t smallBetween1 = .1;
  Double_t smallBetween2 = .1;
  Double_t smallBetween3 = .1;
  Double_t smallBetween4 = .1;

  Double_t small = .00001;

  TLine TLine;
  TLatex *t = new TLatex();
  t->SetTextSize(0.03);
  t->SetTextFont(42);
  t->SetTextAlign(12);
  t->SetTextColor(1);
  t->SetTextFont(12);

  TCanvas *C1;


 smallBetween1 = .125;
  smallBetween2 = .05;
  smallBetween3 = .05;
  smallBetween4 = .125;
  TString cleg = "MyPicture";
    C1 = new TCanvas(cleg, cleg, 10, 10, 800, 800);
    gPad->SetLeftMargin(smallBetween1);
    gPad->SetRightMargin(smallBetween2);
    gPad->SetTopMargin(smallBetween3);
    gPad->SetBottomMargin(smallBetween4);
}
