#include "TCanvas.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TH1F.h"
#include <string>
#include "TGraphAsymmErrors.h"

#ifndef SimpleCanvas_h
#define SimpleCanvas_h
// --- Class SimpleCanvas
class SimpleCanvas 
{
  public:
  SimpleCanvas(string, int nPanels = 1);
  TCanvas *can_;
  TPad* can_up_;
  TPad* can_dw_;
  TPaveText *title_; 
  TPaveText *pt1_; 
  TPaveText *pt2_; 
  TPaveText *pt3_; 
  TPaveText *st1_; 
  TPaveText *st2_; 
  TPaveText *st3_; 
  TPaveText *phys14_1_; 
  TPaveText *phys14_2_; 
  TPaveText *phys14_3_; 

  int nPanels_;
  string canName_;
  void Up(){if (nPanels_==2)can_up_->cd();}
  void Dw(){if (nPanels_==2)can_dw_->cd();}
  void SetLogy(){if (nPanels_==2)can_up_->SetLogy(); if(nPanels_ == 1) can_->SetLogy();}
  void SetLogz(){if (nPanels_==2)can_up_->SetLogz(); if(nPanels_ == 1) can_->SetLogz();}
  void Save()
  {
    can_->SaveAs((canName_+".pdf").c_str());
    can_->SaveAs((canName_+".png").c_str());
    can_->SaveAs((canName_+".C").c_str());
  }
  void Save(string outputDir)
  {
    can_->SaveAs((outputDir+canName_+".pdf").c_str());
    can_->SaveAs((outputDir+canName_+".png").c_str());
    can_->SaveAs((outputDir+canName_+".C").c_str());
  }
  void ShapeMeDw(TH1F*);
  void ShapeMeDw(TGraphAsymmErrors*);
  void ShapeMeUp(TH1F*);
  void ShapeMe(TH1*);
  void CMSPre(){can_up_->cd();pt1_->Draw("same");pt2_->Draw("same");pt3_->Draw("same");}
  void CMSSim(){can_up_->cd();st1_->Draw("same");st2_->Draw("same");}
  void CMSPhys14(){ if(nPanels_==2)can_up_->cd(); phys14_1_->Draw("same");phys14_2_->Draw("same");phys14_3_->Draw("same");}
  void SynchUpDw(TH1F*,TGraphAsymmErrors*);
};

void SimpleCanvas::SynchUpDw(TH1F* h1,TGraphAsymmErrors* h2)
{

  float xMin = h1->GetXaxis()->GetXmin();
  float xMax = h1->GetXaxis()->GetXmax();
  h2->GetXaxis()->SetLimits(xMin, xMax);
  h2->GetXaxis()->SetRangeUser(xMin, xMax);
}
void SimpleCanvas::ShapeMe(TH1*hist)
{

//  hist->GetXaxis()->SetTitleFont(43);
//  hist->GetXaxis()->SetTitleSize(28);
//  hist->GetXaxis()->SetLabelFont(43);
//  hist->GetXaxis()->SetLabelSize(26);
//  hist->GetXaxis()->SetTitleOffset(1.0);
//  hist->GetYaxis()->SetTitleFont(43);
//  hist->GetYaxis()->SetTitleSize(28);
//  hist->GetYaxis()->SetLabelFont(43);
//  hist->GetYaxis()->SetLabelSize(26);
  hist->GetYaxis()->SetTitleOffset(1.13);

}

void SimpleCanvas::ShapeMeDw(TH1F*hist)
{

  hist->GetXaxis()->SetTitleFont(43);
  hist->GetXaxis()->SetTitleSize(28);
  hist->GetXaxis()->SetLabelFont(43);
  hist->GetXaxis()->SetLabelSize(26);
  hist->GetXaxis()->SetTitleOffset(3.0);
  hist->GetYaxis()->SetTitleFont(43);
  hist->GetYaxis()->SetTitleSize(28);
  hist->GetYaxis()->SetLabelFont(43);
  hist->GetYaxis()->SetLabelSize(26);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetYaxis()->SetNdivisions(509);

}

void SimpleCanvas::ShapeMeDw(TGraphAsymmErrors *hist)
{
  hist->GetXaxis()->SetTitleFont(43);
  hist->GetXaxis()->SetTitleSize(28);
  hist->GetXaxis()->SetLabelFont(43);
  hist->GetXaxis()->SetLabelSize(26);
  hist->GetXaxis()->SetTitleOffset(3.0);
  hist->GetYaxis()->SetTitleFont(43);
  hist->GetYaxis()->SetTitleSize(28);
  hist->GetYaxis()->SetLabelFont(43);
  hist->GetYaxis()->SetLabelSize(26);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetYaxis()->SetNdivisions(509);
}


void SimpleCanvas::ShapeMeUp(TH1F*hist)
{
  hist->GetXaxis()->SetTitleFont(43);
  hist->GetXaxis()->SetTitleSize(28);
  hist->GetXaxis()->SetLabelFont(43);
  hist->GetXaxis()->SetLabelSize(26);
  hist->GetXaxis()->SetLabelOffset(0.025);
  hist->GetXaxis()->SetTitleOffset(1.35);
  hist->GetYaxis()->SetTitleFont(43);
  hist->GetYaxis()->SetTitleSize(28);
  hist->GetYaxis()->SetLabelFont(43);
  hist->GetYaxis()->SetLabelSize(26);
  hist->GetYaxis()->SetTitleOffset(1.45);
}



SimpleCanvas::SimpleCanvas(string canName, int nPanels)
{
  canName_ = canName;
  nPanels_ = nPanels;

  title_ = new TPaveText(0.20, 0.95, 0.80, 1.0,"blNDC");
  title_->SetFillStyle(0);
  title_->SetFillColor(kWhite);
  title_->SetBorderSize(0);
  //title_->AddText("CMS Preliminary, 19.4 fb^{-1} (8 TeV)"); // README put correct header
//  title_->AddText("ATLAS Preliminary, 20.3 fb^{-1} (8 TeV)"); // README put correct header
  title_->AddText(""); // README put correct header

  // --- define pt1,2,3
  pt1_ = new TPaveText(0.129,0.921,0.258,1,"blNDC");
  pt1_->SetBorderSize(0);
  pt1_->SetFillColor(0);
  pt1_->SetFillStyle(0);
  pt1_->SetTextFont(63);
  pt1_->SetTextSize(30);
  pt1_->AddText("CMS");

  pt2_ = new TPaveText(0.26,0.91,0.51,1,"blNDC");
  pt2_->SetBorderSize(0);
  pt2_->SetFillColor(0);
  pt2_->SetFillStyle(0);
  pt2_->SetTextFont(53);
  pt2_->SetTextSize(25);
  pt2_->AddText("Preliminary");

  pt3_ = new TPaveText(0.60,0.921,0.95,1,"blNDC");
  pt3_->SetBorderSize(0);
  pt3_->SetFillColor(0);
  pt3_->SetFillStyle(0);
  pt3_->SetTextFont(43);
  pt3_->SetTextSize(25);
  pt3_->AddText("19.4 fb^{-1} (8 TeV)");

  // --- define st1,2,3
  st1_ = new TPaveText(0.129,0.921,0.258,1,"blNDC");
  st1_->SetBorderSize(0);
  st1_->SetFillColor(0);
  st1_->SetFillStyle(0);
  st1_->SetTextFont(63);
  st1_->SetTextSize(30);
  st1_->AddText("CMS");

  st2_ = new TPaveText(0.26,0.91,0.51,1,"blNDC");
  st2_->SetBorderSize(0);
  st2_->SetFillColor(0);
  st2_->SetFillStyle(0);
  st2_->SetTextFont(53);
  st2_->SetTextSize(25);
  st2_->AddText("Simulation");

  st3_ = new TPaveText(0.60,0.921,0.95,1,"blNDC");
  st3_->SetBorderSize(0);
  st3_->SetFillColor(0);
  st3_->SetFillStyle(0);
  st3_->SetTextFont(43);
  st3_->SetTextSize(25);
  st3_->AddText("(8 TeV)");

  // --- define phys14_1,2,3
  float pushToRight1 = 0.0;
  float pushToRight2 = 0.03;
  phys14_1_ = new TPaveText(0.129+pushToRight1,0.921,0.258+pushToRight1,1,"blNDC");
  phys14_1_->SetBorderSize(0);
  phys14_1_->SetFillColor(0);
  phys14_1_->SetFillStyle(0);
  phys14_1_->SetTextFont(63);
  phys14_1_->SetTextSize(30);
  phys14_1_->AddText("CMS");

  phys14_2_ = new TPaveText(0.26+pushToRight1,0.93,0.51+pushToRight1,1,"blNDC");
  phys14_2_->SetBorderSize(0);
  phys14_2_->SetFillColor(0);
  phys14_2_->SetFillStyle(0);
  phys14_2_->SetTextFont(53);
  phys14_2_->SetTextSize(25);
  phys14_2_->AddText("Simulation");

  phys14_3_ = new TPaveText(0.60+pushToRight2,0.921,0.95+pushToRight2,1,"blNDC");
  phys14_3_->SetBorderSize(0);
  phys14_3_->SetFillColor(0);
  phys14_3_->SetFillStyle(0);
  phys14_3_->SetTextFont(43);
  phys14_3_->SetTextSize(25);
  phys14_3_->AddText("1 fb^{-1} (13 TeV)");

  if(nPanels_ == 2)
  {
    can_ = new TCanvas(canName_.c_str(), canName_.c_str(), 500,650);
    can_->Divide(1, 2);
    can_up_ = (TPad*) can_->GetListOfPrimitives()->FindObject((canName_+"_1").c_str());
    can_dw_ = (TPad*) can_->GetListOfPrimitives()->FindObject((canName_+"_2").c_str());
    can_up_->SetPad(0., 0.35, 1., 1.);
    can_dw_->SetPad(0., 0., 1., 0.35);
    can_up_->SetFrameFillColor(0);
    can_up_->SetFillColor(0);
    can_dw_->SetFillColor(0);
    can_dw_->SetFrameFillColor(0);
    can_up_->SetTopMargin(0.08);
    can_up_->SetBottomMargin(0.03);
    can_dw_->SetBottomMargin(0.3);
    can_dw_->SetGridy();
    can_dw_->SetGridx();
//    can_up_->SetLeftMargin(0.14);
//    can_dw_->SetLeftMargin(0.14);
    can_up_->cd();
    title_->Draw("");
  }

  if(nPanels_ == 1)
  {
    can_ = new TCanvas(canName_.c_str(), canName_.c_str(), 500, 500);
    can_->SetTopMargin(0.08);
    can_->SetLeftMargin(0.14);
    title_->Draw("");
  }
}
// --- end: Class SimpleCanvas
#endif 
