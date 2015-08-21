#include "TLegend.h"
#include "TLegendEntry.h"
#include "THStack.h"

#include "SimpleDriver.h"
#include "SimpleStack.h"

#ifndef SimpleLegend_h
#define SimpleLegend_h
class SimpleLegend:public TLegend
{
  public:
  SimpleLegend();
  SimpleLegend(string position);
  void FillLegend(SimpleStack *);
  void SetStyle();
};


void SimpleLegend::SetStyle()
{
  this->SetTextSize(24);
  this->SetTextFont(43);
  this->SetFillColor(kWhite);
  this->SetFillStyle(0);
  this->SetBorderSize(0);
}

SimpleLegend::SimpleLegend()
{
  this->SetX1NDC(0.65);
  this->SetY1NDC(0.47);
  this->SetX2NDC(0.95);
  this->SetY2NDC(0.89);
  SetStyle();
}
SimpleLegend::SimpleLegend(string position)
{
  if(position == "TR")
  {
    this->SetX1NDC(0.75);
    this->SetY1NDC(0.47);
    this->SetX2NDC(0.95);
    this->SetY2NDC(0.89);
  }
  if(position == "TL") // 0.15, 0.5, 0.5, 0.9
  {
    this->SetX1NDC(0.15);
    this->SetY1NDC(0.5);
    this->SetX2NDC(0.5);
    this->SetY2NDC(0.9);
  }

  if(position == "TLSF")
  {
    this->SetX1NDC(0.15);
    this->SetY1NDC(0.48);
    this->SetX2NDC(0.50);
    this->SetY2NDC(0.90);
  }

  if(position == "TLOF")
  {
    this->SetX1NDC(0.15);
    this->SetY1NDC(0.48);
    this->SetX2NDC(0.50);
    this->SetY2NDC(0.90);
  }


  if(position == "TRSF")
  {
    this->SetX1NDC(0.63);
    this->SetY1NDC(0.48);
    this->SetX2NDC(0.95);
    this->SetY2NDC(0.90);
  }

  if(position == "TROF")
  {
    this->SetX1NDC(0.60);
    this->SetY1NDC(0.48);
    this->SetX2NDC(0.92);
    this->SetY2NDC(0.90);
  }

  if(position == "TRSFjzb_pos")
  {
    this->SetX1NDC(0.63);
    this->SetY1NDC(0.41);
    this->SetX2NDC(0.95);
    this->SetY2NDC(0.90);
  }

  if(position == "TROFjzb_pos")
  {
    this->SetX1NDC(0.60);
    this->SetY1NDC(0.39);
    this->SetX2NDC(0.91);
    this->SetY2NDC(0.88);
  }

  if(position == "TRSFjzb_neg")
  {
    this->SetX1NDC(0.63);
    this->SetY1NDC(0.41);
    this->SetX2NDC(0.95);
    this->SetY2NDC(0.90);
  }

  if(position == "TROFjzb_neg")
  {
    this->SetX1NDC(0.60);
    this->SetY1NDC(0.39);
    this->SetX2NDC(0.91);
    this->SetY2NDC(0.88);
  }

  if(position == "jzb_exp_mode")
  {
    this->SetX1NDC(0.53);
    this->SetY1NDC(0.68);
    this->SetX2NDC(0.95);
    this->SetY2NDC(0.90);
  }

  if(position == "TTR")
  {
    this->SetX1NDC(0.58);
    this->SetY1NDC(0.77);
    this->SetX2NDC(0.85);
    this->SetY2NDC(0.92);
  }

  SetStyle();
}

void SimpleLegend::FillLegend(SimpleStack *sstack)
{
  int nSamples = (int)sstack->sampleTitles_.size();
  for(int sCounter = nSamples - 1; sCounter >= 0; --sCounter)
  {
    string sampleTitle = sstack->sampleTitles_[sCounter];
    TH1F *histPointer  = sstack->histoPointers_[sCounter];

    bool addEntry = true;        

    if(sCounter < nSamples - 1)  // check if histogram is not to be visible in the TLegend (has been merged with previous having the same fill color)
    {
      TH1F *previousHist = sstack->histoPointers_[sCounter+1];
      if(previousHist->GetFillColor() == histPointer->GetFillColor()){addEntry = false; cout <<"sample " << sampleTitle << " has the same fillcolor as previous one, will not adding it in the TLegend" << endl;  }
    }
  
    if(addEntry && sampleTitle != "Data") this->AddEntry(histPointer, sampleTitle.c_str(), "fl");
    if(addEntry && sampleTitle == "Data") this->AddEntry(histPointer, sampleTitle.c_str(), "pe");
  }
}
#endif
