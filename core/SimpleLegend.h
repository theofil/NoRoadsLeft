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
  this->SetTextFont(63);
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
    this->SetX1NDC(0.65);
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
      if(previousHist->GetFillColor() == histPointer->GetFillColor())addEntry = false;
    }
  
    if(addEntry && sampleTitle != "Data") this->AddEntry(histPointer, sampleTitle.c_str(), "fl");
    if(addEntry && sampleTitle == "Data") this->AddEntry(histPointer, sampleTitle.c_str(), "pe");
  }
}
#endif