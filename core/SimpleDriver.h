#include <functional>
#include <string>
#include "any2string.h"
#include "SimpleSample.h"
#include "SimpleStack.h"

#include "TH1F.h"
#include "THStack.h"
#include "TCut.h"

#ifndef SimpleDriver_h
#define SimpleDriver_h


class SimpleDriver: public std::vector<SimpleSample*>
{
  public: 
  SimpleDriver(){hcounter_ = 0; scounter_ = 0;}

  SimpleStack * getSimpleStackTH1F(string var, string histTitle, int nBins, float xMin, float xMax, TCut myCut)
  {
    string stackName = var + "_s" + any2string(scounter_);  
    simple_stacks_[scounter_] = new SimpleStack(stackName.c_str(), stackName.c_str());
    scounter_++;

    for(size_t sampleIt = 0; sampleIt < this->size(); ++sampleIt)
    {
     string histName = var + "_h" + any2string(hcounter_) + "_" + getHash();
     histos_[hcounter_] = new TH1F(histName.c_str(), histTitle.c_str(), nBins, xMin, xMax);
     histos_[hcounter_]->Sumw2();
     (*this)[sampleIt]->SetHistoStyle(histos_[hcounter_]); 
     hcounter_++;
     string drawCommand = var + ">>" + histName;
     (*this)[sampleIt]-> events_->Draw(drawCommand.c_str(),  myCut*(*this)[sampleIt]->tcut_, "goff", (*this)[sampleIt]->maxN_);

     moveOverflowToLastBin ( histos_[hcounter_-1] );
     moveUnderflowToLastBin( histos_[hcounter_-1] );

     simple_stacks_[scounter_-1]->Add(histos_[hcounter_-1]);
     simple_stacks_[scounter_-1]->sampleTitles_.push_back( (*this)[sampleIt]->title_);
     simple_stacks_[scounter_-1]->histoPointers_.push_back( histos_[hcounter_-1]);
    }
    return  simple_stacks_[scounter_-1];
  }

  THStack * getStackTH1F(string var, string histTitle, int nBins, float xMin, float xMax, TCut myCut)
  {
    string stackName = var + "_s" + any2string(scounter_);  
    stacks_[scounter_] = new THStack(stackName.c_str(), stackName.c_str());
    scounter_++;

    for(size_t sampleIt = 0; sampleIt < this->size(); ++sampleIt)
    {
     string histName = var + "_h" + any2string(hcounter_) + "_" + getHash();
     histos_[hcounter_] = new TH1F(histName.c_str(), histTitle.c_str(), nBins, xMin, xMax);
     histos_[hcounter_]->Sumw2();
     (*this)[sampleIt]->SetHistoStyle(histos_[hcounter_]); 
     hcounter_++;
     string drawCommand = var + ">>" + histName;
     (*this)[sampleIt]-> events_->Draw(drawCommand.c_str(),  myCut*(*this)[sampleIt]->tcut_, "goff", (*this)[sampleIt]->maxN_);

     moveOverflowToLastBin ( histos_[hcounter_-1] );
     moveUnderflowToLastBin( histos_[hcounter_-1] );

     stacks_[scounter_-1]->Add(histos_[hcounter_-1]);
    }
    return  stacks_[scounter_-1];
  }

  TH1F *getHistoTH1F(THStack *hs)
  {
    TList *histos = hs->GetHists();
    TH1F  *sum = (TH1F*) hs->GetStack()->First()->Clone("");
    sum->Reset();
    sum->Merge(histos);
    sum->SetFillColor(kWhite);
    sum->SetLineColor(kBlack);
    return sum;
  }

  TH1F *getHistoTH1F(string var, string histTitle, int nBins, float xMin, float xMax, TCut myCut) // uses getHistoTH1F(THStack *hs) but first creates the stack
  {
    THStack *hs_temp = getStackTH1F(var, histTitle, nBins, xMin, xMax, myCut);
    return getHistoTH1F(hs_temp);
  }

  void moveOverflowToLastBin(TH1* h1)
  {
      Int_t nBins = h1->GetNbinsX();
      Int_t overflowBin = nBins + 1;
      Int_t lastBin = nBins;


      float overflowBinContent = h1->GetBinContent(overflowBin);
      float overflowBinError   = h1->GetBinError  (overflowBin); 
      float lastBinContent     = h1->GetBinContent(lastBin);
      float lastBinError       = h1->GetBinError  (lastBin); 

      float totBinContent      = overflowBinContent + lastBinContent;
      float totBinError        = sqrt(overflowBinError*overflowBinError + lastBinError*lastBinError);
     
      h1->SetBinContent(lastBin, totBinContent);
      h1->SetBinError(lastBin, totBinError);
  }

  void moveUnderflowToLastBin(TH1* h1)
  {
      Int_t underflowBin = 0; // underflow bin
      Int_t firstBin     = 1; // first bin with low-edge xlow INCLUDED

      float underflowBinContent = h1->GetBinContent(underflowBin);
      float underflowBinError   = h1->GetBinError  (underflowBin); 
      float firstBinContent     = h1->GetBinContent(firstBin);
      float firstBinError       = h1->GetBinError  (firstBin); 

      float totBinContent      = underflowBinContent + firstBinContent;
      float totBinError        = sqrt(underflowBinError*underflowBinError + firstBinError*firstBinError);
     
      h1->SetBinContent(firstBin, totBinContent);
      h1->SetBinError(firstBin, totBinError);
  }

  string getHash()
  {
    string hash = "";
    for(size_t sampleIt = 0; sampleIt < this->size(); ++sampleIt)
    {
      hash += (*this)[sampleIt]->title_ + (*this)[sampleIt]->tcut_ + (*this)[sampleIt]->filename_;
    }
    std::hash<std::string> str_hash;
    hash = any2string(str_hash(hash));
    return hash;
  }
  SimpleStack *simple_stacks_[30];
  THStack *stacks_[30];
  TH1F *histos_[100];
  int scounter_ ;
  int hcounter_ ;
};
#endif
