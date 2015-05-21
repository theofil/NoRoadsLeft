#include <functional>
#include <string>
#include "any2string.h"
#include "SimpleSample.h"
#include "SimpleStack.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "THStack.h"
#include "TCut.h"

#ifndef SimpleDriver_h
#define SimpleDriver_h


class SimpleDriver: public std::vector<SimpleSample*>
{
  public: 
  SimpleDriver(){h1counter_ = 0; scounter_ = 0; h2counter_ = 0; pXcounter_ = 0;}

  SimpleStack * getSimpleStackTH1F(string var, string histTitle, int nBinsX, float minX, float maxX, TCut myCut)
  {
    string stackName = var + "_s" + any2string(scounter_);  
    simple_stacks_[scounter_] = new SimpleStack(stackName.c_str(), stackName.c_str());
    scounter_++;

    for(size_t sampleIt = 0; sampleIt < this->size(); ++sampleIt)
    {
     string histName = var + "_h" + any2string(h1counter_) + "_" + getHash();
     h1_[h1counter_] = new TH1F(histName.c_str(), histTitle.c_str(), nBinsX, minX, maxX);
     h1_[h1counter_]->Sumw2();
     (*this)[sampleIt]->SetHistoStyle(h1_[h1counter_]); 
     h1counter_++;
     string drawCommand = var + ">>" + histName;
     (*this)[sampleIt]-> events_->Draw(drawCommand.c_str(),  myCut*(*this)[sampleIt]->tcut_, "goff", (*this)[sampleIt]->maxN_);

     moveOverflowToLastBin ( h1_[h1counter_-1] );
     moveUnderflowToLastBin( h1_[h1counter_-1] );

     simple_stacks_[scounter_-1]->Add(h1_[h1counter_-1]);
     simple_stacks_[scounter_-1]->sampleTitles_.push_back( (*this)[sampleIt]->title_);
     simple_stacks_[scounter_-1]->histoPointers_.push_back( h1_[h1counter_-1]);
    }
    return  simple_stacks_[scounter_-1];
  }

  THStack * getStackTH1F(string var, string histTitle, int nBinsX, float minX, float maxX, TCut myCut)
  {
    string stackName = var + "_s" + any2string(scounter_);  
    stacks_[scounter_] = new THStack(stackName.c_str(), stackName.c_str());
    scounter_++;

    for(size_t sampleIt = 0; sampleIt < this->size(); ++sampleIt)
    {
     string histName = var + "_h_" + any2string(h1counter_) + "_" + getHash();
     h1_[h1counter_] = new TH1F(histName.c_str(), histTitle.c_str(), nBinsX, minX, maxX);
     h1_[h1counter_]->Sumw2();
     (*this)[sampleIt]->SetHistoStyle(h1_[h1counter_]); 
     h1counter_++;
     string drawCommand = var + ">>" + histName;
     (*this)[sampleIt]-> events_->Draw(drawCommand.c_str(),  myCut*(*this)[sampleIt]->tcut_, "goff", (*this)[sampleIt]->maxN_);

     moveOverflowToLastBin ( h1_[h1counter_-1] );
     moveUnderflowToLastBin( h1_[h1counter_-1] );

     stacks_[scounter_-1]->Add(h1_[h1counter_-1]);
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

  TH2F * getHistoTH2F(string var, string histTitle, int nBinsX, float minX, float maxX, int nBinsY, float minY, float maxY, TCut myCut)
  {
    string histName = "h2_" + any2string(h2counter_) + "_" + getHash(var+any2string(nBinsX*(minX+maxX)*nBinsY*(minY+maxY))); 
    h2_[h2counter_] = new TH2F(histName.c_str(), histTitle.c_str(), nBinsX, minX, maxX, nBinsY, minY, maxY); 
    h2_[h2counter_]->Sumw2();
    h2counter_++;
    for(size_t sampleIt = 0; sampleIt < this->size(); ++sampleIt)
    {
     string drawCommand = var + ">>+" + histName;
     (*this)[sampleIt]-> events_->Draw(drawCommand.c_str(),  myCut*(*this)[sampleIt]->tcut_, "goff", (*this)[sampleIt]->maxN_);
    }
    return  h2_[h2counter_-1];
  }

  TH1F *getHistoTH1F(string var, string histTitle, int nBinsX, float minX, float maxX, TCut myCut) // uses getHistoTH1F(THStack *hs) but first creates the stack
  {
    THStack *hs_temp = getStackTH1F(var, histTitle, nBinsX, minX, maxX, myCut);
    return getHistoTH1F(hs_temp);
  }
  //	TProfile(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup, Option_t* option = "")

  TProfile *getProfileX(string var, string histTitle, std::vector<float> xBins, TCut myCut)
  {
        Int_t nBins = xBins.size();
        Float_t localXBins[nBins];

        for(int i =0 ; i<nBins; i++) localXBins[i] = xBins[i];
        string histName = "pX_" + any2string(pXcounter_) + "_" + getHash(var+any2string(nBins)); 
	pX_[pXcounter_] = new TProfile(histName.c_str(), histTitle.c_str(), nBins-1, localXBins);

        pXcounter_++;
  	for(size_t sampleIt = 0; sampleIt < this->size(); ++sampleIt)
      	{
	    string drawCommand = var + ">>+" + histName;
	    (*this)[sampleIt]-> events_->Draw(drawCommand.c_str(),  myCut*(*this)[sampleIt]->tcut_, "goff", (*this)[sampleIt]->maxN_);
	}
	return  pX_[pXcounter_-1];
  }



  void moveOverflowToLastBin(TH1* h1)
  {
      Int_t nBinsX = h1->GetNbinsX();
      Int_t overflowBin = nBinsX + 1;
      Int_t lastBin = nBinsX;


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

  string getHash(string input)
  {
    string hash = "";
    for(size_t sampleIt = 0; sampleIt < this->size(); ++sampleIt)
    {
      hash += (*this)[sampleIt]->title_ + (*this)[sampleIt]->tcut_ + (*this)[sampleIt]->filename_ + input;
    }
    std::hash<std::string> str_hash;
    hash = any2string(str_hash(hash));
    return hash;
  }

  SimpleStack *simple_stacks_[30];
  THStack *stacks_[30];
  TH1F *h1_[50];
  TH2F *h2_[50];
  TProfile *pX_[50];
  int scounter_ ;
  int h1counter_ ;
  int h2counter_ ;
  int pXcounter_ ;
};
#endif
