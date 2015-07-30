#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include <string>
#include "../core/any2string.h"
#include "TSystem.h"

#include "TH1.h"

#ifndef SimpleSample_h
#define SimpleSample_h
// Class SimpleSample, with members a TFile, a TTree, a Title and a TCut(for the x-section)
class SimpleSample 
{
  public:
  SimpleSample(string filename, string title,  TCut tcut, bool fast = true, int fillColor = kWhite, int lineColor = kBlack,  int lineStyle = 1, int lineWidth = 2):
  title_(title),tcut_(tcut), fillColor_(fillColor), lineColor_(lineColor), lineStyle_(lineStyle), lineWidth_(lineWidth), filename_(filename)
  {
    fp_ = new TFile(filename.c_str(),"OPEN");
    events_ = (TTree*)fp_->Get("demo/events");
    maxN_ = events_->GetEntries();
    if(fast) maxN_ = ULong64_t(0.01*maxN_); // don't change this 0.01 unless you change also totEveW*100 when fast = true
    
    int filenameSize = filename.size();
    string friendTreePossiblePath = filename.substr(0,filenameSize-5) + ".aux.root";

    // "checking the existance of friend tree file in: friendTreePossiblePath 
    FileStat_t myfilestat;
    bool auxFileExists = !gSystem->GetPathInfo(friendTreePossiblePath.c_str(), myfilestat) ;
    if(auxFileExists){
        cout << "### a good friend was found in: " << friendTreePossiblePath.c_str() << endl;
	events_->AddFriend("events",friendTreePossiblePath.c_str());
    }

    if(tcut == TCut("0")) maxN_ = 0; // that's special for the null trivial sample
    //if(title != "Data" && !auxFileExists)tcut_ = tcut_*TCut( ("1./"+any2string(maxN_)).c_str()); // if aux file doesn't exists use 1/maxN_ for the MC normalization
    //if(title != "Data" && auxFileExists)tcut_ = tcut_*TCut("totEveW"); // if aux file exists use multiply with totEveW (takes into account 1/Ntot, PU etc)
    if(title != "Data")tcut_ = tcut_*TCut("totEveW"); // if aux file exists use multiply with totEveW (takes into account 1/Ntot, PU etc)

    // --- set aliases
    events_->SetAlias("jzb","t1vHT-l1l2Pt");
    events_->SetAlias("rawjzb","vHT-l1l2Pt");

    cout << "opening "<< filename << " ; title = "<< title_ << " ; TCut = " << tcut_.GetTitle() << " ; maxN = " << maxN_ << " ;  events_->GetEntries() = " << events_->GetEntries() << endl;
  } 

  void SetStyle(int fillColor, int lineStyle, int lineColor, int lineWidth)
  {
    fillColor_ = fillColor;
    lineStyle_ = lineStyle;
    lineColor_ = lineColor;
    lineWidth_ = lineWidth;
  }
 
  void SetHistoStyle(TH1*h1)
  {
      h1->SetFillColor(fillColor_);
      h1->SetLineStyle(lineStyle_);
      h1->SetLineColor(lineColor_);
      h1->SetLineWidth(lineWidth_);
      h1->GetXaxis()->SetNdivisions(509);
      h1->GetYaxis()->SetNdivisions(509);
  }
 
  TFile *fp_;
  TTree *events_;
  ULong64_t maxN_;
  TCut  tcut_;
  int  fillColor_;
  int  lineStyle_;
  int  lineColor_;
  int  lineWidth_;
  string title_;
  string filename_;
  TH1I *mcWeights_;
};
#endif
