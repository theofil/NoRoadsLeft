#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include <vector>
#include <string>
#include "TEfficiency.h"


#ifndef autoRebin_h
#define autoRebin_h
vector <float> getMyBins(TH1F *histo_in);                                // returns a vector with the low edges + the terminal high edge as last element
TH1F *reBinHisto(TH1F*histo_in, vector<float> myBins, bool divideByBinWidth = true);                   // returns a rebinned histogram for myBins
vector <float> findNiceBins(TH1F *histo_in, float thresh = 25);          // finds nice bins
void reBinMe(TH1F *&histo_in, vector<float> myBins, bool divideByBinWidth = true);                     // rebins a histo for a given array of bins
void autoRebin(TH1F *&histo_in, float thresh = 25, vector<float> myBins = vector<float>());
void autoRebin(TH1F *&h1,TH1F *&h2, float thresh = 25);
TH1F *doRatio(TH1F *h1,TH1F *h2, float thresh = 25, bool addOverflowBin = true);
TGraphAsymmErrors *doEffRatio(TH1F *h1,TH1F *h2, float thresh = 25);


TGraphAsymmErrors *doEffRatio(TH1F *h1,TH1F *h2, float thresh)
{
  TH1F *h1_clone = (TH1F*)h1->Clone("");
  TH1F *h2_clone = (TH1F*)h2->Clone("");
  vector<float> myBins = findNiceBins(h1_clone, thresh);
  reBinMe(h1_clone, myBins, false); // false is for not dividing by bin width in eff calculations
  reBinMe(h2_clone, myBins, false);
  TEfficiency *myeff =  new TEfficiency(*h1_clone, *h2_clone);
  TGraphAsymmErrors *mygraph = myeff->CreateGraph();
  mygraph->GetXaxis()->SetRangeUser(h1_clone->GetXaxis()->GetXmin(), h1_clone->GetXaxis()->GetXmax()); // TEfficiency x-axis is larger than the original histo
  return mygraph;
}

TH1F *doRatio(TH1F *h1,TH1F *h2, float thresh, bool addOverflowBin)
{
  TH1F *h1_clone = (TH1F*)h1->Clone("");
  TH1F *h2_clone = (TH1F*)h2->Clone("");
  autoRebin(h1_clone, thresh);
  reBinMe(h2_clone, getMyBins(h1_clone));
  h1_clone->Divide(h2_clone); 
  vector<float> myRatioBins = getMyBins(h1_clone);
  
  // if addOverflowBin = true  move overflow to the last bin of a new histogram 
  // --- first get overflow bin x-cordinates
  if(addOverflowBin) 
  { 
    Int_t nBins = h1->GetNbinsX();
    Int_t nOverflowBin = nBins + 1;
    float overflowBinWidth = h1->GetBinWidth(nOverflowBin);
    float overflowBinLowEdge = h2->GetBinLowEdge(nOverflowBin);
    float overflowBinHighEdge = overflowBinLowEdge + overflowBinWidth;
    myRatioBins.push_back(overflowBinHighEdge); // add overflow bin as last bin in the ratio
  }
  
  int myRatioBinsSize = myRatioBins.size();
  float ratioHistBins[myRatioBinsSize];
  for(int ii =0; ii<int(myRatioBins.size()) ; ii++) // copy myRatioBins inside the float array ratioHistBins
  {
    ratioHistBins[ii] = myRatioBins[ii];
   // cout << myRatioBins[ii] << endl;
  }

  string histName = string(h1->GetName()); 
  histName.pop_back(); histName.pop_back(); histName.pop_back();
  //histName += "_ratio";
  histName = "";

  TH1F *h_ratioHist = new TH1F(histName.c_str(),h1_clone->GetTitle(),myRatioBinsSize-1, ratioHistBins);
  h_ratioHist->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
  int h_ratioHistNBins = h_ratioHist->GetNbinsX();
  
  for(int ii =0; ii <= h_ratioHistNBins ; ii++) // start from overflow bin up to last bin of "h_ratioHist" which is overflow of h1_clone
  {
    float binContent = h1_clone->GetBinContent(ii);
    float binError   = h1_clone->GetBinError(ii);
    h_ratioHist->SetBinContent(ii, binContent);
    h_ratioHist->SetBinError(ii, binError);
  }
 

  delete h2_clone; 
  delete h1_clone;
  return h_ratioHist;
}

void autoRebin(TH1F *&h1,TH1F *&h2, float thresh)
{
  autoRebin(h1, thresh);
  reBinMe(h2, getMyBins(h1));
}

void autoRebin(TH1F *&histo_in, float thresh, vector<float> myBins)
{
  if (myBins == vector<float>()) // no bins are provided
  {
    myBins = findNiceBins(histo_in, thresh);
  }
  if (myBins != vector<float>()) // bins are provided
  {
    reBinMe(histo_in,myBins);
  }
}

vector <float> getMyBins(TH1F *histo_in)
{
  vector <float> myBins;

  int nBins = histo_in->GetNbinsX() ;
  float xMin = histo_in->GetXaxis()->GetXmin();
  float xMax = histo_in->GetXaxis()->GetXmax();
  
  for(int ii=1; ii<=nBins ; ii++) // get or lower edges
  {
    myBins.push_back(histo_in->GetBinLowEdge(ii));
  }
  myBins.push_back(histo_in->GetBinLowEdge(nBins) + histo_in->GetBinWidth(nBins)); // put the terminal "high" edge
  
  if(int(myBins.size())==0) cout << "myBins.size() = " << myBins.size() << " nBins = " << nBins << endl;
  if(myBins[nBins]!=xMax) cout << "segmentation problem myBins[nBins] = " << myBins[nBins] << " xMax = " << xMax << endl;
  if(myBins[0]!=xMin) cout << "segmentation problem myBins[0] = " << myBins[0] << " xMin = " << xMin << endl;
  
  return myBins;
}

void reBinMe(TH1F *&histo_in, vector<float> myBins, bool divideByBinWidth)
{
  TH1F *tmp2 = histo_in;
  TH1F *tmp = reBinHisto(histo_in, myBins, divideByBinWidth);
  histo_in = tmp;
  delete tmp2;
}

TH1F *reBinHisto(TH1F*histo_in, vector<float> myBins, bool divideByBinWidth)
{
   TH1F *histo_out;
   int myBinsArraySize = myBins.size();
   Double_t xbins[myBins.size()];
   for(int ii=0; ii<myBinsArraySize; ii++)xbins[ii] = myBins[ii];
   string histname = string(histo_in->GetName()) + "rbnd";
   histo_out = (TH1F*)histo_in->Rebin(myBinsArraySize-1,histname.c_str(),xbins);
   if(divideByBinWidth)histo_out->Scale(1,"width");
   return histo_out;
}

vector <float> findNiceBins(TH1F *histo_in, float thresh)
{
  bool debug = false;
  bool failed = false;
  vector <float> myNiceBins;
  vector <float> myTmpNiceBins;
  TH1F *hist_clone = (TH1F*)histo_in->Clone();

  int nBins = histo_in->GetNbinsX() ;
  float xMin = histo_in->GetXaxis()->GetXmin();
  float xMax = histo_in->GetXaxis()->GetXmax();
  if(debug) cout << "xMin = " << xMin << " xMax = " << xMax << " nBins = "<< nBins << endl;
 
  int firstNonZeroBin = hist_clone->FindFirstBinAbove(0.01);
  int lastNonZeroBin  = hist_clone->FindLastBinAbove(0.01);
  if(debug) cout << "firstNonZeroBin = " << firstNonZeroBin << " lastNonZeroBin = " << lastNonZeroBin << endl;

  float firstNonZeroBinX = hist_clone->GetBinLowEdge(firstNonZeroBin);
  float lastNonZeroBinX  = hist_clone->GetBinLowEdge(lastNonZeroBin) + hist_clone->GetBinWidth(lastNonZeroBin);
  if(debug) cout << "firstNonZeroBinX = " << firstNonZeroBinX << " lastNonZeroBin = " << lastNonZeroBinX << endl;
  
  // --- first check if the histo is already OK 
  string histname = string(histo_in->GetName());
  string forbidenKey = "rbnd";
  std::size_t found = histname.find(forbidenKey);
  if(found!=std::string::npos) 
  {
   cout << "forbidenKey found " << endl;
   return myNiceBins;
  }

  //--- start filling a tmp array of the non-empty part of the histo
  myTmpNiceBins.push_back(firstNonZeroBinX);
  float binSum=0;

  for(int ii=firstNonZeroBin; ii <= lastNonZeroBin ; ii++ )
  {
    float highEdge = histo_in->GetBinLowEdge(ii) + histo_in->GetBinWidth(ii);

    binSum += histo_in->GetBinContent(ii);

    if(binSum >= thresh) 
    {
      myTmpNiceBins.push_back(highEdge);     
      binSum =0;
    }

    if(ii == lastNonZeroBin && binSum<thresh) // if last bin reach and still below thresh merge last bin
    {
      myTmpNiceBins.pop_back();
      myTmpNiceBins.push_back(highEdge);     
    }        
  }
  if(debug) cout << "myTmpNiceBins.size() = " << myTmpNiceBins.size() << endl;
  if(debug) {for(unsigned int i=0; i<myTmpNiceBins.size() ; i++ )cout << myTmpNiceBins[i] << endl;}
 

  // instert as many bins are needed in "front" of myTmpNiceBins in order to start from the xMin
  float myCurrentXMin = myTmpNiceBins.front();
  int nTrials = 0;
  while( fabs(myCurrentXMin-xMin)>1.e-5)
  { 
    myCurrentXMin = myCurrentXMin - (float) hist_clone->GetBinWidth(firstNonZeroBin); 
    myTmpNiceBins.insert(myTmpNiceBins.begin(),myCurrentXMin);
    nTrials++;
    if(nTrials == 50) {cout << "WARNING: failed to find nice bins in autorebin, breaking here." << endl; failed = true; break;}
  }

  // --- now fill the final returned bins "myNiceBins"
  //int myTmpNiceBinsSize = myTmpNiceBins.size();

  for(int ii=0; ii<int(myTmpNiceBins.size()) ; ii++)
  {
    myNiceBins.push_back(myTmpNiceBins[ii]);
  }

  // add as many bins of original bin size are needed in order the end to be reached
  float originalBinWidth = histo_in->GetBinWidth(nBins);
  float currentEndEdge = myNiceBins.back(); 
  nTrials = 0;
  while(fabs(xMax- currentEndEdge)>1.e-5)
  {
    myNiceBins.push_back(currentEndEdge + originalBinWidth);
    currentEndEdge = myNiceBins.back();
    nTrials++;
    if(nTrials == 50) {cout << "WARNING: failed to find nice bins in autorebin, breaking here." << endl;failed = true; break;}
  }

  if(debug) cout << "printing out final myNiceBins" << endl;
  if(debug) for(unsigned int i=0; i<myNiceBins.size() ; i++ ) cout << myNiceBins[i] << endl;
  if(failed) {cout << "failed: returning original bins" << endl;return getMyBins(histo_in);}
  return myNiceBins;
}
#endif
