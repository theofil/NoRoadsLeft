#ifndef utilities_h
#define utilities_h

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <map>
#include <set>

#include "TFile.h"
#include "TLine.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TF1.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TBox.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TPolyLine.h"
#include "setTDRStyle.C" // README we need this in the same dir  "./" with the root macro or make c++ to know where to find it
#include "pnumber.h"
#include "TRandom3.h"
#include "TChain.h"
#include "THStack.h"
#include "TEntryList.h"

#include "any2string.h"

using namespace std;

float ratioError(float a,float b);
float ratioError(float a,float b, float da , float db);
float productError(float a,float b);
float productError(float a,float b, float da , float db);
float histIntegral(TH1F *hist, float minX){return hist->Integral(hist->FindBin(double(minX)),hist->GetNbinsX()+1);}
float histIntegral(TH1F *hist, float minX, float maxX){return hist->Integral(hist->FindBin(double(minX)),hist->FindBin(double(maxX)));}
float histIntegralAndError(TH1F *hist, float minX, float & error);
float histIntegralAndError(TH1F *hist, float minX, float maxX, float & error); 
pnumber histIntegralPN(TH1F*hist, float minX, float maxX);
void  errorBand(TH1F* hist, float sysUncert = 1);
TH1F *squareRootHist(TH1F*);
float significance(float,float);
float effRatioError(float A, float B);
float statErrorN(float x){return x - 0.5*TMath::ChisquareQuantile(0.3173/2,2*x);}
float statErrorP(float x){return 0.5*TMath::ChisquareQuantile(1-0.3173/2,2*(x+1))-x;}
pnumber nEveW(TH1F* h1);
pnumber nEve(TH1F* h1);
double getNSigma(double prob);
TLegend *simpleLegend(float x1ndc = 0.53, float x2ndc = 0.65,  float y1ndc = 0.89, float y2ndc =0.89);
TF1 * chi2Dist; // chi2Dist = new TF1("fgamma", "TMath::GammaDist(x, [0], [1], [2])", 0, 10); chi2Dist->SetParameters(ndf/2, 0 , 2); // execute this in main scope


float significance(float n, float b)
{
   return (n-b)/sqrt(n + b);
}

void errorBand(TH1F*hist, float sysUncert)
{
  for(int i =0 ; i <= hist->GetNbinsX()+1;i++)
  {
    float sysError  =  hist->GetBinContent(i)*sysUncert; 
    float statError =  sqrt(hist->GetBinContent(i));
    float totError  =  sqrt(sysError*sysError + statError*statError); 
    hist->SetBinError(i,totError);
  } 
}

TH1F  *squareRootHist(TH1F *hist)
{
  TH1F *newHist = (TH1F*)hist->Clone();
  for(int i =0 ; i <= newHist->GetNbinsX()+1;i++)
  {
    float binContent = newHist->GetBinContent(i);
    float binError  = newHist->GetBinError(i);
   
    float myBinContent=0; 
    float myBinError=0;
    if(binContent>0) 
    {
      myBinContent = sqrt(binContent);
      myBinError = 0.5*binError*(1/myBinContent);
    }
    newHist->SetBinError(i,myBinError);
    newHist->SetBinContent(i,myBinContent);
  } 
  return newHist;
}

void normHist(TH1F *hist)
{
 hist->Scale(1.0/hist->Integral());
 hist->GetYaxis()->SetTitle("PDF");
}

TH1F  *normHistPointer(TH1F *hist)
{
 TH1F *newHist = (TH1F*)hist->Clone();
 newHist->Scale(1.0/hist->Integral());
 newHist->GetYaxis()->SetTitle("PDF");
 return newHist;
}


float ratioError(float a,float b)
{
  return sqrt(a/(a*a) + b/(b*b))*(a/b);
}

float ratioError(float a,float b, float da , float db)
{
  return sqrt(pow(da/a,2) + pow(db/b,2))*(a/b);
}

float productError(float a,float b)
{
  return sqrt(a/(a*a) + b/(b*b))*(a*b);
}

float productError(float a,float b, float da , float db)
{
  return sqrt(pow(da/a,2) + pow(db/b,2))*(a*b);
}

float effRatioError(float A, float B)
{
   float eff = A/B;
   float tmpError2 = eff*(1-eff)/B;
   return sqrt(tmpError2);
} 

float histIntegralAndError(TH1F *hist, float minX, float &error)
{
  //      bin = 0;       underflow bin
  //      bin = 1;       first bin with low-edge xlow INCLUDED
  //      bin = nbins;   last bin with upper-edge xup EXCLUDED
  //      bin = nbins+1; overflow bin
  float val = hist->Integral(hist->FindBin(double(minX)),hist->GetNbinsX()+1);
  double valError=0.;
  hist->IntegralAndError(hist->FindBin(double(minX)),hist->GetNbinsX()+1,valError);
  error = float(valError);
  return val;
}


float histIntegralAndError(TH1F *hist, float minX, float maxX, float &error)
{
  float val = hist->Integral  (hist->FindBin(double(minX)),hist->FindBin(double(maxX)));
  double valError=0.;
  hist->IntegralAndError(hist->FindBin(double(minX)),hist->FindBin(double(maxX)),valError);
  error = float(valError);
  return val;
}

pnumber histIntegralPN(TH1F*hist, float minX, float maxX)
{
  float integral      = 0.;
  float integralError = 0.; 
  integral  =  histIntegralAndError(hist, minX, maxX, integralError);
  return pnumber(integral, integralError) ;
}


pnumber nEveW(TH1F* h1)
{ 
  pnumber val;
  int nBins = h1->GetNbinsX();
  val.x = h1->Integral(0, nBins+1);
 
  Double_t error =0 ; h1->IntegralAndError(0,nBins+1,error);
  val.xE = float(error);

  return val;
}

pnumber nEve(TH1F* h1)
{ 
  pnumber val(h1->GetEntries());
  return val;
}


double getNSigma(double prob) // with prob to be the one-sided Gaussian probability: p = 1-CDF
{
  return sqrt(2)*TMath::ErfInverse(1-2*prob);
}


TH1F  *fractionalUncertTH1F(TH1F *hist) // changes a histogram to a fraction uncertainty band [error/binContent]
{
  TH1F *newHist = (TH1F*)hist->Clone();

  for(int i =0 ; i <= newHist->GetNbinsX()+1;i++)
  {
    float binContent = newHist->GetBinContent(i);
    float binError  = newHist->GetBinError(i);
   
    float myBinContent=0; 
    float myBinError=0;
    myBinContent = 1.0;

    if(binContent != 0) myBinError = sqrt(binError)/binContent;
    if(binContent == 0) myBinError = 0;
  
    if (binContent != 0) 
    { 
      newHist->SetBinError(i,myBinError);
      newHist->SetBinContent(i,myBinContent);
    }
  } 
  return newHist;
}

TH1F  *ScaleTH1F(TH1F *hist, pnumber scale) 
{
  TH1F *newHist = (TH1F*)hist->Clone();

  for(int i =0 ; i <= newHist->GetNbinsX()+1;i++)
  {
    float binContent = newHist->GetBinContent(i);
    float binError  = newHist->GetBinError(i);
   
    pnumber newBinContent(binContent,binError);
    newBinContent = newBinContent*scale; 
 
    newHist->SetBinContent(i, newBinContent.x);
    newHist->SetBinError  (i, newBinContent.xE);
  } 
  return newHist;
}


float coolme(float x, float b){return round(x/b)*b;}

#endif

