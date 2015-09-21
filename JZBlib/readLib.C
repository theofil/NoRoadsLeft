#include "../core/utilities.h"
#include "../core/simpleROOT_cuts.h"
#include "../core/autoRebin.h"
#include "../core/SimpleCanvas.h"
#include "../core/SimpleSample.h"
#include "../core/SimpleDriver.h"
#include "../core/SimpleLegend.h"
#include "../core/SimpleStack.h"
#include "../core/CrossSections13TeV.h"
#include "../core/LocalFilePath.h"
#include "../core/SimplePaveText.h"
#include "../example/v1_drivers.C"

using namespace std;
bool goFast(false);

class PtAbsEtaBin
{
    public:
    PtAbsEtaBin(float PtMin, float PtMax, float AbsEtaMin, float AbsEtaMax):PtMin_(PtMin), PtMax_(PtMax), AbsEtaMin_(AbsEtaMin), AbsEtaMax_(AbsEtaMax)
    {
      string tcut_string_ = "jetPt[0] > " + any2string(PtMin_) + " && " + "jetPt[0] < " + any2string(PtMax_) + " && " + " abs(jetEta[0]) > " + any2string(AbsEtaMin_);
      tcut_string_ +=  " && abs(jetEta[0]) < " + any2string(AbsEtaMax_);
      tcut_ = TCut(tcut_string_.c_str()); 
    }
    string getTitle(){return any2string(PtMin_) + " < P_{T} (jet) < " + any2string(PtMax_) + " , " + any2string(AbsEtaMin_) + " < |#eta (jet)| < " + any2string(AbsEtaMax_); }
    TCut   getTCut() {return tcut_;}

    private:
    float PtMin_;
    float PtMax_;
    float AbsEtaMin_;
    float AbsEtaMax_;
    string tcut_string_;
    TCut  tcut_; 
};


TFile *fp_in;
TH1F *histo[30];
TF1  *fgauss[30];


void readLib()
{ 
  setTDRStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.5);

  vector<PtAbsEtaBin> myPtEtaBins;
  myPtEtaBins.push_back(PtAbsEtaBin(100, 6500, 1.4, 2.4));
  myPtEtaBins.push_back(PtAbsEtaBin(75, 100, 1.4, 2.4));
  myPtEtaBins.push_back(PtAbsEtaBin(50, 75, 1.4, 2.4));
  myPtEtaBins.push_back(PtAbsEtaBin(40, 50, 1.4, 2.4));
  myPtEtaBins.push_back(PtAbsEtaBin(30, 40, 1.4, 2.4));

  myPtEtaBins.push_back(PtAbsEtaBin(100, 6500, 0, 1.4));
  myPtEtaBins.push_back(PtAbsEtaBin(75, 100, 0, 1.4));
  myPtEtaBins.push_back(PtAbsEtaBin(50, 75, 0, 1.4));
  myPtEtaBins.push_back(PtAbsEtaBin(40, 50, 0, 1.4));
  myPtEtaBins.push_back(PtAbsEtaBin(30, 40, 0, 1.4));

  v1_drivers(goFast);

  cout << "..:: makeLib log ::..." << endl;
  cout << endl;

  fp_in = new TFile("../data/JZBlib_mc.root", "READ");

  for(size_t bin = 0; bin < myPtEtaBins.size(); ++bin )
  {
    PtAbsEtaBin myPtEtaBin = myPtEtaBins[bin];
    cout << myPtEtaBin.getTitle() << endl;
    cout << myPtEtaBin.getTCut().GetTitle() << endl;
  
    string v_met = "met";
    string v_hr  = "vHT";
    string jzbvar        = (v_hr+"-l1l2Pt").c_str();

    string plotTitle     = "bin_"+any2string(bin);
    string ds            = " (mc)";
    SimpleDriver myDriver = mcDriver;
  
    histo[bin] = (TH1F*)fp_in->Get((plotTitle+"_TH1F").c_str());
//    histo[bin]->Sumw2(false);
//    histo[bin]->SetBinErrorOption(TH1::kPoisson);

/*
    fgauss[bin] = new TF1((plotTitle+"_TH1F"+"_TF1").c_str(), "gaus", -200, 200);
    fgauss[bin]->SetLineColor(kBlue);
    fgauss[bin]->SetLineWidth(2);
    histo[bin]->Fit(fgauss[bin],"R");
*/
    
//    h1_tmp->SetTitle((myPtEtaBin.getTitle() + ds).c_str());
//    h1_tmp->SetName((plotTitle+"_TH1F").c_str());
  }
}



