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
    string getTitle(){return any2string(PtMin_) + " < P_T (jet) < " + any2string(PtMax_) + " and " + any2string(AbsEtaMin_) + " < |eta(jet)| < " + any2string(AbsEtaMax_); }
    TCut   getTCut() {return tcut_;}

    private:
    float PtMin_;
    float PtMax_;
    float AbsEtaMin_;
    float AbsEtaMax_;
    string tcut_string_;
    TCut  tcut_; 
};

void makeLib()
{ 


  setTDRStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.5);

  vector<PtAbsEtaBin> myPtEtaBins;
  myPtEtaBins.push_back(PtAbsEtaBin(100, 6500, 1.4, 2.4));
  
  v1_drivers(goFast);

  cout << "..:: makeLib log ::..." << endl;
  cout << endl;

  for(size_t bin = 0; bin < myPtEtaBins.size(); ++bin )
  {
    PtAbsEtaBin myPtEtaBin = myPtEtaBins[bin];
    cout << myPtEtaBin.getTitle() << endl;
    cout << myPtEtaBin.getTCut().GetTitle() << endl;
  
    string v_met = "met";
    string v_hr  = "vHT";
    string jzbvar        = (v_hr+"-l1l2Pt").c_str();

    string plotTitle     = "bin_"+any2string(bin);
    SimpleDriver myDriver = dataDriver;
  
    TCut sel_magic       = TCut("((lepID[0]*lepID[1] == -11*13) ? -1 : 0) + ((lepID[0]*lepID[1] == -11*11) ? 1 : 0) + ((lepID[0]*lepID[1] == -13*13) ? 1 : 0)");
    TCut sel_cut         = (sel_basic && sel_M81101 && sel_ej1 && myPtEtaBin.getTCut())*sel_magic;
  
    // --- plot JZB
    {
        string title = "jzb"+plotTitle;
        SimpleStack * sstack    = myDriver.getSimpleStackTH1F(jzbvar,";JZB [GeV]; events / 5 GeV;", 80, -200, 200, sel_cut);
        TH1F *sstack_hall       = myDriver.getHistoTH1F(sstack);	
        SimpleLegend *sleg = new SimpleLegend("TLSF");
        SimpleCanvas *simpleCan = new SimpleCanvas(title.c_str(), 1);
        sstack->Draw("hist");
        sstack_hall->Draw("hist same");
        simpleCan->CMSPre();
        simpleCan->SetLogy();
        sleg->FillLegend(sstack);
        sleg->Draw("same");
        simpleCan->Save(outputDir);
    }
  }

}



