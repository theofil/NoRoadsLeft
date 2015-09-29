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


TFile *fp_out;


void makeLib()
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

  fp_out = new TFile("/tmp/theofil/JZBlib_mc.root", "RECREATE");
  fp_out->cd();

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
  
    TCut sel_magic       = TCut("((lepID[0]*lepID[1] == -11*13) ? -1 : 0) + ((lepID[0]*lepID[1] == -11*11) ? 1 : 0) + ((lepID[0]*lepID[1] == -13*13) ? 1 : 0)");
    TCut sel_cut         = (sel_basic && sel_M81101 && sel_ej1 && myPtEtaBin.getTCut())*sel_magic;

    SimpleStack * ss_tmp    = myDriver.getSimpleStackTH1F(jzbvar,";JZB [GeV]; events / 5 GeV;", 80, -200, 200, sel_cut);
    TH1F *h1_tmp    = myDriver.getHistoTH1F(ss_tmp);	
    SimpleLegend *sleg = new SimpleLegend("TLSF");
    SimpleCanvas *simpleCan = new SimpleCanvas(plotTitle, 1);
    ss_tmp->Draw("hist"); 
    ss_tmp->GetXaxis()->SetTitle(h1_tmp->GetXaxis()->GetTitle());
    h1_tmp->Draw("hist same");
    h1_tmp->Draw("axis same");
    sleg->FillLegend(ss_tmp);
    sleg->Draw("same");
    TPaveText *pave = new TPaveText(0.14, 0.94, 0.94, 0.98,"blNDC");
    pave->SetBorderSize(0);
    pave->SetFillColor(0);
    pave->SetFillStyle(0);
    pave->SetTextFont(43);
    pave->SetTextColor(kBlue);
    pave->SetTextSize(18);
    pave->AddText((myPtEtaBin.getTitle() + ds).c_str());
    pave->Draw("same");
    simpleCan->Save(outputDir);

    h1_tmp->SetTitle((myPtEtaBin.getTitle() + ds).c_str());
    h1_tmp->SetName((plotTitle+"_TH1F").c_str());
    simpleCan->can_->Write();
    h1_tmp->Write();
  }

  // here start the event loop 

}



