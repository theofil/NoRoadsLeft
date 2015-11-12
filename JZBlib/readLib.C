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
#include "../friends/readNTPL.h"

#include "TVector2.h"

using namespace std;
bool goFast(false);

TH1F *h_jzb_true_ej1;    
TH1F *h_jzb_smeared_ej1;

TH1F *h_jzb_true_ij2;    
TH1F *h_jzb_smeared_ij2;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
 
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
    bool   isTheBin(float pt, float eta){bool res = false; if(pt >= PtMin_ && pt <= PtMax_ && fabs(eta) >= AbsEtaMin_ && fabs(eta) <= AbsEtaMax_) res = true; return res;}

    private:
    float PtMin_;
    float PtMax_;
    float AbsEtaMin_;
    float AbsEtaMax_;
    string tcut_string_;
    TCut  tcut_; 
};

int findPtEtaBin(vector<PtAbsEtaBin> myPtEtaBins, float pt, float eta)
{
    int res=-1;

    for(int bin = 0 ; bin < int(myPtEtaBins.size()) ; bin++)
    {
	if(myPtEtaBins[bin].isTheBin(pt, eta)) res = bin;	
    }
 
    if(res == -1) cout << "(" << pt << "," << eta << ") bin was not found " << endl;
    return res;	    
}

TFile *fp_in_local;
TH1F *histo[30];
TF1  *fgauss[30];

double fnc_dscb(double*xx,double*pp);

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


  cout << "..:: makeLib log ::..." << endl;
  cout << endl;

  fp_in_local = new TFile("../data/JZBlib_mc.root", "READ");

  TF1* fdscb[myPtEtaBins.size()];
  TF1* fgaus[myPtEtaBins.size()];

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
  
    histo[bin] = (TH1F*)fp_in_local->Get((plotTitle+"_TH1F").c_str());

    fdscb[bin] = new TF1(("fdscb_b"+any2string(bin)).c_str(),fnc_dscb, -200,200,7);
    fgaus[bin] = new TF1(("fgaus_b"+any2string(bin)).c_str(),"gaus",-10,10);
    fgaus[bin]->SetLineColor(kBlue);

//    TCanvas *cantmp = new TCanvas();
//    cantmp->SetLogy();
    SimpleCanvas *simpleCan = new SimpleCanvas(plotTitle, 1);
    simpleCan->SetLogy();
    histo[bin]->Fit(fgaus[bin]);
    fdscb[bin]->FixParameter(0, fgaus[bin]->GetParameter(0));
    fdscb[bin]->FixParameter(1, fgaus[bin]->GetParameter(1));
    fdscb[bin]->FixParameter(2, fgaus[bin]->GetParameter(2));
    fdscb[bin]->SetParameter(3, 1);
    fdscb[bin]->SetParameter(4, 0.5);
    fdscb[bin]->SetParameter(5, 1);
    fdscb[bin]->SetParameter(6, 0.5);
    histo[bin]->Fit(fdscb[bin]);
    histo[bin]->Draw();
   // fgaus[bin]->Draw("same");
    fdscb[bin]->Draw("same");
    simpleCan->Save(outputDir);
  }

  
  // try first to make sense in the first jet bin
  TCut sel_magic       = TCut("((lepID[0]*lepID[1] == -11*13) ? -1 : 0) + ((lepID[0]*lepID[1] == -11*11) ? 1 : 0) + ((lepID[0]*lepID[1] == -13*13) ? 1 : 0)");
  TCut sel_cut         = (sel_basic && sel_M81101 && sel_ej1)*sel_magic;

  string v_met = "met";
  string v_hr  = "vHT";
  string jzbvar        = (v_hr+"-l1l2Pt").c_str();

 
  v1_drivers(goFast);

  SimpleDriver myDriver = mcDriver;


  TChain *ch = myDriver.getTChain();
  init(ch);
  float totEveW_;
  float xsection_;
  ch->SetBranchAddress("totEveW", &totEveW_);   
  ch->SetBranchAddress("xsection", &xsection_);

  ULong64_t chainEntries = ch->GetEntries();
  if(goFast) chainEntries = 0.01*chainEntries;

//    sel_basic && sel_M81101 && sel_ej1; 

  h_jzb_true_ej1    = new TH1F("h_jzb_true_ej1", ";JZB [GeV]; events / 5 GeV", 80, -200, 200);
  h_jzb_smeared_ej1 = new TH1F("h_jzb_smeared_ej1", ";JZB [GeV]; events / 5 GeV", 80, -200, 200);
  h_jzb_true_ij2    = new TH1F("h_jzb_true_ij2", ";JZB [GeV]; events / 5 GeV", 80, -200, 200);
  h_jzb_smeared_ij2 = new TH1F("h_jzb_smeared_ij2", ";JZB [GeV]; events / 5 GeV", 80, -200, 200);

  for(ULong64_t ii=0; ii < chainEntries;ii++)
  {
    ch->GetEntry(ii);
    bool bool_sel_M81101 = l1l2M_>81 && l1l2M_<101; 
    float float_sel_magic = ((lepID_[0]*lepID_[1] == -11*13) ? -1 : 0) + ((lepID_[0]*lepID_[1] == -11*11) ? 1 : 0) + ((lepID_[0]*lepID_[1] == -13*13) ? 1 : 0);
    bool bool_basic = goodVtx_ && lepPt_[0]>25 && lepPt_[1]>20 && l1l2DR_>0.3 && !(fabs(lepEta_[0]) > 1.4 && fabs(lepEta_[0]) < 1.6) && !(fabs(lepEta_[1]) > 1.4 && fabs(lepEta_[1]) < 1.6);
    float eventW = totEveW_*xsection_*42;
    bool bool_sel_ej1 = njets_==1;
    bool bool_sel_ij2 = njets_>=2;
   
    if(bool_sel_M81101 && bool_basic && bool_sel_ej1)
    {
	float jzbvar_ = vHT_ - l1l2Pt_;
        h_jzb_true_ej1->Fill(jzbvar_, float_sel_magic*eventW);
    }

    if(bool_sel_M81101 && bool_basic && bool_sel_ij2)
    {
	float jzbvar_ = vHT_ - l1l2Pt_;
        h_jzb_true_ij2->Fill(jzbvar_, float_sel_magic*eventW);
    }

    if(bool_sel_M81101 && bool_basic && bool_sel_ej1 && met_ < 20)
    {
/*
	float jzbvar_ = vHT_ - l1l2Pt_;
        
        for(int jj =0 ; jj < 1; jj++){
	    float tmp_jetPt = jetPt_[0];    
	    float tmp_jetEta = jetEta_[0];    
	    jzbvar_ =   fdscb[findPtEtaBin(myPtEtaBins, tmp_jetPt, tmp_jetEta)]->GetRandom();
	    h_jzb_smeared_ej1->Fill(jzbvar_, float_sel_magic*eventW);
        }
*/  

    	TVector2 jzbVec(0,0);
      	int nGen = 2;
	for(int jj =0 ; jj < nGen; jj++)
	{
	    float jzbvar_ = vHT_ - l1l2Pt_;
	    for(int nj =0 ; nj < njets_; nj++)
	    {
		float tmp_jetPt  = jetPt_[nj];    
	    	float tmp_jetEta = jetEta_[nj];
	    	float tmp_jetPhi = jetPhi_[nj];
	      	float jzb_mag = fdscb[findPtEtaBin(myPtEtaBins, tmp_jetPt, tmp_jetEta)]->GetRandom();
		if(jzb_mag < 0) {jzb_mag = -jzb_mag; tmp_jetPhi = -tmp_jetPhi;}
                
		jzbVec += TVector2(fdscb[findPtEtaBin(myPtEtaBins, tmp_jetPt, tmp_jetEta)]->GetRandom(), tmp_jetPhi);
	    }
            
	    jzbvar_ = jzbVec.Mod();
            
	    h_jzb_smeared_ej1->Fill(jzbvar_, float_sel_magic*eventW);
	}
   }

/*
    if(bool_sel_M81101 && bool_basic && bool_sel_ij2 && met_ < 20)
    {
	float jzbvar_ = vHT_ - l1l2Pt_;
        TVector2 jzbVec(0,0);
        int nGen = 2;
        for(int jj =0 ; jj < nGen; jj++)
	{
            for(int nj =0 ; nj < njets_; nj++)
	    {
                float tmp_jetPhi = jetPhi_[nj];
                float jzb_mag = fdscb[findPtEtaBin(myPtEtaBins, tmp_jetPt, tmp_jetEta)]->GetRandom();
                if(jzb_mag < 0) {jzb_mag = -jzb_mag; tmp_jetPhi = tmp_jetPhi + 3.14159;}

                jzbVec += TVector2(fdscb[findPtEtaBin(myPtEtaBins, tmp_jetPt, tmp_jetEta)]->GetRandom(), tmp_jetPhi);
	    }
	    jzbvar_ = jzbVec.Mod();
	    h_jzb_smeared_ij2->Fill(jzbvar_, float_sel_magic*eventW);
        }
    }
*/
    
  } 


  {
    string plotTitle     = "ej1_closure";
    SimpleStack * ss_tmp    = myDriver.getSimpleStackTH1F(jzbvar,";JZB [GeV]; events / 5 GeV;", 80, -200, 200, sel_cut);
    TH1F *h1_tmp    = myDriver.getHistoTH1F(ss_tmp);	
    SimpleLegend *sleg = new SimpleLegend("TLSF");
    SimpleCanvas *simpleCan = new SimpleCanvas(plotTitle, 1);
    ss_tmp->Draw("hist"); 
    ss_tmp->GetXaxis()->SetTitle(h1_tmp->GetXaxis()->GetTitle());
    h1_tmp->Draw("hist same");
    h1_tmp->Draw("axis same");
    h_jzb_true_ej1->Draw("e1 same");
    sleg->FillLegend(ss_tmp);
    sleg->Draw("same");
    simpleCan->Save(outputDir);
  }
  
  TH1F *h_jzb_smeared_ej1_clone = (TH1F*)h_jzb_smeared_ej1->Clone();
  TH1F *h_jzb_true_ej1_clone = (TH1F*)h_jzb_true_ej1->Clone();
  normHist(h_jzb_smeared_ej1_clone);
  normHist(h_jzb_true_ej1_clone);

  h_jzb_smeared_ej1_clone ->SetLineColor(kRed);

  {
    string plotTitle     = "ej1_smeared_closure";
    SimpleCanvas *simpleCan = new SimpleCanvas(plotTitle, 2);
    simpleCan->Up();
    simpleCan->ShapeMeUp(h_jzb_true_ej1_clone);
    simpleCan->ShapeMeUp(h_jzb_smeared_ej1_clone);
    h_jzb_true_ej1_clone->DrawNormalized("e1");
    h_jzb_smeared_ej1_clone->DrawNormalized("hist same");
    simpleCan->SetLogy();
    simpleCan->Dw();
    TH1F *hratio = doRatio(h_jzb_smeared_ej1_clone, h_jzb_true_ej1_clone);
    simpleCan->ShapeMeDw(hratio);
    hratio->GetYaxis()->SetTitle("pred. / true");
    hratio->Draw("e1");
    hratio->GetYaxis()->SetRangeUser(0.0, 2.0);
    hratio->GetYaxis()->SetNdivisions(507);
    simpleCan->Save(outputDir);
  }

  

  


}

// from https://github.com/cms-analysis/JetMETAnalysis-JetAnalyzers/blob/master/bin/jet_response_fitter_x.cc

double fnc_dscb(double*xx,double*pp)
{
    double x   = xx[0];
    // gaussian core
    double N   = pp[0];//norm
    double mu  = pp[1];//mean
    double sig = pp[2];//variance

    // transition parameters
    double a1  = pp[3];
    double p1  = pp[4];
    double a2  = pp[5];
    double p2  = pp[6];
   
    double u   = (x-mu)/sig;
    double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
    double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
    double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
    double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);
    
    double result(N);
    if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
    else if (u<a2)  result *= TMath::Exp(-u*u/2);
    else            result *= A2*TMath::Power(B2+u,-p2);
   
    return result;
}

