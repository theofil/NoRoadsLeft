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

bool goFast(false);

SimpleCanvas *simpleCan;

SimpleStack *hs_mc, *hs_data;
TH2F        *h2_mc, *h2_data;

void basicDataMCPlot(string title, SimpleStack *hs_data, SimpleStack *hs_mc, bool setLog = true , string legPos = "TR");
void basicDataMCPlotFSSub(string title, SimpleStack *hs_data, SimpleStack *hs_mc, bool setLog = true , string legPos = "TR");


void exZZ()
{
  setTDRStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.5);

  v1_drivers(goFast);
  
  TCut mySelection_SF = sel_basic && sel_SF && sel_SF_trig && sel_ej0 && TCut("l1l2Pt > 50") && TCut("nleps==2") && TCut("met > 40") && TCut("l1l2M>83.5 && l1l2M<98.5");

  hs_data = dataDriver.getSimpleStackTH1F("met", ";MET [GeV]; events / 10 GeV", 35, 0, 350, mySelection_SF ); 
  hs_mc   = mcDriver.getSimpleStackTH1F  ("met", ";MET [GeV]; events / 10 GeV", 35, 0, 350, mySelection_SF ); 
  basicDataMCPlot("met_ZZsel_SF", hs_data, hs_mc, true, "TR");

  hs_data = dataDriver.getSimpleStackTH1F("l1l2M", ";dilepton mass [GeV]; events / 5 GeV", 80, 20, 420, mySelection_SF ); 
  hs_mc   = mcDriver.getSimpleStackTH1F  ("l1l2M", ";dilepton mass [GeV]; events / 5 GeV", 80, 20, 420, mySelection_SF ); 
  basicDataMCPlot("l1l2M_ZZsel_SF", hs_data, hs_mc, true, "TR");

  hs_data = dataDriver.getSimpleStackTH1F("jzb", ";JZB [GeV]; events / 10 GeV", 60, -300, 300, mySelection_SF ); 
  hs_mc   = mcDriver.getSimpleStackTH1F  ("jzb", ";JZB [GeV]; events / 10 GeV", 60, -300, 300, mySelection_SF ); 
  basicDataMCPlot("jzb_ZZsel_SF", hs_data, hs_mc, true, "TL");

}


void basicDataMCPlot(string title, SimpleStack *hs_data, SimpleStack *hs_mc, bool setLog, string legPos)
{
  cout << "basicDataMCPlot: "+title << endl;
  SimpleCanvas *simpleCan = new SimpleCanvas(title.c_str(), 2);
  simpleCan->Up();
  TH1F *hall_mc   = mcDriver.getHistoTH1F(hs_mc);
  TH1F *hall_data = dataDriver.getHistoTH1F(hs_data);
  simpleCan->ShapeMeUp(hall_data);
  simpleCan->ShapeMeUp(hall_mc);
  hall_data->Draw("e1");
  hall_mc->Draw("hist same");
  hs_mc  ->Draw("hist same");
  hall_mc->Draw("axis same");
  hall_mc->Draw("hist same");
  hs_data->Draw("same e1");
  simpleCan->CMSPre();
  if(setLog) simpleCan->SetLogy();
  SimpleLegend *sleg = new SimpleLegend(legPos.c_str());
  sleg->FillLegend(hs_data);
  sleg->FillLegend(hs_mc);
  sleg->Draw("same");
  simpleCan->Dw();
  TH1F *hratio = doRatio(hall_data, hall_mc);
  simpleCan->ShapeMeDw(hratio);
  hratio->GetYaxis()->SetTitle("Data / MC");
  hratio->Draw("e1");
  hratio->GetYaxis()->SetRangeUser(0.0, 2.0);
  hratio->GetYaxis()->SetNdivisions(507);
  simpleCan->Save(outputDir);
}

/*

  h2_mc   = mcDriver.getHistoTH2F  ("l1l2Pt:l1l2M", ";dilepton mass [GeV];  dilepton p_{T} [GeV]", 28, 20, 160, 50, 0, 500, mySelection); 
  string title = "l1l2Pt_l1l2M_CE_SF_FSSub";
  cout << "basic 1 panel plot:"+title << endl;
  SimpleCanvas *simpleCan = new SimpleCanvas(title.c_str(), 1);
  h2_mc->Draw("colz");
  simpleCan->ShapeMe(h2_mc);
  simpleCan->SetLogy();
  simpleCan->CMSPhys14();
  simpleCan->Save(outputDir);
*/
