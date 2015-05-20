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

bool goFast(true);
unsigned int totEvents;

SimpleCanvas *simpleCan;

SimpleDriver mcDriver, dataDriver;
SimpleStack *hs_mc, *hs_data;

void basicDataMCPlot(string title, SimpleStack *hs_data, SimpleStack *hs_mc, bool setLog = true , string legPos = "TR");

void example()
{
  setTDRStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.5);

  dataDriver.push_back(new SimpleSample(fp_DYJetsToLL, "Data", TCut("0"), goFast));   // dummy file for pre-DATA period

  TCut Lumi("1000"); // 1000 pb-1

  mcDriver.push_back(new SimpleSample(fp_ZZ4L          , "ZZ(4l)"        , xs_ZZ4L*Lumi                                          ,goFast, kGray,      kGray));
  mcDriver.push_back(new SimpleSample(fp_WZJetsTo3LNu  , "WZ"            , xs_WZJetsTo3LNu*Lumi                                  ,goFast, 32,            32));
  mcDriver.push_back(new SimpleSample(fp_TTJets        , "t#bar{t}"      , xs_TTJets*Lumi                                        ,goFast, 40,            40)); 
  mcDriver.push_back(new SimpleSample(fp_DYJetsToLL    , "DY(#tau#tau)"  , TCut("isDYTauTau ? 1:0")*xs_DYJetsToLL*Lumi           ,goFast, 47,            47)); 
  mcDriver.push_back(new SimpleSample(fp_DYJetsToLL    , "DY(#mu#mu,ee)" , TCut("!isDYTauTau ? 1:0")*xs_DYJetsToLL*Lumi          ,goFast, kWhite          )); 
  
  TCut mySelection = sel_basic && sel_CE && sel_SF;

  hs_data = dataDriver.getSimpleStackTH1F("l1l2M", ";dilepton mass [GeV]; events / 5 GeV", 28, 20, 160, mySelection ); 
  hs_mc   = mcDriver.getSimpleStackTH1F  ("l1l2M", ";dilepton mass [GeV]; events / 5 GeV", 28, 20, 160, mySelection ); 
  basicDataMCPlot("l1l2M_basic_ce_sf", hs_data, hs_mc, true, "TR");

  hs_data = dataDriver.getSimpleStackTH1F("jzb", ";JZB [GeV]; events / 10 GeV", 60, -300, 300, mySelection ); 
  hs_mc   = mcDriver.getSimpleStackTH1F  ("jzb", ";JZB [GeV]; events / 10 GeV", 60, -300, 300, mySelection ); 
  basicDataMCPlot("jzb_basic_ce_sf", hs_data, hs_mc, true, "TL");

  hs_data = dataDriver.getSimpleStackTH1F("t1met", ";MET [GeV]; events / 10 GeV", 25, 0, 250, mySelection ); 
  hs_mc   = mcDriver.getSimpleStackTH1F  ("t1met", ";MET [GeV]; events / 10 GeV", 25, 0, 250, mySelection ); 
  basicDataMCPlot("t1met_basic_ce_sf", hs_data, hs_mc, true, "TR");
}


void basicDataMCPlot(string title, SimpleStack *hs_data, SimpleStack *hs_mc, bool setLog, string legPos)
{
  cout << "basicDataMCPlot: "+title << endl;
  SimpleCanvas *simpleCan = new SimpleCanvas(title.c_str(), 2);
  simpleCan->Up();
  TH1F *hall_mc   = mcDriver.getHistoTH1F(hs_mc);
  TH1F *hall_data = dataDriver.getHistoTH1F(hs_data);
  simpleCan->ShapeMeUp(hall_mc);
  hall_mc->Draw("hist");
  hs_mc  ->Draw("hist same");
  hall_mc->Draw("axis same");
  hall_mc->Draw("hist same");
  hs_data->Draw("same e1");
  simpleCan->CMSPhys14();
  if(setLog) simpleCan->SetLogy();
  SimpleLegend *sleg = new SimpleLegend(legPos.c_str());
  sleg->FillLegend(hs_data);
  sleg->FillLegend(hs_mc);
  sleg->Draw("same");
  simpleCan->Dw();
  TH1F *hratio = doRatio(hall_mc, hall_mc);
  simpleCan->ShapeMeDw(hratio);
  hratio->GetYaxis()->SetTitle("Data / MC");
  hratio->Draw("axis");
  hratio->GetYaxis()->SetRangeUser(0.6,1.4);
  hratio->GetYaxis()->SetNdivisions(507);
  simpleCan->Save(outputDir);
}
