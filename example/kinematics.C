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

void basicDataMCPlot(string title, SimpleStack *hs_data, SimpleStack *hs_mc, bool setLog = true , string legPos = "TR", string headerTitle = "header");


void kinematics()
{
  setTDRStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.5);

  v1_drivers(goFast);
  
  TCut mySelection_SF   = sel_basic && sel_SF && sel_SF_trig;
  TCut mySelection_ee   = sel_basic && sel_ee && sel_ee_trig;
  TCut mySelection_mumu = sel_basic && sel_mumu && sel_mumu_trig;
  TCut mySelection_emu  = sel_basic && sel_emu && sel_emu_trig;

/*
  // mySelection_ee
  {
    string headerTitle = "e^{+}e^{-}";
    hs_data = dataDriver.getSimpleStackTH1F("lepPt[0]", ";P_{T}(el1) [GeV]; events / 10 GeV", 35, 0, 350, mySelection_ee ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepPt[0]", ";P_{T}(el1) [GeV]; events / 10 GeV", 35, 0, 350, mySelection_ee ); 
    basicDataMCPlot("lepPt1_basic_ee", hs_data, hs_mc, true, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("lepPt[1]", ";P_{T}(el2) [GeV]; events / 10 GeV", 35, 0, 350, mySelection_ee ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepPt[1]", ";P_{T}(el2) [GeV]; events / 10 GeV", 35, 0, 350, mySelection_ee ); 
    basicDataMCPlot("lepPt2_basic_ee", hs_data, hs_mc, true, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("lepEta[0]", ";#eta(el1); events / 0.1", 50, -2.5, 2.5, mySelection_ee ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepEta[0]", ";#eta(el1); events / 0.1", 50, -2.5, 2.5, mySelection_ee ); 
    basicDataMCPlot("lepEta1_basic_ee", hs_data, hs_mc, false, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("lepEta[1]", ";#eta(el2); events / 0.1", 50, -2.5, 2.5, mySelection_ee ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepEta[1]", ";#eta(el2); events / 0.1", 50, -2.5, 2.5, mySelection_ee ); 
    basicDataMCPlot("lepEta2_basic_ee", hs_data, hs_mc, false, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("lepPhi[0]", ";#phi(el1); events / 0.2", 32, -3.2, 3.2, mySelection_ee ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepPhi[0]", ";#phi(el1); events / 0.2", 32, -3.2, 3.2, mySelection_ee ); 
    basicDataMCPlot("lepPhi1_basic_ee", hs_data, hs_mc, false, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("lepPhi[1]", ";#phi(el2); events / 0.2", 32, -3.2, 3.2, mySelection_ee ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepPhi[1]", ";#phi(el2); events / 0.2", 32, -3.2, 3.2, mySelection_ee ); 
    basicDataMCPlot("lepPhi2_basic_ee", hs_data, hs_mc, false, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("met", ";MET [GeV]; events / 10 GeV", 35, 0, 350, mySelection_ee ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("met", ";MET [GeV]; events / 10 GeV", 35, 0, 350, mySelection_ee ); 
    basicDataMCPlot("met_basic_ee", hs_data, hs_mc, true, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("l1l2M", ";mass(el1, el2) [GeV]; events / 5 GeV", 60, 20, 320, mySelection_ee ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("l1l2M", ";mass(el1, el2) [GeV]; events / 5 GeV", 60, 20, 320, mySelection_ee ); 
    basicDataMCPlot("l1l2M_basic_ee", hs_data, hs_mc, true, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("l1l2Pt", ";P_{T}(el1, el2) [GeV]; events / 10 GeV", 50, 0, 500, mySelection_ee ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("l1l2Pt", ";P_{T}(el1, el2) [GeV]; events / 10 GeV", 50, 0, 500, mySelection_ee ); 
    basicDataMCPlot("l1l2Pt_basic_ee", hs_data, hs_mc, true, "TR", headerTitle );
  }

  // mySelection_mumu
  {
    string headerTitle = "#mu^{+}#mu^{-}";
    hs_data = dataDriver.getSimpleStackTH1F("lepPt[0]", ";P_{T}(mu1) [GeV]; events / 10 GeV", 35, 0, 350, mySelection_mumu ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepPt[0]", ";P_{T}(mu1) [GeV]; events / 10 GeV", 35, 0, 350, mySelection_mumu ); 
    basicDataMCPlot("lepPt1_basic_mumu", hs_data, hs_mc, true, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("lepPt[1]", ";P_{T}(mu2) [GeV]; events / 10 GeV", 35, 0, 350, mySelection_mumu ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepPt[1]", ";P_{T}(mu2) [GeV]; events / 10 GeV", 35, 0, 350, mySelection_mumu ); 
    basicDataMCPlot("lepPt2_basic_mumu", hs_data, hs_mc, true, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("lepEta[0]", ";#eta(mu1); events / 0.1", 50, -2.5, 2.5, mySelection_mumu ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepEta[0]", ";#eta(mu1); events / 0.1", 50, -2.5, 2.5, mySelection_mumu ); 
    basicDataMCPlot("lepEta1_basic_mumu", hs_data, hs_mc, false, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("lepEta[1]", ";#eta(mu2); events / 0.1", 50, -2.5, 2.5, mySelection_mumu ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepEta[1]", ";#eta(mu2); events / 0.1", 50, -2.5, 2.5, mySelection_mumu ); 
    basicDataMCPlot("lepEta2_basic_mumu", hs_data, hs_mc, false, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("lepPhi[0]", ";#phi(mu1); events / 0.2", 32, -3.2, 3.2, mySelection_mumu ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepPhi[0]", ";#phi(mu1); events / 0.2", 32, -3.2, 3.2, mySelection_mumu ); 
    basicDataMCPlot("lepPhi1_basic_mumu", hs_data, hs_mc, false, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("lepPhi[1]", ";#phi(mu2); events / 0.2", 32, -3.2, 3.2, mySelection_mumu ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepPhi[1]", ";#phi(mu2); events / 0.2", 32, -3.2, 3.2, mySelection_mumu ); 
    basicDataMCPlot("lepPhi2_basic_mumu", hs_data, hs_mc, false, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("met", ";MET [GeV]; events / 10 GeV", 35, 0, 350, mySelection_mumu ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("met", ";MET [GeV]; events / 10 GeV", 35, 0, 350, mySelection_mumu ); 
    basicDataMCPlot("met_basic_mumu", hs_data, hs_mc, true, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("l1l2M", ";mass(mu1, mu2) [GeV]; events / 5 GeV", 60, 20, 320, mySelection_mumu ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("l1l2M", ";mass(mu1, mu2) [GeV]; events / 5 GeV", 60, 20, 320, mySelection_mumu ); 
    basicDataMCPlot("l1l2M_basic_mumu", hs_data, hs_mc, true, "TR", headerTitle );
  
    hs_data = dataDriver.getSimpleStackTH1F("l1l2Pt", ";P_{T}(mu1, mu2) [GeV]; events / 10 GeV", 50, 0, 500, mySelection_mumu ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("l1l2Pt", ";P_{T}(mu1, mu2) [GeV]; events / 10 GeV", 50, 0, 500, mySelection_mumu ); 
    basicDataMCPlot("l1l2Pt_basic_mumu", hs_data, hs_mc, true, "TR", headerTitle );
  }
*/
  // mySelection_emu
  {
    string headerTitle = "e^{#pm}#mu^{#mp}";
    hs_data = dataDriver.getSimpleStackTH1F("lepPt[0]", ";P_{T}(l1) [GeV]; events / 10 GeV", 35, 0, 350, mySelection_emu ); 
    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepPt[0]", ";P_{T}(l1) [GeV]; events / 10 GeV", 35, 0, 350, mySelection_emu ); 
    basicDataMCPlot("lepPt1_basic_emu", hs_data, hs_mc, true, "TR", headerTitle );
  
//    hs_data = dataDriver.getSimpleStackTH1F("lepPt[1]", ";P_{T}(l2) [GeV]; events / 10 GeV", 35, 0, 350, mySelection_emu ); 
//    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepPt[1]", ";P_{T}(l2) [GeV]; events / 10 GeV", 35, 0, 350, mySelection_emu ); 
//    basicDataMCPlot("lepPt2_basic_emu", hs_data, hs_mc, true, "TR", headerTitle );
//  
//    hs_data = dataDriver.getSimpleStackTH1F("lepEta[0]", ";#eta(l1); events / 0.1", 50, -2.5, 2.5, mySelection_emu ); 
//    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepEta[0]", ";#eta(l1); events / 0.1", 50, -2.5, 2.5, mySelection_emu ); 
//    basicDataMCPlot("lepEta1_basic_emu", hs_data, hs_mc, false, "TR", headerTitle );
//  
//    hs_data = dataDriver.getSimpleStackTH1F("lepEta[1]", ";#eta(l2); events / 0.1", 50, -2.5, 2.5, mySelection_emu ); 
//    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepEta[1]", ";#eta(l2); events / 0.1", 50, -2.5, 2.5, mySelection_emu ); 
//    basicDataMCPlot("lepEta2_basic_emu", hs_data, hs_mc, false, "TR", headerTitle );
//  
//    hs_data = dataDriver.getSimpleStackTH1F("lepPhi[0]", ";#phi(l1); events / 0.2", 32, -3.2, 3.2, mySelection_emu ); 
//    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepPhi[0]", ";#phi(l1); events / 0.2", 32, -3.2, 3.2, mySelection_emu ); 
//    basicDataMCPlot("lepPhi1_basic_emu", hs_data, hs_mc, false, "TR", headerTitle );
//  
//    hs_data = dataDriver.getSimpleStackTH1F("lepPhi[1]", ";#phi(l2); events / 0.2", 32, -3.2, 3.2, mySelection_emu ); 
//    hs_mc   = mcDriver.getSimpleStackTH1F  ("lepPhi[1]", ";#phi(l2); events / 0.2", 32, -3.2, 3.2, mySelection_emu ); 
//    basicDataMCPlot("lepPhi2_basic_emu", hs_data, hs_mc, false, "TR", headerTitle );
//  
//    hs_data = dataDriver.getSimpleStackTH1F("met", ";MET [GeV]; events / 10 GeV", 35, 0, 350, mySelection_emu ); 
//    hs_mc   = mcDriver.getSimpleStackTH1F  ("met", ";MET [GeV]; events / 10 GeV", 35, 0, 350, mySelection_emu ); 
//    basicDataMCPlot("met_basic_emu", hs_data, hs_mc, true, "TR", headerTitle );
//  
//    hs_data = dataDriver.getSimpleStackTH1F("l1l2M", ";mass(l1, l2) [GeV]; events / 5 GeV", 60, 20, 320, mySelection_emu ); 
//    hs_mc   = mcDriver.getSimpleStackTH1F  ("l1l2M", ";mass(l1, l2) [GeV]; events / 5 GeV", 60, 20, 320, mySelection_emu ); 
//    basicDataMCPlot("l1l2M_basic_emu", hs_data, hs_mc, true, "TR", headerTitle );
//  
//    hs_data = dataDriver.getSimpleStackTH1F("l1l2Pt", ";P_{T}(l1, l2) [GeV]; events / 10 GeV", 50, 0, 500, mySelection_emu ); 
//    hs_mc   = mcDriver.getSimpleStackTH1F  ("l1l2Pt", ";P_{T}(l1, l2) [GeV]; events / 10 GeV", 50, 0, 500, mySelection_emu ); 
//    basicDataMCPlot("l1l2Pt_basic_emu", hs_data, hs_mc, true, "TR", headerTitle );
  }


}


void basicDataMCPlot(string title, SimpleStack *hs_data, SimpleStack *hs_mc, bool setLog, string legPos, string headerTitle)
{
  cout << "basicDataMCPlot: "+title << endl;
  SimpleCanvas *simpleCan = new SimpleCanvas(title.c_str(), 2);
  simpleCan->Up();
  TH1F *hall_mc   = mcDriver.getHistoTH1F(hs_mc);
  TH1F *hall_data = dataDriver.getHistoTH1F(hs_data);
  simpleCan->ShapeMeUp(hall_data);
  simpleCan->ShapeMeUp(hall_mc);
  hall_data->Sumw2(kFALSE); hall_data->SetBinErrorOption(TH1::kPoisson);
  hall_data->Draw("e1");
  hall_mc->Draw("hist same");
  hs_mc  ->Draw("hist same");
  hall_mc->Draw("axis same");
  hall_mc->Draw("hist same");
  hall_data->Draw("same e1");
  simpleCan->CMSPre();
  if(setLog) simpleCan->SetLogy();
  SimpleLegend *sleg = new SimpleLegend(legPos.c_str());
  sleg->SetHeader(headerTitle.c_str());
//  sleg->SetHeader("dielectron");
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

