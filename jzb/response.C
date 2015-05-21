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

bool goFast(false);
unsigned int totEvents;

SimpleCanvas *simpleCan, *sc1, *sc2;
SimpleLegend *simpleLeg, *sl1, *sl2;
string can_title;

SimpleDriver mcDriver, dataDriver;

enum jetBins{ij0_CE, ij0_FW, ij0, ej0, ej1, ej2, ij3}myJetBins;
TProfile *pro[10];
vector<float> myBinsX[10];

TFile *fpout;

void response()
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
  

  string filename = (outputDir + "/jzb_resp_mc_Phys14.root").c_str();
  fpout = new TFile(filename.c_str() ,"RECREATE");


  TCut sel_massBin = sel_M81101 && TCut("t1met<50"); 

  // --- build first canvas
  myBinsX[ij0_CE] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,100,200,500};
  myBinsX[ij0_FW] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,100,200,500};
  pro[ij0_CE] = mcDriver.getProfileX("t1vHT/l1l2Pt:l1l2Pt", ";dilepton pt [GeV]; jzb response", myBinsX[ij0_CE], (sel_basic && sel_massBin && sel_ij0 && sel_CE)*sel_FSSub);
  pro[ij0_FW] = mcDriver.getProfileX("t1vHT/l1l2Pt:l1l2Pt", ";dilepton pt [GeV]; jzb response", myBinsX[ij0_FW], (sel_basic && sel_massBin && sel_ij0 && sel_FW)*sel_FSSub);
  pro[ij0_FW]->SetLineColor(kRed);
  pro[ij0_FW]->SetLineWidth(2);

  can_title = "resp_summary_ij0_CE_FW";
  sc1 = new SimpleCanvas(can_title.c_str(), 1);
  pro[ij0_CE]->Draw("e1");
  pro[ij0_FW]->Draw("hist same");
  pro[ij0_CE]->Draw("e1 same");
  sc1->ShapeMe(pro[ij0_CE]);
  sc1->ShapeMe(pro[ij0_FW]);
  sl1 = new SimpleLegend();
  sl1->AddEntry(pro[ij0_CE],"ij0_CE","ple");
  sl1->AddEntry(pro[ij0_FW],"ij0_FW","l");
  sl1->Draw("same");
  sc1->SetLogy();
  sc1->CMSPhys14();
  sc1->Save(outputDir);
  
  pro[ij0_CE]->Write();
  pro[ij0_FW]->Write();
  sc1->can_->Write();
  // --- build second canvas 
  myBinsX[ij0] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,100,200,500};
  myBinsX[ej0] = {0,5,10,15,20,25,30,35,40,70,100,500};
  myBinsX[ej1] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,100,200,500};
  myBinsX[ej2] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,100,200,500};
  myBinsX[ij3] = {0,20,40,50,60,65,70,75,80,85,90,100,200,500};

  pro[ij0] = mcDriver.getProfileX("t1vHT/l1l2Pt:l1l2Pt", ";dilepton pt [GeV]; jzb response", myBinsX[ij0], (sel_basic && sel_massBin && sel_ij0)*sel_FSSub);
  pro[ej0] = mcDriver.getProfileX("t1vHT/l1l2Pt:l1l2Pt", ";dilepton pt [GeV]; jzb response", myBinsX[ej0], (sel_basic && sel_massBin && sel_ej0)*sel_FSSub);
  pro[ej1] = mcDriver.getProfileX("t1vHT/l1l2Pt:l1l2Pt", ";dilepton pt [GeV]; jzb response", myBinsX[ej1], (sel_basic && sel_massBin && sel_ej1)*sel_FSSub);
  pro[ej2] = mcDriver.getProfileX("t1vHT/l1l2Pt:l1l2Pt", ";dilepton pt [GeV]; jzb response", myBinsX[ej2], (sel_basic && sel_massBin && sel_ej2)*sel_FSSub);
  pro[ij3] = mcDriver.getProfileX("t1vHT/l1l2Pt:l1l2Pt", ";dilepton pt [GeV]; jzb response", myBinsX[ij3], (sel_basic && sel_massBin && sel_ij3)*sel_FSSub);

  pro[ej0]->SetLineColor(kRed);
  pro[ej0]->SetLineWidth(2);

  pro[ej1]->SetLineColor(39);
  pro[ej1]->SetLineWidth(2);

  pro[ej2]->SetLineColor(kGreen-2);
  pro[ej2]->SetLineWidth(2);

  pro[ij3]->SetLineColor(9);
  pro[ij3]->SetLineWidth(2);

  can_title = "resp_summary_njets";
  sc2 = new SimpleCanvas(can_title.c_str(), 1);
  pro[ij0]->GetXaxis()->SetNdivisions(509);
  pro[ij0]->Draw("e1");
  pro[ej0]->Draw("hist same");
  pro[ej1]->Draw("hist same");
  pro[ej2]->Draw("hist same");
  pro[ij3]->Draw("hist same");
  pro[ij0]->Draw("e1 same");
  sc2->ShapeMe(pro[ij0]);
  sl2 = new SimpleLegend();
  sl2->AddEntry(pro[ij0],"ij0","ple");
  sl2->AddEntry(pro[ej0],"ej0","l");
  sl2->AddEntry(pro[ej1],"ej1","l");
  sl2->AddEntry(pro[ej2],"ej2","l");
  sl2->AddEntry(pro[ij3],"ij3","l");
  sl2->Draw("same");
  sc2->SetLogy();
  sc2->CMSPhys14();
  sc2->Save(outputDir);

  pro[ij0]->Write();
  pro[ej0]->Write();
  pro[ej1]->Write();
  pro[ej2]->Write();
  pro[ij3]->Write();

  sc2->can_->Write();
  fpout->Write();
}


