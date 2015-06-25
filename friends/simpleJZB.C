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

unsigned int totEvents;


SimpleDriver mcDriver, dataDriver;

using namespace std;
bool goFast(false);


void simpleJZB()
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

  // needs to become a function
  TCut sel_cut_SF      = sel_basic && sel_M81101 && sel_ej2 && sel_CE && sel_SF;
  TCut sel_cut_OF      = sel_basic && sel_M81101 && sel_ej2 && sel_CE && sel_OF;
  TCut sel_cut_jzb_pos = TCut("t1vHT-l1l2Pt>0");  
  TCut sel_cut_jzb_neg = TCut("t1vHT-l1l2Pt<0");  
  string metvar        = "t1met";
  string plotTitle     = "_ej2_CE_M81101";
  pnumber OF_scale     = pnumber(1, 0.05);

/*
  // --- plot MET all
  {
      char SF_header[] = "Same Flavor";
      char OF_header[] = "Opposite Flavor";
      string title     = "";
       
      title = metvar+plotTitle+"_SF";
      SimpleStack * sstack_SF    = mcDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",50,0,500, sel_cut_SF);
      TH1F *sstack_SF_hall       = mcDriver.getHistoTH1F(sstack_SF);	
      SimpleLegend *sleg_SF = new SimpleLegend("TRSF");
      SimpleCanvas *simpleCan_SF = new SimpleCanvas(title.c_str(), 1);
      sstack_SF_hall->Draw("hist");
      sstack_SF->Draw("hist same");
      sstack_SF_hall->Draw("hist same");
      sstack_SF_hall->Draw("axis same");
      simpleCan_SF->CMSPhys14();
      simpleCan_SF->SetLogy();
      sleg_SF->SetHeader(SF_header);
      sleg_SF->FillLegend(sstack_SF);
      sleg_SF->Draw("same");
      simpleCan_SF->Save(outputDir);

      title = metvar+plotTitle+"_OF";
      SimpleStack * sstack_OF    = mcDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",50,0,500, sel_cut_OF);
      TH1F *sstack_OF_hall       = mcDriver.getHistoTH1F(sstack_OF);	
      SimpleLegend *sleg_OF = new SimpleLegend("TROF");
      SimpleCanvas *simpleCan_OF = new SimpleCanvas(title.c_str(), 1);
      sstack_SF_hall->Draw("axis");
      sstack_OF_hall->Draw("hist same");
      sstack_OF->Draw("hist same");
      sstack_OF_hall->Draw("hist same");
      sstack_OF_hall->Draw("axis same");
      simpleCan_OF->CMSPhys14();
      simpleCan_OF->SetLogy();
      sleg_OF->SetHeader(OF_header);
      sleg_OF->FillLegend(sstack_OF);
      sleg_OF->Draw("same");
      simpleCan_OF->Save(outputDir);
  }
*/

  // --- plot MET for jzb_pos and jzb_neg
  {
      char SF_jzb_pos_header[]    = "#splitline{Same Flavor}{JZB > 0}";
      char OF_jzb_pos_header[]    = "#splitline{Opposite Flavor}{JZB > 0}";
      char SF_jzb_neg_header[]    = "#splitline{Same Flavor}{JZB < 0}";
      char OF_jzb_neg_header[]    = "#splitline{Opposite Flavor}{JZB < 0}";
      string title        = "";
      
      title = metvar+plotTitle+"_SF"+"_jzb_pos";
      SimpleStack * sstack_SF_jzb_pos    = mcDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",50,0,500, sel_cut_SF && sel_cut_jzb_pos);
      TH1F *sstack_SF_jzb_pos_hall       = mcDriver.getHistoTH1F(sstack_SF_jzb_pos);	
      SimpleLegend *sleg_SF_jzb_pos = new SimpleLegend("TRSFjzb_pos");
      SimpleCanvas *simpleCan_SF_jzb_pos = new SimpleCanvas(title.c_str(), 1);
      sstack_SF_jzb_pos_hall->Draw("hist");
      sstack_SF_jzb_pos->Draw("hist same");
      sstack_SF_jzb_pos_hall->Draw("hist same");
      sstack_SF_jzb_pos_hall->Draw("axis same");
      simpleCan_SF_jzb_pos->CMSPhys14();
      simpleCan_SF_jzb_pos->SetLogy();
      sleg_SF_jzb_pos->SetHeader(SF_jzb_pos_header);
      sleg_SF_jzb_pos->FillLegend(sstack_SF_jzb_pos);
      sleg_SF_jzb_pos->Draw("same");
      simpleCan_SF_jzb_pos->Save(outputDir);

      title = metvar+plotTitle+"_OF"+"_jzb_pos";
      SimpleStack * sstack_OF_jzb_pos    = mcDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",50,0,500, sel_cut_OF && sel_cut_jzb_pos);
      TH1F *sstack_OF_jzb_pos_hall       = mcDriver.getHistoTH1F(sstack_OF_jzb_pos);	
      SimpleLegend *sleg_OF_jzb_pos = new SimpleLegend("TROFjzb_pos");
      SimpleCanvas *simpleCan_OF_jzb_pos = new SimpleCanvas(title.c_str(), 1);
      sstack_SF_jzb_pos_hall->Draw("axis");
      sstack_OF_jzb_pos_hall->Draw("hist same");
      sstack_OF_jzb_pos->Draw("hist same");
      sstack_OF_jzb_pos_hall->Draw("hist same");
      sstack_OF_jzb_pos_hall->Draw("axis same");
      simpleCan_OF_jzb_pos->CMSPhys14();
      simpleCan_OF_jzb_pos->SetLogy();
      sleg_OF_jzb_pos->SetHeader(OF_jzb_pos_header);
      sleg_OF_jzb_pos->FillLegend(sstack_OF_jzb_pos);
      sleg_OF_jzb_pos->Draw("same");
      simpleCan_OF_jzb_pos->Save(outputDir);


      title = metvar+plotTitle+"_SF"+"_jzb_neg";
      SimpleStack * sstack_SF_jzb_neg    = mcDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",50,0,500, sel_cut_SF && sel_cut_jzb_neg);
      TH1F *sstack_SF_jzb_neg_hall       = mcDriver.getHistoTH1F(sstack_SF_jzb_neg);	
      SimpleLegend *sleg_SF_jzb_neg = new SimpleLegend("TRSFjzb_neg");
      SimpleCanvas *simpleCan_SF_jzb_neg = new SimpleCanvas(title.c_str(), 1);
      sstack_SF_jzb_neg_hall->Draw("hist");
      sstack_SF_jzb_neg->Draw("hist same");
      sstack_SF_jzb_neg_hall->Draw("hist same");
      sstack_SF_jzb_neg_hall->Draw("axis same");
      simpleCan_SF_jzb_neg->CMSPhys14();
      simpleCan_SF_jzb_neg->SetLogy();
      sleg_SF_jzb_neg->SetHeader(SF_jzb_neg_header);
      sleg_SF_jzb_neg->FillLegend(sstack_SF_jzb_neg);
      sleg_SF_jzb_neg->Draw("same");
      simpleCan_SF_jzb_neg->Save(outputDir);

      title = metvar+plotTitle+"_OF"+"_jzb_neg";
      SimpleStack * sstack_OF_jzb_neg    = mcDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",50,0,500, sel_cut_OF && sel_cut_jzb_neg);
      TH1F *sstack_OF_jzb_neg_hall       = mcDriver.getHistoTH1F(sstack_OF_jzb_neg);	
      SimpleLegend *sleg_OF_jzb_neg = new SimpleLegend("TROFjzb_neg");
      SimpleCanvas *simpleCan_OF_jzb_neg = new SimpleCanvas(title.c_str(), 1);
      sstack_SF_jzb_neg_hall->Draw("axis");
      sstack_OF_jzb_neg_hall->Draw("hist same");
      sstack_OF_jzb_neg->Draw("hist same");
      sstack_OF_jzb_neg_hall->Draw("hist same");
      sstack_OF_jzb_neg_hall->Draw("axis same");
      simpleCan_OF_jzb_neg->CMSPhys14();
      simpleCan_OF_jzb_neg->SetLogy();
      sleg_OF_jzb_neg->SetHeader(OF_jzb_neg_header);
      sleg_OF_jzb_neg->FillLegend(sstack_OF_jzb_neg);
      sleg_OF_jzb_neg->Draw("same");
      simpleCan_OF_jzb_neg->Save(outputDir);

      pnumber jzb_pos_0_50_SF = histIntegralPN(sstack_SF_jzb_pos_hall, 0 ,50);
      pnumber jzb_neg_0_50_SF = histIntegralPN(sstack_SF_jzb_neg_hall, 0 ,50);
      pnumber jzb_pos_0_50_OF = histIntegralPN(sstack_OF_jzb_pos_hall, 0 ,50);
      pnumber jzb_neg_0_50_OF = histIntegralPN(sstack_OF_jzb_neg_hall, 0 ,50);

      pnumber jzb_norm = (jzb_pos_0_50_SF-jzb_pos_0_50_OF)/(jzb_neg_0_50_SF-jzb_neg_0_50_OF);
      cout << "jzb_pos_0_50_SF/jzb_neg_0_50_SF = " << jzb_pos_0_50_SF/jzb_neg_0_50_SF << endl;
      cout << "jzb_norm = " << jzb_norm << endl;

      vector<TH1F*> SF_jzb_pos_TH1F_vec = sstack_SF_jzb_pos->histoPointers_;
      vector<TH1F*> SF_jzb_neg_TH1F_vec = sstack_SF_jzb_neg->histoPointers_;
      vector<TH1F*> OF_jzb_pos_TH1F_vec = sstack_OF_jzb_pos->histoPointers_;
      vector<TH1F*> OF_jzb_neg_TH1F_vec = sstack_OF_jzb_neg->histoPointers_;
    
      vector<string> sampleTitles = sstack_SF_jzb_pos->sampleTitles_;

      for(size_t ii = 0; ii < sampleTitles.size(); ++ii)
      {
        cout <<"processing "<< sampleTitles[ii] << endl;

       // --- assumes the OF_scale is uncorrelated in JZB > 0, JZB < 0
       TH1F *SF_jzb_pos = (TH1F*)SF_jzb_pos_TH1F_vec[ii]->Clone(); 
       TH1F *SF_jzb_neg = (TH1F*)SF_jzb_neg_TH1F_vec[ii]->Clone();
       TH1F *OF_jzb_neg = ScaleTH1F(OF_jzb_neg_TH1F_vec[ii], OF_scale);
       TH1F *OF_jzb_pos = ScaleTH1F(OF_jzb_pos_TH1F_vec[ii], OF_scale);
       TH1F *jzb_only = (TH1F*)SF_jzb_neg->Clone();
       jzb_only ->Add(OF_jzb_neg,-1);
       TH1F *jzb_scaled = ScaleTH1F(jzb_only, jzb_norm);
       TH1F *pred = jzb_scaled;
       pred->Add(OF_jzb_pos);

        // --- assumes that all scales are uncorrelated
	//TH1F *SF_jzb_pos = ScaleTH1F(SF_jzb_pos_TH1F_vec[ii], pnumber(1,0));
	//TH1F *SF_jzb_neg = ScaleTH1F(SF_jzb_neg_TH1F_vec[ii], pnumber(1,0)*jzb_norm);
	//TH1F *OF_jzb_neg = ScaleTH1F(OF_jzb_neg_TH1F_vec[ii], OF_scale*pnumber(-1,0)*jzb_norm);
    	//TH1F *OF_jzb_pos = ScaleTH1F(OF_jzb_pos_TH1F_vec[ii], OF_scale);
        //TH1F *pred = (TH1F*)SF_jzb_neg->Clone();
        //pred->Add(OF_jzb_neg);
        //pred->Add(OF_jzb_pos);

        SF_jzb_pos->SetLineColor(kBlack);
        SF_jzb_pos->SetFillColor(kWhite);
        pred->SetLineColor(kRed);
        pred->SetFillColor(kWhite);
        pred->SetLineWidth(2);
        pred->SetLineStyle(2);

        title = metvar+plotTitle+"_pred"+any2string(ii);  
        SimpleLegend *sleg      = new SimpleLegend("jzb_exp_mode");
        SimpleCanvas *simpleCan = new SimpleCanvas(title.c_str(), 2);
        simpleCan->CMSPhys14();
	simpleCan->SetLogy();
        simpleCan->Up();
        simpleCan->ShapeMeUp(SF_jzb_pos); 
        simpleCan->ShapeMeUp(pred); 
        SF_jzb_pos->Draw("hist");
        pred->Draw("hist same"); 
        sleg->SetHeader((sampleTitles[ii]+" [SF]").c_str());
        sleg->AddEntry(SF_jzb_pos,"JZB > 0","FL");
        sleg->AddEntry(pred, "prediction", "FL");
        sleg->Draw("same");
        simpleCan->Dw();
        TH1F *hratio = doRatio(SF_jzb_pos, pred);
        hratio->GetYaxis()->SetTitle("ratio");
        hratio->GetYaxis()->CenterTitle();
        simpleCan->ShapeMeDw(hratio);
        hratio->Draw("e1");
        hratio->GetYaxis()->SetRangeUser(0.6,1.4);
        hratio->GetYaxis()->SetNdivisions(507);
        simpleCan->Save(outputDir);

      }	
  }

}


