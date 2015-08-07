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

using namespace std;
bool goFast(false);

unsigned int totEvents;
SimpleDriver mcDriver, dataDriver;


void plotMET_datamc();


// needs to become a function
string v_met = "met";
string v_hr  = "vHT";
string plotTitle     = "_ej0_CE_M81101_"+v_met+"_"+v_hr;
TCut sel_cut_gauge   = TCut("!isDYTauTau") && TCut("vHT>0 && l1l2Pt>0");
TCut sel_cut_SF      = sel_basic && sel_M81101 && sel_ej0 && sel_CE && sel_SF_trig && sel_cut_gauge;
TCut sel_cut_OF      = sel_basic && sel_M81101 && sel_ej0 && sel_CE && sel_OF_trig && sel_cut_gauge;
TCut sel_cut_jzb_pos = TCut((v_hr+"-l1l2Pt>0").c_str());  
TCut sel_cut_jzb_neg = TCut((v_hr+"-l1l2Pt<0").c_str());  
string metvar        = v_met.c_str();
string jzbvar        = (v_hr+"-l1l2Pt").c_str();
pnumber OF_scale     = pnumber(1, 0.05);

// HACK no scale factor between positive and negative


void experimentalJZB()
{ 

  setTDRStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.5);

  TCut Lumi("40.24"); // 1000 pb-1
  if(goFast) Lumi = Lumi*TCut("0.01");
  cout << "loading v1_drivers" << endl;
  cout << "goFast = " << goFast << endl;
  cout << "Lumi = " << Lumi.GetTitle() << endl;

  string myLocalPath       = "/afs/cern.ch/work/t/theofil/public/ntuples/data_v1/dileptonSkim/";
//  string myLocalPath       = "/afs/cern.ch/work/t/theofil/public/ntuples/data_v1/dileptonSkim_ej1/";
//  string myLocalPath       = "/afs/cern.ch/work/t/theofil/public/ntuples/data_v1/dileptonSkim_ij2/";

  mcDriver.push_back(new SimpleSample(myLocalPath+fname_DYJetsM50    , "DY"                    , Lumi                           ,goFast, kWhite , kRed+1         )); 

  cout << "..:: SimpleJZB log ::..." << endl;
  cout << "sel_cut_SF = " << sel_cut_SF.GetTitle() << endl;
  cout << "sel_cut_SF = " << sel_cut_OF.GetTitle() << endl;
  cout << "sel_cut_jzb_pos = " << sel_cut_jzb_pos.GetTitle() << endl;
  cout << "sel_cut_jzb_neg = " << sel_cut_jzb_neg.GetTitle() << endl;
  cout << "metvar = " << metvar << endl;
  cout << "jzbvar = " << jzbvar << endl;
  cout << "plotTitle = "<< plotTitle << endl;
  cout << "OF_scale = " << OF_scale << endl;
  cout << endl;

  // --- plot JZB all for mc
  {
      char SF_header[] = "Same Flavor";
      char OF_header[] = "Opposite Flavor";
      string title     = "";
       
      title = "jzb"+plotTitle+"_SF";
      SimpleStack * sstack_mc_SF    = mcDriver.getSimpleStackTH1F(jzbvar,";JZB [GeV]; events / 2 GeV;",60,-60,60, sel_cut_SF);
      TH1F *sstack_mc_SF_hall       = mcDriver.getHistoTH1F(sstack_mc_SF);	
      SimpleLegend *sleg_SF = new SimpleLegend("TLSF");
      SimpleCanvas *simpleCan_SF = new SimpleCanvas(title.c_str(), 1);
      sstack_mc_SF_hall->Draw("hist");
      sstack_mc_SF_hall->Draw("hist same");
      sstack_mc_SF->Draw("hist same");
      sstack_mc_SF_hall->Draw("hist same");
      sstack_mc_SF_hall->Draw("axis same");
      simpleCan_SF->CMSPre();
      simpleCan_SF->SetLogy();
      sleg_SF->SetHeader(SF_header);
      sleg_SF->FillLegend(sstack_mc_SF);
      sleg_SF->Draw("same");
      simpleCan_SF->Save(outputDir);

  }

  // --- contrast hadronic recoil versus dilepton pt
  {
      char SF_header[] = "Same Flavor";
      char OF_header[] = "Opposite Flavor";
      string title     = "";
       
      title = "hr_l1l2Pt"+plotTitle+"_SF";
      SimpleStack * sstack_hr_mc_SF    = mcDriver.getSimpleStackTH1F(v_hr,";P_{T} [GeV]; events / 5 GeV;",24,0,120, sel_cut_SF);
      TH1F *sstack_hr_mc_SF_hall       = mcDriver.getHistoTH1F(sstack_hr_mc_SF);	
      SimpleStack * sstack_l1l2Pt_mc_SF    = mcDriver.getSimpleStackTH1F("l1l2Pt",";P_{T} [GeV]; events / 5 GeV;",24,0,120, sel_cut_SF);
      TH1F *sstack_l1l2Pt_mc_SF_hall       = mcDriver.getHistoTH1F(sstack_l1l2Pt_mc_SF);	
      SimpleLegend *sleg = new SimpleLegend("TTR");
      SimpleCanvas *simpleCan = new SimpleCanvas(title.c_str(), 2);
      simpleCan->Up();
      simpleCan->ShapeMeUp(sstack_l1l2Pt_mc_SF_hall); 
      simpleCan->ShapeMeUp(sstack_hr_mc_SF_hall); 
      sstack_l1l2Pt_mc_SF_hall->SetLineColor(kBlue);
      sstack_l1l2Pt_mc_SF_hall->Draw("hist");
      sstack_hr_mc_SF_hall->Draw("hist same");
      simpleCan->CMSPre();
      simpleCan->SetLogy();
      sleg->AddEntry(sstack_l1l2Pt_mc_SF_hall, "P_{T} (l1,l2)", "F");
      sleg->AddEntry(sstack_hr_mc_SF_hall, "P_{T} (had recoil)", "F");
      sleg->Draw("same");
      simpleCan->Dw();
      //TH1F *hratio = doSimpleRatio(sstack_hr_mc_SF_hall, sstack_l1l2Pt_mc_SF_hall);
      TH1F *hratio = doRatio(sstack_hr_mc_SF_hall, sstack_l1l2Pt_mc_SF_hall, 3);
      hratio->GetYaxis()->SetTitle("ratio");
      hratio->GetYaxis()->CenterTitle();
      simpleCan->ShapeMeDw(hratio);
      hratio->SetLineColor(kBlack);
      hratio->Draw("e1");
      hratio->GetYaxis()->SetRangeUser(0.0,2.0);
      hratio->GetYaxis()->SetNdivisions(507);
      simpleCan->Save(outputDir);

  }

  // --- calculate MET for jzb_pos and jzb_neg
  {
      char SF_jzb_pos_header[]    = "#splitline{Same Flavor}{JZB > 0}";
      char SF_jzb_neg_header[]    = "#splitline{Same Flavor}{JZB < 0}";
      string title        = "";
      
      title = metvar+plotTitle+"_SF"+"_jzb_pos";
      SimpleStack * sstack_mc_SF_jzb_pos    = mcDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",30,0,300, sel_cut_SF && sel_cut_jzb_pos);
      TH1F *sstack_mc_SF_jzb_pos_hall       = mcDriver.getHistoTH1F(sstack_mc_SF_jzb_pos);	

      title = metvar+plotTitle+"_SF"+"_jzb_neg";
      SimpleStack * sstack_mc_SF_jzb_neg    = mcDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",30,0,300, sel_cut_SF && sel_cut_jzb_neg);
      TH1F *sstack_mc_SF_jzb_neg_hall       = mcDriver.getHistoTH1F(sstack_mc_SF_jzb_neg);	


      //-------------------------------------------
      //--- jzb core for mc
      //-------------------------------------------
      
      pnumber jzb_pos_0_50_mc_SF = histIntegralPN(sstack_mc_SF_jzb_pos_hall, 0 ,50);
      pnumber jzb_neg_0_50_mc_SF = histIntegralPN(sstack_mc_SF_jzb_neg_hall, 0 ,50);
      pnumber jzb_pos_0_50_mc_OF = pnumber(0,0);
      pnumber jzb_neg_0_50_mc_OF = pnumber(0,0);

      pnumber jzb_norm_mc = (jzb_pos_0_50_mc_SF-jzb_pos_0_50_mc_OF)/(jzb_neg_0_50_mc_SF-jzb_neg_0_50_mc_OF);
      cout << "jzb_pos_0_50_mc_SF/jzb_neg_0_50_mc_SF = " << jzb_pos_0_50_mc_SF/jzb_neg_0_50_mc_SF << endl;
      cout << "jzb_norm_mc = " << jzb_norm_mc << endl;

      vector<TH1F*> SF_jzb_pos_TH1F_vec_mc = sstack_mc_SF_jzb_pos->histoPointers_;
      vector<TH1F*> SF_jzb_neg_TH1F_vec_mc = sstack_mc_SF_jzb_neg->histoPointers_;
    
      vector<string> sampleTitles_mc = sstack_mc_SF_jzb_pos->sampleTitles_;
      vector<TH1F*> jzb_predictions_mc;


      TH1F *pred_all_mc = (TH1F*)sstack_mc_SF_jzb_neg_hall->Clone();
     // TH1F *jzb_only = (TH1F*)sstack_mc_SF_jzb_neg_hall->Clone();
     // TH1F *jzb_scaled = ScaleTH1F(jzb_only, jzb_norm_mc); // HACK

      {
  	title = metvar+plotTitle+"_pred_all_mc";  
  	SimpleLegend *sleg      = new SimpleLegend("jzb_exp_mode");
       	SimpleCanvas *simpleCan = new SimpleCanvas(title.c_str(), 2);
        simpleCan->CMSPre();
        simpleCan->SetLogy();
        simpleCan->Up();
        simpleCan->ShapeMeUp(sstack_mc_SF_jzb_pos_hall); 
        simpleCan->ShapeMeUp(pred_all_mc); 
        sstack_mc_SF_jzb_pos_hall->SetLineColor(kBlue);
        sstack_mc_SF_jzb_pos_hall->Draw("hist");
        pred_all_mc->Draw("hist same"); 
        sleg->SetHeader("Same Flavor [MC]");
        sleg->AddEntry(sstack_mc_SF_jzb_pos_hall,"JZB > 0","FL");
        sleg->AddEntry(pred_all_mc, "prediction", "FL");
        sleg->Draw("same");
        simpleCan->Dw();
        //TH1F *hratio = doSimpleRatio(sstack_mc_SF_jzb_pos_hall, pred_all_mc);
        TH1F *hratio = doRatio(sstack_mc_SF_jzb_pos_hall, pred_all_mc, 3);
        hratio->GetYaxis()->SetTitle("ratio");
        hratio->GetYaxis()->CenterTitle();
        simpleCan->ShapeMeDw(hratio);
        hratio->SetLineColor(kBlack);
        hratio->Draw("e1");
        hratio->GetYaxis()->SetRangeUser(0.0,2.0);
        hratio->GetYaxis()->SetNdivisions(507);
        simpleCan->Save(outputDir);

        // --- printout
        pnumber obs_50_100 = histIntegralPN(sstack_mc_SF_jzb_pos_hall, 50, 100);
        pnumber pred_50_100 = histIntegralPN(pred_all_mc, 50, 100);
        cout << "mc: [50, 100] "<< "all" << " obs = " << obs_50_100 << " pred = " << pred_50_100 << " r = " << obs_50_100/pred_50_100 << endl;

        pnumber obs_100_1000 = histIntegralPN(sstack_mc_SF_jzb_pos_hall, 100, 1000);
        pnumber pred_100_1000 = histIntegralPN(pred_all_mc, 100, 1000);
        cout << "mc: [100, 1000] "<< "all" << " obs = " << obs_100_1000 << " pred = " << pred_100_1000 << " r = " << obs_100_1000/pred_100_1000 <<  endl;
        // end:printout
      }


  }

}



void plotMET_datamc()
{
  // --- plot MET all
  {
      char SF_header[] = "Same Flavor";
      char OF_header[] = "Opposite Flavor";
      string title     = "";
       
      title = metvar+plotTitle+"_SF";
      SimpleStack * sstack_data_SF    = dataDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",30,0,300, sel_cut_SF);
      TH1F *sstack_data_SF_hall       = dataDriver.getHistoTH1F(sstack_data_SF);	
      SimpleStack * sstack_mc_SF    = mcDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",30,0,300, sel_cut_SF);
      TH1F *sstack_mc_SF_hall       = mcDriver.getHistoTH1F(sstack_mc_SF);	
      SimpleLegend *sleg_SF = new SimpleLegend("TRSF");
      SimpleCanvas *simpleCan_SF = new SimpleCanvas(title.c_str(), 1);
      sstack_data_SF_hall->Draw("e1");
      sstack_mc_SF_hall->Draw("hist same");
      sstack_mc_SF->Draw("hist same");
      sstack_mc_SF_hall->Draw("hist same");
      sstack_mc_SF_hall->Draw("axis same");
      simpleCan_SF->CMSPre();
      simpleCan_SF->SetLogy();
      sleg_SF->SetHeader(SF_header);
      sleg_SF->FillLegend(sstack_data_SF);
      sleg_SF->FillLegend(sstack_mc_SF);
      sleg_SF->Draw("same");
      sstack_data_SF_hall->Draw("e1 same");
      simpleCan_SF->Save(outputDir);


      title = metvar+plotTitle+"_OF";
      SimpleStack * sstack_data_OF    = dataDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",30,0,300, sel_cut_OF);
      TH1F *sstack_data_OF_hall       = dataDriver.getHistoTH1F(sstack_data_OF);	
      SimpleStack * sstack_mc_OF    = mcDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",30,0,300, sel_cut_OF);
      TH1F *sstack_mc_OF_hall       = mcDriver.getHistoTH1F(sstack_mc_OF);	
      SimpleLegend *sleg_OF = new SimpleLegend("TROF");
      SimpleCanvas *simpleCan_OF = new SimpleCanvas(title.c_str(), 1);
      sstack_data_SF_hall->Draw("axis");
      sstack_data_OF_hall->Draw("e1 same");
      sstack_mc_OF_hall->Draw("hist same");
      sstack_mc_OF->Draw("hist same");
      sstack_mc_OF_hall->Draw("hist same");
      sstack_mc_OF_hall->Draw("axis same");
      simpleCan_OF->CMSPre();
      simpleCan_OF->SetLogy();
      sleg_OF->SetHeader(OF_header);
      sleg_OF->FillLegend(sstack_data_OF);
      sleg_OF->FillLegend(sstack_mc_OF);
      sleg_OF->Draw("same");
      sstack_data_OF_hall->Draw("e1 same");
      simpleCan_OF->Save(outputDir);
  }


}
