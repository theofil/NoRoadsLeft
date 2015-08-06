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

void simpleJZB()
{ 


  setTDRStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.5);

  v1_drivers(goFast);
  bool doJZBnorm(false);

  // needs to become a function
  string v_met = "met";
  string v_hr  = "vHT";
  string plotTitle     = "_ij2_CE_M81101_"+v_met+"_"+v_hr;
  TCut sel_cut_gauge   = TCut("!isDYTauTau") && TCut("vHT>0 && l1l2Pt>0 && jetPt[1]>50 ");
  TCut sel_cut_SF      = sel_basic && sel_M81101 && sel_ij2 && sel_CE && sel_SF_trig && sel_cut_gauge;
  TCut sel_cut_OF      = sel_basic && sel_M81101 && sel_ij2 && sel_CE && sel_OF_trig && sel_cut_gauge;
  TCut sel_cut_jzb_pos = TCut((v_hr+"-l1l2Pt>0").c_str());  
  TCut sel_cut_jzb_neg = TCut((v_hr+"-l1l2Pt<0").c_str());  
  string metvar        = v_met.c_str();
  string jzbvar        = (v_hr+"-l1l2Pt").c_str();
  pnumber OF_scale     = pnumber(1, 0.05);


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

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // --- plot JZB all for mc
  {
      char SF_header[] = "Same Flavor";
      char OF_header[] = "Opposite Flavor";
      string title     = "";
       
      title = "jzb"+plotTitle+"_SF";
      SimpleStack * sstack_data_SF    = dataDriver.getSimpleStackTH1F(jzbvar,";JZB [GeV]; events / 10 GeV;",40,-200,200, sel_cut_SF);
      TH1F *sstack_data_SF_hall       = dataDriver.getHistoTH1F(sstack_data_SF);	
      SimpleStack * sstack_mc_SF    = mcDriver.getSimpleStackTH1F(jzbvar,";JZB [GeV]; events / 10 GeV;",40,-200,200, sel_cut_SF);
      TH1F *sstack_mc_SF_hall       = mcDriver.getHistoTH1F(sstack_mc_SF);	
      SimpleLegend *sleg_SF = new SimpleLegend("TLSF");
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

      title = "jzb"+plotTitle+"_OF";
      SimpleStack * sstack_data_OF    = dataDriver.getSimpleStackTH1F(jzbvar,";JZB [GeV]; events / 10 GeV;",40,-200,200, sel_cut_OF);
      TH1F *sstack_data_OF_hall       = dataDriver.getHistoTH1F(sstack_data_OF);	
      SimpleStack * sstack_mc_OF    = mcDriver.getSimpleStackTH1F(jzbvar,";JZB [GeV]; events / 10 GeV;",40,-200,200, sel_cut_OF);
      TH1F *sstack_mc_OF_hall       = mcDriver.getHistoTH1F(sstack_mc_OF);	
      SimpleLegend *sleg_OF = new SimpleLegend("TLOF");
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

  // --- plot MET for jzb_pos and jzb_neg
  {
      char SF_jzb_pos_header[]    = "#splitline{Same Flavor}{JZB > 0}";
      char OF_jzb_pos_header[]    = "#splitline{Opposite Flavor}{JZB > 0}";
      char SF_jzb_neg_header[]    = "#splitline{Same Flavor}{JZB < 0}";
      char OF_jzb_neg_header[]    = "#splitline{Opposite Flavor}{JZB < 0}";
      string title        = "";
      
      title = metvar+plotTitle+"_SF"+"_jzb_pos";
      SimpleStack * sstack_data_SF_jzb_pos    = dataDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",30,0,300, sel_cut_SF && sel_cut_jzb_pos);
      TH1F *sstack_data_SF_jzb_pos_hall       = dataDriver.getHistoTH1F(sstack_data_SF_jzb_pos);	
      SimpleStack * sstack_mc_SF_jzb_pos    = mcDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",30,0,300, sel_cut_SF && sel_cut_jzb_pos);
      TH1F *sstack_mc_SF_jzb_pos_hall       = mcDriver.getHistoTH1F(sstack_mc_SF_jzb_pos);	
      SimpleLegend *sleg_SF_jzb_pos = new SimpleLegend("TRSFjzb_pos");
      SimpleCanvas *simpleCan_SF_jzb_pos = new SimpleCanvas(title.c_str(), 1);
      sstack_data_SF_jzb_pos_hall->Draw("e1");
      sstack_mc_SF_jzb_pos_hall->Draw("hist same");
      sstack_mc_SF_jzb_pos->Draw("hist same");
      sstack_mc_SF_jzb_pos_hall->Draw("hist same");
      sstack_mc_SF_jzb_pos_hall->Draw("axis same");
      simpleCan_SF_jzb_pos->CMSPre();
      simpleCan_SF_jzb_pos->SetLogy();
      sleg_SF_jzb_pos->SetHeader(SF_jzb_pos_header);
      sleg_SF_jzb_pos->FillLegend(sstack_data_SF_jzb_pos);
      sleg_SF_jzb_pos->FillLegend(sstack_mc_SF_jzb_pos);
      sleg_SF_jzb_pos->Draw("same");
      sstack_data_SF_jzb_pos_hall->Draw("e1 same");
      simpleCan_SF_jzb_pos->Save(outputDir);

      title = metvar+plotTitle+"_OF"+"_jzb_pos";
      SimpleStack * sstack_data_OF_jzb_pos    = dataDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",30,0,300, sel_cut_OF && sel_cut_jzb_pos);
      TH1F *sstack_data_OF_jzb_pos_hall       = dataDriver.getHistoTH1F(sstack_data_OF_jzb_pos);	
      SimpleStack * sstack_mc_OF_jzb_pos    = mcDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",30,0,300, sel_cut_OF && sel_cut_jzb_pos);
      TH1F *sstack_mc_OF_jzb_pos_hall       = mcDriver.getHistoTH1F(sstack_mc_OF_jzb_pos);	
      SimpleLegend *sleg_OF_jzb_pos = new SimpleLegend("TROFjzb_pos");
      SimpleCanvas *simpleCan_OF_jzb_pos = new SimpleCanvas(title.c_str(), 1);
      sstack_data_SF_jzb_pos_hall->Draw("axis");
      sstack_data_OF_jzb_pos_hall->Draw("e1 same");
      sstack_mc_OF_jzb_pos_hall->Draw("hist same");
      sstack_mc_OF_jzb_pos->Draw("hist same");
      sstack_mc_OF_jzb_pos_hall->Draw("hist same");
      sstack_mc_OF_jzb_pos_hall->Draw("axis same");
      simpleCan_OF_jzb_pos->CMSPre();
      simpleCan_OF_jzb_pos->SetLogy();
      sleg_OF_jzb_pos->SetHeader(OF_jzb_pos_header);
      sleg_OF_jzb_pos->FillLegend(sstack_data_OF_jzb_pos);
      sleg_OF_jzb_pos->FillLegend(sstack_mc_OF_jzb_pos);
      sleg_OF_jzb_pos->Draw("same");
      sstack_data_OF_jzb_pos_hall->Draw("e1 same");
      simpleCan_OF_jzb_pos->Save(outputDir);


      title = metvar+plotTitle+"_SF"+"_jzb_neg";
      SimpleStack * sstack_data_SF_jzb_neg    = dataDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",30,0,300, sel_cut_SF && sel_cut_jzb_neg);
      TH1F *sstack_data_SF_jzb_neg_hall       = dataDriver.getHistoTH1F(sstack_data_SF_jzb_neg);	
      SimpleStack * sstack_mc_SF_jzb_neg    = mcDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",30,0,300, sel_cut_SF && sel_cut_jzb_neg);
      TH1F *sstack_mc_SF_jzb_neg_hall       = mcDriver.getHistoTH1F(sstack_mc_SF_jzb_neg);	
      SimpleLegend *sleg_SF_jzb_neg = new SimpleLegend("TRSFjzb_neg");
      SimpleCanvas *simpleCan_SF_jzb_neg = new SimpleCanvas(title.c_str(), 1);
      sstack_data_SF_jzb_neg_hall->Draw("e1");
      sstack_mc_SF_jzb_neg_hall->Draw("hist same");
      sstack_mc_SF_jzb_neg->Draw("hist same");
      sstack_mc_SF_jzb_neg_hall->Draw("hist same");
      sstack_mc_SF_jzb_neg_hall->Draw("axis same");
      simpleCan_SF_jzb_neg->CMSPre();
      simpleCan_SF_jzb_neg->SetLogy();
      sleg_SF_jzb_neg->SetHeader(SF_jzb_neg_header);
      sleg_SF_jzb_neg->FillLegend(sstack_data_SF_jzb_neg);
      sleg_SF_jzb_neg->FillLegend(sstack_mc_SF_jzb_neg);
      sleg_SF_jzb_neg->Draw("same");
      sstack_data_SF_jzb_neg_hall->Draw("e1 same");
      simpleCan_SF_jzb_neg->Save(outputDir);

      title = metvar+plotTitle+"_OF"+"_jzb_neg";
      SimpleStack * sstack_data_OF_jzb_neg    = dataDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",30,0,300, sel_cut_OF && sel_cut_jzb_neg);
      TH1F *sstack_data_OF_jzb_neg_hall       = dataDriver.getHistoTH1F(sstack_data_OF_jzb_neg);	
      SimpleStack * sstack_mc_OF_jzb_neg    = mcDriver.getSimpleStackTH1F(metvar,";MET [GeV]; events / 10 GeV;",30,0,300, sel_cut_OF && sel_cut_jzb_neg);
      TH1F *sstack_mc_OF_jzb_neg_hall       = mcDriver.getHistoTH1F(sstack_mc_OF_jzb_neg);	
      SimpleLegend *sleg_OF_jzb_neg = new SimpleLegend("TROFjzb_neg");
      SimpleCanvas *simpleCan_OF_jzb_neg = new SimpleCanvas(title.c_str(), 1);
      sstack_data_SF_jzb_neg_hall->Draw("axis");
      sstack_data_OF_jzb_neg_hall->Draw("e1");
      sstack_mc_OF_jzb_neg_hall->Draw("hist same");
      sstack_mc_OF_jzb_neg->Draw("hist same");
      sstack_mc_OF_jzb_neg_hall->Draw("hist same");
      sstack_mc_OF_jzb_neg_hall->Draw("axis same");
      simpleCan_OF_jzb_neg->CMSPre();
      simpleCan_OF_jzb_neg->SetLogy();
      sleg_OF_jzb_neg->SetHeader(OF_jzb_neg_header);
      sleg_OF_jzb_neg->FillLegend(sstack_data_OF_jzb_neg);
      sleg_OF_jzb_neg->FillLegend(sstack_mc_OF_jzb_neg);
      sleg_OF_jzb_neg->Draw("same");
      sstack_data_OF_jzb_neg_hall->Draw("e1 same");
      simpleCan_OF_jzb_neg->Save(outputDir);

      //-------------------------------------------
      //--- jzb core for mc
      //-------------------------------------------
      
      pnumber jzb_pos_0_50_mc_SF = histIntegralPN(sstack_mc_SF_jzb_pos_hall, 0 ,50);
      pnumber jzb_neg_0_50_mc_SF = histIntegralPN(sstack_mc_SF_jzb_neg_hall, 0 ,50);
      pnumber jzb_pos_0_50_mc_OF = histIntegralPN(sstack_mc_OF_jzb_pos_hall, 0 ,50);
      pnumber jzb_neg_0_50_mc_OF = histIntegralPN(sstack_mc_OF_jzb_neg_hall, 0 ,50);

      pnumber jzb_norm_mc = (jzb_pos_0_50_mc_SF-jzb_pos_0_50_mc_OF)/(jzb_neg_0_50_mc_SF-jzb_neg_0_50_mc_OF);
      cout << "jzb_pos_0_50_mc_SF/jzb_neg_0_50_mc_SF = " << jzb_pos_0_50_mc_SF/jzb_neg_0_50_mc_SF << endl;
      if(!doJZBnorm) jzb_norm_mc = pnumber(1,0);
      cout << "jzb_norm_mc = " << jzb_norm_mc << endl;

      vector<TH1F*> SF_jzb_pos_TH1F_vec_mc = sstack_mc_SF_jzb_pos->histoPointers_;
      vector<TH1F*> SF_jzb_neg_TH1F_vec_mc = sstack_mc_SF_jzb_neg->histoPointers_;
      vector<TH1F*> OF_jzb_pos_TH1F_vec_mc = sstack_mc_OF_jzb_pos->histoPointers_;
      vector<TH1F*> OF_jzb_neg_TH1F_vec_mc = sstack_mc_OF_jzb_neg->histoPointers_;
    
      vector<string> sampleTitles_mc = sstack_mc_SF_jzb_pos->sampleTitles_;
      vector<TH1F*> jzb_predictions_mc;

      for(size_t ii = 0; ii < sampleTitles_mc.size(); ++ii)
      {
       	cout <<"processing "<< sampleTitles_mc[ii] << endl;

       	// --- assumes the OF_scale is uncorrelated in JZB > 0, JZB < 0
        TH1F *SF_jzb_pos = (TH1F*)SF_jzb_pos_TH1F_vec_mc[ii]->Clone(); 
        TH1F *SF_jzb_neg = (TH1F*)SF_jzb_neg_TH1F_vec_mc[ii]->Clone();
        TH1F *OF_jzb_neg = ScaleTH1F(OF_jzb_neg_TH1F_vec_mc[ii], OF_scale);
        TH1F *OF_jzb_pos = ScaleTH1F(OF_jzb_pos_TH1F_vec_mc[ii], OF_scale);
        TH1F *jzb_only = (TH1F*)SF_jzb_neg->Clone();
        jzb_only ->Add(OF_jzb_neg,-1);
        TH1F *jzb_scaled = ScaleTH1F(jzb_only, jzb_norm_mc);
        TH1F *pred = jzb_scaled;
        pred->Add(OF_jzb_pos);

        jzb_predictions_mc.push_back(pred);

        SF_jzb_pos->SetLineColor(kBlack);
        SF_jzb_pos->SetFillColor(kWhite);
        pred->SetLineColor(kRed);
        pred->SetFillColor(kWhite);
        pred->SetLineWidth(2);
        pred->SetLineStyle(2);

        title = metvar+plotTitle+"_mc"+"_pred"+any2string(ii);  
        SimpleLegend *sleg      = new SimpleLegend("jzb_exp_mode");
        SimpleCanvas *simpleCan = new SimpleCanvas(title.c_str(), 2);
        simpleCan->CMSPre();
        simpleCan->SetLogy();
        simpleCan->Up();
        simpleCan->ShapeMeUp(SF_jzb_pos); 
        simpleCan->ShapeMeUp(pred); 
        SF_jzb_pos->Draw("hist");
        pred->Draw("hist same"); 
        sleg->SetHeader((sampleTitles_mc[ii]+" [SF]").c_str());
        sleg->AddEntry(SF_jzb_pos,"JZB > 0","FL");
        sleg->AddEntry(pred, "prediction", "FL");
        sleg->Draw("same");
        simpleCan->Dw();
        TH1F *hratio = doSimpleRatio(SF_jzb_pos, pred);
        hratio->GetYaxis()->SetTitle("ratio");
        hratio->GetYaxis()->CenterTitle();
        simpleCan->ShapeMeDw(hratio);
        hratio->Draw("e1");
        hratio->GetYaxis()->SetRangeUser(0.0,2.0);
        hratio->GetYaxis()->SetNdivisions(507);
        simpleCan->Save(outputDir);
      }	

      TH1F *pred_all_mc;
      for(size_t ii = 0; ii < sampleTitles_mc.size(); ++ii)
      {
	 if(ii==0) pred_all_mc = (TH1F*)jzb_predictions_mc[ii]->Clone("pred_all_mc");
         if(ii!=0) pred_all_mc -> Add(jzb_predictions_mc[ii]); 

         // --- printout
         pnumber obs_50_100 = histIntegralPN(SF_jzb_pos_TH1F_vec_mc[ii], 50, 100);
         pnumber pred_50_100 = histIntegralPN(jzb_predictions_mc[ii], 50, 100);
         cout << "[50, 100] "<< sampleTitles_mc[ii] << " obs = " << obs_50_100 << " pred = " << pred_50_100 << " r = " << obs_50_100/pred_50_100 << endl;

         pnumber obs_100_1000 = histIntegralPN(SF_jzb_pos_TH1F_vec_mc[ii], 100, 1000);
         pnumber pred_100_1000 = histIntegralPN(jzb_predictions_mc[ii], 100, 1000);
         cout << "[100, 1000] "<< sampleTitles_mc[ii] << " obs = " << obs_100_1000 << " pred = " << pred_100_1000 << " r = " << obs_50_100/pred_50_100 <<  endl;
         // end:printout
      }

      {
        TH1F *sstack_mc_SF_jzb_pos_hall_clone = (TH1F*)sstack_mc_SF_jzb_pos_hall->Clone();
        sstack_mc_SF_jzb_pos_hall_clone->SetLineColor(kBlack);
  	title = metvar+plotTitle+"_pred_all_mc";  
  	SimpleLegend *sleg      = new SimpleLegend("jzb_exp_mode");
       	SimpleCanvas *simpleCan = new SimpleCanvas(title.c_str(), 2);
        simpleCan->CMSPre();
        simpleCan->SetLogy();
        simpleCan->Up();
        simpleCan->ShapeMeUp(sstack_mc_SF_jzb_pos_hall_clone); 
        simpleCan->ShapeMeUp(pred_all_mc); 
        sstack_mc_SF_jzb_pos_hall_clone->Draw("hist");
        pred_all_mc->Draw("hist same"); 
        sleg->SetHeader("Same Flavor [MC]");
        sleg->AddEntry(sstack_mc_SF_jzb_pos_hall_clone,"JZB > 0","FL");
        sleg->AddEntry(pred_all_mc, "prediction", "FL");
        sleg->Draw("same");
        simpleCan->Dw();
        TH1F *hratio = doSimpleRatio(sstack_mc_SF_jzb_pos_hall_clone, pred_all_mc);
        hratio->GetYaxis()->SetTitle("ratio");
        hratio->GetYaxis()->CenterTitle();
        simpleCan->ShapeMeDw(hratio);
        hratio->Draw("e1");
        hratio->GetYaxis()->SetRangeUser(0.0,2.0);
        hratio->GetYaxis()->SetNdivisions(507);
        simpleCan->Save(outputDir);

        // --- printout
        pnumber obs_50_100 = histIntegralPN(sstack_mc_SF_jzb_pos_hall_clone, 50, 100);
        pnumber pred_50_100 = histIntegralPN(pred_all_mc, 50, 100);
        cout << "mc: [50, 100] "<< "all" << " obs = " << obs_50_100 << " pred = " << pred_50_100 << " r = " << obs_50_100/pred_50_100 << endl;

        pnumber obs_100_1000 = histIntegralPN(sstack_mc_SF_jzb_pos_hall_clone, 100, 1000);
        pnumber pred_100_1000 = histIntegralPN(pred_all_mc, 100, 1000);
        cout << "mc: [100, 1000] "<< "all" << " obs = " << obs_100_1000 << " pred = " << pred_100_1000 << " r = " << obs_100_1000/pred_100_1000 <<  endl;
        // end:printout
      }

      //-------------------------------------------
      //--- jzb core for data
      //-------------------------------------------

      pnumber jzb_pos_0_50_data_SF = histIntegralPN(sstack_data_SF_jzb_pos_hall, 0 ,50);
      pnumber jzb_neg_0_50_data_SF = histIntegralPN(sstack_data_SF_jzb_neg_hall, 0 ,50);
      pnumber jzb_pos_0_50_data_OF = histIntegralPN(sstack_data_OF_jzb_pos_hall, 0 ,50);
      pnumber jzb_neg_0_50_data_OF = histIntegralPN(sstack_data_OF_jzb_neg_hall, 0 ,50);

      pnumber jzb_norm_data = (jzb_pos_0_50_data_SF-jzb_pos_0_50_data_OF)/(jzb_neg_0_50_data_SF-jzb_neg_0_50_data_OF);
      cout << "jzb_pos_0_50_data_SF/jzb_neg_0_50_data_SF = " << jzb_pos_0_50_data_SF/jzb_neg_0_50_data_SF << endl;
      if(!doJZBnorm)jzb_norm_data = pnumber(1,0);
      cout << "jzb_norm_data = " << jzb_norm_data << endl;

      vector<TH1F*> SF_jzb_pos_TH1F_vec_data = sstack_data_SF_jzb_pos->histoPointers_;
      vector<TH1F*> SF_jzb_neg_TH1F_vec_data = sstack_data_SF_jzb_neg->histoPointers_;
      vector<TH1F*> OF_jzb_pos_TH1F_vec_data = sstack_data_OF_jzb_pos->histoPointers_;
      vector<TH1F*> OF_jzb_neg_TH1F_vec_data = sstack_data_OF_jzb_neg->histoPointers_;
    
      vector<string> sampleTitles_data = sstack_data_SF_jzb_pos->sampleTitles_;
      vector<TH1F*> jzb_predictions_data;

      for(size_t ii = 0; ii < sampleTitles_data.size(); ++ii)
      {
       	cout <<"processing "<< sampleTitles_data[ii] << endl;

       	// --- assumes the OF_scale is uncorrelated in JZB > 0, JZB < 0
        TH1F *SF_jzb_pos = (TH1F*)SF_jzb_pos_TH1F_vec_data[ii]->Clone(); 
        TH1F *SF_jzb_neg = (TH1F*)SF_jzb_neg_TH1F_vec_data[ii]->Clone();
        TH1F *OF_jzb_neg = ScaleTH1F(OF_jzb_neg_TH1F_vec_data[ii], OF_scale);
        TH1F *OF_jzb_pos = ScaleTH1F(OF_jzb_pos_TH1F_vec_data[ii], OF_scale);
        TH1F *jzb_only = (TH1F*)SF_jzb_neg->Clone();
        jzb_only ->Add(OF_jzb_neg,-1);
        TH1F *jzb_scaled = ScaleTH1F(jzb_only, jzb_norm_data);
        TH1F *pred = jzb_scaled;
        pred->Add(OF_jzb_pos);

        jzb_predictions_data.push_back(pred);

        SF_jzb_pos->SetLineColor(kBlack);
        SF_jzb_pos->SetFillColor(kWhite);
        pred->SetLineColor(kRed);
        pred->SetFillColor(kWhite);
        pred->SetLineWidth(2);
        pred->SetLineStyle(2);

//        title = metvar+plotTitle+"_data"+"_pred"+any2string(ii);  
//        SimpleLegend *sleg      = new SimpleLegend("jzb_exp_mode");
//        SimpleCanvas *simpleCan = new SimpleCanvas(title.c_str(), 2);
//        simpleCan->CMSPre();
//        simpleCan->SetLogy();
//        simpleCan->Up();
//        simpleCan->ShapeMeUp(SF_jzb_pos); 
//        simpleCan->ShapeMeUp(pred); 
//        SF_jzb_pos->Draw("hist");
//        pred->Draw("hist same"); 
//        sleg->SetHeader((sampleTitles_data[ii]+" [SF]").c_str());
//        sleg->AddEntry(SF_jzb_pos,"JZB > 0","FL");
//        sleg->AddEntry(pred, "prediction", "FL");
//        sleg->Draw("same");
//        simpleCan->Dw();
//        TH1F *hratio = doRatio(SF_jzb_pos, pred);
//        hratio->GetYaxis()->SetTitle("ratio");
//        hratio->GetYaxis()->CenterTitle();
//        simpleCan->ShapeMeDw(hratio);
//        hratio->Draw("e1");
//        hratio->GetYaxis()->SetRangeUser(0.0,2.0);
//        hratio->GetYaxis()->SetNdivisions(507);
//        simpleCan->Save(outputDir);
      }	

      TH1F *pred_all_data;
      for(size_t ii = 0; ii < sampleTitles_data.size(); ++ii)
      {
	 if(ii==0) pred_all_data = (TH1F*)jzb_predictions_data[ii]->Clone("pred_all_data");
         if(ii!=0) pred_all_data -> Add(jzb_predictions_data[ii]); 

         // --- printout
         pnumber obs_50_100 = histIntegralPN(SF_jzb_pos_TH1F_vec_data[ii], 50, 100);
         pnumber pred_50_100 = histIntegralPN(jzb_predictions_data[ii], 50, 100);
         cout << "[50, 100] "<< sampleTitles_data[ii] << " obs = " << obs_50_100 << " pred = " << pred_50_100 << " r = " << obs_50_100/pred_50_100 << endl;

         pnumber obs_100_1000 = histIntegralPN(SF_jzb_pos_TH1F_vec_data[ii], 100, 1000);
         pnumber pred_100_1000 = histIntegralPN(jzb_predictions_data[ii], 100, 1000);
         cout << "[100, 1000] "<< sampleTitles_data[ii] << " obs = " << obs_100_1000 << " pred = " << pred_100_1000 << " r = " << obs_50_100/pred_50_100 <<  endl;
         // end:printout
      }

      {
  	title = metvar+plotTitle+"_pred_all_data";  
  	SimpleLegend *sleg      = new SimpleLegend("jzb_exp_mode");
       	SimpleCanvas *simpleCan = new SimpleCanvas(title.c_str(), 2);
        simpleCan->CMSPre();
        simpleCan->SetLogy();
        simpleCan->Up();
        simpleCan->ShapeMeUp(sstack_data_SF_jzb_pos_hall); 
        simpleCan->ShapeMeUp(pred_all_data); 
        sstack_data_SF_jzb_pos_hall->Draw("e1");
        pred_all_data->Draw("hist same"); 
        sleg->SetHeader("Same Flavor [Data]");
        sleg->AddEntry(sstack_data_SF_jzb_pos_hall,"JZB > 0","FL");
        sleg->AddEntry(pred_all_data, "prediction", "FL");
        sleg->Draw("same");
        sstack_data_SF_jzb_pos_hall->Draw("e1 same");
        simpleCan->Dw();
        TH1F *hratio = doRatio(sstack_data_SF_jzb_pos_hall, pred_all_data);
        hratio->GetYaxis()->SetTitle("ratio");
        hratio->GetYaxis()->CenterTitle();
        simpleCan->ShapeMeDw(hratio);
        hratio->Draw("e1");
        hratio->GetYaxis()->SetRangeUser(0.0,2.0);
        hratio->GetYaxis()->SetNdivisions(507);
        simpleCan->Save(outputDir);

        // --- printout
        pnumber obs_50_100 = histIntegralPN(sstack_data_SF_jzb_pos_hall, 50, 100);
        pnumber pred_50_100 = histIntegralPN(pred_all_data, 50, 100);
        cout << "data: [50, 100] "<< "all" << " obs = " << obs_50_100 << " pred = " << pred_50_100 << " r = " << obs_50_100/pred_50_100 << endl;

        pnumber obs_100_1000 = histIntegralPN(sstack_data_SF_jzb_pos_hall, 100, 1000);
        pnumber pred_100_1000 = histIntegralPN(pred_all_data, 100, 1000);
        cout << "data: [100, 1000] "<< "all" << " obs = " << obs_100_1000 << " pred = " << pred_100_1000 << " r = " << obs_100_1000/pred_100_1000 <<  endl;
        // end:printout
      }

  }

}



