#include "../core/utilities.h"
#include "../core/LocalFilePath.h"
#include "../core/simpleROOT_cuts.h"
#include "../friends/readNTPL.h"

using namespace std;
bool goFast(false);


TFile *fp_out;
TTree *events_out;

void makeAuxFile()
{ 
  string filename_in = fp_DYJetsM50;
  string filename_out = outputDir + "dym50.aux.root";
  bool hasNegWeights = true;

//  string filename_in = fp_TTJets;
//  string filename_out = outputDir + "ttbar.aux.root";
//  bool hasNegWeights = true;

  cout << "opening " << filename_in << endl;
  cout << "will write output to " << filename_out << endl;

  fp_in      = new TFile ( filename_in.c_str()  ,  "OPEN");
  events_in  = (TTree*)fp_in->Get("demo/events");
  init();
  
  events_in->SetBranchStatus("*",0); //disable all branches
  events_in->SetBranchStatus("genWeight", 1);

  //setting events_out, fp_out
  fp_out = new TFile(filename_out.c_str(), "RECREATE");
  events_out = new TTree("events","events"); 

  // setting variables to be stored
  float  totEveW;
  events_out->Branch("totEveW",&totEveW, "totEveW/F"); 

  // calculating Neff here
  float Neff = 0;
  
  float tot_entries = events_in ->GetEntries();

  mcWeights_ = (TH1I*)fp_in->Get("demo/mcWeights");
  float neg_eve = mcWeights_->GetBinContent(2);
  float pos_eve = mcWeights_->GetBinContent(3);
  float frac_neg_eve = neg_eve/tot_entries;
  float tot_eff_entries = tot_entries*(1-2*frac_neg_eve);

  cout << "neg_eve + pos_eve = " << neg_eve << " + " << pos_eve << " = " << neg_eve + pos_eve << " ; tot_entries = " << tot_entries << endl;
  cout << "frac_neg_eve = neg_eve/tot_entries = " << frac_neg_eve << endl;  
  cout << "tot_eff_entries = tot_entries*(1-2*frac_neg_eve)  = " << tot_eff_entries << endl;
  Neff = tot_eff_entries;

 
  float NeffInv = 1/Neff;
  cout << "NeffInv = " << NeffInv << endl;

  // --- event loop
  ULong64_t nentries = events_in->GetEntries();
  if(goFast)nentries = nentries*0.01;

  // --- calculate sum of weights
  double sumW =0 ;
  for (ULong64_t ii=0; ii<nentries;ii++)
  {
    events_in->GetEntry(ii);
    sumW += genWeight_;
  }
  cout << "sumW = " << sumW << endl;

  double init_test = 0;
  // --- standard loop
  for (ULong64_t ii=0; ii<nentries;ii++)
  {
     events_in->GetEntry(ii);

     float  genWeightSign = 0 ;
     if(genWeight_ > 0) genWeightSign = +1;  
     if(genWeight_ < 0) genWeightSign = -1;
     if(genWeight_ == 0) cout << "zero weight was found for entry: " << ii << endl;

     totEveW = NeffInv*genWeightSign ;
     init_test += totEveW;

     if(ii < 30) cout << setprecision(10)<< "genWeight_ =" << genWeight_ << " ; totEveW = " << totEveW << "; genWeight_/sumW = " << genWeight_/float(sumW) << endl;
     if(fabs(totEveW - genWeight_/float(sumW))>0.001) cout << "### genWeight_ =" << genWeight_ << " ; totEveW = " << totEveW << "; genWeight_/sumW = " << genWeight_/float(sumW) << endl;


     events_out->Fill();
  }


  cout << "init_test = " << init_test << endl;
  cout << "writting the output file: " << filename_out << endl;

  fp_out->cd();
  events_out->Write();
  fp_out->Close();
}


/*
     sel_basic = false;
     if(nleps_>=2 && lepPt_[0]>25 && lepPt_[1]>20 && l1l2M_>81 && l1l2M_<101 && (lepID_[0]*lepID_[1]==-11*11 || lepID_[0]*lepID_[1]==-13*13))
     if(!isDYTauTau_)
     sel_basic = true; 

     TLorentzVector l1(0,0,0,0);
     TLorentzVector l2(0,0,0,0);
     if(nleps_>=2)
     {
        l1.SetPtEtaPhiM(lepPt_[0], lepEta_[0], lepPhi_[0], lepM_[0]);
        l2.SetPtEtaPhiM(lepPt_[1], lepEta_[1], lepPhi_[1], lepM_[1]);
     } 
*/

