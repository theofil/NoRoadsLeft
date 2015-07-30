#include "../core/utilities.h"
#include "../core/LocalFilePath.h"
#include "../core/simpleROOT_cuts.h"
#include "../core/ActiveSamples.h"
#include "TSystem.h"


TFile *fp_in;
TTree *events_in;

bool goFast(false);

void makeSkims()
{

  ActiveSamples();

  for(auto it: myActiveSamples) 
  {
   float xsection      = std::stof(string(it.sampleXS_.GetTitle())) ;
   string filename     = it.sampleName_;
   string filename_in  = localPath + filename;
   string filename_out = outputDir + filename;
 
   cout << "reading: " << filename_in << " with xsection = "<< xsection << endl;
   cout << "writing: " << filename_out << endl;
 
   //setting fp_in, events    
   fp_in      = new TFile ( filename_in.c_str()  ,  "OPEN");
   events_in     = (TTree*)fp_in->Get("demo/events");

   // loading friend TTree
   TTree *events_friend;
   int filenameSize = filename_in.size();
   string friendTreePossiblePath = filename_in.substr(0,filenameSize-5) + ".aux.root";
   FileStat_t myfilestat;
   bool auxFileExists = !gSystem->GetPathInfo(friendTreePossiblePath.c_str(), myfilestat) ;
   if(!auxFileExists){cout << "### ERROR: can't find friend tree " << friendTreePossiblePath.c_str() << endl;}
   TFile *fp_friend = new TFile(friendTreePossiblePath.c_str(), "OPEN");
   events_friend = (TTree*)fp_friend->Get("events");
   float totEveW;
   events_friend->SetBranchAddress("totEveW", &totEveW);

   // set branches that will take part in the event selection
   unsigned short nleps = 0;
   events_in->SetBranchAddress("nleps", &nleps);
   unsigned short njets = 0;
   events_in->SetBranchAddress("njets", &njets);
 
   //setting fp_out
   TFile * fp_out = new TFile(filename_out.c_str(), "RECREATE");
   fp_out->cd();
   fp_out->mkdir("demo");
   fp_out->cd("demo");

   TTree *events_out = events_in->CloneTree(0);
   events_out->Branch("xsection",&xsection, "xsection/F");
   events_out->Branch("totEveW",&totEveW, "totEveW/F");
 
   Long64_t nentries = events_in->GetEntries(); 
   if(goFast) nentries = 100;   

   for (Long64_t ii = 0 ; ii < nentries; ii++) 
   {
      events_in->GetEntry(ii);
      events_friend->GetEntry(ii); // you definetly need this!
      if(nleps==2 && njets==1)events_out->Fill();
   }
   
   events_out->AutoSave();
   cout << "closing: " << filename_out << endl;
   fp_out->Close();
  }
}

