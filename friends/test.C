#include "../core/utilities.h"
#include "../core/LocalFilePath.h"
#include "../core/simpleROOT_cuts.h"


using namespace std;
bool goFast(false);

TFile *fp_in;
TTree *events_in;

TFile *fp_out;
TTree *events_out;

TFile *fp_calib;
TProfile *pro_ij0;

void init();

// --- tree variables
#define njetsMax 30 
#define nrjetsMax 10 
#define nlepsMax 10
#define ngenlepsMax 10
#define ngenpartsMax 10

bool goodVtx_; 
unsigned short nVtx_; 

unsigned long long eventNum_;
//unsigned long eventNum_;
unsigned int runNum_;
unsigned int lumi_;

float l1l2DPhi_;
float l1l2DR_;
float l1l2Pt_;
float l1l2M_;
float l1l2Eta_;
float l1l2Phi_;

unsigned short nleps_;
float lepPt_                      [nlepsMax]; 
float lepEta_                     [nlepsMax]; 
float lepPhi_                     [nlepsMax]; 
float lepM_                       [nlepsMax];   
float lepIso_                     [nlepsMax];
float lepPtRel_                   [nlepsMax];
short lepID_                      [nlepsMax];
short lepGenMatchIndex_           [nlepsMax];   // index in the stored genLep[] collection
bool  lepTriggerMatch_            [nlepsMax];

unsigned short njets_;
float jetPt_                      [njetsMax]; 
float jetEta_                     [njetsMax]; 
float jetPhi_                     [njetsMax]; 
float jetM_                       [njetsMax]; 
float jetBTag_                    [njetsMax];

unsigned short nrjets_;          
float rjetPt_                     [nrjetsMax]; 
float rjetEta_                    [nrjetsMax]; 
float rjetPhi_                    [nrjetsMax]; 
float rjetM_                      [nrjetsMax]; 
float rjetBTag_                   [nrjetsMax];

unsigned short ngenleps_;
float genlepPt_                   [ngenlepsMax]; 
float genlepEta_                  [ngenlepsMax]; 
float genlepPhi_                  [ngenlepsMax]; 
float genlepM_                    [ngenlepsMax];   
short genlepID_                   [ngenlepsMax];
int   genlepMID_                  [ngenlepsMax]; // mother 
int   genlepGMID_                 [ngenlepsMax]; // grand mother
int   genlepGGMID_                [ngenlepsMax]; // grand grand mother

float genl1l2DPhi_;
float genl1l2DR_;
float genl1l2Pt_;
float genl1l2M_;
float genl1l2Eta_;
float genl1l2Phi_;

unsigned short ngenparts_;
float genpartPt_                   [ngenpartsMax]; 
float genpartEta_                  [ngenpartsMax]; 
float genpartPhi_                  [ngenpartsMax]; 
float genpartM_                    [ngenpartsMax];   
int   genpartID_                   [ngenpartsMax]; 
int   genpartDID1_                 [ngenpartsMax]; // daugther 1
int   genpartDID2_                 [ngenpartsMax]; // daugther 2
short genpartDRMI1_                [ngenpartsMax]; // daugther reco match index (genleptons to leptons, genjet to jet)
short genpartDRMI2_                [ngenpartsMax]; // daugther reco match index 

float met_;          // raw pf-met
float metPhi_;
float t1met_;
float t1metPhi_;
float genmet_;          
float sumEt_;
float t1sumEt_;

float rho_;
float nPU_;
float nPUTrue_;

float vHT_;           // pt of the recoil vector of all objects excluding the dilepton [= -met - l1l2] 
float t1vHT_;         // as vHT but with t1 correction
float jvHT_;          // recoil of hard jets and subleading leptons

bool HLT_e1e2_; 
bool HLT_mu1mu2_; 
bool HLT_mu1e2_; 
bool HLT_e1mu2_; 
bool HLT_pfmet_; 
bool HLT_pfmetCSV_; 

bool Flag_trackingFailureFilter_;		        
bool Flag_goodVertices_;			 
bool Flag_CSCTightHaloFilter_;		 
bool Flag_trkPOGFilters_;			 
bool Flag_trkPOG_logErrorTooManyClusters_;	 
bool Flag_EcalDeadCellTriggerPrimitiveFilter_; 
bool Flag_ecalLaserCorrFilter_;		 
bool Flag_trkPOG_manystripclus53X_;		 
bool Flag_eeBadScFilter_;			 
bool Flag_METFilters_;			 
bool Flag_HBHENoiseFilter_;			 
bool Flag_trkPOG_toomanystripclus53X_;	 
bool Flag_hcalLaserEventFilter_;		 

bool isDYTauTau_; 

unsigned short nphos_;
// --- end of ntuple variables

void makejzb()
{ 
  string filename_in = fp_DYJetsToLL;
  string filename_out = outputDir + "dy.aux.root";

//  string filename_in = fp_TTJets;
//  string filename_out = outputDir + "ttbar.aux.root";

  string calibration_file = "data/jzb_resp_mc_Phys14_2015_05_21.root";

  cout << "opening " << filename_in << endl;
  cout << "will write output to " << filename_out << endl;

  fp_in      = new TFile ( filename_in.c_str()  ,  "OPEN");
  events_in  = (TTree*)fp_in->Get("demo/events");
  init();

  // opening calibration file
  cout << "opening calibration file: " << calibration_file << endl;
  fp_calib = new TFile(calibration_file.c_str(), "READ");
  pro_ij0 = (TProfile*)fp_calib->Get("pro_ij0");

  //setting events_out, fp_out
  fp_out = new TFile(filename_out.c_str(), "RECREATE");
  events_out = new TTree("events","events"); 
  float  cjzb_resp;
  unsigned short    njets;
  float  l1l2Pt;
  float  l1l2M;
  float  t1vHT;
  float  resp;
  float  jzbraw;
  float  t1met;
  float  l1jl2jPt;
  float  l1jl2jM;
  float  recoilmet_lep;
  float  recoiljzb_lep;
  float  MHT;
  float  rvHT;
  bool   sel_basic;
  
  events_out->Branch("cjzb_resp",&cjzb_resp, "cjzb_resp/F"); // jzb response corrected from TProfile
  events_out->Branch("njets",&njets, "njets/s");             // njets
  events_out->Branch("l1l2Pt",&l1l2Pt, "l1l2Pt/F");          
  events_out->Branch("l1l2M",&l1l2M, "l1l2M/F");          
  events_out->Branch("t1vHT",&t1vHT, "t1vHT/F");
  events_out->Branch("t1met",&t1met, "t1met/F");
  events_out->Branch("resp",&resp, "resp/F");
  events_out->Branch("jzbraw",&jzbraw, "jzbraw/F");
  events_out->Branch("sel_basic",&sel_basic, "sel_basic/O");
  events_out->Branch("l1jl2jPt",&l1jl2jPt, "l1jl2jPt/F");
  events_out->Branch("l1jl2jM",&l1jl2jM, "l1jl2jM/F");
  events_out->Branch("recoilmet_lep",&recoilmet_lep, "recoilmet_lep/F");
  events_out->Branch("recoiljzb_lep",&recoiljzb_lep, "recoiljzb_lep/F");
  events_out->Branch("MHT",&MHT, "MHT/F");
  events_out->Branch("rvHT",&rvHT, "rvHT/F");

  // --- event loop
  ULong64_t nentries = events_in->GetEntries();

  if(goFast)nentries = nentries*0.01;
  for (ULong64_t ii=0; ii<nentries;ii++)
  {
     events_in->GetEntry(ii);

     l1l2Pt = l1l2Pt_; 
     l1l2M = l1l2M_; 

     resp = pro_ij0->Interpolate(l1l2Pt);      
     
     sel_basic = false;
     if(nleps_>=2 && lepPt_[0]>25 && lepPt_[1]>20 && l1l2M_>81 && l1l2M_<101 && (lepID_[0]*lepID_[1]==-11*11 || lepID_[0]*lepID_[1]==-13*13))
     if(!isDYTauTau_)
     sel_basic = true; 

     TLorentzVector l1(0,0,0,0);
     TLorentzVector l2(0,0,0,0);
     TLorentzVector l1j(0,0,0,0); 
     TLorentzVector l2j(0,0,0,0);
     TLorentzVector hadronic_jet_vector(0,0,0,0);

     if(nleps_>=2)
     {
        l1.SetPtEtaPhiM(lepPt_[0], lepEta_[0], lepPhi_[0], lepM_[0]);
        l2.SetPtEtaPhiM(lepPt_[1], lepEta_[1], lepPhi_[1], lepM_[1]);
     } 

     if(nrjets_ == 2) // needs to be improved
     {
        l1j.SetPtEtaPhiM(rjetPt_[0], rjetEta_[0], rjetPhi_[0], rjetM_[0]);
        l2j.SetPtEtaPhiM(rjetPt_[1], rjetEta_[1], rjetPhi_[1], rjetM_[1]);
     }

     l1jl2jPt = (l1j + l2j).Pt();
     l1jl2jM = (l1j + l2j).M();

     for(int ii=0; ii < njets; ++ii)
     {
	TLorentzVector aJet(0,0,0,0);
        aJet.SetPtEtaPhiM(jetPt_[ii], jetEta_[ii], jetPhi_[ii], jetM_[ii]);
        hadronic_jet_vector += aJet;
     }

     recoilmet_lep = (hadronic_jet_vector + l1 + l2).Pt();
     recoiljzb_lep = (hadronic_jet_vector).Pt() - (l1 + l2).Pt();
     
     rvHT = (hadronic_jet_vector).Pt();
     
     for(int ii=0; ii < nrjets_; ++ii)
     {
	TLorentzVector aJet(0,0,0,0);
        aJet.SetPtEtaPhiM(rjetPt_[ii], rjetEta_[ii], rjetPhi_[ii], rjetM_[ii]);
        hadronic_jet_vector += aJet;
     }
     MHT = hadronic_jet_vector.Pt();


     cjzb_resp = t1vHT_ - resp*l1l2Pt_;
     jzbraw = t1vHT_ - l1l2Pt_;
     resp = resp;
     t1vHT = t1vHT_;
     t1met = t1met_;

     njets = njets_;
     l1l2Pt = l1l2Pt_;     

     events_out->Fill();
  }


  cout << "writting the output file: " << filename_out << endl;

  fp_out->cd();
  events_out->Write();
  fp_out->Close();
}

void init()
{
    events_in->SetBranchAddress("goodVtx"          ,&goodVtx_               );
    events_in->SetBranchAddress("nVtx"             ,&nVtx_                  );
    events_in->SetBranchAddress("eventNum"         ,&eventNum_              );
    events_in->SetBranchAddress("runNum"           ,&runNum_                );
    events_in->SetBranchAddress("lumi"             ,&lumi_                  );

    events_in->SetBranchAddress("l1l2M"            ,&l1l2M_                 );
    events_in->SetBranchAddress("l1l2Pt"           ,&l1l2Pt_                );
    events_in->SetBranchAddress("l1l2Eta"          ,&l1l2Eta_               );
    events_in->SetBranchAddress("l1l2Phi"          ,&l1l2Phi_               );
    events_in->SetBranchAddress("l1l2DPhi"         ,&l1l2DPhi_              );
    events_in->SetBranchAddress("l1l2DR"           ,&l1l2DR_                );

    events_in->SetBranchAddress("nleps"            ,&nleps_                 );
    events_in->SetBranchAddress("lepPt"            ,lepPt_                  );
    events_in->SetBranchAddress("lepEta"           ,lepEta_                 );
    events_in->SetBranchAddress("lepPhi"           ,lepPhi_                 );
    events_in->SetBranchAddress("lepM"             ,lepM_                   );
    events_in->SetBranchAddress("lepIso"           ,lepIso_                 );
    events_in->SetBranchAddress("lepPtRel"         ,lepPtRel_               );
    events_in->SetBranchAddress("lepID"            ,lepID_                  );
    events_in->SetBranchAddress("lepGenMatchIndex" ,lepGenMatchIndex_       );
    events_in->SetBranchAddress("lepTriggerMatch"  ,lepTriggerMatch_        );

    events_in->SetBranchAddress("njets"            ,&njets_                 );
    events_in->SetBranchAddress("jetPt"            ,jetPt_                  );
    events_in->SetBranchAddress("jetEta"           ,jetEta_                 );
    events_in->SetBranchAddress("jetPhi"           ,jetPhi_                 );
    events_in->SetBranchAddress("jetM"             ,jetM_                   );
    events_in->SetBranchAddress("jetBTag"          ,jetBTag_                );

    events_in->SetBranchAddress("nrjets"           ,&nrjets_                );
    events_in->SetBranchAddress("rjetPt"           ,rjetPt_                 );
    events_in->SetBranchAddress("rjetEta"          ,rjetEta_                );
    events_in->SetBranchAddress("rjetPhi"          ,rjetPhi_                );
    events_in->SetBranchAddress("rjetM"            ,rjetM_                  );
    events_in->SetBranchAddress("rjetBTag"         ,rjetBTag_               );

    events_in->SetBranchAddress("ngenleps"         ,&ngenleps_              );
    events_in->SetBranchAddress("genlepPt"         ,genlepPt_               );
    events_in->SetBranchAddress("genlepEta"        ,genlepEta_              );
    events_in->SetBranchAddress("genlepPhi"        ,genlepPhi_              );
    events_in->SetBranchAddress("genlepM"          ,genlepM_                );
    events_in->SetBranchAddress("genlepID"         ,genlepID_               );
    events_in->SetBranchAddress("genlepMID"        ,genlepMID_              );
    events_in->SetBranchAddress("genlepGMID"       ,genlepGMID_             );
    events_in->SetBranchAddress("genlepGGMID"      ,genlepGGMID_            );

    events_in->SetBranchAddress("genl1l2M"         ,&genl1l2M_              );
    events_in->SetBranchAddress("genl1l2Pt"        ,&genl1l2Pt_             );
    events_in->SetBranchAddress("genl1l2Eta"       ,&genl1l2Eta_            );
    events_in->SetBranchAddress("genl1l2Phi"       ,&genl1l2Phi_            );
    events_in->SetBranchAddress("genl1l2DPhi"      ,&genl1l2DPhi_           );
    events_in->SetBranchAddress("genl1l2DR"        ,&genl1l2DR_             );

    events_in->SetBranchAddress("nphos"            ,&nphos_                 );
    events_in->SetBranchAddress("met"              ,&met_                   );
    events_in->SetBranchAddress("metPhi"           ,&metPhi_                );
    events_in->SetBranchAddress("genmet"           ,&genmet_                );
    events_in->SetBranchAddress("t1met"            ,&t1met_                 );
    events_in->SetBranchAddress("t1metPhi"         ,&t1metPhi_              );
    events_in->SetBranchAddress("sumEt"            ,&sumEt_                 );
    events_in->SetBranchAddress("t1sumEt"          ,&t1sumEt_               );

    events_in->SetBranchAddress("rho"              ,&rho_                   );
    events_in->SetBranchAddress("nPU"              ,&nPU_                   );
    events_in->SetBranchAddress("nPUTrue"          ,&nPUTrue_               );

    events_in->SetBranchAddress("vHT"              ,&vHT_                   );
    events_in->SetBranchAddress("t1vHT"            ,&t1vHT_                 );
    events_in->SetBranchAddress("jvHT"             ,&jvHT_                  );

    events_in->SetBranchAddress("HLT_e1e2"         ,&HLT_e1e2_              );
    events_in->SetBranchAddress("HLT_mu1mu2"       ,&HLT_mu1mu2_            );
    events_in->SetBranchAddress("HLT_mu1e2"        ,&HLT_mu1e2_             );
    events_in->SetBranchAddress("HLT_e1mu2"        ,&HLT_e1mu2_             );
    events_in->SetBranchAddress("HLT_pfmet"        ,&HLT_pfmet_             );
    events_in->SetBranchAddress("HLT_pfmetCSV"     ,&HLT_pfmetCSV_          );

    events_in->SetBranchAddress("Flag_trackingFailureFilter"                   ,&Flag_trackingFailureFilter_	                );
    events_in->SetBranchAddress("Flag_goodVertices"                            ,&Flag_goodVertices_			                );
    events_in->SetBranchAddress("Flag_CSCTightHaloFilter"                      ,&Flag_CSCTightHaloFilter_		                );
    events_in->SetBranchAddress("Flag_trkPOGFilters"                           ,&Flag_trkPOGFilters_		                );
    events_in->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters"          ,&Flag_trkPOG_logErrorTooManyClusters_	        );
    events_in->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter"      ,&Flag_EcalDeadCellTriggerPrimitiveFilter_           );
    events_in->SetBranchAddress("Flag_ecalLaserCorrFilter"                     ,&Flag_ecalLaserCorrFilter_		                );
    events_in->SetBranchAddress("Flag_trkPOG_manystripclus53X"                 ,&Flag_trkPOG_manystripclus53X_		        );
    events_in->SetBranchAddress("Flag_eeBadScFilter"                           ,&Flag_eeBadScFilter_		                );
    events_in->SetBranchAddress("Flag_METFilters"                              ,&Flag_METFilters_			                );
    events_in->SetBranchAddress("Flag_HBHENoiseFilter"                         ,&Flag_HBHENoiseFilter_			        );
    events_in->SetBranchAddress("Flag_trkPOG_toomanystripclus53X"              ,&Flag_trkPOG_toomanystripclus53X_	                );
    events_in->SetBranchAddress("Flag_hcalLaserEventFilter"                    ,&Flag_hcalLaserEventFilter_		                );

    events_in->SetBranchAddress("isDYTauTau"         ,&isDYTauTau_              );

    events_in->SetBranchAddress("ngenparts"             ,&ngenparts_                );
    events_in->SetBranchAddress("genpartPt"             ,genpartPt_                 );
    events_in->SetBranchAddress("genpartEta"            ,genpartEta_                );
    events_in->SetBranchAddress("genpartPhi"            ,genpartPhi_                );
    events_in->SetBranchAddress("genpartM"              ,genpartM_                  );
    events_in->SetBranchAddress("genpartID"             ,genpartID_                 );
    events_in->SetBranchAddress("genpartDID1"           ,genpartDID1_               );
    events_in->SetBranchAddress("genpartDID2"           ,genpartDID2_               );
    events_in->SetBranchAddress("genpartDRMI1"          ,genpartDRMI1_              );
    events_in->SetBranchAddress("genpartDRMI2"          ,genpartDRMI2_              ); 
}
