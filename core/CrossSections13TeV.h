// https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
// in synch with CMG https://github.com/CERN-PH-CMG/cmg-cmssw/blob/CMGTools-from-CMSSW_7_2_3/CMGTools/TTHAnalysis/python/samples/samples_13TeV_PHYS14.py
// sync also with https://twiki.cern.ch/twiki/bin/view/CMS/DesyTauAnalysesRun2?rev=25
#ifndef CrossSections13TeV_h
#define CrossSections13TeV_h
TCut xs_Data("1");                                        // trivial
TCut xs_WW("63.21");                                      // WW_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1
TCut xs_WZ("22.82");                                      // WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1
TCut xs_ZZ("10.32");                                      // ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3
TCut xs_WWTo2L2Nu("10.481 ");                             // WWTo2L2Nu_13TeV-powheg/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1

TCut xs_DYJetsToLL("6025.2");                            //DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU4bx50_PHYS14_25_V1-v1/MINIAODSIM
TCut xs_DYJetsToLLM50("6025.2");                         //DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU4bx50_PHYS14_25_V1-v1/MINIAODSIM
TCut xs_DYJetsM1050("18610");                            //DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU4bx50_PHYS14_25_V1-v1/MINIAODSIM
TCut xs_TTJets("831.76");                                //TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU4bx50_PHYS14_25_V1-v1/MINIAODSIM
TCut xs_WZJetsTo3LNu("2.165");                           //WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
TCut xs_ZZ4L("31.8*(3*0.03366*0.03366)");                //ZZTo4L_Tune4C_13TeV-powheg-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIMM

float xs_float_Data(1);                                        // trivial
float xs_float_WW(63.21);                                      // WW_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1
float xs_float_WZ(22.82);                                      // WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1
float xs_float_ZZ(10.32);                                      // ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3
float xs_float_WWTo2L2Nu(10.481 );                             // WWTo2L2Nu_13TeV-powheg/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1

float xs_float_DYJetsToLL(6025.2);                            //DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU4bx50_PHYS14_25_V1-v1/MINIAODSIM
float xs_float_DYJetsToLLM50(6025.2);                         //DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU4bx50_PHYS14_25_V1-v1/MINIAODSIM
float xs_float_DYJetsM1050(18610);                            //DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU4bx50_PHYS14_25_V1-v1/MINIAODSIM
float xs_float_TTJets(831.76);                                //TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU4bx50_PHYS14_25_V1-v1/MINIAODSIM
float xs_float_WZJetsTo3LNu(2.165);                           //WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
float xs_float_ZZ4L(31.8*(3*0.03366*0.03366));                //ZZTo4L_Tune4C_13TeV-powheg-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIMM
#endif
