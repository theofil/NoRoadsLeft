// https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
// in synch with CMG https://github.com/CERN-PH-CMG/cmg-cmssw/blob/CMGTools-from-CMSSW_7_2_3/CMGTools/TTHAnalysis/python/samples/samples_13TeV_PHYS14.py
#ifndef CrossSections13TeV_h
#define CrossSections13TeV_h
TCut xs_DYJetsToLL("2008.*3");                            //DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU4bx50_PHYS14_25_V1-v1/MINIAODSIM
TCut xs_TTJets("809.1");                                  //TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU4bx50_PHYS14_25_V1-v1/MINIAODSIM
TCut xs_WZJetsTo3LNu("2.165");                            //WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
TCut xs_ZZ4L("31.8*(3*0.03366*0.03366)");                 //ZZTo4L_Tune4C_13TeV-powheg-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
#endif
