#include "../core/utilities.h"
#include "../core/SimpleDriver.h"
#include "../core/SimpleStack.h"
#include "../core/CrossSections13TeV.h"
#include "../core/simpleROOT_cuts.h"
#include "../core/LocalFilePath.h"

unsigned int totEvents;

SimpleDriver mcDriver, dataDriver;

void v1_drivers();
void v1_drivers(bool);

void v1_drivers()
{
    v1_drivers(false);
}

void v1_drivers(bool goFast)
{

    dataDriver.push_back(new SimpleSample(fp_Data, "Data", TCut("1"), goFast));   // dummy file for pre-DATA period

    TCut Lumi("40.24"); // 1000 pb-1
    if(goFast) Lumi = Lumi*TCut("0.01");
    cout << "loading v1_drivers" << endl;
    cout << "goFast = " << goFast << endl;
    cout << "Lumi = " << Lumi.GetTitle() << endl;

    mcDriver.push_back(new SimpleSample(fp_TTJets       , "t#bar{t}"      , xs_TTJets*Lumi                               ,goFast, 40,            40)); 
    mcDriver.push_back(new SimpleSample(fp_DYJetsM50    , "DY(#tau#tau)"  , TCut("isDYTauTau ? 1:0")*xs_DYJetsToLL*Lumi  ,goFast, 47,            47)); 
    mcDriver.push_back(new SimpleSample(fp_DYJetsM50    , "DY(#mu#mu,ee)" , TCut("!isDYTauTau ? 1:0")*xs_DYJetsToLL*Lumi ,goFast, kWhite , kRed+1         )); 
}
