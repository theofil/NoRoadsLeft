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

    TCut Lumi("42"); // 1000 pb-1
    if(goFast) Lumi = Lumi*TCut("0.01");
    cout << "loading v1_drivers" << endl;
    cout << "goFast = " << goFast << endl;
    cout << "Lumi = " << Lumi.GetTitle() << endl;

/*
    mcDriver.push_back(new SimpleSample(fp_ZZ           , "ZZ"                    ,  Lumi                           ,goFast, kGray,             kGray)); 
    mcDriver.push_back(new SimpleSample(fp_WZ           , "WZ"                    ,  Lumi                           ,goFast, 37,                   37)); 
    mcDriver.push_back(new SimpleSample(fp_WW2l2nu      , "WW2l2nu"               ,  Lumi                           ,goFast, 40,                   40)); 
    mcDriver.push_back(new SimpleSample(fp_TTJets       , "t#bar{t}/WW"           ,  Lumi                           ,goFast, 40,                   40)); 
    mcDriver.push_back(new SimpleSample(fp_DYJetsM1050  , "DY(#tau#tau)"  , TCut("isDYTauTau ? 1:0")*Lumi   ,goFast, 47,                   47)); 
    mcDriver.push_back(new SimpleSample(fp_DYJetsM50    , "DY(#tau#tau)"  , TCut("isDYTauTau ? 1:0")*Lumi   ,goFast, 47,                   47)); 
    mcDriver.push_back(new SimpleSample(fp_DYJetsM50    , "DY(#mu#mu,ee)" , TCut("!isDYTauTau ? 1:0")*Lumi  ,goFast, kWhite , kRed+1         )); 
    mcDriver.push_back(new SimpleSample(fp_DYJetsM1050  , "DY(#mu#mu,ee)" , TCut("!isDYTauTau ? 1:0")*Lumi  ,goFast, kWhite , kRed+1         )); 
*/

    mcDriver.push_back(new SimpleSample(fp_TTJets       , "t#bar{t}"              , Lumi                           ,goFast, 40,                   40)); 
    mcDriver.push_back(new SimpleSample(fp_WW2l2nu      , "WW"                    , Lumi                           ,goFast, 47,                   47)); 
    mcDriver.push_back(new SimpleSample(fp_WZ           , "WZ"                    , Lumi                           ,goFast, 37,                   37)); 
    mcDriver.push_back(new SimpleSample(fp_ZZ           , "ZZ"                    , Lumi                           ,goFast, kGray,             kGray)); 
    mcDriver.push_back(new SimpleSample(fp_DYJetsM1050  , "DY"                    , Lumi                           ,goFast, kWhite , kWhite        )); 
    mcDriver.push_back(new SimpleSample(fp_DYJetsM50    , "DY"                    , Lumi                           ,goFast, kWhite , kRed+1         )); 

}
