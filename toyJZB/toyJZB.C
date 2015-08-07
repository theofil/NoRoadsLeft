#include "../core/utilities.h"
#include "../toyJZB/initPtSpectrum.C"
#include "../core/autoRebin.h"

bool goFast(true);
unsigned long nToys = 100000000;
TRandom3 qTGenerator(122562123);

TH1F *qTSpectrum;
TH1F *ZpTSpectrum;
TH1F *jpTSpectrum;
TH1F *jzbPos;
TH1F *jzbNeg;

TF1 sigmaOverE("sigmaOverE", "sqrt(([0]**2)/x + [1]**2)",0, 500);
TF1 jres("jres","sqrt([0]*x + ([1]**2)*(x**2))",0, 500);
TF1 lowPtJetBias("lowPtJetBias","[0] + [1]*x + [2]*x**2", 0,55); // valid up to 0.55, after that is 1:1
TF1 highPtJetBias("histPtJetBias", "[0] + [1]*x", 55, 500);

/* lowPtJetBias  
 * p0                        =      21.6225   +/-   0.0965604   
 * p1                        =     0.213766   +/-   0.00949187  
 * p2                        =   0.00661534   +/-   0.000157606 

 * highPtJetBias
 * p0                        =      2.56538   +/-   0.757161    
 * p1                        =      0.92931   +/-   0.00810402  
 */

TF1 *f1;

// --- parameteres 
float stochastic = 1.;          
float constant   = 0.2;
float jetbias    = 1.0; 
float Zres       = 0.25;  // 0.2 is the RMS of the true MC resolution for genl1l2Pt < 10; is not Gaussian, more picky and with far larger tails tune this to 0.25
float ZpT_thresh = 5; // nominal value = 5
float qT_thresh  = 3; // nominal value = 3

bool  applyLowPtJetBias       = false;
bool  applyHighPtJetBias      = false;
bool  applyFlatFractionalBias = false; // 0.94 flat
bool  applyJetPtThreshold     = true;  //
float JpT_thresh              = 5;    //

void toyJZB()
{
    cout << "jzb toys are starting, have fun" << endl;
    setTDRStyle();
    gStyle->SetOptTitle(0);
    initPtSpectrum(); // loads external pt spectrum from file

    if(goFast) {nToys = 1000000; cout << "goFast with nToys = " << nToys << endl;}


    qTSpectrum   = new TH1F("qTSpectrum",";qT [GeV]; Events / GeV", 100, 0, 100);
    ZpTSpectrum  = new TH1F("ZpTSpectrum",";ZpT [GeV]; Events / GeV", 100, 0, 100);
    jpTSpectrum  = new TH1F("jpTSpectrum",";jpT [GeV]; Events / GeV", 100, 0, 100);
    jzbPos       = new TH1F("jzbPos","JZB > 0;|JZB| [GeV]; Events / GeV", 100, 0, 100);
    jzbNeg       = new TH1F("jzbNeg","JZB < 0;|JZB| [GeV]; Events / GeV", 100, 0, 100);

    // --- PDF parameters
    sigmaOverE.SetParameter(0, stochastic);
    sigmaOverE.SetParameter(1, constant);

    jres.SetParameter(0, stochastic);
    jres.SetParameter(1, constant);

    lowPtJetBias.SetParameter(0, 21.6225); 
    lowPtJetBias.SetParameter(1, 0.213766); 
    lowPtJetBias.SetParameter(2, 0.00661534); 

    highPtJetBias.SetParameter(0, 2.56538); 
    highPtJetBias.SetParameter(1, 0.92931); 


    f1 = new TF1("f1","[0]*x**[1]", qT_thresh, 500);
    f1->SetParameter(0 , 1.41883e+00);
    f1->SetParameter(1,  -1.8e+00);

//    h_ptSpectrum->Rebin(10); float qT    =  h_ptSpectrum->GetRandom();
    
    for(unsigned long entry = 0; entry < nToys; entry++)
    {
	float qT    =  f1->GetRandom();
        float sigma = jres.Eval(qT);
        float ZpT   = qTGenerator.Gaus(qT, Zres*qT); // emulates the leptonic turn-on 

        if(qT > 0 && applyLowPtJetBias && qT < 55)  jetbias = lowPtJetBias.Eval(qT)/qT;
        if(qT > 55 && applyHighPtJetBias) jetbias = highPtJetBias.Eval(qT)/qT;
        if(applyFlatFractionalBias) jetbias = 0.94;

        float jpT     = qTGenerator.Gaus(jetbias*qT, sigma);

        if(jpT<0) jpT = 0;
        if(ZpT<0) ZpT = 0;

        float jzb = jpT - ZpT;
        bool selEvent(true); 

        selEvent = false; 
        if(ZpT > ZpT_thresh) selEvent=true;
        if(applyJetPtThreshold && jpT < JpT_thresh) selEvent = false;

        if(selEvent) 
        {
          qTSpectrum   -> Fill(qT);
          ZpTSpectrum  -> Fill(ZpT);
          jpTSpectrum  ->Fill(jpT);
          if(jzb > 0) jzbPos->Fill(fabs(jzb));
          if(jzb < 0) jzbNeg->Fill(fabs(jzb));
	}
    }
    cout << "selection efficiency = " << float(jpTSpectrum->GetEntries())/float(nToys) << endl;

    TCanvas *c_sigmaOverE = new TCanvas("c_sigmaOverE","c_sigmaOverE");
    sigmaOverE.GetXaxis()->SetTitle("qT [GeV]"); 
    sigmaOverE.GetYaxis()->SetTitle("#sigma(qT)/qT"); 
    sigmaOverE.Draw();

    TCanvas *c_qTSpectrum = new TCanvas("c_qTSpectrum","c_qTSpectrum");
    c_qTSpectrum->SetLogy();
    qTSpectrum->Draw("hist");

    TCanvas *c_ZpTSpectrum = new TCanvas("c_ZpTSpectrum","c_ZpTSpectrum");
    c_ZpTSpectrum->SetLogy();
    ZpTSpectrum->Draw("hist");

    TCanvas *c_jpTSpectrum = new TCanvas("c_jpTSpectrum","c_jpTSpectrum");
    c_jpTSpectrum->SetLogy();
    jpTSpectrum->Draw("hist");

    TCanvas *c_jzb = new TCanvas("c_jzb","c_jzb");
    c_jzb->SetLogy();
    jzbPos->SetLineColor(kBlue);
    jzbNeg->SetLineColor(kRed);
    jzbPos->Draw("hist");
    jzbNeg->Draw("hist same");

    gStyle->SetErrorX(0.5);
 
    TCanvas *c_jzbRatio = new TCanvas("c_jzbRatio","c_jzbRatio");
    TH1F *jzbRatio = doRatio(jzbPos, jzbNeg, 100);
//    TH1F *jzbRatio = (TH1F*) jzbPos->Clone();
//    jzbRatio->Divide(jzbNeg);
    jzbRatio->SetLineColor(kBlack); 
    jzbRatio->GetYaxis()->SetRangeUser(0.6, 1.4);
    jzbRatio->Draw();

}

