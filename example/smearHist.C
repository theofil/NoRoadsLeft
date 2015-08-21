#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "../core/SimpleCanvas.h"

TFile *fp;
TH1F *h1, *h2, *h3, *h4, *h5;

#define nMax 1000000

void normHist(TH1F *hist)
{
    hist->Scale(1.0/hist->Integral());
}

void smearHist()
{
    fp = new TFile("tmp.root","OPEN");
    h1 = (TH1F*)fp->Get("h1");
    h2 = (TH1F*)fp->Get("h2");

    h3 = (TH1F*)h1->Clone();
    h3->Reset();

    h4 = (TH1F*)h1->Clone();
    h4->Reset();

   for(int n = 0; n < nMax; n++)
   {
     h3->Fill( (h1->GetRandom() + h1->GetRandom()) );
     h4->Fill( (h1->GetRandom() - h1->GetRandom()) );
   }
   
   normHist(h3);
   normHist(h4);

   h3->SetLineColor(kRed);
   h3->SetLineStyle(2);
   h3->SetLineWidth(2);

   h4->SetLineColor(kBlue);
   h4->SetLineStyle(2);
   h4->SetLineWidth(3);
 //  h4->Fit("gaus");

   h5 = (TH1F*)h4->Clone();
   h5->Add(h3);

   h5->SetLineColor(kGray);
   h5->SetLineStyle(3);
   h5->SetLineWidth(5);

   TCanvas *can1 = new TCanvas("can1","can1");  
   can1->SetLogy();
   h3->DrawNormalized("hist same");  
   h4->DrawNormalized("hist same");  
   h5->DrawNormalized("hist same");  
   h2->Draw("same e1");
 
/*
   hratio = (TH1F*) h2->Clone();
   hratio->Divide(h3);

   TCanvas *can2 = new TCanvas("can2","can2");  
   hratio->GetYaxis()->SetRangeUser(0,2);
   hratio->Draw();  
*/

}
