#include "TH1.h"
#include "TF1.h"

{
    TFile *f1=new TFile("../histogram/efficiency/before/13C_efficiency_single_particle.root");
    
    Double_t temp=0;
    Double_t count1 = 60000;
    Double_t sigma1 = 40;
    Double_t centroid1 = 3323;
    Double_t count2 = 61000;
    Double_t sigma2 = 39;
    Double_t centroid2 = 3434;
    Double_t background = 928;
    
    ADC_eff_1_ch15->Draw();
    
    TRandom3 r;
    TH1F* h1 = new TH1F("h1","test1",4000,0,5000);
    for (int i=0;i<count1;i++){
       h1->Fill(r.Gaus(centroid1,sigma1));
    }
    for (int i=0;i<10000;i++){
        temp=h1->GetBinContent(i);
        h1->SetBinContent(i, temp+background);
    }
    TH1F* h2 = new TH1F("h2","test2",4000,0,5000);
    for (int i=0;i<count2;i++){
        h2->Fill(r.Gaus(centroid2,sigma2));
    }
    for (int i=0;i<10000;i++){
        temp=h2->GetBinContent(i);
        h2->SetBinContent(i, temp+background);
    }
    TH1F* h3 = new TH1F("h3","test3",4000,0,5000);
    for (int i=0;i<10000;i++){
        temp=h3->GetBinContent(i);
        h3->SetBinContent(i, temp-background);
    }
    h3->Add(h1,1.0);
    h3->Add(h2,1.0);
    
    h3->SetLineColor(1);
    h2->SetLineColor(2);
    h1->SetLineColor(3);
    h1->Draw("same");
    h2->Draw("same");
    h3->Draw("same");
}




