void FittingPeak(){
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(000000);
    
    TFile *DataFile = new TFile("Plots/PeakHighStatistic.root","READ");
    TH1D *FCAL_ZeroZDC_hist= (TH1D*)DataFile->Get("FCAL_ZeroZDC_hist");
    TF1 *Fitted = new TF1("Fitted", "gaus(0) + expo(3)", 0.0, 0.5);

    Fitted->SetParameter(0, 25000.0);
    Fitted->SetParLimits(0, 10000.0, 50000.0); 

    Fitted->FixParameter(1, 0.265);      

    Fitted->FixParameter(2, 0.053);

    Fitted->SetParameter(3, 14.9);
    Fitted->SetParLimits(3, 13.0, 16.0);

    Fitted->SetParameter(4, -60.0);
    Fitted->SetParLimits(4, -80.0, -20.0);  

    FCAL_ZeroZDC_hist->Fit(Fitted, "RN"); 

    TF1 *Gauss = new TF1("Gauss", "gaus", 0, 0.5);
    Gauss->FixParameter(0,Fitted->GetParameter(0));
    Gauss->FixParameter(1,Fitted->GetParameter(1));
    Gauss->FixParameter(2,Fitted->GetParameter(2));

    TF1 *Exponent = new TF1("Exponent", "expo", 0, 0.5);
    Exponent->FixParameter(0,Fitted->GetParameter(3));
    Exponent->FixParameter(1,Fitted->GetParameter(4));

    TCanvas *c1=new TCanvas("c1","c1",800,600);
    gPad->SetLogy(1);
    FCAL_ZeroZDC_hist->Draw();
    Gauss->SetLineColor(kBlack);

    Gauss->SetLineColor(kBlack);
    Exponent->SetLineColor(kBlack);
    Fitted->SetLineWidth(3);

    Fitted->Draw("same");
    Gauss->Draw("same");
    Exponent->Draw("same");

    

    c1->SaveAs("Plots/PeakWithFit.png");
    c1->SaveAs("Plots/PeakWithFit.pdf");
    c1->SaveAs("Plots/PeakWithFit.root");
    gPad->SetLogy(0);


    TF1 *Difference = new TF1("Difference", "gaus(0) - expo(3)", 0.1, 0.4);
    Difference->SetParameter(0, Fitted->GetParameter(0));
    Difference->SetParameter(1, Fitted->GetParameter(1));
    Difference->SetParameter(2, Fitted->GetParameter(2));
    Difference->SetParameter(3, Fitted->GetParameter(3));
    Difference->SetParameter(4, Fitted->GetParameter(4));
    double intersection_x = Difference->GetX(0.0, 0.10, 0.15);
    double intersection_y = Gauss->Eval(intersection_x);
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Punkt przeciecia krzywych:" << std::endl;
    std::cout << "X = " << intersection_x << std::endl;
    std::cout << "Y = " << intersection_y << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

}