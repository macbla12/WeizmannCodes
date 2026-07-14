void PlottingZDC(){
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(000000);
    
    TFile *DataFile = new TFile("Plots/FCaL-CentralityHighStat.root","READ");
    TH1D *ZDC_ASide= (TH1D*)DataFile->Get("ZDC_ASide_hist");
    
    TCanvas *c1 = new TCanvas("c1", "ZDC", 800, 600);
    
    // Rysujemy histogram
    ZDC_ASide->SetTitle(";ZDC Energy [PeV];Events"); 
    ZDC_ASide->SetLineWidth(2); 
    ZDC_ASide->GetYaxis()->SetRangeUser(0, 500000);
    ZDC_ASide->Draw("HIST");

    TLine *boundryLine = new TLine(0.001, 0, 0.001, 500000);
    boundryLine->SetLineColor(kRed); // Ciemniejszy szary, żeby nie zniknął na białym tle
    boundryLine->SetLineStyle(3);       // Kod stylu dla linii kropkowanej/przerywanej (dashed/dotted)
    boundryLine->SetLineWidth(3);       // Grubość linii
    boundryLine->Draw("same");


    TLatex *text1 = new TLatex();
        text1->SetNDC();
        text1->SetTextSize(0.04);
        text1->DrawLatex(0.52, 0.84, "ATLAS Internal 2026 O+O");
        text1->DrawLatex(0.52, 0.78, "Run - 501859 #sqrt{s_{NN}} = 5.36 TeV");


    c1->Modified();
    c1->Update();

    c1->SaveAs("Plots/ZDC.png");
    c1->SaveAs("Plots/ZDC.pdf");
    c1->SaveAs("Plots/ZDC.root");

    
    gPad->SetLogy(0);
}