void PlottingBining(){
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(000000);
    
    TFile *DataFile = new TFile("Plots/FCaL-CentralityHighStat.root","READ");
    TH1D *FCAL_ZeroZDC_hist= (TH1D*)DataFile->Get("FCAL_ZeroZDC_hist");
    
    TCanvas *c1 = new TCanvas("c1", "Centrality Bins", 800, 600);
    
    // Rysujemy histogram
    FCAL_ZeroZDC_hist->SetTitle("FCal - zero nuetrons;FCAL [TeV];Events"); 

    FCAL_ZeroZDC_hist->Draw("HIST"); 
    
    gPad->SetLogy(1);
    
    double yMin = FCAL_ZeroZDC_hist->GetMinimum();
    if (yMin <= 0) yMin = 0.5;

    double FcalBoundry = 0.160;
    double Centrality[] = {299.736, 277.788, 262.923, 251.217, 241.358, 232.613, 224.684, 217.391, 210.599, 204.217, 198.159, 192.374, 186.851, 181.527, 176.384, 171.405, 166.613, 161.958, 157.427, 153.033};
    int nLines = sizeof(Centrality) / sizeof(double);
    TLegend *leg1 = new TLegend(0.50, 0.63, 0.68, 0.75);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.05); 
    leg1->SetBorderSize(0);  
    leg1->SetFillStyle(0);   
   
    for(int i = 0; i < nLines; i++) {
        if(FcalBoundry >= Centrality[i]/1000) break;
        double xVal = Centrality[i]/1000;
        
        int bin = FCAL_ZeroZDC_hist->GetXaxis()->FindBin(xVal);
        double yMax = FCAL_ZeroZDC_hist->GetBinContent(bin);
        
        TLine *line = new TLine(xVal, yMin, xVal, yMax);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->Draw("same");
        if(i==0) leg1->AddEntry(line, "Bin boundry - 1%", "l");

    }
    // === TUTAJ DODAJEMY SZARĄ PRZERYWANĄ LINIĘ DLA FcalBoundry ===
    TLine *boundryLine = new TLine(FcalBoundry, 0.2, FcalBoundry, 1.2*FCAL_ZeroZDC_hist->GetMaximum());
    boundryLine->SetLineColor(kGray+2); // Ciemniejszy szary, żeby nie zniknął na białym tle
    boundryLine->SetLineStyle(3);       // Kod stylu dla linii kropkowanej/przerywanej (dashed/dotted)
    boundryLine->SetLineWidth(2);       // Grubość linii
    boundryLine->Draw("same");
    leg1->AddEntry(boundryLine, "Boundry - 0.16 TeV", "l");

    leg1->Draw();

    TLatex *text1 = new TLatex();
        text1->SetNDC();
        text1->SetTextSize(0.04);
        text1->DrawLatex(0.52, 0.84, "ATLAS Internal 2026 O+O");
        text1->DrawLatex(0.52, 0.78, "Run - 501859 #sqrt{s_{NN}} = 5.36 TeV");


    c1->Modified();
    c1->Update();

    c1->SaveAs("Plots/FCALBin.png");
    c1->SaveAs("Plots/FCALBin.pdf");
    c1->SaveAs("Plots/FCALBin.root");

    
    gPad->SetLogy(0);
}