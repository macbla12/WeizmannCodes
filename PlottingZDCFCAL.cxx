void PlottingZDCFCAL(){
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(000000);
    
    // 1. Otwarcie pliku ROOT
    TFile *DataFile = new TFile("Plots/ReproducedPlotHighStat.root", "READ");
    
    // Pobranie 3 histogramów ZDC z pliku
    TH1D *h_ZDC       = (TH1D*)DataFile->Get("ZDC_hist");
    TH1D *h_ZDC_1pr   = (TH1D*)DataFile->Get("ZDC_centrality1pr");
    TH1D *h_ZDC_10pr  = (TH1D*)DataFile->Get("ZDC_centrality10pr");

    // 2. Tworzenie kanwy
    TCanvas *c1 = new TCanvas("c1", "ZDC Spectra", 800, 600);
    
    // Włączenie skali logarytmicznej na osi Y (przydatne przy widmach ZDC/etykietach)
    gPad->SetLogy(1);

    // 3. Stylizacja i rysowanie histogramów
    h_ZDC->SetTitle(";ZDC Energy [TeV];Events"); 
    h_ZDC->SetLineColor(kBlack);
    h_ZDC->SetLineWidth(2);
    h_ZDC->Scale(1./h_ZDC->GetEntries());
    h_ZDC->Draw("HIST");

    h_ZDC_1pr->SetLineColor(kRed);
    h_ZDC_1pr->SetLineWidth(2);
    h_ZDC_1pr->Scale(1./h_ZDC_1pr->GetEntries());
    h_ZDC_1pr->Draw("HIST SAME");

    h_ZDC_10pr->SetLineColor(kBlue);
    h_ZDC_10pr->SetLineWidth(2);
    h_ZDC_10pr->Scale(1./h_ZDC_10pr->GetEntries());
    h_ZDC_10pr->Draw("HIST SAME");

    // 4. Legenda
    TLegend *leg1 = new TLegend(0.62, 0.65, 0.85, 0.82);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.035); 
    leg1->SetBorderSize(0);  
    leg1->SetFillStyle(0);   
    leg1->AddEntry(h_ZDC,      "ZDC Inclusive", "l");
    leg1->AddEntry(h_ZDC_1pr,  "Centrality 1%", "l");
    leg1->AddEntry(h_ZDC_10pr, "Centrality 10%", "l");
    leg1->Draw();

    // 5. Adnotacje ATLAS
    TLatex *text1 = new TLatex();
    text1->SetNDC();
    text1->SetTextSize(0.04);
    text1->DrawLatex(0.20, 0.84, "ATLAS Internal 2026 O+O");
    text1->DrawLatex(0.20, 0.78, "Run - 501859 #sqrt{s_{NN}} = 5.36 TeV");

    c1->Modified();
    c1->Update();

    // Zapis do plików
    c1->SaveAs("Plots/ZDC_Spectra.png");
    c1->SaveAs("Plots/ZDC_Spectra.pdf");
    c1->SaveAs("Plots/ZDC_Spectra.root");

    gPad->SetLogy(0);
}