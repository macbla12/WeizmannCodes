void PlottingFCALRatio(){
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(000000);
    
    // 1. Otwarcie pliku ROOT
    TFile *DataFile = new TFile("Plots/ReproducedPlotHighStat.root", "READ");
    
    // Pobranie histogramów z pliku
    TH1D *h_FCAL      = (TH1D*)DataFile->Get("FCAL_hist");
    TH1D *h_0n0n      = (TH1D*)DataFile->Get("FCAL0n0noverFCALmb");
    TH1D *h_1n0n      = (TH1D*)DataFile->Get("FCAL1n0noverFCALmb");
    TH1D *h_1n1n      = (TH1D*)DataFile->Get("FCAL1n1noverFCALmb");

    // 2. Klonowanie histogramów i dzielenie przez FCAL_hist
    TH1D *r_0n0n = (TH1D*)h_0n0n->Clone("r_0n0n");
    TH1D *r_1n0n = (TH1D*)h_1n0n->Clone("r_1n0n");
    TH1D *r_1n1n = (TH1D*)h_1n1n->Clone("r_1n1n");

    // Poprawne wyliczenie błędów i dzielenie
    r_0n0n->Sumw2();
    r_1n0n->Sumw2();
    r_1n1n->Sumw2();

    r_0n0n->Divide(h_FCAL);
    r_1n0n->Divide(h_FCAL);
    r_1n1n->Divide(h_FCAL);

    // 3. Tworzenie kanwy i rysowanie
    TCanvas *c1 = new TCanvas("c1", "FCAL Ratios", 800, 600);
    
    // Ustawienie formatowania osi i stylów dla pierwszego histogramu
    r_0n0n->SetTitle(";FCal E_{T} [TeV];Ratio to Minimum Bias"); 
    r_0n0n->SetLineColor(kRed);
    r_0n0n->SetLineWidth(2);
    r_0n0n->SetMinimum(0.0); // Opcjonalnie: zmień zakres osi Y
    r_0n0n->SetMaximum(1.2); 
    r_0n0n->GetXaxis()->SetRangeUser(0,0.45);
    r_0n0n->Draw("HIST");

    r_1n0n->SetLineColor(kBlue);
    r_1n0n->SetLineWidth(2);
    r_1n0n->GetXaxis()->SetRangeUser(0,0.45);
    r_1n0n->Draw("HIST SAME");

    r_1n1n->SetLineColor(kGreen+2);
    r_1n1n->SetLineWidth(2);
    r_1n1n->GetXaxis()->SetRangeUser(0,0.45);
    r_1n1n->Draw("HIST SAME");

    // 4. Legenda
    TLegend *leg1 = new TLegend(0.60, 0.65, 0.85, 0.82);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.035); 
    leg1->SetBorderSize(0);  
    leg1->SetFillStyle(0);   
    leg1->AddEntry(r_0n0n, "0n0n / FCal", "l");
    leg1->AddEntry(r_1n0n, "1n0n / FCal", "l");
    leg1->AddEntry(r_1n1n, "1n1n / FCal", "l");
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
    c1->SaveAs("Plots/FCAL_Ratios.png");
    c1->SaveAs("Plots/FCAL_Ratios.pdf");
    c1->SaveAs("Plots/FCAL_Ratios.root");
}