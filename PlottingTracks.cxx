void PlottingTracks(){
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(000000);
    
    TFile *DataFile = new TFile("Plots/FCaL-CentralityHighStat.root","READ");
    TH1D *Track_Full= (TH1D*)DataFile->Get("Tracks_fullspectrum_hist");
    TH1D *Track_Select= (TH1D*)DataFile->Get("Tracks_selectedspectrum_hist");

    
    TCanvas *c1 = new TCanvas("c1", "Tracks", 800, 600);
    
    // Rysujemy histogram
    Track_Full->SetTitle(";Number of Tracks;Events"); 
    Track_Full->SetLineWidth(2); 
    Track_Full->Draw("HIST");

    Track_Select->SetLineColor(kRed); 
    Track_Select->SetLineWidth(2); 
    Track_Select->Draw("HIST SAME"); 

    
    gPad->SetLogy(1);
    
    
    TLegend *leg1 = new TLegend(0.63, 0.63, 0.83, 0.75);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.05); 
    leg1->SetBorderSize(0);  
    leg1->SetFillStyle(0);   
    leg1->AddEntry(Track_Full, "Full Spectrum", "l");
    leg1->AddEntry(Track_Select, "Zero-neutron", "l");

    leg1->Draw();

    TLatex *text1 = new TLatex();
        text1->SetNDC();
        text1->SetTextSize(0.04);
        text1->DrawLatex(0.52, 0.84, "ATLAS Internal 2026 O+O");
        text1->DrawLatex(0.52, 0.78, "Run - 501859 #sqrt{s_{NN}} = 5.36 TeV");

    TLatex *text2 = new TLatex();
        text2->SetNDC();
        text2->SetTextSize(0.04);
        text2->DrawLatex(0.72, 0.59, "p_{T}>300 MeV");
        text2->DrawLatex(0.78, 0.54, "|#eta|<0.5");


    c1->Modified();
    c1->Update();

    c1->SaveAs("Plots/Tracks.png");
    c1->SaveAs("Plots/Tracks.pdf");
    c1->SaveAs("Plots/Tracks.root");

    
    gPad->SetLogy(0);
}