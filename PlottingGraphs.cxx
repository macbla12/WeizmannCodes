
void PlottingGraphs()
{
    TFile *inputFile = TFile::Open("Plots/ReproducedPlot.root");
//Npart
    TGraph *graph = dynamic_cast<TGraph*>(inputFile->Get("g_Tracks_vs_Npart"));
    graph->SetMarkerStyle(20);        
    graph->SetMarkerSize(1.2);       
    graph->SetMarkerColor(kBlue);  
    graph->SetLineColor(kBlue);    
    graph->SetTitle(";<N_{part}>;<1/<N_{part}> dn/d#eta");
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    graph->GetYaxis()->SetRangeUser(0, 8);
    graph->GetXaxis()->SetLimits(0.0, 36.0);
    graph->Draw("APL");
    c1->SetGrid();
    c1->SaveAs("Plots/Tracks_vs_Npart.pdf");


//TAA
    TGraph *graphTAA = dynamic_cast<TGraph*>(inputFile->Get("g_Tracks_vs_TAA"));
    graphTAA->SetMarkerStyle(20);        
    graphTAA->SetMarkerSize(1.2);       
    graphTAA->SetMarkerColor(kBlack);  
    graphTAA->SetLineColor(kBlack);    
    graphTAA->SetTitle(";<T_{AA}>;<1/<T_{AA}> dn/d#eta>");
    TCanvas *c2 = new TCanvas("c2", "Canvas", 800, 600);
    graphTAA->GetYaxis()->SetRangeUser(0, 50);
    c2->SetGrid();
    graphTAA->GetXaxis()->SetLimits(0.0, 0.7);
    graphTAA->Draw("APL");
    c2->SaveAs("Plots/Tracks_vs_TAA.pdf");


    inputFile->Close();
}