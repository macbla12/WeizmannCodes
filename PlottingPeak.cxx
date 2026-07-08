void PlottingPeak()
{

    TString File ="data.root"; 
    TChain *mychain = new TChain("Z_analysis");
        mychain->Add(File);
    TTreeReader tree_reader(mychain);

    TFile *ZOutput = new TFile("Plots/Peak.root","RECREATE");

    TTreeReaderArray<float> ZDC_CSide(tree_reader, "ZDC_energy_side_C");
    TTreeReaderArray<float> ZDC_ASide(tree_reader, "ZDC_energy_side_A");

    TTreeReaderArray<float> FCAL(tree_reader, "FCAL");

    TH1D *ZDC_ASide_hist =new TH1D("ZDC_ASide_hist","ZDC_ASide",100,0,0.01);
    TH1D *ZDC_CSide_hist =new TH1D("ZDC_CSide_hist","ZDC_CSide",100,0,0.01);

    TH1D *FCAL_hist =new TH1D("FCAL_hist","FCAL",100,0,0.6);
    TH1D *FCAL_ZeroZDC_hist =new TH1D("FCAL_ZeroZDC_hist","FCAL_ZeroZDC",100,0,0.6);

    Long64_t total_events = mychain->GetEntries();
    cout<<"============================="<<endl;
    cout<<"Starting Analysis"<<endl;
    int eventID=0;
    
    while(tree_reader.Next()){
        if (eventID % 10000 == 0) {
            std::cout << Form("\rDone: %.2f", (100.0 * eventID / total_events)) <<"%   "<< std::flush;
        }
        eventID++;
        ZDC_CSide_hist->Fill(ZDC_CSide[0]);
        ZDC_ASide_hist->Fill(ZDC_ASide[0]);
        FCAL_hist->Fill(FCAL[0]);
        if(ZDC_CSide[0]<0.001 && ZDC_ASide[0]<0.001)
        {
            FCAL_ZeroZDC_hist->Fill(FCAL[0]);
        }
    }
    cout<<endl;

    ZOutput->cd();
    ZDC_CSide_hist->Write();
    ZDC_ASide_hist->Write();

    FCAL_hist->Write();
    FCAL_ZeroZDC_hist->Write();

    ZOutput->Close();
    cout<<"============================="<<endl;

}