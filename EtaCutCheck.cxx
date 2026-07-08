void EtaCutCheck()
{
    TString File ="Files/dataOOwithTracks.root"; 
    
    TChain *mychain = new TChain("Z_analysis");
        mychain->Add(File);
    TTreeReader tree_reader(mychain);

    TFile *ZOutput = new TFile("Plots/EtaCheck.root","RECREATE");

    TTreeReaderArray<float> TrackEta(tree_reader, "trk_eta");
    TTreeReaderArray<vector<int>> TrackQuality(tree_reader, "trk_qualities");
    
    TH1D *EtaDistribution= new TH1D("EtaDistribution","EtaDistribution",100,-2.8,2.8);
    TH1D *QualityDistribution= new TH1D("QualityDistribution","QualityDistribution",5,-.5,4.5);


    Long64_t total_events = mychain->GetEntries();
    cout<<"============================="<<endl;
    cout<<"Starting Analysis"<<endl;

    int eventID=0;
    double EtaSelection=0, FullEta=0;
    while(tree_reader.Next()){
        eventID++;
        if (eventID % 10000 == 0) {
            std::cout << Form("\rDone: %.2f", (100.0 * eventID / total_events)) <<"%   "<< std::flush;
        }
        for(int track=0; track<TrackEta.GetSize();track++)
        {
            if(abs(TrackEta[track])<2.5 && TrackQuality[track][0]){
                FullEta++;
                EtaDistribution->Fill(TrackEta[track]);
                if(abs(TrackEta[track])<0.5) EtaSelection++; 
            } 
        }      

    }

    EtaDistribution->Write();
    QualityDistribution->Write();
    ZOutput->Close();


    cout<<"FullEta="<<FullEta<<"    EtaSelection="<<EtaSelection<<endl;
    cout<<"EtaCut="<<EtaSelection/FullEta<<endl;
    cout<<"============================="<<endl;

}