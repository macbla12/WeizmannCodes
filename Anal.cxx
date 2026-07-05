void Anal()
{

    double FcalBoundry=0.20;



    Float_t MuonMass=0.1056583755;
    Float_t ElectronMass=0.00051099895;

    TString File ="data.root"; 
    TChain *mychain = new TChain("Z_analysis");
        mychain->Add(File);
    TTreeReader tree_reader(mychain);

    TFile *ZOutput = new TFile("Plots/FCaL-Centrality.root","RECREATE");


    TTreeReaderArray<float> ElRecopT(tree_reader, "e_reco_pT");
    TTreeReaderArray<float> ElRecoEta(tree_reader, "e_reco_eta");
    TTreeReaderArray<float> ElRecoPhi(tree_reader, "e_reco_phi");
    TTreeReaderArray<float> ElRecoQuality(tree_reader, "e_quality");

    TTreeReaderArray<float> MuRecopT(tree_reader, "mu_reco_pT");
    TTreeReaderArray<float> MuRecoEta(tree_reader, "mu_reco_eta");
    TTreeReaderArray<float> MuRecoPhi(tree_reader, "mu_reco_phi");
    TTreeReaderArray<float> MuRecoQuality(tree_reader, "mu_quality");

    TTreeReaderArray<float> ZDC_CSide(tree_reader, "ZDC_energy_side_C");
    TTreeReaderArray<float> ZDC_ASide(tree_reader, "ZDC_energy_side_A");

    TTreeReaderArray<int> afp_n_tracks_station(tree_reader, "afp_n_tracks_station");
    TTreeReaderArray<int> afp_n_clusters_station(tree_reader, "afp_n_clusters_station");
    TTreeReaderArray<int> afp_n_clusters_layer(tree_reader, "afp_n_clusters_layer");
    TTreeReaderArray<float> afp_meanclusterSize_layer(tree_reader, "afp_meanclusterSize_layer");
    TTreeReaderArray<float> afp_meanclusterCharge_layer(tree_reader, "afp_meanclusterCharge_layer");
    TTreeReaderArray<int> afp_n_protons_side(tree_reader, "afp_n_protons_side");
    TTreeReaderArray<float> afp_proton_sum_e_side(tree_reader, "afp_proton_sum_e_side");

    TTreeReaderArray<float> FCAL(tree_reader, "FCAL");
    TTreeReaderArray<float> FCALCside(tree_reader, "FCAL_side_C_from_layers");
    TTreeReaderArray<float> FCALAside(tree_reader, "FCAL_side_A_from_layers");
    
    TTreeReaderArray<int> Nvertex(tree_reader, "nvertex");
    TTreeReaderArray<int> RunNumber(tree_reader, "run_number");

    TH1D *ZDC_ASide_hist =new TH1D("ZDC_ASide_hist","ZDC_ASide",100,0,0.01);
    TH1D *ZDC_CSide_hist =new TH1D("ZDC_CSide_hist","ZDC_CSide",100,0,0.01);


    TH1D *FCAL_hist =new TH1D("FCAL_hist","FCAL",100,0,0.6);
    TH1D *FCAL_ZeroZDC_hist =new TH1D("FCAL_ZeroZDC_hist","FCAL_ZeroZDC",100,0,0.6);
    TH1D *FCAL_Bining =new TH1D("FCAL_Bining","FCAL_Bining",10,-0.5,9.5);
    for (int i = 0; i < 10; i++) {
        FCAL_Bining->GetXaxis()->SetBinLabel(i + 1, Form("%d-%d %%", i, i + 1));
    }



    Long64_t total_events = mychain->GetEntries();
    cout<<"============================="<<endl;
    cout<<"Starting Analysis"<<endl;
    int eventID=0;
    double Centrality[]={299.736, 277.788, 262.923, 251.217, 241.358, 232.613, 224.684, 217.391, 210.599, 204.217, 198.159, 192.374, 186.851, 181.527, 176.384, 171.405, 166.613, 161.958, 157.427, 153.033};
    double Nparticipates[]={27.88, 26.96, 26.37, 25.86, 25.39, 24.93, 24.49, 24.04, 23.60, 23.15, 21.80, 19.52};
    double Scaling[]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,};
    double SumParticipatesAllEvents=0,EventPassed=0;
    int MaxBin=size(Centrality);
    for(int bin=0; bin<size(Centrality);bin++) if(Centrality[bin]<FcalBoundry) { MaxBin=bin; break;}
    while(tree_reader.Next()){
        if (eventID % 10000 == 0) {
            std::cout << Form("\rDone: %.2f", (100.0 * eventID / total_events)) <<"%   "<< std::flush;
        }
        eventID++;
        //if(Nvertex[0]!=1) continue;
        ZDC_CSide_hist->Fill(ZDC_CSide[0]);
        ZDC_ASide_hist->Fill(ZDC_ASide[0]);

        FCAL_hist->Fill(FCAL[0]);

        if(ZDC_CSide[0]<0.001 && ZDC_ASide[0]<0.001 )
        {
            FCAL_ZeroZDC_hist->Fill(FCAL[0]);
            if(FCAL[0]>FcalBoundry)

            for(int bin=0; bin<MaxBin;bin++)
            {
                if(FCAL[0]*1000>Centrality[bin]) 
                {
                    FCAL_Bining->Fill(bin,1);
                    if(bin<10) Scaling[bin]++;
                    else if(bin<15) Scaling[10]++;
                    else if(bin<20) Scaling[11]++;
                    EventPassed++;
                    break;
                }
            }
        }
    }
    cout<<endl;

    double Average=0;
    for(int bin=0;bin<size(Scaling);bin++)
    {
        Scaling[bin]/=EventPassed;
        cout<<Scaling[bin]*100<<"%"<<endl;
        Average+=Scaling[bin]*Nparticipates[bin];
    }
    
    cout<<endl;
    cout<<"End of Analysis"<<endl;
    cout<<"Integral Ratio: "<<FCAL_ZeroZDC_hist->Integral(29,70)/FCAL_hist->Integral(29,70)<<endl;
    cout<<"Average Participates: "<<Average<<endl;

    ZOutput->cd();
    ZDC_CSide_hist->Write();
    ZDC_ASide_hist->Write();


    FCAL_hist->Write();
    FCAL_ZeroZDC_hist->Write();
    FCAL_Bining->Write();


    ZOutput->Close();
    cout<<"============================="<<endl;

}