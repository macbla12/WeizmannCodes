void Anal()
{

    double FcalBoundry=0.175;
    
    TString File ="Files/dataOO.root"; 
    double Centrality[]={299.736, 277.788, 262.923, 251.217, 241.358, 232.613, 224.684, 217.391, 210.599, 204.217, 198.159, 192.374, 186.851, 181.527, 176.384, 171.405, 166.613, 161.958, 157.427, 153.033};
    double Nparticipates[]={27.88, 26.96, 26.37, 25.86, 25.39, 24.93, 24.49, 24.04, 23.60, 23.15, 21.80, 19.52};
    double NDeltaparticipates[]={0.23368, 0.31184, 0.35184, 0.38368, 0.41184, 0.45184, 0.47184, 0.51000, 0.53184, 0.57000, 0.63184,0.71000};

    double Ncollisions[]={45.90, 42.54, 40.67, 39.17, 37.86, 36.65, 35.51, 34.41, 33.35, 32.31, 29.34,24.70};
    double NDeltacollisions[]={4.00736, 3.48704, 3.61552, 3.56184, 3.49368, 3.43552, 3.38184, 3.30736, 3.21656, 3.13576, 2.89968, 2.70496};

    double TAA[]={0.6750, 0.6255, 0.5981, 0.5761, 0.5567, 0.5390, 0.5222, 0.5061, 0.4904, 0.4752, 0.4314, 0.3633};
    double DeltaTAA[]={0.05891, 0.05143, 0.05317, 0.05232, 0.05144, 0.05056, 0.04974, 0.04853, 0.04736, 0.04609, 0.04264,0.03963};
    
    
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

    TTreeReaderArray<int> NLoose(tree_reader, "n_loose");

    TTreeReaderArray<float> FCAL(tree_reader, "FCAL");
    TTreeReaderArray<float> FCALCside(tree_reader, "FCAL_side_C_from_layers");
    TTreeReaderArray<float> FCALAside(tree_reader, "FCAL_side_A_from_layers");
    
    TTreeReaderArray<int> Nvertex(tree_reader, "nvertex");
    TTreeReaderArray<float> ZVertex(tree_reader, "vertex_z");
    TTreeReaderArray<float> ZVertexVar(tree_reader, "vertex_z_var");


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


    double Scaling[]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,};


    int eventID=0;
    double SumParticipatesAllEvents=0,EventPassed=0;
    double HILooseTrackAll=0, HILooseTrackSelection=0;
    int MaxBin=size(Centrality);
    for(int bin=0; bin<size(Centrality);bin++) if(Centrality[bin]<FcalBoundry) { MaxBin=bin; break;}
    while(tree_reader.Next()){
        //if(eventID>1000) break;
        if (eventID % 10000 == 0) {
            std::cout << Form("\rDone: %.2f", (100.0 * eventID / total_events)) <<"%   "<< std::flush;
        }
        eventID++;
        //Vertex Cut
        int NumberOfVertex=0;
        for(int vertex=0; vertex<Nvertex[0]; vertex++)
        {
            //cout<<ZVertexVar[vertex]<<" Variacja"<<endl;
            if(abs(ZVertex[vertex])>150) continue;
            if(ZVertexVar[vertex]<0.02) NumberOfVertex++;
        }
        //cout<<NumberOfVertex<<" Vertex"<<endl;
        if(NumberOfVertex!=1) continue;
        
        
        if(FCAL[0]<0.020983) continue;
        FCAL_hist->Fill(FCAL[0]);
        HILooseTrackAll+=NLoose[0];
        if(ZDC_CSide[0]<0.001 && ZDC_ASide[0]<0.001)
        {
            ZDC_CSide_hist->Fill(ZDC_CSide[0]);
            ZDC_ASide_hist->Fill(ZDC_ASide[0]);
            FCAL_ZeroZDC_hist->Fill(FCAL[0]);
            if(FCAL[0]>FcalBoundry)
            {
                HILooseTrackSelection+=NLoose[0];

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

    }
    cout<<endl;
    cout<<EventPassed<<"==================================================="<<endl;
    double AverageParticipates=0, AverageParticipatesDeltaSq=0,AverageCollisions=0, AverageCollisionsDeltaSq=0,AverageTAA=0, AverageTAADeltaSq=0;
    for(int bin=0;bin<size(Scaling);bin++)
    {
        Scaling[bin]/=EventPassed;
        cout<<Scaling[bin]*100<<"%"<<endl;
        AverageParticipates+=Scaling[bin]*Nparticipates[bin];
        AverageParticipatesDeltaSq += (Scaling[bin] * Scaling[bin]) * (NDeltaparticipates[bin] * NDeltaparticipates[bin]);

        AverageCollisions+=Scaling[bin]*Ncollisions[bin];
        AverageCollisionsDeltaSq += (Scaling[bin] * Scaling[bin]) * (NDeltacollisions[bin] * NDeltacollisions[bin]);

        AverageTAA+=Scaling[bin]*TAA[bin];
        AverageTAADeltaSq += (Scaling[bin] * Scaling[bin]) * (DeltaTAA[bin] * DeltaTAA[bin]);
    }
    double AverageParticpatesDelta = sqrt(AverageParticipatesDeltaSq);
    double AverageCollisionsDelta = sqrt(AverageCollisionsDeltaSq);
    double AverageTAADelta= sqrt(AverageTAADeltaSq);

    cout<<endl;
    cout<<"End of Analysis"<<endl;
    cout<<"Integral Ratio: "<<FCAL_ZeroZDC_hist->Integral(FCAL_ZeroZDC_hist->FindBin(FcalBoundry),70)/FCAL_hist->Integral(FCAL_hist->FindBin(FcalBoundry),70)<<endl;
    cout<<"Average Participates: "<<AverageParticipates<< " +/- " << AverageParticpatesDelta <<endl;
    cout<<"Average Collisions: "<<AverageCollisions<< " +/- " << AverageCollisionsDelta <<endl;
    cout<<"Average TAA: "<<AverageTAA<< " +/- " << AverageTAADelta <<endl;

    double TrackEfficiency=0.5; //=0.8;//
    double EtaCut= 0.21; //abs(eta)<0.5
    double FullSpectrumNparicipates=10.5;

    double All_mean_dn_deta = (HILooseTrackAll*EtaCut) / (eventID*TrackEfficiency); 
    double final_value_All = All_mean_dn_deta / FullSpectrumNparicipates;

    double Selection_mean_dn_deta = (HILooseTrackSelection*EtaCut) / (EventPassed*TrackEfficiency); 
    double final_value_Seletion = Selection_mean_dn_deta / AverageParticipates;

    cout<<"Full spectrum N tracks: "<< final_value_All<<" with "<<FullSpectrumNparicipates<<endl;
    cout<<"Selection spectrum N tracks: "<<final_value_Seletion <<" with "<<AverageParticipates<<endl;

    cout<<"Full spectrum N tracks: "<<HILooseTrackAll <<" Selected spectrum N tracks: "<<HILooseTrackSelection<<endl;




    ZOutput->cd();
    ZDC_CSide_hist->Write();
    ZDC_ASide_hist->Write();


    FCAL_hist->Write();
    FCAL_ZeroZDC_hist->Write();
    FCAL_Bining->Write();


    ZOutput->Close();
    cout<<"============================="<<endl;

}