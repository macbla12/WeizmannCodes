void AnalOO()
{

    double FcalBoundary=0.160;
    
    TString File ="Files/dataOOwithTracks.root"; 
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

    TFile *Output = new TFile("Plots/FCaL-Centrality.root","RECREATE");


    TTreeReaderArray<int> NTracks(tree_reader, "ntracks");
    TTreeReaderArray<float> TrackEta(tree_reader, "trk_eta");
    TTreeReaderArray<vector<int>> TrackQuality(tree_reader, "trk_qualities");
    TTreeReaderArray<float> TrackD0(tree_reader, "trk_d0");
    TTreeReaderArray<float> TrackZ0(tree_reader, "trk_z0");
    TTreeReaderArray<float> TrackTheta(tree_reader, "trk_theta");



    TTreeReaderArray<float> ZDC_CSide(tree_reader, "ZDC_energy_side_C");
    TTreeReaderArray<float> ZDC_ASide(tree_reader, "ZDC_energy_side_A");

    TTreeReaderArray<int> NLoose(tree_reader, "n_loose");

    TTreeReaderArray<float> FCAL(tree_reader, "FCAL");
    TTreeReaderArray<float> FCALCside(tree_reader, "FCAL_side_C_from_layers");
    TTreeReaderArray<float> FCALAside(tree_reader, "FCAL_side_A_from_layers");
    
    //TTreeReaderArray<float> Trigger(tree_reader, "trigger_prescale");

    TTreeReaderArray<int> Nvertex(tree_reader, "nvertex");
    TTreeReaderArray<float> ZVertex(tree_reader, "vertex_z");
    TTreeReaderArray<float> ZVertexVar(tree_reader, "vertex_z_var");

    TTreeReaderArray<bool> MB_topo_flag(tree_reader, "MB_topo_flag");
    TTreeReaderArray<int> RunNumber(tree_reader, "run_number");
    TTreeReaderArray<int> LumiBlock(tree_reader, "lb");

    TH1D *track_theta_hist = new TH1D("track_theta_hist","track_theta_hist",100,-2.5,2.5);
    TH1D *track_pT_hist = new TH1D("track_pT_hist","track_pT_hist",100,-2.5,2.5);


    TH1D *ZDC_ASide_hist =new TH1D("ZDC_ASide_hist","ZDC_ASide",100,0,0.01);
    TH1D *ZDC_CSide_hist =new TH1D("ZDC_CSide_hist","ZDC_CSide",100,0,0.01);

    TH1D *Tracks_withPileUp_hist =new TH1D("Tracks_withPileUp_hist","Tracks_withPileUp_hist",250,0,250);
    TH1D *Tracks_fullspectrum_hist =new TH1D("Tracks_fullspectrum_hist","Tracks_fullspectrum_hist",250,0,250);
    TH1D *Tracks_selectedspectrum_hist =new TH1D("Tracks_selectedspectrum_hist","Tracks_selectedspectrum_hist",250,0,250);

    TH2D *Tracks_withPileUp_histvsFCAL =new TH2D("Tracks_withPileUp_histvsFCAL","Tracks_withPileUp_histvsFCAL",250,0,250,100,0,0.6);

    TH1D *FCAL_hist =new TH1D("FCAL_hist","FCAL",100,0,0.6);
    TH1D *FCAL_ZeroZDC_hist =new TH1D("FCAL_ZeroZDC_hist","FCAL_ZeroZDC",100,0,0.6);
    TH1D *FCAL_Bining =new TH1D("FCAL_Bining","FCAL_Bining",12,-0.5,11.5);
    for (int i = 0; i < 12; i++) {
        if(i<10) FCAL_Bining->GetXaxis()->SetBinLabel(i + 1, Form("%d-%d %%", i, i + 1));
        if(i==10) FCAL_Bining->GetXaxis()->SetBinLabel(i + 1, "10-15%");
        if(i==11) FCAL_Bining->GetXaxis()->SetBinLabel(i + 1, "15-20%");
    }

    Long64_t total_events = mychain->GetEntries();
    cout<<"============================="<<endl;
    cout<<"Starting Analysis"<<endl;
    double Scaling[]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,};
    int eventID=0;
    double SumParticipatesAllEvents=0,EventPassed=0,CountFullsepctrum=0;;
    double HILooseTrackAll=0, HILooseTrackSelection=0;
    int MaxBin=size(Centrality);
    for(int bin=0; bin<size(Centrality);bin++) if(Centrality[bin]<FcalBoundary) { MaxBin=bin; break;}
    while(tree_reader.Next()){
        //if(eventID>1000) break;
        if (eventID % 1000 == 0) {
            cout << Form("\rDone: %.2f", (100.0 * eventID / total_events)) <<"%   "<< flush;
        }
        eventID++;
        //cout<<Trigger[0]<<endl;
        //if(Trigger[0]<=0) continue;
        int NumberOfTracks=0;
        for(int track=0; track<NTracks[0]; track++) if(abs(TrackEta[track])<0.5 && TrackQuality[track][0]) NumberOfTracks++; // && abs(TrackD0[track])<1.5 && abs((TrackZ0[track])*sin(TrackTheta[track]))<1.5
        //if(eventID<100) cout<<NumberOfTracks<<endl;
        //if(FCAL[0] > (0.005 * NumberOfTracks + 0.1)) continue;
        //if(NumberOfTracks<5) continue;

        Tracks_withPileUp_hist->Fill(NumberOfTracks);
        Tracks_withPileUp_histvsFCAL->Fill(NumberOfTracks,FCAL[0]);

        
        //Vertex Cut
        vector<pair<double, int>> vertexOrder;
        vertexOrder.reserve(Nvertex[0]);
        if(Nvertex[0]==0) continue;
        //if(Nvertex[0]>2) continue;

        if(abs(ZVertex[0]) > 150) continue; 

        for (int vertex = 0; vertex < Nvertex[0]; vertex++)
        {
            vertexOrder.push_back({ZVertexVar[vertex], vertex});
        }
        sort(vertexOrder.begin(), vertexOrder.end());
        int NumberOfVertex = 1;
        for (size_t i = 1; i < vertexOrder.size(); i++) if(vertexOrder[i].first < 0.02) NumberOfVertex++;

        if (NumberOfVertex != 1) continue;
        
        if(!MB_topo_flag[0]) continue;

        //if(LumiBlock[0]<1100) continue;
        
        if(FCAL[0]<0.020215) continue;


        FCAL_hist->Fill(FCAL[0]);
        double extrapolation=0;
        if(FCAL[0]<37.123) extrapolation=1.55;
        else extrapolation=1.4;
        CountFullsepctrum++;
        HILooseTrackAll+=NumberOfTracks*extrapolation;
        Tracks_fullspectrum_hist->Fill(NumberOfTracks);

        ZDC_CSide_hist->Fill(ZDC_CSide[0]);
        ZDC_ASide_hist->Fill(ZDC_ASide[0]);


        if(ZDC_CSide[0]<0.001 && ZDC_ASide[0]<0.001)
        {
            
            FCAL_ZeroZDC_hist->Fill(FCAL[0]);

            if(FCAL[0]>FcalBoundary)
            {
                HILooseTrackSelection+=NumberOfTracks*1.4;  
                Tracks_selectedspectrum_hist->Fill(NumberOfTracks);


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
    cout<<"==================================================="<<endl;
    double AverageParticipates=0, AverageParticipatesDeltaSq=0,AverageCollisions=0, AverageCollisionsDeltaSq=0,AverageTAA=0, AverageTAADeltaSq=0;
    for(int bin=0;bin<size(Scaling);bin++)
    {
        Scaling[bin]/=EventPassed;
        cout<<bin<<" Bin with "<<Scaling[bin]*100<<"%"<<endl;
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
    cout<<"Average Participates: "<<AverageParticipates<< " +/- " << AverageParticpatesDelta <<endl;
    cout<<"Average Collisions: "<<AverageCollisions<< " +/- " << AverageCollisionsDelta <<endl;
    cout<<"Average TAA: "<<AverageTAA<< " +/- " << AverageTAADelta <<endl;

    double TrackEfficiency=1;
    //double EtaCut= 0.21; //abs(eta)<0.5
    double FullSpectrumNparicipates=12.5;

    double All_mean_dn_deta = (HILooseTrackAll) / (CountFullsepctrum*TrackEfficiency); 
    double final_value_All = All_mean_dn_deta / FullSpectrumNparicipates;

    double Selection_mean_dn_deta = (HILooseTrackSelection) / (EventPassed*TrackEfficiency); 
    double final_value_Seletion = Selection_mean_dn_deta / (AverageParticipates);

    cout<<"Full spectrum N tracks: "<< final_value_All<<" with "<<FullSpectrumNparicipates<<endl;
    cout<<"Selection spectrum N tracks: "<<final_value_Seletion <<" with "<<AverageParticipates<<endl;
    cout<<"Ratio of tracks: "<<final_value_Seletion/final_value_All<<endl;

    // cout<<"Full spectrum N tracks: "<<HILooseTrackAll <<" Selected spectrum N tracks: "<<HILooseTrackSelection<<endl;
    Output->cd();
    ZDC_CSide_hist->Write();
    ZDC_ASide_hist->Write();
    FCAL_hist->Write();
    FCAL_ZeroZDC_hist->Write();
    FCAL_Bining->Write();
    Tracks_withPileUp_hist->Write();
    Tracks_fullspectrum_hist->Write();
    Tracks_selectedspectrum_hist->Write();
    Tracks_withPileUp_histvsFCAL->Write();
    Output->Close();
    cout<<"============================="<<endl;

}