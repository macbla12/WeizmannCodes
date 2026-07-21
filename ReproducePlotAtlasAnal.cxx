void ReproducePlotAtlasAnal()
{

    
    TString File ="Files/dataOOwithTracks.root"; 
    double Centrality[]={299.736, 277.788, 262.923, 251.217, 241.358, 204.217,176.384,153.033, 132.738,  114.885, 99.082, 85.128, 72.841, 62.061, 52.594,  44.348, 37.123, 30.774, 25.169, 20.215};

    double Nparticipates[]={27.88, 26.96, 26.37, 25.86, 25.39, 24.04, 21.80, 19.52, 17.29, 15.22, 13.34, 11.64, 10.11, 8.75, 7.54, 6.47, 5.52, 4.72, 4.04, 3.50};
    double TAA[]={0.6750, 0.6255, 0.5981, 0.5761, 0.5567, 0.5066, 0.4314, 0.3633, 0.3030, 0.2516, 0.2082, 0.1719, 0.1415, 0.1162, 0.0953, 0.0780, 0.0636, 0.0518, 0.0424, 0.0350};

    double Extrapolation[]={1.364, 1.358, 1.359, 1.367, 1.370, 1.361, 1.379, 1.375, 1.400, 1.378, 1.381, 1.400, 1.398, 1.388, 1.403, 1.469, 1.536, 1.572, 1.611, 1.519};

    TFile *file = TFile::Open("Files/OO_eff005.root");
    TH1F *h_eff = (TH1F*)file->Get("OO_eff005"); 


    TChain *mychain = new TChain("Z_analysis");
        mychain->Add(File);
    TTreeReader tree_reader(mychain);

    TFile *Output = new TFile("Plots/ReproducedPlot.root","RECREATE");

    TTreeReaderArray<int> NTracks(tree_reader, "ntracks");
    TTreeReaderArray<float> TrackEta(tree_reader, "trk_eta");
    TTreeReaderArray<float> TrackPt(tree_reader, "trk_pt");
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

    TF1* pt_dep_d0_func = new TF1("pt_dep_d0", "0.14+ 0.24/x+(1.5e-04)*x", 0, 500); 

    double trkPtBins_wguard[] = {0.1, 0.2, 0.31, 0.407, 0.5, 0.575, 0.66, 0.755, 0.86, 0.999, 1.13, 1.3, 1.5, 1.7, 1.95, 2.25, 2.57, 2.95, 3.35, 3.87, 4.45, 5.1, 5.85, 6.7, 
                                7.65, 8.8, 10, 11.5, 13.2, 15.1, 17.3, 20, 23, 26, 30, 37.5, 47.5, 60, 75, 95, 120, 150, 190, 240, 290, 350};
    TH1D *TracksPtBeforeEff =new TH1D("TracksPtBeforeEff","TracksPtBeforeEff", 45, trkPtBins_wguard);
    TH1D *TracksPtBeforeEff_ZDCcut =new TH1D("TracksPtBeforeEff_ZDCcut","TracksPtBeforeEff_ZDCcut", 45, trkPtBins_wguard);


    TH1D *FCAL_hist =new TH1D("FCAL_hist","FCAL",100,0,0.6);
    TH1D *FCAL_Bining =new TH1D("FCAL_Bining","FCAL_Bining",12,-0.5,11.5);
    for (int i = 0; i < 12; i++) {
        if(i<5) FCAL_Bining->GetXaxis()->SetBinLabel(i + 1, Form("%d-%d %%", i, i + 1));
        else FCAL_Bining->GetXaxis()->SetBinLabel(i - 4, Form("%d-%d %%", i*5, (i + 1)*5));
    }

    TH1D *FCAL0n0noverFCALmb =new TH1D("FCAL0n0noverFCALmb","FCAL0n0noverFCALmb",100,0,0.6);
    TH1D *FCAL1n0noverFCALmb =new TH1D("FCAL1n0noverFCALmb","FCAL1n0noverFCALmb",100,0,0.6);
    TH1D *FCAL1n1noverFCALmb =new TH1D("FCAL1n1noverFCALmb","FCAL1n1noverFCALmb",100,0,0.6);

    TH1D *ZDC_hist =new TH1D("ZDC_hist","ZDC_ASide",100,0,0.01);
    TH1D *ZDC_centrality1pr =new TH1D("ZDC_centrality1pr","ZDC_centrality1pr",100,0,0.01);
    TH1D *ZDC_centrality10pr =new TH1D("ZDC_centrality10pr","ZDC_centrality10pr",100,0,0.01);




    Long64_t total_events = mychain->GetEntries();
    cout<<"============================="<<endl;
    cout<<"Starting Analysis"<<endl;
    int eventID=0;
    double NumberOfTracksBins[20] = {0.0};
    double NumberOfEventsBins[20] = {0.0};
    double  Counts10procent=0,CountsZDC=0;

    
    while(tree_reader.Next()){
        //if(eventID>1000) break;
        if (eventID % 100 == 0) {
            cout << Form("\rDone: %.2f", (100.0 * eventID / total_events)) <<"%   "<< flush;
        }
        eventID++;
        //Vertex Cut
        vector<pair<double, int>> vertexOrder;
        vertexOrder.reserve(Nvertex[0]);
        if(Nvertex[0]==0) continue;
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
        if(FCAL[0]<0.020983) continue;
        
        FCAL_hist->Fill(FCAL[0]);
        if(FCAL[0]*1000>204.217) Counts10procent++;
        if(FCAL[0]*1000>160){ if(ZDC_CSide[0]<0.001 && ZDC_ASide[0]<0.001) CountsZDC++;}


        double NumberOfTracks=0;
        for(int track=0; track<NTracks[0]; track++){   
            double pT=abs(TrackPt[track]);
            int bin = h_eff->FindBin(pT);
            double eff = h_eff->GetBinContent(bin);
            if(abs(TrackEta[track])<0.5 && pT>0.31 && TrackQuality[track][1]&& std::abs(TrackD0[track]) < pt_dep_d0_func->Eval(pT) && eff>0)
            {
                NumberOfTracks+= 1.0/eff;
                if(FCAL[0]*1000>204.217) TracksPtBeforeEff->Fill(pT);
                if(FCAL[0]*1000>160){ if(ZDC_CSide[0]<0.001 && ZDC_ASide[0]<0.001) TracksPtBeforeEff_ZDCcut->Fill(pT);}
            }
        } 
        ZDC_hist->Fill(ZDC_CSide[0]);
        for(int bin=0; bin<20;bin++)
        {
            if(FCAL[0]*1000>Centrality[bin]) 
            {
                if(bin==0) ZDC_centrality1pr->Fill(ZDC_CSide[0]);
                if(bin<6) ZDC_centrality10pr->Fill(ZDC_CSide[0]);

                FCAL_Bining->Fill(bin,1);
                NumberOfTracksBins[bin]+=NumberOfTracks;
                NumberOfEventsBins[bin]++;
                break;
            }
        }

        if(ZDC_CSide[0]<0.001 && ZDC_ASide[0]<0.001)
        {
            FCAL0n0noverFCALmb->Fill(FCAL[0]);
        }
        if(((ZDC_CSide[0]>0.001 && ZDC_CSide[0]<0.004) && ZDC_ASide[0]<0.001) || ((ZDC_ASide[0]>0.001 && ZDC_ASide[0]<0.004) && ZDC_CSide[0]<0.001))
        {
            FCAL1n0noverFCALmb->Fill(FCAL[0]);
        }
        if((ZDC_CSide[0]>0.001 && ZDC_CSide[0]<0.004) && (ZDC_ASide[0]>0.001 && ZDC_ASide[0]<0.004))
        {
            FCAL1n1noverFCALmb->Fill(FCAL[0]);
        }
        

    }
    cout<<endl;
    cout<<"==================================================="<<endl;
    cout<<"End of Analysis"<<endl;

    double TrackEfficiency=1;
    //double EtaCut= 0.21; //abs(eta)<0.5
    double FinalNumberOfTracksBin[20]={0.0};
    double FinalTAANumberOfTracksBin[20]={0.0};
    for(int bin=0; bin<20;bin++)
    {
        FinalNumberOfTracksBin[bin]= NumberOfTracksBins[bin]*(Extrapolation[bin]*1.083)/(NumberOfEventsBins[bin]*Nparticipates[bin]);

        FinalTAANumberOfTracksBin[bin]= NumberOfTracksBins[bin]*(Extrapolation[bin]*1.083)/(NumberOfEventsBins[bin]*TAA[bin]);
        cout<<bin+1<<" Bin: N particles "<<Nparticipates[bin]<<" 1/<Npart>dn/deta="<<FinalNumberOfTracksBin[bin]<<endl;


    }


    double sumTracks_0_5 = 0, sumEvents_0_5 = 0;
    double sumNpart_weighted_0_5 = 0, sumExtrap_weighted_0_5 = 0;

    for(int bin=0; bin<5; bin++)
    {
        sumTracks_0_5 += NumberOfTracksBins[bin];
        sumEvents_0_5 += NumberOfEventsBins[bin];
        sumNpart_weighted_0_5  += Nparticipates[bin]  * NumberOfEventsBins[bin];
        sumExtrap_weighted_0_5 += Extrapolation[bin]  * NumberOfEventsBins[bin];
    }

    double Npart_0_5  = sumNpart_weighted_0_5  / sumEvents_0_5;
    double Extrap_0_5 = sumExtrap_weighted_0_5 / sumEvents_0_5;

    double FinalValue_0_5 = sumTracks_0_5 * (Extrap_0_5 * 1.083) / (sumEvents_0_5 * Npart_0_5);
    cout << "Extrapolation in 0-5% ="<<Extrap_0_5<<endl;

    cout << "0-5%: Npart=" << Npart_0_5 << " 1/<Npart>dn/deta=" << FinalValue_0_5 <<"  dn/deta="<<FinalValue_0_5*Npart_0_5<<endl;

    double sumTracks_0_10 = 0, sumEvents_0_10 = 0;
    double sumNpart_weighted_0_10 = 0, sumExtrap_weighted_0_10 = 0;

    for(int bin=0; bin<6; bin++)
    {
        sumTracks_0_10 += NumberOfTracksBins[bin];
        sumEvents_0_10 += NumberOfEventsBins[bin];
        sumNpart_weighted_0_10  += Nparticipates[bin]  * NumberOfEventsBins[bin];
        sumExtrap_weighted_0_10 += Extrapolation[bin]  * NumberOfEventsBins[bin];
    }

    double Npart_0_10  = sumNpart_weighted_0_10  / sumEvents_0_10;
    double Extrap_0_10 = sumExtrap_weighted_0_10 / sumEvents_0_10;

    double FinalValue_0_10 = sumTracks_0_10 * (Extrap_0_10 * 1.083) / (sumEvents_0_10 * Npart_0_10);
    cout << "Extrapolation in 0-10% ="<<Extrap_0_10* 1.083<<endl;

    cout << "0-10%: Npart=" << Npart_0_10 << " 1/<Npart>dn/deta=" << FinalValue_0_10 <<"  dn/deta="<<FinalValue_0_10*Npart_0_10<<endl;

    double sumTracks_0_80 = 0, sumEvents_0_80 = 0;
    double sumNpart_weighted_0_80 = 0, sumExtrap_weighted_0_80 = 0;

    for(int bin=0; bin<20; bin++)
    {
        sumTracks_0_80 += NumberOfTracksBins[bin];
        sumEvents_0_80 += NumberOfEventsBins[bin];
        sumNpart_weighted_0_80  += Nparticipates[bin]  * NumberOfEventsBins[bin];
        sumExtrap_weighted_0_80 += Extrapolation[bin]  * NumberOfEventsBins[bin];
    }

    double Npart_0_80  = sumNpart_weighted_0_80  / sumEvents_0_80;
    double Extrap_0_80 = sumExtrap_weighted_0_80 / sumEvents_0_80;

    double FinalValue_0_80 = sumTracks_0_80 * (Extrap_0_80 * 1.083) / (sumEvents_0_80 * Npart_0_80);
    cout << "Extrapolation in 0-80% ="<<Extrap_0_80* 1.083<<endl;

    cout << "0-80%: Npart=" << Npart_0_80 << " 1/<Npart>dn/deta=" << FinalValue_0_80 <<"  dn/deta="<<FinalValue_0_80*Npart_0_80<<endl;

    cout<<"Ratio: "<<FinalValue_0_10/FinalValue_0_80<<endl;



    TH1D *TracksPtAfterEff = (TH1D*)TracksPtBeforeEff->Clone("TracksPtAfterEff");
    TracksPtAfterEff->SetTitle("TracksPtAfterEff"); 

    TH1D *TracksPtAfterEff_ZDCcut = (TH1D*)TracksPtBeforeEff_ZDCcut->Clone("TracksPtAfterEff_ZDCcut");
    TracksPtAfterEff_ZDCcut->SetTitle("TracksPtAfterEff_ZDCcut");

    TracksPtAfterEff->Divide(h_eff);
    TracksPtAfterEff_ZDCcut->Divide(h_eff);
    //cout<<Counts10procent<<endl;
    TracksPtBeforeEff->Scale(1/Counts10procent);
    TracksPtAfterEff->Scale(1/Counts10procent);
    TracksPtBeforeEff_ZDCcut->Scale(1/CountsZDC);
    TracksPtAfterEff_ZDCcut->Scale(1/CountsZDC);



    TGraph *graph = new TGraph(20, Nparticipates, FinalNumberOfTracksBin);
    graph->SetName("g_Tracks_vs_Npart");


    TGraph *graphTAA = new TGraph(20, TAA, FinalTAANumberOfTracksBin);
    graphTAA->SetName("g_Tracks_vs_TAA");

    Output->cd();
    
    graph->Write(); 
    graphTAA->Write(); 
    
    TracksPtBeforeEff->Write();
    TracksPtAfterEff->Write();
    TracksPtBeforeEff_ZDCcut->Write();
    TracksPtAfterEff_ZDCcut->Write();
    
    FCAL_hist->Write();
    FCAL0n0noverFCALmb->Write();
    FCAL1n0noverFCALmb->Write();
    FCAL1n1noverFCALmb->Write();

    ZDC_hist->Write();
    ZDC_centrality1pr->Write();
    ZDC_centrality10pr->Write();

    Output->Close();


    cout<<"============================="<<endl;

}