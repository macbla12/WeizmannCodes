void ReproducePlotAnal()
{

    double FcalBoundary=0.200;
    
    TString File ="Files/dataOOwithTracks.root"; 
    double Centrality[]={299.736, 277.788, 262.923, 251.217, 241.358, 204.217,176.384,153.033, 132.738,  114.885, 99.082, 85.128, 72.841, 62.061, 52.594,  44.348, 37.123, 30.774, 25.169, 20.215};

    double Nparticipates[]={27.88, 26.96, 26.37, 25.86, 25.39, 24.04, 21.80, 19.52, 17.29, 15.22, 13.34, 11.64, 10.11, 8.75, 7.54, 6.47, 5.52, 4.72, 4.04, 3.50};
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

    Long64_t total_events = mychain->GetEntries();
    cout<<"============================="<<endl;
    cout<<"Starting Analysis"<<endl;
    int eventID=0;
    double NumberOfTracksBins[20] = {0.0};
    double NumberOfEventsBins[20] = {0.0};
    double  Counts10procent=0,CountsZDC=0;

    
    while(tree_reader.Next()){
        //if(eventID>1000) break;
        if (eventID % 1000 == 0) {
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

        for(int bin=0; bin<20;bin++)
        {
            if(FCAL[0]*1000>Centrality[bin]) 
            {
                FCAL_Bining->Fill(bin,1);
                NumberOfTracksBins[bin]+=NumberOfTracks;
                NumberOfEventsBins[bin]++;
                break;
            }
        }
        
        

    }
    cout<<endl;
    cout<<"==================================================="<<endl;
    cout<<"End of Analysis"<<endl;

    double TrackEfficiency=1;
    //double EtaCut= 0.21; //abs(eta)<0.5
    double FinalNumberOfTracksBin[20]={0.0};
    for(int bin=0; bin<20;bin++)
    {
        FinalNumberOfTracksBin[bin]= NumberOfTracksBins[bin]*(Extrapolation[bin]*1.09)/(NumberOfEventsBins[bin]*Nparticipates[bin]);
        //cout<<bin+1<<" Bin: 1/<Npart>dn/deta="<<FinalNumberOfTracksBin[bin]<<endl;
    }

    TH1D *TracksPtAfterEff = (TH1D*)TracksPtBeforeEff->Clone("TracksPtAfterEff");
    TracksPtAfterEff->SetTitle("TracksPtAfterEff"); 

    TH1D *TracksPtAfterEff_ZDCcut = (TH1D*)TracksPtBeforeEff_ZDCcut->Clone("TracksPtAfterEff_ZDCcut");
    TracksPtAfterEff_ZDCcut->SetTitle("TracksPtAfterEff_ZDCcut");

    TracksPtAfterEff->Divide(h_eff);
    TracksPtAfterEff_ZDCcut->Divide(h_eff);
    cout<<Counts10procent<<endl;
    TracksPtBeforeEff->Scale(1/Counts10procent);
    TracksPtAfterEff->Scale(1/Counts10procent);
    TracksPtBeforeEff_ZDCcut->Scale(1/CountsZDC);
    TracksPtAfterEff_ZDCcut->Scale(1/CountsZDC);

    TGraph *graph = new TGraph(20, Nparticipates, FinalNumberOfTracksBin);
    
    graph->SetMarkerStyle(20);        
    graph->SetMarkerSize(1.2);       
    graph->SetMarkerColor(kRed);  
    graph->SetLineColor(kRed);    
    graph->SetName("g_Tracks_vs_Npart");
    graph->SetTitle(";<N_{part}>;<1/<N_{part}> dn/d#eta>");
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    graph->GetYaxis()->SetRangeUser(0, 8);
    graph->GetXaxis()->SetLimits(0.0, 36.0);

    graph->Draw("APL"); 

    c1->SaveAs("Plots/Tracks_vs_Npart.pdf");
    Output->cd();
    graph->Write(); 
    TracksPtBeforeEff->Write();
    TracksPtAfterEff->Write();
    TracksPtBeforeEff_ZDCcut->Write();
    TracksPtAfterEff_ZDCcut->Write();
    Output->Close();
    cout<<"============================="<<endl;

}