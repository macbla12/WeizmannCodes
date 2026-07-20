double SampleGamma(double k, double theta)
{
    if (k < 1.0) {
        double u = gRandom->Uniform();
        return SampleGamma(k + 1.0, theta) * TMath::Power(u, 1.0/k);
    }
    double d = k - 1.0/3.0;
    double c = 1.0/TMath::Sqrt(9.0*d);
    double x, v, u;
    while (true) {
        do {
            x = gRandom->Gaus(0,1);
            v = 1.0 + c*x;
        } while (v <= 0);
        v = v*v*v;
        u = gRandom->Uniform();
        if (u < 1.0 - 0.0331*x*x*x*x) return d*v*theta;
        if (TMath::Log(u) < 0.5*x*x + d*(1.0 - v + TMath::Log(v))) return d*v*theta;
    }
}



void FCALtest()
{

    TString InFile="Files/glauber_O16.root";

    TChain *mychain = new TChain("nt_O_O");
        mychain->Add(InFile);
    TTreeReader tree_reader(mychain);
    TFile *Output = new TFile("Plots/FCALtest.root","RECREATE");

    TTreeReaderArray<float> Npart(tree_reader, "Npart");
    TTreeReaderArray<float> Ncoll(tree_reader, "Ncoll");
    TTreeReaderArray<float> ImpPara(tree_reader, "B");


    TH1D* NcollHist = new TH1D("Ncoll","Ncoll",100,-0.5,99.5); 
    TH1D* NpartHist = new TH1D("Npart","Npart",37,-0.5,36.5);
    TH1D* NelemHist = new TH1D("Nelem","Nelem",30,0,30); 
    TH1D* SumETHist = new TH1D("SumET","SumET",100,0,500);
    TH1D* GammaHist = new TH1D("GammaHist","GammaHist",50,0,60);
    TH2D* FCALoverNpartvsNpartHist = new TH2D("FCALoverNpartvsNpartHist","FCaLoverNpartvsNpartHist",35,-.5,34.5,100,0,60);
    TH2D* NpartvsCentralityHist = new TH2D("NpartvsCentralityHist","NpartvsCentralityHist",100,0,150,35,-.5,34.5);
    TH2D* FCALvsCentralityHist = new TH2D("FCALvsCentralityHist","FCALvsCentralityHist",100,0,150,100,0,500);



    
    TH1D* Ncoll80100CentralityHist = new TH1D("Ncoll80100CentralityHist","Ncoll80100CentralityHist",71,-0.5,70.5);
    TH1D* Ncoll6080CentralityHist = new TH1D("Ncoll6080Centrality","Ncoll6080Centrality",71,-0.5,70.5);
    TH1D* Ncoll4060CentralityHist = new TH1D("Ncoll4060Centrality","Ncoll4060Centrality",71,-0.5,70.5);
    TH1D* Ncoll2040CentralityHist = new TH1D("Ncoll2040Centrality","Ncoll2040Centrality",71,-0.5,70.5);
    TH1D* Ncoll0020CentralityHist = new TH1D("Ncoll0020Centrality","Ncoll0020Centrality",71,-0.5,70.5);


    double k_pp     = 2.25;   
    double theta_pp = 6.53;  
    double k_oo =1.95;
    double theta_oo = 7.43;
    double sigma_noise = 1.5; 

    int Events=0;

    while(tree_reader.Next()){
        Events++;
        NcollHist->Fill(Ncoll[0]);
        NpartHist->Fill(Npart[0]);

        double x=0.11;
        double Nelem=x*Ncoll[0]+(1-x)*Npart[0]/2;
        NelemHist->Fill(Nelem);
        GammaHist->Fill(SampleGamma(k_pp, theta_pp));

        double ET = SampleGamma(k_oo*Nelem, theta_oo);

        ET += gRandom->Gaus(0, sigma_noise);
        SumETHist->Fill(ET);
        if(ET>20.215) FCALoverNpartvsNpartHist->Fill(Npart[0],ET/Npart[0]);
        double centrality=(ImpPara[0]*ImpPara[0]);
        
        if(ET<20.215) //80%-60% 
        {
            Ncoll80100CentralityHist->Fill(Ncoll[0]);
            NpartvsCentralityHist->Fill(centrality,Npart[0]);
            FCALvsCentralityHist->Fill(centrality,ET);
        }
        else if(ET<44.348 && ET>20.215) //80%-60% 
        {
            Ncoll6080CentralityHist->Fill(Ncoll[0]);
        }
        else if(ET<85.128 && ET>44.348) //60%-40% 
        {
            Ncoll4060CentralityHist->Fill(Ncoll[0]);
        }
        else if(ET<153.033 && ET>85.128) //40%-20% 
        {
            Ncoll2040CentralityHist->Fill(Ncoll[0]);
        }
        else if(ET>153.033) //20%-0% 
        {
            Ncoll0020CentralityHist->Fill(Ncoll[0]);
        }
        
    }
    TProfile *profFCALoverNpart = FCALoverNpartvsNpartHist->ProfileX("profFCALoverNpart");
    TProfile *profCentralityvsNpartHist = NpartvsCentralityHist->ProfileX("profCentralityvsNpartHist");

    NcollHist->Write();
    NpartHist->Write();
    NelemHist->Write();
    SumETHist->Write();
    GammaHist->Write();
    Ncoll80100CentralityHist->Write();
    Ncoll6080CentralityHist->Write();
    Ncoll4060CentralityHist->Write();
    Ncoll2040CentralityHist->Write();
    Ncoll0020CentralityHist->Write();
    FCALoverNpartvsNpartHist->Write();
    profFCALoverNpart->Write();
    NpartvsCentralityHist->Write();
    profCentralityvsNpartHist->Write();
    FCALvsCentralityHist->Write();


    Output->Close();

}