const int N = 4;
const int Nbin = 11;
TFile * f[N];

TString filename[N] = {
	"../result/FairMUanalyzer_18July25_single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF0_output.root",
	"../result/FairMUanalyzer_18July25_single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root",
	"../result/FairMUanalyzer_18July25_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF0_output.root",
	"../result/FairMUanalyzer_18July25_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root"
};

TString label[Nbin] = {"","total","golden","t0mem","t0mee","t0mmm","t0me<m","t1mem","t1mee","t1mmm","t1me<m"};

TH1F * hCaseDist[N];

void ana_CaseDist(){

	for(int i = 0; i<N; i++){

		cout<<filename[i].Data()<<endl;

		f[i] = new TFile(filename[i].Data());
		hCaseDist[i] = (TH1F*) f[i]->Get("hCaseDist");
		hCaseDist[i]->SetName(Form("hCaseDist_%i",i));

		for(int j = 1; j<= Nbin; j++)cout<<label[j]<<": "<<hCaseDist[i]->GetBinContent(j)<<endl;

		cout<<"golden/total: "<<hCaseDist[i]->GetBinContent(2)*1.0/hCaseDist[i]->GetBinContent(1)<<endl;
		
		int t0_allscat = hCaseDist[i]->GetBinContent(3) + hCaseDist[i]->GetBinContent(4) + hCaseDist[i]->GetBinContent(5);
		int t1_allscat = hCaseDist[i]->GetBinContent(7) + hCaseDist[i]->GetBinContent(8) + hCaseDist[i]->GetBinContent(9);

		cout<<"t0_allscat: "<<t0_allscat<<endl;
		cout<<"t1_allscat: "<<t1_allscat<<endl;

		/*
		cout<<"total: "<<hCaseDist[i]->GetBinContent(1)<<endl;
		cout<<"golden: "<<hCaseDist[i]->GetBinContent(2)<<endl;
		cout<<"t0mem: "<<hCaseDist[i]->GetBinContent(3)<<endl;
    	cout<<"t0mee: "<<hCaseDist[i]->GetBinContent(4)<<endl;
    	cout<<"t0mmm: "<<hCaseDist[i]->GetBinContent(5)<<endl;
    	cout<<"t0me<m: "<<hCaseDist[i]->GetBinContent(6)<<endl;
    	cout<<"t1mem: "<<hCaseDist[i]->GetBinContent(7)<<endl;
    	cout<<"t1mee: "<<hCaseDist[i]->GetBinContent(8)<<endl;
    	cout<<"t1mmm: "<<hCaseDist[i]->GetBinContent(9)<<endl;
    	cout<<"t1me<m: "<<hCaseDist[i]->GetBinContent(10)<<endl;
		*/

	}



}