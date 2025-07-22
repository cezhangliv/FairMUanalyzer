const int N = 4;
TFile * f[N];

TSting filename[N] = {
	"../result/FairMUanalyzer_18July25_single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF0_output.root",
	"../result/FairMUanalyzer_18July25_single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root",
	"../result/FairMUanalyzer_18July25_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF0_output.root",
	"../result/FairMUanalyzer_18July25_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root"
}

TH1F * hCaseDist[N];

void ana_CaseDist(){

	for(int i = 0; i<N; i++){

		f[i] = new TFile(filename[i].Data());
		hCaseDist[i] = (TH1F*) f[i]->Get("hCaseDist");
		hCaseDist[i]->SetName(Form("hCaseDist_%i",i));



	}



}