#include <TGraph.h>
#include <iostream>

#include <TGraph.h>
#include <vector>
void TruncateGraphAtX(TGraph* g, double x_max = 0.0315) {
    if (!g) return;
    std::cout << "Total N = " << g->GetN() << std::endl;

    std::vector<double> x_vals, y_vals;
    double x, y;

    int N = g->GetN();
    for (int i = 0; i < N; ++i) {
        g->GetPoint(i, x, y);
        if (x > x_max) break;  // 提前终止
        x_vals.push_back(x);
        y_vals.push_back(y);
    }

    g->Set(0);  // 清空原图
    for (size_t i = 0; i < x_vals.size(); ++i)
        g->SetPoint(i, x_vals[i], y_vals[i]);

    std::cout << "Truncated graph to " << x_vals.size() << " points with x <= " << x_max << std::endl;
    std::cout << "Total N = " << g->GetN() << std::endl;
}


const int N = 6;
const int Nbin = 11;
TFile * f[N];

TString filename[N] = {
	"../result/FairMUanalyzer_18July25_single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF0_output.root",
	"../result/FairMUanalyzer_18July25_single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root",
	"../result/FairMUanalyzer_single_muon_interaction_0_CDbugfix11July25_run8_21July25_muedaq04-1750227094-1750228937_MF1_maxNumberOfSharedHits2_output.root",
	"../result/FairMUanalyzer_18July25_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF0_output.root",
	"../result/FairMUanalyzer_18July25_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root",
	"../result/FairMUanalyzer_single_muon_interaction_1_CDbugfix11July25_run8_21July25_muedaq04-1750227094-1750228937_MF1_maxNumberOfSharedHits2_output.root"
	
};
TString savebasename[N];

TString label[Nbin] = {"total","golden","t0mem","t0mee","t0mmm","t0me<m","t1mem","t1mee","t1mmm","t1me<m",""};

TH1F * hCaseDist[N];
TH2D * h2d[N];
TH2D * h2d_ref[N];
TGraph * g_elastic[N];
TCanvas * c1[N];

void ana_h2d(){

	for(int i = 0; i<N; i++){

		cout<<filename[i].Data()<<endl;

		f[i] = new TFile(filename[i].Data());
		hCaseDist[i] = (TH1F*) f[i]->Get("hCaseDist");
		hCaseDist[i]->SetName(Form("hCaseDist_%i",i));

		for(int j = 1; j<= Nbin; j++)cout<<label[j-1]<<": "<<hCaseDist[i]->GetBinContent(j)<<endl;

		cout<<"golden/total: "<<hCaseDist[i]->GetBinContent(2)*1.0/hCaseDist[i]->GetBinContent(1)<<endl;
		
		int t0_allscat = hCaseDist[i]->GetBinContent(3) + hCaseDist[i]->GetBinContent(4) + hCaseDist[i]->GetBinContent(5);
		int t1_allscat = hCaseDist[i]->GetBinContent(7) + hCaseDist[i]->GetBinContent(8) + hCaseDist[i]->GetBinContent(9);

		cout<<"t0_allscat: "<<t0_allscat<<endl;
		cout<<"t1_allscat: "<<t1_allscat<<endl;

		
		h2d[i] = (TH2D*)f[i]->Get("h_2d");
		h2d[i]->SetName(Form("h_2d_%i",i));
		
		h2d_ref[i] = (TH2D*)f[i]->Get("h_2d_ref");
		h2d_ref[i]->SetName(Form("h_2d_ref_%i",i));
		g_elastic[i] = (TGraph*)f[i]->Get("g_elastic");
		g_elastic[i]->SetName(Form("g_elastic_%i",i));
		TruncateGraphAtX(g_elastic[i]);
		c1[i] = new TCanvas(Form("c1_%i",i),Form("c1_%i",i));
		h2d[i]->Draw("colz");
		h2d_ref[i]->Draw("samecolscat");
		g_elastic[i]->Draw("sameL");
		
		//save the plots

		savebasename[i] = filename[i];
		savebasename[i].ReplaceAll("../result/","");
		savebasename[i].ReplaceAll(".root","");

		c1[i]->SaveAs(Form("h2d_%s.pdf",savebasename[i].Data()));
		

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