#include <TGraph.h>
#include <iostream>
#include <TGraph.h>
#include <vector>
#include <algorithm>
#include "TH1.h"

bool savepdf = 0;

std::pair<double, double> GetYRange(const std::vector<TH1*>& hists) {
    double minY = std::numeric_limits<double>::max();
    double maxY = -std::numeric_limits<double>::max();

    for (auto* hist : hists) {
        if (!hist) continue;

        for (int i = 1; i <= hist->GetNbinsX(); ++i) {
            double y = hist->GetBinContent(i);
            if (y < minY) minY = y;
            if (y > maxY) maxY = y;
        }

        // Also consider underflow and overflow bins if relevant:
        /*
        double under = hist->GetBinContent(0);
        double over  = hist->GetBinContent(hist->GetNbinsX() + 1);
        minY = std::min(minY, under);
        maxY = std::max(maxY, under);
        minY = std::min(minY, over);
        maxY = std::max(maxY, over);
        */
    }

    return {minY, maxY};
}

void TruncateGraphAtX(TGraph* g, double x_max = 0.0315) {
    if (!g) return;
    std::cout << "Total N = " << g->GetN() << std::endl;

    std::vector<double> x_vals, y_vals;
    double x, y;

    int N = g->GetN();
    for (int i = 0; i < N; ++i) {
        g->GetPoint(i, x, y);
        if (x > x_max) break;  
        x_vals.push_back(x);
        y_vals.push_back(y);
    }

    g->Set(0);  
    for (size_t i = 0; i < x_vals.size(); ++i)
        g->SetPoint(i, x_vals[i], y_vals[i]);

    std::cout << "Truncated graph to " << x_vals.size() << " points with x <= " << x_max << std::endl;
    std::cout << "Total N = " << g->GetN() << std::endl;
}



const int N = 2;
const int Nbin = 11;
TFile * f[N];

TString filename[N] = {
	/*
	// 26July25 archive
	"../result/FairMUanalyzer_18July25_single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF0_output.root",
	"../result/FairMUanalyzer_18July25_single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root",
	"../result/FairMUanalyzer_single_muon_interaction_0_CDbugfix11July25_run8_21July25_muedaq04-1750227094-1750228937_MF1_maxNumberOfSharedHits2_output.root",
	"../result/FairMUanalyzer_18July25_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF0_output.root",
	"../result/FairMUanalyzer_18July25_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root",
	"../result/FairMUanalyzer_single_muon_interaction_1_CDbugfix11July25_run8_21July25_muedaq04-1750227094-1750228937_MF1_maxNumberOfSharedHits2_output.root"
	*/

	//9Aug25 backup, before Milosz. All my own generation
	/*
	"../result-26July25/FairMUanalyzer_18July25_single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF0_output.root",
	"../result-26July25/FairMUanalyzer_18July25_single_muon_interaction_0_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root",
	"../result-26July25/FairMUanalyzer_single_muon_interaction_0_CDbugfix11July25_run8_21July25_muedaq04-1750227094-1750228937_MF1_maxNumberOfSharedHits2_output.root",
	"../result-26July25/FairMUanalyzer_18July25_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF0_output.root",
	"../result-06Aug25/FairMUanalyzer_18July25_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_NoTightTrackCutTgt2_output.root",
	"../result-06Aug25/FairMUanalyzer_single_muon_interaction_1_CDbugfix11July25_run8_21July25_muedaq04-1750227094-1750228937_MF1_maxNumberOfSharedHits2_NoTightTrackCutTgt2_output.root"
	*/
	//9Aug25 backup, run17 by Milosz.
	//"../result/FairMUanalyzer_06Aug25_Milosz_run17_reconstruction_0_CDbugfix11July25_MF1_output.root",
	//"../result/FairMUanalyzer_06Aug25_Milosz_run17_reconstruction_1_CDbugfix11July25_MF1_output.root"

	//21Aug first look at run24
	"../result/FairMUanalyzer_21Aug25_run24_reconstruction_0_CDbugfix11July25_MF1_output.root",
	"../result/FairMUanalyzer_21Aug25_run24_reconstruction_1_CDbugfix11July25_MF1_output.root"

	
};
TString savebasename[N];

TString label[Nbin] = {"total","golden","t0mem","t0mee","t0mmm","t0me<m","t1mem","t1mee","t1mmm","t1me<m",""};

TH1F * hCaseDist[N];
TH2D * h2d[N];
TH2D * h2d_ref[N];
TGraph * g_elastic[N];

TH1D* h1_x[N];
TH1D* h1_y[N];

TH1D* ratio_x[N];
TH1D* ratio_y[N];

TLegend *lx[N];
TLegend *ly[N];

TCanvas * c1[N];
TCanvas * c2[N];
TCanvas * c3[N];

void ana_h2d(){

	//gStyle->SetOptStat(0); 
	gStyle->SetOptStat(1110); // entries, mean, RMS, name

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

		cout<<"t0_allscat/total: "<<t0_allscat*1.0/hCaseDist[i]->GetBinContent(1)<<endl;
		cout<<"t1_allscat/total: "<<t1_allscat*1.0/hCaseDist[i]->GetBinContent(1)<<endl;

		
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
		
		gPad->Update();            
		TPaveStats* stats = (TPaveStats*)h2d[i]->FindObject("stats");
		if (stats) {
		    stats->SetName("nostatsname");
		}
		
		//save the plots

		savebasename[i] = filename[i];
		TRegexp re("^\\.\\./result[^/]*/");
		savebasename[i](re) = ""; 
		//savebasename[i].ReplaceAll("../result/","");
		savebasename[i].ReplaceAll(".root","");

		if(savepdf)c1[i]->SaveAs(Form("h2d_%s.pdf",savebasename[i].Data()));
		

		h1_x[i] = h2d[i]->ProjectionX(Form("h1_x_%i",i));
		h1_x[i]->Rebin(5);
		h1_x[i]->SetTitle("Projection to electron Angle");
		h1_y[i] = h2d[i]->ProjectionY(Form("h1_y_%i",i));
		h1_y[i]->Rebin(5);
		h1_y[i]->SetTitle("Projection to muon Angle");

		c2[i] = new TCanvas(Form("c2_%i",i),Form("c2_%i",i),1400*0.8,500*0.8);
		c2[i]->Divide(2,1);
		c2[i]->cd(1);h1_x[i]->Draw();
		c2[i]->cd(2);h1_y[i]->Draw();

		
		if(i%3==2){

			std::vector<TH1*> hvec = {h1_x[i-2], h1_x[i-1], h1_x[i]};
			auto [ymin, ymax] = GetYRange(hvec);
			std::vector<TH1*> hvec2 = {h1_y[i-2], h1_y[i-1], h1_y[i]};
			auto [ymin2, ymax2] = GetYRange(hvec2);

			//gStyle->SetOptStat(1);


			c3[i] = new TCanvas(Form("c3_%i",i),Form("c3_%i",i),1400*0.7,1000*0.7);
			c3[i]->Divide(2,2);
			c3[i]->cd(1);	
			h1_x[i-2]->SetMinimum(ymin);
			h1_x[i-2]->SetMaximum(1.1*ymax);
			h1_x[i-2]->SetLineColor(kBlack);h1_x[i-2]->Draw();
			h1_x[i-1]->SetLineColor(kBlue);h1_x[i-1]->Draw("same");
			h1_x[i]->SetLineColor(kRed);h1_x[i]->Draw("same");

			lx[i] = new TLegend(0.7,0.6,0.9,0.8);
			lx[i]->AddEntry(h1_x[i-2],"MF0");
			lx[i]->AddEntry(h1_x[i-1],"MF1; SH0");
			lx[i]->AddEntry(h1_x[i],"MF1; SH2");
			lx[i]->Draw();

			
			c3[i]->cd(2);
			h1_y[i-2]->SetMinimum(ymin2);
			h1_y[i-2]->SetMaximum(1.1*ymax2);
			h1_y[i-2]->SetLineColor(kBlack);h1_y[i-2]->Draw();
			h1_y[i-1]->SetLineColor(kBlue);h1_y[i-1]->Draw("same");
			h1_y[i]->SetLineColor(kRed);h1_y[i]->Draw("same");

			ly[i] = new TLegend(0.7,0.6,0.9,0.8);
			ly[i]->AddEntry(h1_y[i-2],"MF0");
			ly[i]->AddEntry(h1_y[i-1],"MF1; SH0");
			ly[i]->AddEntry(h1_y[i],"MF1; SH2");
			ly[i]->Draw();

			ratio_x[i] = (TH1D*)h1_x[i]->Clone(Form("ratio_x_%i",i));
			ratio_x[i]->Divide(h1_x[i],h1_x[i-1],1.0,1.0,"B"); //"B" for beysian propagation
			c3[i]->cd(3);
			ratio_x[i]->SetTitle("SH2 / SH0 for projection to electron Angle");
			ratio_x[i]->SetLineColor(kBlack);
			ratio_x[i]->SetMarkerStyle(20);
			ratio_x[i]->SetMarkerSize(0.5);
			ratio_x[i]->Draw("E");

			ratio_y[i] = (TH1D*)h1_y[i]->Clone(Form("ratio_y_%i",i));
			ratio_y[i]->Divide(h1_y[i],h1_y[i-1],1.0,1.0,"B");
			c3[i]->cd(4);
			ratio_y[i]->SetTitle("SH2 / SH0 for projection to muon Angle");
			ratio_y[i]->SetLineColor(kBlack);
			ratio_y[i]->SetMarkerStyle(20);
			ratio_y[i]->SetMarkerSize(0.5);
			ratio_y[i]->Draw("E");
		}

		

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