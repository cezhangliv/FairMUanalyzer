#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <vector>

void draw_three_files() {
    
    std::vector<TString> filenames = {
        "../result/FairMUanalyzer_04Sep25_GAtest502_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root",
        "../result/FairMUanalyzer_04Sep25_GAtest603_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root",
        "../result/FairMUanalyzer_04Sep25_GAtest604_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root"
    };

    
    std::vector<TString> histonames = {
        "h2d_t1all",
        "h2d_t1mem",
        "h2d_t1mmm",
        "h2d_t1mee",
        "h2d_bstvtx_t1all",
        "h2d_bstvtx_t1mem",
        "h2d_bstvtx_t1mmm",
        "h2d_bstvtx_t1mee"
    };

    for (size_t i = 0; i < filenames.size(); i++) {
        TFile *f = TFile::Open(filenames[i], "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Error: cannot open " << filenames[i] << std::endl;
            continue;
        }

        TRegexp re("GAtest[0-9]+");
        TString match;
        if (filenames[i].Index(re) != kNPOS) {
            match = filenames[i](re);  // 结果是 "GAtest603"
            match.ReplaceAll("GAtest", ""); // 结果就是 "603"
        }

        TString cname = TString::Format("c%s", match.Data());
        TCanvas *c = new TCanvas(cname, filenames[i], 1600, 1000);
        c->Divide(4, 2);

        for (size_t j = 0; j < histonames.size(); j++) {
            TH2D *h = (TH2D*)f->Get(histonames[j]);
            if (!h) {
                std::cerr << "Warning: histogram " << histonames[j] 
                          << " not found in " << filenames[i] << std::endl;
                continue;
            }
            c->cd(j+1);
            h->Draw("COLZ");
        }

        c->SaveAs(TString::Format("%s.pdf", cname.Data()));
        f->Close();
    }
}
