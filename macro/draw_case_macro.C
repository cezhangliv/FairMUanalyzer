
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TRegexp.h>
#include <TString.h>
#include <TStyle.h>

#include <iostream>
#include <vector>
#include <string>

TH1* GetHistWithSemiColon(TFile* f, const TString& searchName) {
    TIter nextkey(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)nextkey())) {
        TString kname = key->GetName();
        if (kname == searchName) {
            return dynamic_cast<TH1*>(key->ReadObj());
        }
    }
    return nullptr;
}

std::vector<TString> getSingleHistNames(const TString& type) {
    std::vector<TString> names;

    names.push_back("case_h2d_" + type);
    names.push_back("case_h2d_bstvtx_" + type);

    names.push_back("case_h1d_vertex_" + type);
    names.push_back("case_h1d_Vtxchi2_" + type);
    names.push_back("case_h1d_aco_" + type);
    names.push_back("case_h1d_Leftoverhits0_" + type + ";Hits");
    names.push_back("case_h1d_Leftoverhits1_" + type + ";Hits");
    names.push_back("case_h1d_Leftoverhits2_" + type + ";Hits");

    return names;
}


std::vector<TString> getPerModuleNames(const TString& type, const TString& statTag) {
    std::vector<TString> names;
    for (int j = 0; j < 6; j++) {
        names.push_back(Form("case_h1d_LeftoverhitsPerModule_%s_Module_%d_%s;Hits",
                             statTag.Data(), j, type.Data()));
    }
    return names;
}


void draw_case_macro(const std::vector<TString>& targetTypes = {"Total","golden","t1mem","t1mee","t1mmm","t1mmmband","t1mmmoutofband"}) {


    std::vector<TString> filenames = {
        "../result/FairMUanalyzer_13Sep25_GAtest902_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root"
    };



    for (size_t i = 0; i < filenames.size(); i++) {
        TFile *f = TFile::Open(filenames[i], "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Error: cannot open " << filenames[i] << std::endl;
            continue;
        }
        /*
        f->cd();
        TIter nextkey(f->GetListOfKeys());
        TKey *key;
        while ((key = (TKey*)nextkey())) {
            std::cout << "Key name = [" << key->GetName() << "]" << std::endl;
        }
        */

        TRegexp re("GAtest[0-9]+");
        TString match;
        if (filenames[i].Index(re) != kNPOS) {
            match = filenames[i](re);  
            match.ReplaceAll("GAtest", ""); 
        } else {
            match = Form("f%d", (int)i);
        }

        gSystem->Exec(Form("mkdir -p %s", match.Data()));

        for (auto& type : targetTypes) {

            
            auto histNames = getSingleHistNames(type);
            for (auto& hname : histNames) {
                //TObject* obj = f->Get(hname);
                TObject* obj = GetHistWithSemiColon(f, hname);

                if (!obj) {
                    std::cerr << "Warning: histogram " << hname 
                              << " not found in " << filenames[i] << std::endl;
                    continue;
                }

                TCanvas *c = new TCanvas(Form("c_%s_%s_%s", hname.Data(), match.Data(), type.Data()),
                                         hname, 1000, 800);
                if (hname.BeginsWith("case_h2d")) {
                    ((TH2*)obj)->Draw("COLZ");
                } else {
                    ((TH1*)obj)->Draw();
                }
                c->SaveAs(Form("%s/%s_%s_%s.pdf", match.Data(), hname.Data(), match.Data(), type.Data()));
                c->SaveAs(Form("%s/%s_%s_%s.png", match.Data(), hname.Data(), match.Data(), type.Data()));

                
            }

            
            std::vector<TString> statTags = {"stat0","stat1","stat2"};
            for (auto& stat : statTags) {
                auto modNames = getPerModuleNames(type, stat);
                TCanvas *c = new TCanvas(Form("c_%s_%s_%s", stat.Data(), match.Data(), type.Data()),
                                         Form("perModule_%s_%s", stat.Data(), type.Data()),
                                         1400, 1000);
                c->Divide(3,2);

                for (int j = 0; j < 6; j++) {
                    //TObject* obj = f->Get(modNames[j]);
                    TObject* obj = GetHistWithSemiColon(f, modNames[j]);

                    if (!obj) {
                        std::cerr << "Warning: histogram " << modNames[j] 
                                  << " not found in " << filenames[i] << std::endl;
                        continue;
                    }
                    c->cd(j+1);
                    ((TH1*)obj)->Draw();
                }

                c->SaveAs(Form("%s/perModule_%s_%s_%s.pdf", match.Data(), stat.Data(), match.Data(), type.Data()));
                c->SaveAs(Form("%s/perModule_%s_%s_%s.png", match.Data(), stat.Data(), match.Data(), type.Data()));
            }

        } 

        f->Close();
    }
}

