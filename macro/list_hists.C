#include <TFile.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <iostream>

void list_hists(const char* filename = "/Users/zhangce/WorkArea/MUonE/TB25/FairMUanalyzer/result/FairMUanalyzer_v0.17.6_Data2025_run32_singleMu1_sharedHit0_recov0.17.6_output.root") {
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::cout << "=== Histograms in file: " << filename << " ===" << std::endl;

    TIter nextkey(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)nextkey())) {
        TObject *obj = key->ReadObj();
        if (!obj) continue;

        // 判断是否是 histogram
        if (obj->InheritsFrom("TH1")) {
            std::cout << "Name: " << obj->GetName()
                      << " | Class: " << obj->ClassName();

            TH1* h = (TH1*)obj;
            std::cout << " | Nbins: " << h->GetNbinsX();

            if (obj->InheritsFrom("TH2")) {
                TH2* h2 = (TH2*)obj;
                std::cout << " x " << h2->GetNbinsY();
            }

            std::cout << std::endl;
        }
    }

    f->Close();
}
