// File: add_TH2D.C
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"

void add_TH2D() {
    // 输入 ROOT 文件
    TFile *f = TFile::Open("../result/FairMUanalyzer_29Aug25_GAtest204_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root", "READ");
    if (!f || f->IsZombie()) {
        printf("Error: cannot open file!\n");
        return;
    }

    // 取出三张二维直方图
    TH2D *h_mmm = (TH2D*)f->Get("h2d_t1mmm");
    TH2D *h_mem = (TH2D*)f->Get("h2d_t1mem");
    TH2D *h_mee = (TH2D*)f->Get("h2d_t1mee");

    if (!h_mmm || !h_mem || !h_mee) {
        printf("Error: one or more histograms not found!\n");
        return;
    }

    // 复制第一张图作为总图
    TH2D *h_all = (TH2D*)h_mmm->Clone("h2d_all");
    h_all->SetTitle("mem+mmm+mee");

    // 累加另外两张
    h_all->Add(h_mem);
    h_all->Add(h_mee);

    // 画出来
    TCanvas *c1 = new TCanvas("c1", "h2d_all", 800, 600);
    h_all->Draw("COLZ");

    // 保存结果到新的 root 文件
    //TFile *fout = new TFile("combined_output.root", "RECREATE");
    //h_all->Write();
    //fout->Close();

    //printf("Done! Combined histogram saved to combined_output.root\n");
}
