// draw_case_macro.C
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

//------------------------------------------------------------
// 生成目标 hist 名称列表
//------------------------------------------------------------
std::vector<TString> getHistNames(const TString& type) {
    std::vector<TString> names;

    // 2D
    names.push_back("case_h2d_" + type);
    names.push_back("case_h2d_bstvtx_" + type);

    // 1D
    names.push_back("case_h1d_vertex_" + type);
    names.push_back("case_h1d_Vtxchi2_" + type);
    names.push_back("case_h1d_aco_" + type);
    names.push_back("case_h1d_Leftoverhits0_" + type);
    names.push_back("case_h1d_Leftoverhits1_" + type);
    names.push_back("case_h1d_Leftoverhits2_" + type);

    // perModule (6 modules)
    for (int j = 0; j < 6; j++) {
        names.push_back(Form("case_h1d_LeftoverhitsPerModule_stat0_Module_%d_%s", j, type.Data()));
        names.push_back(Form("case_h1d_LeftoverhitsPerModule_stat1_Module_%d_%s", j, type.Data()));
        names.push_back(Form("case_h1d_LeftoverhitsPerModule_stat2_Module_%d_%s", j, type.Data()));
    }

    return names;
}

//------------------------------------------------------------
// 主函数：读文件 + 画图
//------------------------------------------------------------
void draw_case_macro(const TString& targetType="t1mmm") {

    std::vector<TString> filenames = {
        "../result/FairMUanalyzer_13Sep25_GAtesttest_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root"
    };

    // 需要画的 hist 名称
    auto histNames = getHistNames(targetType);

    // 遍历文件
    for (size_t i = 0; i < filenames.size(); i++) {

        TFile *f = TFile::Open(filenames[i], "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Error: cannot open " << filenames[i] << std::endl;
            continue;
        }

        // 从文件名提取 GA 编号
        TRegexp re("GAtest[0-9]+");
        TString match;
        if (filenames[i].Index(re) != kNPOS) {
            match = filenames[i](re);  // "GAtest603"
            match.ReplaceAll("GAtest", ""); // "603"
        } else {
            match = Form("f%d", (int)i);
        }

        // 每个 hist 单独一个 canvas
        for (size_t j = 0; j < histNames.size(); j++) {
            TString hname = histNames[j];
            TObject* obj = f->Get(hname);
            if (!obj) {
                std::cerr << "Warning: histogram " << hname 
                          << " not found in " << filenames[i] << std::endl;
                continue;
            }

            TCanvas *c = new TCanvas(Form("c_%s_%s", hname.Data(), match.Data()),
                                     hname, 1000, 800);

            if (hname.BeginsWith("case_h2d")) {
                ((TH2*)obj)->Draw("COLZ");
            } else {
                ((TH1*)obj)->Draw();
            }

            c->SaveAs(Form("%s_%s.pdf", hname.Data(), match.Data()));
            c->SaveAs(Form("%s_%s.png", hname.Data(), match.Data()));
        }

        f->Close();
    }
}
