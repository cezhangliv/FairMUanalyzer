#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <sys/stat.h>
#include <iostream>
#include <vector>

void compare_hists(
    //const char* file1 = "../result/25Nov_FairMUanalyzer_TEST_v0.17.2_Data2025_run32_passingMu_sharedHit0_recoTEST_v0.17.2_output.root",
    //const char* file2 = "../result/25Nov_FairMUanalyzer_v0.17.6_Data2025_run32_passingMu_sharedHit0_recov0.17.6_output.root"

    //const char* file1 = "../result/25Nov_FairMUanalyzer_TEST_v0.17.2_Data2025_run32_singleMu0_sharedHit2_recoTEST_v0.17.2_output.root",
    //const char* file2 = "../result/25Nov_FairMUanalyzer_v0.17.6_Data2025_run32_singleMu0_sharedHit2_recov0.17.6_output.root"

    const char* file1 = "../result/25Nov_FairMUanalyzer_TEST_v0.17.2_Data2025_run32_singleMu1_sharedHit2_recoTEST_v0.17.2_output.root",
    const char* file2 = "../result/25Nov_FairMUanalyzer_v0.17.6_Data2025_run32_singleMu1_sharedHit2_recov0.17.6_output.root"
) {
    
    std::vector<std::string> targetHists = {
        /*
        "h_res_allTrack_st0Module0",
        "h_res_allTrack_st0Module1",
        "h_res_allTrack_st0Module2",
        "h_res_allTrack_st0Module3",
        "h_res_allTrack_st1Module0",
        "h_res_allTrack_st1Module1",
        "h_res_allTrack_st1Module2",
        "h_res_allTrack_st1Module3",
        "h_res_allTrack_st2Module0",
        "h_res_allTrack_st2Module1",
        "h_res_allTrack_st2Module2",
        "h_res_allTrack_st2Module3",
        "h_hits_muon_1mu",
        "h_hitsModuleID_zcut_1mu",
        "h_hits_zcut (MF hits)",
        "h_Ntracks",
        "h_isMuon"
        */

        /*
        "h_Ntracks",
        "h_Nhits0",
        "h_Nhits1",
        "h_Nhits2",

        "case_h2d_t0all",
        "case_h2d_t0me<m",
        "case_h2d_t0mee",
        "case_h2d_t0mem",
        "case_h2d_t0mmm",

        "case_h1d_vertex_golden",
        "case_h1d_vertex_t0all",
        "case_h1d_vertex_t0me<m",
        "case_h1d_vertex_t0mee",
        "case_h1d_vertex_t0mem",
        "case_h1d_vertex_t0mmm"
        */

        "h_Ntracks",
        "h_Nhits0",
        "h_Nhits1",
        "h_Nhits2",

        "case_h2d_t1all",
        "case_h2d_t1me<m",
        "case_h2d_t1mee",
        "case_h2d_t1mem",
        "case_h2d_t1mmm",

        "case_h1d_vertex_golden",
        "case_h1d_vertex_t1all",
        "case_h1d_vertex_t1me<m",
        "case_h1d_vertex_t1mee",
        "case_h1d_vertex_t1mem",
        "case_h1d_vertex_t1mmm"


    };

    // 创建输出目录 ./compare/
    system("mkdir -p compare");

    TFile* f1 = TFile::Open(file1);
    TFile* f2 = TFile::Open(file2);
    if (!f1 || !f2 || f1->IsZombie() || f2->IsZombie()) {
        std::cerr << "Error opening files." << std::endl;
        return;
    }

    for (const auto& name : targetHists) {
        TObject* obj1 = f1->Get(name.c_str());
        TObject* obj2 = f2->Get(name.c_str());

        if (!obj1 || !obj2) {
            std::cerr << "[WARNING] " << name << " not found in both files." << std::endl;
            continue;
        }

        if (!obj1->InheritsFrom("TH1") || !obj2->InheritsFrom("TH1")) {
            std::cerr << "[WARNING] " << name << " is not TH1 type." << std::endl;
            continue;
        }

        TH1* h1 = (TH1*)obj1;
        TH1* h2 = (TH1*)obj2;

        // 创建 canvas
        TCanvas* c = new TCanvas(("c_"+name).c_str(), name.c_str(), 500, 600);
        c->Divide(1,2);

        // pad1: overlay
        c->cd(1);
        h1->SetLineColor(kRed);
        h2->SetLineColor(kBlue);
        h1->Draw();
        

        TLegend* leg = new TLegend(0.15,0.75,0.38,0.88);
        //leg->AddEntry(h1, (std::string(file1)+" : "+name).c_str());
        //leg->AddEntry(h2, (std::string(file2)+" : "+name).c_str());
        //leg->AddEntry(h1, ("v0.17.2 passing muon: "+name).c_str());
        //leg->AddEntry(h2, ("v0.17.6 passing muon: "+name).c_str());
        leg->AddEntry(h1, "v0.17.2 passing muon");
        leg->AddEntry(h2, "v0.17.6 passing muon");
        leg->Draw();

        // pad2: difference: h1-h2
        c->cd(2);
        h2->Draw("SAME");
        //TH1* hdiff = (TH1*)h1->Clone(("diff_"+name).c_str());
        //hdiff->Add(h2, -1);
        //hdiff->SetLineColor(kBlack);
        //hdiff->SetTitle((name + " difference (A-B)").c_str());
        //hdiff->Draw();

        // 保存图像
        c->SaveAs(("compare/" + name + ".png").c_str());

        delete c;
    }

    f1->Close();
    f2->Close();
}
