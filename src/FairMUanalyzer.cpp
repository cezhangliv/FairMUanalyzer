#include "FairMUanalyzer.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include <iostream>
#include <set>
#include <algorithm>
#include <cmath>

FairMUanalyzer::FairMUanalyzer() : inputFile_(nullptr), cbmsim_(nullptr), reco_(nullptr), MuonFilterHits_(3), outputPrefix_("result/FairMUanalyzer"), savepdf_(true) {
    
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(1111);

    for (int i = 0; i < 3; ++i) {
        h_residual_hitOnTrack[i] = new TH1F(Form("h_res_onTrack_st%d", i), Form("Residual (OnTrack) - Station %d;Distance [cm];Entries", i), 200, -2, 2);
        h_residual_hitOffTrack[i] = new TH1F(Form("h_res_offTrack_st%d", i), Form("Residual (OffTrack) - Station %d;Distance [cm];Entries", i), 200, -2, 2);
        for (int j = 0; j < 4; ++j) {
            h_residual_hitOnTrackModule[i][j] = new TH1F(Form("h_res_onTrack_st%dModule%d", i, j), Form("Residual (OnTrack) - Station %d, Module %d;Distance [cm];Entries", i, j), 200, -2, 2);
            h_residual_hitOffTrackModule[i][j] = new TH1F(Form("h_res_offTrack_st%dModule%d", i, j), Form("Residual (OffTrack) - Station %d, Module %d;Distance [cm];Entries", i, j), 200, -2, 2);
            h_residual_hitAllTrackModule[i][j] = new TH1F(Form("h_res_allTrack_st%dModule%d", i, j), Form("Residual (AllTracks) - Station %d, Module %d;Distance [cm];Entries", i, j), 200, -2, 2);
        }
    }
}

FairMUanalyzer::~FairMUanalyzer() {
    if (inputFile_) inputFile_->Close();
}

void FairMUanalyzer::SetInputFile(const std::string& path) {
    inputFilePath_ = path;
}

void FairMUanalyzer::SetOutputPrefix(const std::string& prefix) {
    outputPrefix_ = prefix.c_str();
}

void FairMUanalyzer::SetMuonFilterHits(int val) {
    MuonFilterHits_ = val;
}

void FairMUanalyzer::Run() {
    Init();
    Analyze();
    SaveResults();
}

void FairMUanalyzer::Init() {
    inputFile_ = new TFile(inputFilePath_.c_str());
    if (!inputFile_ || inputFile_->IsZombie()) {
        std::cerr << "Failed to open input file: " << inputFilePath_ << std::endl;
        return;
    }

    cbmsim_ = (TTree*) inputFile_->Get("cbmsim");
    if (!cbmsim_) {
        std::cerr << "Cannot find cbmsim tree in file: " << inputFilePath_ << std::endl;
        return;
    }

    cbmsim_->SetBranchAddress("ReconstructionOutput", &reco_);
}

TVector3 FairMUanalyzer::getXYfromHitMF(const MUonERecoOutputHitAnalysis& hit) {
    double pos = hit.positionPerpendicular();
    double x = 0, y = 0;
    int moduleID = hit.moduleID();

    switch(moduleID) {
        case 0:
        case 2:
            x = pos; y = 0; break;
        case 1:
        case 3:
            x = 0; y = pos; break;
        default:
            std::cerr << "Unknown moduleID: " << moduleID << std::endl;
            break;
    }

    return TVector3(x, y, hit.z());
}

double FairMUanalyzer::computeSigned2DResidualMF(const TVector3& p3D, const TVector3& x0, const TVector3& h, int moduleID) {
    TVector3 p = p3D;
    TVector3 delta = h - x0;
    TVector3 normal;

    switch (moduleID) {
        case 0:
        case 2:
            p.SetY(0); delta.SetY(0); normal = TVector3(0, 1, 0); break;
        case 1:
        case 3:
            p.SetX(0); delta.SetX(0); normal = TVector3(1, 0, 0); break;
        default:
            std::cerr << "Unknown moduleID!" << std::endl; return 0;
    }

    TVector3 cross = delta.Cross(p);
    return cross.Dot(normal) / p.Mag();
}

void FairMUanalyzer::Analyze() {
    AnalyzeMF();
}

void FairMUanalyzer::SaveResults() {
    TFile* fout = new TFile(Form("%s_output.root", outputPrefix_.Data()), "RECREATE");
    for (int i = 0; i < 3; ++i) {
        h_residual_hitOnTrack[i]->Write();
        h_residual_hitOffTrack[i]->Write();
        for (int j = 0; j < 4; ++j) {
            h_residual_hitOnTrackModule[i][j]->Write();
            h_residual_hitOffTrackModule[i][j]->Write();
            h_residual_hitAllTrackModule[i][j]->Write();
        }
    }
    fout->Close();

    if(!savepdf_)return;

    TCanvas* c3 = new TCanvas("c3", "Residuals by station", 1800, 400);
    c3->Divide(3, 1);
    for (int i = 0; i < 3; ++i) {
        c3->cd(i + 1);
        h_residual_hitOnTrack[i]->SetLineColor(kRed);
        h_residual_hitOffTrack[i]->SetLineColor(kBlue);
        h_residual_hitOnTrack[i]->Draw();
        h_residual_hitOffTrack[i]->Draw("sames");

        auto legend = new TLegend(0.7, 0.6, 0.9, 0.7);
        legend->AddEntry(h_residual_hitOnTrack[i], "On Track", "l");
        legend->AddEntry(h_residual_hitOffTrack[i], "Off Track", "l");
        legend->Draw();
    }
    c3->SaveAs(Form("%s_residuals_station.pdf", outputPrefix_.Data()));

    TCanvas* c4 = new TCanvas("c4", "Residuals by station/module", 2400, 1200);
    c4->Divide(4, 3);
    for (int i = 0; i < 12; ++i) {
        int ii = i % 4; // module
        int jj = i / 4; // station
        c4->cd(i + 1);
        h_residual_hitOnTrackModule[jj][ii]->SetLineColor(kRed);
        h_residual_hitOffTrackModule[jj][ii]->SetLineColor(kBlue);
        h_residual_hitOnTrackModule[jj][ii]->Draw();
        h_residual_hitOffTrackModule[jj][ii]->Draw("sames");

        auto legend = new TLegend(0.7, 0.6, 0.9, 0.7);
        legend->AddEntry(h_residual_hitOnTrackModule[jj][ii], "On Track", "l");
        legend->AddEntry(h_residual_hitOffTrackModule[jj][ii], "Off Track", "l");
        legend->Draw();
    }
    c4->SaveAs(Form("%s_residuals_station_module.pdf", outputPrefix_.Data()));
}
