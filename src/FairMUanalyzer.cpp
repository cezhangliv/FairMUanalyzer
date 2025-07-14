#include "FairMUanalyzer.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include <iostream>
#include <set>
#include <algorithm>
#include <cmath>
/*
double me=0.51099906e-3;
double mu=105.65836900e-3;
double Ebeam=160;

double Eevsth(double *x, double *par){
    double th = *x;
    double Emu = par[0];
    double r = sqrt(Emu*Emu-mu*mu)/(Emu+me);
    return  me*(1.+r*r*cos(th)*cos(th))/(1-r*r*cos(th)*cos(th));
}


double thmu_vs_the(double *x, double *par){
    //double th = *x/1000.;
    double th = *x/1.;
    double Emu = par[0];
    double pmu = sqrt(Emu*Emu-mu*mu);
    double Ee = Eevsth(&th, par);
    double pe = sqrt(Ee*Ee-me*me);
    //return 1000.*acos((pmu-pe*cos(th))/sqrt(pe*pe+pmu*pmu-2*pmu*pe*cos(th)));
    return 1.*acos((pmu-pe*cos(th))/sqrt(pe*pe+pmu*pmu-2*pmu*pe*cos(th)));
}
*/

FairMUanalyzer::FairMUanalyzer() : inputFile_(nullptr), cbmsim_(nullptr), reco_(nullptr), MuonFilterHits_(3), outputPrefix_("result/FairMUanalyzer"), savepdf_(true) {
    
    goldenevents_ = 0;

    f_elastic = new TF1("f_elastic", thmu_vs_the,  0., 0.032,1);
    f_elastic->SetParameter(0, Ebeam_);
    f_elastic->SetNpx(1e5);

    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(1111);

    h_hits_zcut = new TH1F("h_hits_zcut", "Number of hits in MF; Hits;Entries", 20, 0, 20);
    h_isMuon = new TH1F("h_isMuon", "Number of 'Muon' tracks; isMuon tracks per event;Entries", 20, 0, 20);
    h_Ntracks = new TH1F("h_Ntracks", "Tracks multiplicity; tracks per event;Entries", 20, 0, 20);
    
    h_hitsModuleID_zcut[5];
    h_hitsPerMuonTrack_zcut[5];

    for (int i = 1; i <= 4; i++) {
        h_hitsModuleID_zcut[i] = new TH1F(Form("h_hitsModuleID_zcut_%imu", i),
                                              Form("MF hits module ID related to muon track (%i muti trks);Module ID;Entries", i),
                                              5, 0, 5);
        h_hitsPerMuonTrack_zcut[i] = new TH1F(Form("h_hits_muon_%imu", i),
                                              Form("MF hits related to muon track (%i muti trks);N hits;Entries", i),
                                              10, 0, 10);
    }

    h_goldenMuon_isMuon[3];
    for (int i = 0; i < 3; i++) {
        h_goldenMuon_isMuon[i] = new TH1F( Form("h_goldenMuon_isMuon_atStat%i",i), Form("golden muon isMuon? at station %i;isMuon;Entries",i), 2, 0, 2);
    }

    for (int i = 0; i < 3; ++i) {
        h_residual_hitOnTrack[i] = new TH1F(Form("h_res_onTrack_st%d", i), Form("Residual (OnTrack) - Station %d;Distance [cm];Entries", i), 200, -2, 2);
        h_residual_hitOffTrack[i] = new TH1F(Form("h_res_offTrack_st%d", i), Form("Residual (OffTrack) - Station %d;Distance [cm];Entries", i), 200, -2, 2);
        for (int j = 0; j < 4; ++j) {
            h_residual_hitOnTrackModule[i][j] = new TH1F(Form("h_res_onTrack_st%dModule%d", i, j), Form("Residual (OnTrack) - Station %d, Module %d;Distance [cm];Entries", i, j), 200, -2, 2);
            h_residual_hitOffTrackModule[i][j] = new TH1F(Form("h_res_offTrack_st%dModule%d", i, j), Form("Residual (OffTrack) - Station %d, Module %d;Distance [cm];Entries", i, j), 200, -2, 2);
            h_residual_hitAllTrackModule[i][j] = new TH1F(Form("h_res_allTrack_st%dModule%d", i, j), Form("Residual (AllTracks) - Station %d, Module %d;Distance [cm];Entries", i, j), 200, -2, 2);
        }
    }

    h_2d = new TH2D("h_2d","Electron VS Muon angle; Electron [rad]; Muon [rad]" ,500,0.,0.032,500,0.,0.005);
    h_2d_ref = new TH2D("h_2d_ref","Electron VS Muon angle; Electron [rad]; Muon [rad]" ,500,0.,0.032,500,0.,0.005);
}

double FairMUanalyzer::Eevsth(double* x, double* par) {
    double th = *x;
    double Emu = par[0];
    double r = sqrt(Emu*Emu - mu_*mu_) / (Emu + me_);
    return me_ * (1. + r*r * cos(th) * cos(th)) / (1. - r*r * cos(th) * cos(th));
}

double FairMUanalyzer::thmu_vs_the(double* x, double* par) {
    double th = *x;  // in rad
    double Emu = par[0];
    double pmu = sqrt(Emu*Emu - mu_*mu_);
    double Ee = Eevsth(&th, par);
    double pe = sqrt(Ee*Ee - me_*me_);
    return acos((pmu - pe * cos(th)) / sqrt(pe*pe + pmu*pmu - 2 * pmu * pe * cos(th)));
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
    //AnalyzeMF();
    AnalyzeTRK();
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


    TCanvas* c10 = new TCanvas(Form("c10_%s", outputPrefix_.Data()), "Theta_e_theta_mu", 600, 400);
    h_2d->Draw("colz");
    double k = 1; 
    double b = 0.0; 

    for (int i = 1; i <= h_2d_ref->GetNbinsX(); ++i) {
        double x = h_2d_ref->GetXaxis()->GetBinCenter(i);
        double y = k * x + b;
        
        int j = h_2d_ref->GetYaxis()->FindBin(y);
        if (j >= 1 && j <= h_2d_ref->GetNbinsY()) {
            h_2d_ref->SetBinContent(i, j, 1.0); 
        }
    }
    h_2d_ref->Draw("same");
    f_elastic->SetLineWidth(1);
    f_elastic->SetLineColor(kRed);
    f_elastic->Draw("same");

    c10->SaveAs(Form("%s_theta_e_theta_mu.pdf", outputPrefix_.Data()));
    c10->SaveAs(Form("%s_theta_e_theta_mu.root", outputPrefix_.Data()));

    TCanvas* c0 = new TCanvas(Form("c0_%s", outputPrefix_.Data()), "Tracks multiplicity", 600, 400);
    h_Ntracks->GetYaxis()->SetRangeUser(0, h_Ntracks->GetMaximum() * 1.2);
    h_Ntracks->Draw();
    auto legend = new TLegend(0.4, 0.7, 0.9, 0.9);
    legend->AddEntry(h_Ntracks, "(all) tracks multiplicity", "l");
    legend->Draw();
    c0->SaveAs(Form("%s_all_tracksMulti.pdf", outputPrefix_.Data()));

    

    TCanvas* c11 = new TCanvas(Form("c11_%s", outputPrefix_.Data()), "Muon tracks multiplicity", 1200, 400);
    c11->Divide(2,1);
    c11->cd(1);
    h_hits_zcut->Draw();
    c11->cd(2);
    h_isMuon->GetYaxis()->SetRangeUser(0, h_isMuon->GetMaximum() * 1.2);
    h_isMuon->SetLineColor(kRed);
    h_isMuon->Draw("");
    auto legend11 = new TLegend(0.4, 0.7, 0.9, 0.9);
    legend11->AddEntry(h_isMuon, "muon tracks multiplicity by MF", "l");
    legend11->Draw();
    std::cout<<outputPrefix_.Data()<<", isMuon==0: "<<h_isMuon->GetBinContent(1)<<std::endl;
    c11->SaveAs(Form("%s_all_tracksMulti.pdf", outputPrefix_.Data()));

    
    TCanvas* c1 = new TCanvas(Form("c1_%s", outputPrefix_.Data()), "Hit distributions", 2400, 800);
    c1->Divide(4, 2);
    for (int i = 1; i <= 4; i++) {
        c1->cd(i);
        h_hitsPerMuonTrack_zcut[i]->SetLineColor(kRed);
        h_hitsPerMuonTrack_zcut[i]->Draw();
    }
    for (int i = 5; i <= 8; i++) {
        c1->cd(i);
        h_hitsModuleID_zcut[i-4]->SetLineColor(kRed);
        h_hitsModuleID_zcut[i-4]->Draw();
    }
    c1->SaveAs(Form("%s_hit_module.pdf", outputPrefix_.Data()));

    
    TCanvas* c2 = new TCanvas(Form("c2_%s", outputPrefix_.Data()), "Golden muon identified", 1800, 400);
    c2->Divide(3,1);
    for (int i = 1; i <= 3; i++) {
        c2->cd(i);
        h_goldenMuon_isMuon[i-1]->GetYaxis()->SetRangeUser(0, h_goldenMuon_isMuon[i-1]->GetMaximum() * 1.2);
        h_goldenMuon_isMuon[i-1]->Draw();
        std::cout<<outputPrefix_.Data()<<", "<<i<<" station golden muon failed: "<<h_goldenMuon_isMuon[i-1]->GetBinContent(1)<<std::endl;
    }
    c2->SaveAs(Form("%s_golden_muon.pdf", outputPrefix_.Data()));

    TCanvas* c3 = new TCanvas(Form("c3_%s", outputPrefix_.Data()), "Residuals by station", 1800, 400);
    c3->Divide(3, 1);
    for (int i = 0; i < 3; ++i) {
        c3->cd(i + 1);
        h_residual_hitOnTrack[i]->SetLineColor(kRed);
        h_residual_hitOffTrack[i]->SetLineColor(kBlue);
        h_residual_hitOnTrack[i]->Draw();
        h_residual_hitOffTrack[i]->Draw("sames");

        c3->Update();  
        TPaveStats* stats1 = (TPaveStats*)h_residual_hitOnTrack[i]->FindObject("stats");
        TPaveStats* stats2 = (TPaveStats*)h_residual_hitOffTrack[i]->FindObject("stats");

        if (stats1) {
            stats1->SetTextColor(kRed);
            stats1->SetX1NDC(0.70); stats1->SetX2NDC(0.85);
            stats1->SetY1NDC(0.75); stats1->SetY2NDC(0.90);
        }

        if (stats2) {
            stats2->SetTextColor(kBlue);
            stats2->SetX1NDC(0.85); stats2->SetX2NDC(1.00);
            stats2->SetY1NDC(0.75); stats2->SetY2NDC(0.90);
        }

        auto legend2 = new TLegend(0.7, 0.6, 0.9, 0.7);
        legend2->AddEntry(h_residual_hitOnTrack[i], "On Track", "l");
        legend2->AddEntry(h_residual_hitOffTrack[i], "Off Track", "l");
        legend2->Draw();
    }
    c3->SaveAs(Form("%s_residuals_station.pdf", outputPrefix_.Data()));

    TCanvas* c4 = new TCanvas(Form("c4_%s", outputPrefix_.Data()), "Residuals by station/module", 2400, 1200);
    c4->Divide(4, 3);
    for (int i = 0; i < 12; ++i) {
        int ii = i % 4; // module
        int jj = i / 4; // station
        c4->cd(i + 1);
        h_residual_hitOnTrackModule[jj][ii]->SetLineColor(kRed);
        h_residual_hitOffTrackModule[jj][ii]->SetLineColor(kBlue);
        h_residual_hitOnTrackModule[jj][ii]->Draw();
        h_residual_hitOffTrackModule[jj][ii]->Draw("sames");

        c4->Update();  
        TPaveStats* stats1 = (TPaveStats*)h_residual_hitOnTrackModule[jj][ii]->FindObject("stats");
        TPaveStats* stats2 = (TPaveStats*)h_residual_hitOffTrackModule[jj][ii]->FindObject("stats");

        if (stats1) {
            stats1->SetTextColor(kRed);
            stats1->SetX1NDC(0.70); stats1->SetX2NDC(0.85);
            stats1->SetY1NDC(0.75); stats1->SetY2NDC(0.90);
        }

        if (stats2) {
            stats2->SetTextColor(kBlue);
            stats2->SetX1NDC(0.85); stats2->SetX2NDC(1.00);
            stats2->SetY1NDC(0.75); stats2->SetY2NDC(0.90);
        }

        auto legend3 = new TLegend(0.7, 0.6, 0.9, 0.7);
        legend3->AddEntry(h_residual_hitOnTrackModule[jj][ii], "On Track", "l");
        legend3->AddEntry(h_residual_hitOffTrackModule[jj][ii], "Off Track", "l");
        legend3->Draw();
    }
    c4->SaveAs(Form("%s_residuals_station_module.pdf", outputPrefix_.Data()));

    TCanvas* c5 = new TCanvas(Form("c5_%s", outputPrefix_.Data()), "Residuals by station", 2400, 1200);
    c5->Divide(4,3);

    for (int i = 0; i < 12; ++i) {
        c5->cd(i + 1);
        
        int ii = i%4;//module
        int jj = i/4;//station
        h_residual_hitAllTrackModule[jj][ii]->SetName(Form("hAll_%d", i));
        h_residual_hitAllTrackModule[jj][ii]->SetLineColor(kGreen+1);
        h_residual_hitAllTrackModule[jj][ii]->Draw();

        c5->Update();  
        TPaveStats* stats1 = (TPaveStats*)h_residual_hitAllTrackModule[jj][ii]->FindObject("stats");
        if (stats1) {
            stats1->SetTextColor(kBlack);
            stats1->SetX1NDC(0.60); stats1->SetX2NDC(0.85);
            stats1->SetY1NDC(0.65); stats1->SetY2NDC(0.90);
        }

        auto legend4 = new TLegend(0.7, 0.5, 0.9, 0.6);
        legend4->AddEntry(h_residual_hitAllTrackModule[jj][ii], "All Track", "l");
        legend4->Draw();
    }
    c5->SaveAs(Form("%s_alltracks_residuals_station_module.pdf", outputPrefix_.Data()));
}
