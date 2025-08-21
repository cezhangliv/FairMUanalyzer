// draw_thmu_vs_the.C
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include <cmath>

static constexpr double Ebeam_ = 4.0; // Beam energy [GeV]
static constexpr double mu_ = 0.1056583745; // muon mass [GeV/c^2]
static constexpr double me_ = 0.0005109989461; // electron mass [GeV/c^2]

double Eevsth(double* x, double* par) {
    double th = *x;     // electron scattering angle [rad]
    double Emu = par[0];
    double r = sqrt(Emu*Emu - mu_*mu_) / (Emu + me_);
    return me_ * (1. + r*r * cos(th) * cos(th)) /
                 (1. - r*r * cos(th) * cos(th));
}

double thmu_vs_the(double* x, double* par) {
    double th = *x;  // electron scattering angle [rad]
    double Emu = par[0];
    double pmu = sqrt(Emu*Emu - mu_*mu_);
    double Ee = Eevsth(&th, par);
    double pe = sqrt(Ee*Ee - me_*me_);
    return acos((pmu - pe * cos(th)) /
                sqrt(pe*pe + pmu*pmu - 2 * pmu * pe * cos(th)));
}

void draw_thmu_vs_the2() {
    
    //TF1* f_elastic = new TF1("f_elastic", thmu_vs_the, 0., 0.032, 1);
    TF1* f_elastic = new TF1("f_elastic", thmu_vs_the, 0., 0.1, 1);
    f_elastic->SetParameter(0, Ebeam_);
    f_elastic->SetNpx(100000);
    f_elastic->SetTitle("Theoretical muon scattering angle vs electron angle;Electron angle [rad];Muon angle [rad]");

    TCanvas* c1 = new TCanvas("c1", "Muon vs Electron scattering angle", 800, 600);
    f_elastic->Draw();
}

void draw_thmu_vs_the() {
    
    TCanvas* c1 = new TCanvas("c1", "Muon vs Electron scattering angle", 800, 600);

    
    double energies[3] = {2.0, 4.0, 10.0}; // GeV
    int colors[3] = {kRed, kBlue, kGreen+2};
    TF1* funcs[3];

    
    for (int i = 0; i < 3; i++) {
        funcs[i] = new TF1(Form("f_elastic_%.0fGeV", energies[i]), thmu_vs_the, 0., 0.1, 1);
        //funcs[i] = new TF1(Form("f_elastic_%.0fGeV", energies[i]), thmu_vs_the, 0., 1.0, 1);
        funcs[i]->SetParameter(0, energies[i]);
        funcs[i]->SetNpx(100000);
        funcs[i]->SetLineColor(colors[i]);
        funcs[i]->SetLineWidth(2);
        funcs[i]->SetTitle("Theoretical muon scattering angle vs electron angle;Electron angle [rad];Muon angle [rad]");
        if (i == 0)
            funcs[i]->Draw();
        else
            funcs[i]->Draw("SAME");
    }

    TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    for (int i = 0; i < 3; i++) {
        leg->AddEntry(funcs[i], Form("E_{beam} = %.0f GeV", energies[i]), "l");
    }
    leg->Draw();
}
