#ifndef FAIRMUANALYZER_H
#define FAIRMUANALYZER_H

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "MUonERecoOutputAnalysis.h"


#include <vector>
#include <string>

class FairMUanalyzer {
public:
    FairMUanalyzer();
    ~FairMUanalyzer();

    void SetInputFile(const std::string& path);
    void SetOutputPrefix(const std::string& prefix);
    void SetMuonFilterHits(int val);
    void Run();

    static double Eevsth(double* x, double* par);
    static double thmu_vs_the(double* x, double* par);


private:
    std::string inputFilePath_;
    TString outputPrefix_;
    int MuonFilterHits_;
    int goldenevents_;
    bool savepdf_;

    TFile* inputFile_;
    TTree* cbmsim_;
    MUonERecoOutputAnalysis* reco_;

    // Histograms
    TH1F* h_residual_hitOnTrack[3];
    TH1F* h_residual_hitOffTrack[3];
    TH1F* h_residual_hitOnTrackModule[3][4];
    TH1F* h_residual_hitOffTrackModule[3][4];
    TH1F* h_residual_hitAllTrackModule[3][4];

    TH1F* h_hits_zcut;
    TH1F* h_hitsModuleID_zcut[5];
    TH1F* h_hitsPerMuonTrack_zcut[5];
    TH1F* h_isMuon;
    TH1F* h_Ntracks;
    TH1F* h_goldenMuon_isMuon[3];
    TH2D *h_2d;
    TH2D *h_2d_ref;
    TF1 *f_elastic;

    void Init();
    void Analyze();
    void AnalyzeTRK();
    void AnalyzeMF();
    void SaveResults();

    static constexpr double me_ = 0.51099906e-3;   // Electron mass [GeV]
    static constexpr double mu_ = 105.65836900e-3; // Muon mass [GeV]
    static constexpr double Ebeam_ = 160;          // Beam energy [GeV]




    TVector3 getXYfromHitMF(const MUonERecoOutputHitAnalysis& hit);
    double computeSigned2DResidualMF(const TVector3& p3D, const TVector3& x0, const TVector3& h, int moduleID);
    
};

#endif
