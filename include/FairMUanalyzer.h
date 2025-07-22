#ifndef FAIRMUANALYZER_H
#define FAIRMUANALYZER_H

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TGraph.h"
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
    
    void SetSavepdf(bool val);
    void SetMf(bool val);
    void SetRunN(Long64_t val);
    void SetTgt(int val);

    void SetMuonFilterHits(int val);//not very useful
    
    void Run();

    static double Eevsth(double* x, double* par);
    static double thmu_vs_the(double* x, double* par);


private:
    std::string inputFilePath_;
    TString outputPrefix_;
    
    bool mf_;
    bool savepdf_;
    Long64_t runN_;
    int tgt_;
    int MuonFilterHits_;

    int goldenevents_;

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

    TH1F* h_goldenMuon_isMuon[3];

    TH1F* h_isMuon;
    TH1F* h_Ntracks;
    
    TH2D *h_2d;
    TH2D *h_2d_ref;
    TF1 *f_elastic;
    TGraph* g_elastic;

    TH1I* hCaseDist;
    std::map<std::string, int> case_counts;


    void Init();
    void Analyze();
    void AnalyzeTRK();
    void AnalyzeMF();
    void SaveResults();

    static constexpr double me_ = 0.51099906e-3;   // Electron mass [GeV]
    static constexpr double mu_ = 105.65836900e-3; // Muon mass [GeV]
    static constexpr double Ebeam_ = 160;          // Beam energy [GeV]

    static constexpr double z_tgt1_ = 664.6;          // origin: 667.3 ; targetRelativePosition: -2.7
    static constexpr double z_tgt2_ = 780.7;        //cm,  origin: 784.6 ; targetRelativePosition: -3.9

    static constexpr int maxNhitInStat_ = 30; // Giovanni A suggested to set a max hit cut in a station

    double intersecX_ = -1;


    TVector3 getXYfromHitMF(const MUonERecoOutputHitAnalysis& hit);
    double computeSigned2DResidualMF(const TVector3& p3D, const TVector3& x0, const TVector3& h, int moduleID);

    double acoplanarity(const TVector3 in, const TVector3 out1, const TVector3 out2);
    
};

#endif
