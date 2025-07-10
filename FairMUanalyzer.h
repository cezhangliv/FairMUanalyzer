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

private:
    std::string inputFilePath_;
    TString outputPrefix_;
    int MuonFilterHits_;

    TFile* inputFile_;
    TTree* cbmsim_;
    MUonERecoOutputAnalysis* reco_;

    // Histograms
    TH1F* h_residual_hitOnTrack[3];
    TH1F* h_residual_hitOffTrack[3];
    TH1F* h_residual_hitOnTrackModule[3][4];
    TH1F* h_residual_hitOffTrackModule[3][4];
    TH1F* h_residual_hitAllTrackModule[3][4];

    void Init();
    void Analyze();
    void SaveResults();

    TVector3 getXYfromHitMF(const MUonERecoOutputHitAnalysis& hit);
    double computeSigned2DResidualMF(const TVector3& p3D, const TVector3& x0, const TVector3& h, int moduleID);
};

#endif
