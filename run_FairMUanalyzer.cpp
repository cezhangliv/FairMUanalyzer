#include "FairMUanalyzer.h"
#include <iostream>

//g++ run_FairMUanalyzer.cpp FairMUanalyzer.cpp $(root-config --cflags --libs) -o run_FairMUanalyzer

//g++ -O2 -Wall -std=c++17 \
//    -I/cvmfs/fairsoft.gsi.de/centos8/fairsoft/nov22p1/include \
//    -I/afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_8July25/include \
//    FairMUanalyzer.cpp run_FairMUanalyzer.cpp \
//    -L/afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_8July25/lib64 \
//    -lMUonERecoOutput \
//    $(root-config --libs) \
//    -o run_FairMUanalyzer


//./run_FairMUanalyzer root/myfile.root result/outputName 3


int main(int argc, char** argv) {
    
    if (argc < 1) {
        std::cerr << "Usage: " << argv[0] << " input.root [output_prefix] [MuonFilterHits]" << std::endl;
        return 1;
    }

    int MuonFilterHits=3;
    TString singlefile = "passing_muon_muedaq04-1750192160_MuonFilterHits";

    std::string inputFile = (argc > 2) ? argv[1]: 
        Form("root/%s%i_CDbugfix8July25.root",singlefile.Data(), MuonFilterHits);
    std::string outputPrefix = (argc > 2) ? argv[2] : Form("result/FairMUanalyzer_%s",singlefile.Data());
    int muonFilterHits = (argc > 3) ? std::stoi(argv[3]) : 3;

    FairMUanalyzer analyzer;
    analyzer.SetInputFile(inputFile);
    analyzer.SetOutputPrefix(outputPrefix);
    analyzer.SetMuonFilterHits(muonFilterHits);
    analyzer.Run();

    return 0;
}
