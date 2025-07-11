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
        std::cerr << "Usage: " << argv[0] << " [input.root] [MuonFilterHits] [output_prefix]" << std::endl;
        return 1;
    }

    TString singlefile = (argc > 1) ? argv[1]:"passing_muon_muedaq04-1750192160_MuonFilterHits";
    int muonFilterHits = (argc > 2) ? std::stoi(argv[2]) : 3;
    
    std::string inputFile = Form("root/%s%i_CDbugfix8July25.root",singlefile.Data(), muonFilterHits);
    std::string outputPrefix = (argc > 3) ? argv[3] : Form("result/FairMUanalyzer_%s%i",singlefile.Data(),muonFilterHits);

    FairMUanalyzer analyzer;
    analyzer.SetInputFile(inputFile);
    analyzer.SetOutputPrefix(outputPrefix);
    analyzer.SetMuonFilterHits(muonFilterHits);
    analyzer.Run();

    return 0;
}
