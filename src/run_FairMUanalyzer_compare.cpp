#include "FairMUanalyzer.h"
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <string>


TChain* LoadChain(const std::string& input) {
    TChain* chain = new TChain("cbmsim");

    if (input.size() > 4 && input.substr(input.size() - 4) == ".txt") {
        std::ifstream infile(input);
        std::string line;
        while (std::getline(infile, line)) {
            if (!line.empty()) {
                chain->Add(line.c_str());
            }
        }
    } else {
        chain->Add(input.c_str());
    }

    return chain;
}


bool CompareEvents(FairMUanalyzer& ana1, FairMUanalyzer& ana2, Long64_t i) {
    bool g1 = ana1.IsGoldenEvent();
    bool g2 = ana2.IsGoldenEvent();

    if (g1 != g2) {
        std::cout << "Event " << i << " mismatch in golden status: "
                  << g1 << " vs " << g2 << std::endl;
        return false;
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: ./run_FairMUanalyzer_compare input1.root|.txt input2.root|.txt" << std::endl;
        return 1;
    }

    TChain* chain1 = LoadChain(argv[1]);
    TChain* chain2 = LoadChain(argv[2]);
    std::cout<<chain1->GetEntries()<<" "<<chain2->GetEntries()<<std::endl;

    FairMUanalyzer analyzer1, analyzer2;

    
    analyzer1.SetInputTree(chain1);
    analyzer2.SetInputTree(chain2);

    //analyzer.SetTgt(1);
    //analyzer.SetMf();


    Long64_t nEntries = std::min(chain1->GetEntries(), chain2->GetEntries());


    
    int mismatch = 0;
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain1->GetEntry(i);
        chain2->GetEntry(i);

        analyzer1.ProcessEvent(i);
        analyzer2.ProcessEvent(i);

        if(!CompareEvents(analyzer1, analyzer2, i)){
            mismatch++;
            if(!analyzer1.IsGoldenEvent())continue;
            std::cout<<"HS0"<<std::endl;
            analyzer1.ProcessEvent(i,1);
            std::cout<<"HS2"<<std::endl;
            analyzer2.ProcessEvent(i,1);
        }
        if(mismatch>1)break;
    }
    

    return 0;
}
