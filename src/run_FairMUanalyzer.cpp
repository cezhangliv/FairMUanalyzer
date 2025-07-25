#include "FairMUanalyzer.h"
#include <iostream>
#include <string>
#include <TSystem.h>

//./run_FairMUanalyzer [input.root] [output_prefix] [RunN] [Tgt] 
// more see run_FairMUanalyzer.sh

int main(int argc, char** argv) {
    
    std::string inputFile;
    if (argc < 2) {
        
        inputFile = "root/single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root";
        //std::cerr << "Usage: " << argv[0] << " input.root [output_prefix]" << std::endl;
        //return 1;
    }
    else inputFile = argv[1];

    
    TString base = gSystem->BaseName(inputFile.c_str());
    base.ReplaceAll(".root", "");
    std::string outputPrefix = (argc > 2) ? argv[2] : "result/FairMUanalyzer_" + std::string(base.Data());

    std::cout<<"inputFile: "<<inputFile<<std::endl;
    std::cout<<"outputPrefix: "<<outputPrefix<<std::endl;


    FairMUanalyzer analyzer;
    analyzer.SetInputFile(inputFile);
    analyzer.SetOutputPrefix(outputPrefix);
    
    //// configurations

    analyzer.SetSavepdf(true);
    analyzer.SetRunN((argc > 3) ? std::stoll(argv[3]) : -1);
	analyzer.SetTgt((argc > 4) ? std::stoi(argv[4]) : 1);
    
    bool mf_flag = true;
	if (argc > 5) {
    	std::string arg = argv[5];
    	std::transform(arg.begin(), arg.end(), arg.begin(), ::tolower);
    	mf_flag = (arg == "true" || arg == "1");
	}
	analyzer.SetMf(mf_flag);
    analyzer.SetUseTightTrackCut(false,1);//default is true: only for tgt1 Ntrack==4, otherwise Ntrack>=4 for golden muon first step selection
    //analyzer.SetMuonFilterHits(3);

    ////////////////////

    analyzer.Run();

    return 0;
}
