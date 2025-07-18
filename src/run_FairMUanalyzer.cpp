#include "FairMUanalyzer.h"
#include <iostream>
#include <TSystem.h>

//./run_FairMUanalyzer [input.root] [output_prefix]
//./run_FairMUanalyzer root/event123.root custom/output/path/myresult
// output to custom/output/path/myresult.root



int main(int argc, char** argv) {
    
    std::string inputFile;
    if (argc < 2) {
        
        //inputFile = "root/single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1.root";
        //inputFile = "root/single_muon_interaction_1_CDbugfix11July25_run11_muedaq04-1750493893-1750494687_MF1.root";
        //inputFile = "root/single_muon_interaction_1_CDbugfix11July25_run12_muedaq04-1752206496-1752209387_MF1.root";
		inputFile = "root/single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_maxNumberOfSharedHits2.root";
        
        
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
    analyzer.Run();

    return 0;
}
