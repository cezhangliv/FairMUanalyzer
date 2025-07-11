#include "FairMUanalyzer.h"
#include <iostream>
#include <TSystem.h>

//./run_FairMUanalyzer <input.root> [output_prefix]
// ./run_FairMUanalyzer root/event123.root
// output: result/FairMUanalyzer_event123.root
//./run_FairMUanalyzer root/event123.root custom/output/path/myresult
// output to custom/output/path/myresult.root

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input.root [output_prefix]" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    TString base = gSystem->BaseName(inputFile.c_str());
    base.ReplaceAll(".root", "");
    std::string outputPrefix = (argc > 2) ? argv[2] : "result/FairMUanalyzer_" + std::string(base.Data());

    FairMUanalyzer analyzer;
    analyzer.SetInputFile(inputFile);
    analyzer.SetOutputPrefix(outputPrefix);
    analyzer.Run();

    return 0;
}
