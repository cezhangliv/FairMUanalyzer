#include "FairMUanalyzer.h"
#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <TSystem.h>

//./batch_run root/prefix

namespace fs = std::filesystem;

void processFile(const std::string& infile) {
    TString base = gSystem->BaseName(infile.c_str());
    base.ReplaceAll(".root", "");
    std::string outdir = "result/" + std::string(base.Data());

    std::error_code ec;
    fs::create_directories(outdir, ec);
    if (ec) {
        std::cerr << "! Failed to create directory: " << outdir << " (" << ec.message() << ")" << std::endl;
        return;
    }

    FairMUanalyzer analyzer;
    analyzer.SetInputFile(infile);
    analyzer.SetOutputPrefix(outdir + "/FairMUanalyzer_" + std::string(base.Data()));
    
    analyzer.SetSavepdf(true);
    analyzer.SetRunN(-1);
    analyzer.SetTgt(0);//0/1
    analyzer.SetMf(true);

    //analyzer.SetUseTightTrackCut(false,1);//default is true for tgt1(1 in braket) and flase for tgt0(0 in braket): only for tgt1 Ntrack==4, otherwise Ntrack>=4 for golden muon first step selection
    //analyzer.SetMuonFilterHits(3);

    analyzer.Run();

    std::cout << "Finished processing: " << infile << std::endl;
    std::cout << "-----" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_dir_or_prefix>" << std::endl;
        return 1;
    }

    std::string inputArg = argv[1];
    std::string searchDir = inputArg;
    std::string prefix;

    if (!fs::is_directory(inputArg)) {
        fs::path full(inputArg);
        searchDir = full.parent_path().string();
        prefix = full.filename().string();
    }

    std::vector<std::string> files;
    for (const auto& entry : fs::directory_iterator(searchDir)) {
        if (entry.path().extension() == ".root") {
            std::string filename = entry.path().filename().string();
            if (prefix.empty() || filename.rfind(prefix, 0) == 0) {
                files.push_back(entry.path().string());
            }
        }
    }

    std::cout << "Found " << files.size() << " ROOT files in " << inputArg << std::endl;

    for (const auto& file : files) {
        processFile(file);  
    }

    return 0;
}
