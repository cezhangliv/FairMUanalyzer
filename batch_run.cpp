#include "FairMUanalyzer.h"
#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <thread>
#include <mutex>
#include <TSystem.h>

//./batch_run <input_dir> [num_threads=4]
//./batch_run root/
//./batch_run root/ 8

namespace fs = std::filesystem;
std::mutex printMutex;

void processFile(const std::string& infile) {
    TString base = gSystem->BaseName(infile.c_str());
    base.ReplaceAll(".root", "");
    std::string outdir = "result/" + std::string(base.Data());

    // 创建输出目录
    std::error_code ec;
    std::filesystem::create_directories(outdir, ec);
    if (ec) {
        std::lock_guard<std::mutex> lock(printMutex);
        std::cerr << "[!] Failed to create directory: " << outdir << " (" << ec.message() << ")" << std::endl;
        return;
    }

    FairMUanalyzer analyzer;
    analyzer.SetInputFile(infile);
    analyzer.SetOutputPrefix(outdir + "/FairMUanalyzer_" + std::string(base.Data()));
    analyzer.Run();

    std::lock_guard<std::mutex> lock(printMutex);
    std::cout << "[✓] Finished processing: " << infile << std::endl;
}

#ifndef __CINT__
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_dir> [threads=4]" << std::endl;
        return 1;
    }

    std::string inputDir = argv[1];
    int numThreads = (argc > 2) ? std::stoi(argv[2]) : 4;

    std::vector<std::string> files;
    for (const auto& entry : fs::directory_iterator(inputDir)) {
        if (entry.path().extension() == ".root") {
            files.push_back(entry.path().string());
        }
    }

    std::cout << "Found " << files.size() << " ROOT files in " << inputDir << std::endl;

    std::vector<std::thread> workers;
    int i = 0;
    while (i < files.size()) {
        while (workers.size() < numThreads && i < files.size()) {
            workers.emplace_back(processFile, files[i]);
            ++i;
        }
        for (auto& t : workers) {
            if (t.joinable()) t.join();
        }
        workers.clear();
    }

    return 0;
}
#endif
