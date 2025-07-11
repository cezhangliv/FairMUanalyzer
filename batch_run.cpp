//make all
//make batch

//运行示例：处理某个路径下的所有 ROOT 文件，用 4 个线程，并设置 MuonFilterHits = 3
//./batch_run root/ 3 4

//or

//root -l -q 'batch_run.C+'


#include "FairMUanalyzer.h"
#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <thread>
#include <mutex>
#include <TSystem.h>

namespace fs = std::filesystem;
std::mutex printMutex;

void processFile(const std::string& infile, int muonFilterHits) {
    TString base = gSystem->BaseName(infile.c_str());  // e.g. input1.root
    base.ReplaceAll(".root", "");                      // → input1
    std::string outdir = "result/" + std::string(base.Data());

    // 创建 result/input1/ 子目录
    std::error_code ec;
    std::filesystem::create_directories(outdir, ec);
    if (ec) {
        std::lock_guard<std::mutex> lock(printMutex);
        std::cerr << "[!] Failed to create directory: " << outdir << " (" << ec.message() << ")" << std::endl;
        return;
    }

    FairMUanalyzer analyzer;
    analyzer.SetInputFile(infile);
    analyzer.SetOutputPrefix(outdir + "/" + base.Data());  // output path with prefix
    analyzer.SetMuonFilterHits(muonFilterHits);
    analyzer.Run();

    std::lock_guard<std::mutex> lock(printMutex);
    std::cout << "[✓] Finished processing: " << infile << std::endl;
}


// ----------------------------
// 1. C++ main entry
// ----------------------------
#ifndef __CINT__
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_dir> [MuonFilterHits=3] [threads=4]" << std::endl;
        return 1;
    }

    std::string inputDir = argv[1];
    int muonFilterHits = (argc > 2) ? std::stoi(argv[2]) : 3;
    int numThreads = (argc > 3) ? std::stoi(argv[3]) : 4;

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
            workers.emplace_back(processFile, files[i], muonFilterHits);
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


