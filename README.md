# FairMUanalyzer

A C++ analysis tool designed to process MUonE `cbmsim` ROOT files and extract hit residuals from golden muon tracks. Built for use within the MUonE environment (FairRoot + MUonEReco* classes), it supports both single-file and multi-file batch processing with multithreading.

---

## 🔧 Environment Setup

You must have the MUonE and FairRoot environments loaded. Typically:

```bash
export SIMPATH=/cvmfs/fairsoft.gsi.de/centos8/fairsoft/nov22p1
source /afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_8July25/bin/FairRootConfig.sh
```

Ensure `root-config` and MUonE shared libraries are available.

---

## ⚙️ Build Instructions

To build everything:

```bash
make
```

This creates:

- `run_FairMUanalyzer`: for analyzing one file at a time
- `batch_run`: for multithreaded batch processing over many `.root` files in a directory

If needed, clean with:

```bash
make clean
```

---

## 🚀 Usage

### 1. Single-file analysis

```bash
./run_FairMUanalyzer input.root result/output_name 3
```

- `input.root`: path to ROOT file with `cbmsim` TTree
- `result/output_name`: output prefix, creates `output_name_output.root` and `.pdf` in `result/`
- `3`: (optional) required number of hits to identify muon (default: 3)

### 2. Batch processing with multithreading

```bash
./batch_run root_directory 3 4
```

- `root_directory`: folder containing multiple `.root` files
- `3`: required number of hits for muon
- `4`: number of threads to use

Outputs will be saved under `result/<basename>/`.

---

## 📁 Output Contents

Each processed file produces:

- `*_output.root`: contains all histograms
- `*_residuals_station.pdf`: comparison of residuals for each station
- `*_residuals_station_module.pdf`: per-module residuals (on/off track)

These are useful for studying tracking residuals and module performance.

---

## 📚 Class Overview: `FairMUanalyzer`

Key methods:

```cpp
void SetInputFile(const std::string& path);
void SetOutputPrefix(const std::string& prefix);
void SetMuonFilterHits(int val);  // default = 3
void Run();  // calls Init(), Analyze(), SaveResults()
```

---

## 🧪 Dependencies

- ROOT (via `root-config`)
- MUonE reconstruction libraries:
  - `libMUonEReconstruction.so`
  - `libMUonEReconstructedEventsFilter.so`
  - `libMUonEReconstructionConfiguration.so`
- C++17 (`std::filesystem` and threads)

Linking requires:

```bash
-lstdc++fs  # for GCC < 10
```

---

## Drafts

make         # 编译 run_FairMUanalyzer 和 batch_run
运行单个文件：
./run_FairMUanalyzer input.root result/outputname 3
多线程处理整个目录（自动建 result 子目录）：
./batch_run root/ 3 4

如你升级到 GCC ≥ 10，则 std::filesystem 无需 -lstdc++fs，但 CERN 环境多为 8.x/9.x，建议始终保留以确保兼容性。

g++ -O2 -Wall -std=c++17 -pthread \
  -I$(root-config --incdir) \
  -I/cvmfs/fairsoft.gsi.de/centos8/fairsoft/nov22p1/include \
  -I/afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_8July25/include \
  FairMUanalyzer.cpp batch_run.cpp \
  -L/afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_8July25/lib64 \
  -lMUonEReconstruction -lMUonEReconstructedEventsFilter -lMUonEReconstructionConfiguration \
  $(root-config --libs) \
  -lstdc++fs \
  -o batch_run


g++ -O2 -Wall -std=c++17 \
  -I$(root-config --incdir) \
  -I/cvmfs/fairsoft.gsi.de/centos8/fairsoft/nov22p1/include \
  -I/afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_8July25/include \
  FairMUanalyzer.cpp run_FairMUanalyzer.cpp \
  -L/afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_8July25/lib64 \
  -lMUonEReconstruction -lMUonEReconstructedEventsFilter -lMUonEReconstructionConfiguration \
  $(root-config --libs) \
  -lstdc++fs \
  -o run_FairMUanalyzer

-std=c++17  支持 std::filesystem 和 modern C++
-I$(root-config --incdir) 引入 ROOT 头文件，例如 TFile.h, TTree.h 等
-I/cvmfs/...  引入 FairSoft 相关头文件
-I/afs/...  引入你的 MUonE 安装头文件
-lMUonEReconstruction 等 你的项目所需的 MUonE .so
-lstdc++fs  GCC <10 时需要手动链接 std::filesystem
$(root-config --libs) 链接所有标准 ROOT 库

运行示例

./batch_run root/ 3 4   # 目录 + muon filter hits + 线程数
./run_FairMUanalyzer root/single.root result/output 3




## 📝 License

Internal use only — maintained by Ce Zhang at Liverpool / CERN.
