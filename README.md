# FairMUanalyzer

A modular C++ analysis framework for MUonE ROOT files using FairRoot infrastructure. Developed by C Zhang.

---

## Requirements

- ROOT (via `root-config`)
- FairSoft (nov22p1): [CVMFS path] `/cvmfs/fairsoft.gsi.de/centos8/fairsoft/nov22p1`
- FairRoot (custom build): [set your path] e.g. `/afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_11July25`

Before building or running, source your environment:
```bash
source /afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_11July25/bin/FairRootConfig.sh
```

---

## Directory Structure

```
FairMUanalyzer/
├── include/                        # Contains FairMUanalyzer.h
├── src/                            # Contains core analysis code
│   ├── FairMUanalyzer.cpp
│   ├── FairMUanalyzer_MF/TRK/OTHERS.cpp
│   └── run_FairMUanalyzer.cpp
├── batch_run.cpp                   # Multi-threaded runner
├── Makefile
└── README.md
```

---

## Build

Compile both single-run and batch-run modes:

```bash
make clean
make -j
```

This produces:
- `run_FairMUanalyzer`: run a single ROOT file
- `batch_run`: multi-threaded batch processor over directory

---

## Usage

### Analyze a single file:
```bash
./run_FairMUanalyzer input.root
```

or you can set the file path inside run_FairMUanalyzer (hard coded), then just run
```bash
./run_FairMUanalyzer 
```

### Analyze all `.root` files in a directory with 4 threads:
```bash
./batch_run root/  4
```
- `root/` is the input directory
- `4` is the number of threads

Each output is saved in:
```
result/<basename>/
```

---

## Output Contents

Each processed file produces:

- `*_output.root`: contains all histograms
- `*_residuals_station.pdf`: comparison of residuals for each station
- `*_residuals_station_module.pdf`: per-module residuals (on/off track)

These are useful for studying tracking residuals and module performance.

---

## Notes

- You can add more `AnalyzeXXX()` modules inside `FairMUanalyzer_SOMETHING.cpp`.
- `.DS_Store`, build artifacts, and result folders are ignored by default via `.gitignore`.

---

## Debugging

If compilation fails:
- Check `Makefile` paths to `FairRoot` and `FairSoft`
- Ensure `MUonERecoOutputAnalysis.h` is found under `include/`
- Ensure ROOT is properly initialized with `root-config`

---

---

## Dependencies

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

## g++
```bash
source /afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_11July25/bin/FairRootConfig.sh
```

```bash
g++ -std=c++17 -O2 -Wall \
  -Iinclude \
  -I/cvmfs/fairsoft.gsi.de/centos8/fairsoft/nov22p1/include \
  -I/afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_11July25/include \
  $(root-config --cflags) \
  src/FairMUanalyzer.cpp \
  src/FairMUanalyzer_Analyze.cpp \
  src/run_FairMUanalyzer.cpp \
  -o run_FairMUanalyzer \
  $(root-config --libs) \
  -L/afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_11July25/lib64 \
  -lMUonEReconstruction -lMUonEReconstructedEventsFilter -lMUonEReconstructionConfiguration \
  -lstdc++fs
```
```bash
g++ -std=c++17 -O2 -Wall -pthread \
  -Iinclude \
  -I/cvmfs/fairsoft.gsi.de/centos8/fairsoft/nov22p1/include \
  -I/afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_11July25/include \
  $(root-config --cflags) \
  src/FairMUanalyzer.cpp \
  src/FairMUanalyzer_Analyze.cpp \
  batch_run.cpp \
  -o batch_run \
  $(root-config --libs) \
  -L/afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_11July25/lib64 \
  -lMUonEReconstruction -lMUonEReconstructedEventsFilter -lMUonEReconstructionConfiguration \
  -lstdc++fs
```


```bash
-std=c++17  支持 std::filesystem 和 modern C++
-I$(root-config --incdir) 引入 ROOT 头文件，例如 TFile.h, TTree.h 等
-I/cvmfs/...  引入 FairSoft 相关头文件
-I/afs/...  引入你的 MUonE 安装头文件
-lMUonEReconstruction 等 你的项目所需的 MUonE .so
-lstdc++fs  GCC <10 时需要手动链接 std::filesystem
$(root-config --libs) 链接所有标准 ROOT 库
```

---

## Drafts

make         # 编译 run_FairMUanalyzer 和 batch_run
运行单个文件：
./run_FairMUanalyzer input.root result/outputname 3
多线程处理整个目录（自动建 result 子目录）：
./batch_run root/ 3 4

如你升级到 GCC ≥ 10，则 std::filesystem 无需 -lstdc++fs，但 CERN 环境多为 8.x/9.x，建议始终保留以确保兼容性。

运行示例见单独文件


