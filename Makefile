# === FairSoft and FairRoot paths ===
FAIRSOFT_BASE  = /cvmfs/fairsoft.gsi.de/centos8/fairsoft/nov22p1
#FAIRROOT_BASE  = /afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_11July25
FAIRROOT_BASE  = /afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_09Oct25
#FAIRROOT_BASE  = /afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_21Nov25_v0176

# === Compiler and Flags ===
CXX = g++
CXXFLAGS = -O2 -Wall -std=c++17 $(shell root-config --cflags) \
           -I$(FAIRSOFT_BASE)/include \
           -I$(FAIRROOT_BASE)/include \
           -Iinclude
CXXFLAGS += -I$(shell root-config --incdir)

LDFLAGS = $(shell root-config --libs) \
  -L$(FAIRROOT_BASE)/lib64 \
  -lMUonEReconstruction -lMUonEReconstructedEventsFilter -lMUonEReconstructionConfiguration \
  -lstdc++fs

# === Source Files ===
FAIRMU_SRC = \
  build/FairMUanalyzer.o \
  build/FairMUanalyzer_ProcessEvent.o \
  build/FairMUanalyzer_MF.o \
  build/FairMUanalyzer_TRK.o

# === Targets ===
TARGETS = run_FairMUanalyzer run_FairMUanalyzer_compare batch_run

all: $(TARGETS)

run_FairMUanalyzer: $(FAIRMU_SRC) build/run_FairMUanalyzer.o
	$(CXX) -o $@ $^ $(LDFLAGS)

run_FairMUanalyzer_compare: $(FAIRMU_SRC) build/run_FairMUanalyzer_compare.o
	$(CXX) -o $@ $^ $(LDFLAGS)

batch_run: $(FAIRMU_SRC) build/batch_run.o
	$(CXX) -o $@ $^ $(LDFLAGS) -pthread

# === Compilation Rules ===
build/%.o: src/%.cpp include/FairMUanalyzer.h | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/batch_run.o: batch_run.cpp include/FairMUanalyzer.h | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

# === Create build directory if missing ===
build:
	mkdir -p build

# === Clean ===
clean:
	rm -rf build $(TARGETS)
