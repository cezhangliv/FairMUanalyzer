
# === FairSoft and FairRoot paths ===
FAIRSOFT_BASE  = /cvmfs/fairsoft.gsi.de/centos8/fairsoft/nov22p1
FAIRROOT_BASE  = /afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_11July25

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

# === Targets ===
TARGETS = run_FairMUanalyzer run_FairMUanalyzer_compare batch_run

all: $(TARGETS)

run_FairMUanalyzer: build/FairMUanalyzer.o build/FairMUanalyzer_MF.o build/FairMUanalyzer_TRK.o build/run_FairMUanalyzer.o
	$(CXX) -o $@ $^ $(LDFLAGS)

run_FairMUanalyzer_compare: build/FairMUanalyzer.o build/FairMUanalyzer_MF.o build/FairMUanalyzer_TRK.o build/run_FairMUanalyzer_compare.o
	$(CXX) -o $@ $^ $(LDFLAGS)

batch_run: build/FairMUanalyzer.o build/FairMUanalyzer_MF.o build/FairMUanalyzer_TRK.o build/batch_run.o
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
