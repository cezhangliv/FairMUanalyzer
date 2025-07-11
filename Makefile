# === FairSoft and FairRoot paths ===
FAIRSOFT_BASE  = /cvmfs/fairsoft.gsi.de/centos8/fairsoft/nov22p1
FAIRROOT_BASE  = /afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_8July25

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

# === Sources and Objects ===
SRCS = $(wildcard src/*.cpp)
OBJS = $(SRCS:src/%.cpp=build/%.o)

BATCH_OBJS = build/FairMUanalyzer.o build/FairMUanalyzer_MF.o build/batch_run.o

TARGET = run_FairMUanalyzer
BATCH  = batch_run

# === Targets ===
all: $(TARGET) $(BATCH)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BATCH): $(BATCH_OBJS)
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
	rm -rf build $(TARGET) $(BATCH)
