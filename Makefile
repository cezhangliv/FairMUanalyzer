# === FairSoft and FairRoot paths ===
FAIRSOFT_BASE  = /cvmfs/fairsoft.gsi.de/centos8/fairsoft/nov22p1
FAIRROOT_BASE  = /afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_8July25

# === Compiler and Flags ===
CXX = g++
CXXFLAGS = -O2 -Wall -std=c++17 $(shell root-config --cflags) \
           -I$(FAIRSOFT_BASE)/include \
           -I$(FAIRROOT_BASE)/include
CXXFLAGS += -I$(shell root-config --incdir)

LDFLAGS = $(shell root-config --libs) \
  -L/afs/cern.ch/user/c/cez/eos/Soft/fair_install/FairRoot/install_8July25/lib64 \
  -lMUonEReconstruction -lMUonEReconstructedEventsFilter -lMUonEReconstructionConfiguration \
  -lstdc++fs


# === File Lists ===
HEADERS = FairMUanalyzer.h
SRCS    = FairMUanalyzer.cpp run_FairMUanalyzer.cpp
BATCH_SRCS = FairMUanalyzer.cpp batch_run.cpp

OBJS       = $(SRCS:.cpp=.o)
BATCH_OBJS = $(BATCH_SRCS:.cpp=.o)

TARGET = run_FairMUanalyzer
BATCH  = batch_run

# === Default build ===
all: $(TARGET) $(BATCH)

# === Main Analyzer Executable ===
$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

# === Batch Mode Executable (multi-threaded) ===
$(BATCH): $(BATCH_OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS) -pthread

# === Generic Compilation Rule ===
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# === Clean Target ===
clean:
	rm -f *.o $(TARGET) $(BATCH)
