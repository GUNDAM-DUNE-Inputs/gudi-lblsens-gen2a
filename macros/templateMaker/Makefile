# ==========================================================
# Project Settings
# ==========================================================

TARGET      := analysis
CXX         := g++
STD         := -std=c++17

# ==========================================================
# Build Type (default: release)
#   make BUILD=debug
# ==========================================================

BUILD ?= release

ifeq ($(BUILD),debug)
    OPTFLAGS := -O0 -g -DDEBUG
else
    OPTFLAGS := -O3 -march=native -flto
endif

WARNFLAGS   := -Wall -Wextra -Wpedantic

# ==========================================================
# ROOT Configuration
# ==========================================================

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLIBS    := $(shell root-config --libs)

# ==========================================================
# yaml-cpp Configuration
# ==========================================================

YAMLCFLAGS  := $(shell pkg-config --cflags yaml-cpp)
YAMLLIBS    := $(shell pkg-config --libs yaml-cpp)

# ==========================================================
# Compiler Flags
# ==========================================================

CXXFLAGS := $(STD) $(OPTFLAGS) $(WARNFLAGS) \
            $(ROOTCFLAGS) $(YAMLCFLAGS)

LDFLAGS  := $(ROOTLIBS) $(YAMLLIBS)

# ==========================================================
# Source & Object Files
# ==========================================================

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(SRCS:.cpp=.d)

# ==========================================================
# Default Target
# ==========================================================

all: $(TARGET)

# ==========================================================
# Link Step
# ==========================================================

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) -o $@

# ==========================================================
# Compile Step (with auto-dependency generation)
# ==========================================================

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# ==========================================================
# Clean
# ==========================================================

clean:
	rm -f *.o *.d $(TARGET)

