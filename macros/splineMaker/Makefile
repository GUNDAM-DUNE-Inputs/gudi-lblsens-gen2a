CXX      := c++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

YAML_CFLAGS := $(shell pkg-config --cflags yaml-cpp)
YAML_LIBS   := $(shell pkg-config --libs yaml-cpp)

TARGET := MakeSplines
SRC    := MakeFDESSplines.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(YAML_CFLAGS) $< \
	       $(ROOTLIBS) $(YAML_LIBS) -o $@

clean:
	rm -f $(TARGET) *.o *.so *.pcm *.d *.rootmap

.PHONY: all clean

