CXX      := c++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

YAMLLIBS := -lyaml-cpp

TARGET := Converter
SRC    := Converter.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $< $(ROOTLIBS) $(YAMLLIBS) -o $@


clean:
	rm -f $(TARGET) *.o *.so *.pcm *.d *.rootmap

.PHONY: all clean

