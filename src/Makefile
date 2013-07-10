CXX = g++
CXXFLAGS = -O3 -Wall -DNDEBUG -Wno-deprecated -pedantic -finline-functions -foptimize-sibling-calls -Wcast-qual -Wwrite-strings -Wsign-promo -Wcast-align -Wno-long-long -fexpensive-optimizations -funroll-all-loops -ffast-math -fomit-frame-pointer -pipe
LINKFLAGS = -lm

ifeq ($(OS),Windows_NT)
    CCFLAGS += -D WIN32
    ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
        CCFLAGS += -std=c++0x -D AMD64
    endif
    ifeq ($(PROCESSOR_ARCHITECTURE),x86)
        CCFLAGS += -std=c++0x -D IA32
    endif
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        CXXFLAGS += -std=c++0x
    endif
    ifeq ($(UNAME_S),Darwin)

    endif
    UNAME_P := $(shell uname -p)
    ifeq ($(UNAME_P),x86_64)
        %CXXFLAGS += -D AMD64
    endif
    ifneq ($(filter %86,$(UNAME_P)),)
        %CXXFLAGS += -D IA32
    endif
    ifneq ($(filter arm%,$(UNAME_P)),)
        %CXXFLAGS += -D ARM
    endif
endif

SRCS1 = \
	bit_array.cpp \
	wat_array.cpp \
	fmIndex.cpp \

SRCS2 = \
	randal_index.cpp \

SRCS3 = \
	randal_search.cpp \
	randal_align.cpp \

OBJS1 = $(SRCS1:%.cpp=%.o)
OBJS2 = $(SRCS2:%.cpp=%.o)
OBJS3 = $(SRCS3:%.cpp=%.o)

all: index align

index: $(OBJS1) $(OBJS2)
	$(CXX) $(CXXFLAGS) $(LINKFLAGS) $(OBJS1) $(OBJS2) -o randal-index

align: $(OBJS1) $(OBJS3)
	$(CXX) $(CXXFLAGS) $(LINKFLAGS) $(OBJS1) $(OBJS3) -o randal-align

debug:
	make all CXXFLAGS="-ggdb -W -Wall -pedantic"

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f randal-index randal-align *.o *~