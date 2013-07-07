CXX = g++

CXXFLAGS = -O3 -Wall -DNDEBUG -Wno-deprecated -pedantic -finline-functions -foptimize-sibling-calls -Wcast-qual -Wwrite-strings -Wsign-promo -Wcast-align -Wno-long-long -fexpensive-optimizations -funroll-all-loops -ffast-math -fomit-frame-pointer -pipe

%CXXFLAGS = -O3 -Wall -DNDEBUG -Wno-deprecated -pedantic -finline-functions -foptimize-sibling-calls -Wcast-qual -Wwrite-strings -Wsign-promo -Wcast-align -Wno-long-long -fexpensive-optimizations -funroll-all-loops -ffast-math -fomit-frame-pointer -pipe -std=c++0x
%CXXFLAGS = -O3 -Wall -g -DNDEBUG -Wno-deprecated -pedantic -ansi -finline-functions -foptimize-sibling-calls -Wcast-qual -Wwrite-strings -Wsign-promo -Wcast-align -Wno-long-long -fexpensive-optimizations -funroll-all-loops -ffast-math -pipe -I./ -L./ -std=c++0x
%CXXFLAGS = -O3 -DNDEBUG -W -Wall -Wno-deprecated -pedantic -ansi -finline-functions -foptimize-sibling-calls -Wcast-qual -Wwrite-strings -Wsign-promo -Wcast-align -Wno-long-long -fexpensive-optimizations -funroll-all-loops -ffast-math -fomit-frame-pointer -pipe -I./ -L./ -std=c++0x
%CXXFLAGS = -g -W -Wall -Wno-deprecated

LINKFLAGS = -lm
OTHERFLAGS = -pg

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