DIRS := . ../src

SOURCES := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.cpp))
OBJS := $(patsubst %.cpp, %.o, $(SOURCES))

CFLAGS := -std=c++17 -ltensorflow -lpthread -no-pie -fmax-errors=5 -I ~/libtensorflow2/include/ -I libs/cppflow-master/include/
CXX ?= g++
LIBS :=
INCLUDES :=
LIBDIR :=  libs/blossom5-v2.05.src/PMduals.o libs/blossom5-v2.05.src/PMexpand.o  \
		libs/blossom5-v2.05.src/PMinit.o libs/blossom5-v2.05.src/PMinterface.o \
		libs/blossom5-v2.05.src/PMmain.o libs/blossom5-v2.05.src/PMrepair.o \
		libs/blossom5-v2.05.src/PMshrink.o libs/blossom5-v2.05.src/misc.o \
		libs/blossom5-v2.05.src/MinCost/MinCost.o libs/blossom5-v2.05.src/GEOM/GPMinit.o \
		libs/blossom5-v2.05.src/GEOM/GPMinterface.o libs/blossom5-v2.05.src/GEOM/GPMkdtree.o \
		libs/blossom5-v2.05.src/GEOM/GPMmain.o -O3 \


all: simulate

simulate: ${OBJS}
	$(CXX) $(CFLAGS) ${LIBDIR} -o $@ ${OBJS} ${LIBS}

.cpp.o:
	$(CXX) $(CFLAGS) ${INCLUDES} $< -c -o $@ | echo $<

clean:
	rm -f ${OBJS} simulate

		
