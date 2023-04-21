#definition of a recursive wildcard "ref: https://stackoverflow.com/questions/2483182/recursive-wildcards-in-gnu-make"
rwildcard=$(foreach d,$(wildcard $(1:=/*)),$(call rwildcard,$d,$2) $(filter $(subst *,%,$2),$d))

DIRS := . ../src
SOURCES := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.cpp))
#search for all object files in blossom5 folder
BLOSSOM_OBJECTS:= $(call rwildcard, ../src/libs/blossom5-v2.05.src,*.o)
OBJS := $(patsubst %.cpp, %.o, $(SOURCES))

CFLAGS := -std=c++17 -ltensorflow -lpthread -no-pie -fmax-errors=5 -O3 \
		-I ../src \
		-I ../src/libs/ \
		-I ../src/libs/blossom5-v2.05.src \
		-I ~/libtensorflow2/include/ \
		-I ../src/libs/cppflow-master/include \
		-I ../src/libs/cppflow-master/include/cppflow
CXX ?= g++
INCLUDES :=
LIBDIR :=

all: simulate

simulate: ${OBJS}
	export LIBRARY_PATH=$(LIBRARY_PATH):~/libtensorflow2/lib:/usr/local/Cellar/libtensorflow/lib; \
	export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH):~/libtensorflow2/lib:/usr/local/Cellar/libtensorflow/lib; \
	$(CXX) $(CFLAGS) ${BLOSSOM_OBJECTS} -o $@ ${OBJS} ${LIBS} \

.cpp.o:
	$(CXX) $(CFLAGS) ${INCLUDES} $< -c -o $@ | echo $<

clean:
	rm -f ${OBJS} simulate

		
