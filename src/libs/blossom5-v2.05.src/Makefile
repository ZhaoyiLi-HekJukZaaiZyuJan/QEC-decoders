DIRS := . MinCost GEOM

SOURCES := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.cpp))
OBJS := $(patsubst %.cpp, %.o, $(SOURCES))

CFLAGS := -O3 -D_NDEBUG
CXX ?= c++
LIBS :=
INCLUDES :=
LIBDIR :=

# Add librt if the target platform is not Darwin (OS X)
ifneq ($(shell uname -s),Darwin)
    LIBS += -lrt
endif

all: ${OBJS} #only compile, don't link

# all: blossom5 ##previous version

# blossom5: ${OBJS}
# 	$(CXX) $(CFLAGS) ${LIBDIR} -o $@ ${OBJS} ${LIBS}

.cpp.o:
	$(CXX) $(CFLAGS) ${INCLUDES} $< -c -o $@

clean:
	rm -f ${OBJS} blossom5