# ---- sources/objects -------------------------------------------------
rwildcard = $(foreach d,$(wildcard $(1:=/*)),$(call rwildcard,$d,$2) $(filter $(subst *,%,$2),$d))

DIRS     := . ../src
SOURCES  := $(foreach dir,$(DIRS),$(wildcard $(dir)/*.cpp))
OBJS     := $(patsubst %.cpp,%.o,$(SOURCES))

# prebuilt blossom .o files (if that’s really what you have)
BLOSSOM_OBJECTS := $(call rwildcard,../src/libs/blossom5-v2.05.src,*.o)

# ---- toolchain -------------------------------------------------------
CXX      ?= g++
# or on macOS you might prefer:
# CXX    ?= clang++

# ---- include & compile flags ----------------------------------------
CXXFLAGS := -std=c++17 -fmax-errors=5 -O3 \
            -I/usr/local/include \
            -I../src \
            -I../src/libs \
            -I../src/libs/blossom5-v2.05.src \
            -I../src/libs/cppflow-master/include \
            -I../src/libs/cppflow-master/include/cppflow

# ---- link flags & libs ----------------------------------------------
# Tell the linker where to find libtensorflow, and embed an rpath so the
# runtime loader (dyld) finds it without DYLD_LIBRARY_PATH.
LDFLAGS  := -L/usr/local/lib -Wl,-rpath,/usr/local/lib
LDLIBS   := -ltensorflow -lpthread

# (macOS doesn’t need -no-pie; remove it)

# ---- targets ---------------------------------------------------------
all: simulate

simulate: $(OBJS) $(BLOSSOM_OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(BLOSSOM_OBJECTS) $(LDLIBS)

# compile rule (no link flags here)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@echo "compiled $<"

clean:
	rm -f $(OBJS) simulate