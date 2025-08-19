# ---- sources/objects -------------------------------------------------
# recursive wildcard (ref: https://stackoverflow.com/q/2483182)
rwildcard = $(foreach d,$(wildcard $(1)/*),$(call rwildcard,$d,$2)) \
            $(filter $(subst *,%,$2),$(wildcard $(1)/$2))

DIRS     := . ../src
# (non-recursive sources, same as your original)
SOURCES  := $(foreach dir,$(DIRS),$(wildcard $(dir)/*.cpp))
OBJS     := $(patsubst %.cpp,%.o,$(SOURCES))

# search for all object files in blossom5 folder (recursive)
BLOSSOM_OBJECTS := $(call rwildcard,../src/libs/blossom5-v2.05.src,*.o)

# ---- toolchain -------------------------------------------------------
CXX ?= g++          # on macOS you can also use: clang++

# ---- include & compile flags ----------------------------------------
CXXFLAGS := -std=c++17 -O3 \
            -I/usr/local/include \
            -I../src \
            -I../src/libs \
            -I../src/libs/blossom5-v2.05.src \
            -I../src/libs/cppflow-master/include \
            -I../src/libs/cppflow-master/include/cppflow
# clang warns -fmax-errors is unused; keep or drop as you like:
# CXXFLAGS += -fmax-errors=5

# ---- link flags & libs ----------------------------------------------
# Tell the linker where to find libtensorflow, and embed an rpath so dyld
# finds it without DYLD_LIBRARY_PATH (macOS).
LDFLAGS  := -L/usr/local/lib -Wl,-rpath,/usr/local/lib
LDLIBS   := -ltensorflow -lpthread

# ---- targets ---------------------------------------------------------
.PHONY: all clean
all: simulate

simulate: $(OBJS) $(BLOSSOM_OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(BLOSSOM_OBJECTS) $(LDLIBS)

# compile rule (no link flags here)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@echo "compiled $<"

clean:
	rm -f $(OBJS) simulate