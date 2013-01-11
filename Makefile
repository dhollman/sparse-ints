sources := $(wildcard *.cc)
hfiles := $(wildcard *.h)
objfiles := $(sources:.cc=.o)
executable := sparse_ints

# Change this to the path to your installed sc-config script.
SCCONFIG = /Users/dhollman/Projects/mpqc/install/latest_repo/mpi/bin/sc-config

# Get the rest of the details using the sc-config script
CXX := $(shell $(SCCONFIG) --cxx)
CXXFLAGS := $(shell $(SCCONFIG) --cxxflags)
CC := $(shell $(SCCONFIG) --cc)
CCFLAGS := $(shell $(SCCONFIG) --cflags)
F77 := $(shell $(SCCONFIG) --f77)
F90 := $(F77)
FFLAGS := $(shell $(SCCONFIG) --f77flags)
F90FLAGS := $(FFLAGS)
CPPFLAGS := $(shell $(SCCONFIG) --cppflags)
LIBS := $(shell $(SCCONFIG) --libs)
LIBDIR  := $(shell $(SCCONFIG) --libdir)
LTCOMP := $(shell $(SCCONFIG) --ltcomp)
LTLINK := $(shell $(SCCONFIG) --ltlink)
LTLINKBINOPTS := $(shell $(SCCONFIG) --ltlinkbinopts)

# Ignore c++11 warnings
CXXFLAGS += -Wno-c++11-extensions

all: $(executable)

$(executable): $(objfiles)
	$(LTLINK) $(CXX) $(CXXFLAGS) -o $@ $^ -L$(LIBDIR) $(LIBS) $(LTLINKBINOPTS)
	
$(objfiles): %.o: %.cc $(hfiles)
	$(COMPILE.cc) $<
	
clean:
	-rm -f $(executable) $(objfiles) 
