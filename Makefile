ccfiles := $(wildcard *.cc)
hfiles := $(wildcard *.h)
objfiles := $(ccfiles:.cc=.o)
executables := $(ccfiles:%.cc=%)

# Change this to the path to your installed sc-config script.
SCCONFIG = /Users/dhollman/Projects/mpqc/install/latest_repo/mpi/bin/sc-config

# This doesn't work
# $(ccfiles:%.cc=%_nompi) : SCCONFIG = /Users/dhollman/Projects/mpqc/install/latest_repo/debug/bin/sc-config

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


CXXFLAGS += -Wno-c++11-extensions

$(objfiles) : $(hfiles)

$(executables): %: %.o
	$(LTLINK) $(CXX) $(CXXFLAGS) -o $@ $^ -L$(LIBDIR) $(LIBS) $(LTLINKBINOPTS)
	
#%: %.o 
#$(LTLINK) $(CXX) $(CXXFLAGS) -o $@ $^ -L$(LIBDIR) $(LIBS) $(LTLINKBINOPTS)
	
#%.o: %.cc $(hfiles)
	
	
clean:
	-rm -f $(executables) $(objfiles) 
