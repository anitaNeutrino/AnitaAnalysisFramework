include Makefile.arch
include Makefile.config

CXXFLAGS+= -g

CXXFLAGS     += $(ROOTCFLAGS) $(SYSINCLUDES) -I$(ANITA_UTIL_INSTALL_DIR)/include
LDFLAGS      += $(ROOTLDFLAGS)  -L$(ANITA_UTIL_INSTALL_DIR)/lib
#LIBS          = $(ROOTLIBS) -g -Wl,-z,defs -lMathMore -lRootFftwWrapper -lAnitaEvent
LIBS          = $(ROOTLIBS) -g -lMathMore -lRootFftwWrapper -lAnitaEvent
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)
LIBDIR=lib
BUILDDIR=build
INCLUDEDIR=include
BINDIR=bin

.PHONY: clean install all

OBJS := $(addprefix $(BUILDDIR)/, FilteredAnitaEvent.o FilterOperation.o FilterStrategy.o AnalysisWaveform.o AnitaEventSummary.o dict.o)

#BINARIES := $(addprefix $(BINDIR)/, binary);

# INCLUDES := $(addprefix $(INCLUDEDIR)/, $(shell ls includes/*.h ))
INCLUDES := $(shell ls $(INCLUDEDIR)/*.h )

LIBNAME = $(LIBDIR)/libAnitaAnalysis.${DllSuf}
LINKLIBNAME=AnitaAnalysis

all: $(LIBNAME) $(BINARIES) 

### probably need some magic for Mac OS X here? 
$(LIBNAME): $(OBJS) | $(LIBDIR)
	@echo Building shared library $@
	@$(LD) $(SOFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) $(GLIBS) -shared -o $@


$(OBJS): | $(BUILDDIR)

$(BUILDDIR): 
	mkdir -p $(BUILDDIR)

$(BINDIR): 
	mkdir -p $(BINDIR)

$(LIBDIR): 
	mkdir -p $(LIBDIR)


$(BUILDDIR)/%.o: src/%.cc $(INCLUDES) Makefile | $(BUILDDIR) $(VECTORIZE)
	@echo Compiling  $< 
	$(LD)  -I./include $(CXXFLAGS) -o $@ -c $< 

$(BUILDDIR)/%.o: build/%.cc $(INCLUDES) Makefile | $(BUILDDIR) 
	@echo Compiling  $< 
	$(LD)  -I../include -I./ $(CXXFLAGS) -o $@ -c $< 


$(BINDIR)/%: %.cc $(INCLUDES) Makefile $(LIBNAME) | $(BINDIR)
	@echo Compiling $<
	$(LD)  -I./include -I./ $(CXXFLAGS) -o $@ $(LDFLAGS) -L./$(LIBDIR) -lAnitaAnalysis  $< 

$(BUILDDIR)/dict.cc: $(INCLUDES) LinkDef.h | $(BUILDDIR)
	@echo Running rootcint
	rootcint  -f $@ -c -p -I$(ANITA_UTIL_INSTALL_DIR)/include $(INCLUDES) LinkDef.h

install: 
ifndef ANITA_UTIL_INSTALL_DIR 
	$(error Please define ANITA_UTIL_INSTALL_DIR)
endif 
	install -d $(ANITA_UTIL_INSTALL_DIR)/lib 
	install -d $(ANITA_UTIL_INSTALL_DIR)/include 
	install -c -m 755 $(LIBNAME) $(ANITA_UTIL_INSTALL_DIR)/lib
	install -c -m 644 $(INCLUDES) $(ANITA_UTIL_INSTALL_DIR)/include 

clean: 
	rm -rf build
	rm -rf bin
	rm -rf lib
