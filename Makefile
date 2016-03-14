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

OBJS := $(addprefix $(BUILDDIR)/, FilteredAnitaEvent.o FilterOperation.o FilterStrategy.o AnalysisWaveform.o AnitaEventSummary.o analysisDict.o)

#BINARIES := $(addprefix $(BINDIR)/, binary);

# INCLUDES := $(addprefix $(INCLUDEDIR)/, $(shell ls includes/*.h ))
INCLUDES := $(shell ls $(INCLUDEDIR)/*.h )

ROOT_LIBRARY = $(LIBDIR)/libAnitaAnalysis.${DllSuf}

all: $(ROOT_LIBRARY) $(BINARIES) 

## probably need some magic for Mac OS X here? Yes you do, and it is a truly mysterious dark art
$(ROOT_LIBRARY) : $(OBJS) | $(LIBDIR) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(LIBS) $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $(LIBS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
#	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $(OBJS) -o $@
# $(LIBNAME): $(OBJS) | $(LIBDIR)
# 	@echo Building shared library $@
	$(LD) $(SOFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) $(GLIBS) -shared -o $@
endif




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


$(BINDIR)/%: %.cc $(INCLUDES) Makefile $(ROOT_LIBRARY) | $(BINDIR)
	@echo Compiling $<
	$(LD)  -I./include -I./ $(CXXFLAGS) -o $@ $(LDFLAGS) -L./$(LIBDIR) -lAnitaAnalysis  $< 

$(BUILDDIR)/analysisDict.cc: $(INCLUDES) LinkDef.h | $(BUILDDIR)
	@echo Running rootcint
	rootcint  -f $@ -c -p -I$(ANITA_UTIL_INSTALL_DIR)/include $(INCLUDES) LinkDef.h

install: $(ROOT_LIBRARY)
ifndef ANITA_UTIL_INSTALL_DIR 
	$(error Please define ANITA_UTIL_INSTALL_DIR)
endif 
	install -d $(ANITA_UTIL_INSTALL_DIR)/lib 
	install -d $(ANITA_UTIL_INSTALL_DIR)/include 
	install -c -m 755 $(ROOT_LIBRARY) $(ANITA_UTIL_INSTALL_DIR)/lib
	install -c -m 644 $(INCLUDES) $(ANITA_UTIL_INSTALL_DIR)/include 
	if [ -e $(BUILDDIR)/analysisDict_rdict.pcm ];  then install -c -m 755 $(BUILDDIR)/analysisDict_rdict.pcm $(ANITA_UTIL_INSTALL_DIR)/lib; fi; 
clean: 
	rm -rf build
	rm -rf bin
	rm -rf lib
