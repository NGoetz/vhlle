
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)

CXX           = g++
CXXFLAGS      = -Wall -fPIC -O1 -g
LD            = g++
LDFLAGS       = -O1 -g

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)

vpath %.cpp src
objdir     = obj

SRC        = cll.cpp eos.cpp eo3.cpp eo1.cpp eoChiral.cpp eoCMF.cpp eoHadron.cpp eoAZH.cpp eoSmash.cpp trancoeff.cpp fld.cpp hdo.cpp s95p.cpp icurqmd.cpp ic.cpp ickw.cpp icPartUrqmd.cpp icPartSMASH.cpp main.cpp rmn.cpp cornelius.cpp \
             icGlauber.cpp icGubser.cpp icGlissando.cpp icTrento.cpp icIPG.cpp
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC))

TARGET	   = hlle_visc
#------------------------------------------------------------------------------
$(TARGET):       $(OBJS)
		$(LD)  $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(TARGET)

$(OBJS): | $(objdir)

$(objdir):
	@mkdir -p $(objdir)

obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
