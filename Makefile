INSTALLDIR    = /usr/local/bin

ARCH         := $(shell root-config --arch)

CXX           =
IncSuf        = h
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTLIBDIR   := $(shell root-config --libdir)


ifeq ($(ARCH),linux)
# Linux with egcs, gcc 2.9x, gcc 3.x (>= RedHat 5.2)
CXX           = g++
CXXFLAGS      = -std=c++11 -g -O -Wall -fPIC # -Weffc++
LD            = g++
LDFLAGS       = -g -O
SOFLAGS       = -shared
endif

ifeq (,$(findstring $(ARCH),x8664))
# Linux with egcs, gcc 2.9x, gcc 3.x (>= RedHat 5.2)
CXX           = g++
CXXFLAGS      = -std=c++11 -g -O -Wall -fPIC # -Weffc++
LD            = g++
LDFLAGS       = -g -O
SOFLAGS       = -shared
endif

ifeq ($(CXX),)
$(error $(ARCH) invalid architecture)
endif

ifeq ($(USENCURSES),1)
CXXFLAGS     += -DUSENCURSES
SYSLIBS       = -lmenu -lncurses
endif

#QCINC      = -I$(HOME)/alice/sw/ubuntu1804_x86-64/QualityControl/latest/include
#COMO2INC   = -I$(HOME)/alice/sw/ubuntu1804_x86-64/Common-O2/latest/include
#O2INC      = -I$(HOME)/alice/sw/ubuntu1804_x86-64/O2/latest/include
#BOOSTINC   = -I$(HOME)/alice/sw/ubuntu1804_x86-64/boost/latest/include
#FAIRLOGINC = -I$(HOME)/alice/sw/ubuntu1804_x86-64/FairLogger/latest/include
#FAIRMQINC  = -I$(HOME)/alice/sw/ubuntu1804_x86-64/FairMQ/latest/include/fairmq
#FAIRMQINC2 = -I$(HOME)/alice/sw/ubuntu1804_x86-64/FairMQ/latest/include/
#JALIENINC  = -I$(HOME)/alice/sw/ubuntu1804_x86-64/JAliEn-ROOT/latest/include
#FMTINC     = -I$(HOME)/alice/sw/ubuntu1804_x86-64/fmt/latest/include/

QCINC      = -I/data/dataalice/cdemart/QC/QCUtils/sw/slc7_x86-64/QualityControl/latest/include
COMO2INC   = -I/data/dataalice/cdemart/QC/QCUtils/sw/slc7_x86-64/Common-O2/latest/include
O2INC      = -I/data/dataalice/cdemart/QC/QCUtils/sw/slc7_x86-64/O2/latest/include
BOOSTINC   = -I/data/dataalice/cdemart/QC/QCUtils/sw/slc7_x86-64/boost/latest/include
FAIRLOGINC = -I/data/dataalice/cdemart/QC/QCUtils/sw/slc7_x86-64/FairLogger/latest/include
FAIRMQINC  = -I/data/dataalice/cdemart/QC/QCUtils/sw/slc7_x86-64/FairMQ/latest/include/fairmq
FAIRMQINC2 = -I/data/dataalice/cdemart/QC/QCUtils/sw/slc7_x86-64/FairMQ/latest/include/
JALIENINC  = -I/data/dataalice/cdemart/QC/QCUtils/sw/slc7_x86-64/JAliEn-ROOT/latest/include
FMTINC     = -I/data/dataalice/cdemart/QC/QCUtils/sw/slc7_x86-64/fmt/latest/include/

ROOTINC    = -I$(ROOTSYS)/include
#INCLUDES   = $(QCINC) $(COMO2INC) $(O2INC) $(FAIRMQINC) $(FAIRLOGINC) $(MSGSLINC) $(BOOSTINC) $(ARROWINC) $(INFOLOGINC) $(ROOTINC)
INCLUDES   =  $(BOOSTINC) $(QCINC) $(COMO2INC) $(O2INC) $(ROOTINC) $(FAIRLOGINC) $(FAIRMQINC) $(FAIRMQINC2) $(FMTINC) $(JALIENINC)
CXXFLAGS   += $(ROOTCFLAGS) $(INCLUDES)

#QCLIBS     = -L$(HOME)/alice/sw/ubuntu1804_x86-64/QualityControl/latest/lib -lO2QualityControl
#O2LIBS     = -L$(HOME)/alice/sw/ubuntu1804_x86-64/O2/latest/lib -lO2CCDB
QCLIBS     = -L/data/dataalice/cdemart/QC/QCUtils/sw/slc7_x86-64/QualityControl/latest/lib -lO2QualityControl
O2LIBS     = -L/data/dataalice/cdemart/QC/QCUtils/sw/slc7_x86-64/O2/latest/lib -lO2CCDB

LIBS       = $(QCLIBS) $(O2LIBS) $(ROOTLIBS) $(SYSLIBS)

#------------------------------------------------------------------------------
THISPROG     := getObject
##THISPROGSRC  := $(wildcard *.$(SrcSuf)) DataCompDict.$(SrcSuf)
##THISPROGSRC  := $(wildcard *.$(SrcSuf))
THISPROGSRC  := $(THISPROG).$(SrcSuf)
THISPROGOBJ  := $(THISPROGSRC:.$(SrcSuf)=.$(ObjSuf))
THISPROGEXE  := $(THISPROG)$(ExeSuf)
THISPROGSO   := $(THISPROG).$(DllSuf)

OBJS          = $(THISPROGOBJ)

#PROGRAMS      = $(THISPROGEXE) $(THISPROGSO)
PROGRAMS      = $(THISPROGEXE)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf) .$(IncSuf)

.PHONY: all install clean someclean package

all:            $(PROGRAMS)

$(THISPROGEXE):  $(THISPROGOBJ)
		$(LD) $(LDFLAGS) $^ $(LIBS) -lEG -lHtml -lGeom -lThread -lRMySQL $(OutPutOpt)$@
		@echo "$@ done"

$(THISPROGSO):  $(THISPROGOBJ)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(EXPLLINKLIBS) $(OutPutOpt)$@
		@echo "$@ done"

install: all
		[ -d $(INSTALLDIR) ] || mkdir -p $(INSTALLDIR)
		@cp -p jtag $(INSTALLDIR)


someclean:
		@rm -f $(OBJS) core

clean:      someclean
		@rm -f $(PROGRAMS) *Dict.* *.def *.exp \
		   *.geom *.root *.ps *.so .def so_locations
		@rm -rf cxx_repository

package:
		@rm -f $(THISPROG).tar.gz
		@tar zcpvf $(THISPROG)_`date +%Y%m%d_%H%M%S`.tar.gz *.$(SrcSuf) *.$(IncSuf) Makefile


##DataCompDict.$(SrcSuf):  DigiScanLib.h
##		@echo "****** Generating dictionary $@... ******"
##		@rootcint -f $@ -c $^



###

.$(SrcSuf).$(ObjSuf):
		@echo "****** Compiling " $< " ******"
		$(CXX) $(CXXFLAGS) -c $<
