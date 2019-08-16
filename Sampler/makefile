#******************************************************************************#
# Copyright (c) 2018
#
#******************************************************************************#

OBJDIR = ../.o
DEPDIR = ../.d
$(shell mkdir -p $(OBJDIR) >/dev/null)
$(shell mkdir -p $(DEPDIR) >/dev/null)

#CXX = g++-7 # XPS
CXX = mpicxx # XPS - mpi wrapper
#CXX = g++ # Compute servers
#CXX = cc # Archer

CXXFLAGS  = -g -std=c++11 -O3 -fopenmp -fopenmp-simd $(INCLUDE) -lstdc++ -ftree-vectorize -funroll-loops -ffast-math
CPP_FILES:= $(shell find ../src -type f -name '*.cpp')
# GSLFLAGS  = -lgsl -lgslcblas -lm
HH_FILES := $(shell find ../src -type f -name '*.hpp')
INCLUDE  := $(sort $(dir $(CPP_FILES)+$(HH_FILES)))
INCLUDE  := $(patsubst %, -I%, $(INCLUDE))
OBJS     := $(patsubst ../%.cpp, $(OBJDIR)/%.o, $(CPP_FILES))
BIN       = ../bin

DEPS:= $(patsubst $(OBJDIR)/%.o,$(DEPDIR)/%.d,$(OBJS))

CUROBJDIR = .o
CURDEPDIR = .d
$(shell mkdir -p $(CUROBJDIR) >/dev/null)
$(shell mkdir -p $(CURDEPDIR) >/dev/null)
CUR_CPP_FILES :=$(shell find . -type f -name '*.cpp')
CUR_HH_FILES  :=$(shell find . -type f -name '*.hh')
CUR_INCLUDE   := $(sort $(dir $(CUR_CPP_FILES)+$(CUR_HH_FILES)))
CUR_INCLUDE   := $(patsubst %, -I%, $(INCLUDE))
CUR_OBJS      := $(patsubst %.cpp, $(CUROBJDIR)/%.o, $(CUR_CPP_FILES))

CURDEPS:= $(patsubst $(CUROBJDIR)/%.o,$(CURDEPDIR)/%.d,$(CUR_OBJS))

#MAIN=LJcluster
MAIN=DoubleWell
#MAIN=MercedesBenz

#MercBenz: MercedesBenz.cpp ${OBJS} ${CUR_OBJS}
#LJcluster: LJcluster.cpp ${OBJS} ${CUR_OBJS}
$(MAIN): $(MAIN).cpp ${OBJS} ${CUR_OBJS}
	@mkdir -p $(BIN)
	@echo ${OBJS}
	# $(CXX) $(CXXFLAGS) -o $(BIN)/MercBenz ${OBJS} ${CUR_OBJS} -lstdc++ -lgsl -lgslcblas
	# $(CXX) $(CXXFLAGS) -o $(BIN)/LJcluster ${OBJS} ${CUR_OBJS} -lstdc++ -lgsl -lgslcblas
	$(CXX) $(CXXFLAGS) -o $(BIN)/$(MAIN) ${OBJS} ${CUR_OBJS} -lstdc++ -lgsl -lgslcblas

-include $(DEPS)
-include $(CURDEPS)

${CUROBJDIR}/%.o : %.cpp $(wildcard %.hpp)
	@mkdir -p $(dir $@)
	@mkdir -p $(CURDEPDIR)/$(dir $*.cpp)
	$(CXX) $(CXXFLAGS) -c -o ${CUROBJDIR}/$*.o $*.cpp
	$(CXX) $(CXXFLAGS) -MM $*.cpp > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|${CUROBJDIR}/$*.o:|' < $*.d.tmp > $*.d
	@mv $*.d ${CURDEPDIR}/$*.d
	@rm -f $*.d.tmp

${OBJDIR}/%.o : ../%.cpp $(wildcard ../%.hpp)
	@mkdir -p $(dir $@)
	@mkdir -p $(DEPDIR)/$(dir $*.cpp)
	$(CXX) $(CXXFLAGS) -c -o ${OBJDIR}/$*.o ../$*.cpp
	$(CXX) $(CXXFLAGS) -MM ../$*.cpp > ../$*.d
	@mv -f ../$*.d ../$*.d.tmp
	@sed -e 's|.*:|${OBJDIR}/$*.o:|' < ../$*.d.tmp > ../$*.d
	@mv ../$*.d ${DEPDIR}/$*.d
	@rm -f ../$*.d.tmp

debug: CXXFLAGS += -fsanitize=address -fno-omit-frame-pointer \
                   -fsanitize-recover=address
debug: Sampler

clean:
	rm -rf $(OBJDIR) $(DEPDIR) $(CUROBJDIR) $(CURDEPDIR) $(BIN)/Sampler
