include $(SKOFL_ROOT)/config.gmk  # pulls in libskroot.so as well
PWD=`pwd`

# C++ compiler flags - XXX config.gmk sets this already, so APPEND ONLY XXX
CXXFLAGS += -fPIC -O3 -g -std=c++11 -fdiagnostics-color=always -Wno-reorder -Wno-sign-compare -Wno-unused-variable -Wno-unused-but-set-variable -Werror=array-bounds -lgfortran # -Wpadded -Wpacked -malign-double -mpreferred-stack-boundary=8  # -Wpedantic << too many pybind warnings?

# debug mode: disable the try{}-catch{} around all Tool methods.
# Combine with -lSegFault to cause exceptions to invoke a segfault, printing a backtrace.
ifeq ($(MAKECMDGOALS),debug)
CXXFLAGS+= -O0 -g -lSegFault -rdynamic -DDEBUG
endif

# flags required for gprof profiling
#CXXFLAGS    += -g -pg -ggdb3

# Fortran compiler flags. APPEND ONLY probably
FCFLAGS += -w -fPIC -lstdc++ -fimplicit-none # -falign-commons

# SK Offline Library (SKOFL) Headers & Libraries
# n.b. no spaces between -I and following path allowed here (rootcint doesn't like it)
SKOFLINCLUDE = -I$(SKOFL_ROOT)/include -I$(SKOFL_ROOT)/inc -I$(SKOFL_ROOT)/include/lowe -I$(SKOFL_ROOT)/inc/lowe

# lowe libraries - some of these may not be required in this list
SKOFLLIB = -L $(SKOFL_ROOT)/lib -lbonsai_3.3 -lsklowe_7.0 -lwtlib_5.1 -lsollib_4.0 -lgeom -lskrd -lastro -lzbs -lgeom -lsklib -llibrary -liolib -lrfa -lskroot -lDataDefinition -ltqrealroot -lloweroot -latmpdroot -lmcinfo -lsofttrgroot -lidod_xtlk_root -lConnectionTableReader

# Atmospheric, Muon and Proton Decay libraries (ATMPD) Headers & Libraries
OLD_NTAG_GD_ROOT = $(ATMPD_ROOT)/src/analysis/neutron/ntag_gd
ATMPDINCLUDE = -I $(ATMPD_ROOT)/include -I $(OLD_NTAG_GD_ROOT) -I $(ATMPD_ROOT)/src/recon/fitqun
ATMPDLIB = -L $(ATMPD_ROOT)/lib -lapdrlib -laplib -lringlib -ltp -ltf -lringlib -laplib -lmsfit -lmslib -lseplib -lmsfit -lprtlib -lmuelib -lffit -lodlib -lstmu -laplowe -laplib -lfiTQun -ltf -lmslib -llelib -lntuple_t2k

# not all fortran routines are built into libraries as part of compiling SKOFL & ATMPD.
# figure out why standalones don't need to specify a full path when listing in dependencies of a target....?
LOCAL_OBJS = $(SKOFL_ROOT)/examples/root2zbs/fort_fopen.o $(ATMPD_ROOT)/src/analysis/neutron/merge/zbsinit.o
#$(ATMPD_ROOT)/src/analysis/neutron/ntag/bonsai.o
#$(RELICWORKDIR)/data_reduc/neutron_tagging/src/bonsai.o
# ATMP version calls lfallfit_sk4_data, which invokes bonsaifit.
# relic version (copied to SK2p2MeV Tool directory for now) invokes bonsaifit, then lfneweff_sk4.
# TODO Make a bonsai tool? Though it seems most scripts invoke bonsai via lfallfit_* or similar...
# Right now it's unclear what actions each of these routines are doing, and their overlaps/differences, though.
LOCAL_LIBS += $(LOCAL_OBJS)

# CERNLib Headers & Libraries
CERNLIB = `cernlib jetset74 graflib grafX11 packlib mathlib kernlib lapack3 blas`

# These MUST be in LDLIBS for things (e.g. cerlib) to work properly!
LD_RUN_PATH=$(SKOFL_ROOT)/lib:$(ATMPD_ROOT)/lib
LDLIBS += $(SKOFLLIB)
LDLIBS += $(ATMPDLIB)
LDLIBS += $(CERNLIB)

# ROOT Headers & Libraries
ROOTINCLUDE= `root-config --cflags`
ROOTLIB = `root-config --libs --evelibs --glibs` -lMinuit -lXMLIO -lMLP
ROOTSTLLIBS = -L${HOME}/stllibs -lRootStl

TMVASYS = $(Dependencies)/TMVA
TMVAINCLUDE = -I $(TMVASYS)/include
TMVALIB = -L $(TMVASYS)/lib -lTMVA.1

# Python Headers & Libraries
HASPYTHON3 := $(shell python3-config --cflags >/dev/null 2>&1; echo $$?)
ifeq ($(HASPYTHON3),0)
    CXXFLAGS+= -DPYTHON=1 
    PythonInclude = `python3-config --cflags`
    PythonLib = `python3-config --ldflags --libs`
    # note that pybind11 documentation says the following
    # | it's better to (intentionally) *not* link against libpython.
    # | The symbols will be resolved when the extension library is loaded into a Python binary.
    # | This is preferable because you might have several different installations of a given Python version
    # | (e.g. the system-provided Python, and one that ships with a piece of commercial software).
    # | In this way, the plugin will work with both versions, instead of possibly importing a second 
    # | Python library into a process that already contains one (which will lead to a segfault).
    # ...?? both `python3-config --ldflags` and `python3-config --libs` both pull in `libpython3.6m`,
    # so maybe it means "don't manually add `-lpython`"?
    
    # Pybind11
    PybindInclude=`python3 -m pybind11 --includes`
    # we also need this flag
#    CXXFLAGS+= -fvisibility=hidden
    # see the following for details
    # pybind11.readthedocs.io/en/stable/faq.html#someclass-declared-with-greater-visibility-than-the-type-of-its-field-someclass-member-wattributes
    # but it can only be set for building DataModel.lib or we get all sorts of undefined references
    # *** IMPORTANT *** #
    # The pybind header includes...
    #| #include <pybind11/pybind11.h>
    #| #include <pybind11/embed.h>
    #| #include <pybind11/numpy.h>
    # ... MUST be placed at the top of the DataModel.h file.
    # Dunno why, but otherwise it complains "invalid conversion from 'int' to 'const char*'"
else
    $(warn "Did not find python3-config...")
endif

# support compression in BStores
# must be either defined or not for both this and ToolFrameworkCore in Dependencies!
#ZBLIB= -lz
#ZLIBFLAG= -DZLIB
#CXXFLAGS += $(ZLIBFLAG) $(ZLIB)

# Third Reduction Library (part of SRN analysis)
THIRDREDLIB = -L${HOME}/relic_sk4_ana/relic_work_dir/data_reduc/third/lib -lthirdredvars

# all user classes that the user may wish to write to ROOT files require a dictionary.
# TODO maybe we should put these in a separate directory or something so they don't need to be listed explicitly
ROOTCLASSES = DataModel/Candidate.h DataModel/Cluster.h DataModel/EventCandidates.h DataModel/EventParticles.h DataModel/EventTrueCaptures.h DataModel/Particle.h DataModel/PMTHitCluster.h DataModel/PMTHit.h DataModel/TrueCapture.h
#DataModel/Calculator.h

# Combine all external libraries and headers needed by the DataModel
DataModelInclude = $(ROOTINCLUDE) $(SKOFLINCLUDE) $(ATMPDINCLUDE) $(PybindInclude)
DataModelLib = $(ROOTLIB)

# Combine all external libraries and headers needed by user Tools
MyToolsInclude = $(SKOFLINCLUDE) $(ATMPDINCLUDE) $(PythonInclude) $(TMVAINCLUDE)
MyToolsLib = $(LDFLAGS) $(LDLIBS) $(PythonLib) $(THIRDREDLIB) $(TMVALIB) $(ROOTSTLLIBS)

# To add user classes:
# 1. Add a rule to build them into a shared library (see libMyClass.so example)
#    (the ouput library must be built into libs/libXXX.so)
# 2. Add the library here
UserLibs = lib/libMyClass.so

# these can't be put in rules, they won't be executed there.
USERLIBS1=$(patsubst %.so,%,$(UserLibs))
USERLIBS2=$(patsubst lib/lib%,-l%,$(USERLIBS1))

Dependencies=Dependencies

debug: all

all: lib/libMyTools.so lib/libToolChain.so lib/libStore.so include/Tool.h lib/libDataModel.so lib/libLogging.so main

main: src/main.cpp lib/libStore.so lib/libLogging.so lib/libToolChain.so | lib/libMyTools.so lib/libDataModel.so  lib/liblowfit_sk4_stripped.so lib/libRootDict.so lib/libBStore_RootDict.so $(UserLibs)
	@echo -e "\e[38;5;214m\n*************** Making " $@ "****************\e[0m"
	g++ $(CXXFLAGS) -L lib -llowfit_sk4_stripped -I include $(DataModelInclude) $(MyToolsInclude) src/main.cpp -o $@ $(DataModelLib) $(MyToolsLib) -L lib -lStore -lMyTools -lToolChain -lDataModel -lLogging -lpthread $(ROOTLIB) $(ATMPDLIB) $(SKOFLLIB) $(CERNLIB) -lRootDict $(USERLIBS2)

lib/libStore.so: $(Dependencies)/ToolFrameworkCore/src/Store/*
	cd $(Dependencies)/ToolFrameworkCore && $(MAKE) lib/libStore.so
	@echo -e "\e[38;5;118m\n*************** Copying " $@ "****************\e[0m"
	cp $(Dependencies)/ToolFrameworkCore/src/Store/*.h include/
	cp $(Dependencies)/ToolFrameworkCore/lib/libStore.so lib/

# note: these MUST use the same `-DZLIB` and `-lz` options as libStore and main!
BStore_RootDict.cxx: include/BStore.h BStore_Linkdef.hh | UserTools/Factory/Factory.o
	@echo "making $@"
	rootcint -f $@ -c -p $(ZLIBFLAG) -fPIC `root-config --cflags` -I./include $^

lib/libBStore_RootDict.so: BStore_RootDict.cxx lib/libStore.so lib/libRootStl.so $(UserLibs)
	@echo "making $@"
	# pull in libraries defined in BStore_RootDict_libs.txt (drop first comment line)
	g++ $(CXXFLAGS) -shared `root-config --cflags` -I./include $< -L lib -lStore $(USERLIBS2) `root-config --libs` -lRootStl -o $@

lib/libRootStl.so: $(Dependencies)/RootStl/libRootStl.so
	cp $^ $@

include/Tool.h:  $(Dependencies)/ToolFrameworkCore/src/Tool/Tool.h
	@echo -e "\e[38;5;118m\n*************** Copying " $@ "****************\e[0m"
	cp $(Dependencies)/ToolFrameworkCore/src/Tool/Tool.h include/

include/BStore.h:  $(Dependencies)/ToolFrameworkCore/src/Store/Store.h
	@echo -e "\e[38;5;118m\n*************** Copying " $@ "****************\e[0m"
	cp $(Dependencies)/ToolFrameworkCore/src/Store/Store.h include/

lib/libToolChain.so: $(Dependencies)/ToolFrameworkCore/src/ToolChain/*  lib/libLogging.so lib/libStore.so | lib/libMyTools.so lib/libDataModel.so
	@echo -e "\e[38;5;214m\n*************** Making " $@ "****************\e[0m"
	cp $(Dependencies)/ToolFrameworkCore/src/ToolChain/*.h include/
	g++ $(CXXFLAGS) -shared $(Dependencies)/ToolFrameworkCore/src/ToolChain/ToolChain.cpp -I include -lpthread -L lib -lStore -lDataModel -lLogging -lMyTools -o lib/libToolChain.so $(DataModelInclude) $(DataModelLib) $(MyToolsInclude) $(MyToolsLib)

clean: 
	@echo -e "\e[38;5;201m\n*************** Cleaning up ****************\e[0m"
	rm -f include/*.h
	rm -f lib/*.so
	rm -f main
	rm -f UserTools/*/*.o
	rm -f DataModel/*.o

lib/libDataModel.so: DataModel/* lib/libLogging.so lib/libStore.so  $(patsubst DataModel/%.cpp, DataModel/%.o, $(wildcard DataModel/*.cpp)) lib/libRootDict.so
	@echo -e "\e[38;5;214m\n*************** Making " $@ "****************\e[0m"
	g++ $(CXXFLAGS) -fvisibility=hidden -shared DataModel/*.o -I include -L lib -lStore -lLogging -o lib/libDataModel.so $(DataModelInclude) $(DataModelLib) -lRootDict


lib/libMyTools.so: UserTools/*/* UserTools/* lib/libStore.so include/Tool.h lib/libLogging.so UserTools/Factory/Factory.o | lib/libDataModel.so 
	@echo -e "\e[38;5;214m\n*************** Making " $@ "****************\e[0m"
	g++ $(CXXFLAGS) -shared UserTools/*/*.o -I include -L lib -lStore -lDataModel -lLogging -o lib/libMyTools.so $(MyToolsInclude) $(DataModelInclude) $(MyToolsLib) $(DataModelLib)

lib/libLogging.so:  $(Dependencies)/ToolFrameworkCore/src/Logging/* | lib/libStore.so
	cd $(Dependencies)/ToolFrameworkCore && $(MAKE) lib/libLogging.so
	@echo -e "\e[38;5;118m\n*************** Copying " $@ "****************\e[0m"
	cp $(Dependencies)/ToolFrameworkCore/src/Logging/Logging.h include/
	cp $(Dependencies)/ToolFrameworkCore/lib/libLogging.so lib/

UserTools/Factory/Factory.o: UserTools/Factory/Factory.cpp lib/libStore.so include/Tool.h lib/libLogging.so lib/libDataModel.so  $(filter-out UserTools/Factory/Factory.o, $(patsubst UserTools/%.cpp, UserTools/%.o, $(wildcard UserTools/*/*.cpp)) $(patsubst UserTools/%.cc, UserTools/%.o, $(wildcard UserTools/*/*.cc)) $(patsubst UserTools/%.F, UserTools/%.o, $(wildcard UserTools/*/*.F))) | include/Tool.h
	@echo -e "\e[38;5;214m\n*************** Making " $@ "****************\e[0m"
	cp UserTools/Factory/Factory.h include
	cp UserTools/Unity.h include
	-g++ $(CXXFLAGS) -c -o $@ $< -I include -L lib -lStore -lDataModel -lLogging $(MyToolsInclude) $(MyToolsLib) $(DataModelInclude) $(DataModelib)

update:
	@echo -e "\e[38;5;51m\n*************** Updating ****************\e[0m"
	cd $(Dependencies)/ToolFrameworkCore; git pull
	git pull

UserTools/%.o: UserTools/%.cpp lib/libStore.so include/Tool.h lib/libLogging.so lib/libDataModel.so | include/Tool.h
	@echo -e "\e[38;5;214m\n*************** Making " $@ "****************\e[0m"
	-cp $(shell dirname $<)/*.h include
	-g++ $(CXXFLAGS) -c -o $@ $< -I include -L lib -lStore -lDataModel -lLogging $(MyToolsInclude) $(MyToolsLib) $(DataModelInclude) $(DataModelLib)

UserTools/%.o: UserTools/%.cc lib/libStore.so include/Tool.h lib/libLogging.so lib/libDataModel.so
	@echo -e "\e[38;5;214m\n*************** Making c++ Tool " $@ "****************\e[0m"
	-cp $(shell dirname $<)/*.h include 2>/dev/null || :
	-g++ $(CXXFLAGS) -c -o $@ $< -I include -L lib -lStore -lDataModel -lLogging $(MyToolsInclude) $(MyToolsLib) $(DataModelInclude) $(DataModelLib)

UserTools/%.o: UserTools/%.F lib/libStore.so include/Tool.h lib/libLogging.so lib/libDataModel.so
	@echo -e "\e[38;5;214m\n*************** Making fortran Tool " $@ "****************\e[0m"
	-cp $(shell dirname $<)/*.h include 2>/dev/null || :
	-gfortran $(FCFLAGS) -c -o $@ $< $(SKOFLINCLUDE) -I $(ATMPD_ROOT)/inc

target: remove $(patsubst %.cpp, %.o, $(wildcard UserTools/$(TOOL)/*.cpp)) $(patsubst %.cc, %.o, $(wildcard UserTools/$(TOOL)/*.cc)) $(patsubst UserTools/%.F, UserTools/%.o, $(wildcard UserTools/*/*.F))

remove:
	echo "removing"
	-rm UserTools/$(TOOL)/*.o

DataModel/%.o: DataModel/%.cpp lib/libLogging.so lib/libStore.so  
	@echo -e "\e[38;5;214m\n*************** Making " $@ "****************\e[0m"
	cp $(shell dirname $<)/*.h include
	-g++ $(CXXFLAGS) -c -o $@ $< -I include -L lib -lStore -lLogging  $(DataModelInclude) $(DataModelLib)

Docs:
	doxygen Doxyfile

# this is the 'hack' library that must be linked into the final executable,
# which together with certain orderings of arguments (place the library EARLY)
# somehow avoids cernlib 64-bit address issues. Do not edit this or the `main` target method!
lib/liblowfit_sk4_stripped.so: DataModel/lowfit_sk4_stripped.cc DataModel/lowfit_sk4_stripped.h
	@echo -e "\n*************** Making " $@ "****************"
	-g++ $(CXXFLAGS) -shared $(SKOFLINCLUDE) -L${CERN_ROOT}/lib $< -o $@

# TODO maybe we should have some way to distinguish DataModel classes that we wish to have a root dictionary for
DataModel/NTagDataModelDict.cxx: $(ROOTCLASSES) DataModel/NTagLinkDef.hh
	@echo "making $@"
	# the dictionary sourcefiles need to be built within the directory they're to reside in
	# otherwise the paths encoded into the dictionary don't match those when the object is built.
	# why doesn't this complain about c++11 stuff?
	# maybe because he uses rootcint to build the dictionary as an object
	# then builds a lib from all the objects using g++ with -std=c++11 (in CXXFLAGS) ??
	# but then the libRootDict.so doesn't use the object, just uses the dictionary cxx.
	# i think the latter `root-config --cxx ... `command here isn't required... 
	cd DataModel && \
	rootcint -f NTagDataModelDict.cxx -c -I./include -I./ -I$(shell root-config --incdir) $(SKOFLINCLUDE) $(notdir $^)

lib/libRootDict.rootmap:
	@echo "Library.EventParticles: $(PWD)/lib/libRootDict.so" >> $@
	@echo "Library.Particle: $(PWD)/lib/libRootDict.so" >> $@
	@echo "Library.EventTrueCaptures: $(PWD)/lib/libRootDict.so" >> $@
	@echo "Library.TrueCapture: $(PWD)/lib/libRootDict.so" >> $@
	@echo "Library.Cluster: $(PWD)/lib/libRootDict.so" >> $@
	@#echo "Library.Cluster<Particle>: $(PWD)/lib/libRootDict.so" >> $@
	@#echo "Library.Cluster<TrueCapture>: $(PWD)/lib/libRootDict.so" >> $@

lib/libCalculator.so: DataModel/Calculator.cpp DataModel/Calculator.h
	-g++ $(CXXFLAGS) -shared -fPIC `root-config --cflags` $< -o $@ `root-config --libs`

lib/libRootDict.so: DataModel/NTagDataModelDict.cxx lib/libRootDict.rootmap lib/libCalculator.so $(ROOTCLASSES:%.h=%.o)
	@echo "making lib/libRootDict.so"
	`root-config --cxx --cflags` $(CXXFLAGS) -W -Wall -fPIC -shared -o $@ $(ROOTCLASSES:%.h=%.o) $(SKOFLINCLUDE) $< `root-config --glibs` $(SKOFLLIB) -L lib/ -lCalculator

#################################################################################
## Adding user classes:                                                        ##
## Follow the below template to build a ROOT dictionary for your class         ##
## Then add the library and its directory to BStore_RootDict_libs.txt          ##
#################################################################################
lib/libMyClass.so: UserTools/TestTool/MyClass.cpp UserTools/TestTool/MyClass.h UserTools/TestTool/MyClass_Linkdef.hh
	@echo "making $@"
	cd UserTools/TestTool && \
	rootcint -f MyClass_RootDict.cxx -c -p -I./ -fPIC `root-config --cflags` MyClass.h MyClass_Linkdef.hh && \
	cd ../../ && \
	g++ -shared -fPIC `root-config --cflags` UserTools/TestTool/MyClass_RootDict.cxx $^ -L./ `root-config --libs` -o $@
