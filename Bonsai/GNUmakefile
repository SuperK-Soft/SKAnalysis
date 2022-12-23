#
#  In case of building locally, set SKOFL_ROOT variable 
#      setenv SKOFL_ROOT   ... 
#  or directly set in this makefile 
# SKOFL_ROOT = /skofl
#
ifndef SKOFL_ROOT
  SKOFL_ROOT = ../..
endif

include $(SKOFL_ROOT)/config.gmk

CXXOPTIONS += -D ConstDirC="\"$(CONSTDIR):$(CONSTDIR)/lowe:$(SKAMCONSTDIR):\""

#
#  Objects
#

LIBNAME = bonsai_3.3

OBJS=	fit_param.o pmt_geometry.o binfile.o bonsai.o bonsaifit.o plato.o \
	searchgrid.o fourhitgrid.o centroid.o hits.o hitsel.o goodness.o \
	timefit.o likelihood.o vertex.o pot.o tree.o vertexfit.o \
	bscalls.o ariadne.o

#
#  Rules for building library 
#

.PHONY:  lib$(LIBNAME).a $(LIBDIR)lib$(LIBNAME).a

lib$(LIBNAME).a : $(OBJS)
	$(RM) $@
	$(AR) $@ $(OBJS) 
	$(RANLIB) $@

$(LIBDIR)lib$(LIBNAME).a : lib$(LIBNAME).a
	$(RM) $@
	$(INSTALL_LIB) $< $@


#
#  Targets 
#

.PHONY:  clean setup includes install.includes depend lib install.lib exec install.exec


emptyrule:: lib


clean::
	$(RM) *.o *~ *.a *.so core* 

setup::

includes::

install.includes::

depend::

lib:: lib$(LIBNAME).a

install.lib:: $(LIBDIR)lib$(LIBNAME).a

exec::

install.exec:: 
 
