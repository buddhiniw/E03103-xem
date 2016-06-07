#
# NOTE: in order to use this for the ntuples created on general, recon.f MUST
# be compiled on general!!!!!!!
#


#buddhini 02-03-2015
my_objs = dummy_ana.o hbook_init.o makedathistos.o makesimhistos.o get_physics.o getcharge.o getkine.o getgenvol.o sig_xiscale.o calc_xsec_from_delta.o calc_xsec_from_x.o dis_smear.o xsechd.o nform.o y_calc.o fy.o target_info.o r1990.o calc_eloss.o enerloss_new.o gauss1.o mt19937.o total_eloss.o h_targ_trans.o h_targ_trans_init.o sig_bar_df.o 

#INC = c_ntuple.cmn hms.cmn sos.cmn physics.cmn hms_recon_elements.cmn \
#      sos_recon_elemnts.cmn
RM        = rm -f 
SHELL     = /bin/sh

#ifndef MYOS  # New tcsh defines this
  MYOS := $(subst -,,$(shell uname))
#endif

#CERNLIBS = -lgeant$(GEANTVER) -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib
#OURLIBS = $(UTILLIB) $(GMCLIB)

CERN_ROOT = /apps/cernlib/i386_rhel3/2003
CERNLIBS = -L$(CERN_ROOT)/lib -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib

ifeq ($(MYOS),HPUX)
  FFLAGS=+U77 +ppu -C +es -O +Obb1000 +FPVZOU
  LDFLAGS=-Wl,-a archive
  OTHERLIBS =
endif

ifeq ($(MYOS),ULTRIX)
  FFLAGS=-check_bounds
  LDFLAGS=
  OTHERLIBS =
endif

ifeq ($(MYOS),OSF1)
  FFLAGS=-O -C -extend_source -fpe -warn unused
  LDFLAGS=
  OTHERLIBS =
endif

#ifeq ($(MYOS),Linux)  
#   F77=g77
#   FFLAGS=-g -ffixed-line-length-132
##   FFLAGS=-g -I$(Csoft)/INCLUDE -ffixed-line-length-132
#   DISPFLAGS=$(FFLAGS)
#   OTHERLIBS += -lc -lm -lnsl
#   OURLIBS := $(OURGENLIBS) $(LIBROOT)/libport.a
#endif

ifeq ($(MYOS),Linux)
# Uncomment these two lines for redhat 7.2
#  ABSOFT=/apps/absoft/PRO/usr/absoft
#  CERN_ROOT=/apps/cernlib/i386_redhat72/2001
# Uncomment these two lines for Enterprise Linux 3.
#  ABSOFT=/apps/absoft/PRO/opt/absoft
  ABSOFT=/apps/absoft/absoft-8.2/opt/absoft
  CERN_ROOT = /apps/cernlib/i386_rhel3/2003
#  FABSFLAGS=-O -V -W -f -s -N1 -B108 -B100 -N90 -N22 -N2 -N113
  FABSFLAGS=-O -V -W -f -s -N1 -B108 -B100 -N90 -N22 -N113
  INCLUDES=-I.,..,MEC/ 
  EXTRAFLAGS=-DABSOFTFORTRAN
  FFLAGS= $(INCLUDES) $(FABSFLAGS) $(EXTRAFLAGS)
  FFLAG1=$(FFLAGS) -c
  OTHERLIBS = -L$(LIBROOT) -lctp \
        -L$(CERN_ROOT)/lib $(CERNLIBS) -lV77 -lU77 -lg2c -lc -lm \
#	-lnsl -lcrypt
# Uncomment the last line above (-lnsl -lcrypt) if you are using cernlib 2001.
  FC  := $(ABSOFT)/bin/f77
  F77 :=$(ABSOFT)/bin/f77
endif

none: dummy_ana

all: dummy_ana

dummy_ana: $(my_objs) Makefile $(OURLIBS)
	$(F77) -g -o $@ $(FFLAGS) $(my_objs) $(OTHERLIBS) $(CERNLIBS)

clean: 
	$(RM) *.o MEC/*.o dummy_ana
