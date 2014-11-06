#########################
#		        #
#  Makefile for PFem3d  #
#                       #
#########################

####################################################################
# Directories 
#

include ../Makefile.path


####################################################################
# Compiler parameters
#

CC = mpicc
LINK = mpicxx
DFLAGS		= #-DNDEBUG

# BG/Q machines
# CC = mpixlc_r
# LINK = mpixlf77_r
# DFLAGS	+= -DARCH_BGQ

#TAU_OPTS = -optKeepFiles -tau_options=-optTauSelectFile=tau_instrument.txt
#CC = tau_cc.sh $(TAU_OPTS)
#LINK = tau_cxx.sh $(TAU_OPTS)

# Un/comment the appropriate CFLAGS variable for using intel or GNU
# compiler for correct strict warnings
CFLAGS_INTEL	= -Wcheck
#CFLAGS_GCC	= -pedantic -Werror
CFLAGS          = -Wall -Werror -std=c99 -O3 -g $(CFLAGS_INTEL) $(CFLAGS_GCC)

#DFLAGS		+= -DNO_VTK_LIB=1
#DFLAGS		+= -DDISP_DEBUG=1
#DFLAGS		+= -DVD_INT_ALG_DEBUG=1
#DFLAGS		+= -DSTAB_DEBUG=1
#DFLAGS		+= -DUL_DEBUG=1
#DFLAGS		+= -DVD_DEBUG=1
#DFLAGS		+= -DARC_DEBUG=1
#DFLAGS		+= -DDOF_DEBUG=1
#DFLAGS		+= -DUTILS_DEBUG=1
#DFLAGS		+= -DNR_UPDATE=1
#DFLAGS		+= -DARC_UPDATE=1
#DFLAGS		+= -DNO_BUBBLE=1
#DFLAGS		+= -DMINI_DEBUG=1
#DFLAGS		+= -DMINI_P1_P1=1
#DFLAGS		+= -DMINI_3F_DEBUG=1
#DFLAGS		+= -DMINI_3F_P1_P1=1
#DFLAGS		+= -DHW_DEBUG=1
#DFLAGS		+= -DSTAB_DEBUG=1
#DFLAGS		+= -DPFEM_DEBUG_GRN=1
#DFLAGS		+= -DNEW_ASSEM=1
#DFLAGS		+= -DPFEM_TIMERS=1
#DFLAGS		+= -DPFEM_TRACE=1
#DFLAGS		+= -DORIGINAL_STIFFMAT_FD=1
#DFLAGS		+= -DPFEM_DEBUG=1
#DFLAGS		+= -DPFEM_DEBUG_ALL=1
#DFLAGS		+= -DDEBUG_STIFFNESS=1
#DFLAGS		+= -DPFEM_PRINT=1
#DFLAGS		+= -DPFEM_PRINT_MATRIX=1

CFLAGS		+= $(DFLAGS)

LDFLAGS		=
LDOPTIONS       = $(CFLAGS) $(LDFLAGS)

####################################################################
# Includes
#
INCLUDES        = -I$(pfem_incl_dir) \
		  $(RENUMBERING_incl) \
		  $(VTK_IO_incl) \
		  $(PARMETIS_incl) \
		  $(SUITESPARSE_incl) \
		  $(BLOCKSOLVE_incl) \
		  $(HYPRE_incl) \
		  $(VTK_incl) \
		  $(MKL_incl)

####################################################################
# Libraries
#
LIBS	= $(RENUMBERING_lib) \
	  $(VTK_IO_lib) \
	  $(PARMETIS_lib) \
	  $(SUITESPARSE_lib) \
	  $(BLOCKSOLVE_lib) \
	  $(HYPRE_lib) \
	  $(VTK_lib) \
	  $(MKL_lib) \
	  $(MATH_lib)

# uncomment to link with totalview
LIBS	+= $(TOTALVIEW_lib)

####################################################################
# List of object files
#

include Makefile.objects
SRCS = $(OBJS:.o=.c)

####################################################################
# Probably will not need to change anything below here
RM              = rm -rf

vpath %.c $(srcdir)

####################################################################
#  Targets
#

all: PGFem3D

PGFem3D: $(OBJS)
	$(RM) $(pfem_bin_dir)/$@
	$(LINK) -o $(pfem_bin_dir)/$@ $(LDOPTIONS) $(OBJS) $(LIBS)

objects: $(OBJS)

%.o: %.c
	$(CC) -c $(INCLUDES) $(CFLAGS) $<

clean:
	$(RM)  *.o *~

realclean: clean
	$(RM) $(pfem_bin_dir)/PGFem3D

docs:
	doxygen docs.config

tags:
	etags $(SRCS)

depend:
	makedepend -f- -- $(INCLUDES) $(LOC_INCL) $(CFLAGS) -- $(SRCS) > Makefile.depend

# DO NOT DELETE
-include Makefile.depend
