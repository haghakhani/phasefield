#
# automake template for src/tecplot directory
#
# $Id: Makefile.am 211 2009-06-16 20:02:10Z dkumar $
#

noinst_LIBRARIES = libtecplot.a

libtecplot_a_SOURCES = outsum.C tecplot.C grassout.C

if HAVE_HDF5
   libtecplot_a_SOURCES += xdfm_write.C GMFG_hdfconstants.h hd5calls.c
endif

INCLUDES=-I../header/
