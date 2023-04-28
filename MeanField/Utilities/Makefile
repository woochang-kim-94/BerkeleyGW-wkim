# BerkleyGW: Meanfield/Utilities Makefile
#
# Most rules are defined in Common/common-rules.mk, and dependencies in
# make_deps.mk. Rerun dev-scripts/make_deps.py if you change any dependency!

PREFIX=../..
include $(PREFIX)/Common/common-rules.mk

INTERNAL =
ifeq ($(findstring -DHDF5,$(MATHFLAG)),-DHDF5)
  INTERNAL += wfn2hdf.x hdf2wfn.x
endif
ifeq ($(TYPEFLAG),-DCPLX)
  INTERNAL += wfnmix_spinor.x wfnmix_spinor_rho_vxc.x
endif
#BEGIN_INTERNAL_ONLY
INTERNAL += analyzebz.x wfnreduce.x wfnmix_QSGW$(FLAVOR).x real2complex.x
ifeq ($(findstring -DHDF5,$(MATHFLAG)),-DHDF5)
INTERNAL += compute_vcoul.x
endif
ifeq ($(findstring -DHDF5,$(MATHFLAG)),-DHDF5)
  INTERNAL += extract_vcoul.x
endif
#END_INTERNAL_ONLY

LIBS = $(TILE_LIBS_EXT) $(SCALAPACKLIB) $(LAPACKLIB) $(FFTWLIB) $(HDF5LIB) $(PERFORMANCE)

default: all

all: shell_scripts mf_convert.x degeneracy_check.x wfnmerge.x \
     wfn_rho_vxc_info.x scissors2eqp.x wfn_dotproduct.x fix_occ.x \
     wfn_time_reversal.x wfn_modify_spinor.x  $(INTERNAL)

# FHJ: Extra dependency
fix_occ.x: fix_occ_cp.o

# FHJ: NOTE: do not put any explicit dependecy here. Simply run make_deps.py
# and commit the newly generated make_deps.mk. End users will never have to
# do this step of generating make_deps.mk
include make_deps.mk
