# Install script for directory: /home/anhvt89/Documents/c++lib/cxsc-2-5-4

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cxsc_mpicomm.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/RtsTyp.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cxsc_mpicomm_templ.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/compiler.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rtsrmath.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/RtsFunc.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/r_ari.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/ivecrmat.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/sparsevector.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lx_interval.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cmatrix.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/liveclrmat.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/imatrix.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cxscmatr.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/dot.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_rmatrix.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/dot_defs.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/sparsematrix.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/scvector.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cvector.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rmatrix.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/sparseidot.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/civector.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cinterval.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/except.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_rvector.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/util.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cxsc_mpicomm.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/ioflags.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/intmatrix.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lrvecimat.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/srvector.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/testcomp.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/interval.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/livecimat.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/sparsedot.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rmath.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lx_complex.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lx_ivector.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/civecrmat.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/livecrmat.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/complex.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rvector.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_rmath.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts_real.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_ivector.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lrvecrmat.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/srmatrix.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/test.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_cimath.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/iveclrmat.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_imath.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cmatimat.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lx_civector.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cimath.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/scivector.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_cinterval.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_complex.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/testintv.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lrvecivec.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/docu.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/iveccmat.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/real.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/vector.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_interval.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lx_real.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/simatrix.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/idot.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/testdot.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/civeccmat.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cdot.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_cmath.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cvecimat.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lx_cinterval.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cxsc_blas.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_real.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/testvect.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/civecimat.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/iveccvec.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/testclss.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/testmatr.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/ivector.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cimatrix.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/imath.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/scimatrix.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/xscclass.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cvecrmat.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_imatrix.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cxscvect.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/sivector.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/sparsecidot.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/scmatrix.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/matrix.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/intvector.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/testsklr.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/sparsecdot.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cidot.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lrmatimat.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lx_real.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/civecrmat.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/vector_friend_declarations.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cinterval.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/complex.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_ivector.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/liveclrmat.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/sparsevector.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lrvecimat.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cdot.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_rmath.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/idotk.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/sparsedot.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rmatrix.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/iveccvec.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/real.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/livecimat.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cmatrix.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cmatimat.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lx_complex.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/imath.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/matrix_friend_declarations.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_cinterval.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_complex.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cvecrmat.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/dotk.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cvector.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lrvecrmat.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/civector.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cxsc_blas.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lrvecivec.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_rvector.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/matrix.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/intvector.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_interval.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cidotk.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_rmatrix.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/intmatrix.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/civecimat.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/livecrmat.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lrmatimat.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/ivector.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cvecimat.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/imatrix.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lx_ivector.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rmath.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cimatrix.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/sparsecdot.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/civeccmat.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/vector.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lx_cinterval.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/sparseidot.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/iveccmat.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cidot.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/cdotk.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lx_interval.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_cmath.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/sparsematrix.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rvector.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/l_imatrix.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/ivecrmat.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/idot.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/iveclrmat.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/lx_civector.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/interval.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/b_64bt.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/o_defs.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/b_defs.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/s_fcth.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/t_exc_.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/o_revs.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/t_ieee.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/y_fcth.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/l_fcth.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/o_spec.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/o_syst.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/a_defs.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/p88rts.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/addbody.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/subbody.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/b_fcth.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/o_slct.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/o_type.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/body.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/f_fcth.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/t_cnst.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/t_drea.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/d_fcth.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/t_cond.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/divtrap.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/r_fcth.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/o_msg1.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/divbody.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/o_name.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/b_lpi_.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/d_defs.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/t_name.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/t_fcth.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/l_defs.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/a_fcth.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/r_defs.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/e_fcth.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/o_fcth.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/t_defs.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/mulbody.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/t_ddev.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/b_lari.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/e_defs.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/rts/t_mach.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/asm/r_ari.h"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/fi_lib/fi_lib_consts.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/src/fi_lib/fi_lib.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/linsys.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/lss.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/cipoly.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/ddf_ari.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/clss.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/rev_simp.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/i_util.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/mvi_util.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/gop.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/rpeval.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/expreval.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/r_util.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/stacksz.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/matinv_aprx.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/cpzero.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/lst_ari.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/gop1.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/fastlss.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/xi_ari.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/nlfzero.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/lop.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/lst1_ari.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/lop_ari.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/grad_ari.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/ci_util.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/rpoly.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/ilss.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/mv_util.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/cilss.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/cpoly.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/hess_ari.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/nlinsys.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/matinv.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/set_ari.hpp"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/fastlss.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Modules/ddf_ari.inl"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscconf.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/doc" TYPE FILE FILES
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/changelog"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/README"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcxsc.so.2.5.4"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcxsc.so.2"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcxsc.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/libcxsc.so.2.5.4"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/libcxsc.so.2"
    "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/libcxsc.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcxsc.so.2.5.4"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcxsc.so.2"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libcxsc.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH ":::::::::::::::::::::::::::::::::::::::::::::::::::::::"
           NEW_RPATH "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/cxscbuild/lib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/anhvt89/Documents/c++lib/cxsc-2-5-4/examples/cmake_install.cmake")
  include("/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/Programs/cmake_install.cmake")
  include("/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/test/cmake_install.cmake")
  include("/home/anhvt89/Documents/c++lib/cxsc-2-5-4/CToolbox/test/test_programms/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/anhvt89/Documents/c++lib/cxsc-2-5-4/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
