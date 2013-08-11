#ifndef _MHD_H_
#define _MHD_H_

#if 0
!##############################################################################
!# ****************************************************************************
!# <name> mhd </name>
!# ****************************************************************************
!#
!# <purpose>
!# This is the magnetohydrodynamic header file
!#
!# </purpose>
!##############################################################################
#endif

#include "flagship.h"
#include "kernel/System/fmath.h"
#include "kernel/System/idxmanager.h"


#if 0
!##############################################################################
! Number of variables
!##############################################################################
#endif

#define NVAR1D 7
#define NVAR2D 8
#define NVAR3D 8


#if 0
!##############################################################################
! Global type of system format
!##############################################################################
#endif

#define SYSTEM_SCALAR 0
#define SYSTEM_BLOCK  1


#if 0
!##############################################################################
! Global type of system coupling approach
!##############################################################################
#endif

#define SYSTEM_SEGREGATED 0
#define SYSTEM_ALLCOUPLED 1


#if 0
!##############################################################################
! Global type of dissipation
!##############################################################################
#endif

#define DISSIPATION_ZERO            0
#define DISSIPATION_SCALAR          1
#define DISSIPATION_ROE             2
#define DISSIPATION_RUSANOV         3
#define DISSIPATION_SCALAR_DSPLIT  -DISSIPATION_SCALAR
#define DISSIPATION_ROE_DSPLIT     -DISSIPATION_ROE
#define DISSIPATION_RUSANOV_DSPLIT -DISSIPATION_RUSANOV


#if 0
!##############################################################################
! Position of variables in vector of conservative variables
!
! By setting MHD_VARPOS_EXTERNAL it is possible to define the position of
! the conservative variables externally, e.g., in a larger solution vector.
!##############################################################################
#endif

#ifndef MHD_VARPOS_EXTERNAL

#define _MHD_DENSITY_1D_     1
#define _MHD_XMOMENTUM_1D_   2
#define _MHD_YMOMENTUM_1D_   3
#define _MHD_ZMOMENTUM_1D_   4
#undef  _MHD_XMAGFIELD_1D_
#define _MHD_YMAGFIELD_1D_   5
#define _MHD_ZMAGFIELD_1D_   6
#define _MHD_TOTALENERGY_1D_ 7

#define _MHD_DENSITY_2D_     1
#define _MHD_XMOMENTUM_2D_   2
#define _MHD_YMOMENTUM_2D_   3
#define _MHD_ZMOMENTUM_2D_   4
#define _MHD_XMAGFIELD_2D_   5
#define _MHD_YMAGFIELD_2D_   6
#define _MHD_ZMAGFIELD_2D_   7
#define _MHD_TOTALENERGY_2D_ 8

#define _MHD_DENSITY_3D_     1
#define _MHD_XMOMENTUM_3D_   2
#define _MHD_YMOMENTUM_3D_   3
#define _MHD_ZMOMENTUM_3D_   4
#define _MHD_XMAGFIELD_3D_   5
#define _MHD_YMAGFIELD_3D_   6
#define _MHD_ZMAGFIELD_3D_   7
#define _MHD_TOTALENERGY_3D_ 8

#endif


#if 0
! In one dimension, the x-component of the magnetic field is constant.
#endif
#ifndef MHD_XMAGFIELD_CONST
#define MHD_XMAGFIELD_CONST 1.0
#endif


#if 0
!##############################################################################
! Conservative variables
!##############################################################################
#endif

#define DENSITY1_1D(U,IDX)                  IDX(U,_MHD_DENSITY_1D_)
#define DENSITY2_1D(U,IDX,i,n1,n2)          IDX(U,_MHD_DENSITY_1D_,i,n1,n2)
#define DENSITY3_1D(U,IDX,i,j,n1,n2,n3)     IDX(U,_MHD_DENSITY_1D_,i,j,n1,n2,n3)

#define DENSITY1_2D(U,IDX)                  IDX(U,_MHD_DENSITY_2D_)
#define DENSITY2_2D(U,IDX,i,n1,n2)          IDX(U,_MHD_DENSITY_2D_,i,n1,n2)
#define DENSITY3_2D(U,IDX,i,j,n1,n2,n3)     IDX(U,_MHD_DENSITY_2D_,i,j,n1,n2,n3)

#define DENSITY1_3D(U,IDX)                  IDX(U,_MHD_DENSITY_3D_)
#define DENSITY2_3D(U,IDX,i,n1,n2)          IDX(U,_MHD_DENSITY_3D_,i,n1,n2)
#define DENSITY3_3D(U,IDX,i,j,n1,n2,n3)     IDX(U,_MHD_DENSITY_3D_,i,j,n1,n2,n3)

#define XMOMENTUM1_1D(U,IDX)                IDX(U,_MHD_XMOMENTUM_1D_)
#define XMOMENTUM2_1D(U,IDX,i,n1,n2)        IDX(U,_MHD_XMOMENTUM_1D_,i,n1,n2)
#define XMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_XMOMENTUM_1D_,i,j,n1,n2,n3)

#define XMOMENTUM1_2D(U,IDX)                IDX(U,_MHD_XMOMENTUM_2D_)
#define XMOMENTUM2_2D(U,IDX,i,n1,n2)        IDX(U,_MHD_XMOMENTUM_2D_,i,n1,n2)
#define XMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_XMOMENTUM_2D_,i,j,n1,n2,n3)

#define XMOMENTUM1_3D(U,IDX)                IDX(U,_MHD_XMOMENTUM_3D_)
#define XMOMENTUM2_3D(U,IDX,i,n1,n2)        IDX(U,_MHD_XMOMENTUM_3D_,i,n1,n2)
#define XMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_XMOMENTUM_3D_,i,j,n1,n2,n3)
      
#define YMOMENTUM1_1D(U,IDX)                IDX(U,_MHD_YMOMENTUM_1D_)
#define YMOMENTUM2_1D(U,IDX,i,n1,n2)        IDX(U,_MHD_YMOMENTUM_1D_,i,n1,n2)
#define YMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_YMOMENTUM_1D_,i,j,n1,n2,n3)

#define YMOMENTUM1_2D(U,IDX)                IDX(U,_MHD_YMOMENTUM_2D_)
#define YMOMENTUM2_2D(U,IDX,i,n1,n2)        IDX(U,_MHD_YMOMENTUM_2D_,i,n1,n2)
#define YMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_YMOMENTUM_2D_,i,j,n1,n2,n3)

#define YMOMENTUM1_3D(U,IDX)                IDX(U,_MHD_YMOMENTUM_3D_)
#define YMOMENTUM2_3D(U,IDX,i,n1,n2)        IDX(U,_MHD_YMOMENTUM_3D_,i,n1,n2)
#define YMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_YMOMENTUM_3D_,i,j,n1,n2,n3)

#define ZMOMENTUM1_1D(U,IDX)                IDX(U,_MHD_ZMOMENTUM_1D_)
#define ZMOMENTUM2_1D(U,IDX,i,n1,n2)        IDX(U,_MHD_ZMOMENTUM_1D_,i,n1,n2)
#define ZMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_ZMOMENTUM_1D_,i,j,n1,n2,n3)

#define ZMOMENTUM1_2D(U,IDX)                IDX(U,_MHD_ZMOMENTUM_2D_)
#define ZMOMENTUM2_2D(U,IDX,i,n1,n2)        IDX(U,_MHD_ZMOMENTUM_2D_,i,n1,n2)
#define ZMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_ZMOMENTUM_2D_,i,j,n1,n2,n3)

#define ZMOMENTUM1_3D(U,IDX)                IDX(U,_MHD_ZMOMENTUM_3D_)
#define ZMOMENTUM2_3D(U,IDX,i,n1,n2)        IDX(U,_MHD_ZMOMENTUM_3D_,i,n1,n2)
#define ZMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_ZMOMENTUM_3D_,i,j,n1,n2,n3)

#define XMAGFIELD1_1D(U,IDX)                (DCONST(MHD_XMAGFIELD_CONST))
#define XMAGFIELD2_1D(U,IDX,i,n1,n2)        (DCONST(MHD_XMAGFIELD_CONST))
#define XMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3)   (DCONST(MHD_XMAGFIELD_CONST))

#define XMAGFIELD1_2D(U,IDX)                IDX(U,_MHD_XMAGFIELD_2D_)
#define XMAGFIELD2_2D(U,IDX,i,n1,n2)        IDX(U,_MHD_XMAGFIELD_2D_,i,n1,n2)
#define XMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_XMAGFIELD_2D_,i,j,n1,n2,n3)

#define XMAGFIELD1_3D(U,IDX)                IDX(U,_MHD_XMAGFIELD_3D_)
#define XMAGFIELD2_3D(U,IDX,i,n1,n2)        IDX(U,_MHD_XMAGFIELD_3D_,i,n1,n2)
#define XMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_XMAGFIELD_3D_,i,j,n1,n2,n3)

#define YMAGFIELD1_1D(U,IDX)                IDX(U,_MHD_YMAGFIELD_1D_)
#define YMAGFIELD2_1D(U,IDX,i,n1,n2)        IDX(U,_MHD_YMAGFIELD_1D_,i,n1,n2)
#define YMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_YMAGFIELD_1D_,i,j,n1,n2,n3)

#define YMAGFIELD1_2D(U,IDX)                IDX(U,_MHD_YMAGFIELD_2D_)
#define YMAGFIELD2_2D(U,IDX,i,n1,n2)        IDX(U,_MHD_YMAGFIELD_2D_,i,n1,n2)
#define YMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_YMAGFIELD_2D_,i,j,n1,n2,n3)

#define YMAGFIELD1_3D(U,IDX)                IDX(U,_MHD_YMAGFIELD_3D_)
#define YMAGFIELD2_3D(U,IDX,i,n1,n2)        IDX(U,_MHD_YMAGFIELD_3D_,i,n1,n2)
#define YMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_YMAGFIELD_3D_,i,j,n1,n2,n3)

#define ZMAGFIELD1_1D(U,IDX)                IDX(U,_MHD_ZMAGFIELD_1D_)
#define ZMAGFIELD2_1D(U,IDX,i,n1,n2)        IDX(U,_MHD_ZMAGFIELD_1D_,i,n1,n2)
#define ZMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_ZMAGFIELD_1D_,i,j,n1,n2,n3)

#define ZMAGFIELD1_2D(U,IDX)                IDX(U,_MHD_ZMAGFIELD_2D_)
#define ZMAGFIELD2_2D(U,IDX,i,n1,n2)        IDX(U,_MHD_ZMAGFIELD_2D_,i,n1,n2)
#define ZMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_ZMAGFIELD_2D_,i,j,n1,n2,n3)

#define ZMAGFIELD1_3D(U,IDX)                IDX(U,_MHD_ZMAGFIELD_3D_)
#define ZMAGFIELD2_3D(U,IDX,i,n1,n2)        IDX(U,_MHD_ZMAGFIELD_3D_,i,n1,n2)
#define ZMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3)   IDX(U,_MHD_ZMAGFIELD_3D_,i,j,n1,n2,n3)

#define TOTALENERGY1_1D(U,IDX)              IDX(U,_MHD_TOTALENERGY_1D_)
#define TOTALENERGY2_1D(U,IDX,i,n1,n2)      IDX(U,_MHD_TOTALENERGY_1D_,i,n1,n2)
#define TOTALENERGY3_1D(U,IDX,i,j,n1,n2,n3) IDX(U,_MHD_TOTALENERGY_1D_,i,j,n1,n2,n3)

#define TOTALENERGY1_2D(U,IDX)              IDX(U,_MHD_TOTALENERGY_2D_)
#define TOTALENERGY2_2D(U,IDX,i,n1,n2)      IDX(U,_MHD_TOTALENERGY_2D_,i,n1,n2)
#define TOTALENERGY3_2D(U,IDX,i,j,n1,n2,n3) IDX(U,_MHD_TOTALENERGY_2D_,i,j,n1,n2,n3)

#define TOTALENERGY1_3D(U,IDX)              IDX(U,_MHD_TOTALENERGY_3D_)
#define TOTALENERGY2_3D(U,IDX,i,n1,n2)      IDX(U,_MHD_TOTALENERGY_3D_,i,n1,n2)
#define TOTALENERGY3_3D(U,IDX,i,j,n1,n2,n3) IDX(U,_MHD_TOTALENERGY_3D_,i,j,n1,n2,n3)


#if 0
!##############################################################################
! Roe average values
!##############################################################################
#endif

#define ROE_MEAN_RATIO(ul,ur)       (sqrt(ul/ur))
#define ROE_MEAN_VALUE(ul,ur,ratio) ((ratio*ul+MYNEWLINE \
      ur)/MYNEWLINE \
      (ratio+DCONST(1.0)))


#if 0
!##############################################################################
! Auxiliary quantities
!##############################################################################
#endif

#define MAG_DOT_VEL1_1D(U,IDX,u,v,w)              (XMAGFIELD1_1D(U,IDX)*u+MYNEWLINE \
      YMAGFIELD1_1D(U,IDX)*v+MYNEWLINE \
      ZMAGFIELD1_1D(U,IDX)*w)
#define MAG_DOT_VEL2_1D(U,IDX,i,n1,n2,u,v,w)      (XMAGFIELD2_1D(U,IDX,i,n1,n2)*u+MYNEWLINE \
      YMAGFIELD2_1D(U,IDX,i,n1,n2)*v+MYNEWLINE \
      ZMAGFIELD2_1D(U,IDX,i,n1,n2)*w)
#define MAG_DOT_VEL3_1D(U,IDX,i,j,n1,n2,n3,u,v,w) (XMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3)*u+MYNEWLINE \
      YMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3)*v+MYNEWLINE \
      ZMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3)*w)

#define MAG_DOT_VEL1_2D(U,IDX,u,v,w)              (XMAGFIELD1_2D(U,IDX)*u+MYNEWLINE \
      YMAGFIELD1_2D(U,IDX)*v+MYNEWLINE \
      ZMAGFIELD1_2D(U,IDX)*w)
#define MAG_DOT_VEL2_2D(U,IDX,i,n1,n2,u,v,w)      (XMAGFIELD2_2D(U,IDX,i,n1,n2)*u+MYNEWLINE \
      YMAGFIELD2_2D(U,IDX,i,n1,n2)*v+MYNEWLINE \
      ZMAGFIELD2_2D(U,IDX,i,n1,n2)*w)
#define MAG_DOT_VEL3_2D(U,IDX,i,j,n1,n2,n3,u,v,w) (XMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*u+MYNEWLINE \
      YMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*v+MYNEWLINE \
      ZMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*w)

#define MAG_DOT_VEL1_3D(U,IDX,u,v,w)              (XMAGFIELD1_3D(U,IDX)*u+MYNEWLINE \
      YMAGFIELD1_3D(U,IDX)*v+MYNEWLINE \
      ZMAGFIELD1_3D(U,IDX)*w)
#define MAG_DOT_VEL2_3D(U,IDX,i,n1,n2,u,v,w)      (XMAGFIELD2_3D(U,IDX,i,n1,n2)*u+MYNEWLINE \
      YMAGFIELD2_3D(U,IDX,i,n1,n2)*v+MYNEWLINE \
      ZMAGFIELD2_3D(U,IDX,i,n1,n2)*w)
#define MAG_DOT_VEL3_3D(U,IDX,i,j,n1,n2,n3,u,v,w) (XMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3)*u+MYNEWLINE \
      YMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3)*v+MYNEWLINE \
      ZMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3)*w)


#if 0
!##############################################################################
! Include magnetohydrodynamic header file
!##############################################################################
#endif

#if 0
! Specify ratio of specific heats
#endif
#ifdef MHD_GAMMA
#define MAGNETOHYDRODYN_GAMMA MHD_GAMMA
#endif

#if 0
! Specify Vacuum permittivity
#endif
#ifdef MHD_VACUUM_PERM
#define MAGNETOHYDRODYN_VACUUM_PERM MHD_VACUUM_PERM
#endif

#include "kernel/magnetohydrodynamics.h"


#if 0
!##############################################################################
! Fluxes and matrices in 1D
!##############################################################################
#endif

#if 0
! Flux in x-direction for ideal MHD in 1D
#endif

#define INVISCIDFLUX1_XDIR1_1D(U,IDX,u,v,w,p,q)  XMOMENTUM1_1D(U,IDX)
#define INVISCIDFLUX2_XDIR1_1D(U,IDX,u,v,w,p,q) (XMOMENTUM1_1D(U,IDX)*u+p)
#define INVISCIDFLUX3_XDIR1_1D(U,IDX,u,v,w,p,q) (YMOMENTUM1_1D(U,IDX)*u-MYNEWLINE \
      XMAGFIELD1_1D(U,IDX)*MYNEWLINE \
      YMAGFIELD(U,IDX))
#define INVISCIDFLUX4_XDIR1_1D(U,IDX,u,v,w,p,q) (ZMOMENTUM1_1D(U,IDX)*u-MYNEWLINE \
      XMAGFIELD1_1D(U,IDX)*MYNEWLINE \
      ZMAGFIELD1_1D(U,IDX))
#define INVISCIDFLUX5_XDIR1_1D(U,IDX,u,v,w,p,q) (YMAGFIELD1_1D(U,IDX)*u-MYNEWLINE \
      XMAGFIELD1_1D(U,IDX)*v)
#define INVISCIDFLUX6_XDIR1_1D(U,IDX,u,v,w,p,q) (ZMAGFIELD1_1D(U,IDX)*u-MYNEWLINE \
      XMAGFIELD1_1D(U,IDX)*w)
#define INVISCIDFLUX7_XDIR1_1D(U,IDX,u,v,w,p,q) (TOTALENERGY1_1D(U,IDX)*u+p*u-MYNEWLINE \
      XMAGFIELD1_1D(U,IDX)*q)

#define INVISCIDFLUX1_XDIR2_1D(U,IDX,i,n1,n2,u,v,w,p,q)  XMOMENTUM2_1D(U,IDX,i,n1,n2)
#define INVISCIDFLUX2_XDIR2_1D(U,IDX,i,n1,n2,u,v,w,p,q) (XMOMENTUM2_1D(U,IDX,i,n1,n2)*u+p)
#define INVISCIDFLUX3_XDIR2_1D(U,IDX,i,n1,n2,u,v,w,p,q) (YMOMENTUM2_1D(U,IDX,i,n1,n2)*u-MYNEWLINE \
      XMAGFIELD2_1D(U,IDX,i,n1,n2)*MYNEWLINE \
      YMAGFIELD2_1D(U,IDX,i,n1,n2))
#define INVISCIDFLUX4_XDIR2_1D(U,IDX,i,n1,n2,u,v,w,p,q) (ZMOMENTUM2_1D(U,IDX,i,n1,n2)*u-MYNEWLINE \
      XMAGFIELD2_1D(U,IDX,i,n1,n2)*MYNEWLINE \
      ZMAGFIELD2_1D(U,IDX,i,n1,n2))
#define INVISCIDFLUX5_XDIR2_1D(U,IDX,i,n1,n2,u,v,w,p,q) (YMAGFIELD2_1D(U,IDX,i,n1,n2)*u-MYNEWLINE \
      XMAGFIELD2_1D(U,IDX,i,n1,n2)*v)
#define INVISCIDFLUX6_XDIR2_1D(U,IDX,i,n1,n2,u,v,w,p,q) (ZMAGFIELD2_1D(U,IDX,i,n1,n2)*u-MYNEWLINE \
      XMAGFIELD2_1D(U,IDX,i,n1,n2)*w)
#define INVISCIDFLUX7_XDIR2_1D(U,IDX,i,n1,n2,u,v,w,p,q) (TOTALENERGY2_1D(U,IDX,i,n1,n2)*u+p*u-MYNEWLINE \
      XMAGFIELD2_1D(U,IDX,i,n1,n2)*q)

#define INVISCIDFLUX1_XDIR3_1D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q)  XMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3)
#define INVISCIDFLUX2_XDIR3_1D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (XMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3)*u+p)
#define INVISCIDFLUX3_XDIR3_1D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (YMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3)*u-MYNEWLINE \
      XMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3)*MYNEWLINE \
      YMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3))
#define INVISCIDFLUX4_XDIR3_1D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (ZMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3)*u-MYNEWLINE \
      XMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3)*MYNEWLINE \
      ZMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3))
#define INVISCIDFLUX5_XDIR3_1D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (YMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3)*u-MYNEWLINE \
      XMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3)*v)
#define INVISCIDFLUX6_XDIR3_1D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (ZMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3)*u-MYNEWLINE \
      XMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3)*w)
#define INVISCIDFLUX7_XDIR3_1D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (TOTALENERGY3_1D(U,IDX,i,j,n1,n2,n3)*u+p*u-MYNEWLINE \
      XMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3)*q)

#if 0
! Flux Jacobian matrix for ideal MHD in 1D
#endif

#if 0
! Flux Jacobian matrix in x-direction for ideal MHD in 1D
#endif


#if 0
!##############################################################################
! Fluxes and matrices in 2D
!##############################################################################
#endif

#if 0
! Flux in x-direction for ideal MHD in 2D
#endif

#define INVISCIDFLUX1_XDIR1_2D(U,IDX,u,v,w,p,q)  XMOMENTUM1_2D(U,IDX)
#define INVISCIDFLUX2_XDIR1_2D(U,IDX,u,v,w,p,q) (XMOMENTUM1_2D(U,IDX)*u+p-MYNEWLINE \
      POW(XMAGFIELD1_2D(U,IDX),2))
#define INVISCIDFLUX3_XDIR1_2D(U,IDX,u,v,w,p,q) (YMOMENTUM1_2D(U,IDX)*u-MYNEWLINE \
      XMAGFIELD1_2D(U,IDX)*MYNEWLINE \
      YMAGFIELD1_2D(U,IDX))
#define INVISCIDFLUX4_XDIR1_2D(U,IDX,u,v,w,p,q) (ZMOMENTUM1_2D(U,IDX)*u-MYNEWLINE \
      XMAGFIELD1_2D(U,IDX)*MYNEWLINE \
      ZMAGFIELD1_2D(U,IDX))
#define INVISCIDFLUX5_XDIR1_2D(U,IDX,u,v,w,p,q)  DCONST(0.0)
#define INVISCIDFLUX6_XDIR1_2D(U,IDX,u,v,w,p,q) (YMAGFIELD1_2D(U,IDX)*u-MYNEWLINE \
      XMAGFIELD1_2D(U,IDX)*v)
#define INVISCIDFLUX7_XDIR1_2D(U,IDX,u,v,w,p,q) (ZMAGFIELD1_2D(U,IDX)*u-MYNEWLINE \
      XMAGFIELD1_2D(U,IDX)*w)
#define INVISCIDFLUX8_XDIR1_2D(U,IDX,u,v,w,p,q) (TOTALENERGY1_2D(U,IDX)*u+p*u-MYNEWLINE \
      XMAGFIELD1_2D(U,IDX)*q)

#define INVISCIDFLUX1_XDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q)  XMOMENTUM2_2D(U,IDX,i,n1,n2)
#define INVISCIDFLUX2_XDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q) (XMOMENTUM2_2D(U,IDX,i,n1,n2)*u+p-MYNEWLINE \
      POW(XMAGFIELD2_2D(U,IDX,i,n1,n2),2))
#define INVISCIDFLUX3_XDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q) (YMOMENTUM2_2D(U,IDX,i,n1,n2)*u-MYNEWLINE \
      XMAGFIELD2_2D(U,IDX,i,n1,n2)*MYNEWLINE \
      YMAGFIELD2_2D(U,IDX,i,n1,n2))
#define INVISCIDFLUX4_XDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q) (ZMOMENTUM2_2D(U,IDX,i,n1,n2)*u-MYNEWLINE \
      XMAGFIELD2_2D(U,IDX,i,n1,n2)*MYNEWLINE \
      ZMAGFIELD2_2D(U,IDX,i,n1,n2))
#define INVISCIDFLUX5_XDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q)  DCONST(0.0)
#define INVISCIDFLUX6_XDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q) (YMAGFIELD2_2D(U,IDX,i,n1,n2)*u-MYNEWLINE \
      XMAGFIELD2_2D(U,IDX,i,n1,n2)*v)
#define INVISCIDFLUX7_XDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q) (ZMAGFIELD2_2D(U,IDX,i,n1,n2)*u-MYNEWLINE \
      XMAGFIELD2_2D(U,IDX,i,n1,n2)*w)
#define INVISCIDFLUX8_XDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q) (TOTALENERGY2_2D(U,IDX,i,n1,n2)*u+p*u-MYNEWLINE \
      XMAGFIELD2_2D(U,IDX,i,n1,n2)*q)

#define INVISCIDFLUX1_XDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q)  XMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)
#define INVISCIDFLUX2_XDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (XMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)*u+p-MYNEWLINE \
      POW(XMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3),2))
#define INVISCIDFLUX3_XDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (YMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)*u-MYNEWLINE \
      XMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*MYNEWLINE \
      YMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3))
#define INVISCIDFLUX4_XDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (ZMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)*u-MYNEWLINE \
      XMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*MYNEWLINE \
      ZMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3))
#define INVISCIDFLUX5_XDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q)  DCONST(0.0)
#define INVISCIDFLUX6_XDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (YMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*u-MYNEWLINE \
      XMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*v)
#define INVISCIDFLUX7_XDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (ZMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*u-MYNEWLINE \
      XMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*w)
#define INVISCIDFLUX8_XDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (TOTALENERGY3_2D(U,IDX,i,j,n1,n2,n3)*u+p*u-MYNEWLINE \
      XMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*q)


#if 0
! Flux in y-direction for ideal MHD in 2D
#endif

#define INVISCIDFLUX1_YDIR1_2D(U,IDX,u,v,w,p,q)  YMOMENTUM1_2D(U,IDX)
#define INVISCIDFLUX2_YDIR1_2D(U,IDX,u,v,w,p,q) (XMOMENTUM1_2D(U,IDX)*v-MYNEWLINE \
      YMAGFIELD1_2D(U,IDX)*MYNEWLINE \
      XMAGFIELD1_2D(U,IDX))
#define INVISCIDFLUX3_YDIR1_2D(U,IDX,u,v,w,p,q) (YMOMENTUM1_2D(U,IDX)*v+p-MYNEWLINE \
      POW(YMAGFIELD1_2D(U,IDX),2))
#define INVISCIDFLUX4_YDIR1_2D(U,IDX,u,v,w,p,q) (ZMOMENTUM1_2D(U,IDX)*v-MYNEWLINE \
      YMAGFIELD1_2D(U,IDX)*MYNEWLINE \
      ZMAGFIELD1_2D(U,IDX))
#define INVISCIDFLUX5_YDIR1_2D(U,IDX,u,v,w,p,q) (XMAGFIELD1_2D(U,IDX)*v-MYNEWLINE \
      YMAGFIELD1_2D(U,IDX)*u)
#define INVISCIDFLUX6_YDIR1_2D(U,IDX,u,v,w,p,q)  DCONST(0.0)
#define INVISCIDFLUX7_YDIR1_2D(U,IDX,u,v,w,p,q) (ZMAGFIELD1_2D(U,IDX)*v-MYNEWLINE \
      YMAGFIELD1_2D(U,IDX)*w)
#define INVISCIDFLUX8_YDIR1_2D(U,IDX,u,v,w,p,q) (TOTALENERGY1_2D(U,IDX)*v+p*v-MYNEWLINE \
      YMAGFIELD1_2D(U,IDX)*q)

#define INVISCIDFLUX1_YDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q)  YMOMENTUM2_2D(U,IDX,i,n1,n2)
#define INVISCIDFLUX2_YDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q) (XMOMENTUM2_2D(U,IDX,i,n1,n2)*v-MYNEWLINE \
      YMAGFIELD2_2D(U,IDX,i,n1,n2)*MYNEWLINE \
      XMAGFIELD2_2D(U,IDX,i,n1,n2))
#define INVISCIDFLUX3_YDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q) (YMOMENTUM2_2D(U,IDX,i,n1,n2)*v+p-MYNEWLINE \
      POW(YMAGFIELD2_2D(U,IDX,i,n1,n2),2))
#define INVISCIDFLUX4_YDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q) (ZMOMENTUM2_2D(U,IDX,i,n1,n2)*v-MYNEWLINE \
      YMAGFIELD2_2D(U,IDX,i,n1,n2)*MYNEWLINE \
      ZMAGFIELD2_2D(U,IDX,i,n1,n2))
#define INVISCIDFLUX5_YDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q) (XMAGFIELD2_2D(U,IDX,i,n1,n2)*v-MYNEWLINE \
      YMAGFIELD2_2D(U,IDX,i,n1,n2)*u)
#define INVISCIDFLUX6_YDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q)  DCONST(0.0)
#define INVISCIDFLUX7_YDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q) (ZMAGFIELD2_2D(U,IDX,i,n1,n2)*v-MYNEWLINE \
      YMAGFIELD2_2D(U,IDX,i,n1,n2)*w)
#define INVISCIDFLUX8_YDIR2_2D(U,IDX,i,n1,n2,u,v,w,p,q) (TOTALENERGY2_2D(U,IDX,i,n1,n2)*v+p*v-MYNEWLINE \
      YMAGFIELD2_2D(U,IDX,i,n1,n2)*q)

#define INVISCIDFLUX1_YDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q)  YMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)
#define INVISCIDFLUX2_YDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (XMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)*v-MYNEWLINE \
      YMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*MYNEWLINE \
      XMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3))
#define INVISCIDFLUX3_YDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (YMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)*v+p-MYNEWLINE \
      POW(YMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3),2))
#define INVISCIDFLUX4_YDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (ZMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)*v-MYNEWLINE \
      YMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*MYNEWLINE \
      ZMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3))
#define INVISCIDFLUX5_YDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (XMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*v-MYNEWLINE \
      YMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*u)
#define INVISCIDFLUX6_YDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q)  DCONST(0.0)
#define INVISCIDFLUX7_YDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (ZMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*v-MYNEWLINE \
      YMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*w)
#define INVISCIDFLUX8_YDIR3_2D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (TOTALENERGY3_2D(U,IDX,i,j,n1,n2,n3)*v+p*v-MYNEWLINE \
      YMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3)*q)


#if 0
!##############################################################################
! Fluxes and matrices in 3D
!##############################################################################
#endif

#if 0
  ! Flux in x-direction for ideal MHD in 3D (same as in 2D)
#endif

#define INVISCIDFLUX1_XDIR1_3D INVISCIDFLUX1_XDIR1_2D
#define INVISCIDFLUX2_XDIR1_3D INVISCIDFLUX2_XDIR1_2D
#define INVISCIDFLUX3_XDIR1_3D INVISCIDFLUX3_XDIR1_2D
#define INVISCIDFLUX4_XDIR1_3D INVISCIDFLUX4_XDIR1_2D
#define INVISCIDFLUX5_XDIR1_3D INVISCIDFLUX5_XDIR1_2D
#define INVISCIDFLUX6_XDIR1_3D INVISCIDFLUX6_XDIR1_2D
#define INVISCIDFLUX7_XDIR1_3D INVISCIDFLUX7_XDIR1_2D
#define INVISCIDFLUX8_XDIR1_3D INVISCIDFLUX8_XDIR1_2D

#define INVISCIDFLUX1_XDIR2_3D INVISCIDFLUX1_XDIR2_2D
#define INVISCIDFLUX2_XDIR2_3D INVISCIDFLUX2_XDIR2_2D
#define INVISCIDFLUX3_XDIR2_3D INVISCIDFLUX3_XDIR2_2D
#define INVISCIDFLUX4_XDIR2_3D INVISCIDFLUX4_XDIR2_2D
#define INVISCIDFLUX5_XDIR2_3D INVISCIDFLUX5_XDIR2_2D
#define INVISCIDFLUX6_XDIR2_3D INVISCIDFLUX6_XDIR2_2D
#define INVISCIDFLUX7_XDIR2_3D INVISCIDFLUX7_XDIR2_2D
#define INVISCIDFLUX8_XDIR2_3D INVISCIDFLUX8_XDIR2_2D

#define INVISCIDFLUX1_XDIR3_3D INVISCIDFLUX1_XDIR3_2D
#define INVISCIDFLUX2_XDIR3_3D INVISCIDFLUX2_XDIR3_2D
#define INVISCIDFLUX3_XDIR3_3D INVISCIDFLUX3_XDIR3_2D
#define INVISCIDFLUX4_XDIR3_3D INVISCIDFLUX4_XDIR3_2D
#define INVISCIDFLUX5_XDIR3_3D INVISCIDFLUX5_XDIR3_2D
#define INVISCIDFLUX6_XDIR3_3D INVISCIDFLUX6_XDIR3_2D
#define INVISCIDFLUX7_XDIR3_3D INVISCIDFLUX7_XDIR3_2D
#define INVISCIDFLUX8_XDIR3_3D INVISCIDFLUX8_XDIR3_2D

#if 0
  ! Flux in y-direction for ideal MHD in 3D (same as in 2D)
#endif

#define INVISCIDFLUX1_YDIR1_3D INVISCIDFLUX1_YDIR1_2D
#define INVISCIDFLUX2_YDIR1_3D INVISCIDFLUX2_YDIR1_2D
#define INVISCIDFLUX3_YDIR1_3D INVISCIDFLUX3_YDIR1_2D
#define INVISCIDFLUX4_YDIR1_3D INVISCIDFLUX4_YDIR1_2D
#define INVISCIDFLUX5_YDIR1_3D INVISCIDFLUX5_YDIR1_2D
#define INVISCIDFLUX6_YDIR1_3D INVISCIDFLUX6_YDIR1_2D
#define INVISCIDFLUX7_YDIR1_3D INVISCIDFLUX7_YDIR1_2D
#define INVISCIDFLUX8_YDIR1_3D INVISCIDFLUX8_YDIR1_2D

#define INVISCIDFLUX1_YDIR2_3D INVISCIDFLUX1_YDIR2_2D
#define INVISCIDFLUX2_YDIR2_3D INVISCIDFLUX2_YDIR2_2D
#define INVISCIDFLUX3_YDIR2_3D INVISCIDFLUX3_YDIR2_2D
#define INVISCIDFLUX4_YDIR2_3D INVISCIDFLUX4_YDIR2_2D
#define INVISCIDFLUX5_YDIR2_3D INVISCIDFLUX5_YDIR2_2D
#define INVISCIDFLUX6_YDIR2_3D INVISCIDFLUX6_YDIR2_2D
#define INVISCIDFLUX7_YDIR2_3D INVISCIDFLUX7_YDIR2_2D
#define INVISCIDFLUX8_YDIR2_3D INVISCIDFLUX8_YDIR2_2D

#define INVISCIDFLUX1_YDIR3_3D INVISCIDFLUX1_YDIR3_2D
#define INVISCIDFLUX2_YDIR3_3D INVISCIDFLUX2_YDIR3_2D
#define INVISCIDFLUX3_YDIR3_3D INVISCIDFLUX3_YDIR3_2D
#define INVISCIDFLUX4_YDIR3_3D INVISCIDFLUX4_YDIR3_2D
#define INVISCIDFLUX5_YDIR3_3D INVISCIDFLUX5_YDIR3_2D
#define INVISCIDFLUX6_YDIR3_3D INVISCIDFLUX6_YDIR3_2D
#define INVISCIDFLUX7_YDIR3_3D INVISCIDFLUX7_YDIR3_2D
#define INVISCIDFLUX8_YDIR3_3D INVISCIDFLUX8_YDIR3_2D

#if 0
  ! Flux in z-direction for ideal MHD in 3D
#endif

#define INVISCIDFLUX1_ZDIR1_3D(U,IDX,u,v,w,p,q)  ZMOMENTUM1_3D(U,IDX)
#define INVISCIDFLUX2_ZDIR1_3D(U,IDX,u,v,w,p,q) (XMOMENTUM1_3D(U,IDX)*w-MYNEWLINE \
      ZMAGFIELD1_3D(U,IDX)*MYNEWLINE \
      XMAGFIELD1_3D(U,IDX))
#define INVISCIDFLUX3_ZDIR1_3D(U,IDX,u,v,w,p,q) (YMOMENTUM1_3D(U,IDX)*w-MYNEWLINE \
      ZMAGFIELD1_3D(U,IDX)*MYNEWLINE \
      YMAGFIELD1_3D(U,IDX))
#define INVISCIDFLUX4_ZDIR1_3D(U,IDX,u,v,w,p,q) (ZMOMENTUM1_3D(U,IDX)*w+p-MYNEWLINE \
      POW(ZMAGFIELD1_3D(U,IDX),2))
#define INVISCIDFLUX5_ZDIR1_3D(U,IDX,u,v,w,p,q) (XMAGFIELD1_3D(U,IDX)*w-MYNEWLINE \
      ZMAGFIELD1_3D(U,IDX)*u)
#define INVISCIDFLUX6_ZDIR1_3D(U,IDX,u,v,w,p,q) (YMAGFIELD1_3D(U,IDX)*w-MYNEWLINE \
      ZMAGFIELD1_3D(U,IDX)*v)
#define INVISCIDFLUX7_ZDIR1_3D(U,IDX,u,v,w,p,q)  DCONST(0.0)
#define INVISCIDFLUX8_ZDIR1_3D(U,IDX,u,v,w,p,q) (TOTALENERGY1_3D(U,IDX)*w+p*w-MYNEWLINE \
      ZMAGFIELD1_3D(U,IDX)*q)

#define INVISCIDFLUX1_ZDIR2_3D(U,IDX,i,n1,n2,u,v,w,p,q)  ZMOMENTUM2_3D(U,IDX,i,n1,n2)
#define INVISCIDFLUX2_ZDIR2_3D(U,IDX,i,n1,n2,u,v,w,p,q) (XMOMENTUM2_3D(U,IDX,i,n1,n2)*w-MYNEWLINE \
      ZMAGFIELD2_3D(U,IDX,i,n1,n2)*MYNEWLINE \
      XMAGFIELD2_3D(U,IDX,i,n1,n2))
#define INVISCIDFLUX3_ZDIR2_3D(U,IDX,i,n1,n2,u,v,w,p,q) (YMOMENTUM2_3D(U,IDX,i,n1,n2)*w-MYNEWLINE \
      ZMAGFIELD2_3D(U,IDX,i,n1,n2)*MYNEWLINE \
      YMAGFIELD2_3D(U,IDX,i,n1,n2))
#define INVISCIDFLUX4_ZDIR2_3D(U,IDX,i,n1,n2,u,v,w,p,q) (ZMOMENTUM2_3D(U,IDX,i,n1,n2)*w+p-MYNEWLINE \
      POW(ZMAGFIELD2_3D(U,IDX,i,n1,n2),2))
#define INVISCIDFLUX5_ZDIR2_3D(U,IDX,i,n1,n2,u,v,w,p,q) (XMAGFIELD2_3D(U,IDX,i,n1,n2)*w-MYNEWLINE \
      ZMAGFIELD2_3D(U,IDX,i,n1,n2)*u)
#define INVISCIDFLUX6_ZDIR2_3D(U,IDX,i,n1,n2,u,v,w,p,q) (YMAGFIELD2_3D(U,IDX,i,n1,n2)*w-MYNEWLINE \
      ZMAGFIELD2_3D(U,IDX,i,n1,n2)*v)
#define INVISCIDFLUX7_ZDIR2_3D(U,IDX,i,n1,n2,u,v,w,p,q)  DCONST(0.0)
#define INVISCIDFLUX8_ZDIR2_3D(U,IDX,i,n1,n2,u,v,w,p,q) (TOTALENERGY2_3D(U,IDX,i,n1,n2)*w+p*w-MYNEWLINE \
      ZMAGFIELD2_3D(U,IDX,i,n1,n2)*q)

#define INVISCIDFLUX1_ZDIR3_3D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q)  ZMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)
#define INVISCIDFLUX2_ZDIR3_3D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (XMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)*w-MYNEWLINE \
      ZMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3)*MYNEWLINE \
      XMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3))
#define INVISCIDFLUX3_ZDIR3_3D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (YMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)*w-MYNEWLINE \
      ZMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3)*MYNEWLINE \
      YMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3))
#define INVISCIDFLUX4_ZDIR3_3D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (ZMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)*w+p-MYNEWLINE \
      POW(ZMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3),2))
#define INVISCIDFLUX5_ZDIR3_3D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (XMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3)*w-MYNEWLINE \
      ZMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3)*u)
#define INVISCIDFLUX6_ZDIR3_3D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (YMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3)*w-MYNEWLINE \
      ZMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3)*v)
#define INVISCIDFLUX7_ZDIR3_3D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q)  DCONST(0.0)
#define INVISCIDFLUX8_ZDIR3_3D(U,IDX,i,j,n1,n2,n3,u,v,w,p,q) (TOTALENERGY3_3D(U,IDX,i,j,n1,n2,n3)*w+p*w-MYNEWLINE \
      ZMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3)*q)

#endif
