#ifndef _HYDRO_H_
#define _HYDRO_H_

#if 0
!##############################################################################
!# ****************************************************************************
!# <name> hydro </name>
!# ****************************************************************************
!#
!# <purpose>
!# This is the hydrodynamic header file
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

#define NVAR1D 3
#define NVAR2D 4
#define NVAR3D 5


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
! By setting HYDRO_VARPOS_EXTERNAL it is possible to define the position of
! the conservative variables externally, e.g., in a larger solution vector.
!##############################################################################
#endif

#ifndef HYDRO_VARPOS_EXTERNAL

#define _HYDRO_DENSITY_1D_     1
#define _HYDRO_XMOMENTUM_1D_   2
#undef  _HYDRO_YMOMENTUM_1D_
#undef  _HYDRO_ZMOMENTUM_1D_
#define _HYDRO_TOTALENERGY_1D_ 3
      
#define _HYDRO_DENSITY_2D_     1
#define _HYDRO_XMOMENTUM_2D_   2
#define _HYDRO_YMOMENTUM_2D_   3
#undef  _HYDRO_ZMOMENTUM_2D_
#define _HYDRO_TOTALENERGY_2D_ 4
      
#define _HYDRO_DENSITY_3D_     1
#define _HYDRO_XMOMENTUM_3D_   2
#define _HYDRO_YMOMENTUM_3D_   3
#define _HYDRO_ZMOMENTUM_3D_   4
#define _HYDRO_TOTALENERGY_3D_ 5
      
#endif

#if 0
!##############################################################################
! Conservative variables
!##############################################################################
#endif

#define DENSITY1_1D(U,IDX)                  IDX(U,_HYDRO_DENSITY_1D_)
#define DENSITY2_1D(U,IDX,i,n1,n2)          IDX(U,_HYDRO_DENSITY_1D_,i,n1,n2)
#define DENSITY3_1D(U,IDX,i,j,n1,n2,n3)     IDX(U,_HYDRO_DENSITY_1D_,i,j,n1,n2,n3)

#define DENSITY1_2D(U,IDX)                  IDX(U,_HYDRO_DENSITY_2D_)
#define DENSITY2_2D(U,IDX,i,n1,n2)          IDX(U,_HYDRO_DENSITY_2D_,i,n1,n2)
#define DENSITY3_2D(U,IDX,i,j,n1,n2,n3)     IDX(U,_HYDRO_DENSITY_2D_,i,j,n1,n2,n3)

#define DENSITY1_3D(U,IDX)                  IDX(U,_HYDRO_DENSITY_3D_)
#define DENSITY2_3D(U,IDX,i,n1,n2)          IDX(U,_HYDRO_DENSITY_3D_,i,n1,n2)
#define DENSITY3_3D(U,IDX,i,j,n1,n2,n3)     IDX(U,_HYDRO_DENSITY_3D_,i,j,n1,n2,n3)
      
#define XMOMENTUM1_1D(U,IDX)                IDX(U,_HYDRO_XMOMENTUM_1D_)
#define XMOMENTUM2_1D(U,IDX,i,n1,n2)        IDX(U,_HYDRO_XMOMENTUM_1D_,i,n1,n2)
#define XMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3)   IDX(U,_HYDRO_XMOMENTUM_1D_,i,j,n1,n2,n3)

#define XMOMENTUM1_2D(U,IDX)                IDX(U,_HYDRO_XMOMENTUM_2D_)
#define XMOMENTUM2_2D(U,IDX,i,n1,n2)        IDX(U,_HYDRO_XMOMENTUM_2D_,i,n1,n2)
#define XMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)   IDX(U,_HYDRO_XMOMENTUM_2D_,i,j,n1,n2,n3)

#define XMOMENTUM1_3D(U,IDX)                IDX(U,_HYDRO_XMOMENTUM_3D_)
#define XMOMENTUM2_3D(U,IDX,i,n1,n2)        IDX(U,_HYDRO_XMOMENTUM_3D_,i,n1,n2)
#define XMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)   IDX(U,_HYDRO_XMOMENTUM_3D_,i,j,n1,n2,n3)

#define YMOMENTUM1_2D(U,IDX)                IDX(U,_HYDRO_YMOMENTUM_2D_)
#define YMOMENTUM2_2D(U,IDX,i,n1,n2)        IDX(U,_HYDRO_YMOMENTUM_2D_,i,n1,n2)
#define YMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)   IDX(U,_HYDRO_YMOMENTUM_2D_,i,j,n1,n2,n3)

#define YMOMENTUM1_3D(U,IDX)                IDX(U,_HYDRO_YMOMENTUM_3D_)
#define YMOMENTUM2_3D(U,IDX,i,n1,n2)        IDX(U,_HYDRO_YMOMENTUM_3D_,i,n1,n2)
#define YMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)   IDX(U,_HYDRO_YMOMENTUM_3D_,i,j,n1,n2,n3)

#define ZMOMENTUM1_3D(U,IDX)                IDX(U,_HYDRO_ZMOMENTUM_3D_)
#define ZMOMENTUM2_3D(U,IDX,i,n1,n2)        IDX(U,_HYDRO_ZMOMENTUM_3D_,i,n1,n2)
#define ZMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)   IDX(U,_HYDRO_ZMOMENTUM_3D_,i,j,n1,n2,n3)

#define TOTALENERGY1_1D(U,IDX)              IDX(U,_HYDRO_TOTALENERGY_1D_)
#define TOTALENERGY2_1D(U,IDX,i,n1,n2)      IDX(U,_HYDRO_TOTALENERGY_1D_,i,n1,n2)
#define TOTALENERGY3_1D(U,IDX,i,j,n1,n2,n3) IDX(U,_HYDRO_TOTALENERGY_1D_,i,j,n1,n2,n3)

#define TOTALENERGY1_2D(U,IDX)              IDX(U,_HYDRO_TOTALENERGY_2D_)
#define TOTALENERGY2_2D(U,IDX,i,n1,n2)      IDX(U,_HYDRO_TOTALENERGY_2D_,i,n1,n2)
#define TOTALENERGY3_2D(U,IDX,i,j,n1,n2,n3) IDX(U,_HYDRO_TOTALENERGY_2D_,i,j,n1,n2,n3)

#define TOTALENERGY1_3D(U,IDX)              IDX(U,_HYDRO_TOTALENERGY_3D_)
#define TOTALENERGY2_3D(U,IDX,i,n1,n2)      IDX(U,_HYDRO_TOTALENERGY_3D_,i,n1,n2)
#define TOTALENERGY3_3D(U,IDX,i,j,n1,n2,n3) IDX(U,_HYDRO_TOTALENERGY_3D_,i,j,n1,n2,n3)


#if 0
!##############################################################################
! Roe average values
!##############################################################################
#endif

#define ROE_MEAN_RATIO(ul,ur)       (sqrt(ul/ur))
#define ROE_MEAN_VALUE(ul,ur,ratio) ((ratio*ul+ur)/(ratio+DCONST(1.0)))


#if 0
!##############################################################################
! Include thermodynamic header file
!##############################################################################
#endif

      
#if 0
! Specify ratio of specific heats
#endif
#ifdef HYDRO_GAMMA
#define THERMODYN_GAMMA HYDRO_GAMMA
#endif

#include "kernel/thermodynamics.h"


#if 0
!##############################################################################
! Fluxes and matrices in 1D
!##############################################################################
#endif

#if 0
! Flux in x-direction for inviscid hydrodynamics in 1D
#endif

#define INVISCIDFLUX1_XDIR1_1D(U,IDX,u,p)  XMOMENTUM1_1D(U,IDX)
#define INVISCIDFLUX2_XDIR1_1D(U,IDX,u,p) (XMOMENTUM1_1D(U,IDX)*u+p)
#define INVISCIDFLUX3_XDIR1_1D(U,IDX,u,p) (TOTALENERGY1_1D(U,IDX)*u+p*u)

#define INVISCIDFLUX1_XDIR2_1D(U,IDX,i,n1,n2,u,p)  XMOMENTUM2_1D(U,IDX,i,n1,n2)
#define INVISCIDFLUX2_XDIR2_1D(U,IDX,i,n1,n2,u,p) (XMOMENTUM2_1D(U,IDX,i,n1,n2)*u+p)
#define INVISCIDFLUX3_XDIR2_1D(U,IDX,i,n1,n2,u,p) (TOTALENERGY2_1D(U,IDX,i,n1,n2)*u+p*u)

#define INVISCIDFLUX1_XDIR3_1D(U,IDX,i,j,n1,n2,n3,u,p)  XMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3)
#define INVISCIDFLUX2_XDIR3_1D(U,IDX,i,j,n1,n2,n3,u,p) (XMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3)*u+p)
#define INVISCIDFLUX3_XDIR3_1D(U,IDX,i,j,n1,n2,n3,u,p) (TOTALENERGY3_1D(U,IDX,i,j,n1,n2,n3)*u+p*u)

#if 0
! Flux Jacobian matrix for inviscid hydrodynamics in 1D
#endif

#define INVISCIDFLUXJACOBIMATRIX11_1D(dscale,CX,u,E) DCONST(0.0)
#define INVISCIDFLUXJACOBIMATRIX12_1D(dscale,CX,u,E) dscale*CX
#define INVISCIDFLUXJACOBIMATRIX13_1D(dscale,CX,u,E) DCONST(0.0)
#define INVISCIDFLUXJACOBIMATRIX21_1D(dscale,CX,u,E) dscale*((DCONST(HYDRO_GAMMA))-DCONST(3.0))/DCONST(2.0)*u*u*CX
#define INVISCIDFLUXJACOBIMATRIX22_1D(dscale,CX,u,E) dscale*(DCONST(3.0)-(DCONST(HYDRO_GAMMA)))*u*CX
#define INVISCIDFLUXJACOBIMATRIX23_1D(dscale,CX,u,E) dscale*((DCONST(HYDRO_GAMMA))-DCONST(1.0))*CX
#define INVISCIDFLUXJACOBIMATRIX31_1D(dscale,CX,u,E) dscale*(((DCONST(HYDRO_GAMMA))-DCONST(1.0))*u*u-(DCONST(HYDRO_GAMMA))*E)*u*CX
#define INVISCIDFLUXJACOBIMATRIX32_1D(dscale,CX,u,E) dscale*((DCONST(HYDRO_GAMMA))*E-DCONST(3.0)*((DCONST(HYDRO_GAMMA))-DCONST(1.0))/DCONST(2.0)*u*u)*CX
#define INVISCIDFLUXJACOBIMATRIX33_1D(dscale,CX,u,E) dscale*(DCONST(HYDRO_GAMMA))*u*CX

#if 0
! Flux Jacobian matrix in x-direction for inviscid hydrodynamics in 1D
#endif

#define INVISCIDFLUXJACOBIMATRIX11_XDIR_1D(dscale,u,E) INVISCIDFLUXJACOBIMATRIX11_1D(dscale,DCONST(1.0),u,E)
#define INVISCIDFLUXJACOBIMATRIX12_XDIR_1D(dscale,u,E) INVISCIDFLUXJACOBIMATRIX12_1D(dscale,DCONST(1.0),u,E)
#define INVISCIDFLUXJACOBIMATRIX13_XDIR_1D(dscale,u,E) INVISCIDFLUXJACOBIMATRIX13_1D(dscale,DCONST(1.0),u,E)
#define INVISCIDFLUXJACOBIMATRIX21_XDIR_1D(dscale,u,E) INVISCIDFLUXJACOBIMATRIX21_1D(dscale,DCONST(1.0),u,E)
#define INVISCIDFLUXJACOBIMATRIX22_XDIR_1D(dscale,u,E) INVISCIDFLUXJACOBIMATRIX22_1D(dscale,DCONST(1.0),u,E)
#define INVISCIDFLUXJACOBIMATRIX23_XDIR_1D(dscale,u,E) INVISCIDFLUXJACOBIMATRIX23_1D(dscale,DCONST(1.0),u,E)
#define INVISCIDFLUXJACOBIMATRIX31_XDIR_1D(dscale,u,E) INVISCIDFLUXJACOBIMATRIX31_1D(dscale,DCONST(1.0),u,E)
#define INVISCIDFLUXJACOBIMATRIX32_XDIR_1D(dscale,u,E) INVISCIDFLUXJACOBIMATRIX32_1D(dscale,DCONST(1.0),u,E)
#define INVISCIDFLUXJACOBIMATRIX33_XDIR_1D(dscale,u,E) INVISCIDFLUXJACOBIMATRIX33_1D(dscale,DCONST(1.0),u,E)


#if 0
!##############################################################################
! Fluxes and matrices in 2D
!##############################################################################
#endif

#if 0
! Flux in x-direction for inviscid hydrodynamics in 2D
#endif

#define INVISCIDFLUX1_XDIR1_2D(U,IDX,u,p)  XMOMENTUM1_2D(U,IDX)
#define INVISCIDFLUX2_XDIR1_2D(U,IDX,u,p) (XMOMENTUM1_2D(U,IDX)*u+p)
#define INVISCIDFLUX3_XDIR1_2D(U,IDX,u,p) (YMOMENTUM1_2D(U,IDX)*u)
#define INVISCIDFLUX4_XDIR1_2D(U,IDX,u,p) (TOTALENERGY1_2D(U,IDX)*u+p*u)

#define INVISCIDFLUX1_XDIR2_2D(U,IDX,i,n1,n2,u,p)  XMOMENTUM2_2D(U,IDX,i,n1,n2)
#define INVISCIDFLUX2_XDIR2_2D(U,IDX,i,n1,n2,u,p) (XMOMENTUM2_2D(U,IDX,i,n1,n2)*u+p)
#define INVISCIDFLUX3_XDIR2_2D(U,IDX,i,n1,n2,u,p) (YMOMENTUM2_2D(U,IDX,i,n1,n2)*u)
#define INVISCIDFLUX4_XDIR2_2D(U,IDX,i,n1,n2,u,p) (TOTALENERGY2_2D(U,IDX,i,n1,n2)*u+p*u)

#define INVISCIDFLUX1_XDIR3_2D(U,IDX,i,j,n1,n2,n3,u,p)  XMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)
#define INVISCIDFLUX2_XDIR3_2D(U,IDX,i,j,n1,n2,n3,u,p) (XMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)*u+p)
#define INVISCIDFLUX3_XDIR3_2D(U,IDX,i,j,n1,n2,n3,u,p) (YMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)*u)
#define INVISCIDFLUX4_XDIR3_2D(U,IDX,i,j,n1,n2,n3,u,p) (TOTALENERGY3_2D(U,IDX,i,j,n1,n2,n3)*u+p*u)

#if 0
! Flux in y-direction for inviscid hydrodynamics in 2D
#endif

#define INVISCIDFLUX1_YDIR1_2D(U,IDX,v,p)  YMOMENTUM1_2D(U,IDX)
#define INVISCIDFLUX2_YDIR1_2D(U,IDX,v,p) (XMOMENTUM1_2D(U,IDX)*v)
#define INVISCIDFLUX3_YDIR1_2D(U,IDX,v,p) (YMOMENTUM1_2D(U,IDX)*v+p)
#define INVISCIDFLUX4_YDIR1_2D(U,IDX,v,p) (TOTALENERGY1_2D(U,IDX)*v+p*v)

#define INVISCIDFLUX1_YDIR2_2D(U,IDX,i,n1,n2,v,p)  YMOMENTUM2_2D(U,IDX,i,n1,n2)
#define INVISCIDFLUX2_YDIR2_2D(U,IDX,i,n1,n2,v,p) (XMOMENTUM2_2D(U,IDX,i,n1,n2)*v)
#define INVISCIDFLUX3_YDIR2_2D(U,IDX,i,n1,n2,v,p) (YMOMENTUM2_2D(U,IDX,i,n1,n2)*v+p)
#define INVISCIDFLUX4_YDIR2_2D(U,IDX,i,n1,n2,v,p) (TOTALENERGY2_2D(U,IDX,i,n1,n2)*v+p*v)

#define INVISCIDFLUX1_YDIR3_2D(U,IDX,i,j,n1,n2,n3,v,p)  YMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)
#define INVISCIDFLUX2_YDIR3_2D(U,IDX,i,j,n1,n2,n3,v,p) (XMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)*v)
#define INVISCIDFLUX3_YDIR3_2D(U,IDX,i,j,n1,n2,n3,v,p) (YMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)*v+p)
#define INVISCIDFLUX4_YDIR3_2D(U,IDX,i,j,n1,n2,n3,v,p) (TOTALENERGY3_2D(U,IDX,i,j,n1,n2,n3)*v+p*v)

#if 0
! Flux Jacobian matrix for inviscid hydrodynamics in 2D
#endif

#define INVISCIDFLUXJACOBIMATRIX11_2D(dscale,CX,CY,u,v,E) DCONST(0.0)
#define INVISCIDFLUXJACOBIMATRIX12_2D(dscale,CX,CY,u,v,E) dscale*CX
#define INVISCIDFLUXJACOBIMATRIX13_2D(dscale,CX,CY,u,v,E) dscale*CY
#define INVISCIDFLUXJACOBIMATRIX14_2D(dscale,CX,CY,u,v,E) DCONST(0.0)
#define INVISCIDFLUXJACOBIMATRIX21_2D(dscale,CX,CY,u,v,E) dscale*((((DCONST(HYDRO_GAMMA))-DCONST(1.0))/DCONST(2.0)*(u*u+v*v)-MYNEWLINE \
                                                                   u*u)*CX-u*v*CY)
#define INVISCIDFLUXJACOBIMATRIX22_2D(dscale,CX,CY,u,v,E) dscale*((DCONST(3.0)-(DCONST(HYDRO_GAMMA)))*u*CX+v*CY)
#define INVISCIDFLUXJACOBIMATRIX23_2D(dscale,CX,CY,u,v,E) dscale*(u*CY-((DCONST(HYDRO_GAMMA))-DCONST(1.0))*v*CX)
#define INVISCIDFLUXJACOBIMATRIX24_2D(dscale,CX,CY,u,v,E) dscale*((DCONST(HYDRO_GAMMA))-DCONST(1.0))*CX
#define INVISCIDFLUXJACOBIMATRIX31_2D(dscale,CX,CY,u,v,E) dscale*((((DCONST(HYDRO_GAMMA))-DCONST(1.0))/DCONST(2.0)*(u*u+v*v)-v*v)*CY-MYNEWLINE \
                                                                  u*v*CX)
#define INVISCIDFLUXJACOBIMATRIX32_2D(dscale,CX,CY,u,v,E) dscale*(v*CX-((DCONST(HYDRO_GAMMA))-DCONST(1.0))*u*CY)
#define INVISCIDFLUXJACOBIMATRIX33_2D(dscale,CX,CY,u,v,E) dscale*(u*CX+(DCONST(3.0)-(DCONST(HYDRO_GAMMA)))*v*CY)
#define INVISCIDFLUXJACOBIMATRIX34_2D(dscale,CX,CY,u,v,E) dscale*((DCONST(HYDRO_GAMMA))-DCONST(1.0))*CY
#define INVISCIDFLUXJACOBIMATRIX41_2D(dscale,CX,CY,u,v,E) dscale*(((DCONST(HYDRO_GAMMA))-DCONST(1.0))*(u*u+v*v)-(DCONST(HYDRO_GAMMA))*E)*(u*CX+v*CY)
#define INVISCIDFLUXJACOBIMATRIX42_2D(dscale,CX,CY,u,v,E) dscale*(((DCONST(HYDRO_GAMMA))*E-((DCONST(HYDRO_GAMMA))-DCONST(1.0))/DCONST(2.0)*(u*u+v*v))*CX-MYNEWLINE \
                                                                  ((DCONST(HYDRO_GAMMA))-DCONST(1.0))*u*(u*CX+v*CY))
#define INVISCIDFLUXJACOBIMATRIX43_2D(dscale,CX,CY,u,v,E) dscale*(((DCONST(HYDRO_GAMMA))*E-((DCONST(HYDRO_GAMMA))-DCONST(1.0))/DCONST(2.0)*(u*u+v*v))*CY-MYNEWLINE \
                                                                  ((DCONST(HYDRO_GAMMA))-DCONST(1.0))*v*(u*CX+v*CY))
#define INVISCIDFLUXJACOBIMATRIX44_2D(dscale,CX,CY,u,v,E) dscale*((DCONST(HYDRO_GAMMA))*(u*CX+v*CY))

#if 0
! Flux Jacobian matrix in x-direction for inviscid hydrodynamics in 2D
#endif

#define INVISCIDFLUXJACOBIMATRIX11_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX11_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX12_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX12_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX13_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX13_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX14_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX14_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX21_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX21_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX22_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX22_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX23_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX23_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX24_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX24_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX31_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX31_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX32_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX32_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX33_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX33_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX34_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX34_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX41_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX41_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX42_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX42_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX43_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX43_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX44_XDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX44_2D(dscale,DCONST(1.0),DCONST(0.0),u,v,E)

#if 0
! Flux Jacobian matrix in y-direction for inviscid hydrodynamics in 2D
#endif

#define INVISCIDFLUXJACOBIMATRIX11_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX11_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX12_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX12_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX13_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX13_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX14_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX14_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX21_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX21_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX22_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX22_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX23_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX23_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX24_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX24_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX31_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX31_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX32_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX32_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX33_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX33_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX34_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX34_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX41_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX41_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX42_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX42_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX43_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX43_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)
#define INVISCIDFLUXJACOBIMATRIX44_YDIR_2D(dscale,u,v,E) INVISCIDFLUXJACOBIMATRIX44_2D(dscale,DCONST(0.0),DCONST(1.0),u,v,E)


#if 0
!##############################################################################
! Fluxes and matrices in 3D
!##############################################################################
#endif

#if 0
! Flux in x-direction for inviscid hydrodynamics in 3D
#endif

#define INVISCIDFLUX1_XDIR1_3D(U,IDX,u,p)  XMOMENTUM1_3D(U,IDX)
#define INVISCIDFLUX2_XDIR1_3D(U,IDX,u,p) (XMOMENTUM1_3D(U,IDX)*u+p)
#define INVISCIDFLUX3_XDIR1_3D(U,IDX,u,p) (YMOMENTUM1_3D(U,IDX)*u)
#define INVISCIDFLUX4_XDIR1_3D(U,IDX,u,p) (ZMOMENTUM1_3D(U,IDX)*u)
#define INVISCIDFLUX5_XDIR1_3D(U,IDX,u,p) (TOTALENERGY1_3D(U,IDX)*u+p*u)

#define INVISCIDFLUX1_XDIR2_3D(U,IDX,i,n1,n2,u,p)  XMOMENTUM2_3D(U,IDX,i,n1,n2)
#define INVISCIDFLUX2_XDIR2_3D(U,IDX,i,n1,n2,u,p) (XMOMENTUM2_3D(U,IDX,i,n1,n2)*u+p)
#define INVISCIDFLUX3_XDIR2_3D(U,IDX,i,n1,n2,u,p) (YMOMENTUM2_3D(U,IDX,i,n1,n2)*u)
#define INVISCIDFLUX4_XDIR2_3D(U,IDX,i,n1,n2,u,p) (ZMOMENTUM2_3D(U,IDX,i,n1,n2)*u)
#define INVISCIDFLUX5_XDIR2_3D(U,IDX,i,n1,n2,u,p) (TOTALENERGY2_3D(U,IDX,i,n1,n2)*u+p*u)

#define INVISCIDFLUX1_XDIR3_3D(U,IDX,i,j,n1,n2,n3,u,p)  XMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)
#define INVISCIDFLUX2_XDIR3_3D(U,IDX,i,j,n1,n2,n3,u,p) (XMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)*u+p)
#define INVISCIDFLUX3_XDIR3_3D(U,IDX,i,j,n1,n2,n3,u,p) (YMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)*u)
#define INVISCIDFLUX4_XDIR3_3D(U,IDX,i,j,n1,n2,n3,u,p) (ZMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)*u)
#define INVISCIDFLUX5_XDIR3_3D(U,IDX,i,j,n1,n2,n3,u,p) (TOTALENERGY3_3D(U,IDX,i,j,n1,n2,n3)*u+p*u)

#if 0
! Flux in y-direction for inviscid hydrodynamics in 3D
#endif

#define INVISCIDFLUX1_YDIR1_3D(U,IDX,v,p)  YMOMENTUM1_3D(U,IDX)
#define INVISCIDFLUX2_YDIR1_3D(U,IDX,v,p) (XMOMENTUM1_3D(U,IDX)*v)
#define INVISCIDFLUX3_YDIR1_3D(U,IDX,v,p) (YMOMENTUM1_3D(U,IDX)*v+p)
#define INVISCIDFLUX4_YDIR1_3D(U,IDX,v,p) (ZMOMENTUM1_3D(U,IDX)*v)
#define INVISCIDFLUX5_YDIR1_3D(U,IDX,v,p) (TOTALENERGY1_3D(U,IDX)*v+p*v)

#define INVISCIDFLUX1_YDIR2_3D(U,IDX,i,n1,n2,v,p)  YMOMENTUM2_3D(U,IDX,i,n1,n2)
#define INVISCIDFLUX2_YDIR2_3D(U,IDX,i,n1,n2,v,p) (XMOMENTUM2_3D(U,IDX,i,n1,n2)*v)
#define INVISCIDFLUX3_YDIR2_3D(U,IDX,i,n1,n2,v,p) (YMOMENTUM2_3D(U,IDX,i,n1,n2)*v+p)
#define INVISCIDFLUX4_YDIR2_3D(U,IDX,i,n1,n2,v,p) (ZMOMENTUM2_3D(U,IDX,i,n1,n2)*v)
#define INVISCIDFLUX5_YDIR2_3D(U,IDX,i,n1,n2,v,p) (TOTALENERGY2_3D(U,IDX,i,n1,n2)*v+p*v)

#define INVISCIDFLUX1_YDIR3_3D(U,IDX,i,j,n1,n2,n3,v,p)  YMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)
#define INVISCIDFLUX2_YDIR3_3D(U,IDX,i,j,n1,n2,n3,v,p) (XMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)*v)
#define INVISCIDFLUX3_YDIR3_3D(U,IDX,i,j,n1,n2,n3,v,p) (YMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)*v+p)
#define INVISCIDFLUX4_YDIR3_3D(U,IDX,i,j,n1,n2,n3,v,p) (ZMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)*v)
#define INVISCIDFLUX5_YDIR3_3D(U,IDX,i,j,n1,n2,n3,v,p) (TOTALENERGY3_3D(U,IDX,i,j,n1,n2,n3)*v+p*v)

#if 0
! Flux in z-direction for inviscid hydrodynamics in 3D
#endif

#define INVISCIDFLUX1_ZDIR1_3D(U,IDX,w,p)  ZMOMENTUM1_3D(U,IDX)
#define INVISCIDFLUX2_ZDIR1_3D(U,IDX,w,p) (XMOMENTUM1_3D(U,IDX)*w)
#define INVISCIDFLUX3_ZDIR1_3D(U,IDX,w,p) (YMOMENTUM1_3D(U,IDX)*w)
#define INVISCIDFLUX4_ZDIR1_3D(U,IDX,w,p) (ZMOMENTUM1_3D(U,IDX)*w+p)
#define INVISCIDFLUX5_ZDIR1_3D(U,IDX,w,p) (TOTALENERGY1_3D(U,IDX)*w+p*w)

#define INVISCIDFLUX1_ZDIR2_3D(U,IDX,i,n1,n2,w,p)  ZMOMENTUM2_3D(U,IDX,i,n1,n2)
#define INVISCIDFLUX2_ZDIR2_3D(U,IDX,i,n1,n2,w,p) (XMOMENTUM2_3D(U,IDX,i,n1,n2)*w)
#define INVISCIDFLUX3_ZDIR2_3D(U,IDX,i,n1,n2,w,p) (YMOMENTUM2_3D(U,IDX,i,n1,n2)*w)
#define INVISCIDFLUX4_ZDIR2_3D(U,IDX,i,n1,n2,w,p) (ZMOMENTUM2_3D(U,IDX,i,n1,n2)*w+p)
#define INVISCIDFLUX5_ZDIR2_3D(U,IDX,i,n1,n2,w,p) (TOTALENERGY2_3D(U,IDX,i,n1,n2)*w+p*w)

#define INVISCIDFLUX1_ZDIR3_3D(U,IDX,i,j,n1,n2,n3,w,p)  ZMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)
#define INVISCIDFLUX2_ZDIR3_3D(U,IDX,i,j,n1,n2,n3,w,p) (XMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)*w)
#define INVISCIDFLUX3_ZDIR3_3D(U,IDX,i,j,n1,n2,n3,w,p) (YMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)*w)
#define INVISCIDFLUX4_ZDIR3_3D(U,IDX,i,j,n1,n2,n3,w,p) (ZMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)*w+p)
#define INVISCIDFLUX5_ZDIR3_3D(U,IDX,i,j,n1,n2,n3,w,p) (TOTALENERGY3_3D(U,IDX,i,j,n1,n2,n3)*w+p*w)

#if 0
! Flux Jacobian matrix for inviscid hydrodynamics in 3D
#endif

#define INVISCIDFLUXJACOBIMATRIX11_3D(dscale,CX,CY,CZ,u,v,w,E) DCONST(0.0)
#define INVISCIDFLUXJACOBIMATRIX12_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*CX
#define INVISCIDFLUXJACOBIMATRIX13_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*CY
#define INVISCIDFLUXJACOBIMATRIX14_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*CZ
#define INVISCIDFLUXJACOBIMATRIX15_3D(dscale,CX,CY,CZ,u,v,w,E) DCONST(0.0)
#define INVISCIDFLUXJACOBIMATRIX21_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*((((DCONST(HYDRO_GAMMA))-DCONST(1.0))/DCONST(2.0)*(ui*ui+vi*vi+wi*wi)-ui*ui)*CX-MYNEWLINE \
                                                                       ui*vi*CY-MYNEWLINE \
                                                                       ui*wi*CZ)
#define INVISCIDFLUXJACOBIMATRIX22_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*((DCONST(3.0)-(DCONST(HYDRO_GAMMA)))*ui*CX+MYNEWLINE \
                                                                       vi*CY+MYNEWLINE \
                                                                       wi*CZ)
#define INVISCIDFLUXJACOBIMATRIX23_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*(ui*CY-MYNEWLINE \
                                                                       ((DCONST(HYDRO_GAMMA))-DCONST(1.0))*vi*CX)
#define INVISCIDFLUXJACOBIMATRIX24_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*(ui*CZ-MYNEWLINE \
                                                                       ((DCONST(HYDRO_GAMMA))-DCONST(1.0))*wi*CX)
#define INVISCIDFLUXJACOBIMATRIX25_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*((DCONST(HYDRO_GAMMA))-DCONST(1.0))*CX
#define INVISCIDFLUXJACOBIMATRIX31_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*((((DCONST(HYDRO_GAMMA))-DCONST(1.0))/DCONST(2.0)*(ui*ui+vi*vi+wi*wi)-vi*vi)*CY-MYNEWLINE \
                                                                       ui*vi*CX-MYNEWLINE \
                                                                       vi*wi*CZ)
#define INVISCIDFLUXJACOBIMATRIX32_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*(vi*CX-MYNEWLINE \
                                                                       ((DCONST(HYDRO_GAMMA))-DCONST(1.0))*ui*CY)
#define INVISCIDFLUXJACOBIMATRIX33_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*((DCONST(3.0)-(DCONST(HYDRO_GAMMA)))*vi*CY+MYNEWLINE \
                                                                       ui*CX+MYNEWLINE \
                                                                       wi*CZ)
#define INVISCIDFLUXJACOBIMATRIX34_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*(vi*CZ-MYNEWLINE \
                                                                       ((DCONST(HYDRO_GAMMA))-DCONST(1.0))*wi*CY)
#define INVISCIDFLUXJACOBIMATRIX35_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*((DCONST(HYDRO_GAMMA))-DCONST(1.0))*CY
#define INVISCIDFLUXJACOBIMATRIX41_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*((((DCONST(HYDRO_GAMMA))-DCONST(1.0))/DCONST(2.0)*(ui*ui+vi*vi+wi*wi)-wi*wi)*CZ-MYNEWLINE \
                                                                       ui*wi*CX-MYNEWLINE \
                                                                       vi*wi*CY)
#define INVISCIDFLUXJACOBIMATRIX42_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*(wi*CX-MYNEWLINE \
                                                                       ((DCONST(HYDRO_GAMMA))-DCONST(1.0))*ui*CZ)
#define INVISCIDFLUXJACOBIMATRIX43_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*(wi*CY-MYNEWLINE \
                                                                       ((DCONST(HYDRO_GAMMA))-DCONST(1.0))*vi*CZ)
#define INVISCIDFLUXJACOBIMATRIX44_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*((DCONST(3.0)-(DCONST(HYDRO_GAMMA)))*wi*CZ+MYNEWLINE \
                                                                       ui*CX+MYNEWLINE \
                                                                       vi*CY)
#define INVISCIDFLUXJACOBIMATRIX45_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*((DCONST(HYDRO_GAMMA))-DCONST(1.0))*CZ
#define INVISCIDFLUXJACOBIMATRIX51_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*(((DCONST(HYDRO_GAMMA))-DCONST(1.0))*(ui*ui+vi*vi+wi*wi)-MYNEWLINE \
                                                                       (DCONST(HYDRO_GAMMA))*Ei)*(ui*CX+MYNEWLINE \
                                                                                                  vi*CY+MYNEWLINE \
                                                                                                  wi*CZ)
#define INVISCIDFLUXJACOBIMATRIX52_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*(((DCONST(HYDRO_GAMMA))*Ei-((DCONST(HYDRO_GAMMA))-DCONST(1.0))/DCONST(2.0)*(ui*ui+vi*vi+wi*wi))*CX-MYNEWLINE \
                                                                       ((DCONST(HYDRO_GAMMA))-DCONST(1.0))*ui*(ui*CX+MYNEWLINE \
                                                                                                               vi*CY+MYNEWLINE \
                                                                                                               wi*CZ))
#define INVISCIDFLUXJACOBIMATRIX53_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*(((DCONST(HYDRO_GAMMA))*Ei-((DCONST(HYDRO_GAMMA))-DCONST(1.0))/DCONST(2.0)*(ui*ui+vi*vi+wi*wi))*CY-MYNEWLINE \
                                                                       ((DCONST(HYDRO_GAMMA))-DCONST(1.0))*vi*(ui*CX+MYNEWLINE \
                                                                                                               vi*CY+MYNEWLINE \
                                                                                                               wi*CZ))
#define INVISCIDFLUXJACOBIMATRIX54_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*(((DCONST(HYDRO_GAMMA))*Ei-((DCONST(HYDRO_GAMMA))-DCONST(1.0))/DCONST(2.0)*(ui*ui+vi*vi+wi*wi))*CZ-MYNEWLINE \
                                                                       ((DCONST(HYDRO_GAMMA))-DCONST(1.0))*wi*(ui*CX+MYNEWLINE \
                                                                                                               vi*CY+MYNEWLINE \
                                                                                                               wi*CZ))
#define INVISCIDFLUXJACOBIMATRIX55_3D(dscale,CX,CY,CZ,u,v,w,E) dscale*((DCONST(HYDRO_GAMMA))*(ui*CX+MYNEWLINE \
                                                                                              vi*CY+MYNEWLINE \
                                                                                              wi*CZ))
      
#if 0
! Flux Jacobian matrix in x-direction for inviscid hydrodynamics in 2D
#endif

#define INVISCIDFLUXJACOBIMATRIX11_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX11_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX12_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX12_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX13_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX13_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX14_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX14_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX15_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX15_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX21_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX21_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX22_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX22_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX23_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX23_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX24_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX24_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX25_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX25_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX31_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX31_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX32_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX32_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX33_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX33_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX34_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX34_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX35_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX35_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX41_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX41_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX42_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX42_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX43_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX43_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX44_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX44_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX45_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX45_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX51_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX51_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX52_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX52_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX53_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX53_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX54_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX54_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX55_XDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX55_3D(dscale,DCONST(1.0),DCONST(0.0),DCONST(0.0),u,v,w,E)

#if 0
! Flux Jacobian matrix in y-direction for inviscid hydrodynamics in 2D
#endif

#define INVISCIDFLUXJACOBIMATRIX11_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX11_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX12_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX12_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX13_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX13_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX14_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX14_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX15_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX15_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX21_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX21_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX22_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX22_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX23_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX23_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX24_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX24_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX25_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX25_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX31_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX31_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX32_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX32_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX33_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX33_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX34_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX34_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX35_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX35_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX41_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX41_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX42_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX42_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX43_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX43_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX44_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX44_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX45_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX45_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX51_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX51_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX52_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX52_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX53_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX53_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX54_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX54_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX55_YDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX55_3D(dscale,DCONST(0.0),DCONST(1.0),DCONST(0.0),u,v,w,E)

#if 0
! Flux Jacobian matrix in z-direction for inviscid hydrodynamics in 2D
#endif

#define INVISCIDFLUXJACOBIMATRIX11_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX11_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX12_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX12_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX13_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX13_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX14_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX14_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX15_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX15_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX21_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX21_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX22_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX22_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX23_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX23_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX24_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX24_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX25_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX25_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX31_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX31_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX32_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX32_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX33_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX33_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX34_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX34_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX35_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX35_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX41_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX41_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX42_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX42_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX43_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX43_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX44_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX44_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX45_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX45_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX51_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX51_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX52_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX52_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX53_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX53_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX54_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX54_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)
#define INVISCIDFLUXJACOBIMATRIX55_ZDIR_3D(dscale,u,v,w,E) INVISCIDFLUXJACOBIMATRIX55_3D(dscale,DCONST(0.0),DCONST(0.0),DCONST(1.0),u,v,w,E)

#endif
