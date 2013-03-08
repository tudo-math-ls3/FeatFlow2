#include "../feat2constants.h"

! ------------------------------------------------------------------------------
! Convention: all internal identifiers have two leading and two
! trailing underscores and are undefined at the end of this file
! ------------------------------------------------------------------------------

#ifndef MatDT
#error "Definition of MatDT is missing!"
#else
#if   MatDT == SINGLE_PREC
#define __MatName__ SP
#define __MatType__ SP
#define __MatOne__  1.0_SP
#define __MatZero__ 0.0_SP
#elif MatDT == DOUBLE_PREC
#define __MatName__ DP
#define __MatType__ DP
#define __MatOne__  1.0_DP
#define __MatZero__ 0.0_DP
#elif MatDT == QUAD_PREC
#define __MatName__ QP
#define __MatType__ QP
#define __MatOne__  1.0_DP
#define __MatZero__ 0.0_QP
#else
#error "Unsupported matrix datatype!"
#endif
#endif

#ifndef VecDT
#error "Definition of VecDT is missing!"
#else
#if   VecDT == SINGLE_PREC
#define __VecName__ SP
#define __VecType__ SP
#define __VecOne__  1.0_SP
#define __VecZero__ 0.0_SP
#elif VecDT == DOUBLE_PREC
#define __VecName__ DP
#define __VecType__ DP
#define __VecOne__  1.0_DP
#define __VecZero__ 0.0_DP
#elif VecDT == QUAD_PREC
#define __VecName__ QP
#define __VecType__ QP
#define __VecOne__  1.0_QP
#define __VecZero__ 0.0_QP
#else
#error "Unsupported vector datatype!"
#endif
#endif

#define template_(name,mat,vec) name##mat##vec
#define template(name,mat,vec) template_(name,mat,vec)

  !**************************************************************
  ! Format 7 and Format 9 multiplication
  ! matrix dataype    : real(__MatType__)
  ! vectors data type : real(__VecType__)

  subroutine template(lsyssc_LAX79,__MatName__,__VecName__)&
      (Kld,Kcol,Da,Dx,Dy,cx,cy,NEQ,NVAR,rperfconfig)

    real(__MatType__), dimension(:), intent(in) :: Da
    real(__VecType__), dimension(:), intent(in) :: Dx
    real(__VecType__), dimension(:), intent(inout) :: Dy
    integer, dimension(:), intent(in) :: Kld
    integer, dimension(:), intent(in) :: Kcol
    real(__MatType__), intent(in) :: cx
    real(__VecType__), intent(in) :: cy
    integer, intent(in) :: NEQ,NVAR

    type(t_perfconfig), intent(in) :: rperfconfig

    integer :: ia,icol,irow,ivar
    real(__MatType__), dimension(NVAR) :: Ddtmp
    real(__MatType__) :: dtmp

    ! --------------------------------------------------------------------------
    ! Explicit instantiation for NVAR = 1,...,8
    ! --------------------------------------------------------------------------

    select case(NVAR)

    case(1)
#ifdef __NVAR__
#undef __NVAR__
#endif
#include "lsyssc_scalarMatVec_lax79.h"

    case(2)
#define __NVAR__ 2
#include "lsyssc_scalarMatVec_lax79.h"
#undef __NVAR__

    case(3)
#define __NVAR__ 3
#include "lsyssc_scalarMatVec_lax79.h"
#undef __NVAR__

    case(4)
#define __NVAR__ 4
#include "lsyssc_scalarMatVec_lax79.h"
#undef __NVAR__

    case(5)
#define __NVAR__ 5
#include "lsyssc_scalarMatVec_lax79.h"
#undef __NVAR__

    case(6)
#define __NVAR__ 6
#include "lsyssc_scalarMatVec_lax79.h"
#undef __NVAR__

    case(7)
#define __NVAR__ 7
#include "lsyssc_scalarMatVec_lax79.h"
#undef __NVAR__

    case(8)
#define __NVAR__ 8
#include "lsyssc_scalarMatVec_lax79.h"
#undef __NVAR__

    case default
#define __NVAR__ NVAR
#include "lsyssc_scalarMatVec_lax79.h"
#undef __NVAR__

    end select

  end subroutine

  !**************************************************************
  ! Format 9 row-compressed multiplication
  ! matrix dataype    : real(__MatType__)
  ! vectors data type : real(__VecType__)

  subroutine template(lsyssc_LAX9rowc,__MatName__,__VecName__)&
      (Kld,Kcol,KrowIdx,Da,Dx,Dy,cx,cy,NEQ,NVAR,rperfconfig)

    real(__MatType__), dimension(:), intent(in) :: Da
    real(__VecType__), dimension(:), intent(in) :: Dx
    real(__VecType__), dimension(:), intent(inout) :: Dy
    integer, dimension(:), intent(in) :: Kld
    integer, dimension(:), intent(in) :: Kcol
    integer, dimension(:), intent(in) :: KrowIdx
    real(__MatType__), intent(in) :: cx
    real(__VecType__), intent(in) :: cy
    integer, intent(in) :: NEQ,NVAR

    type(t_perfconfig), intent(in) :: rperfconfig

    integer :: ia,icol,irow,ivar
    real(__MatType__), dimension(NVAR) :: Ddtmp
    real(__MatType__) :: dtmp

    ! --------------------------------------------------------------------------
    ! Explicit instantiation for NVAR = 1,...,8
    ! --------------------------------------------------------------------------

    select case(NVAR)

    case(1)
#ifdef __NVAR__
#undef __NVAR__
#endif
#include "lsyssc_scalarMatVec_lax9rowc.h"

    case(2)
#define __NVAR__ 2
#include "lsyssc_scalarMatVec_lax9rowc.h"
#undef __NVAR__

    case(3)
#define __NVAR__ 3
#include "lsyssc_scalarMatVec_lax9rowc.h"
#undef __NVAR__

    case(4)
#define __NVAR__ 4
#include "lsyssc_scalarMatVec_lax9rowc.h"
#undef __NVAR__

    case(5)
#define __NVAR__ 5
#include "lsyssc_scalarMatVec_lax9rowc.h"
#undef __NVAR__

    case(6)
#define __NVAR__ 6
#include "lsyssc_scalarMatVec_lax9rowc.h"
#undef __NVAR__

    case(7)
#define __NVAR__ 7
#include "lsyssc_scalarMatVec_lax9rowc.h"
#undef __NVAR__

    case(8)
#define __NVAR__ 8
#include "lsyssc_scalarMatVec_lax9rowc.h"
#undef __NVAR__

    case default
#define __NVAR__ NVAR
#include "lsyssc_scalarMatVec_lax9rowc.h"
#undef __NVAR__

    end select

  end subroutine

  !**************************************************************
  ! Format 7 and Format 9 full interleaved multiplication
  ! matrix dataype    : real(__MatType__)
  ! vectors data type : real(__VecType__)

  subroutine template(lsyssc_LAX79INTL1,__MatName__,__VecName__)&
      (Kld,Kcol,Da,Dx,Dy,cx,cy,NEQ,NVAR,rperfconfig)

    real(__MatType__), dimension(:), intent(in) :: Da
    real(__VecType__), dimension(:), intent(in) :: Dx
    real(__VecType__), dimension(:), intent(inout) :: Dy
    integer, dimension(:), intent(in) :: Kld
    integer, dimension(:), intent(in) :: Kcol
    real(__MatType__), intent(in) :: cx
    real(__VecType__), intent(in) :: cy
    integer, intent(in) :: NEQ,NVAR

    type(t_perfconfig), intent(in) :: rperfconfig

    integer :: ia,icol,irow,ivar,jvar
    real(__MatType__), dimension(NVAR) :: Ddtmp
    real(__MatType__) :: dtmp

    ! --------------------------------------------------------------------------
    ! Explicit instantiation for NVAR = 1,...,8
    ! --------------------------------------------------------------------------

    select case(NVAR)

    case(1)
#ifdef __NVAR__
#undef __NVAR__
#endif
#include "lsyssc_scalarMatVec_lax79.h"

    case(2)
#define __NVAR__ 2
#include "lsyssc_scalarMatVec_lax79intl1.h"
#undef __NVAR__

    case(3)
#define __NVAR__ 3
#include "lsyssc_scalarMatVec_lax79intl1.h"
#undef __NVAR__

    case(4)
#define __NVAR__ 4
#include "lsyssc_scalarMatVec_lax79intl1.h"
#undef __NVAR__

    case(5)
#define __NVAR__ 5
#include "lsyssc_scalarMatVec_lax79intl1.h"
#undef __NVAR__

    case(6)
#define __NVAR__ 6
#include "lsyssc_scalarMatVec_lax79intl1.h"
#undef __NVAR__

    case(7)
#define __NVAR__ 7
#include "lsyssc_scalarMatVec_lax79intl1.h"
#undef __NVAR__

    case(8)
#define __NVAR__ 8
#include "lsyssc_scalarMatVec_lax79intl1.h"
#undef __NVAR__

    case default
#define __NVAR__ NVAR
#include "lsyssc_scalarMatVec_lax79intl1.h"
#undef __NVAR__

    end select

  end subroutine

  !**************************************************************
  ! Format 7 and Format 9 diagonal interleaved multiplication
  ! matrix dataype    : real(__MatType__)
  ! vectors data type : real(__VecType__)

  subroutine template(lsyssc_LAX79INTLD,__MatName__,__VecName__)&
      (Kld,Kcol,Da,Dx,Dy,cx,cy,NEQ,NVAR,rperfconfig)

    real(__MatType__), dimension(:), intent(in) :: Da
    real(__VecType__), dimension(:), intent(in) :: Dx
    real(__VecType__), dimension(:), intent(inout) :: Dy
    integer, dimension(:), intent(in) :: Kld
    integer, dimension(:), intent(in) :: Kcol
    real(__MatType__), intent(in) :: cx
    real(__VecType__), intent(in) :: cy
    integer, intent(in) :: NEQ,NVAR

    type(t_perfconfig), intent(in) :: rperfconfig

    integer :: ia,icol,irow,ivar
    real(__MatType__), dimension(NVAR) :: Ddtmp
    real(__MatType__) :: dtmp

    ! --------------------------------------------------------------------------
    ! Explicit instantiation for NVAR = 1,...,8
    ! --------------------------------------------------------------------------

    select case(NVAR)

    case(1)
#ifdef __NVAR__
#undef __NVAR__
#endif
#include "lsyssc_scalarMatVec_lax79.h"

    case(2)
#define __NVAR__ 2
#include "lsyssc_scalarMatVec_lax79intld.h"
#undef __NVAR__

    case(3)
#define __NVAR__ 3
#include "lsyssc_scalarMatVec_lax79intld.h"
#undef __NVAR__

    case(4)
#define __NVAR__ 4
#include "lsyssc_scalarMatVec_lax79intld.h"
#undef __NVAR__

    case(5)
#define __NVAR__ 5
#include "lsyssc_scalarMatVec_lax79intld.h"
#undef __NVAR__

    case(6)
#define __NVAR__ 6
#include "lsyssc_scalarMatVec_lax79intld.h"
#undef __NVAR__

    case(7)
#define __NVAR__ 7
#include "lsyssc_scalarMatVec_lax79intld.h"
#undef __NVAR__

    case(8)
#define __NVAR__ 8
#include "lsyssc_scalarMatVec_lax79intld.h"
#undef __NVAR__

    case default
#define __NVAR__ NVAR
#include "lsyssc_scalarMatVec_lax79intld.h"
#undef __NVAR__

    end select

  end subroutine

  !**************************************************************
  ! Format D (diagonal matrix) multiplication
  ! matrix dataype    : real(__MatType__)
  ! vectors data type : real(__VecType__)

  subroutine template(lsyssc_LATXD,__MatName__,__VecName__)&
      (Da,Dx,Dy,cx,cy,NEQ,NVAR,rperfconfig)

    real(__MatType__), dimension(:), intent(in) :: Da
    real(__VecType__), dimension(:), intent(in) :: Dx
    real(__VecType__), dimension(:), intent(inout) :: Dy
    real(__MatType__), intent(in) :: cx
    real(__VecType__), intent(in) :: cy
    integer, intent(in) :: NEQ,NVAR

    type(t_perfconfig), intent(in) :: rperfconfig    

    integer :: irow,ivar

select case(NVAR)

    case(1)
#ifdef __NVAR__
#undef __NVAR__
#endif
#include "lsyssc_scalarMatVec_latxd.h"

    case(2)
#define __NVAR__ 2
#include "lsyssc_scalarMatVec_latxd.h"
#undef __NVAR__

    case(3)
#define __NVAR__ 3
#include "lsyssc_scalarMatVec_latxd.h"
#undef __NVAR__

    case(4)
#define __NVAR__ 4
#include "lsyssc_scalarMatVec_latxd.h"
#undef __NVAR__

    case(5)
#define __NVAR__ 5
#include "lsyssc_scalarMatVec_latxd.h"
#undef __NVAR__

    case(6)
#define __NVAR__ 6
#include "lsyssc_scalarMatVec_latxd.h"
#undef __NVAR__

    case(7)
#define __NVAR__ 7
#include "lsyssc_scalarMatVec_latxd.h"
#undef __NVAR__

    case(8)
#define __NVAR__ 8
#include "lsyssc_scalarMatVec_latxd.h"
#undef __NVAR__

    case default
#define __NVAR__ NVAR
#include "lsyssc_scalarMatVec_latxd.h"
#undef __NVAR__

    end select

  end subroutine

  !**************************************************************
  ! Format 7 and Format 9 multiplication, transposed matrix
  ! matrix dataype    : real(__MatType__)
  ! vectors data type : real(__VecType__)

  subroutine template(lsyssc_LTX79,__MatName__,__VecName__)&
      (Kld,Kcol,Da,Dx,Dy,cx,cy,NEQ,NVAR,rperfconfig)

    real(__MatType__), dimension(:), intent(in) :: Da
    real(__VecType__), dimension(:), intent(in) :: Dx
    real(__VecType__), dimension(:), intent(inout) :: Dy
    integer, dimension(:), intent(in) :: Kld
    integer, dimension(:), intent(in) :: Kcol
    real(__MatType__), intent(in) :: cx
    real(__VecType__), intent(in) :: cy
    integer, intent(in) :: NEQ,NVAR

    type(t_perfconfig), intent(in) :: rperfconfig

    integer :: ia,icol,irow,ivar
    real(__MatType__), dimension(NVAR) :: Ddtmp
    real(__MatType__) :: dtmp

    ! --------------------------------------------------------------------------
    ! Explicit instantiation for NVAR = 1,...,8
    ! --------------------------------------------------------------------------

    select case(NVAR)

    case(1)
#ifdef __NVAR__
#undef __NVAR__
#endif
#include "lsyssc_scalarMatVec_ltx79.h"

    case(2)
#define __NVAR__ 2
#include "lsyssc_scalarMatVec_ltx79.h"
#undef __NVAR__

    case(3)
#define __NVAR__ 3
#include "lsyssc_scalarMatVec_ltx79.h"
#undef __NVAR__

    case(4)
#define __NVAR__ 4
#include "lsyssc_scalarMatVec_ltx79.h"
#undef __NVAR__

    case(5)
#define __NVAR__ 5
#include "lsyssc_scalarMatVec_ltx79.h"
#undef __NVAR__

    case(6)
#define __NVAR__ 6
#include "lsyssc_scalarMatVec_ltx79.h"
#undef __NVAR__

    case(7)
#define __NVAR__ 7
#include "lsyssc_scalarMatVec_ltx79.h"
#undef __NVAR__

    case(8)
#define __NVAR__ 8
#include "lsyssc_scalarMatVec_ltx79.h"
#undef __NVAR__

    case default
#define __NVAR__ NVAR
#include "lsyssc_scalarMatVec_ltx79.h"
#undef __NVAR__

    end select

  end subroutine

! Undefine all internal identifiers
#undef __MatName__
#undef __MatType__
#undef __MatOne__
#undef __MatZero__

#undef __VecName__
#undef __VecType__
#undef __VecOne__
#undef __VecZero__

#undef template_
#undef template
