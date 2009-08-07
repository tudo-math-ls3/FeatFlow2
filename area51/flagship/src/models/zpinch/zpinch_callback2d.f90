!##############################################################################
!# ****************************************************************************
!# <name> euler_callback2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the simplified MHD equations in 2D.
!#
!# The following callback functions are available:
!#
!# 1.) zpinch_setVariable2d = zpinch_setVariableScalar2d /
!#                            zpinch_setVariableBlock2d
!#     -> Sets global variables for external data, e.g., velocity fields in 2D
!#
!# 2.) zpinch_hadaptCallbackScalar2d
!#     -> Performs application specific tasks in the adaptation
!#        algorithm in 2D, whereby the vector is stored in interleave format
!#
!# 3.) zpinch_hadaptCallbackBlock2d
!#     -> Performs application specific tasks in the adaptation
!#        algorithm in 2D, whereby the vector is stored in block format
!#
!# 4.) zpinch_calcMatRusConvectionP2d
!#     -> Calculates the transport coefficients for linear convection in 2D
!#
!# 5.) zpinch_calcMatRusConvectionD2d
!#     -> Calculates the transport coefficients for linear convection in 2D
!#
!# </purpose>
!##############################################################################

module zpinch_callback2d

  use collection
  use euler_basic
  use flagship_callback
  use fsystem
  use genoutput
  use hadaptaux
  use linearsystemblock
  use linearsystemscalar
  use storage

  implicit none
  
  private
  public :: zpinch_setVariable2d
  public :: zpinch_hadaptCallbackScalar2d
  public :: zpinch_hadaptCallbackBlock2d
  public :: zpinch_calcMatRusConvectionP2d
  public :: zpinch_calcMatRusConvectionD2d

  interface zpinch_setVariable2d
    module procedure zpinch_setVariableScalar2d
    module procedure zpinch_setVariableBlock2d
  end interface

!<globals>

  !*****************************************************************
  ! Pointers to external data vectors.
  !
  ! Using global variables is not good programming style but it is the
  ! only way to allow for an efficient access to the velocity data
  ! from within the callback routines which are called repeatedly

  real(DP), dimension(:), pointer, save :: p_Dvariable1 => null()
  real(DP), dimension(:), pointer, save :: p_Dvariable2 => null()
  real(DP), dimension(:), pointer, save :: p_Dvariable3 => null()
  real(DP), dimension(:), pointer, save :: p_Dvariable4 => null()
  real(DP), dimension(:), pointer, save :: p_Dvariable5 => null()

!</globals>

contains

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_setVariableScalar2d(rvector, ivariable)

!<description>
    ! This subroutine sets one of the the global pointers to
    ! the given scalar vector.
!</description>

!<input>
    ! scalar vector
    type(t_vectorScalar), intent(in) :: rvector

    ! variable number
    integer, intent(in) :: ivariable
!</input>
!</subroutine>

    select case(ivariable)
    case (1)
      call lsyssc_getbase_double(rvector, p_Dvariable1)
    case (2)
      call lsyssc_getbase_double(rvector, p_Dvariable2)
    case (3)
      call lsyssc_getbase_double(rvector, p_Dvariable3)
    case (4)
      call lsyssc_getbase_double(rvector, p_Dvariable4)
    case (5)
      call lsyssc_getbase_double(rvector, p_Dvariable5)
    case DEFAULT
      call output_line('Invalid variable number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'zpinch_setVariableScalar2d')
      call sys_halt()
    end select
    
  end subroutine zpinch_setVariableScalar2d

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_setVariableBlock2d(rvector, ivariable)

!<description>
    ! This subroutine sets one of the the global pointers to
    ! the given block vector.
!</description>

!<input>
    ! block vector
    type(t_vectorBlock), intent(in) :: rvector

    ! variable number
    integer, intent(in) :: ivariable
!</input>
!</subroutine>

    select case(ivariable)
    case (1)
      call lsysbl_getbase_double(rvector, p_Dvariable1)
    case (2)
      call lsysbl_getbase_double(rvector, p_Dvariable2)
    case (3)
      call lsysbl_getbase_double(rvector, p_Dvariable3)
    case (4)
      call lsysbl_getbase_double(rvector, p_Dvariable4)
    case (5)
      call lsysbl_getbase_double(rvector, p_Dvariable5)
    case DEFAULT
      call output_line('Invalid variable number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'zpinch_setVariableBlock2d')
      call sys_halt()
    end select
    
  end subroutine zpinch_setVariableBlock2d

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_hadaptCallbackScalar2d(rcollection, iOperation,&
      Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D. The solution vector is assumed
    ! to be store in scalar interleave format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation

    ! Array of vertices involved in the adaptivity step
    integer, dimension(:), intent(in) :: Ivertices

    ! Array of elements involved in the adaptivity step
    integer, dimension(:), intent(in) :: Ielements
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolutionEuler, rsolutionTransport
    real(DP), dimension(:), pointer, save :: p_DsolutionEuler, p_DsolutionTransport
    integer :: ivar


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve solution vectors from colletion and set pointer
      rsolutionEuler     => rcollection%p_rvectorQuickAccess1
      rsolutionTransport => rcollection%p_rvectorQuickAccess2
     
      ! Check if solution is stored in interleave format
      if (rsolutionEuler%nblocks .ne. 1) then
        call output_line('Vector is not in interleave format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_hadaptCallbackScalar2d')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation,&
          Ivertices, Ielements)
      
      
    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vectors
      nullify(rsolutionEuler, p_DsolutionEuler)
      nullify(rsolutionTransport, p_DsolutionTransport)
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation,&
          Ivertices, Ielements)
      

    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector for the Euler model
      if (rsolutionEuler%NEQ .ne. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionEuler, NVAR2D&
            *Ivertices(1), .false., .true.)
        call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      end if

      ! Resize solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .ne. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector for the Euler model
      if (rsolutionEuler%NEQ .lt. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionEuler, NVAR2D&
            *Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      end if
      do ivar = 1, NVAR2D
        p_DsolutionEuler((Ivertices(1)-1)*NVAR2D+ivar) = &
            0.5_DP*(p_DsolutionEuler((Ivertices(2)-1)*NVAR2D+ivar)+&
                    p_DsolutionEuler((Ivertices(3)-1)*NVAR2D+ivar))
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(Ivertices(1)) =&
          0.5_DP*(p_DsolutionTransport(Ivertices(2))+&    
                  p_DsolutionTransport(Ivertices(3)))

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation,&
          Ivertices, Ielements)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector for the Euler model
      if (rsolutionEuler%NEQ .lt. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionEuler, NVAR2D&
            *Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      end if
      do ivar = 1, NVAR2D
        p_DsolutionEuler((Ivertices(1)-1)*NVAR2D+ivar) = &
            0.25_DP*(p_DsolutionEuler((Ivertices(2)-1)*NVAR2D+ivar)+&
                     p_DsolutionEuler((Ivertices(3)-1)*NVAR2D+ivar)+&
                     p_DsolutionEuler((Ivertices(4)-1)*NVAR2D+ivar)+&
                     p_DsolutionEuler((Ivertices(5)-1)*NVAR2D+ivar))
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(Ivertices(1)) =&
          0.25_DP*(p_DsolutionTransport(Ivertices(2))+&
                   p_DsolutionTransport(Ivertices(3))+&
                   p_DsolutionTransport(Ivertices(4))+&
                   p_DsolutionTransport(Ivertices(5)))    
      

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation,&
          Ivertices, Ielements)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution for the Euler model
      if (Ivertices(2) .ne. 0) then
        do ivar = 1, NVAR2D
          p_DsolutionEuler((Ivertices(1)-1)*NVAR2D+ivar) = &
              p_DsolutionEuler((Ivertices(2)-1)*NVAR2D+ivar)
        end do
      else
        do ivar = 1, NVAR2D
          p_DsolutionEuler((Ivertices(1)-1)*NVAR2D+ivar) = 0.0_DP
        end do
      end if

      ! Remove vertex from solution for the scalar transport model
      if (Ivertices(2) .ne. 0) then
        p_DsolutionTransport(Ivertices(1)) =&
            p_DsolutionTransport(Ivertices(2))
      else
        p_DsolutionTransport(Ivertices(1)) = 0.0_DP
      end if
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation,&
          Ivertices, Ielements)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation,&
          Ivertices, Ielements)

    end select
    
  end subroutine zpinch_hadaptCallbackScalar2d

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_hadaptCallbackBlock2d(rcollection, iOperation,&
      Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D. The solution vector is assumed
    ! to be store in block format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation

    ! Array of vertices involved in the adaptivity step
    integer, dimension(:), intent(in) :: Ivertices

    ! Array of elements involved in the adaptivity step
    integer, dimension(:), intent(in) :: Ielements
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolutionEuler, rsolutionTransport
    real(DP), dimension(:), pointer, save :: p_DsolutionEuler, p_DsolutionTransport
    integer :: ivar,neq


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve solution vectors from colletion and set pointer
      rsolutionEuler     => rcollection%p_rvectorQuickAccess1
      rsolutionTransport => rcollection%p_rvectorQuickAccess2
      
      ! Check if solution is stored in interleave format
      if (rsolutionEuler%nblocks .ne. NVAR2D) then
        call output_line('Vector is not in block format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_hadaptCallbackBlock2d')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation,&
          Ivertices, Ielements)
      
      
    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vectors
      nullify(rsolutionEuler, p_DsolutionEuler)
      nullify(rsolutionTransport, p_DsolutionTransport)
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation,&
          Ivertices, Ielements)
      

    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector for the Euler model
      if (rsolutionEuler%NEQ .ne. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionEuler,&
            NVAR2D*Ivertices(1), .false., .true.)
        call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      end if

      ! Resize solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .ne. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector for the Euler model
      if (rsolutionEuler%NEQ .lt. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionEuler, NVAR2D&
            *Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      end if
      neq = rsolutionEuler%NEQ/NVAR2D
      do ivar = 1, NVAR2D
        p_DsolutionEuler((ivar-1)*neq+Ivertices(1)) = &
            0.5_DP*(p_DsolutionEuler((ivar-1)*neq+Ivertices(2))+&
                    p_DsolutionEuler((ivar-1)*neq+Ivertices(3)) )
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(Ivertices(1)) =&
          0.5_DP*(p_DsolutionTransport(Ivertices(2))+&    
                  p_DsolutionTransport(Ivertices(3)))

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation,&
          Ivertices, Ielements)

      
    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector for the Euler model
      if (rsolutionEuler%NEQ .lt. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionEuler, NVAR2D&
            *Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      end if
      neq = rsolutionEuler%NEQ/NVAR2D
      do ivar = 1, NVAR2D
        p_DsolutionEuler((ivar-1)*neq+Ivertices(1)) =&
            0.25_DP*(p_DsolutionEuler((ivar-1)*neq+Ivertices(2))+&
                     p_DsolutionEuler((ivar-1)*neq+Ivertices(3))+&
                     p_DsolutionEuler((ivar-1)*neq+Ivertices(4))+&
                     p_DsolutionEuler((ivar-1)*neq+Ivertices(5)) )
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(Ivertices(1)) =&
          0.25_DP*(p_DsolutionTransport(Ivertices(2))+&
                   p_DsolutionTransport(Ivertices(3))+&
                   p_DsolutionTransport(Ivertices(4))+&
                   p_DsolutionTransport(Ivertices(5)))

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation,&
          Ivertices, Ielements)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution for the Euler model
      if (Ivertices(2) .ne. 0) then
        neq = rsolutionEuler%NEQ/NVAR2D
        do ivar = 1, NVAR2D
          p_DsolutionEuler((ivar-1)*neq+Ivertices(1)) = &
              p_DsolutionEuler((ivar-1)*neq+Ivertices(2))
        end do
      else
        neq = rsolutionEuler%NEQ/NVAR2D
        do ivar = 1, NVAR2D
          p_DsolutionEuler((ivar-1)*neq+Ivertices(1)) = 0.0_DP
        end do
      end if

      ! Remove vertex from solution for the scalar transport model
      if (Ivertices(2) .ne. 0) then
        p_DsolutionTransport(Ivertices(1)) =&
            p_DsolutionTransport(Ivertices(2))
      else
        p_DsolutionTransport(Ivertices(1)) = 0.0_DP
      end if
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation,&
          Ivertices, Ielements)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation,&
          Ivertices, Ielements)

    end select
    
  end subroutine zpinch_hadaptCallbackBlock2d

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcMatRusConvectionP2d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 2D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>
    
!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij,k_ji,d_ij
!</output>
!</subroutine>

    ! local parameters
    real(DP), parameter :: GAMMA = 1.4_DP

    ! local variables
    real(DP) :: hi,hj,Ei,Ej,ui,uj,vi,vj,ci,cj, d2_ij
    integer :: idx, jdx

    ! Compute convective coefficients
    k_ij = -p_Dvariable1(j)*C_ij(1)-p_Dvariable2(j)*C_ij(2)
    k_ji = -p_Dvariable1(i)*C_ji(1)-p_Dvariable2(i)*C_ji(2)

!!$    ! Compute artificial diffusion coefficient
!!$    d_ij = max( abs(k_ij), abs(k_ji) )
    
    ! Compute base indices
    idx = 4*(i-1); jdx = 4*(j-1)
    
    ! Compute auxiliary variables
    ui = p_Dvariable3(idx+2)/p_Dvariable3(idx+1)
    vi = p_Dvariable3(idx+3)/p_Dvariable3(idx+1)
    Ei = p_Dvariable3(idx+4)/p_Dvariable3(idx+1)
    uj = p_Dvariable3(jdx+2)/p_Dvariable3(jdx+1)
    vj = p_Dvariable3(jdx+3)/p_Dvariable3(jdx+1)
    Ej = p_Dvariable3(jdx+4)/p_Dvariable3(jdx+1)

    ! Compute auxiliary quantities
    hi = GAMMA*Ei+(1-GAMMA)*0.5*(ui*ui+vi*vi)
    hj = GAMMA*Ej+(1-GAMMA)*0.5*(uj*uj+vj*vj)

    ci = sqrt(max((GAMMA-1)*(hi-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
    cj = sqrt(max((GAMMA-1)*(hj-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))

    ! Compute dissipation tensor D_ij
    d_ij = max( abs(C_ij(1)*uj+C_ij(2)*vj) + sqrt(C_ij(1)**2+C_ij(2)**2)*cj,&
                abs(C_ji(1)*ui+C_ji(2)*vi) + sqrt(C_ji(1)**2+C_ji(2)**2)*ci )
    
  end subroutine zpinch_calcMatRusConvectionP2d

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_calcMatRusConvectionD2d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 2D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>
    
!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij,k_ji,d_ij
!</output>
!</subroutine>

    ! Compute convective coefficients
    k_ij = p_Dvariable1(j)*C_ij(1)+p_Dvariable2(j)*C_ij(2)
    k_ji = p_Dvariable1(i)*C_ji(1)+p_Dvariable2(i)*C_ji(2)

    ! Compute artificial diffusion coefficient
    d_ij = max( abs(k_ij), abs(k_ji) )
    
  end subroutine zpinch_calcMatRusConvectionD2d
  
end module zpinch_callback2d
