!##############################################################################
!# ****************************************************************************
!# <name> transport_callback3d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve scalar conservation laws in 3D.
!#
!# The following routines are available:
!#
!# 1.) transp_setVariable1d
!#     -> Sets global variables for external data, e.g., velocity fields in 3D
!#
!# 2.) transp_calcMatrixPrimalConst3d
!#     -> Calculates the transport coefficients for linear convection in 3D
!#
!# 3.) transp_calcMatrixDualConst3d
!#     -> Calculates the transport coefficients for linear convection in 3D
!#
!# 4.) transp_hadaptCallback3d
!#      -> Performs application specific tasks in the adaptation algorithm in 3D
!#!# </purpose>
!##############################################################################

module transport_callback3d

  use collection
  use flagship_callback
  use fsystem
  use genoutput
  use hadaptaux
  use linearsystemblock
  use storage

  implicit none

  private
  public :: transp_setVariable3d
  public :: transp_calcMatrixPrimalConst3d
  public :: transp_calcMatrixDualConst3d
  public :: transp_hadaptCallback3d

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

  subroutine transp_setVariable3d(rvector, ivariable)

!<description>
    ! This subroutine sets one of the the global pointers to the given vector.
!</description>

!<input>
    ! scalar vector
    type(t_vectorScalar), intent(IN) :: rvector

    ! variable number
    integer, intent(IN) :: ivariable
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
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_setVariable3d')
      call sys_halt()
    end select
    
  end subroutine transp_setVariable3d

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatrixPrimalConst3d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 1D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>
    
!<input>
    ! solution vector
    real(DP), intent(IN) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij, C_ji

    ! nodal indices
    integer, intent(IN) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(OUT) :: k_ij, k_ji, d_ij
!</output>
!</subroutine>

    ! Compute convective coefficients
    k_ij = -p_Dvariable1(j)*C_ij(1)-p_Dvariable2(j)*C_ij(2)-p_Dvariable3(j)*C_ij(3)
    k_ji = -p_Dvariable1(i)*C_ji(1)-p_Dvariable2(i)*C_ji(2)-p_Dvariable3(i)*C_ji(3)

    ! Compute artificial diffusion coefficient
    d_ij = max(-k_ij, 0.0_DP, -k_ji)

  end subroutine transp_calcMatrixPrimalConst3d

  !*****************************************************************************

!<subroutine>

  pure subroutine transp_calcMatrixDualConst3d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 1D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>
    
!<input>
    ! solution vector
    real(DP), intent(IN) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij, C_ji

    ! nodal indices
    integer, intent(IN) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(OUT) :: k_ij, k_ji, d_ij
!</output>
!</subroutine>

    ! Compute convective coefficients
    k_ij = p_Dvariable1(j)*C_ij(1)+p_Dvariable2(j)*C_ij(2)+p_Dvariable3(j)*C_ij(3)
    k_ji = p_Dvariable1(i)*C_ji(1)+p_Dvariable2(i)*C_ji(2)+p_Dvariable3(i)*C_ji(3)

    ! Compute artificial diffusion coefficient
    d_ij = max(-k_ij, 0.0_DP, -k_ji)
    
  end subroutine transp_calcMatrixDualConst3d

  !*****************************************************************************

!<subroutine>

  subroutine transp_hadaptCallback3d(rcollection, iOperation, Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 3D.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(IN) :: iOperation

    ! Array of vertices involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ivertices

    ! Array of elements involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ielements
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolution
    real(DP), dimension(:), pointer, save :: p_Dsolution


    ! What operation should be performed
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! This subroutine assumes that the name of the solution vector
      ! is stored in the second quick access string.

      ! Retrieve solution vector from colletion and set pointer
      rsolution => collct_getvalue_vec(rcollection,&
                                        trim(rcollection%SquickAccess(2)))
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation, Ivertices, Ielements)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)
      
      ! Call the general callback function
      call flagship_hadaptCallback1d(rcollection, iOperation, Ivertices, Ielements)
      
      
    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      p_Dsolution(Ivertices(1)) = 0.5_DP*(p_Dsolution(Ivertices(2))+&    
                                           p_Dsolution(Ivertices(3)))

      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation, Ivertices, Ielements)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      p_Dsolution(Ivertices(1)) = 0.25_DP*(p_Dsolution(Ivertices(2))+&
                                            p_Dsolution(Ivertices(3))+&
                                            p_Dsolution(Ivertices(4))+&
                                            p_Dsolution(Ivertices(5)))

      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation, Ivertices, Ielements)

    
    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (Ivertices(2) .ne. 0) then
        p_Dsolution(Ivertices(1)) = p_Dsolution(Ivertices(2))
      else
        p_Dsolution(Ivertices(1)) = 0.0_DP
      end if

      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation, Ivertices, Ielements)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback3d(rcollection, iOperation, Ivertices, Ielements)
    end select
    
  end subroutine transp_hadaptCallback3d

end module transport_callback3d
