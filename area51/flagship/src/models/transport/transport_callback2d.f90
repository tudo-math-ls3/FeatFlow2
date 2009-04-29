!##############################################################################
!# ****************************************************************************
!# <name> transport_callback2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve scalar conservation laws in 2D.
!#
!# The following routines are available:
!#
!# 1.) transp_setVariable1d
!#     -> Sets global variables for external data, e.g., velocity fields in 2D
!#
!# 2.) transp_calcMatrixPrimalConst2d
!#     -> Calculates the transport coefficients for linear convection in 2D
!#
!# 3.) transp_calcMatrixDualConst2d
!#     -> Calculates the transport coefficients for linear convection in 2D
!#
!# 4.) transp_calcMatrixPrimalBurgersSpT2d
!#     -> Calculates the transport coefficients for Burgers' equation in space-time
!#
!# 5.) transp_calcMatrixPrimalBuckLevSpT2d
!#     -> Calculates the transport coefficients for Buckley-Leverett equation in space-time
!#
!# 6.) transp_calcMatrixPrimalBurgers2d
!#     -> Calculates the transport coefficients for Burgers' equation in 2D
!#
!# 7.) transp_hadaptCallback2d
!#     -> Performs application specific tasks in the adaptation algorithm in 2D
!#
!# </purpose>
!##############################################################################

module transport_callback2d

  use collection
  use flagship_callback
  use fsystem
  use genoutput
  use hadaptaux
  use linearsystemblock
  use storage

  implicit none

  private
  public :: transp_setVariable2d
  public :: transp_calcMatrixPrimalConst2d
  public :: transp_calcMatrixDualConst2d
  public :: transp_calcMatrixPrimalBurgersSpT2d
  public :: transp_calcMatrixPrimalBuckLevSpT2d
  public :: transp_calcMatrixPrimalBurgers2d
  public :: transp_hadaptCallback2d

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

  subroutine transp_setVariable2d(rvector, ivariable)

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
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_setVariable2d')
      call sys_halt()
    end select
    
  end subroutine transp_setVariable2d

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatrixPrimalConst2d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 2D.
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
    real(DP), intent(OUT) :: k_ij,k_ji,d_ij
!</output>
!</subroutine>

    ! local variables
    real(DP) :: hi,hj,Ei,Ej,ui,uj,vi,vj,ci,cj
    integer :: idx,jdx

    ! Compute convective coefficients
    k_ij = -p_Dvariable1(j)*C_ij(1)-p_Dvariable2(j)*C_ij(2)
    k_ji = -p_Dvariable1(i)*C_ji(1)-p_Dvariable2(i)*C_ji(2)

!!$    ! Compute artificial diffusion coefficient
!!$    d_ij = max(-k_ij, 0.0_DP, -k_ji)

    ! Compute index positions
    idx = 4*(i-1)
    jdx = 4*(j-1)

    ! Compute velocities and energy
    ui = p_Dvariable3(idx+2)/p_Dvariable3(idx+1)
    vi = p_Dvariable3(idx+3)/p_Dvariable3(idx+1)
    Ei = p_Dvariable3(idx+4)/p_Dvariable3(idx+1)

    uj = p_Dvariable3(jdx+2)/p_Dvariable3(jdx+1)
    vj = p_Dvariable3(jdx+3)/p_Dvariable3(jdx+1)
    Ej = p_Dvariable3(jdx+4)/p_Dvariable3(jdx+1)

    ! Compute enthalpy
    hi = 1.4*Ei + (1-1.4)*0.5*(ui*ui+vi*vi)
    hj = 1.4*Ej + (1-1.4)*0.5*(uj*uj+vj*vj)

    ! Compute speed of sound
    ci = sqrt(max((1.4-1)*(hi-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
    cj = sqrt(max((1.4-1)*(hj-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))

    d_ij = max( abs(C_ij(1)*uj+C_ij(2)*vj) + sqrt(C_ij(1)*C_ij(1)+C_ij(2)*C_ij(2))*cj,&
                abs(C_ji(1)*ui+C_ji(2)*vi) + sqrt(C_ji(1)*C_ji(1)+C_ji(2)*C_ji(2))*ci )

  end subroutine transp_calcMatrixPrimalConst2d

  !*****************************************************************************

!<subroutine>

  pure subroutine transp_calcMatrixDualConst2d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 2D.
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
    real(DP), intent(OUT) :: k_ij,k_ji,d_ij
!</output>
!</subroutine>

    ! Compute convective coefficients
    k_ij = p_Dvariable1(j)*C_ij(1)+p_Dvariable2(j)*C_ij(2)
    k_ji = p_Dvariable1(i)*C_ji(1)+p_Dvariable2(i)*C_ji(2)

    ! Compute artificial diffusion coefficient
    d_ij = max(-k_ij, 0.0_DP, -k_ji)
    
  end subroutine transp_calcMatrixDualConst2d

  !*****************************************************************************
    
!<subroutine>

  pure subroutine transp_calcMatrixPrimalBurgersSpT2d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the 
    ! one-dimensional Burgers equation $du/dt+df(u)/dx=0$, whereby
    ! the flux function is given by $f(u)=0.5*u^2$.
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
    k_ij = -0.5_DP*(u_i+u_j)*C_ij(1)-C_ij(2)
    k_ji = -0.5_DP*(u_i+u_j)*C_ji(1)-C_ji(2)

    ! Compute artificial diffusion coefficient
    d_ij = max(-k_ij, 0.0_DP, -k_ji)

  end subroutine transp_calcMatrixPrimalBurgersSpT2d

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatrixPrimalBuckLevSpT2d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the 
    ! Buckley-Leverett equation $du/dt+df(u)/dx=0$, whereby the
    ! flux function is given by $f(u)=u^2/(u^2+0.5*(1-u)^2)$
    !
    ! Here, the characteristic velocity $a(u)=f^\prime(u)$ is given
    ! by $a(u)=\frac{4u(1-u)}{(3u^2-2u+1)^2}$.
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

    ! local variables
    real(DP) :: v_i,v_j
    
    ! Compute velocities
    v_i = 4*u_i*(1-u_i)/(3*u_i*u_i-2*u_i+1)**2
    v_j = 4*u_j*(1-u_j)/(3*u_j*u_j-2*u_j+1)**2

    ! Compute convective coefficients
    k_ij = -v_j*C_ij(1)-C_ij(2)
    k_ji = -v_i*C_ji(1)-C_ji(2)

    ! Compute artificial diffusion coefficient
    d_ij = max(-k_ij, 0.0_DP, -k_ji)
        
  end subroutine transp_calcMatrixPrimalBuckLevSpT2d

  !*****************************************************************************
    
!<subroutine>

  pure subroutine transp_calcMatrixPrimalBurgers2d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for Burgers' equation in 2D.
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
    k_ij = -0.5_DP*(u_i+u_j)*(C_ij(1)+C_ij(2))
    k_ji = -0.5_DP*(u_i+u_j)*(C_ji(1)+C_ji(2))

    ! Compute artificial diffusion coefficient
    d_ij = max(-k_ij, 0.0_DP, -k_ji)
    
  end subroutine transp_calcMatrixPrimalBurgers2d
  
  !*****************************************************************************

!<subroutine>

  pure subroutine transp_calcMatrixZPinchInterleave2d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for the externally prescribed velocity field resulting from the a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 2D.
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
    real(DP), intent(OUT) :: k_ij,k_ji,d_ij
!</output>
!</subroutine>

    ! local variables
    real(DP) :: hi,hj,Ei,Ej,ui,uj,vi,vj,ci,cj
    integer :: idx,jdx

    ! Compute convective coefficients
    k_ij = -p_Dvariable1(j)*C_ij(1)-p_Dvariable2(j)*C_ij(2)
    k_ji = -p_Dvariable1(i)*C_ji(1)-p_Dvariable2(i)*C_ji(2)

    ! TEST CODE >>>
    
    ! Compute index positions
    idx = 4*(i-1)
    jdx = 4*(j-1)

    ! Compute velocities and energy
    ui = p_Dvariable1(idx+2)/p_Dvariable1(idx+1)
    vi = p_Dvariable1(idx+3)/p_Dvariable1(idx+1)
    Ei = p_Dvariable1(idx+4)/p_Dvariable1(idx+1)

    uj = p_Dvariable1(jdx+2)/p_Dvariable1(jdx+1)
    vj = p_Dvariable1(jdx+3)/p_Dvariable1(jdx+1)
    Ej = p_Dvariable1(jdx+4)/p_Dvariable1(jdx+1)

    ! Compute enthalpy
    hi = 1.4*Ei + (1-1.4)*0.5*(ui*ui+vi*vi)
    hj = 1.4*Ej + (1-1.4)*0.5*(uj*uj+vj*vj)

    ! Compute speed of sound
    ci = sqrt(max((1.4-1)*(hi-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
    cj = sqrt(max((1.4-1)*(hj-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))

    d_ij = max( abs(C_ij(1)*uj+C_ij(2)*vj) + sqrt(C_ij(1)*C_ij(1)+C_ij(2)*C_ij(2))*cj,&
                abs(C_ji(1)*ui+C_ji(2)*vi) + sqrt(C_ji(1)*C_ji(1)+C_ji(2)*C_ji(2))*ci )

    ! TEST CODE <<<

  end subroutine transp_calcMatrixZPinchInterleave2d
  
  !*****************************************************************************

!<subroutine>

  subroutine transp_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D.
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
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)

      
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
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)


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
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)

      
    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (Ivertices(2) .ne. 0) then
        p_Dsolution(Ivertices(1)) = p_Dsolution(Ivertices(2))
      else
        p_Dsolution(Ivertices(1)) = 0.0_DP
      end if

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)

    
    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)

    end select

  end subroutine transp_hadaptCallback2d

end module transport_callback2d
