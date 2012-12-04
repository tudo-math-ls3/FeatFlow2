!##############################################################################
!# ****************************************************************************
!# <name> extconvection </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# Implements extended matrix/vector assembly methods for the Navier-Stokes
!# operator.
!#
!# 1.) conv_strdiffOptC2dinitasm
!#     -> Initialises the assembly
!#
!# 2.) conv_strdiffOptC2ddoneasm
!#     -> Cleans up after an assembly
!#
!# 3.) conv_strdiffOptC2dgetMatrix
!#     -> Calculate the matrix of the operator.
!#
!# 4.) conv_strdiffOptC2dgetDerMatrix
!#     -> Calculate the matrix of the derivative of the operator.
!#
!# 5.) conv_strdiffOptC2dgetDefect
!#     -> Calculate the defect.
!#
!# </purpose>
!##############################################################################

module extconvection

  use fsystem
  use storage
  use linearalgebra
  use derivatives
  use element
  use cubature
  use basicgeometry
  use dofmapping
  use triangulation
  use transformation
  use elementpreprocessing
  use spatialdiscretisation
  use bilinearformevaluation
  use linearsystemscalar
  use linearsystemblock
  use collection
  
  implicit none
  
  private
  
  public :: t_navStOperator
  public :: t_navStAssemblyInfo
  public :: conv_strdiffOptC2dinitasm
  public :: conv_strdiffOptC2ddoneasm
  public :: conv_strdiffOptC2dgetMatrix
  public :: conv_strdiffOptC2dgetDerMatrix
  public :: conv_strdiffOptC2dgetDefect

  ! Defines the terms in the Optimal-Control operator to compute.

  type t_navStOperator
  
    ! Stabilisation parameter for the primal equation.
    ! Note: A value of 0.0_DP sets up the convection part without
    ! any stabilisation (central-difference like discretisation).
    real(DP) :: dupsam = 0.0_DP

    ! Model for the viscosity.
    ! =0: Constant viscosity.
    ! =1: Power law: nu = nu_0 * z^(dviscoexponent/2 - 1), nu_0 = 1/RE, z=||D(u)||^2+dviscoEps
    ! =2: Bingham fluid: nu = nu_0 + dviscoyield / sqrt(|D(u)||^2+dviscoEps^2), nu_0 = 1/RE
    integer :: cviscoModel = 0
    
    ! Viscosity parameter $\nu = 1/Re$ if viscosity is constant
    real(dp) :: dnu = 0.0_DP
    
    ! Coefficient in front of the Mass matrix in the primal equation
    real(dp) :: dalpha = 0.0_DP
  
    ! Coefficient in front of the Stokes operator nu*Laplace in the primal equation
    real(dp) :: dbeta = 0.0_DP
    
    ! Coefficient in front of the convection operator \grad(.) y in the
    ! primal equation.
    real(dp) :: ddelta = 0.0_DP

    ! Coefficient in front of the transposed convection operator \grad(.)^t y
    ! in the primal equation.
    real(dp) :: ddeltaTrans = 0.0_DP
  
    ! Coefficient in front of the Newton operator \grad(y)(.) in the
    ! primal equation.
    real(dp) :: dnewton = 0.0_DP
  
    ! Coefficient in front of the Newton operator \grad(y)^t(.) in the
    ! primal equation.
    real(dp) :: dnewtonTrans = 0.0_DP

    ! Calculation of local H.
    ! In 3D, there are 2 different methods to calculate the local H which
    ! is needed for the assembly of the convective part.
    ! =0: Use the cube-root of the volume of the hexahedron as local H
    ! =1: Use the length of the way that a particle travels through
    !     the hexahedron in direction of the flow
    integer :: clocalH = 1

  end type

  ! Assembly structure that saves information about the current
  ! status during the assembly of the local/global matrices.

  type t_navStAssemblyInfo
  
    ! Block discretisation structure that defines the discretisation
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    
    ! Id of the element distribution to be processed.
    integer :: ielemDistribution
    
    ! Pointer to the current element distribution
    type(t_elementDistribution), pointer :: p_relementDistribution

    ! Derivative specifier for evaluating elements
    logical, dimension(EL_MAXNDER) :: Bder
    
    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), pointer :: Idofs
    
    ! Allocateable arrays for the values of the basis functions - 
    ! for test and trial spaces.
    real(DP), dimension(:,:,:,:), pointer :: Dbas

    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef

    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(CUB_MAXCUBP) :: Domega

    ! Type of transformation from the reference to the real element 
    integer :: ctrafoType

    ! An element evaluation set for evaluating elements.
    type(t_evalElementSet) :: revalElementSet
    
    ! Specifier to signal when cubature points are initialised.
    logical :: bcubPtsInitialised
    
    ! Number of local DOF's
    integer :: indof
    
    ! Number of cubature points on the reference element
    integer :: ncubp
    
    ! Maximum number of simultaneously processed elements per element set
    integer :: nelementsPerBlock
    
    ! Local matrices, used during the assembly.
    ! Values and positions of values in the global matrix.
    integer, dimension(:,:,:), pointer :: Kentry
    
    ! Additional contributions for the submatrices A11, A12, A21, A22 stemming from Newton.
    integer, dimension(:,:,:), pointer :: Kentry12
    
    ! Maximum velocity magnitude
    real(DP) :: dumax
    
    ! An array with local DELTA's, each DELTA for one element, primal equation
    real(DP), dimension(:), pointer :: DlocalDelta

    ! An array with local DELTA's, each DELTA for one element, dual equation
    real(DP), dimension(:), pointer :: DlocalDeltaDual
    
    ! A pointer to an element-number list
    integer(I32), dimension(:), pointer :: p_IelementList

    ! Pointer to the primal/dual velocity field in the cubature points.
    real(DP), dimension(:,:,:), pointer :: Dvel
    
    ! Pointer to the velocity X- and Y-derivative of the velocity
    ! in the cubature points
    real(DP), dimension(:,:,:), pointer :: DvelXderiv,DvelYderiv
    
    ! Pointer to the primal velocity DOF's on the elements where to
    ! evaluate the matrices.
    real(DP), dimension(:,:,:), pointer :: DvelDofs

    ! Pointer to the primal velocity DOF's on the elements which are
    ! multiplied to the local matrices to get the derivative of the
    ! operator.
    real(DP), dimension(:,:,:), pointer :: DvelDofsAlt
    
    ! Viscosity in all points on all elements.
    real(dp), dimension(:,:), pointer :: Dnu

    ! Temp memory for ||D(u)||^2 if needed.
    real(dp), dimension(:,:), pointer :: DoperatorNorm
    
  end type

contains

  ! ***************************************************************************

  !<subroutine>

  subroutine computeLocalNu (Dvel,DvelXderiv,DvelYderiv,Dnu,roperator,rassemblyInfo,rcollection)
      
  !<description>
    ! Computes the local viscolity parameter.
    !
    ! Note: For constant viscosity, the routine should do nothing as Dnu is
    ! already initialised with the default nu.
  !</description>
      
  !<input>      
    ! Temporary arrays for the velocity and their
    ! derivatives in the cubature points.
    real(DP), dimension(:,:,:), intent(inout)    :: Dvel
    real(DP), dimension(:,:,:), intent(inout)    :: DvelXderiv
    real(DP), dimension(:,:,:), intent(inout)    :: DvelYderiv

    ! Configuration of the operator
    type(t_navStOperator), intent(in)          :: roperator
    
    ! Collection structure with application settings.
    type(t_collection), intent(inout), optional :: rcollection
  !</input>
  
  !<inputoutput>
    ! Assembly information structure
    type(t_navStAssemblyInfo), intent(inout)   :: rassemblyInfo

    ! Viscosity in the cubature points. To be computed.
    real(DP), dimension(:,:), intent(inout)      :: Dnu
  !</inputoutput>
    
  !<subroutine>
    ! local variables
    integer :: cviscoModel,i,j,isubEquation
    real(DP) :: dnu0,dviscoexponent,dviscoEps,dviscoYield
    type(t_vectorBlock), pointer :: p_rvector
    integer, dimension(2) :: Ibounds
    real(DP), dimension(:,:,:), allocatable :: Ddata
    
    if (.not. present(rcollection)) then
      ! Leave Dnu as it is.
      return
    end if
    
    ! Get the parameters from the collection,
    ! as specified in the call to the SD assembly routine.
    cviscoModel = rcollection%IquickAccess(1)
    isubEquation = rcollection%IquickAccess(2)
    dnu0 = rcollection%DquickAccess(1)
    dviscoexponent = rcollection%DquickAccess(2) 
    dviscoEps = rcollection%DquickAccess(3)
    dviscoYield = rcollection%DquickAccess(4)
    
    ! p_rvector may point to NULL if there is no nonlinearity
    
    ! If we have a nonlinear viscosity, calculate 
    !    z := D(u):D(u) = ||D(u)||^2
    ! The isubEquation defines the shape of the tensor, which may me
    !    D(u) = grad(u) 
    ! or D(u) = 1/2 ( grad(u) + grad(u)^T )
    select case (cviscoModel)
    case (0)
      ! Constant viscosity. There's nothing to do as Dnu is already initialised.
      
    case (1)
    
      ! Calculate ||D(u)||^2 to DoperatorNorm(:,:):
      select case (isubequation)
      case (0)
        ! D(u) = grad(u)
        do i=1,Ibounds(2)
          do j=1,Ibounds(1)
            rassemblyInfo%DoperatorNorm(j,i) = &
                DvelXderiv(1,j,i)**2 + DvelXderiv(2,j,i)**2 + &
                DvelYderiv(1,j,i)**2 + DvelYderiv(2,j,i)**2
          end do
        end do
        
      case (1)
        ! D(u) = 1/2 ( grad(u) + grad(u)^T )
        do i=1,Ibounds(2)
          do j=1,Ibounds(1)
            rassemblyInfo%DoperatorNorm(j,i) = (DvelXderiv(1,j,i)**2 + &
                            0.5_DP * (DvelXderiv(2,j,i) + DvelYderiv(1,j,i))**2 + &
                            DvelYderiv(2,j,i)**2)
          end do
        end do
        
      case default
      
        Ddata(:,:,:) = 0.0_DP
        
      end select
      
      ! Calculate the viscosity
      select case (cviscoModel)
      case (1)
        ! Power law: nu = nu_0 * z^(dviscoexponent/2 - 1), nu_0 = 1/RE, z=||D(u)||^2+dviscoEps
        Dnu(:,:) = dnu0 * (rassemblyInfo%DoperatorNorm(j,i)+dviscoEps)**(0.5_DP*dviscoexponent-1.0_DP)

      case (2)
        ! Bingham fluid: nu = nu_0 + dviscoyield / sqrt(|D(u)||^2+dviscoEps^2), nu_0 = 1/RE
        Dnu(:,:) = dnu0 + dviscoyield / sqrt(rassemblyInfo%DoperatorNorm(j,i)+dviscoEps**2)

      end select
      
      ! Deallocate needed memory.
      deallocate(Ddata)

!    case default
!    
!      ! Viscosity specified by the callback function getNonconstantViscosity.
!      ! Call it and pass the user-defined collection.
!      call getNonconstantViscosity (cterm,rdiscretisation, &
!          nelements,npointsPerElement,Dpoints, &
!          IdofsTest,rdomainIntSubset, &
!          Dnu,p_rvector,rcollection%p_rnextCollection)
    
    end select

  end subroutine
    
  ! ***************************************************************************

  !<subroutine>

  subroutine computeLocalNavStMatrices (Dbas,Domega,Ddetj,ndof,ncubp,&
      Ielements,DvelDofs,&
      Dvel,DvelXderiv,DvelYderiv,Dnu,&
      DentryA11,DentryA22,DentryA12,DentryA21,&
      DlocalDelta,roperator,rassemblyInfo,rcollection)
      
  !<description>
    ! Computes a local matrix of the optimal control operator
    ! which is to be incorporated into the global matrix.
  !</description>
      
  !<input>      
  
    ! Basis functions in all cubature points
    real(DP), dimension(:,:,:,:), intent(in)  :: Dbas
    
    ! Cubature weights
    real(DP), dimension(:), intent(in)        :: Domega
    
    ! Jacobian determinants in all cubature points
    real(DP), dimension(:,:), intent(in)      :: Ddetj
    
    ! Number of local DOF's in each element
    integer, intent(in)                       :: ndof
    
    ! Number of cubature points on each element
    integer, intent(in)                       :: ncubp
    
    ! List of currently handled elements
    integer, dimension(:), intent(in)         :: Ielements
    
    ! DOF's in the velocity
    real(DP), dimension(:,:,:), intent(in)    :: DvelDofs

    ! Temporary arrays for the velocity and their
    ! derivatives in the cubature points.
    real(DP), dimension(:,:,:), intent(inout)    :: Dvel
    real(DP), dimension(:,:,:), intent(inout)    :: DvelXderiv
    real(DP), dimension(:,:,:), intent(inout)    :: DvelYderiv

    ! Configuration of the operator
    type(t_navStOperator), intent(in)          :: roperator
    
    ! Assembly information structure
    type(t_navStAssemblyInfo), intent(inout)   :: rassemblyInfo
    
    ! Collection structure with application settings.
    type(t_collection), intent(inout), optional :: rcollection
  !</input>
    
  !<output>
    
    ! Viscosity in the cubature points.
    real(DP), dimension(:,:), intent(inout)      :: Dnu
    
    ! Temporary space for local delta of the SD method, primal equation
    real(DP), dimension(:), intent(inout)        :: DlocalDelta

    ! Entries in the matrix.
    ! All local matrices are saved transposed, i.e. we have
    ! DentryIJ(column,row,element). This has to be taken care of
    ! during the assembly and was introduced to gain higher speed
    ! as it avoids extensive jumps in the memory access.
    real(DP), dimension(:,:,:), intent(inout) :: DentryA11
    real(DP), dimension(:,:,:), intent(inout) :: DentryA22
    
    real(DP), dimension(:,:,:), intent(inout) :: DentryA12
    real(DP), dimension(:,:,:), intent(inout) :: DentryA21
  !</output>   
  
  !<subroutine>
    
    ! local variables
    integer :: iel, icubp, idofe, jdofe, NEL
    real(DP) :: OM,dtemp,dsumI,dsumJ
    real(DP) :: du1p,du2p,du1d,du2d, db,dbx,dby
    real(DP) :: du1locxp,du1locyp,du2locxp,du2locyp
    real(DP) :: du1locxd,du1locyd,du2locxd,du2locyd
    real(DP) :: dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY
    real(DP) :: AH11,AH22,AH44,AH55,AH12,AH21,AH45,AH54
    real(DP) :: AH41,AH52,AH42,AH51,AH14,AH25,AH24,AH15
    real(dp) :: dtemp1,dtemp2

    ! Compute the local NU in all cubature points.
    call computeLocalNu (Dvel,DvelXderiv,DvelYderiv,Dnu,roperator,rassemblyInfo,rcollection)

    NEL = size(Ielements)
    
    ! ===========================================
    ! Solution evaluation phase
    !
    ! Evaluate the solutions for the nonlinearity
    ! ===========================================

    ! Calculate local DELTA's for streamline diffusion method.
    ! (cf. p. 121 in Turek's CFD book).
    ! For every element, we need a local DELTA.
    ! Every local delta is weighted by the global "ddelta".
    ! If ddelta=0, we don't do anything as this disables the
    ! nonlinear term.
    ! If UPSAM=0.0, we have a central-difference like discretisation, which
    ! is one can see as the local stabilisation weight Delta is also = 0.0.
    ! In this case, we even switch of the calculation of the local Delta,
    ! as it is always =0.0, so we save a little bit time.
    if (roperator%dupsam .ne. 0.0_DP) then
      call getLocalDeltaQuad (DvelDofs,rassemblyInfo%p_rdiscretisation%p_rtriangulation,&
          Ielements,rassemblyInfo%dumax,roperator%dupsam,&
          rassemblyInfo%Dnu,DlocalDelta)
    end if

    ! Loop over all elements in the current set
    do iel=1,NEL
    
      ! Loop over all cubature points on the current element
      do icubp = 1, ncubp
      
        ! Primal/dual velocity
        du1p = 0.0_DP
        du2p = 0.0_DP
        
        ! X/Y-derivative of the primal/dual velocity
        du1locxp = 0.0_DP
        du1locyp = 0.0_DP
        du2locxp = 0.0_DP
        du2locyp = 0.0_DP

        ! Perform a loop through the trial DOF's.
        do jdofe=1,ndof

          ! Get the value of the (test) basis function 
          ! phi_i (our "O") in the cubature point:
          
          db =  Dbas(jdofe,1,icubp,iel)
          dbx = Dbas(jdofe,DER_DERIV2D_X,icubp,iel)
          dby = Dbas(jdofe,DER_DERIV2D_Y,icubp,iel)
          
          ! Sum up to the value in the cubature point
          
          du1p = du1p + DvelDofs(jdofe,iel,1)*db
          du2p = du2p + DvelDofs(jdofe,iel,2)*db
          
          du1locxp = du1locxp + DvelDofs(jdofe,iel,1)*dbx
          du1locyp = du1locyp + DvelDofs(jdofe,iel,1)*dby
          du2locxp = du2locxp + DvelDofs(jdofe,iel,2)*dbx
          du2locyp = du2locyp + DvelDofs(jdofe,iel,2)*dby
                                                    
        end do ! jdofe
        
        ! Save the computed velocities and derivatives
        Dvel(1,icubp,iel) = du1p
        Dvel(2,icubp,iel) = du2p

        DvelXderiv(1,icubp,iel) = du1locxp
        DvelXderiv(2,icubp,iel) = du1locyp
        DvelYderiv(1,icubp,iel) = du2locxp
        DvelYderiv(2,icubp,iel) = du2locyp
      
      end do ! icubp
      
    end do ! iel
    
    ! ===========================================
    ! Matrix evaluation phase
    !
    ! Evaluate the matrices using the solutions
    ! from above.
    ! ===========================================
      
    AH11 = 0.0_DP
    AH22 = 0.0_DP
    
    AH12 = 0.0_DP
    AH21 = 0.0_DP
    
    ! Compute mass matrix part?
    if (roperator%dalpha .ne. 0.0_DP) then
    
      ! Loop over the elements in the current set.
      do iel=1,NEL

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant 
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

          OM = Domega(icubp)*Ddetj(icubp,iel)

          ! Outer loop over the DOF's i=1..indof on our current element, 
          ! which corresponds to the basis functions Phi_i:

          do idofe=1,ndof
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the 
            ! current DOF into HBASIy:
          
            dbasI = Dbas(idofe,1,icubp,iel)
            
            ! Inner loop over the DOF's j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do jdofe=1,ndof
              
              ! Fetch the contributions of the (trial) basis function Phi_j
              ! (out "X") for function value and first derivatives for the 
              ! current DOF into HBASJy:
            
              dbasJ = Dbas(jdofe,1,icubp,iel)

              ! Contribution to the mass matrix -- weighted by the cubature weight
              dtemp = OM * dbasI * dbasJ

              ! Weighten the calculated value AHxy by the cubature
              ! weight OM and add it to the local matrices. After the
              ! loop over all DOF's is finished, each entry contains
              ! the calculated integral.

              DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel) + &
                                           roperator%dalpha * dtemp
                                           
              DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel) + &
                                           roperator%dalpha * dtemp
                                           
            end do ! idofe
            
          end do ! jdofe

        end do ! icubp 
      
      end do ! iel
    end if
    
    ! Laplace matrix
    if (roperator%dbeta .ne. 0.0_DP) then
    
      ! Loop over the elements in the current set.
      do iel=1,NEL

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant 
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

          OM = Domega(icubp)*Ddetj(icubp,iel)

          ! Outer loop over the DOF's i=1..indof on our current element, 
          ! which corresponds to the basis functions Phi_i:

          do idofe=1,ndof
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the 
            ! current DOF into HBASIy:
          
            dbasIX = Dbas(idofe,2,icubp,iel)
            dbasIY = Dbas(idofe,3,icubp,iel)
            
            ! Inner loop over the DOF's j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do jdofe=1,ndof
              
              ! Fetch the contributions of the (trial) basis function Phi_j
              ! (out "X") for function value and first derivatives for the 
              ! current DOF into HBASJy:
            
              dbasJX = Dbas(jdofe,2,icubp,iel)
              dbasJY = Dbas(jdofe,3,icubp,iel)

              ! dny*(grad(phi_j,grad(phi_i))
              dtemp = rassemblyInfo%Dnu(icubp,iel) * (dbasIX*dbasJX + dbasIY*dbasJY)

              ! Weighten the calculated value AHxy by the cubature
              ! weight OM and add it to the local matrices. After the
              ! loop over all DOF's is finished, each entry contains
              ! the calculated integral.

              DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel) + &
                  OM * roperator%dbeta * dtemp
                  
              DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel) + &
                  OM * roperator%dbeta * dtemp
                  
            end do ! idofe
            
          end do ! jdofe

        end do ! icubp 
      
      end do ! iel
    end if    

    ! Convection operator grad(.)*y.
    ! This operator implies the transposed Newton \grad(y)^t*(.) in the
    ! dual equation.
    if (roperator%ddelta .ne. 0.0_DP) then
    
      ! Loop over the elements in the current set.
      do iel=1,NEL

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant 
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

          OM = Domega(icubp)*Ddetj(icubp,iel)
          
          du1p = Dvel(1,icubp,iel)
          du2p = Dvel(2,icubp,iel)
          du1locxp = DvelXderiv(1,icubp,iel)
          du1locyp = DvelXderiv(2,icubp,iel)
          du2locxp = DvelYderiv(1,icubp,iel)
          du2locyp = DvelYderiv(2,icubp,iel)

          ! Outer loop over the DOF's i=1..indof on our current element, 
          ! which corresponds to the basis functions Phi_i:

          do idofe=1,ndof
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the 
            ! current DOF into HBASIy:
          
            dbasI = Dbas(idofe,1,icubp,iel)
            dbasIX = Dbas(idofe,2,icubp,iel)
            dbasIY = Dbas(idofe,3,icubp,iel)
            
            dsumI = dbasIX*du1p + dbasIY*du2p
            
            ! Inner loop over the DOF's j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do jdofe=1,ndof
              
              ! Fetch the contributions of the (trial) basis function Phi_j
              ! (out "X") for function value and first derivatives for the 
              ! current DOF into HBASJy:
            
              dbasJ = Dbas(jdofe,1,icubp,iel)
              dbasJX = Dbas(jdofe,2,icubp,iel)
              dbasJY = Dbas(jdofe,3,icubp,iel)

              dsumJ = dbasJX*du1p + dbasJY*du2p

              ! Transposed Newton in the dual equation.            
              dtemp = dbasJ*dbasI
            
              ! Delta*(U*grad(Phi_j), U*grad(Phi_i)) + (U*grad(Phi_j),Phi_i)
              dtemp1 = dsumJ * (DlocalDelta(iel)*dsumI + dbasI)
              
              ! Weighten the calculated value AHxy by the cubature
              ! weight OM and add it to the local matrices. After the
              ! loop over all DOF's is finished, each entry contains
              ! the calculated integral.

              DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel) + &
                  OM * roperator%ddelta * dtemp1
                  
              DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel) + &
                  OM * roperator%ddelta * dtemp1
                  
            end do ! idofe
            
          end do ! jdofe

        end do ! icubp 
      
      end do ! iel
    end if    
    
    ! Newton operator grad(y)*(.).
    ! This operator implies the reactive operator (.)*grad(\lambda)+(\grad(.)^t)\lambda
    ! in the dual equation.
    if (roperator%dnewton .ne. 0.0_DP) then
    
      ! Loop over the elements in the current set.
      do iel=1,NEL

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant 
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

          OM = Domega(icubp)*Ddetj(icubp,iel)
          
          du1locxp = DvelXderiv(1,icubp,iel)
          du1locyp = DvelXderiv(2,icubp,iel)
          du2locxp = DvelYderiv(1,icubp,iel)
          du2locyp = DvelYderiv(2,icubp,iel)
          
          ! Outer loop over the DOF's i=1..indof on our current element, 
          ! which corresponds to the basis functions Phi_i:

          do idofe=1,ndof
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the 
            ! current DOF into HBASIy:
          
            dbasI = Dbas(idofe,1,icubp,iel)
            dbasIX = Dbas(idofe,2,icubp,iel)
            dbasIY = Dbas(idofe,3,icubp,iel)
            
            ! Inner loop over the DOF's j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do jdofe=1,ndof
              
              ! Fetch the contributions of the (trial) basis function Phi_j
              ! (out "X") for function value and first derivatives for the 
              ! current DOF into HBASJy:
            
              dbasJ = Dbas(jdofe,1,icubp,iel)
              dbasJX = Dbas(jdofe,2,icubp,iel)
              dbasJY = Dbas(jdofe,3,icubp,iel)

              ! Delta*(U*grad(Phi_j), U*grad(Phi_i)) + (U*grad(Phi_j),Phi_i)
              dtemp = dbasJ*dbasI
              dtemp1 = dbasJX*dbasI
              dtemp2 = dbasJY*dbasI

              ! Weighten the calculated value AHxy by the cubature
              ! weight OM and add it to the local matrices. After the
              ! loop over all DOF's is finished, each entry contains
              ! the calculated integral.

              ! Newton operator in the primal
              DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel) + &
                  OM * roperator%dnewton * du1locxp * dtemp
                  
              DentryA12(jdofe,idofe,iel) = DentryA12(jdofe,idofe,iel) + &
                  OM * roperator%dnewton * du1locyp * dtemp
                  
              DentryA21(jdofe,idofe,iel) = DentryA21(jdofe,idofe,iel) + &
                  OM * roperator%dnewton * du2locxp * dtemp
                  
              DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel) + &
                  OM * roperator%dnewton * du2locyp * dtemp

            end do ! idofe
            
          end do ! jdofe

        end do ! icubp 
      
      end do ! iel
    end if    
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_strdiffOptC2dinitasm (rdiscretisation,ielemDistribution,&
      roperator,rassemblyInfo,balternativeVelocity)
  
!<description>
  ! Initialise the rassemblyInfo structure for the assembly of
  ! local matrices. Initialises the Dnu array with the default viscosity.
!</description>

!<input>
  ! Block discretisation structure that defines the discretisation
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisation
  
  ! Id of the element distribution to be processed.
  integer, intent(in) :: ielemDistribution
  
  ! Configuration of the operator
  type(t_navStOperator), intent(in) :: roperator

  ! OPTIONAL: Create arrays for an additional velocity, independent of the
  ! evaluation point of the matrices. If not specified, FALSE is assumed.
  logical, intent(in), optional :: balternativeVelocity
!</input>
  
!<output>
  ! Structure collecting information about the assembly.
  type(t_navStAssemblyInfo), intent(out) :: rassemblyInfo
!</output>

!</subroutine>

    ! local variables
    type(t_elementDistribution), pointer :: p_relementDistribution
    integer :: i,k
    
    ! Cubature point coordinates on the reference element
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi

    ! Get the element distribution. Here, we can find information about
    ! the cubature formula etc...
    p_relementDistribution => rdiscretisation%RelementDistr(ielemDistribution)
    
    ! Remember the discretisation and the element distribution
    rassemblyInfo%p_rdiscretisation => rdiscretisation
    rassemblyInfo%p_relementDistribution => p_relementDistribution
    rassemblyInfo%ielemDistribution = ielemDistribution
    
    ! Get the number of local DOF's for trial/test functions.
    ! We assume trial and test functions to be the same.
    rassemblyInfo%indof = elem_igetNDofLoc(p_relementDistribution%celement)

    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    rassemblyInfo%nelementsPerBlock = min(BILF_NELEMSIM,rdiscretisation%p_rtriangulation%NEL)
    
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    rassemblyInfo%ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)
    
    ! Allocate some memory to hold the cubature points on the reference element
    allocate(rassemblyInfo%p_DcubPtsRef(&
        trafo_igetReferenceDimension(rassemblyInfo%ctrafoType),CUB_MAXCUBP))

    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference element
    call cub_getCubPoints(rassemblyInfo%p_relementDistribution%ccubTypeBilForm, &
        rassemblyInfo%ncubp, Dxi, rassemblyInfo%Domega)
    
    ! Reformat the cubature points; they are in the wrong shape!
    do i=1,rassemblyInfo%ncubp
      do k=1,ubound(rassemblyInfo%p_DcubPtsRef,1)
        rassemblyInfo%p_DcubPtsRef(k,i) = Dxi(i,k)
      end do
    end do
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  allocate(Dbas(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    
    allocate(rassemblyInfo%Dbas(&
        rassemblyInfo%indof,elem_getMaxDerivative(p_relementDistribution%celement), &
        rassemblyInfo%ncubp,rassemblyInfo%nelementsPerBlock))
        
    ! Allocate memory for the DOF's of all the elements.
    allocate(rassemblyInfo%Idofs(rassemblyInfo%indof,rassemblyInfo%nelementsPerBlock))
    
    ! Allocate memory for array with local DELTA's and nu's.
    allocate(rassemblyInfo%DlocalDelta(rassemblyInfo%nelementsPerBlock))
    allocate(rassemblyInfo%DoperatorNorm(rassemblyInfo%ncubp,rassemblyInfo%nelementsPerBlock))
    allocate(rassemblyInfo%Dnu(rassemblyInfo%ncubp,rassemblyInfo%nelementsPerBlock))
    
    ! Initialise nu with the default viscosity.
    rassemblyInfo%Dnu(:,:) = roperator%dnu

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
    !
    ! Kentry (:,:,:) defines the positions of the local matrices
    ! in the submatrices A11 and A22.
    ! KentryA12 (:,:,:) defines the positions of the local matrices
    ! in the submatrices A12 and A21.
    allocate(rassemblyInfo%Kentry(&
        rassemblyInfo%indof,rassemblyInfo%indof,rassemblyInfo%nelementsPerBlock))
    
    allocate(rassemblyInfo%Kentry12(&
        rassemblyInfo%indof,rassemblyInfo%indof,rassemblyInfo%nelementsPerBlock))
    
    ! Allocate memory for the velocites in the cubature points.
    allocate(rassemblyInfo%Dvel(NDIM2D,rassemblyInfo%ncubp,&
        rassemblyInfo%nelementsPerBlock))
    allocate(rassemblyInfo%DvelXderiv(NDIM2D,rassemblyInfo%ncubp,&
        rassemblyInfo%nelementsPerBlock))
    allocate(rassemblyInfo%DvelYderiv(NDIM2D,rassemblyInfo%ncubp,&
        rassemblyInfo%nelementsPerBlock))
    
    ! Allocate memory for the DOF's on the elements.
    allocate (rassemblyInfo%DvelDofs(&
        rassemblyInfo%indof,rassemblyInfo%nelementsPerBlock,NDIM2D))
       
    ! The alternative velocity coincides with the standard velocity.
    rassemblyInfo%DvelDofsAlt => rassemblyInfo%DvelDofs

    if (present(balternativeVelocity)) then
      if (balternativeVelocity) then
        ! Create an additional evaluation point for the velocity.
        allocate (rassemblyInfo%DvelDofsAlt(&
            rassemblyInfo%indof,rassemblyInfo%nelementsPerBlock,NDIM2D))
      end if
    end if
    
    ! Indicate that cubature points must still be initialised in the element set.
    rassemblyInfo%bcubPtsInitialised = .false.
    
    ! Initialise the derivative flags
    rassemblyInfo%Bder = .false.
    rassemblyInfo%Bder(DER_FUNC) = .true.
    rassemblyInfo%Bder(DER_DERIV2D_X) = .true.
    rassemblyInfo%Bder(DER_DERIV2D_Y) = .true.
    
    ! Initialisation of the element set.
    call elprep_init(rassemblyInfo%revalElementSet)

    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_relementDistribution%h_IelementList, &
                              rassemblyInfo%p_IelementList)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_strdiffOptC2dinitelemset (rvelmatrix,rvelmatrixoffdiag,&
      rvelocityVector,roperator,rassemblyInfo,istartElement,iendElement,&
      rvelocityVectorAlt)
  
!<description>
  ! Initialise the rassemblyInfo structure for the assembly of
  ! local matrices for the elements istartElement..istartElement+ielCountMax-1 (of the
  ! current element set). 
!</description>

!<input>
  ! Template FE matrix specifying the connectivity of the primal/dual
  ! velocity space for diagonal blocks.
  type(t_matrixScalar), intent(in) :: rvelmatrix

  ! Template FE matrix specifying the connectivity of the primal/dual
  ! velocity space for offdiagonal blocks.
  type(t_matrixScalar), intent(in) :: rvelmatrixoffdiag
  
  ! Velocity vector, primal equation (Block 1/2) for the evaluation of the matrices
  type(t_vectorBlock), intent(in) :: rvelocityVector

  ! Structure defining the operator to set up.
  type(t_navStOperator), intent(in) :: roperator

  ! Index of the start element of the current element set in rassemblyInfo.
  integer, intent(in) :: istartElement
  
  ! Index of the last element.
  integer, intent(in) :: iendElement

  ! OPTIONAL: Alternative Velocity vector (Block 1/2).
  type(t_vectorBlock), intent(in), optional :: rvelocityVectorAlt

!</input>

!<inputoutput>
  ! Structure collecting information about the assembly.
  type(t_navStAssemblyInfo), intent(inout) :: rassemblyInfo
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel, idofe, jdofe, jcol0, jcol, jdfg
    integer, dimension(:), pointer :: p_Kcol
    integer, dimension(:), pointer :: p_Kld
    integer, dimension(:), pointer :: p_Kcol12
    integer, dimension(:), pointer :: p_Kld12

    ! Primal and dual velocity
    real(DP), dimension(:), pointer :: p_DvelocityX, p_DvelocityY
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag
    
    ! Get pointers to the matrix structure(s).
    call lsyssc_getbase_Kcol (rvelmatrix,p_Kcol)
    call lsyssc_getbase_Kld (rvelmatrix,p_Kld)
    
    call lsyssc_getbase_Kcol (rvelmatrixoffdiag,p_Kcol12)
    call lsyssc_getbase_Kld (rvelmatrixoffdiag,p_Kld12)
   
    ! The outstanding feature with finite elements is: A basis
    ! function for a DOF on one element has common support only
    ! with the DOF's on the same element! E.g. for Q1:
    !
    !        #. . .#. . .#. . .#
    !        .     .     .     .
    !        .  *  .  *  .  *  .
    !        #-----O-----O. . .#
    !        |     |     |     .
    !        |     | iel |  *  .
    !        #-----X-----O. . .#
    !        |     |     |     .
    !        |     |     |  *  .
    !        #-----#-----#. . .#
    !
    ! --> On element iel, the basis function at "X" only interacts
    !     with the basis functions in "O". Elements in the 
    !     neighbourhood ("*") have no support, therefore we only have
    !     to collect all "O" DOF's.
    !
    ! Calculate the global DOF's into IdofsTrial / IdofsTest.
    !
    ! More exactly, we call dof_locGlobMapping_mult to calculate all the
    ! global DOF's of our BILF_NELEMSIM elements simultaneously.
    call dof_locGlobMapping_mult(rassemblyInfo%p_rdiscretisation, &
        rassemblyInfo%p_IelementList(istartElement:iendElement), rassemblyInfo%Idofs)
                                
    ! For the assembly of the global matrix, we use a "local"
    ! approach. At first we build a "local" system matrix according
    ! to the current element. This contains all additive
    ! contributions of element iel, which are later added at the
    ! right positions to the elements in the global system matrix.
    !
    ! We have indofTrial trial DOF's per element and
    ! indofTest test DOF's per element. Therefore there are
    ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j) 
    ! "active" (i.e. have common support) on our current element, each 
    ! giving an additive contribution to the system matrix.
    !
    ! We build a quadratic indofTrial*indofTest local matrix:
    ! Kentry(1..indofTrial,1..indofTest) receives the position 
    !   in the global system matrix, where the corresponding value 
    !   has to be added to.
    ! (The corresponding contrbutions can be saved separately, 
    !  but we directly add them to the global matrix in this 
    !  approach.)
    !
    ! We build local matrices for all our elements 
    ! in the set simultaneously.
    ! Loop through elements in the set and for each element,
    ! loop through the local matrices to initialise them:
    do iel=1,iendelement-istartElement+1
    
      ! For building the local matrices, we have first to
      ! loop through the test functions (the "O"'s), as these
      ! define the rows in the matrix.
      do idofe=1,rassemblyInfo%indof
      
        ! Row idofe of the local matrix corresponds 
        ! to row=global DOF KDFG(idofe) in the global matrix.
        ! This is one of the the "O"'s in the above picture.
        ! Get the starting position of the corresponding row
        ! to jcol0:

        jcol0=p_KLD(rassemblyInfo%Idofs(idofe,iel))
        
        ! Now we loop through the other DOF's on the current element
        ! (the "O"'s).
        ! All these have common support with our current basis function
        ! and will therefore give an additive value to the global
        ! matrix.
        
        do jdofe=1,rassemblyInfo%indof
          
          ! Get the global DOF of the "X" which interacts with 
          ! our "O".
          
          jdfg=rassemblyInfo%Idofs(jdofe,iel)
          
          ! Starting in jcol0 (which points to the beginning of
          ! the line initially), loop through the elements in
          ! the row to find the position of column IDFG.
          ! Jump out of the do loop if we find the column.
          
          do jcol=jcol0,rvelmatrix%NA
            if (p_KCOL(jcol) .eq. jdfg) exit
          end do
          
          ! Because columns in the global matrix are sorted 
          ! ascendingly (except for the diagonal element),
          ! the next search can start after the column we just found.
          
          ! jcol0=jcol+1
          
          ! Save the position of the matrix entry into the local
          ! matrix.
          ! Note that a column in Kentry corresponds to a row in
          ! the real matrix. We aligned Kentry/DENTRY this way to get
          ! higher speed of the assembly routine, since this leads
          ! to better data locality.
          
          rassemblyInfo%Kentry(jdofe,idofe,iel)=jcol
          
        end do ! idofe
        
      end do ! jdofe
      
    end do ! iel
    
    ! If the Newton part is to be calculated, we also need the matrix positions
    ! in A12 and A21. We can skip this part if the column structure is
    ! exactly the same!
    if (associated(p_Kcol,p_Kcol12)) then
    
      rassemblyInfo%Kentry12(:,:,:) = rassemblyInfo%Kentry(:,:,:)
      
    else

      do iel=1,iendelement-istartElement+1
      
        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"'s), as these
        ! define the rows in the matrix.
        do idofe=1,rassemblyInfo%indof
        
          ! Row idofe of the local matrix corresponds 
          ! to row=global DOF KDFG(idofe) in the global matrix.
          ! This is one of the the "O"'s in the above picture.
          ! Get the starting position of the corresponding row
          ! to jcol0:

          jcol0=p_KLD12(rassemblyInfo%Idofs(idofe,iel))
          
          ! Now we loop through the other DOF's on the current element
          ! (the "O"'s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.
          
          do jdofe=1,rassemblyInfo%indof
            
            ! Get the global DOF of the "X" which interacts with 
            ! our "O".
            
            jdfg=rassemblyInfo%Idofs(jdofe,iel)
            
            ! Starting in jcol0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the do loop if we find the column.
            
            do jcol=jcol0,rvelmatrixoffdiag%NA
              if (p_KCOL12(jcol) .eq. jdfg) exit
            end do

            ! Because columns in the global matrix are sorted 
            ! ascendingly (except for the diagonal element),
            ! the next search can start after the column we just found.
            
            ! jcol0=jcol+1
            
            ! Save the position of the matrix entry into the local
            ! matrix.
            ! Note that a column in Kentry corresponds to a row in
            ! the real matrix. We aligned Kentry/DENTRY this way to get
            ! higher speed of the assembly routine, since this leads
            ! to better data locality.
            
            rassemblyInfo%Kentry12(jdofe,idofe,iel)=jcol
            
          end do ! idofe
          
        end do ! jdofe
        
      end do ! iel
      
    end if
    
    ! Ok, we found the positions of the local matrix entries
    ! that we have to change.
    ! To calculate the matrix contributions, we have to evaluate
    ! the elements to give us the values of the basis functions
    ! in all the DOF's in all the elements in our set.

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag. 
    cevaluationTag = elem_getEvaluationTag(rassemblyInfo%p_relementDistribution%celement)
                    
    ! In the first loop, calculate the coordinates on the reference element.
    ! In all later loops, use the precalculated information.
    !
    ! Note: Why not using
    !   if (ielset .EQ. 1) then
    ! here, but this strange concept with the boolean variable?
    ! Because the if-command does not work with OpenMP! bcubPtsInitialised
    ! is a local variable and will therefore ensure that every thread
    ! is initialising its local set of cubature points!
    if (.not. rassemblyInfo%bcubPtsInitialised) then
      rassemblyInfo%bcubPtsInitialised = .true.
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)
    else
      cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
    end if
    
    ! Calculate all information that is necessary to evaluate the finite element
    ! on all cells of our subset. This includes the coordinates of the points
    ! on the cells.
    call elprep_prepareSetForEvaluation (rassemblyInfo%revalElementSet,&
        cevaluationTag, rassemblyInfo%p_rdiscretisation%p_rtriangulation, &
        rassemblyInfo%p_IelementList(istartElement:iendelement), &
        rassemblyInfo%ctrafoType, rassemblyInfo%p_DcubPtsRef(:,1:rassemblyInfo%ncubp))

    ! Calculate the values of the basis functions.
    ! Pass p_DcubPts as point coordinates, which point either to the
    ! coordinates on the reference element (the same for all elements)
    ! or on the real element - depending on whether this is a 
    ! parametric or nonparametric element.
    call elem_generic_sim2 (rassemblyInfo%p_relementDistribution%celement, &
        rassemblyInfo%revalElementSet, rassemblyInfo%Bder, rassemblyInfo%Dbas)
          
    ! Calculate the primal and dual velocity and its derivatives
    ! in all the cubature points on all the elements.
    ! The calculation routine might need them.
    
    ! Get pointers to the velocity to be evaluated.
    call lsyssc_getbase_double (rvelocityVector%RvectorBlock(1),p_DvelocityX)
    call lsyssc_getbase_double (rvelocityVector%RvectorBlock(2),p_DvelocityY)
    
    ! Extract the velocity DOF's omn the elements.
    do iel=1,iendelement-istartElement+1
      do idofe = 1,rassemblyInfo%indof
        jdofe = rassemblyInfo%Idofs(idofe,iel)
        rassemblyInfo%DvelDofs(idofe,iel,1) = p_DvelocityX(jdofe)
        rassemblyInfo%DvelDofs(idofe,iel,2) = p_DvelocityY(jdofe)
      end do
    end do

    if (present(rvelocityVectorAlt)) then
      ! Also fetch the alternative velocity.
      call lsyssc_getbase_double (rvelocityVectorAlt%RvectorBlock(1),p_DvelocityX)
      call lsyssc_getbase_double (rvelocityVectorAlt%RvectorBlock(2),p_DvelocityY)

      ! Extract the velocity DOF's omn the elements.
      do iel=1,iendelement-istartElement+1
        do idofe = 1,rassemblyInfo%indof
          jdofe = rassemblyInfo%Idofs(idofe,iel)
          rassemblyInfo%DvelDofsAlt(idofe,iel,1) = p_DvelocityX(jdofe)
          rassemblyInfo%DvelDofsAlt(idofe,iel,2) = p_DvelocityY(jdofe)
        end do
      end do
    end if
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine conv_strdiffOptC2ddoneasm (rassemblyInfo)
  
!<description>
  ! Clean up the rassemblyInfo structure.
!</description>

!<inputoutput>
  ! Structure collecting information about the assembly.
  type(t_navStAssemblyInfo), intent(inout) :: rassemblyInfo
!</inputoutput>

!</subroutine>

    deallocate(rassemblyInfo%Dbas)
    deallocate(rassemblyInfo%Idofs)
    deallocate(rassemblyInfo%DlocalDelta)
    deallocate(rassemblyInfo%DlocalDeltaDual)
    deallocate(rassemblyInfo%Kentry)
    deallocate(rassemblyInfo%Kentry12)
    
    ! Allocate memory for the velocites in the cubature points.
    deallocate(rassemblyInfo%Dvel)
    deallocate(rassemblyInfo%DvelXderiv)
    deallocate(rassemblyInfo%DvelYderiv)
    
    if (associated(rassemblyInfo%DvelDofsAlt,rassemblyInfo%DvelDofs)) then
      deallocate(rassemblyInfo%DvelDofsAlt)
    else
      nullify(rassemblyInfo%DvelDofsAlt)
    end if
    
    deallocate(rassemblyInfo%DvelDofs)
    
    deallocate(rassemblyInfo%Dnu)
    deallocate(rassemblyInfo%DoperatorNorm)
    deallocate(rassemblyInfo%p_DcubPtsRef)

    ! Clean up of the element set.
    call elprep_releaseElementSet(rassemblyInfo%revalElementSet)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  subroutine conv_strdiffOptC2dgetMatrix (rmatrix,roperator,dweight,&
      rvelocityVector,rcollection)
!<description>
  ! Calculate the matrix of the nonlinear operator:
  !   dweight*A(rvelocityVector)
!</description>

!<input>

  ! Structure defining the operator to set up.
  type(t_navStOperator), intent(in) :: roperator
  
  ! Weight of the operator. Standard = 1.0
  real(dp), intent(in) :: dweight
      
  ! Velocity vector for the nonlinearity.
  ! The first blocks 1/2 in this vector define the evaluation
  ! point (velocity).
  type(t_vectorBlock), intent(in) :: rvelocityVector

  ! Collection structure with application settings.
  type(t_collection), intent(inout), optional :: rcollection
!</input>

!<inputoutput>
  ! The system matrix. The submatrices for the velocity must be in block
  ! A11, A12, A21 and A22 and must be in matrix format 7 or 9.
  ! A11 and A22 must have the same structure. A12 and A21 must have
  ! the same structure.
  type(t_matrixBlock), intent(inout), target :: rmatrix
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: ielset,ielmax
  
  ! Matrix structure arrays
  real(DP), dimension(:), pointer :: p_Da11,p_Da22
  real(DP), dimension(:), pointer :: p_Da12,p_Da21
  
  ! Local matrices
  real(DP), dimension(:,:,:), allocatable :: DentryA11
  real(DP), dimension(:,:,:), allocatable :: DentryA12
  real(DP), dimension(:,:,:), allocatable :: DentryA21
  real(DP), dimension(:,:,:), allocatable :: DentryA22
  
  ! Temp vector
  type(t_vectorBlock) :: rvectorBlock
  
  ! Assembly status structure
  type(t_navStAssemblyInfo) :: rassemblyInfo

    if (roperator%dnu .eq. 0.0_DP) then
      print *,'SD: NU=0 not allowed! Set dbeta=0 to prevent Stokes operator'// &
              ' from being build!'
      call sys_halt()
    end if
    
    ! Get pointers to the matrix content 
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1),p_Da11)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,2),p_Da22)

    nullify(p_Da12,p_Da21)
    if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,2),p_Da12)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,1),p_Da21)
    end if

    ! Initialise the asembly of the local matrices.
    call conv_strdiffOptC2dinitasm (rvelocityVector%p_rblockDiscr%RspatialDiscr(1),&
        1,roperator,rassemblyInfo)

    ! Calculate the maximum norm of the actual velocity field
    ! Round up the norm to 1D-8 if it's too small...
    call lsysbl_deriveSubvector(rvelocityVector,rvectorBlock,1,2,.true.)
    call lsysbl_getVectorMagnitude (rvectorBlock,dumax=rassemblyInfo%dumax)
    call lsysbl_releaseVector (rvectorBlock)
  
    if (rassemblyInfo%dumax .lt. 1E-8_DP) rassemblyInfo%dumax=1E-8_DP

    ! Allocate memory for the local matrices
    allocate(DentryA11(rassemblyInfo%indof,rassemblyInfo%indof,rassemblyInfo%nelementsPerBlock))
    allocate(DentryA12(rassemblyInfo%indof,rassemblyInfo%indof,rassemblyInfo%nelementsPerBlock))
    allocate(DentryA21(rassemblyInfo%indof,rassemblyInfo%indof,rassemblyInfo%nelementsPerBlock))
    allocate(DentryA22(rassemblyInfo%indof,rassemblyInfo%indof,rassemblyInfo%nelementsPerBlock))

    ! Initialise the array with the local Delta values for the stabilisation
    call lalg_clearVectorDble (rassemblyInfo%DlocalDelta)
    call lalg_clearVectorDble (rassemblyInfo%DlocalDeltaDual)

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so BILF_NELEMSIM local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !%OMP do SCHEDULE(dynamic,1)
    do ielset = 1, size(rassemblyInfo%p_IelementList), BILF_NELEMSIM

      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      ielmax = min(size(rassemblyInfo%p_IelementList),ielset-1+BILF_NELEMSIM)
    
      ! Initialise the element set, compute the basis functions in the
      ! cubature points.
      if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
        call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
            rmatrix%RmatrixBlock(1,2),rvelocityVector,&
            roperator,rassemblyInfo,ielset,ielmax)
      else
        call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
            rmatrix%RmatrixBlock(1,1),rvelocityVector,&
            roperator,rassemblyInfo,ielset,ielmax)
      end if

      ! Clear the local matrices. If the Newton part is to be calculated,
      ! we must clear everything, otherwise only Dentry.
      DentryA11 = 0.0_DP
      DentryA12 = 0.0_DP
      DentryA21 = 0.0_DP
      DentryA22 = 0.0_DP
      
      ! Now calculate the local matrices on all the elements.
      call computeLocalNavStMatrices (rassemblyInfo%Dbas,&
          rassemblyInfo%Domega,rassemblyInfo%revalElementSet%p_Ddetj,&
          rassemblyInfo%indof,rassemblyInfo%ncubp,&
          rassemblyInfo%p_IelementList(ielset:ielmax),&
          rassemblyInfo%DvelDofs,rassemblyInfo%Dvel,rassemblyInfo%DvelXderiv,&
          rassemblyInfo%DvelYderiv,rassemblyInfo%Dnu,&
          DentryA11,DentryA22,DentryA12,DentryA21,&
          rassemblyInfo%DlocalDelta,roperator,rassemblyInfo,rcollection)
      
      ! Incorporate the computed local matrices into the global matrix.
      !%OMP CRITICAL
      if (associated(p_Da11)) &
        call incorporateLocalMatrix (DentryA11,p_Da11,rassemblyInfo%Kentry,dweight)

      if (associated(p_Da12)) &
        call incorporateLocalMatrix (DentryA12,p_Da12,rassemblyInfo%Kentry12,dweight)

      if (associated(p_Da21)) &
        call incorporateLocalMatrix (DentryA21,p_Da21,rassemblyInfo%Kentry12,dweight)

      if (associated(p_Da22) .and. .not. associated (p_Da11,p_Da22)) &
        call incorporateLocalMatrix (DentryA22,p_Da22,rassemblyInfo%Kentry,dweight)

    end do ! ielset
    !%OMP end do 
    
    ! Release memory
    deallocate(DentryA11)
    deallocate(DentryA22)
    if (allocated(DentryA12)) deallocate(DentryA12)
    if (allocated(DentryA21)) deallocate(DentryA21)
               
  contains
    
    subroutine incorporateLocalMatrix (DaLocal,Da,Kentry,dweight)
    
    ! Incorporate a local matrix into a global one.
    
    ! The local matrix to incorporate.
    real(dp), dimension(:,:,:), intent(in) :: DaLocal
    
    ! The global matrix
    real(dp), dimension(:), intent(inout) :: Da
    
    ! Positions of the local matrix in the global one.
    integer, dimension(:,:,:), intent(in) :: Kentry
    
    ! Weight of the local matrix
    real(dp) :: dweight

      ! local variables
      integer :: iel,idofe,jdofe
    
      do iel=1,ubound(DaLocal,3)
        do idofe=1,ubound(DaLocal,2)
          do jdofe=1,ubound(DaLocal,1)
            Da(Kentry(jdofe,idofe,iel)) = Da(Kentry(jdofe,idofe,iel)) + &
                dweight * DaLocal(jdofe,idofe,iel) 
          end do
        end do
      end do
    
    end subroutine
    
  end subroutine
     
  ! ***************************************************************************

!<subroutine>
  subroutine conv_strdiffOptC2dgetDerMatrix (rmatrix,roperator,dweight,&
      rvelocityVector,dh,rvelocityVectorAlt,rcollection)
!<description>
  ! Calculate the derivative matrix of the nonlinear operator:
  !   dweight*A(rvelocityVector)
!</description>

!<remarks>
  ! rvelocityVectorAlt Alt allows to specify an
  ! alternative velocity in the following sense:
  ! If NOT specified, the routine assumes that an operator
  !   R(x) = A(x)x
  ! is given and calculates the derivative with respect to x. On the other hand,
  ! if y:= rvelocityVectorAlt is specified,
  ! the routine assumes that the operator is given in the form
  !   R(x) = R(x,y) = A(x)y
  ! and calculates the derivative only with respect to x.
  ! This can be used to calculate the derivative matrix of subblocks in
  ! a block matrix, where the block A(x) changes while the corresponding
  ! evaluation point y stays fixed.
!</remarks>

!<input>

  ! Structure defining the operator to set up.
  type(t_navStOperator), intent(in) :: roperator
  
  ! Weight of the operator. Standard = 1.0
  real(dp), intent(in) :: dweight
      
  ! Velocity vector for the nonlinearity.
  ! The first blocks 1/2 in this vector define the evaluation
  ! point (primal velocity).
  type(t_vectorBlock), intent(in) :: rvelocityVector

  ! Step length for the discrete derivative.
  real(dp), intent(in) :: dh

  ! Alternative Velocity vector for the nonlinearity to be multiplied to the
  ! matrix. The first blocks 1/2 in this vector define the evaluation
  ! point (velocity).
  ! If not present, this coincides with rvelocityVector
  type(t_vectorBlock), intent(in), optional :: rvelocityVectorAlt

  ! Collection structure with application settings.
  type(t_collection), intent(inout), optional :: rcollection
!</input>

!<inputoutput>
  ! The system matrix. The submatrices for the velocity must be in block
  ! A11, A12, A21 and A22 and must be in matrix format 7 or 9.
  ! A11 and A22 must have the same structure. A12 and A21 must have
  ! the same structure.
  type(t_matrixBlock), intent(inout), target :: rmatrix
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: ielset,ielmax,idof,iel
  real(dp) :: dweight2
  
  ! Matrix structure arrays
  real(DP), dimension(:), pointer :: p_Da11,p_Da22
  real(DP), dimension(:), pointer :: p_Da12,p_Da21
  
  ! Local matrices
  real(DP), dimension(:,:,:), allocatable :: DentryA11
  real(DP), dimension(:,:,:), allocatable :: DentryA12
  real(DP), dimension(:,:,:), allocatable :: DentryA21
  real(DP), dimension(:,:,:), allocatable :: DentryA22
                                                      
  real(DP), dimension(:), allocatable :: Dtemp
  
  ! Temp vector
  type(t_vectorBlock) :: rvectorBlock
  
  ! Assembly status structure
  type(t_navStAssemblyInfo) :: rassemblyInfo

    if (roperator%dnu .eq. 0.0_DP) then
      print *,'SD: NU=0 not allowed! Set dbeta=0 to prevent Stokes operator'// &
              ' from being build!'
      call sys_halt()
    end if
    
    ! Get pointers to the matrix content 
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1),p_Da11)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,2),p_Da22)

    nullify(p_Da12,p_Da21)
    if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,2),p_Da12)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,1),p_Da21)
    end if

    ! Initialise the asembly of the local matrices.
    call conv_strdiffOptC2dinitasm (rvelocityVector%p_rblockDiscr%RspatialDiscr(1),&
        1,roperator,rassemblyInfo,present(rvelocityVectorAlt))

    ! Calculate the maximum norm of the actual velocity field
    ! Round up the norm to 1D-8 if it's too small...
    call lsysbl_deriveSubvector(rvelocityVector,rvectorBlock,1,2,.true.)
    call lsysbl_getVectorMagnitude (rvectorBlock,dumax=rassemblyInfo%dumax)
    call lsysbl_releaseVector (rvectorBlock)
  
    if (rassemblyInfo%dumax .lt. 1E-8_DP) rassemblyInfo%dumax=1E-8_DP

    ! Allocate memory for the local matrices
    allocate(DentryA11(rassemblyInfo%indof,rassemblyInfo%indof,&
        rassemblyInfo%nelementsPerBlock))
    allocate(DentryA12(rassemblyInfo%indof,rassemblyInfo%indof,&
        rassemblyInfo%nelementsPerBlock))
    allocate(DentryA21(rassemblyInfo%indof,rassemblyInfo%indof,&
        rassemblyInfo%nelementsPerBlock))
    allocate(DentryA22(rassemblyInfo%indof,rassemblyInfo%indof,&
        rassemblyInfo%nelementsPerBlock))

    ! Temp array
    allocate(Dtemp(rassemblyInfo%indof))

    ! Initialise the array with the local Delta values for the stabilisation
    call lalg_clearVectorDble (rassemblyInfo%DlocalDelta)
    call lalg_clearVectorDble (rassemblyInfo%DlocalDeltaDual)

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so BILF_NELEMSIM local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !%OMP do SCHEDULE(dynamic,1)
    do ielset = 1, size(rassemblyInfo%p_IelementList), BILF_NELEMSIM

      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      ielmax = min(size(rassemblyInfo%p_IelementList),ielset-1+BILF_NELEMSIM)
    
      ! To assemble the discrete Newton matrix, we have to
      ! * Assemble the matrix when modifying each DOF on an element
      ! * Multiplication with the corresponding solution vector gives one
      !   column of the Newton matrix
      !
      ! The global matrix structure of the main matrix looks like:
      !
      !   (A11     ) 
      !   (    A22 ) 
      !
      ! The corresponding Newton matrix has the shape
      !
      !   (A11 A12)
      !   (A21 A22)
      !
      ! We calculate the matrix column-wise and element based:
      ! Let Bij := 1/2h * (Aij(u+h*e_k)) and
      !     Cij := 1/2h * (Aij(u-h*e_k)), then we have with
      ! (pdofp1,pdofp2,ddofp1,ddofp2) = u+h*e_k
      ! (pdofn1,pdofn2,ddofn1,ddofn2) = u-h*e_k
      !
      !   (B11    ) (pdofp1) - (C11    ) (pdofn1)  ->  (    X          )
      !   (    B22) (pdofp2)   (    C22) (pdofn2)      (    X          )
      !                                                     ^k
    
      ! Loop through the DOF's. All DOF's must be once increased and once decreased
      ! by h.
      dweight2 = dweight*0.5_DP/dh
      do idof = 1,rassemblyInfo%indof
    
        ! Initialise the element set, compute the basis functions in the
        ! cubature points.
        if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,2),rvelocityVector,&
              roperator,rassemblyInfo,ielset,ielmax,&
              rvelocityVectorAlt)
        else
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,1),rvelocityVector,&
              roperator,rassemblyInfo,ielset,ielmax,&
              rvelocityVectorAlt)
        end if
        
        ! Increase the idof'th entry in the velocity evaluation point by h.
        do iel = 1,ielmax
          rassemblyInfo%DvelDofs(idof,iel,1) = rassemblyInfo%DvelDofs(idof,iel,1) + dh
        end do

        ! Clear the local matrices. If the Newton part is to be calculated,
        ! we must clear everything, otherwise only Dentry.
        DentryA11 = 0.0_DP
        DentryA12 = 0.0_DP
        DentryA21 = 0.0_DP
        DentryA22 = 0.0_DP
        
        ! Now calculate the local matrices on all the elements.
        call computeLocalNavStMatrices (rassemblyInfo%Dbas,&
            rassemblyInfo%Domega,rassemblyInfo%revalElementSet%p_Ddetj,&
            rassemblyInfo%indof,rassemblyInfo%ncubp,&
            rassemblyInfo%p_IelementList(ielset:ielmax),&
            rassemblyInfo%DvelDofs,rassemblyInfo%Dvel,rassemblyInfo%DvelXderiv,&
            rassemblyInfo%DvelYderiv,rassemblyInfo%Dnu,&
            DentryA11,DentryA22,DentryA12,DentryA21,&
            rassemblyInfo%DlocalDelta,roperator,rassemblyInfo,rcollection)
        
        ! Incorporate the computed local matrices into the global matrix.
        if (associated(p_Da11)) then
          call incorporateLocalMatVecCol (DentryA11,p_Da11,rassemblyInfo%Kentry,&
                                          rassemblyInfo%DvelDofsAlt,1,idof,dweight2,Dtemp)
        end if

        if (associated(p_Da21)) then
          call incorporateLocalMatVecCol (DentryA22,p_Da21,rassemblyInfo%Kentry,&
                                          rassemblyInfo%DvelDofsAlt,2,idof,dweight2,Dtemp)
        end if

      end do
      
      dweight2 = -dweight*0.5_DP/dh
      do idof = 1,rassemblyInfo%indof
    
        ! Initialise the element set, compute the basis functions in the
        ! cubature points.
        if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,2),rvelocityVector,&
              roperator,rassemblyInfo,ielset,ielmax,&
              rvelocityVectorAlt)
        else
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,1),rvelocityVector,&
              roperator,rassemblyInfo,ielset,ielmax,&
              rvelocityVectorAlt)
        end if
        
        ! Increase the idof'th entry in the velocity evaluation point by h.
        do iel = 1,ielmax
          rassemblyInfo%DvelDofs(idof,iel,1) = rassemblyInfo%DvelDofs(idof,iel,1) - dh
        end do

        ! Clear the local matrices. If the Newton part is to be calculated,
        ! we must clear everything, otherwise only Dentry.
        DentryA11 = 0.0_DP
        DentryA12 = 0.0_DP
        DentryA21 = 0.0_DP
        DentryA22 = 0.0_DP
        
        ! Now calculate the local matrices on all the elements.
        call computeLocalNavStMatrices (rassemblyInfo%Dbas,&
            rassemblyInfo%Domega,rassemblyInfo%revalElementSet%p_Ddetj,&
            rassemblyInfo%indof,rassemblyInfo%ncubp,&
            rassemblyInfo%p_IelementList(ielset:ielmax),&
            rassemblyInfo%DvelDofs,rassemblyInfo%Dvel,rassemblyInfo%DvelXderiv,&
            rassemblyInfo%DvelYderiv,rassemblyInfo%Dnu,&
            DentryA11,DentryA22,DentryA12,DentryA21,&
            rassemblyInfo%DlocalDelta,roperator,rassemblyInfo,rcollection)
        
        ! Incorporate the computed local matrices into the global matrix.
        if (associated(p_Da11)) then
          call incorporateLocalMatVecCol (DentryA11,p_Da11,rassemblyInfo%Kentry,&
                                          rassemblyInfo%DvelDofsAlt,1,idof,dweight2,Dtemp)
        end if

        if (associated(p_Da21)) then
          call incorporateLocalMatVecCol (DentryA22,p_Da21,rassemblyInfo%Kentry,&
                                          rassemblyInfo%DvelDofsAlt,2,idof,dweight2,Dtemp)
        end if

      end do
      
      ! Loop through the DOF's. All DOF's must be once increased and once decreased
      ! by h.
      dweight2 = dweight*0.5_DP/dh
      do idof = 1,rassemblyInfo%indof
    
        ! Initialise the element set, compute the basis functions in the
        ! cubature points.
        if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,2),rvelocityVector,&
              roperator,rassemblyInfo,ielset,ielmax,&
              rvelocityVectorAlt)
        else
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,1),rvelocityVector,&
              roperator,rassemblyInfo,ielset,ielmax,&
              rvelocityVectorAlt)
        end if
        
        ! Increase the idof'th entry in the velocity evaluation point by h.
        do iel = 1,ielmax
          rassemblyInfo%DvelDofs(idof,iel,2) = rassemblyInfo%DvelDofs(idof,iel,2) + dh
        end do

        ! Clear the local matrices. If the Newton part is to be calculated,
        ! we must clear everything, otherwise only Dentry.
        DentryA11 = 0.0_DP
        DentryA12 = 0.0_DP
        DentryA21 = 0.0_DP
        DentryA22 = 0.0_DP
        
        ! Now calculate the local matrices on all the elements.
        call computeLocalNavStMatrices (rassemblyInfo%Dbas,&
            rassemblyInfo%Domega,rassemblyInfo%revalElementSet%p_Ddetj,&
            rassemblyInfo%indof,rassemblyInfo%ncubp,&
            rassemblyInfo%p_IelementList(ielset:ielmax),&
            rassemblyInfo%DvelDofs,rassemblyInfo%Dvel,rassemblyInfo%DvelXderiv,&
            rassemblyInfo%DvelYderiv,rassemblyInfo%Dnu,&
            DentryA11,DentryA22,DentryA12,DentryA21,&
            rassemblyInfo%DlocalDelta,roperator,rassemblyInfo,rcollection)
        
        ! Incorporate the computed local matrices into the global matrix.
        if (associated(p_Da12)) then
          call incorporateLocalMatVecCol (DentryA11,p_Da12,rassemblyInfo%Kentry,&
                                          rassemblyInfo%DvelDofsAlt,1,idof,dweight2,Dtemp)
        end if

        if (associated(p_Da22)) then
          call incorporateLocalMatVecCol (DentryA22,p_Da22,rassemblyInfo%Kentry,&
                                          rassemblyInfo%DvelDofsAlt,2,idof,dweight2,Dtemp)
        end if

      end do

      dweight2 = -dweight*0.5_DP/dh
      do idof = 1,rassemblyInfo%indof
    
        ! Initialise the element set, compute the basis functions in the
        ! cubature points.
        if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,2),rvelocityVector,&
              roperator,rassemblyInfo,ielset,ielmax,&
              rvelocityVectorAlt)
        else
          call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
              rmatrix%RmatrixBlock(1,1),rvelocityVector,&
              roperator,rassemblyInfo,ielset,ielmax,&
              rvelocityVectorAlt)
        end if
        
        ! Increase the idof'th entry in the velocity evaluation point by h.
        do iel = 1,ielmax
          rassemblyInfo%DvelDofs(idof,iel,2) = rassemblyInfo%DvelDofs(idof,iel,2) - dh
        end do

        ! Clear the local matrices. If the Newton part is to be calculated,
        ! we must clear everything, otherwise only Dentry.
        DentryA11 = 0.0_DP
        DentryA12 = 0.0_DP
        DentryA21 = 0.0_DP
        DentryA22 = 0.0_DP
        
        ! Now calculate the local matrices on all the elements.
        call computeLocalNavStMatrices (rassemblyInfo%Dbas,&
            rassemblyInfo%Domega,rassemblyInfo%revalElementSet%p_Ddetj,&
            rassemblyInfo%indof,rassemblyInfo%ncubp,&
            rassemblyInfo%p_IelementList(ielset:ielmax),&
            rassemblyInfo%DvelDofs,rassemblyInfo%Dvel,rassemblyInfo%DvelXderiv,&
            rassemblyInfo%DvelYderiv,rassemblyInfo%Dnu,&
            DentryA11,DentryA22,DentryA12,DentryA21,&
            rassemblyInfo%DlocalDelta,roperator,rassemblyInfo,rcollection)
        
        ! Incorporate the computed local matrices into the global matrix.
        if (associated(p_Da12)) then
          call incorporateLocalMatVecCol (DentryA11,p_Da12,rassemblyInfo%Kentry,&
                                          rassemblyInfo%DvelDofsAlt,1,idof,dweight2,Dtemp)
        end if

        if (associated(p_Da22)) then
          call incorporateLocalMatVecCol (DentryA22,p_Da22,rassemblyInfo%Kentry,&
                                          rassemblyInfo%DvelDofsAlt,2,idof,dweight2,Dtemp)
        end if
            
      end do
     
    end do ! ielset
    !%OMP end do 
    
    deallocate(Dtemp)
    
    ! Release memory
    deallocate(DentryA11)
    deallocate(DentryA22)
    if (allocated(DentryA12)) deallocate(DentryA12)
    if (allocated(DentryA21)) deallocate(DentryA21)
               
  contains
  
    subroutine incorporateLocalMatVecCol (DaLocal,Da,Kentry,Dvelocity,ivelcomp,idof,dweight,Dtemp,&
        dmin,dmax)
    
    ! Incorporate a local matrix into a global one.
    
    ! The local matrix to incorporate.
    real(dp), dimension(:,:,:), intent(in) :: DaLocal
    
    ! The global matrix
    real(dp), dimension(:), intent(inout) :: Da
    
    ! Positions of the local matrix in the global one.
    integer, dimension(:,:,:), intent(in) :: Kentry
    
    ! Velocity vector on the elements
    real(dp), dimension(:,:,:), intent(inout) :: Dvelocity
    
    ! Current velocity component
    integer, intent(in) :: ivelcomp
    
    ! Current DOF. Specifies the target column where to write data to.
    integer, intent(in) :: idof
    
    ! Weight of the local matrix
    real(dp) :: dweight
    
    ! Temporary vector
    real(dp), dimension(:), intent(inout) :: Dtemp

    ! Min/Max value for the velocity. If present, the velocity will be
    ! projected to this range.
    real(dp), intent(in), optional :: dmin,dmax

      ! local variables
      integer :: iel,idofe,jdofe
    
      if (.not. present(dmin)) then
      
        do iel=1,ubound(DaLocal,3)
          ! Apply a local matrix-vector multiplication to get the column
          ! which is to be incorporated into the matrix.
          Dtemp(:) = 0.0_DP
          do idofe = 1,ubound(DaLocal,2)
            do jdofe = 1,ubound(DaLocal,1)
              Dtemp(idofe) = Dtemp(idofe) + DaLocal(jdofe,idofe,iel)*Dvelocity(jdofe,iel,ivelcomp)
            end do
          end do
        
          ! Incorporate the vector to the global matrix.
          do idofe=1,ubound(DaLocal,2)
            Da(Kentry(idof,idofe,iel)) = Da(Kentry(idof,idofe,iel)) + &
                dweight * Dtemp(idofe)
          end do
        end do
        
      else
      
        do iel=1,ubound(DaLocal,3)
          ! Apply a local matrix-vector multiplication to get the column
          ! which is to be incorporated into the matrix.
          Dtemp(:) = 0.0_DP
          do idofe = 1,ubound(DaLocal,2)
            do jdofe = 1,ubound(DaLocal,1)
              Dtemp(idofe) = Dtemp(idofe) + &
                  DaLocal(jdofe,idofe,iel)*min(max(Dvelocity(jdofe,iel,ivelcomp),dmin),dmax)
            end do
          end do
        
          ! Incorporate the vector to the global matrix.
          do idofe=1,ubound(DaLocal,2)
            Da(Kentry(idof,idofe,iel)) = Da(Kentry(idof,idofe,iel)) + &
                dweight * Dtemp(idofe)
          end do
        end do
      
      end if
    
    end subroutine
    
  end subroutine
     
  ! ***************************************************************************

!<subroutine>
  subroutine conv_strdiffOptC2dgetDefect (rvelMatrix,roperator,&
      rvelocityVector,dweight,rx,rd,rcollection)
      
!<description>
  ! Calculate the defect of the nonlinear operator:
  !   rd = rd - dweight*A(rvelocityVector)rx 
!</description>

!<input>
  ! Template FE matrix specifying the connectivity in the FE space
  ! of the primal/dual velocity.
  type(t_matrixScalar), intent(in) :: rvelMatrix

  ! Structure defining the operator to set up.
  type(t_navStOperator), intent(in) :: roperator
  
  ! Velocity vector for the nonlinearity.
  ! The first blocks 1/2 in this vector define the evaluation
  ! point (primal velocity).
  type(t_vectorBlock), intent(in) :: rvelocityVector

  ! Solution vector, to be multiplied with the matrix.
  type(t_vectorBlock), intent(in) :: rx
  
  ! Weight for the solution
  real(DP), intent(in) :: dweight

  ! Collection structure with application settings.
  type(t_collection), intent(inout), optional :: rcollection
!</input>
  
!<inputoutput>  
  ! On entry: RHS. On exit: Defect.
  type(t_vectorBlock), intent(inout) :: rd
!</inputoutput>  

!</subroutine>

  ! local variables
  integer :: ielset,ielmax
  
  ! Local matrices
  real(DP), dimension(:,:,:), allocatable :: DentryA11
  real(DP), dimension(:,:,:), allocatable :: DentryA12
  real(DP), dimension(:,:,:), allocatable :: DentryA21
  real(DP), dimension(:,:,:), allocatable :: DentryA22
  
  ! Temp vector
  type(t_vectorBlock) :: rvectorBlock
  
  ! Pointer to the subvectors
  real(dp), dimension(:), pointer :: p_Dx1,p_Dx2,p_Dd1,p_Dd2
  
  ! Assembly status structure
  type(t_navStAssemblyInfo) :: rassemblyInfo

    if (roperator%dnu .eq. 0.0_DP) then
      print *,'SD: NU=0 not allowed! Set dbeta=0 to prevent Stokes operator'// &
              ' from being build!'
      call sys_halt()
    end if
    
    ! Initialise the asembly of the local matrices.
    call conv_strdiffOptC2dinitasm (rvelocityVector%p_rblockDiscr%RspatialDiscr(1),&
        1,roperator,rassemblyInfo)

    ! Calculate the maximum norm of the actual velocity field
    ! Round up the norm to 1D-8 if it's too small...
    call lsysbl_deriveSubvector(rvelocityVector,rvectorBlock,1,2,.true.)
    call lsysbl_getVectorMagnitude (rvectorBlock,dumax=rassemblyInfo%dumax)
    call lsysbl_releaseVector (rvectorBlock)
  
    if (rassemblyInfo%dumax .lt. 1E-8_DP) rassemblyInfo%dumax=1E-8_DP

    ! Allocate memory for the local matrices
    allocate(DentryA11(rassemblyInfo%indof,rassemblyInfo%indof,rassemblyInfo%nelementsPerBlock))
    allocate(DentryA12(rassemblyInfo%indof,rassemblyInfo%indof,rassemblyInfo%nelementsPerBlock))
    allocate(DentryA21(rassemblyInfo%indof,rassemblyInfo%indof,rassemblyInfo%nelementsPerBlock))
    allocate(DentryA22(rassemblyInfo%indof,rassemblyInfo%indof,rassemblyInfo%nelementsPerBlock))

    ! Initialise the array with the local Delta values for the stabilisation
    call lalg_clearVectorDble (rassemblyInfo%DlocalDelta)
    
    ! Get pointers to the subvectors
    call lsyssc_getbase_double (rx%RvectorBlock(1),p_Dx1)
    call lsyssc_getbase_double (rx%RvectorBlock(2),p_Dx2)
    call lsyssc_getbase_double (rd%RvectorBlock(1),p_Dd1)
    call lsyssc_getbase_double (rd%RvectorBlock(2),p_Dd2)

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so BILF_NELEMSIM local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !%OMP do SCHEDULE(dynamic,1)
    do ielset = 1, size(rassemblyInfo%p_IelementList), BILF_NELEMSIM

      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      ielmax = min(size(rassemblyInfo%p_IelementList),ielset-1+BILF_NELEMSIM)
    
      ! Initialise the element set, compute the basis functions in the
      ! cubature points.
      call conv_strdiffOptC2dinitelemset (rvelMatrix,rvelMatrix,&
          rvelocityVector,roperator,&
          rassemblyInfo,ielset,ielmax)

      ! Clear the local matrices. If the Newton part is to be calculated,
      ! we must clear everything, otherwise only Dentry.
      DentryA11 = 0.0_DP
      DentryA12 = 0.0_DP
      DentryA21 = 0.0_DP
      DentryA22 = 0.0_DP
      
      ! Now calculate the local matrices on all the elements.
      call computeLocalNavStMatrices (rassemblyInfo%Dbas,&
          rassemblyInfo%Domega,rassemblyInfo%revalElementSet%p_Ddetj,&
          rassemblyInfo%indof,rassemblyInfo%ncubp,&
          rassemblyInfo%p_IelementList(ielset:ielmax),&
          rassemblyInfo%DvelDofs,rassemblyInfo%Dvel,rassemblyInfo%DvelXderiv,&
          rassemblyInfo%DvelYderiv,rassemblyInfo%Dnu,&
          DentryA11,DentryA22,DentryA12,DentryA21,&
          rassemblyInfo%DlocalDelta,roperator,rassemblyInfo,rcollection)
      
      ! Perform a local matrix-vector multiplication with the local matrices.
      ! Incorporate the computed local matrices into the global matrix.
      call calcDefect (DentryA11,rassemblyInfo%Idofs,p_Dx1,p_Dd1,dweight)
      call calcDefect (DentryA12,rassemblyInfo%Idofs,p_Dx2,p_Dd1,dweight)
      call calcDefect (DentryA21,rassemblyInfo%Idofs,p_Dx1,p_Dd2,dweight)
      call calcDefect (DentryA22,rassemblyInfo%Idofs,p_Dx2,p_Dd2,dweight)
      
    end do ! ielset
    !%OMP end do 
    
    ! Release memory
    deallocate(DentryA11)
    deallocate(DentryA22)
    deallocate(DentryA12)
    deallocate(DentryA21)
                              
  contains
    
    subroutine calcDefect (DaLocal,KcolLocal,Dx,Dd,dweight,dmin,dmax)
    
    ! Calculate the local defect:
    !   Dd = Dd - DaLocal*Dx
    
    ! The local matrix to incorporate.
    real(dp), dimension(:,:,:), intent(in) :: DaLocal
    
    ! Global column numbers of the local matrix columns on each element.
    integer, dimension(:,:), intent(in) :: KcolLocal
    
    ! Solution vector
    real(dp), dimension(:), intent(in) :: Dx
    
    ! RHS vector; receives the local defect.
    real(dp), dimension(:), intent(inout) :: Dd
    
    ! Weight for the solution vectot
    real(dp), intent(in) :: dweight
    
    ! Min/Max value for the velocity. If present, the velocity will be
    ! projected to this range.
    real(dp), intent(in), optional :: dmin,dmax
    
      ! local variables
      integer :: iel,idofe,jdofe
      real(dp) :: dval
    
      if (.not. present(dmin)) then
        do iel=1,ubound(DaLocal,3)
          do idofe=1,ubound(DaLocal,1)
            dval = 0.0_DP
            ! DaLocal*Dx. Note that DaLocal is saved transposed for quicker access,
            ! so jdofe in the first index describes the column!
            do jdofe=1,ubound(DaLocal,2)
              dval = dval + DaLocal(jdofe,idofe,iel) * Dx(KcolLocal(jdofe,iel))
            end do
            Dd(KcolLocal(idofe,iel)) = Dd(KcolLocal(idofe,iel)) - dweight*dval
          end do
        end do
      else
        do iel=1,ubound(DaLocal,3)
          do idofe=1,ubound(DaLocal,1)
            dval = 0.0_DP
            ! DaLocal*Dx. Note that DaLocal is saved transposed for quicker access,
            ! so jdofe in the first index describes the column!
            do jdofe=1,ubound(DaLocal,2)
              dval = dval + DaLocal(jdofe,idofe,iel) * &
                  min(max(Dx(KcolLocal(jdofe,iel)),dmin),dmax)
            end do
            Dd(KcolLocal(idofe,iel)) = Dd(KcolLocal(idofe,iel)) - dweight*dval
          end do
        end do
      end if
    
    end subroutine
    
  end subroutine
     
  ! ----------------------------------------------------------------------

  subroutine getLocalDeltaQuad (Du,rtriangulation,Ielements,&
      dumax,dupsam,Dnu,DlocalDelta)

  ! This routine calculates a local ddelta=DELTA_T for a finite element
  ! T=IEL. This can be used by the streamline diffusion stabilisation
  ! technique as a multiplier of the (local) bilinear form.
  
  ! Velocity DOF's on all elements.
  ! Du(#dofs per element,#elements,ndim2d)
  real(DP), dimension(:,:,:), intent(in) :: Du
  
  ! Underlying triangulation
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! List of elements to process
  integer, dimension(:), intent(in) :: Ielements
  
  ! Maximum norm of velocity in the domain:
  ! duMaxR = ||u||_Omega
  real(DP), intent(IN) :: dumax
  
  ! Viscosity parameters in the cubature points
  real(DP), dimension(:,:), intent(in) :: Dnu
  
  ! user defined parameter for configuring the streamline diffusion.
  ! < 0: Simple calculation of ddelta, using 
  !      ddelta = |UPSAM| * h_T.
  ! > 0: usually UPSAM = 0.1 .. 2; Samarskji-like calculation of ddelta using:
  !      ddelta = UPSAM * h_t/||u||_T * 2*Re_T/(1+Re_T)
  real(DP), intent(IN) :: dupsam
  
  ! local local delta on all elements in the list
  real(DP), dimension(:), intent(out) :: DlocalDelta

    ! local variables
    real(DP) :: dlocalH,du1,du2,dunorm,dreLocal,dumaxR,dlocalnu
    integer :: iel,i
    integer :: idof
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    call storage_getbase_int2d(rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    call storage_getbase_double2d(rtriangulation%h_DvertexCoords,p_DvertexCoords)

    ! We later need the reciprocals of dumax
    dumaxR = 1.0_DP/dumax

    ! Loop over all elements
    do iel = 1,size(Ielements)
    
      ! Calculate the mean velocity on the current element
      ! and using this the norm of the local velocity.
      du1=0.0_dp
      du2=0.0_dp
      do idof=1,ubound(Du,1)
        du1=du1+Du(idof,iel,1)
        du2=du2+Du(idof,iel,2)
      end do
    
      dunorm = sqrt(du1**2+du2**2) / dble(ubound(Du,1))
    
      ! If the norm of the velocity is small, we choose ddelta = 0,
      ! which results in central difference in the streamline diffusion
      ! matrix assembling:

      if (dunorm .le. 1D-8) then
      
        DlocalDelta(iel) = 0.0_DP

      else

        ! u_T defines the "slope" of the velocity through
        ! the element T. At next, calculate the local mesh width
        ! dlocalH = h = h_T on our element T=IEL:

        call getLocalMeshWidthQuad (dlocalH,dunorm, du1, du2, iel, &
            p_IverticesAtElement,p_DvertexCoords)

        ! Calculate ddelta... (cf. p. 121 in Turek's CFD book)

        if (dupsam .lt. 0.0_DP) then

          ! For UPSAM<0, we use simple calculation of ddelta:        
        
          DlocalDelta(iel) = abs(dupsam)*dlocalH
          
        else
        
          ! For UPSAM >= 0, we use standard Samarskji-like calculation
          ! of ddelta. At first calculate the local Reynolds number
          ! RELOC = Re_T = ||u||_T * h_T / NU
          !
          ! To calculate the local nu, take the average nu over all 
          ! cubature points.
          
          dlocalnu = 0
          do i=1,ubound(Dnu,1)
            dlocalnu = dlocalnu+Dnu(i,iel)
          end do
          
          dreLocal = dunorm * dlocalH * real(ubound(Dnu,1),dp)/dlocalnu
          
          ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)
          
          DlocalDelta(iel) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLocal/(1.0_DP+dreLocal))
          
        end if ! (UPSAM.LT.0.0)
        
      end if ! (dunorm.LE.1D-8)

    end do

  end subroutine

  ! ----------------------------------------------------------------------

  pure subroutine getLocalMeshWidthQuad (dlocalH, dunorm,  XBETA1, &
                      XBETA2, JEL,Kvert,Dcorvg)
  
  ! Determine the local mesh width for an element JEL of a 
  ! triangulation.
  
  ! Element where the local h should be calculated
  integer, intent(IN)               :: JEL
  
  integer, dimension(TRIA_MAXNVE2D,*), intent(IN) :: Kvert
  real(DP), dimension(NDIM2D,*), intent(IN)          :: Dcorvg
  
  ! norm ||u||_T = mean velocity through element T=JEL
  real(DP), intent(IN)  :: dunorm
  
  ! mean velocity u_T = (xbeta1,xbeta2) through element T=JEL
  real(DP), intent(IN)  :: XBETA1, XBETA2
  
  ! local mesh width
  real(DP), intent(OUT) :: dlocalH
  
  ! local variables
  real(DP) :: dlambda
  integer :: NECK1,NECK2,NECK3,NECK4
  real(DP) :: X1,Y1,X2,Y2,X3,Y3,X4,Y4
  real(DP) :: dalphaMax, dalpha

    ! Fetch the numbers of the four corners of element JEL

    neck1=Kvert(1,JEL)
    neck2=Kvert(2,JEL)
    neck3=Kvert(3,JEL)
    neck4=Kvert(4,JEL)

    ! Fetch the coordinates of these corners

    x1=Dcorvg(1,neck1)
    y1=Dcorvg(2,neck1)
    x2=Dcorvg(1,neck2)
    y2=Dcorvg(2,neck2)
    x3=Dcorvg(1,neck3)
    y3=Dcorvg(2,neck3)
    x4=Dcorvg(1,neck4)
    y4=Dcorvg(2,neck4)

    ! Scale: (deactivated)

    !  dsp=max(xbeta1,xbeta2)

    !  xbeta1=xbeta1
    !  xbeta2=xbeta2

    dalphaMax=0.0_DP
    
    ! In the next step, we calculate the 'maximum possible mesh with
    ! in direction of the flow'; this is the maximum possible length
    ! that a particle can cross in the current element.
    ! The picture in mind is the following:
    !
    !          G3
    !   +-------------X-------+
    !   |            /        |
    !   |           /         |
    !   |          /          |
    !   |         /           |
    !   |        /            |
    ! G4|       /             | G2
    !   |      ^ (beta1,beta2)|
    !   |     /               |
    !   |    /                |
    !   |   /                 |
    !   |  /                  |
    !   | /                   |
    !   |/                    |
    !   O---------------------+
    !            G1
    !
    ! The vector (beta1,beta2) gives the direction of the flow.
    ! A particle starting in point O and moves at most up to point X.
    ! The length of the line (O,X) is the local mesh with h.
    !
    ! Loop through the four corners of element JEL and check
    ! of a line with slope BETA=(xbeta1,xbeta2) starting in this
    ! corner really intersects with one of the edges of the element.
    ! Remark that we only have to check the two opposite edges
    ! to the current corner!

    ! -----------------------------------------------------------------
    ! Check the first corner:

    call intersectLines2D(X1,Y1,dalpha,XBETA1,XBETA2, &
                X3,Y3,dlambda,X2,Y2)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X1,Y1,dalpha,XBETA1,XBETA2, &
                X3,Y3,dlambda,X4,Y4)
    dalphaMax=max(dalpha,dalphaMax)
    
    ! -----------------------------------------------------------------
    ! The second one...
    
    call intersectLines2D(X2,Y2,dalpha,XBETA1,XBETA2, &
                X4,Y4,dlambda,X1,Y1)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X2,Y2,dalpha,XBETA1,XBETA2, &
                X4,Y4,dlambda,X3,Y3)
    dalphaMax=max(dalpha,dalphaMax)
    
    ! -----------------------------------------------------------------
    ! The third one...
    
    call intersectLines2D(X3,Y3,dalpha,XBETA1,XBETA2, &
                X1,Y1,dlambda,X2,Y2)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X3,Y3,dalpha,XBETA1,XBETA2, &
                X1,Y1,dlambda,X4,Y4)
    dalphaMax=max(dalpha,dalphaMax)
    
    ! -----------------------------------------------------------------
    ! And the fourth=last one...
    
    call intersectLines2D(X4,Y4,dalpha,XBETA1,XBETA2, &
                X2,Y2,dlambda,X1,Y1)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X4,Y4,dalpha,XBETA1,XBETA2, &
                X2,Y2,dlambda,X3,Y3)
    dalphaMax=max(dalpha,dalphaMax)

    ! -----------------------------------------------------------------
    ! finally determine the local h=h_T
    !
    ! dalphaMax is the maximum alpha, normalised as 'parameter value',
    ! i.e. dalphaMax=1.0 corresponds to a vector 1.0*(dbeta1,dbeta2).
    ! We multiply with dunorm=|(dbeta1,dbeta2)| to get the actual length
    ! of the vector which can be placed inside of the element.
    !
    ! Furthermore, we multiply with an additional weight 4. (why ?!?)

    dlocalH=dalphaMax*4.0_DP*dunorm

  end subroutine
  
  ! ----------------------------------------------------------------------

  pure subroutine intersectLines2D (XO,YO,dalpha,BETA1,BETA2, &
                      XA,YA,dlambda,XB,YB)

  ! Intersect two lines in R^2

  ! Origin of line 1
  real(DP), intent(IN) :: XO,YO
  
  ! Direction of line 1
  real(DP), intent(IN) :: BETA1,BETA2
  
  ! One point on the second line
  real(DP), intent(IN) :: XA,YA
  
  ! Another point on the second line
  real(DP), intent(IN) :: XB,YB
  
  ! Parameter value of the intersection point on line 1.
  ! =0.0, if there is no intersection point
  real(DP), intent(OUT) :: dalpha
  
  real(DP), intent(OUT) :: dlambda
  
  ! local variables
  double precision :: dsp

    ! Scalar product of the line (xa,ya)->(xb,yb) with the
    ! counterclockwise normal n1 of (beta1,beta2)
    dsp=BETA2*(XB-XA)-BETA1*(YB-YA)
    
    if (dsp.eq.0.0_DP) then
    
      ! beta and the vector are parallel
      dalpha=0.0_DP
      
    else  

      ! Scalar product of (beta1,beta2) with the (inner) normal vector n2
      ! of the line (xo,yo)->(xa,ya). 
      dlambda=(BETA1*(YA-YO)-BETA2*(XA-XO))/dsp

      !                    (xb,yb)
      !   +-----------------+
      !   |                 |
      !   |                 |
      !   ^ n2              |
      !   !                 |
      !   !  (beta1,beta2)  |    (beta1,beta2)
      !   !    ^            |    ^
      !   !   /  ^__ n1     |   /
      !   !  /      \__     |  /
      !   ! /          \__  | /
      !   !/              \_|/ 
      !   +-----------------+
      ! (xo,yo)            (xa,ya)
      !
      ! (What is this? Documentation incomplete. Has someone a good
      ! reference?)

      ! is the intersection point inside of the element?
      if ((dlambda.ge.-1E-1_DP).and.(dlambda.le.1.11E0_DP)) then
        if (BETA1 .ne. 0.0_DP) then
          dalpha=((XA-XO)+dlambda*(XB-XA))/BETA1
        else
          if (BETA2 .ne. 0.0_DP) then
            dalpha=((YA-YO)+dlambda*(YB-YA))/BETA2
          else
            dalpha=0.0_DP
          end if
        end if
      else
        dalpha=0.0_DP
      end if
      
    end if

  end subroutine

end module
