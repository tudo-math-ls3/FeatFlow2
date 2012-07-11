!##############################################################################
!# ****************************************************************************
!# <name> blockmatassembly </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises an assembly routine for block matrices and block 
!# vectors. The assembly involves a callback routine which calculates the 
!# entries in local element matrices/vectors. These matrices are in a second
!# step incorporated into a global matrix / block vector.
!#
!# To support nonlinear matrices and nonconstant coefficient in vectors, 
!# a list of finite element vectors can be provided to the evaluation routine.
!# These vectors are automatically evaluated in all cubature points. 
!# That way, the callback routine for the evaluation of the local matrices 
!# can directly access the data of FEM vectors in all cubature points up 
!# to a specified derivative.
!#
!# The following routines can be found here:
!#
!# 1.) bma_buildMatrix
!#     -> Assembles a block matrix
!#
!# 2.) bma_buildVector
!#     -> Assembles a block vector
!#
!# 2.) bma_buildIntegral
!#     -> Calculates the value of an integral
!#
!# Furthermore, there are a set of auxiliary routines which can be used
!# to customize the matrix assembly:
!#
!# 1.) bma_initMatAssembly
!#     -> Initialise the matrix assembly
!#
!# 2.) bma_assembleSubmeshMatrix
!#     -> Assemble a part of the matrix based on a set of cells.
!#        All cells must have the same element type, the same cubature formula
!#        and the same transformation
!#
!# 3.) bma_doneMatAssembly
!#     -> Clean up the matrix assembly.
!#
!# 4.) bma_initVecAssembly
!#     -> Initialise the vector assembly
!#
!# 5.) bma_assembleSubmeshVector
!#     -> Assemble a part of a vector based on a set of cells.
!#        All cells must have the same element type, the same cubature formula
!#        and the same transformation
!#
!# 6.) bma_doneVecAssembly
!#     -> Clean up the vector assembly.
!#
!# 7.) bma_initIntAssembly
!#     -> Initialise the integral assembly
!#
!# 8.) bma_assembleSubmeshIntegral
!#     -> Assemble a part of an integral based on a set of cells.
!#        All cells must have the same element type, the same cubature formula
!#        and the same transformation
!#
!# 9.) bma_doneIntAssembly
!#     -> Clean up the integral assembly.
!#
!# The following set of callback routines realise standard operators.
!# They can directly be used with bma_buildMatrix, bma_buildVector and 
!# bma_buildIntegral:
!#
!# 1.) bma_fcalc_mass
!#     -> Used with bma_buildMatrix, this calculates mass matrices in all
!#        diagonal blocks of a block matrix.
!#
!# 2.) bma_fcalc_laplace
!#     -> Used with bma_buildMatrix, this calculates Laplace matrices matrices
!#        in all diagonal blocks of a block matrix.
!#
!# 3.) bma_fcalc_convection
!#     -> Used with bma_buildMatrix, this routine can be used to assemble
!#        a standard convection operator in the first diagonal blocks
!#        of a block matrix
!#
!# 4.) bma_fcalc_rhsOne
!#     -> Calculates the RHS vector based on the function f=1.
!#
!# 5.) bma_fcalc_rhsBubble
!#     -> Calculates the RHS vector based on the function f=32*y*(1-y)+32*x*(1-x)
!#        which is the RHS for u=16*x*(1-x)*y*(1-y) in the Laplace
!#        equartion -Laplace(u)=f.
!#
!# 6.) bma_fcalc_rhsBubblePlusFE
!#     -> Calculates the RHS vector based on the function 
!#        f=32*y*(1-y)+32*x*(1-x) + v(x,y)
!#        with v(x,y) being a finite element function passed via parameters.
!#
!# 7.) bma_fcalc_integralOne
!#     -> Calculates the integral of the function v=1 (which results in the
!#        size of the domain).
!#
!# 8.) bma_fcalc_integralFE
!#     -> Calculates the integral of an arbitrary FEM function.
!#
!# 9.) bma_fcalc_bubbleL2error
!#     -> Calculates the squared L2 error of a FEM function to a 
!#        bubble function
!#        
!# 10.) bma_fcalc_bubbleH1error
!#     -> Calculates the squared H1 error of a FEM function to a 
!#        bubble function
!#
!# 11.) bma_fcalc_L2norm
!#     -> Calculates the squared L2 norm of a FEM function
!#        
!# 12.) bma_fcalc_H1norm
!#     -> Calculates the squared H1 (semi-)norm of a FEM function
!# </purpose>
!##############################################################################

module blockmatassembly

  use fsystem
  use storage
  use collection, only: t_collection
  use basicgeometry
  use boundaryaux
  use cubature
  use element
  use elementpreprocessing
  use domainintegration
  use genoutput
  use linearalgebra
  use scalarpde
  use spatialdiscretisation
  use triangulation
  use perfconfig

  use linearsystemscalar
  use linearsystemblock

  use feevaluation2
  use blockmatassemblybase

  implicit none

  private

  !************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: bma_perfconfig

  !************************************************************************

  public :: bma_fcalc_mass
  public :: bma_fcalc_laplace
  public :: bma_fcalc_convection
  
  public :: bma_fcalc_rhsOne
  public :: bma_fcalc_rhsBubble
  public :: bma_fcalc_rhsBubblePlusFE
  
  public :: bma_fcalc_integralOne
  public :: bma_fcalc_integralFE
  public :: bma_fcalc_bubbleL2error
  public :: bma_fcalc_bubbleH1error
  public :: bma_fcalc_L2norm
  public :: bma_fcalc_H1norm

  public :: bma_initMatAssembly
  public :: bma_doneMatAssembly
  public :: bma_assembleSubmeshMatrix
  public :: bma_buildMatrix

  public :: bma_initVecAssembly
  public :: bma_doneVecAssembly
  public :: bma_assembleSubmeshVector
  public :: bma_buildVector

  public :: bma_initIntAssembly
  public :: bma_doneIntAssembly
  public :: bma_assembleSubmeshIntegral
  public :: bma_buildIntegral

contains

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_mass(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the Mass operator in all diagonal matrices.
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! in all diagonal blocks with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    ! rcollection%DquickAccess(1) = multiplier in front of the mass matrix.
    ! rcollection%IquickAccess(1) = First diagonal block or =0 to 
    !                               start from block 1.
    ! rcollection%IquickAccess(2) = Last diagonal block or =0 to calculate 
    !                               up to the last.
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<subroutine>

    ! Local variables
    real(DP) :: dbasI, dbasJ
    integer :: iel, icubp, idofe, jdofe, i, ivar, nvar
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix
    real(DP), dimension(:,:,:,:), pointer :: p_DlocalMatrixIntl
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData

    integer :: istart, iend
    real(DP) :: dscale

    ! First/last block, multiplier
    istart = 1
    iend = min(ubound(RmatrixData,1),ubound(RmatrixData,2))
    dscale = 1.0_DP
    if (present(rcollection)) then
      if (rcollection%IquickAccess(1) .ne. 0) then
        istart = rcollection%IquickAccess(1)
      end if
      if (rcollection%IquickAccess(2) .ne. 0) then
        iend = rcollection%IquickAccess(2)
      end if
      dscale = rcollection%DquickAccess(1)
    end if

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Loop over the diagonal blocks
    do i = istart,iend

      ! Get local data
      p_rmatrixData => RmatrixData(i,i)
      p_DbasTrial => RmatrixData(i,i)%p_DbasTrial
      p_DbasTest => RmatrixData(i,i)%p_DbasTest

      ! Interleaved matrix?
      if (.not. p_rmatrixData%bisInterleaved) then

        ! Get the matrix data      
        p_DlocalMatrix => RmatrixData(i,i)%p_Dentry

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
              do jdofe=1,p_rmatrixData%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)

                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * dbasJ*dbasI

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel

      else

        ! Get the matrix data      
        p_DlocalMatrixIntl => RmatrixData(i,i)%p_DentryIntl

        ! Interleave-data
        nvar = RmatrixData(i,i)%nvar

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
              do jdofe=1,p_rmatrixData%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)

                do ivar = 1,nvar

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                      p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                          dscale * p_DcubWeight(icubp,iel) * dbasJ*dbasI

                end do ! ivar

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel

      end if

    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_laplace(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the Laplace operator in all diagonal matrices
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! in all diagonal blocks with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    ! rcollection%DquickAccess(1) = multiplier in front of the mass matrix.
    ! rcollection%IquickAccess(1) = First diagonal block or =0 to 
    !                               start from block 1.
    ! rcollection%IquickAccess(2) = Last diagonal block or =0 to calculate 
    !                               up to the last.
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<subroutine>

    real(DP) :: dbasIx, dbasJx, dbasIy, dbasJy, dbasIz, dbasJz
    integer :: iel, icubp, idofe, jdofe, i, ivar, nvar
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix
    real(DP), dimension(:,:,:,:), pointer :: p_DlocalMatrixIntl
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData

    integer :: istart, iend
    real(DP) :: dscale

    ! First/last block, multiplier
    istart = 1
    iend = min(ubound(RmatrixData,1),ubound(RmatrixData,2))
    dscale = 1.0_DP
    if (present(rcollection)) then
      if (rcollection%IquickAccess(1) .ne. 0) then
        istart = rcollection%IquickAccess(1)
      end if
      if (rcollection%IquickAccess(2) .ne. 0) then
        iend = rcollection%IquickAccess(2)
      end if
      dscale = rcollection%DquickAccess(1)
    end if

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Loop through all diagonal blocks   
    do i = istart, iend

      ! Get local data
      p_rmatrixData => RmatrixData(i,i)
      p_DbasTrial => RmatrixData(i,i)%p_DbasTrial
      p_DbasTest => RmatrixData(i,i)%p_DbasTest

      ! Interleaved matrix?
      if (.not. p_rmatrixData%bisInterleaved) then

        ! Get the matrix data
        p_DlocalMatrix => RmatrixData(i,i)%p_Dentry

        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)

          ! 1D Laplace matrix

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Phi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV1D_X,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV1D_X,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * dbasJx*dbasIx

                end do ! idofe

              end do ! jdofe

            end do ! icubp

          end do ! iel

        case (NDIM2D)

          ! 2D Laplace matrix

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Phi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
                  dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * ( dbasJx*dbasIx + dbasJy*dbasIy )

                end do ! idofe

              end do ! jdofe

            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! 3D Laplace matrix

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Phi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV3D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV3D_Y,icubp,iel)
                dbasIz = p_DbasTest(idofe,DER_DERIV3D_Z,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                  dbasJy = p_DbasTrial(jdofe,DER_DERIV3D_Y,icubp,iel)
                  dbasJz = p_DbasTrial(jdofe,DER_DERIV3D_Z,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * &
                              ( dbasJx*dbasIx + dbasJy*dbasIy + dbasJz*dbasIz )

                end do ! idofe

              end do ! jdofe

            end do ! icubp

          end do ! iel

        end select

      else

        ! Get the matrix data      
        p_DlocalMatrixIntl => RmatrixData(i,i)%p_DentryIntl

        ! Interleave-data
        nvar = RmatrixData(i,i)%nvar

        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)

          ! 1D Laplace matrix

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Phi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV1D_X,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV1D_X,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                          dscale * p_DcubWeight(icubp,iel) * dbasJx*dbasIx

                  end do ! ivar

                end do ! idofe

              end do ! jdofe

            end do ! icubp

          end do ! iel

        case (NDIM2D)

          ! 2D Laplace matrix

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Phi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
                  dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                            dscale * p_DcubWeight(icubp,iel) * ( dbasJx*dbasIx + dbasJy*dbasIy )

                  end do ! ivar

                end do ! idofe

              end do ! jdofe

            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! 3D Laplace matrix

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Phi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV3D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV3D_Y,icubp,iel)
                dbasIz = p_DbasTest(idofe,DER_DERIV3D_Z,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                  dbasJy = p_DbasTrial(jdofe,DER_DERIV3D_Y,icubp,iel)
                  dbasJz = p_DbasTrial(jdofe,DER_DERIV3D_Z,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                            dscale * p_DcubWeight(icubp,iel) * &
                                    ( dbasJx*dbasIx + dbasJy*dbasIy + dbasJz*dbasIz )

                  end do ! ivar

                end do ! idofe

              end do ! jdofe

            end do ! icubp

          end do ! iel

        end select

      end if

    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_convection(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a convection operator "(u grad) ." into the top/left position
    ! of a block matrix (position 1/1).
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! in all diagonal blocks with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    !
    ! rcollection%DquickAccess(1) = multiplier in front of the operator.
    ! rcollection%IquickAccess(1) = 0, if the convection is a constant vector field.
    !                                  in this case:
    !                                  1D: rcollection%DquickAccess(2)   = x-velocity
    !                                  2D: rcollection%DquickAccess(2:3) = x/y-velocity
    !                                  3D: rcollection%DquickAccess(2:4) = x/y/z-velocity
    !                             = 1, if the convection is specified by a
    !                                  finite element velocity field. In this case,
    !                                  a finite element velocity field must be specified
    !                                  as parameter revalVectors to the call of 
    !                                  bma_buildMatrix. The first vector must be the
    !                                  X-velocity, the 2nd the Y-velocity and 
    !                                  the third the Z-velocity.
!</description>

!<remarks>
    ! Remark 1: 
    ! Using the routines from feevaluation2, it is possible to specify
    ! a nonconstant velocity field. The corresponding code looks as follows:
    !
    ! <verb>
    !     use feevaluation2
    !
    !     ...
    !
    !     type(t_collection) :: rcollection
    !     type(t_scalarCubatureInfo) :: rcubatureInfo   ! Cubature formula
    !     type(t_matrixBlock) :: rmatrix                ! Matrix to be calculated
    !
    !     type(t_vectorBlock) :: rvelocity        ! The velocity field
    !     type(t_fev2Vectors) :: revalVectors     ! Collection of vectors to evaluate
    !
    !     ! Prepare the cubature formula
    !     call spdiscr_createDefCubStructure (..., rcubatureInfo, CUB_GEN_AUTO)
    !
    !     ...
    !
    !     rcollection%IquickAccess(1) = 1          ! Nonconstant viscosity
    !     rcollection%DquickAccess(1) = 1.0_DP     ! Scaling
    !
    !     ! Add the X-, Y- and Z-velocity to revalVectors
    !     call fev2_addVectorToEvalList(revalVectors,rvelocity%RvectorBlock(1),0)
    !     call fev2_addVectorToEvalList(revalVectors,rvelocity%RvectorBlock(2),0)
    !     call fev2_addVectorToEvalList(revalVectors,rvelocity%RvectorBlock(3),0)
    !
    !     ! Set up the matrix
    !     call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD, bma_fcalc_convection, &
    !         rcollection, revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
    !
    !     ! Release the vector structure
    !     call fev2_releaseVectorList(revalVectors)
    !
    !     ...
    !
    !     ! Release the cubature formula
    !     call spdiscr_releaseCubStructure (rcubatureInfo)
    ! </verb>
    !
    ! Remark 2: 
    ! The routine <verb>fev2_addVectorToEvalList</verb> allows to define
    ! the evaluation of derivatives as well. In 3D, e.g., one may apply
    !
    ! <verb>
    !     call fev2_addVectorToEvalList(revalVectors,rvelocity%RvectorBlock(1),1)   ! u1
    !     call fev2_addVectorToEvalList(revalVectors,rvelocity%RvectorBlock(2),1)   ! u2
    !     call fev2_addVectorToEvalList(revalVectors,rvelocity%RvectorBlock(3),1)   ! u3
    ! </verb>
    !
    ! which calculates function values as well as 1st derivatives of the complete
    ! vector field (due to the "1" at the end). The calculated values in the
    ! cubature points can then be found in the "p_Ddata" elements of revalVectors:
    !
    ! <verb>
    !   revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)      = u1
    !   revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV3D_X) = d/dx u1
    !   revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV3D_Y) = d/dy u1
    !   revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV3D_Z) = d/dz u1
    !
    !   revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_FUNC)      = u2
    !   revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_DERIV3D_X) = d/dx u2
    !   revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_DERIV3D_Y) = d/dy u2
    !   revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_DERIV3D_Z) = d/dz u2
    !
    !   revalVectors%p_RvectorData(3)%p_Ddata(:,:,DER_FUNC)      = u3
    !   revalVectors%p_RvectorData(3)%p_Ddata(:,:,DER_DERIV3D_X) = d/dx u3
    !   revalVectors%p_RvectorData(3)%p_Ddata(:,:,DER_DERIV3D_Y) = d/dy u3
    !   revalVectors%p_RvectorData(3)%p_Ddata(:,:,DER_DERIV3D_Z) = d/dz u3
    ! </verb>
    !
    ! in all cubature points on all elements. The vector data in
    ! revalVectors%p_RvectorData appears exactly in the order, the vectors
    ! are added to revalVectors by fev2_addVectorToEvalList.
    !
    ! Remark 3:
    ! The routine currently assumes that all velocity components are discretised
    ! with the same FEM space.
    !
    ! Remark 4:
    ! The routine currently assumes that all velocity matrices are independent.
    ! Matrices sharing data are not supported. This cound be realised by
    ! taking care of the flags RmatrixData(:,:)%bsharedMatrixData which indicate
    ! which matrix data is shared.
    !
    ! Remark 4:
    ! Interleaved matrices are currently not supported.
!</remarks>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<subroutine>

    ! Local variables
    real(DP) :: dbasI, dbasJx, dbasJy, dbasJz
    integer :: iel, icubp, idofe, jdofe
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix11,p_DlocalMatrix22,p_DlocalMatrix33
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData11,p_rmatrixData22,p_rmatrixData33
    real(DP), dimension(:,:,:), pointer :: p_Du1,p_Du2,p_Du3

    integer :: ndim
    real(DP) :: dscale
    real(DP) :: dvelX, dvelY, dvelZ
    logical :: bvelConst

    ! Dimension of the underlying space
    ndim = rmatrixAssembly%p_rtriangulation%ndim
    
    ! Get parameters
    dscale = 1.0_DP
    dvelX = 0.0_DP
    dvelY = 0.0_DP
    dvelZ = 0.0_DP
    bvelConst = .true.
    
    if (present(rcollection)) then
      dscale = rcollection%DquickAccess(1)
      
      ! Constant velocity?
      bvelConst = (rcollection%IquickAccess(1) .eq. 0)
      
      if (bvelConst) then
        ! Get the constant velocity
        dvelX = rcollection%DquickAccess(2)
        if (ndim .ge. 2) dvelY = rcollection%DquickAccess(3)
        if (ndim .ge. 3) dvelZ = rcollection%DquickAccess(4)
      end if
      
    end if
    
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get local data
    p_DbasTrial => RmatrixData(1,1)%p_DbasTrial
    p_DbasTest => RmatrixData(1,1)%p_DbasTest

    ! Set up the local matrix of the convection.
    select case (ndim)
    case (NDIM1D)
      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(1,1)
      
      p_DlocalMatrix11 => RmatrixData(1,1)%p_Dentry
      
      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_convection")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJx = p_DbasTrial(jdofe,DER_DERIV1D_X,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI )            ! ( u1 phi_x , psi )

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel
        
      else

        if (revalVectors%ncount .lt. 1) then
          call output_line ("FEM function missing.",&
              OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_convection")
          call sys_halt()
        end if

        ! Set up the matrix for nonconstant velocity.
        !
        ! Get the velocity field from the parameters
        p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX = p_Du1(icubp,iel,DER_FUNC)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJx = p_DbasTrial(jdofe,DER_DERIV1D_X,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI )            ! ( u1 phi_x , psi )

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel
        
      end if
      
    case (NDIM2D)

      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(1,1)
      p_rmatrixData22 => RmatrixData(2,2)
      
      p_DlocalMatrix11 => RmatrixData(1,1)%p_Dentry
      p_DlocalMatrix22 => RmatrixData(2,2)%p_Dentry

      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_convection")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
                dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI )            ! ( u2 phi_y , psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI )            ! ( u2 phi_y , psi )

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel
        
      else

        if (revalVectors%ncount .lt. 2) then
          call output_line ("FEM function missing.",&
              OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_convection")
          call sys_halt()
        end if

        ! Set up the matrix for nonconstant velocity.
        !
        ! Get the velocity field from the parameters
        p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
        p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX = p_Du1(icubp,iel,DER_FUNC)
            dvelY = p_Du2(icubp,iel,DER_FUNC)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
                dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI )            ! ( u2 phi_y , psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI )            ! ( u2 phi_y , psi )

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel
        
      end if

    case (NDIM3D)

      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(1,1)
      p_rmatrixData22 => RmatrixData(2,2)
      p_rmatrixData33 => RmatrixData(3,3)

      p_DlocalMatrix11 => RmatrixData(1,1)%p_Dentry
      p_DlocalMatrix22 => RmatrixData(2,2)%p_Dentry
      p_DlocalMatrix33 => RmatrixData(3,3)%p_Dentry

      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_convection")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJx = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                dbasJy = p_DbasTrial(jdofe,DER_DERIV3D_Y,icubp,iel)
                dbasJz = p_DbasTrial(jdofe,DER_DERIV3D_Z,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI + &          ! ( u2 phi_y , psi )
                      dvelZ * dbasJz * dbasI )            ! ( u3 phi_z , psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI + &          ! ( u2 phi_y , psi )
                      dvelZ * dbasJz * dbasI )            ! ( u3 phi_z , psi )

                p_DlocalMatrix33(jdofe,idofe,iel) = p_DlocalMatrix33(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI + &          ! ( u2 phi_y , psi )
                      dvelZ * dbasJz * dbasI )            ! ( u3 phi_z , psi )

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel
        
      else

        if (revalVectors%ncount .lt. 3) then
          call output_line ("FEM function missing.",&
              OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_bubbleH1error")
          call sys_halt()
        end if

        ! Set up the matrix for nonconstant velocity.
        !
        ! Get the velocity field from the parameters
        p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
        p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata
        p_Du3 => revalVectors%p_RvectorData(3)%p_Ddata
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX = p_Du1(icubp,iel,DER_FUNC)
            dvelY = p_Du2(icubp,iel,DER_FUNC)
            dvelZ = p_Du3(icubp,iel,DER_FUNC)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJx = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                dbasJy = p_DbasTrial(jdofe,DER_DERIV3D_Y,icubp,iel)
                dbasJz = p_DbasTrial(jdofe,DER_DERIV3D_Z,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI + &          ! ( u2 phi_y , psi )
                      dvelZ * dbasJz * dbasI )            ! ( u3 phi_z , psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI + &          ! ( u2 phi_y , psi )
                      dvelZ * dbasJz * dbasI )            ! ( u3 phi_z , psi )

                p_DlocalMatrix33(jdofe,idofe,iel) = p_DlocalMatrix33(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI + &          ! ( u2 phi_y , psi )
                      dvelZ * dbasJz * dbasI )            ! ( u3 phi_z , psi )

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel
        
      end if

    end select

  end subroutine

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_rhsOne(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f=1 in all components.
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements
    
    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>
    
!<subroutine>

    ! Local variables
    real(DP) :: dbasI, dval
    integer :: icomp
    integer :: iel, icubp, idofe
    real(DP), dimension(:,:), pointer :: p_DlocalVector
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Loop over the components
    do icomp = 1,size(rvectorData)

      ! Get the data arrays of the subvector
      p_rvectorData => RvectorData(icomp)
      p_DlocalVector => RvectorData(icomp)%p_Dentry
      p_DbasTest => RvectorData(icomp)%p_DbasTest
    
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! The value of the RHS is just the 1.0-function here
          dval = 1.0_DP
          
          ! Outer loop over the DOF's i=1..ndof on our current element,
          ! which corresponds to the (test) basis functions Phi_i:
          do idofe=1,p_rvectorData%ndofTest
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! into dbasI
            dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
            
            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local vectors.
            p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
                p_DcubWeight(icubp,iel) * dval * dbasI
            
          end do ! jdofe

        end do ! icubp
      
      end do ! iel
      
    end do ! icomp
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_rhsBubble(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f=32*y*(1-y)+32*x*(1-x).
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements
    
    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>
    
!<subroutine>

    ! Local variables
    real(DP) :: dbasI, dval, dx, dy
    integer :: icomp
    integer :: iel, icubp, idofe
    real(DP), dimension(:,:), pointer :: p_DlocalVector
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Loop over the components
    do icomp = 1,size(rvectorData)

      ! Get the data arrays of the subvector
      p_rvectorData => RvectorData(icomp)
      p_DlocalVector => RvectorData(icomp)%p_Dentry
      p_DbasTest => RvectorData(icomp)%p_DbasTest
    
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement
        
          ! Get the coordinates of the cubature point.
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)

          ! Calculate the values of the RHS using the coordinates
          ! of the cubature points.
          dval = 32.0_DP*dy*(1.0_DP-dy) + 32_DP*dx*(1.0_DP-dx)
          
          ! Outer loop over the DOF's i=1..ndof on our current element,
          ! which corresponds to the (test) basis functions Phi_i:
          do idofe=1,p_rvectorData%ndofTest
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! into dbasI
            dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
            
            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local vectors.
            p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
                p_DcubWeight(icubp,iel) * dval * dbasI
            
          end do ! jdofe

        end do ! icubp
      
      end do ! iel
      
    end do ! icomp
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_rhsBubblePlusFE(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f = 32*y*(1-y)+32*x*(1-x) + v(x,y)
    ! with v being a FEM function.
    !
    ! The FEM function v(x,y) must be provided in revalVectors.
    ! The routine only supports non-interleaved vectors.
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements
    
    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>
    
!<subroutine>

    ! Local variables
    real(DP) :: dbasI, dval, dx, dy
    integer :: icomp
    integer :: iel, icubp, idofe
    real(DP), dimension(:,:), pointer :: p_DlocalVector
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
    real(DP), dimension(:,:), pointer :: p_Dfunc
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    if (revalVectors%ncount .eq. 0) then
      call output_line ("FEM function missing.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_rhsBubblePlusFE")
      call sys_halt()
    end if
    
    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_Dfunc => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC2D)
    
    ! Loop over the components
    do icomp = 1,size(rvectorData)

      ! Get the data arrays of the subvector
      p_rvectorData => RvectorData(icomp)
      p_DlocalVector => RvectorData(icomp)%p_Dentry
      p_DbasTest => RvectorData(icomp)%p_DbasTest
    
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement
        
          ! Get the coordinates of the cubature point.
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)

          ! Calculate the values of the RHS in the cubature point:
          !     f = 32*y*(1-y)+32*x*(1-x) + v(x,y)
          dval = 32.0_DP*dy*(1.0_DP-dy) + 32_DP*dx*(1.0_DP-dx) + p_Dfunc(icubp,iel)
          
          ! Outer loop over the DOF's i=1..ndof on our current element,
          ! which corresponds to the (test) basis functions Phi_i:
          do idofe=1,p_rvectorData%ndofTest
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! into dbasI
            dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
            
            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local vectors.
            p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
                p_DcubWeight(icubp,iel) * dval * dbasI
            
          end do ! jdofe

        end do ! icubp
      
      end do ! iel
      
    end do ! icomp
    
  end subroutine

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_integralOne(dintvalue,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the integral of the function v=1.
    ! The result is the size of the domain.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rvectorAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!<subroutine>

    ! Local variables
    real(DP) :: dval
    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    dintvalue = 0.0_DP
    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! The value of the function v is just the 1.0-function here
        dval = 1.0_DP
        
        ! Multiply the values by the cubature weight and sum up
        ! into the integral value
        dintvalue = dintvalue + p_DcubWeight(icubp,iel) * dval
          
      end do ! icubp
    
    end do ! iel
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_integralFE(dintvalue,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the integral of an arbitrary finite element function v.
    ! If v has multiple components, the sum of the integrals of all
    ! components is returned.
    !
    ! The FEM function(s) must be provided in revalVectors.
    ! The routine only supports non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rvectorAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!<subroutine>

    ! Local variables
    real(DP) :: dval
    integer :: iel, icubp, ivec
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    dintvalue = 0.0_DP

    ! Loop through all provided FEM functions
    do ivec = 1,revalVectors%ncount
    
      ! Skip interleaved vectors.
      if (revalVectors%p_RvectorData(ivec)%bisInterleaved) cycle

      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_Dfunc => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_FUNC2D)

      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Get the value of the FEM function
          dval = p_Dfunc(icubp,iel)
          
          ! Multiply the values by the cubature weight and sum up
          ! into the integral value
          dintvalue = dintvalue + p_DcubWeight(icubp,iel) * dval
            
        end do ! icubp
      
      end do ! iel
      
    end do ! ivec
    
  end subroutine

  !****************************************************************************

!<--
! The L2-norm of a vector field "rvector" can be calculated, e.g., as follows:
!
!    ! Declare an evaluation structure for the subvectors
!    type(t_fev2Vectors) :: revalVectors
!    ...
!    real(DP) :: dresult
!
!    ! Add the subvectors
!    call fev2_addVectorToEvalList(revalVectors,rvector%RvectorBlock(1),0)
!    call fev2_addVectorToEvalList(revalVectors,rvector%RvectorBlock(2),0)
!    ...
!
!    ! Evaluate the norm to "dresult"
!    call bma_buildIntegral(dresult,BMA_CALC_STANDARD,bma_fcalc_L2norm,&
!        revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
!
!    ! Clean up
!    call fev2_releaseVectorList(revalVectors)
!-->

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_L2norm(dintvalue,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the squared L2-norm of an arbitrary finite element 
    ! function: <tex> ||v||^2 </tex>.
    ! If v has multiple components, the sum of the integrals of all
    ! components is returned.
    !
    ! The FEM function(s) must be provided in revalVectors.
    ! The routine only supports non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rvectorAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!<subroutine>

    ! Local variables
    real(DP) :: dval
    integer :: iel, icubp, ivec
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    dintvalue = 0.0_DP

    ! Loop through all provided FEM functions
    do ivec = 1,revalVectors%ncount
    
      ! Skip interleaved vectors.
      if (revalVectors%p_RvectorData(ivec)%bisInterleaved) cycle

      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_Dfunc => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_FUNC)

      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Get the value of the FEM function
          dval = p_Dfunc(icubp,iel)
          
          ! Multiply the values by the cubature weight and sum up
          ! into the integral value
          dintvalue = dintvalue + p_DcubWeight(icubp,iel) * dval**2
            
        end do ! icubp
      
      end do ! iel
      
    end do ! ivec
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_H1norm(dintvalue,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the squared H1 (semi-)norm of an arbitrary finite element 
    ! function <tex> ||v||^2_{H^1} </tex>.
    ! If v has multiple components, the sum of the integrals of all
    ! components is returned.
    !
    ! The FEM function(s) must be provided in revalVectors.
    ! The routine only supports non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rvectorAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!<subroutine>

    ! Local variables
    real(DP) :: dderivX,dderivY,dderivZ
    integer :: iel, icubp, ivec,ndim
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_DderivX,p_DderivY,p_DderivZ
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    dintvalue = 0.0_DP

    ! Loop through all provided FEM functions
    do ivec = 1,revalVectors%ncount
    
      ! Skip interleaved vectors.
      if (revalVectors%p_RvectorData(ivec)%bisInterleaved) cycle

      ! Dimension?
      ndim = revalVectors%p_RvectorData(ivec)%p_rvector%p_rspatialDiscr%p_rtriangulation%ndim
      select case (ndim)
      
      ! ----------------
      ! 1D
      ! ----------------
      case (NDIM1D)
      
        ! Get the data array with the values of the FEM function
        ! in the cubature points
        p_DderivX => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_DERIV1D_X)

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Get the value of the FEM function
            dderivX = p_DderivX(icubp,iel)
            
            ! Multiply the values by the cubature weight and sum up
            ! into the integral value
            dintvalue = dintvalue + p_DcubWeight(icubp,iel) * dderivX**2
              
          end do ! icubp
        
        end do ! iel

      ! ----------------
      ! 2D
      ! ----------------
      case (NDIM2D)
      
        ! Get the data array with the values of the FEM function
        ! in the cubature points
        p_DderivX => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_DERIV2D_X)
        p_DderivY => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_DERIV2D_Y)

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Get the value of the FEM function
            dderivX = p_DderivX(icubp,iel)
            dderivY = p_DderivY(icubp,iel)
            
            ! Multiply the values by the cubature weight and sum up
            ! into the integral value
            dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
                ( dderivX**2 + dderivY**2 )
              
          end do ! icubp
        
        end do ! iel

      ! ----------------
      ! 3D
      ! ----------------
      case (NDIM3D)
      
        ! Get the data array with the values of the FEM function
        ! in the cubature points
        p_DderivX => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_DERIV3D_X)
        p_DderivY => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_DERIV3D_Y)
        p_DderivZ => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_DERIV3D_Z)

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Get the value of the FEM function
            dderivX = p_DderivX(icubp,iel)
            dderivY = p_DderivY(icubp,iel)
            dderivZ = p_DderivZ(icubp,iel)
            
            ! Multiply the values by the cubature weight and sum up
            ! into the integral value
            dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
                ( dderivX**2 + dderivY**2 + dderivZ**2 )
              
          end do ! icubp
        
        end do ! iel
        
      end select
      
    end do ! ivec
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_bubbleL2error(dintvalue,rassemblyData,rintegralAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the (squared) L2 error of an arbitrary FEM function v
    ! to the bubble function u=16x(1-x)y(1-y):
    !
    !  <tex>  || v - u ||^2_{L2}  </tex>
    !
    ! The FEM function must be provided in revalVectors.
    ! The routine only supports non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintegralAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!<subroutine>

    ! Local variables
    real(DP) :: dval1,dval2,dx,dy
    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
    ! Cancel if no FEM function is given.
    if (revalVectors%ncount .eq. 0) then
      call output_line ("FEM function missing.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_bubbleL2error")
      call sys_halt()
    end if

    ! Skip interleaved vectors.
    if (revalVectors%p_RvectorData(1)%bisInterleaved) return

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_Dfunc => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC2D)

    dintvalue = 0.0_DP
    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Get the value of the bubble function
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)
        
        dval1 = 16.0_DP * dx * (1.0_DP-dx) * dy * (1.0_DP-dy)

        ! Get the error of the FEM function to the bubble function
        dval2 = p_Dfunc(icubp,iel)
        
        ! Multiply the values by the cubature weight and sum up
        ! into the (squared) L2 error:
        dintvalue = dintvalue + &
            p_DcubWeight(icubp,iel) * (dval1 - dval2)**2
          
      end do ! icubp
    
    end do ! iel
      
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_bubbleH1error(dintvalue,rassemblyData,rintegralAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the (squared) H1 error of an arbitrary FEM function v
    ! to the bubble function u=16x(1-x)y(1-y)
    ! (based on the H1 semi-norm).
    !
    !  <tex>  | v - u |^2_{H1}  </tex>
    !
    ! The FEM function must be provided in revalVectors.
    ! The routine only supports non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintegralAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!<subroutine>

    ! Local variables
    real(DP) :: dderivX1,dderivY1,dderivX2,dderivY2,dx,dy
    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_DderivX,p_DderivY
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
    ! Calcel if no FEM function is given.
    if (revalVectors%ncount .eq. 0) then
      call output_line ("FEM function missing.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_bubbleH1error")
      call sys_halt()
    end if

    ! Skip interleaved vectors.
    if (revalVectors%p_RvectorData(1)%bisInterleaved) return

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_DderivX => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV2D_X)
    p_DderivY => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV2D_Y)

    dintvalue = 0.0_DP
    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Get the derivatives of the bubble function in the cubature point
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)
        
        dderivX1 = 16.0_DP*dy*(dy-1.0_DP)*(2.0_DP*dx-1.0_DP)
        dderivY1 = 16.0_DP*dx*(dx-1.0_DP)*(2.0_DP*dy-1.0_DP)

        ! Get the error of the FEM function derivatives of the bubble function
        ! in the cubature point
        dderivX2 = p_DderivX(icubp,iel)
        dderivY2 = p_DderivY(icubp,iel)
        
        ! Multiply the values by the cubature weight and sum up
        ! into the (squared) H1 error:
        dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
            ( (dderivX1 - dderivX2)**2 + (dderivY1 - dderivY2)**2 )
          
      end do ! icubp
    
    end do ! iel
      
  end subroutine

  !****************************************************************************
  !****************************************************************************
  ! Assembly of block matrices
  !****************************************************************************
  !****************************************************************************

  !<subroutine>

  subroutine bma_getLocalMatrixIndices (rmatrix,Irows,Icolumns,Kentry,&
      irowsPerElement,icolsPerElement,nelements)

  !<description>

  ! Calculates index positions of local matrices in a global matrix.
  ! For a set of elements, Icolumns and Irows define the row and column indices
  ! of local matrices which have to be accessed in a global matrix rmatrix.
  ! The routine then calculates the positions of all the corresponding matrix
  ! entries in the data array of the matrix rmatrix and saves the result
  ! to the array Kentry.
  !
  ! IMPORTANT: For performance reasons, the columns and rows in Kentry
  ! are saved *transposed* in comparison to the matrix rmatrix!
  ! That means that
  !    Kentry(j,i,:) = position of element (Irows(i,:),Icolumns(j,:))
  ! holds!

  !</description>

  !<input>

  ! The global matrix which has to be accessed.
  type(t_matrixScalar), intent(in) :: rmatrix

  ! Array identifying all rows in the global matrix which have to be
  ! accessed.
  ! DIMENSION(#rows per element, #elements).
  integer, dimension(:,:), intent(in) :: Irows

  ! Array identifying all columns in the global matrix which have to be
  ! accessed.
  ! DIMENSION(#columns per element, #elements).
  integer, dimension(:,:), intent(in) :: Icolumns

  ! Number of rows per element / in the local matrix
  integer, intent(in) :: irowsPerElement

  ! Number of columns per element / in the local matrix
  integer, intent(in) :: icolsPerElement

  ! Number of elements.
  integer, intent(in) :: nelements

  !</input>

  !<output>

  ! Array receiving the positions of the local matrices in the global matrix.
  ! DIMENSION(#columns per element,#rows per element,#elements).
  ! Saved in a transposed way:
  !    Kentry(j,i,:) = position of element (Irows(i,:),Icolumns(j,:))
  integer, dimension(:,:,:), intent(out) :: Kentry

  !</output>

  !</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Kcol, p_Kld, p_KrowIdx
    integer :: na,iel,idofe,jdofe,ndofTest,ndofTrial,jcol0,jdfg,jcol,nnzrows

    ndofTrial = icolsPerElement
    ndofTest = irowsPerElement

    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX1)

      ! That is easy, we can directly calculate the positions
      do iel = 1,nelements
        do idofe = 1,ndofTest
          do jdofe = 1,ndofTrial
            Kentry(jdofe,idofe,iel) = &
                Irows(idofe,iel) * rmatrix%NCOLS + Icolumns(jdofe,iel)
          end do
        end do
      end do

    case (LSYSSC_MATRIX7,LSYSSC_MATRIX9)

      ! Get pointers to the row/column structure of the matrix
      call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      call lsyssc_getbase_Kld (rmatrix,p_Kld)
      na = rmatrix%NA

      ! We build a quadratic ndofTrial*ndofTest local matrix:
      ! Kentry(1..ndofTrial,1..ndofTest) receives the position
      !   in the global system matrix
      !
      ! Loop through elements in the set and for each element,
      ! loop through the local matrices to initialise them:
      do iel = 1,nelements

        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"`s), as these
        ! define the rows in the matrix.
        do idofe = 1,ndofTest

          ! Row IDOFE of the local matrix corresponds
          ! to row=global DOF KDFG(IDOFE) in the global matrix.
          ! This is one of the the "O"`s in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          jcol0=p_Kld(Irows(idofe,iel))

          ! Now we loop through the other DOF`s on the current element
          ! (the "O"`s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.

          do jdofe = 1,ndofTrial

            ! Get the global DOF of the "X" which interacts with
            ! our "O".

            jdfg=Icolumns(jdofe,iel)

            ! Starting in JCOL0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the DO loop if we find the column.

            do jcol = jcol0,na
              if (p_Kcol(jcol) .eq. jdfg) exit
            end do

            ! Because columns in the global matrix are sorted
            ! ascendingly (except for the diagonal element),
            ! the next search can start after the column we just found.

            ! JCOL0=JCOL+1

            ! Save the position of the matrix entry into the local
            ! matrix.
            ! Note that a column in Kentry corresponds to a row in
            ! the real matrix. We aligned Kentry this way to get
            ! higher speed of the assembly routine, since this leads
            ! to better data locality.

            Kentry(jdofe,idofe,iel)=jcol

          end do ! IDOFE

        end do ! JDOFE

      end do ! IEL

    case (LSYSSC_MATRIX9ROWC)

      ! Get pointers to the row/column structure of the matrix
      call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      call lsyssc_getbase_Kld (rmatrix,p_Kld)
      call lsyssc_getbase_KrowIdx (rmatrix,p_KrowIdx)
      na = rmatrix%NA
      nnzrows = rmatrix%NNZROWS

      ! We build a quadratic ndofTrial*ndofTest local matrix:
      ! Kentry(1..ndofTrial,1..ndofTest) receives the position
      !   in the global system matrix
      !
      ! Loop through elements in the set and for each element,
      ! loop through the local matrices to initialise them:
      do iel = 1,nelements

        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"`s), as these
        ! define the rows in the matrix.
        do idofe = 1,ndofTest

          ! Row IDOFE of the local matrix corresponds
          ! to row=global DOF KDFG(IDOFE) in the global matrix.
          ! This is one of the the "O"`s in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          jcol0=p_Kld(p_KrowIdx(nnzrows+Irows(idofe,iel)))

          ! Now we loop through the other DOF`s on the current element
          ! (the "O"`s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.

          do jdofe = 1,ndofTrial

            ! Get the global DOF of the "X" which interacts with
            ! our "O".

            jdfg=Icolumns(jdofe,iel)

            ! Starting in JCOL0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the DO loop if we find the column.

            do jcol = jcol0,na
              if (p_Kcol(jcol) .eq. jdfg) exit
            end do

            ! Because columns in the global matrix are sorted
            ! ascendingly (except for the diagonal element),
            ! the next search can start after the column we just found.

            ! JCOL0=JCOL+1

            ! Save the position of the matrix entry into the local
            ! matrix.
            ! Note that a column in Kentry corresponds to a row in
            ! the real matrix. We aligned Kentry this way to get
            ! higher speed of the assembly routine, since this leads
            ! to better data locality.
            ! Subtract the offset to get the offset in the compressed matrix.

            Kentry(jdofe,idofe,iel)=jcol

          end do ! IDOFE

        end do ! JDOFE

      end do ! IEL

    end select

  end subroutine

  !****************************************************************************

  subroutine addInt (Ilist,nentries,ientry)
  ! Adds ientry to a list represented by an array
  integer, dimension(:), intent(inout) :: Ilist
  integer, intent(in), target :: ientry
  integer, intent(inout) :: nentries
    nentries = nentries + 1
    Ilist(nentries) = ientry
  end subroutine

  integer function containsInt (Ilist,nentries,ientry)
  ! returns <> 0 (the index) if the list Ilist contains Ientry
  integer, dimension(:), intent(in) :: Ilist
  integer, intent(in), target :: ientry
  integer, intent(in) :: nentries

    integer :: i

    do i=1,nentries
      if (Ilist(i) .eq. ientry) then
        containsInt = i
        return
      end if
    end do

    containsInt = 0

  end function

  !****************************************************************************
  
  integer function containsDiscr (rfemDataBlocks,rentry)
  ! returns <> 0 (the index) if the list Rlist contains rentry
  type(t_fev2FemDataBlocks), intent(in) :: rfemDataBlocks
  type(t_spatialDiscretisation), intent(in), target :: rentry
  
    integer :: i
    
    do i=1,rfemDataBlocks%ncount
      if (associated(rfemDataBlocks%p_RfemData(i)%p_rdiscr,rentry)) then
        containsDiscr = i
        return
      end if
    end do
    
    containsDiscr = 0
    
  end function
    

  !****************************************************************************

!<subroutine>

  subroutine bma_prepareMatrixData(&
      rmatrix,RmatrixData,rfemDataBlocks,ielementDistr)

!<description>
  ! Initialise a matrix assembly structure for assembling the submatrices
  ! of a block matrix.
  ! No memory is allocated.
!</description>

!<input>
  ! The matrix which is going to be assembled.
  type(t_matrixBlock), intent(in) :: rmatrix

  ! ID of the active element distribution. Identifies the FEM space
  ! in the matrix which is to be assembled.
  integer, intent(in) :: ielementDistr

  ! Data of all involved FEM spaces.
  type(t_fev2FemDataBlocks) , intent(in) :: rfemDataBlocks
!</input>

!<output>
  ! An array of block assembly structures for all the submatrices
  ! in the block matrix.
  type(t_bmaMatrixData), dimension(:,:), intent(out) :: RmatrixData
!</output>

!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscrTrial,p_rdiscrTest
    integer :: i,j

    ! List of used matrix handles for data and structure
    integer, dimension(:), pointer :: p_ImatrixDataHandles
    integer :: nmatrixDataHandles

    integer, dimension(:), pointer :: p_ImatrixStrucHandles
    integer :: nmatrixStrucHandles

    ! Allocate memory for existence checks
    i = size(rmatrix%RmatrixBlock)
    allocate (p_ImatrixDataHandles(i))
    allocate (p_ImatrixStrucHandles(i))
    nmatrixDataHandles = 0
    nmatrixStrucHandles = 0

    ! Basic initialisation of the structure.
    ! Loop over all blocks. Figure out which blocks have data.
    ! Get the data arrays for these blocks.
    do j=1,ubound(RmatrixData,2)

      do i=1,ubound(RmatrixData,1)

        RmatrixData(i,j)%bhasData = lsysbl_isSubmatrixPresent (rmatrix,i,j)

        if (RmatrixData(i,j)%bhasData) then

          ! Remember the matrix
          RmatrixData(i,j)%p_rmatrix => rmatrix%RmatrixBlock(i,j)

          ! Get the pointers to the matrix arrays
          select case (rmatrix%RmatrixBlock(i,j)%cmatrixFormat)
          case (LSYSSC_MATRIX9,LSYSSC_MATRIX7,LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7INTL)

            if (rmatrix%RmatrixBlock(i,j)%cdataType .ne. ST_DOUBLE) then
              call output_line ("Matrix data type not supported",&
                  OU_CLASS_ERROR,OU_MODE_STD,"bma_prepareMatrixData")
              call sys_halt()
            end if

            if (.not. lsyssc_hasMatrixContent(rmatrix%RmatrixBlock(i,j))) then
              call output_line ("Empty submatrix, no data allocated!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"bma_prepareMatrixData")
              call sys_halt()
            end if

            ! Get a pointer to the matrix data for faster access
            call lsyssc_getbase_double (rmatrix%RmatrixBlock(i,j),RmatrixData(i,j)%p_Da)

            ! Check if the matrix data is shared with another matrix.
            ! If not, remember it.
            RmatrixData(i,j)%bsharedMatrixData = 0 .ne. &
              containsInt (p_ImatrixDataHandles,nmatrixDataHandles,&
                           rmatrix%RmatrixBlock(i,j)%h_Da)

            if (.not. RmatrixData(i,j)%bsharedMatrixData) then
              ! Remember the matrix data for later checks
              call addInt (p_ImatrixDataHandles,nmatrixDataHandles,&
                      rmatrix%RmatrixBlock(i,j)%h_Da)
            end if

            ! Check if the structure is shared based on Kcol
            RmatrixData(i,j)%bsharedMatrixStructure = 0 .ne. &
              containsInt (p_ImatrixStrucHandles,nmatrixStrucHandles,&
                           rmatrix%RmatrixBlock(i,j)%h_Kcol)

            if (.not. RmatrixData(i,j)%bsharedMatrixStructure) then
              ! Remember the matrix data for later checks
              call addInt (p_ImatrixStrucHandles,nmatrixStrucHandles,&
                      rmatrix%RmatrixBlock(i,j)%h_Kcol)
            end if

            ! Matrix subtype: Interleaved matrices.
            select case (rmatrix%RmatrixBlock(i,j)%cmatrixFormat)
            case (LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7INTL)
              ! Support for interleaved matrices
              RmatrixData(i,j)%nvar = rmatrix%RmatrixBlock(i,j)%nvar
              RmatrixData(i,j)%bisInterleaved = .true.
            end select

          case default
            call output_line ("Matrix format not supported",&
                OU_CLASS_ERROR,OU_MODE_STD,"bma_prepareFemData")
            call sys_halt()

          end select

          ! Get the information about the trial and test spaces
          ! of this block.
          p_rdiscrTrial => rmatrix%RmatrixBlock(i,j)%p_rspatialDiscrTrial
          p_rdiscrTest  => rmatrix%RmatrixBlock(i,j)%p_rspatialDiscrTest

          ! Some basic checks...
          if ((p_rdiscrTrial%inumFESpaces .lt. ielementDistr) .or. &
              (p_rdiscrTrial%inumFESpaces .lt. ielementDistr)) then
            call output_line ("Element distribution does not exist in the discretisation.",&
                OU_CLASS_ERROR,OU_MODE_STD,"bma_prepareMatrixData")
            call sys_halt()
          end if

          ! Get the indices of the trial and the test space.
          RmatrixData(i,j)%iidxFemDataTrial = containsDiscr(rfemDataBlocks,p_rdiscrTrial)
          RmatrixData(i,j)%iidxFemDataTest = containsDiscr (rfemDataBlocks,p_rdiscrTest)
          RmatrixData(i,j)%bIdenticalTrialAndTest = associated(p_rdiscrTrial,p_rdiscrTest)

          ! Element type          
          RmatrixData(i,j)%celementTrial = &
              p_rdiscrTrial%RelementDistr(ielementDistr)%celement

          RmatrixData(i,j)%celementTest = &
              p_rdiscrTest%RelementDistr(ielementDistr)%celement

          ! Get the number of local DOF`s for trial and test functions
          RmatrixData(i,j)%ndofTrial = elem_igetNDofLoc(RmatrixData(i,j)%celementTrial)
          RmatrixData(i,j)%ndofTest = elem_igetNDofLoc(RmatrixData(i,j)%celementTest)

        end if

      end do

    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_initMatAssemblyData(rassemblyData,&
      rmatrix,ielementDistr,RmaxDerivativeTest,RmaxDerivativeTrial,revalVectors)

!<description>
  ! Initialises the main parameters in the structure containing the
  ! assembly data (cubature weights, etc.).
  !
  ! The returned rassemblyData structure can be used as template
  ! for bma_createMatAssemblyData.
!</description>

!<input>
  ! The matrix which is going to be assembled.
  type(t_matrixBlock), intent(in) :: rmatrix

  ! ID of the active element distribution. Identifies the FEM space
  ! in the matrix which is to be assembled.
  integer, intent(in) :: ielementDistr

  ! OPTIONAL: For every block in the matrix, maximum
  ! derivative of the basis functions to be computed. If not
  ! specified or an entry is =-1, the maximum available derivative for 
  ! each FEM space is the default.
  integer, dimension(:,:), intent(in), optional :: RmaxDerivativeTest
  integer, dimension(:,:), intent(in), optional :: RmaxDerivativeTrial

  ! OPTIONAL: A set of FEM functions to be evaluated during the vector
  ! assembly.
  type(t_fev2Vectors), intent(in), optional :: revalVectors
!</input>

!<output>
  ! The matrix assembly data structure to be initialised.
  type(t_bmaMatrixAssemblyData), intent(out) :: rassemblyData
!</output>

!</subroutine>

    ! Initialise the element set; can be simultaneously used for all
    ! FEM spaces, it is element independent.
    call elprep_init(rassemblyData%revalElementSet)

    ! Element sets are not yet initialised.  
    rassemblyData%ninitialisedElements = 0

    ! Prepare FEM data
    call fev2_prepareFemDataBMat(rmatrix,rassemblyData%rfemDataBlocks,&
        ielementDistr,RmaxDerivativeTest,RmaxDerivativeTrial)

    ! Prepare FEM data for the vectors to be evaluated.
    if (present(revalVectors)) then
      call fev2_prepareFemDataVecEval(revalVectors,rassemblyData%rfemDataBlocks,&
          ielementDistr)
    end if

    ! Initialise the blocks according to the structure of the matrix
    ! and the FEM spaces used.
    allocate(rassemblyData%p_RmatrixData( &
        ubound(rmatrix%RmatrixBlock,1),ubound(rmatrix%RmatrixBlock,2)))

    call bma_prepareMatrixData(&
        rmatrix,rassemblyData%p_RmatrixData,rassemblyData%rfemDataBlocks,ielementDistr)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_doneMatAssemblyData(rassemblyData)

!<description>
  ! Releases allocated memory from the assembly data structure.
!</description>

!<inputoutput>
  ! Assembly data structure to be cleaned up.
  type(t_bmaMatrixAssemblyData), intent(inout) :: rassemblyData
!</inputoutput>

!</subroutine>

    ! Release the matrix and FEM data array
    call fev2_cleanupFemData(rassemblyData%rfemDataBlocks)
    deallocate(rassemblyData%p_RmatrixData)

    ! Element sets are not yet initialised.  
    rassemblyData%ninitialisedElements = 0

  end subroutine


  !****************************************************************************

!<subroutine>

  subroutine bma_createMatAssemblyData(rmatrixAssembly,&
      rassemblyData,rassemblyDataTemplate,&
      revalVectors,revalVectorsTemplate)

!<description>
  ! Based on a template structure, this creates an assembly data structure
  ! for the matrix assembly with all data arrays initialised and allocated.
!</description>

!<input>
  ! Template assembly data structure. rassemblyData is created based
  ! on the data in this structure.
  type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyDataTemplate

  ! A matrix assembly structure.
  type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

  ! Template vector data structure. rvectorData is created based
  ! on the data in this structure.
  type(t_fev2Vectors), intent(in) :: revalVectorsTemplate
!</input>

!<output>
  ! The matrix assembly data structure to be initialised.
  type(t_bmaMatrixAssemblyData), intent(out) :: rassemblyData

  ! Vector data structure to be initialised according to revalVectorsTemplate
  type(t_fev2Vectors), intent(out) :: revalVectors
!</output>

!</subroutine>

    ! local variables
    integer :: i,j
    type(t_bmaMatrixData), pointer :: p_rmatrixData

    ! At first, just copy the content from the matrix template structure
    rassemblyData = rassemblyDataTemplate

    ! Now initialise the dynamic content.
    !
    ! Copy the FEM data and the matrix data
    allocate (rassemblyData%p_RmatrixData(&
        ubound(rassemblyDataTemplate%p_RmatrixData,1),&
        ubound(rassemblyDataTemplate%p_RmatrixData,2)))
    rassemblyData%p_RmatrixData(:,:) = rassemblyDataTemplate%p_RmatrixData(:,:)

    call fev2_copyFemData(rassemblyData%rfemDataBlocks,rassemblyDataTemplate%rfemDataBlocks)

    ! Initialise the FEM evaluation structures.
    call fev2_createFemData(rassemblyData%rfemDataBlocks,&
        rmatrixAssembly%ncubp,rmatrixAssembly%nelementsPerBlock)

    ! Loop through the matrices. Whereever there is unshared data,
    ! there is something to compute, so allocate memory for the matrix
    ! entries in this block.
    do j=1,ubound(rassemblyData%p_RmatrixData,2)
      do i=1,ubound(rassemblyData%p_RmatrixData,1)

        p_rmatrixData => rassemblyData%p_RmatrixData(i,j)

        ! Unshared data in that block
        if (p_rmatrixData%bhasData) then 

          if (.not. p_rmatrixData%bsharedMatrixData) then

            ! Allocate memory for the martix data and positions in the matrix
            if (.not. p_rmatrixData%bisInterleaved) then
              ! non-interleaved matrix
              allocate(p_rmatrixData%p_Dentry(&
                      p_rmatrixData%ndofTrial,p_rmatrixData%ndofTest,&
                      rmatrixAssembly%nelementsPerBlock))
            else
              ! interleaved matrix
              allocate(p_rmatrixData%p_DentryIntl(p_rmatrixData%nvar,&
                      p_rmatrixData%ndofTrial,p_rmatrixData%ndofTest,&
                      rmatrixAssembly%nelementsPerBlock))
            end if

            allocate(p_rmatrixData%p_Kentry(&
                    p_rmatrixData%ndofTrial,p_rmatrixData%ndofTest,&
                    rmatrixAssembly%nelementsPerBlock))

          else

            ! Nullify to mark this matrix data as "computed somewhere else".
            nullify(p_rmatrixData%p_Dentry)
            nullify(p_rmatrixData%p_DentryIntl)
            nullify(p_rmatrixData%p_Kentry)

          end if

          ! Map the pointers to the DOFs and the values of the
          ! basis functions into that block
          p_rmatrixData%p_DbasTrial => &
              rassemblyData%rfemDataBlocks%p_RfemData(p_rmatrixData%iidxFemDataTrial)%p_Dbas
          p_rmatrixData%p_DbasTest => &
              rassemblyData%rfemDataBlocks%p_RfemData(p_rmatrixData%iidxFemDataTest)%p_Dbas

          p_rmatrixData%p_IdofsTrial => &
              rassemblyData%rfemDataBlocks%p_RfemData(p_rmatrixData%iidxFemDataTrial)%p_Idofs
          p_rmatrixData%p_IdofsTest => &
              rassemblyData%rfemDataBlocks%p_RfemData(p_rmatrixData%iidxFemDataTest)%p_Idofs

        end if

      end do
    end do

    ! Initialise the vector evaluation structure.
    !
    ! Copy the content from the template.
    revalVectors = revalVectorsTemplate

    ! Re-create the vector array
    if (revalVectors%ncount .eq. 0) then
      nullify(revalVectors%p_RvectorData)
    else
      allocate(revalVectors%p_RvectorData(revalVectorsTemplate%ncount))
      revalVectors%p_RvectorData(1:revalVectorsTemplate%ncount) = &
          revalVectorsTemplate%p_RvectorData(1:revalVectorsTemplate%ncount)

      ! Prepare the evaluation of the FEM functions
      call fev2_initVectorEval(revalVectors,rassemblyData%rfemDataBlocks)
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_releaseMatAssemblyData(rassemblyData,revalVectors)

!<description>
  ! Releases memory of an assembly structure created by bma_createMatAssemblyData.
!</description>

!<inputoutput>
  ! The matrix assembly data structure to be cleaned up
  type(t_bmaMatrixAssemblyData), intent(inout) :: rassemblyData

  ! Vector data structure to be cleaned up
  type(t_fev2Vectors), intent(inout) :: revalVectors
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j
    type(t_bmaMatrixData), pointer :: p_rmatrixData

    ! Release vector evaluzation data
    if (revalVectors%ncount .ne. 0) then
      call fev2_doneVectorEval(revalVectors)
      deallocate(revalVectors%p_RvectorData)
    end if

    ! Release FEM evaluation data
    call fev2_releaseFemData(rassemblyData%rfemDataBlocks)

    ! Release the element set and the cubature weights
    call elprep_releaseElementSet(rassemblyData%revalElementSet)

    if (associated(rassemblyData%p_DcubWeight)) then
      deallocate(rassemblyData%p_DcubWeight)
    end if

    ! Element sets are not initialised anymore
    rassemblyData%ninitialisedElements = 0

    ! Loop through the matrices. Whereever there is unshared data,
    ! there is something to compute, so allocate memory for the matrix
    ! entries in this block.
    do j=1,ubound(rassemblyData%p_RmatrixData,2)
      do i=1,ubound(rassemblyData%p_RmatrixData,1)

        p_rmatrixData => rassemblyData%p_RmatrixData(i,j)

        ! Unshared data in that block
        if (p_rmatrixData%bhasData) then 

          if (.not. p_rmatrixData%bsharedMatrixData) then

            ! Allocate memory for the martix data and positions in the matrix
            if (associated(p_rmatrixData%p_Dentry)) then
              deallocate(p_rmatrixData%p_Dentry)
            end if

            if (associated(p_rmatrixData%p_DentryIntl)) then
              deallocate(p_rmatrixData%p_DentryIntl)
            end if

            deallocate(p_rmatrixData%p_Kentry)

          end if

        end if

      end do
    end do

    ! Release memory
    call fev2_cleanupFemData(rassemblyData%rfemDataBlocks)
    deallocate (rassemblyData%p_RmatrixData)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_initMatAssembly(rmatrixAssembly,ccubType,ctrafoType,cflags,&
      rmatrix,ielementDistr,RmaxDerivativeTest,RmaxDerivativeTrial,&
      revalVectors,rperfconfig)

!<description>
  ! Initialise a matrix assembly structure for assembling a bilinear form
  ! for a certain element distribution.
!</description>

!<input>
  ! Cubature formla to use
  integer(I32), intent(in) :: ccubType
  
  ! Transformation to use
  integer(I32), intent(in) :: ctrafoType

  ! Option field. Combination of BMA_CALC_xxxx flags.
  ! Use BMA_CALC_STANDARD for standard options.
  integer(I32), intent(in) :: cflags

  ! The matrix which is going to be assembled.
  type(t_matrixBlock), intent(in) :: rmatrix

  ! ID of the active element distribution. Identifies the FEM space
  ! in the matrix which is to be assembled.
  integer, intent(in) :: ielementDistr

  ! OPTIONAL: For every block in the matrix, maximum
  ! derivative of the basis functions to be computed. If not
  ! specified or an entry is =0, the maximum available derivative for 
  ! each FEM space is the default.
  integer, dimension(:,:), intent(in), optional :: RmaxDerivativeTest
  integer, dimension(:,:), intent(in), optional :: RmaxDerivativeTrial

  ! OPTIONAL: Set of vectors to be automatically evaluated
  type(t_fev2Vectors), intent(in), optional :: revalVectors

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<output>
  ! A matrix assembly structure.
  type(t_bmaMatrixAssembly), intent(out) :: rmatrixAssembly
!</output>

!</subroutine>

    integer :: i,ncount
    integer(I32) :: cshape
    type(t_fev2FemData), dimension(:), pointer :: p_RfemData

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Remember the performance configuration for later use
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => bma_perfconfig
    end if

    rmatrixAssembly%p_rperfconfig => p_rperfconfig

    ! Get basic data
    rmatrixAssembly%p_rboundary => rmatrix%p_rblockDiscrTest%p_rboundary
    rmatrixAssembly%p_rtriangulation => rmatrix%p_rblockDiscrTest%p_rtriangulation

    ! Save the vector evaluation structure as template.
    ! Make a copy of the stucture, that is simpler here.
    if (present(revalVectors)) then
      rmatrixAssembly%revalVectorsTemplate = revalVectors
    end if

    ! Remember current element distribution
    rmatrixAssembly%ielementDistr = ielementDistr

    ! Initialise the template structure which is used during the
    ! actual matrix assembly.
    call bma_initMatAssemblyData(rmatrixAssembly%rassemblyDataTemplate,&
        rmatrix,ielementDistr,RmaxDerivativeTest,RmaxDerivativeTrial,&
        revalVectors)

    ! Get the cubature formula/transformation on/from the reference element
    rmatrixAssembly%ccubType = ccubType
    rmatrixAssembly%ncubp = cub_igetNumPts(ccubType)
    rmatrixAssembly%ctrafoType = ctrafoType

    ! Allocate two arrays for the points and the weights
    allocate(rmatrixAssembly%p_Domega(rmatrixAssembly%ncubp))
    allocate(rmatrixAssembly%p_DcubPtsRef(&
        trafo_igetReferenceDimension(rmatrixAssembly%ctrafoType),&
        rmatrixAssembly%ncubp))

    ! Get the cubature formula
    call cub_getCubature(ccubType,rmatrixAssembly%p_DcubPtsRef,rmatrixAssembly%p_Domega)

    ! Number of simultaneously processed elements
    rmatrixAssembly%nelementsPerBlock = p_rperfconfig%NELEMSIM

    ! Get the evaluation tag by "OR"-ing all those from the FEM spaces
    ! and the option fields.
    rmatrixAssembly%cevaluationTag = 0

    if (iand(cflags,BMA_CALC_REALCOORDS) .ne. 0) then
      ! Calculate real world coordinates of the cubature points
      rmatrixAssembly%cevaluationTag = &
          ior(rmatrixAssembly%cevaluationTag,EL_EVLTAG_REALPOINTS)
    end if

    if (rmatrixAssembly%rassemblyDataTemplate%rfemDataBlocks%ncount .eq. 0) then
      call output_line ("Matrix is empty, nothing to assemble.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_initMatAssembly")
      call sys_halt()
    end if
    
    ! Get a pointer to the FEM data
    p_RfemData => rmatrixAssembly%rassemblyDataTemplate%rfemDataBlocks%p_RfemData
    ncount = rmatrixAssembly%rassemblyDataTemplate%rfemDataBlocks%ncount

    ! Add the evaluation tags of all FEM spaces to a common evaluation tag.
    do i=1,ncount
      rmatrixAssembly%cevaluationTag = ior(rmatrixAssembly%cevaluationTag,&
          elem_getEvaluationTag(p_RfemData(i)%celement))
    end do

    ! Check that all elements have the same NVE. Otherwise,
    ! something is wrong.
    cshape = elem_igetShape(p_RfemData(1)%celement)

    do i=2,ncount
      if (cshape .ne. elem_igetShape(p_RfemData(i)%celement)) then
        call output_line ("Element spaces incompatible!",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_initMatAssembly")
        call sys_halt()
      end if
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_doneMatAssembly(rmatrixAssembly)

!<description>
  ! Clean up a matrix assembly structure.
!</description>

!<inputoutput>
  ! Matrix assembly structure to clean up
  type(t_bmaMatrixAssembly), intent(inout) :: rmatrixAssembly
!</inputoutput>

!</subroutine>

    ! Release all allocated memory.
    call bma_doneMatAssemblyData(rmatrixAssembly%rassemblyDataTemplate)

    deallocate(rmatrixAssembly%p_DcubPtsRef)
    deallocate(rmatrixAssembly%p_Domega)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_prepareLocalMatrices(rassemblyData,IelementList)

!<description>
  ! Auxiliary subroutine. 
  ! Fetches numbers of the the local degrees of freedon of all FEM spaces.
  ! Calculates the positions of the entries in the local matrices.
!</description>

!<input>
  ! List of elements where to evaluate the FEM basis functions.
  integer, dimension(:), intent(in), target :: IelementList
!</input>

!<inputoutput>
  ! Assembly data structure. The positions of the local matrices in the
  ! global matrix are computed, and the local matrices are cleared.
  type(t_bmaMatrixAssemblyData), intent(inout), target :: rassemblyData
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j,k
    type(t_bmaMatrixData), pointer :: p_rmatrixData

    ! Remember the element list
    rassemblyData%p_IelementList => IelementList
    rassemblyData%nelements = size(IelementList)

    ! Calculate the DOF mapping for the FEM spaces.
    call fev2_calcDofMapping(rassemblyData%rfemDataBlocks,IelementList)

    ! Loop through the matrices. Whereever there is unshared data,
    ! compute the local matrix indices.
    do j=1,ubound(rassemblyData%p_RmatrixData,2)
      do i=1,ubound(rassemblyData%p_RmatrixData,1)

        p_rmatrixData => rassemblyData%p_RmatrixData(i,j)

        ! Unshared data in that block
        if (p_rmatrixData%bhasData) then 

          if (.not. p_rmatrixData%bsharedMatrixData) then

            ! For the assembly of the global matrix, we use a "local"
            ! approach. At first we build a "local" system matrix according
            ! to the current element. This contains all additive
            ! contributions of element iel, which are later added at the
            ! right positions to the elements in the global system matrix.
            !
            ! We have ndofTrial trial DOF`s per element and
            ! ndofTest test DOF`s per element. Therefore there are
            ! ndofTrial*ndofTest tupel of basis-/testfunctions (phi_i,psi_j)
            ! "active" (i.e. have common support) on our current element, each
            ! giving an additive contribution to the system matrix.
            !
            ! We build a quadratic ndofTrial*ndofTest local matrix:
            ! Kentry(1..ndofTrial,1..ndofTest) receives the position
            ! in the global system matrix, where the corresponding value
            ! has to be added to.
            ! (The corresponding contributions can be saved separately,
            ! but we directly add them to the global matrix in this
            ! approach.)
            !
            ! We build local matrices for all our elements
            ! in the set simultaneously. Get the positions of the local matrices
            ! in the global matrix.
            call bma_getLocalMatrixIndices (p_rmatrixData%p_rmatrix,&
                p_rmatrixData%p_IdofsTest,p_rmatrixData%p_IdofsTrial,p_rmatrixData%p_Kentry,&
                ubound(p_rmatrixData%p_IdofsTest,1),ubound(p_rmatrixData%p_IdofsTrial,1),&
                size(IelementList))

            ! Fill Dentry with zero; default values.
            if (.not. p_rmatrixData%bisInterleaved) then
              call lalg_clearVector(p_rmatrixData%p_Dentry)
            else
              do k=1,ubound(p_rmatrixData%p_DentryIntl,4)
                call lalg_clearVector(p_rmatrixData%p_DentryIntl(:,:,:,k))
              end do
            end if

          end if

        end if

      end do
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_evaluateFEMforMat(rassemblyData,rmatrixAssembly)

!<description>
  ! Auxiliary subroutine. 
  ! Evaluates the FEM basis functions in the cubature points.
  ! Updates the element sets according to the coordinates of the
  ! cubature points, Jacobian determinants, etc.
!</description>

!<input>
  ! A matrix assembly structure.
  type(t_bmaMatrixAssembly), intent(in), target :: rmatrixAssembly
!</input>

!<inputoutput>
  ! Assembly data structure receiving the FEM data
  type(t_bmaMatrixAssemblyData), intent(inout) :: rassemblyData
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j
    integer(I32) :: cevaluationTag
    real(DP), dimension(:,:), pointer :: p_Ddetj,p_DcubWeight
    real(DP), dimension(:), pointer :: p_Domega

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag.
    cevaluationTag = rmatrixAssembly%cevaluationTag

    ! In the first loop, calculate the coordinates on the reference element.
    ! In all later loops, use the precalculated information.
    !
    ! If the cubature points are already initialised, do not do it again.
    ! We check this by taking a look to ninitialisedElements which
    ! gives the current maximum of initialised elements.
    if (rassemblyData%nelements .gt. &
        rassemblyData%ninitialisedElements) then

      ! (Re-)initialise!
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)

      ! Remember the new number of initialised elements
      rassemblyData%ninitialisedElements = rassemblyData%nelements

      ! (Re-)allocate the cubature weights      
      if (associated(rassemblyData%p_DcubWeight)) then
        deallocate(rassemblyData%p_DcubWeight)
      end if

      allocate(rassemblyData%p_DcubWeight(&
          rmatrixAssembly%ncubp,rassemblyData%nelements))

    else
      ! No need.
      cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
    end if

    ! Calculate all information that is necessary to evaluate the finite element
    ! on all cells of our subset. This includes the coordinates of the points
    ! on the cells.
    call elprep_prepareSetForEvaluation (rassemblyData%revalElementSet,&
        cevaluationTag, rmatrixAssembly%p_rtriangulation, &
        rassemblyData%p_IelementList, rmatrixAssembly%ctrafoType, &
        rmatrixAssembly%p_DcubPtsRef(:,1:rmatrixAssembly%ncubp), &
        rperfconfig=rmatrixAssembly%p_rperfconfig)

    ! Calculate the cubature weights. THese weights must be
    ! multiplied to the values of the trial/test-functions
    ! for cubature. They calculate from the actual cubature
    ! weights and the corresponding Jacobian determinant.

    p_Ddetj => rassemblyData%revalElementSet%p_Ddetj
    p_Domega => rmatrixAssembly%p_Domega
    p_DcubWeight => rassemblyData%p_DcubWeight

    do j=1,rassemblyData%nelements
      do i=1,rmatrixAssembly%ncubp

        ! Take the absolut value of the determinant of the mapping.
        ! In 2D, the determinant is always positive, whereas in 3D,
        ! the determinant might be negative -- that is normal!

        p_DcubWeight(i,j) = abs(p_Ddetj(i,j))*p_Domega(i)

      end do
    end do

    ! Now, calculate the FEM basis functions in the cubature points
    call fev2_evaluateFemData(rassemblyData%rfemDataBlocks,&
        rassemblyData%revalElementSet)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_incorporateMatToGlobal(rassemblyData,rmatrixAssembly)

!<description>
  ! Auxiliary subroutine. 
  ! Incorporates the entries from the local matrices into the
  ! global matrix.
!</description>

!<input>
  ! Matrix assembly data structure with the local matrices.
  type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData
!</input>

!<inputoutput>
  ! A matrix assembly structure. Contains the global matrix
  ! where the data from the local matrices is incorporated to.
  type(t_bmaMatrixAssembly), intent(inout), target :: rmatrixAssembly
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j,k,l,m,ivar,nvar
    type(t_bmaMatrixData), pointer :: p_rmatrixData
    real(DP), dimension(:), pointer :: p_Da
    real(DP), dimension(:,:,:), pointer :: p_Dentry
    real(DP), dimension(:,:,:,:), pointer :: p_DentryIntl
    integer, dimension(:,:,:), pointer :: p_Kentry

    ! Loop through the vectors to incorporate the data.
    do j=1,ubound(rassemblyData%p_RmatrixData,2)
      do i=1,ubound(rassemblyData%p_RmatrixData,1)

        p_rmatrixData => rassemblyData%p_RmatrixData(i,j)

        if ((p_rmatrixData%bhasData) .and. (.not. p_rmatrixData%bsharedMatrixData)) then 

          ! Unshared data in that block.
          ! Incorporate the local matrices into the global one.

          p_Da => p_rmatrixData%p_Da
          p_Kentry => p_rmatrixData%p_Kentry

          if (.not. p_rmatrixData%bisInterleaved) then

            ! Non-interleaved matrix
            p_Dentry => p_rmatrixData%p_Dentry

            ! The outer-most loop loops only up to the maximum element
            ! number. The allocated memory may be larger.
            do m=1,rassemblyData%nelements
              do l=1,ubound(p_Dentry,2)
                do k=1,ubound(p_Dentry,1)
                  p_Da(p_Kentry(k,l,m)) = p_Da(p_Kentry(k,l,m)) + p_Dentry(k,l,m)
                end do
              end do
            end do

          else

            ! Interleaved matrix

            p_DentryIntl => p_rmatrixData%p_DentryIntl
            nvar = p_rmatrixData%nvar

            ! The outer-most loop loops only up to the maximum element
            ! number. The allocated memory may be larger.
            do m=1,rassemblyData%nelements
              do l=1,ubound(p_DentryIntl,3)
                do k=1,ubound(p_DentryIntl,2)
                  do ivar=1,ubound(p_DentryIntl,1)
                    p_Da(nvar*(p_Kentry(k,l,m)-1)+ivar) = &
                        p_Da(nvar*(p_Kentry(k,l,m)-1)+ivar) + p_DentryIntl(ivar,k,l,m)
                  end do
                end do
              end do
            end do

          end if

        end if

      end do
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_assembleSubmeshMatrix (rmatrixAssembly, IelementList,&
      fcalcLocalMatrices, rcollection)

!<description>
  ! Assembles the matrix entries for a list of elements by integrating
  ! over the domain.
!</description>

!<input>

  ! List of elements where to assemble the bilinear form.
  integer, dimension(:), intent(in), target :: IelementList

  interface

    subroutine fcalcLocalMatrices (RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

      use collection
      use blockmatassemblybase

      ! Calculates the local matrices Dentry in all RmatrixData structures.

      ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
      ! have to be filled with data.
      type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData

      ! Data necessary for the assembly. Contains determinants and
      ! cubature weights for the cubature,...
      type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

      ! Structure with all data about the assembly
      type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

      ! Number of points per element
      integer, intent(in) :: npointsPerElement

      ! Number of elements
      integer, intent(in) :: nelements

      ! Values of FEM functions automatically evaluated in the
      ! cubature points.
      type(t_fev2Vectors), intent(in) :: revalVectors

      ! User defined collection structure
      type(t_collection), intent(inout), target, optional :: rcollection

    end subroutine

  end interface  

!</input>

!<inputoutput>

  ! A matrix assembly structure prepared with bma_initMatAssembly.
  type(t_bmaMatrixAssembly), intent(inout), target :: rmatrixAssembly

  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information.
  type(t_collection), intent(inout), target, optional :: rcollection

!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ielStart,ielMax
    type(t_bmaMatrixAssemblyData) :: rassemblyData
    type(t_fev2Vectors) :: revalVectors

    ! Use the template structure to create a matrix data array.
    ! Allocate memory for the assembly.
    ! NOTE: This can be done in parallel, the information in rassemblyData
    ! is independent!
    call bma_createMatAssemblyData(rmatrixAssembly,&
        rassemblyData,rmatrixAssembly%rassemblyDataTemplate,&
        revalVectors,rmatrixAssembly%revalVectorsTemplate)

    ! Loop blockwise through the element list
    do ielStart = 1,size(IelementList),rmatrixAssembly%nelementsPerBlock

      ! End of the current block
      ielMax = min(ielStart-1+rmatrixAssembly%nelementsPerBlock,size(IelementList))

      ! Calculate the indices of the matrix entries to be modified,
      ! i.e., set up the local matrices.
      call bma_prepareLocalMatrices(rassemblyData,IelementList(ielStart:ielMax))

      ! Calculate the FEM basis functions in the cubature points.
      call bma_evaluateFEMforMat(rassemblyData,rmatrixAssembly)

      ! (Re-)allocate memory for the FEM evaluation if necessary.
      call fev2_prepareVectorEval(revalVectors,rassemblyData%revalElementSet)

      ! Evaluate the attached vectors in the cubature points.
      call fev2_evaluateVectors(revalVectors,rassemblyData%rfemDataBlocks)

      ! Use the callback routine to calculate the local matrix entries.
      call fcalcLocalMatrices(rassemblyData%p_RmatrixData,rassemblyData,rmatrixAssembly,&
          rmatrixAssembly%ncubp,ielMax-ielStart+1,revalVectors,rcollection)

      ! Incorporate the local matrices into the global one.
      ! NOTE: This cannot be done in parallel!
      call bma_incorporateMatToGlobal(rassemblyData,rmatrixAssembly)

    end do

    ! Release memory
    call bma_releaseMatAssemblyData(rassemblyData,revalVectors)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_buildMatrix (rmatrix,cflags,&
      fcalcLocalMatrices, rcollection, &
      RmaxDerivativeTest,RmaxDerivativeTrial,&
      revalVectors,rcubatureInfo,rtrafoInfo,rperfconfig)

!<description>
  ! This subroutine calculates the entries of a block matrix and
  ! adds them to rmatrix. The callback function fcalcLocalMatrices
  ! has to compute the local matrix indices by cubature.
!</description>

!<input>

  ! Option field. Combination of BMA_CALC_xxxx flags.
  ! Use BMA_CALC_STANDARD for standard options.
  integer(I32), intent(in) :: cflags

  interface

    subroutine fcalcLocalMatrices(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

      use collection
      use blockmatassemblybase

      ! Calculates the local matrices Dentry in all RmatrixData structures.

      ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
      ! have to be filled with data.
      type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData

      ! Data necessary for the assembly. Contains determinants and
      ! cubature weights for the cubature,...
      type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

      ! Structure with all data about the assembly
      type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

      ! Number of points per element
      integer, intent(in) :: npointsPerElement

      ! Number of elements
      integer, intent(in) :: nelements

      ! Values of FEM functions automatically evaluated in the
      ! cubature points.
      type(t_fev2Vectors), intent(in) :: revalVectors

      ! User defined collection structure
      type(t_collection), intent(inout), target, optional :: rcollection

    end subroutine

  end interface  

  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information.
  type(t_collection), intent(inout), target, optional :: rcollection

  ! OPTIONAL: For every block in the matrix, maximum
  ! derivative of the basis functions to be computed. If not
  ! specified or an entry is =0, the maximum available derivative for 
  ! each FEM space is the default.
  integer, dimension(:,:), intent(in), optional :: RmaxDerivativeTest
  integer, dimension(:,:), intent(in), optional :: RmaxDerivativeTrial

  ! OPTIONAL: Set of vectors to be automatically evaluated
  type(t_fev2Vectors), intent(in), optional :: revalVectors

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), target, optional :: rcubatureInfo

  ! OPTIONAL: A transformation structure that specifies the transformation
  ! from the reference to the real element(s).
  ! If not specified, default settings are used.
  type(t_scalarTrafoInfo), intent(in), target, optional :: rtrafoInfo

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The FE matrix to calculate
  type(t_matrixBlock), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_scalarCubatureInfo), target :: rtempCubatureInfo
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo
    type(t_bmaMatrixAssembly) :: rmatrixAssembly
    integer :: ielementDistr,icubatureBlock,NEL
    integer, dimension(:), pointer :: p_IelementList
    integer(I32) :: ccubType,celement,ctrafoType
    type(t_spatialDiscretisation), pointer :: p_rdiscr
    type(t_scalarTrafoInfo), pointer :: p_rtrafoInfo
    type(t_scalarTrafoInfo), target :: rtempTrafoInfo

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => bma_perfconfig
    end if

    ! Fetch the first available spatial discretisation.
    ! It defines the way the elements are separated into blocks of elements.
    p_rdiscr => rmatrix%p_rblockDiscrTest%RspatialDiscr(1)

    ! If we do not have it, create a cubature info structure that
    ! defines how to do the assembly.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_createDefCubStructure(p_rdiscr,&
          rtempCubatureInfo,CUB_GEN_DEPR_BILFORM)
      p_rcubatureInfo => rtempCubatureInfo
    else
      p_rcubatureInfo => rcubatureInfo
    end if

    ! Transformation structure?
    if (present(rtrafoInfo)) then
      p_rtrafoInfo => rtrafoInfo
    else
      ! Set up a default transformation structure.
      ! Do we have a discretisation that allows us to set up the structure?
      if (associated(p_rdiscr)) then
        ! We found a discretisation. Set up the transformation
        call spdiscr_createDefTrafoStructure (p_rdiscr, rtempTrafoInfo)
      else
        ! Fallback procedure: No discretisation found.
        ! Create a standard transformation based on the cubature formula.
        call spdiscr_createDefTrafoByCubInfo (rcubatureInfo, rtempTrafoInfo)
      end if
      
      p_rtrafoInfo => rtempTrafoInfo
    end if

    ! Loop over the cubature blocks to discretise
    do icubatureBlock = 1,p_rcubatureInfo%ninfoBlockCount

      ! Get information about that block as well as an appropriate cubature formula
      call spdiscr_getStdDiscrInfo(icubatureBlock,p_rcubatureInfo,&
          p_rdiscr,ielementDistr,celement,ccubType,NEL,p_IelementList,&
          p_rtrafoInfo,ctrafoType=ctrafoType)

      ! Check if element distribution is empty
      if (NEL .le. 0) cycle

      ! Initialise a matrix assembly structure for that element distribution
      call bma_initMatAssembly(rmatrixAssembly,ccubType,ctrafoType,cflags,&
          rmatrix,ielementDistr,&
          RmaxDerivativeTest,RmaxDerivativeTrial,revalVectors,rperfconfig)

      ! Assemble the data for all elements in this element distribution
      call bma_assembleSubmeshMatrix (rmatrixAssembly, p_IelementList(1:NEL),&
          fcalcLocalMatrices, rcollection)

      ! Release the assembly structure.
      call bma_doneMatAssembly(rmatrixAssembly)
    end do

    ! Release the assembly structure if necessary.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_releaseCubStructure(rtempCubatureInfo)
    end if

    ! Release the transformation structure if necessary.
    if (.not. present(rtrafoInfo)) then
      call spdiscr_releaseTrafoStructure (rtempTrafoInfo)
    end if

  end subroutine

  !****************************************************************************
  !****************************************************************************
  ! Assembly of block vectors
  !****************************************************************************
  !****************************************************************************

!<subroutine>

  subroutine bma_prepareVectorData(&
      rvector,RvectorData,rfemDataBlocks,ielementDistr)

!<description>
  ! Initialise a vector assembly structure for assembling the subvectors
  ! of a block vector.
  ! No memory is allocated.
!</description>

!<input>
  ! The vector which is going to be assembled.
  type(t_vectorBlock), intent(in) :: rvector

  ! ID of the active element distribution. Identifies the FEM space
  ! in the vector which is to be assembled.
  integer, intent(in) :: ielementDistr

  ! Structure for all involved FEM spaces.
  type(t_fev2FemDataBlocks), intent(in) :: rfemDataBlocks
!</input>

!<output>
  ! An array of block assembly structures for all the subvectors
  ! in the block vector.
  type(t_bmaVectorData), dimension(:), intent(out) :: RvectorData
!</output>

!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscrTest
    integer :: i

    ! Allocate memory for existence checks
    i = size(rvector%RvectorBlock)

    ! Basic initialisation of the structure.
    ! Loop over all blocks. Figure out which blocks have data.
    ! Get the data arrays for these blocks.
    do i=1,ubound(RvectorData,1)

      ! Remember the vector
      RvectorData(i)%p_rvector => rvector%RvectorBlock(i)

      ! Get a pointer to the vector data for faster access
      call lsyssc_getbase_double (rvector%RvectorBlock(i),RvectorData(i)%p_Ddata)

      ! Check if we assemble an interleaved vector
      RvectorData(i)%nvar = rvector%RvectorBlock(i)%nvar
      RvectorData(i)%bisInterleaved = RvectorData(i)%nvar .gt. 1

      ! Get the information about the trial and test spaces
      ! of this block.
      p_rdiscrTest => rvector%RvectorBlock(i)%p_rspatialDiscr

      ! Get the indices of the trial and the test space.
      RvectorData(i)%iidxFemDataTest = containsDiscr (rfemDataBlocks,p_rdiscrTest)

      ! Element type          
      RvectorData(i)%celementTest = &
          p_rdiscrTest%RelementDistr(ielementDistr)%celement

      ! Get the number of local DOF`s for trial and test functions
      RvectorData(i)%ndofTest = elem_igetNDofLoc(RvectorData(i)%celementTest)

    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_initVecAssemblyData(rassemblyData,&
      rvector,ielementDistr,RmaxDerivativeTest,revalVectors)

!<description>
  ! Initialises the main parameters in the structure containing the
  ! assembly data (cubature weights, etc.).
  !
  ! The returned rassemblyData structure can be used as template
  ! for bma_createVecAssemblyData.
!</description>

!<input>
  ! The vector which is going to be assembled.
  type(t_vectorBlock), intent(in) :: rvector

  ! ID of the active element distribution. Identifies the FEM space
  ! in the vector which is to be assembled.
  integer, intent(in) :: ielementDistr

  ! OPTIONAL: For every block in the vector, maximum
  ! derivative of the basis functions to be computed. If not
  ! specified or an entry is =-1, only the function values for 
  ! each FEM space are computed, derivatives are omitted.
  integer, dimension(:), intent(in), optional :: RmaxDerivativeTest
  
  ! OPTIONAL: A set of FEM functions to be evaluated during the vector
  ! assembly.
  type(t_fev2Vectors), intent(in), optional :: revalVectors
!</input>

!<output>
  ! The vector assembly data structure to be initialised.
  type(t_bmaVectorAssemblyData), intent(out) :: rassemblyData
!</output>

!</subroutine>

    ! Initialise the element set; can be simultaneously used for all
    ! FEM spaces, it is element independent.
    call elprep_init(rassemblyData%revalElementSet)

    ! Element sets are not yet initialised.  
    rassemblyData%ninitialisedElements = 0

    ! Prepare FEM data
    call fev2_prepareFemDataBVec(rvector,rassemblyData%rfemDataBlocks,&
        ielementDistr,RmaxDerivativeTest)
        
    ! Prepare FEM data for the vectors to be evaluated.
    if (present(revalVectors)) then
      call fev2_prepareFemDataVecEval(revalVectors,rassemblyData%rfemDataBlocks,&
          ielementDistr)
    end if

    ! Initialise the blocks according to the structure of the vector
    ! and the FEM spaces used.
    allocate(rassemblyData%p_RvectorData(ubound(rvector%RvectorBlock,1)))

    call bma_prepareVectorData(&
        rvector,rassemblyData%p_RvectorData,rassemblyData%rfemDataBlocks,ielementDistr)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_doneVecAssemblyData(rassemblyData)

!<description>
  ! Releases allocated memory from the assembly data structure.
!</description>

!<inputoutput>
  ! Assembly data structure to be cleaned up.
  type(t_bmaVectorAssemblyData), intent(inout) :: rassemblyData
!</inputoutput>

!</subroutine>

    ! Release the vector and FEM data array
    call fev2_cleanupFemData (rassemblyData%rfemDataBlocks)
    deallocate(rassemblyData%p_RvectorData)

    ! Element sets are not yet initialised.  
    rassemblyData%ninitialisedElements = 0

  end subroutine


  !****************************************************************************

!<subroutine>

  subroutine bma_createVecAssemblyData(rvectorAssembly,&
      rassemblyData,rassemblyDataTemplate,&
      revalVectors,revalVectorsTemplate)

!<description>
  ! Based on a template structure, this creates an assembly data structure
  ! for the vector assembly with all data arrays initialised and allocated.
!</description>

!<input>
  ! Template assembly data structure. rassemblyData is created based
  ! on the data in this structure.
  type(t_bmaVectorAssemblyData), intent(in) :: rassemblyDataTemplate

  ! A vector assembly structure.
  type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly

  ! Template vector data structure. rvectorData is created based
  ! on the data in this structure.
  type(t_fev2Vectors), intent(in) :: revalVectorsTemplate
!</input>

!<output>
  ! The vector assembly data structure to be initialised.
  type(t_bmaVectorAssemblyData), intent(out) :: rassemblyData

  ! Vector data structure to be initialised according to revalVectorsTemplate
  type(t_fev2Vectors), intent(out) :: revalVectors
!</output>

!</subroutine>

    ! local variables
    integer :: i
    type(t_bmaVectorData), pointer :: p_rvectorData

    ! At first, just copy the content from the vector template structure
    rassemblyData = rassemblyDataTemplate

    ! Now initialise the dynamic content.
    !
    ! Copy the FEM data and the vector data
    allocate (rassemblyData%p_RvectorData(&
        ubound(rassemblyDataTemplate%p_RvectorData,1)))
    rassemblyData%p_RvectorData(:) = rassemblyDataTemplate%p_RvectorData(:)

    call fev2_copyFemData(rassemblyData%rfemDataBlocks,rassemblyDataTemplate%rfemDataBlocks)
    
    ! Initialise the FEM evaluation structures.
    call fev2_createFemData(rassemblyData%rfemDataBlocks,&
        rvectorAssembly%ncubp,rvectorAssembly%nelementsPerBlock)

    ! Loop through the subvectors. Aallocate memory for the vector
    ! entries in this block.
    do i=1,ubound(rassemblyData%p_RvectorData,1)

      p_rvectorData => rassemblyData%p_RvectorData(i)

      ! Allocate memory for the vector data and positions in the vector
      if (.not. p_rvectorData%bisInterleaved) then
        ! Non-interleaved vector
        allocate(p_rvectorData%p_Dentry(&
                p_rvectorData%ndofTest,rvectorAssembly%nelementsPerBlock))
      else
        ! Interleaved vector
        allocate(p_rvectorData%p_DentryIntl(p_rvectorData%nvar,&
                p_rvectorData%ndofTest,rvectorAssembly%nelementsPerBlock))
      end if

      ! Map the pointers to the DOFs and the values of the
      ! basis functions into that block
      p_rvectorData%p_DbasTest => &
          rassemblyData%rfemDataBlocks%p_RfemData(p_rvectorData%iidxFemDataTest)%p_Dbas

      p_rvectorData%p_IdofsTest => &
          rassemblyData%rfemDataBlocks%p_RfemData(p_rvectorData%iidxFemDataTest)%p_Idofs

    end do

    ! Initialise the vector evaluation structure.
    !
    ! Copy the content from the template.
    revalVectors = revalVectorsTemplate

    ! Re-create the vector array
    if (revalVectors%ncount .eq. 0) then
      nullify(revalVectors%p_RvectorData)
    else
      allocate(revalVectors%p_RvectorData(revalVectorsTemplate%ncount))
      revalVectors%p_RvectorData(1:revalVectorsTemplate%ncount) = &
          revalVectorsTemplate%p_RvectorData(1:revalVectorsTemplate%ncount)

      ! Prepare the evaluation of the FEM functions
      call fev2_initVectorEval(revalVectors,rassemblyData%rfemDataBlocks)
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_releaseVecAssemblyData(rassemblyData,revalVectors)

!<description>
  ! Releases memory of an assembly structure created by bma_createVecAssemblyData.
!</description>

!<inputoutput>
  ! The vector assembly data structure to be cleaned up
  type(t_bmaVectorAssemblyData), intent(inout) :: rassemblyData

  ! Vector data structure to be cleaned up
  type(t_fev2Vectors), intent(inout) :: revalVectors
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    type(t_bmaVectorData), pointer :: p_rvectorData

    ! Release vector evaluzation data
    if (revalVectors%ncount .ne. 0) then
      call fev2_doneVectorEval(revalVectors)
      deallocate(revalVectors%p_RvectorData)
    end if

    ! Release FEM evaluation data
    call fev2_releaseFemData(rassemblyData%rfemDataBlocks)

    ! Release the element set and the cubature weights
    call elprep_releaseElementSet(rassemblyData%revalElementSet)

    if (associated(rassemblyData%p_DcubWeight)) then
      deallocate(rassemblyData%p_DcubWeight)
    end if

    ! Element sets are not initialised anymore
    rassemblyData%ninitialisedElements = 0

    ! Loop through the subvectors and release memory
    do i=1,ubound(rassemblyData%p_rvectorData,1)

      p_rvectorData => rassemblyData%p_RvectorData(i)

      ! Release memory of the local vector  entries
      if (associated(p_rvectorData%p_Dentry)) then
        deallocate(p_rvectorData%p_Dentry)
      end if

      if (associated(p_rvectorData%p_DentryIntl)) then
        deallocate(p_rvectorData%p_DentryIntl)
      end if

    end do

    ! Release memory
    call fev2_cleanupFemData (rassemblyData%rfemDataBlocks)
    deallocate (rassemblyData%p_rvectorData)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_initVecAssembly(rvectorAssembly,ccubType,ctrafoType,cflags,&
      rvector,ielementDistr,RmaxDerivativeTest,&
      revalVectors,rperfconfig)

!<description>
  ! Initialise a vector assembly structure for assembling a linear form
  ! for a certain element distribution.
!</description>

!<input>
  ! Cubature formla to use
  integer(I32), intent(in) :: ccubType

  ! Transformation formla to use
  integer(I32), intent(in) :: ctrafoType

  ! Option field. Combination of BMA_CALC_xxxx flags.
  ! Use BMA_CALC_STANDARD for standard options.
  integer(I32), intent(in) :: cflags

  ! The vector which is going to be assembled.
  type(t_vectorBlock), intent(in) :: rvector

  ! ID of the active element distribution. Identifies the FEM space
  ! in the vector which is to be assembled.
  integer, intent(in) :: ielementDistr

  ! OPTIONAL: For every block in the vector, maximum
  ! derivative of the basis functions to be computed. If not
  ! specified or an entry is =-1, only the function values for 
  ! each FEM space are computed, derivatives are omitted.
  integer, dimension(:), intent(in), optional :: RmaxDerivativeTest

  ! OPTIONAL: Set of vectors to be automatically evaluated
  type(t_fev2Vectors), intent(in), optional :: revalVectors

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<output>
  ! A vector assembly structure.
  type(t_bmaVectorAssembly), intent(out) :: rvectorAssembly
!</output>

!</subroutine>

    integer :: i,ncount
    integer(I32) :: cshape
    type(t_fev2FemData), dimension(:), pointer :: p_RfemData

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Remember the performance configuration for later use
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => bma_perfconfig
    end if

    rvectorAssembly%p_rperfconfig => p_rperfconfig

    ! Get basic data
    rvectorAssembly%p_rboundary => rvector%p_rblockDiscr%p_rboundary
    rvectorAssembly%p_rtriangulation => rvector%p_rblockDiscr%p_rtriangulation

    ! Save the vector evaluation structure as template.
    ! Make a copy of the stucture, that is simpler here.
    if (present(revalVectors)) then
      rvectorAssembly%revalVectorsTemplate = revalVectors
    end if

    ! Remember current element distribution
    rvectorAssembly%ielementDistr = ielementDistr

    ! Initialise the template structure which is used during the
    ! actual vector assembly.
    call bma_initVecAssemblyData(rvectorAssembly%rassemblyDataTemplate,&
        rvector,ielementDistr,RmaxDerivativeTest,&
        revalVectors)

    ! Get the cubature/transformation formula on/from the reference element
    rvectorAssembly%ccubType = ccubType
    rvectorAssembly%ncubp = cub_igetNumPts(ccubType)
    rvectorAssembly%ctrafoType = ctrafoType

    ! Allocate two arrays for the points and the weights
    allocate(rvectorAssembly%p_Domega(rvectorAssembly%ncubp))
    allocate(rvectorAssembly%p_DcubPtsRef(&
        trafo_igetReferenceDimension(rvectorAssembly%ctrafoType),&
        rvectorAssembly%ncubp))

    ! Get the cubature formula
    call cub_getCubature(ccubType,rvectorAssembly%p_DcubPtsRef,rvectorAssembly%p_Domega)

    ! Number of simultaneously processed elements
    rvectorAssembly%nelementsPerBlock = p_rperfconfig%NELEMSIM

    ! Get the evaluation tag by "OR"-ing all those from the FEM spaces
    ! and the option fields.
    rvectorAssembly%cevaluationTag = 0

    if (iand(cflags,BMA_CALC_REALCOORDS) .ne. 0) then
      ! Calculate real world coordinates of the cubature points
      rvectorAssembly%cevaluationTag = &
          ior(rvectorAssembly%cevaluationTag,EL_EVLTAG_REALPOINTS)
    end if

    if (rvectorAssembly%rassemblyDataTemplate%rfemDataBlocks%ncount .eq. 0) then
      call output_line ("Vector is empty, nothing to assemble.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_initVecAssembly")
      call sys_halt()
    end if

    ! Get a pointer to the FEM data
    p_RfemData => rvectorAssembly%rassemblyDataTemplate%rfemDataBlocks%p_RfemData
    ncount = rvectorAssembly%rassemblyDataTemplate%rfemDataBlocks%ncount

    ! Add the evaluation tags of all FEM spaces to a common evaluation tag.
    do i=1,ncount
      rvectorAssembly%cevaluationTag = ior(rvectorAssembly%cevaluationTag,&
          elem_getEvaluationTag(p_RfemData(i)%celement))
    end do

    ! Check that all elements have the same NVE. Otherwise,
    ! something is wrong.
    cshape = elem_igetShape(p_RfemData(1)%celement)

    do i=2,ncount
      if (cshape .ne. elem_igetShape(p_RfemData(i)%celement)) then
        call output_line ("Element spaces incompatible!",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_initVecAssembly")
        call sys_halt()
      end if
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_doneVecAssembly(rvectorAssembly)

!<description>
  ! Clean up a vector assembly structure.
!</description>

!<inputoutput>
  ! Vector assembly structure to clean up
  type(t_bmaVectorAssembly), intent(inout) :: rvectorAssembly
!</inputoutput>

!</subroutine>

    ! Release all allocated memory.
    call bma_doneVecAssemblyData(rvectorAssembly%rassemblyDataTemplate)

    deallocate(rvectorAssembly%p_DcubPtsRef)
    deallocate(rvectorAssembly%p_Domega)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_prepareLocalVectors(rassemblyData,IelementList)

!<description>
  ! Auxiliary subroutine. 
  ! Initialises the data arrays which hold the vector data.
!</description>

!<input>
  ! List of elements where to evaluate the FEM basis functions.
  integer, dimension(:), intent(in), target :: IelementList
!</input>

!<inputoutput>
  ! Assembly data structure. The DOFs in to be modified in the
  ! block vector are computed, and the local vectors are cleared.
  type(t_bmaVectorAssemblyData), intent(inout), target :: rassemblyData
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    type(t_bmaVectorData), pointer :: p_rvectorData

    ! Remember the element list
    rassemblyData%p_IelementList => IelementList
    rassemblyData%nelements = size(IelementList)

    ! Calculate the DOF mapping for the FEM spaces.
    call fev2_calcDofMapping(rassemblyData%rfemDataBlocks,IelementList)

    ! Loop through the subvectors. Whereever there is unshared data,
    ! Clear the destination arrays.
    do i=1,ubound(rassemblyData%p_RvectorData,1)

      p_rvectorData => rassemblyData%p_RvectorData(i)

      ! Fill Dentry with zero; default values.
      if (.not. p_rvectorData%bisInterleaved) then
        ! Non-Interleaved
        call lalg_clearVector(p_rvectorData%p_Dentry)
      else
        ! Interleaved
        call lalg_clearVector(p_rvectorData%p_DentryIntl)
      end if

    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_evaluateFEMforVec(rassemblyData,rvectorAssembly)

!<description>
  ! Auxiliary subroutine. 
  ! Evaluates the FEM basis functions in the cubature points.
  ! Updates the element sets according to the coordinates of the
  ! cubature points, Jacobian determinants, etc.
!</description>

!<input>
  ! A vector assembly structure.
  type(t_bmaVectorAssembly), intent(in), target :: rvectorAssembly
!</input>

!<inputoutput>
  ! Assembly data structure receiving the FEM data
  type(t_bmaVectorAssemblyData), intent(inout) :: rassemblyData
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j
    integer(I32) :: cevaluationTag
    real(DP), dimension(:,:), pointer :: p_Ddetj,p_DcubWeight
    real(DP), dimension(:), pointer :: p_Domega

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag.
    cevaluationTag = rvectorAssembly%cevaluationTag

    ! In the first loop, calculate the coordinates on the reference element.
    ! In all later loops, use the precalculated information.
    !
    ! If the cubature points are already initialised, do not do it again.
    ! We check this by taking a look to ninitialisedElements which
    ! gives the current maximum of initialised elements.
    if (rassemblyData%nelements .gt. &
        rassemblyData%ninitialisedElements) then

      ! (Re-)initialise!
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)

      ! Remember the new number of initialised elements
      rassemblyData%ninitialisedElements = rassemblyData%nelements

      ! (Re-)allocate the cubature weights      
      if (associated(rassemblyData%p_DcubWeight)) then
        deallocate(rassemblyData%p_DcubWeight)
      end if

      allocate(rassemblyData%p_DcubWeight(&
          rvectorAssembly%ncubp,rassemblyData%nelements))

    else
      ! No need.
      cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
    end if

    ! Calculate all information that is necessary to evaluate the finite element
    ! on all cells of our subset. This includes the coordinates of the points
    ! on the cells.
    call elprep_prepareSetForEvaluation (rassemblyData%revalElementSet,&
        cevaluationTag, rvectorAssembly%p_rtriangulation, &
        rassemblyData%p_IelementList, rvectorAssembly%ctrafoType, &
        rvectorAssembly%p_DcubPtsRef(:,1:rvectorAssembly%ncubp), &
        rperfconfig=rvectorAssembly%p_rperfconfig)

    ! Calculate the cubature weights. THese weights must be
    ! multiplied to the values of the trial/test-functions
    ! for cubature. They calculate from the actual cubature
    ! weights and the corresponding Jacobian determinant.

    p_Ddetj => rassemblyData%revalElementSet%p_Ddetj
    p_Domega => rvectorAssembly%p_Domega
    p_DcubWeight => rassemblyData%p_DcubWeight

    do j=1,rassemblyData%nelements
      do i=1,rvectorAssembly%ncubp

        ! Take the absolut value of the determinant of the mapping.
        ! In 2D, the determinant is always positive, whereas in 3D,
        ! the determinant might be negative -- that is normal!

        p_DcubWeight(i,j) = abs(p_Ddetj(i,j))*p_Domega(i)

      end do
    end do

    ! Now, calculate the FEM basis functions in the cubature points
    call fev2_evaluateFemData(rassemblyData%rfemDataBlocks,&
        rassemblyData%revalElementSet)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_incorporateVecToGlobal(rassemblyData,rvectorAssembly)

!<description>
  ! Auxiliary subroutine. 
  ! Incorporates the entries from the local vectors into the
  ! global vector.
!</description>

!<input>
  ! Vector assembly data structure with the local vectors.
  type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData
!</input>

!<inputoutput>
  ! A vector assembly structure. Contains the global vector
  ! where the data from the local vectors is incorporated to.
  type(t_bmaVectorAssembly), intent(inout), target :: rvectorAssembly
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,k,l,ivar,nvar
    type(t_bmaVectorData), pointer :: p_rvectorData
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(:,:), pointer :: p_Dentry
    real(DP), dimension(:,:,:), pointer :: p_DentryIntl
    integer, dimension(:,:), pointer :: p_IdofsTest

    ! Loop through the vectors to incorporate the data.
    do i=1,ubound(rassemblyData%p_rvectorData,1)

      p_rvectorData => rassemblyData%p_rvectorData(i)
      p_IdofsTest => p_rvectorData%p_IdofsTest

      ! Incorporate the local vectors into the global one.

      p_Ddata => p_rvectorData%p_Ddata

      if (.not. p_rvectorData%bisInterleaved) then

        ! Non-Interleaved
        p_Dentry => p_rvectorData%p_Dentry

        ! The outer-most loop loops only up to the maximum element
        ! number. The allocated memory may be larger.
        do l=1,rassemblyData%nelements
          do k=1,ubound(p_Dentry,1)
            p_Ddata(p_IdofsTest(k,l)) = p_Ddata(p_IdofsTest(k,l)) + p_Dentry(k,l)
          end do
        end do

      else

        ! Interleaved
        p_DentryIntl => p_rvectorData%p_DentryIntl
        nvar = p_rvectorData%nvar

        ! The outer-most loop loops only up to the maximum element
        ! number. The allocated memory may be larger.
        do l=1,rassemblyData%nelements
          do k=1,ubound(p_DentryIntl,2)
            do ivar = 1,ubound(p_DentryIntl,1)
              p_Ddata(nvar*(p_IdofsTest(k,l)-1)+ivar) = &
                  p_Ddata(nvar*(p_IdofsTest(k,l)-1)+ivar) + p_DentryIntl(ivar,k,l)
            end do
          end do
        end do

      end if

    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_assembleSubmeshVector (rvectorAssembly, IelementList,&
      fcalcLocalVectors, rcollection)

!<description>
  ! Assembles the vector entries for a list of elements by integrating
  ! over the domain.
!</description>

!<input>

  ! List of elements where to assemble the bilinear form.
  integer, dimension(:), intent(in), target :: IelementList

  interface

    subroutine fcalcLocalVectors (rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

      use collection
      use blockmatassemblybase

      ! Calculates the local vectors Dentry in all rvectorData structures.

      ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
      ! have to be filled with data.
      type(t_bmaVectorData), dimension(:), intent(inout), target :: rvectorData

      ! Data necessary for the assembly. Contains determinants and
      ! cubature weights for the cubature,...
      type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

      ! Structure with all data about the assembly
      type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly

      ! Number of points per element
      integer, intent(in) :: npointsPerElement

      ! Number of elements
      integer, intent(in) :: nelements

      ! Values of FEM functions automatically evaluated in the
      ! cubature points.
      type(t_fev2Vectors), intent(in) :: revalVectors

      ! User defined collection structure
      type(t_collection), intent(inout), target, optional :: rcollection

    end subroutine

  end interface  

!</input>

!<inputoutput>

  ! A vector assembly structure prepared with bma_initMatAssembly.
  type(t_bmaVectorAssembly), intent(inout), target :: rvectorAssembly

  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information.
  type(t_collection), intent(inout), target, optional :: rcollection

!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ielStart,ielMax
    type(t_bmaVectorAssemblyData) :: rassemblyData
    type(t_fev2Vectors) :: revalVectors

    ! Use the template structure to create a vector data array.
    ! Allocate memory for the assembly.
    ! NOTE: This can be done in parallel, the information in rassemblyData
    ! is independent!
    call bma_createVecAssemblyData(rvectorAssembly,&
        rassemblyData,rvectorAssembly%rassemblyDataTemplate,&
        revalVectors,rvectorAssembly%revalVectorsTemplate)

    ! Loop blockwise through the element list
    do ielStart = 1,size(IelementList),rvectorAssembly%nelementsPerBlock

      ! End of the current block
      ielMax = min(ielStart-1+rvectorAssembly%nelementsPerBlock,size(IelementList))

      ! Calculate the indices of the vector entries to be modified,
      ! i.e., set up the local vectors.
      call bma_prepareLocalVectors(rassemblyData,IelementList(ielStart:ielMax))

      ! Calculate the FEM basis functions in the cubature points.
      call bma_evaluateFEMforVec(rassemblyData,rvectorAssembly)

      ! (Re-)allocate memory for the FEM evaluation if necessary.
      call fev2_prepareVectorEval(revalVectors,rassemblyData%revalElementSet)

      ! Evaluate the attached vectors in the cubature points.
      call fev2_evaluateVectors(revalVectors,rassemblyData%rfemDataBlocks)

      ! Use the callback routine to calculate the local vector entries.
      call fcalcLocalVectors(rassemblyData%p_rvectorData,rassemblyData,rvectorAssembly,&
          rvectorAssembly%ncubp,ielMax-ielStart+1,revalVectors,rcollection)

      ! Incorporate the local vectors into the global one.
      ! NOTE: This cannot be done in parallel!
      call bma_incorporateVecToGlobal(rassemblyData,rvectorAssembly)

    end do

    ! Release memory
    call bma_releaseVecAssemblyData(rassemblyData,revalVectors)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_buildVector (rvector,cflags,&
      fcalcLocalVectors, rcollection, &
      RmaxDerivativeTest,revalVectors,rcubatureInfo,rtrafoInfo,rperfconfig)

!<description>
  ! This subroutine calculates the entries of a block vector and
  ! adds them to rvector. The callback function fcalcLocalVectors
  ! has to compute the local vector indices by cubature.
!</description>

!<input>

  ! Option field. Combination of BMA_CALC_xxxx flags.
  ! Use BMA_CALC_STANDARD for standard options.
  integer(I32), intent(in) :: cflags

  interface

    subroutine fcalcLocalVectors (rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

      use collection
      use blockmatassemblybase

      ! Calculates the local vectors Dentry in all rvectorData structures.

      ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
      ! have to be filled with data.
      type(t_bmaVectorData), dimension(:), intent(inout), target :: rvectorData

      ! Data necessary for the assembly. Contains determinants and
      ! cubature weights for the cubature,...
      type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

      ! Structure with all data about the assembly
      type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly

      ! Number of points per element
      integer, intent(in) :: npointsPerElement

      ! Number of elements
      integer, intent(in) :: nelements

      ! Values of FEM functions automatically evaluated in the
      ! cubature points.
      type(t_fev2Vectors), intent(in) :: revalVectors

      ! User defined collection structure
      type(t_collection), intent(inout), target, optional :: rcollection

    end subroutine

  end interface  

  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information.
  type(t_collection), intent(inout), target, optional :: rcollection

  ! OPTIONAL: For every block in the vector, maximum
  ! derivative of the basis functions to be computed. If not
  ! specified or an entry is =0, the maximum available derivative for 
  ! each FEM space is the default.
  integer, dimension(:), intent(in), optional :: RmaxDerivativeTest

  ! OPTIONAL: Set of vectors to be automatically evaluated
  type(t_fev2Vectors), intent(in), optional :: revalVectors

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), target, optional :: rcubatureInfo

  ! OPTIONAL: A transformation structure that specifies the transformation
  ! from the reference to the real element(s).
  ! If not specified, default settings are used.
  type(t_scalarTrafoInfo), intent(in), target, optional :: rtrafoInfo

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The block vector to calculate
  type(t_vectorBlock), intent(inout) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_scalarCubatureInfo), target :: rtempCubatureInfo
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo
    type(t_bmaVectorAssembly) :: rvectorAssembly
    integer :: ielementDistr,icubatureBlock,NEL
    integer, dimension(:), pointer :: p_IelementList
    integer(I32) :: ccubType,celement,ctrafoType
    type(t_spatialDiscretisation), pointer :: p_rdiscr
    type(t_scalarTrafoInfo), pointer :: p_rtrafoInfo
    type(t_scalarTrafoInfo), target :: rtempTrafoInfo

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => bma_perfconfig
    end if

    ! Fetch the first available spatial discretisation.
    ! It defines the way the elements are separated into blocks of elements.
    p_rdiscr => rvector%p_rblockDiscr%RspatialDiscr(1)

    ! If we do not have it, create a cubature info structure that
    ! defines how to do the assembly.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_createDefCubStructure(p_rdiscr,&
          rtempCubatureInfo,CUB_GEN_DEPR_BILFORM)
      p_rcubatureInfo => rtempCubatureInfo
    else
      p_rcubatureInfo => rcubatureInfo
    end if

    ! Transformation structure?
    if (present(rtrafoInfo)) then
      p_rtrafoInfo => rtrafoInfo
    else
      ! Set up a default transformation structure.
      ! Do we have a discretisation that allows us to set up the structure?
      if (associated(p_rdiscr)) then
        ! We found a discretisation. Set up the transformation
        call spdiscr_createDefTrafoStructure (p_rdiscr, rtempTrafoInfo)
      else
        ! Fallback procedure: No discretisation found.
        ! Create a standard transformation based on the cubature formula.
        call spdiscr_createDefTrafoByCubInfo (rcubatureInfo, rtempTrafoInfo)
      end if
      
      p_rtrafoInfo => rtempTrafoInfo
    end if

    ! Loop over the cubature blocks to discretise
    do icubatureBlock = 1,p_rcubatureInfo%ninfoBlockCount

      ! Get information about that block as well as an appropriate cubature formula
      call spdiscr_getStdDiscrInfo(icubatureBlock,p_rcubatureInfo,&
          p_rdiscr,ielementDistr,celement,ccubType,NEL,p_IelementList,&
          p_rtrafoInfo,ctrafoType=ctrafoType)

      ! Check if element distribution is empty
      if (NEL .le. 0 ) cycle

      ! Initialise a vector assembly structure for that element distribution
      call bma_initVecAssembly(rvectorAssembly,ccubType,ctrafoType,cflags,&
          rvector,ielementDistr,RmaxDerivativeTest,revalVectors,rperfconfig)

      ! Assemble the data for all elements in this element distribution
      call bma_assembleSubmeshVector (rvectorAssembly, p_IelementList(1:NEL),&
          fcalcLocalVectors, rcollection)

      ! Release the assembly structure.
      call bma_doneVecAssembly(rvectorAssembly)
    end do

    ! Release the assembly structure if necessary.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_releaseCubStructure(rtempCubatureInfo)
    end if

    ! Release the transformation structure if necessary.
    if (.not. present(rtrafoInfo)) then
      call spdiscr_releaseTrafoStructure (rtempTrafoInfo)
    end if

  end subroutine

  !****************************************************************************
  !****************************************************************************
  ! Assembly of integrals with block techniques
  !****************************************************************************
  !****************************************************************************

!<subroutine>

  subroutine bma_initIntAssemblyData(rassemblyData,revalVectors,ielementDistr)

!<description>
  ! Initialises the main parameters in the structure containing the
  ! assembly data (cubature weights, etc.).
  !
  ! The returned rassemblyData structure can be used as template
  ! for bma_createIntAssemblyData.
!</description>

!<input>
  ! OPTIONAL: A set of FEM functions to be evaluated during the vector
  ! assembly.
  type(t_fev2Vectors), intent(in), optional :: revalVectors

  ! OPTIONAL: Current element distribution.
  ! Must be specified if revalVectors is given.
  integer, intent(in), optional :: ielementDistr
!</input>

!<output>
  ! The integral assembly data structure to be initialised.
  type(t_bmaIntegralAssemblyData), intent(out) :: rassemblyData
!</output>

!</subroutine>

    ! Initialise the element set; can be simultaneously used for all
    ! FEM spaces, it is element independent.
    call elprep_init(rassemblyData%revalElementSet)

    ! Element sets are not yet initialised.  
    rassemblyData%ninitialisedElements = 0

    ! Prepare FEM data for the vectors to be evaluated.
    if (present(revalVectors)) then
      call fev2_prepareFemDataVecEval(revalVectors,rassemblyData%rfemDataBlocks,&
          ielementDistr)
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_doneIntAssemblyData(rassemblyData)

!<description>
  ! Releases allocated memory from the assembly data structure.
!</description>

!<inputoutput>
  ! Assembly data structure to be cleaned up.
  type(t_bmaIntegralAssemblyData), intent(inout) :: rassemblyData
!</inputoutput>

!</subroutine>

    ! Release the vector and FEM data array
    call fev2_cleanupFemData(rassemblyData%rfemDataBlocks)

    ! Element sets are not yet initialised.  
    rassemblyData%ninitialisedElements = 0

  end subroutine


  !****************************************************************************

!<subroutine>

  subroutine bma_createIntAssemblyData(rvectorAssembly,&
      rassemblyData,rassemblyDataTemplate,&
      revalVectors,revalVectorsTemplate)

!<description>
  ! Based on a template structure, this creates an assembly data structure
  ! for the vector assembly with all data arrays initialised and allocated.
!</description>

!<input>
  ! Template assembly data structure. rassemblyData is created based
  ! on the data in this structure.
  type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyDataTemplate

  ! An integral assembly structure.
  type(t_bmaIntegralAssembly), intent(in) :: rvectorAssembly

  ! Template vector data structure. rvectorData is created based
  ! on the data in this structure.
  type(t_fev2Vectors), intent(in) :: revalVectorsTemplate
!</input>

!<output>
  ! The vector assembly data structure to be initialised.
  type(t_bmaIntegralAssemblyData), intent(out) :: rassemblyData

  ! Vector data structure to be initialised according to revalVectorsTemplate
  type(t_fev2Vectors), intent(out) :: revalVectors
!</output>

!</subroutine>

    ! At first, just copy the content from the vector template structure
    rassemblyData = rassemblyDataTemplate

    ! Now initialise the dynamic content.
    !
    ! Copy the FEM data
    call fev2_copyFemData(rassemblyData%rfemDataBlocks,rassemblyDataTemplate%rfemDataBlocks)

    ! Initialise the FEM evaluation structures.
    call fev2_createFemData(rassemblyData%rfemDataBlocks,&
        rvectorAssembly%ncubp,rvectorAssembly%nelementsPerBlock)

    ! Initialise the vector evaluation structure.
    !
    ! Copy the content from the template.
    revalVectors = revalVectorsTemplate

    ! Re-create the vector array
    if (revalVectors%ncount .eq. 0) then
      nullify(revalVectors%p_RvectorData)
    else
      allocate(revalVectors%p_RvectorData(revalVectorsTemplate%ncount))
      revalVectors%p_RvectorData(1:revalVectorsTemplate%ncount) = &
          revalVectorsTemplate%p_RvectorData(1:revalVectorsTemplate%ncount)

      ! Prepare the evaluation of the FEM functions
      call fev2_initVectorEval(revalVectors,rassemblyData%rfemDataBlocks)
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_releaseIntAssemblyData(rassemblyData,revalVectors)

!<description>
  ! Releases memory of an assembly structure created by bma_createIntAssemblyData.
!</description>

!<inputoutput>
  ! The integral assembly data structure to be cleaned up
  type(t_bmaIntegralAssemblyData), intent(inout) :: rassemblyData

  ! Vector data structure to be cleaned up
  type(t_fev2Vectors), intent(inout) :: revalVectors
!</inputoutput>

!</subroutine>

    ! Release vector evaluzation data
    if (revalVectors%ncount .ne. 0) then
      call fev2_doneVectorEval(revalVectors)
      deallocate(revalVectors%p_RvectorData)
    end if

    ! Release FEM evaluation data
    call fev2_releaseFemData(rassemblyData%rfemDataBlocks)

    ! Release the element set and the cubature weights
    call elprep_releaseElementSet(rassemblyData%revalElementSet)

    if (associated(rassemblyData%p_DcubWeight)) then
      deallocate(rassemblyData%p_DcubWeight)
    end if

    ! Element sets are not initialised anymore
    rassemblyData%ninitialisedElements = 0

    ! Release memory
    call fev2_cleanupFemData (rassemblyData%rfemDataBlocks)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_initIntAssembly(rintegralAssembly,ccubType,ctrafoType,cflags,&
      rtriangulation,rboundary,revalVectors,ielementDistr,rperfconfig)

!<description>
  ! Initialise a vector assembly structure for assembling a linear form
  ! for a certain element distribution.
!</description>

!<input>
  ! Cubature formla to use
  integer(I32), intent(in) :: ccubType

  ! Transformation formla to use
  integer(I32), intent(in) :: ctrafoType

  ! Option field. Combination of BMA_CALC_xxxx flags.
  ! Use BMA_CALC_STANDARD for standard options.
  integer(I32), intent(in) :: cflags
  
  ! Underlying triangulation
  type(t_triangulation), target :: rtriangulation
  
  ! OPTIONAL: Underlying domain definition
  type(t_boundary), target, optional :: rboundary

  ! OPTIONAL: Set of vectors to be automatically evaluated
  type(t_fev2Vectors), intent(in), optional :: revalVectors

  ! OPTIONAL: Current element distribution.
  ! Must be specified if revalVectors is given.
  integer, intent(in), optional :: ielementDistr
  
  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<output>
  ! A vector assembly structure.
  type(t_bmaIntegralAssembly), intent(out) :: rintegralAssembly
!</output>

!</subroutine>

    integer :: i,ncount
    integer(I32) :: cshape
    type(t_fev2FemData), dimension(:), pointer :: p_RfemData

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Remember the performance configuration for later use
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => bma_perfconfig
    end if

    rintegralAssembly%p_rperfconfig => p_rperfconfig

    ! Get basic data
    if (present(rboundary)) then
      rintegralAssembly%p_rboundary => rboundary
    else
      nullify(rintegralAssembly%p_rboundary)
    end if
    rintegralAssembly%p_rtriangulation => rtriangulation

    ! Save the vector evaluation structure as template.
    ! Make a copy of the stucture, that is simpler here.
    if (present(revalVectors)) then
      rintegralAssembly%revalVectorsTemplate = revalVectors
    end if

    ! Initialise the template structure which is used during the
    ! actual vector assembly.
    call bma_initIntAssemblyData(rintegralAssembly%rassemblyDataTemplate,&
        revalVectors,ielementDistr)

    ! Get the cubature/transformation formula on/from the reference element
    rintegralAssembly%ccubType = ccubType
    rintegralAssembly%ncubp = cub_igetNumPts(ccubType)
    rintegralAssembly%ctrafoType = ctrafoType

    ! Allocate two arrays for the points and the weights
    allocate(rintegralAssembly%p_Domega(rintegralAssembly%ncubp))
    allocate(rintegralAssembly%p_DcubPtsRef(&
        trafo_igetReferenceDimension(rintegralAssembly%ctrafoType),&
        rintegralAssembly%ncubp))

    ! Get the cubature formula
    call cub_getCubature(ccubType,&
        rintegralAssembly%p_DcubPtsRef,rintegralAssembly%p_Domega)

    ! Number of simultaneously processed elements
    rintegralAssembly%nelementsPerBlock = p_rperfconfig%NELEMSIM

    ! Transformation type
    rintegralAssembly%ctrafoType = ctrafoType

    ! Get the evaluation tag by "OR"-ing all those from the FEM spaces
    ! and the option fields.
    rintegralAssembly%cevaluationTag = 0

    if (iand(cflags,BMA_CALC_REALCOORDS) .ne. 0) then
      ! Calculate real world coordinates of the cubature points
      rintegralAssembly%cevaluationTag = &
          ior(rintegralAssembly%cevaluationTag,EL_EVLTAG_REALPOINTS)
    end if

    if (rintegralAssembly%rassemblyDataTemplate%rfemDataBlocks%ncount .ne. 0) then

      ! Get a pointer to the FEM data
      p_RfemData => rintegralAssembly%rassemblyDataTemplate%rfemDataBlocks%p_RfemData
      ncount = rintegralAssembly%rassemblyDataTemplate%rfemDataBlocks%ncount

      ! Add the evaluation tags of all FEM spaces to a common evaluation tag.
      do i=1,ncount
        rintegralAssembly%cevaluationTag = ior(rintegralAssembly%cevaluationTag,&
            elem_getEvaluationTag(p_RfemData(i)%celement))
      end do

      ! Check that all elements have the same NVE. Otherwise,
      ! something is wrong.
      cshape = elem_igetShape(p_RfemData(1)%celement)

      do i=2,ncount
        if (cshape .ne. elem_igetShape(p_RfemData(i)%celement)) then
          call output_line ("Element spaces incompatible!",&
              OU_CLASS_ERROR,OU_MODE_STD,"bma_initVecAssembly")
          call sys_halt()
        end if
      end do
      
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_doneIntAssembly(rintegralAssembly)

!<description>
  ! Clean up an interal assembly structure.
!</description>

!<inputoutput>
  ! Vector assembly structure to clean up
  type(t_bmaIntegralAssembly), intent(inout) :: rintegralAssembly
!</inputoutput>

!</subroutine>

    ! Release all allocated memory.
    call bma_doneIntAssemblyData(rintegralAssembly%rassemblyDataTemplate)

    deallocate(rintegralAssembly%p_DcubPtsRef)
    deallocate(rintegralAssembly%p_Domega)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_prepareAssemblyData(rassemblyData,IelementList)

!<description>
  ! Auxiliary subroutine. 
  ! Initialises the data arrays.
!</description>

!<input>
  ! List of elements where to evaluate the FEM basis functions.
  integer, dimension(:), intent(in), target :: IelementList
!</input>

!<inputoutput>
  ! Assembly data structure.
  type(t_bmaIntegralAssemblyData), intent(inout), target :: rassemblyData
!</inputoutput>

!</subroutine>

    ! Remember the element list
    rassemblyData%p_IelementList => IelementList
    rassemblyData%nelements = size(IelementList)

    ! Calculate the DOF mapping for the FEM spaces.
    call fev2_calcDofMapping(rassemblyData%rfemDataBlocks,IelementList)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_evaluateFEMforInt(rassemblyData,rintegralAssembly)

!<description>
  ! Auxiliary subroutine. 
  ! Evaluates the FEM basis functions in the cubature points.
  ! Updates the element sets according to the coordinates of the
  ! cubature points, Jacobian determinants, etc.
!</description>

!<input>
  ! A vector assembly structure.
  type(t_bmaIntegralAssembly), intent(in), target :: rintegralAssembly
!</input>

!<inputoutput>
  ! Assembly data structure receiving the FEM data
  type(t_bmaIntegralAssemblyData), intent(inout) :: rassemblyData
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j
    integer(I32) :: cevaluationTag
    real(DP), dimension(:,:), pointer :: p_Ddetj,p_DcubWeight
    real(DP), dimension(:), pointer :: p_Domega

    if (rassemblyData%rfemDataBlocks%ncount .ne. 0) then
    
      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag.
      cevaluationTag = rintegralAssembly%cevaluationTag

      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      !
      ! If the cubature points are already initialised, do not do it again.
      ! We check this by taking a look to ninitialisedElements which
      ! gives the current maximum of initialised elements.
      if (rassemblyData%nelements .gt. &
          rassemblyData%ninitialisedElements) then

        ! (Re-)initialise!
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)

        ! Remember the new number of initialised elements
        rassemblyData%ninitialisedElements = rassemblyData%nelements

        ! (Re-)allocate the cubature weights      
        if (associated(rassemblyData%p_DcubWeight)) then
          deallocate(rassemblyData%p_DcubWeight)
        end if

        allocate(rassemblyData%p_DcubWeight(&
            rintegralAssembly%ncubp,rassemblyData%nelements))

      else
        ! No need.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (rassemblyData%revalElementSet,&
          cevaluationTag, rintegralAssembly%p_rtriangulation, &
          rassemblyData%p_IelementList, rintegralAssembly%ctrafoType, &
          rintegralAssembly%p_DcubPtsRef(:,1:rintegralAssembly%ncubp), &
          rperfconfig=rintegralAssembly%p_rperfconfig)

      ! Calculate the cubature weights. THese weights must be
      ! multiplied to the values of the trial/test-functions
      ! for cubature. They calculate from the actual cubature
      ! weights and the corresponding Jacobian determinant.

      p_Ddetj => rassemblyData%revalElementSet%p_Ddetj
      p_Domega => rintegralAssembly%p_Domega
      p_DcubWeight => rassemblyData%p_DcubWeight

      do j=1,rassemblyData%nelements
        do i=1,rintegralAssembly%ncubp

          ! Take the absolut value of the determinant of the mapping.
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that is normal!

          p_DcubWeight(i,j) = abs(p_Ddetj(i,j))*p_Domega(i)

        end do
      end do

      ! Now, calculate the FEM basis functions in the cubature points
      call fev2_evaluateFemData(rassemblyData%rfemDataBlocks,rassemblyData%revalElementSet)
      
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_assembleSubmeshIntegral(dintvalue,rintegralAssembly, IelementList,&
      fcalcLocalIntegral, rcollection)

!<description>
  ! Assembles the vector entries for a list of elements by integrating
  ! over the domain.
!</description>

!<input>

  ! List of elements where to assemble the bilinear form.
  integer, dimension(:), intent(in), target :: IelementList

  interface

    subroutine fcalcLocalIntegral (dintvalue,rassemblyData,rintegralAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

      use collection
      use blockmatassemblybase

      ! Calculates the value of the integral for a set of elements.
      
      ! Returns the value of the integral
      real(DP), intent(out) :: dintvalue

      ! Data necessary for the assembly. Contains determinants and
      ! cubature weights for the cubature,...
      type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

      ! Structure with all data about the assembly
      type(t_bmaIntegralAssembly), intent(in) :: rintegralAssembly

      ! Number of points per element
      integer, intent(in) :: npointsPerElement

      ! Number of elements
      integer, intent(in) :: nelements

      ! Values of FEM functions automatically evaluated in the
      ! cubature points.
      type(t_fev2Vectors), intent(in) :: revalVectors

      ! User defined collection structure
      type(t_collection), intent(inout), target, optional :: rcollection

    end subroutine

  end interface  

!</input>

!<inputoutput>

  ! A vector assembly structure prepared with bma_initMatAssembly.
  type(t_bmaIntegralAssembly), intent(inout), target :: rintegralAssembly

  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information.
  type(t_collection), intent(inout), target, optional :: rcollection

!</inputoutput>

!<output>
  ! Value of the integral
  real(DP), intent(out) :: dintvalue
!</output>

!</subroutine>

    ! local variables
    integer :: ielStart,ielMax
    real(DP) :: dinttemp
    type(t_bmaIntegralAssemblyData) :: rassemblyData
    type(t_fev2Vectors) :: revalVectors

    ! Use the template structure to create a vector data array.
    ! Allocate memory for the assembly.
    ! NOTE: This can be done in parallel, the information in rassemblyData
    ! is independent!
    call bma_createIntAssemblyData(rintegralAssembly,&
        rassemblyData,rintegralAssembly%rassemblyDataTemplate,&
        revalVectors,rintegralAssembly%revalVectorsTemplate)

    ! Loop blockwise through the element list
    dintvalue = 0.0_DP
    do ielStart = 1,size(IelementList),rintegralAssembly%nelementsPerBlock

      ! End of the current block
      ielMax = min(ielStart-1+rintegralAssembly%nelementsPerBlock,size(IelementList))

      ! Final preparation of the assembly data
      call bma_prepareAssemblyData(rassemblyData,IelementList(ielStart:ielMax))

      ! Calculate the FEM basis functions in the cubature points
      ! (for FEM functions being evaluated)
      call bma_evaluateFEMforInt(rassemblyData,rintegralAssembly)

      ! (Re-)allocate memory for the FEM evaluation if necessary.
      call fev2_prepareVectorEval(revalVectors,rassemblyData%revalElementSet)

      ! Evaluate the attached vectors in the cubature points.
      call fev2_evaluateVectors(revalVectors,rassemblyData%rfemDataBlocks)

      ! Use the callback routine to calculate the local vector entries.
      dinttemp = 0.0_DP
      call fcalcLocalIntegral(dinttemp,rassemblyData,rintegralAssembly,&
          rintegralAssembly%ncubp,ielMax-ielStart+1,revalVectors,rcollection)
          
      ! Sum up the integral
      dintvalue = dintvalue + dinttemp

    end do

    ! Release memory
    call bma_releaseIntAssemblyData(rassemblyData,revalVectors)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_buildIntegral (dintvalue,cflags,&
      fcalcLocalIntegral,rtriangulation,rboundary,rcollection, &
      revalVectors,rcubatureInfo,rtrafoInfo,rperfconfig)

!<description>
  ! This subroutine calculates the entries of a block vector and
  ! adds them to rvector. The callback function fcalcLocalVectors
  ! has to compute the local vector indices by cubature.
!</description>

!<input>

  ! Option field. Combination of BMA_CALC_xxxx flags.
  ! Use BMA_CALC_STANDARD for standard options.
  integer(I32), intent(in) :: cflags

  interface

    subroutine fcalcLocalIntegral (dintvalue,rassemblyData,rintegralAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

      use collection
      use blockmatassemblybase

      ! Calculates the value of the integral for a set of elements.
      
      ! Returns the value of the integral
      real(DP), intent(out) :: dintvalue

      ! Data necessary for the assembly. Contains determinants and
      ! cubature weights for the cubature,...
      type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

      ! Structure with all data about the assembly
      type(t_bmaIntegralAssembly), intent(in) :: rintegralAssembly

      ! Number of points per element
      integer, intent(in) :: npointsPerElement

      ! Number of elements
      integer, intent(in) :: nelements

      ! Values of FEM functions automatically evaluated in the
      ! cubature points.
      type(t_fev2Vectors), intent(in) :: revalVectors

      ! User defined collection structure
      type(t_collection), intent(inout), target, optional :: rcollection

    end subroutine

  end interface  

  ! OPTIONAL: Underlying triangulation.
  ! Cal be omitted if a vector is passed via revalVectors.
  type(t_triangulation), target, optional :: rtriangulation
  
  ! OPTIONAL: Underlying domain definition
  type(t_boundary), target, optional :: rboundary

  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information.
  type(t_collection), intent(inout), target, optional :: rcollection

  ! OPTIONAL: Set of vectors to be automatically evaluated
  type(t_fev2Vectors), intent(in), optional :: revalVectors

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), target, optional :: rcubatureInfo

  ! OPTIONAL: A transformation structure that specifies the transformation
  ! from the reference to the real element(s).
  ! If not specified, default settings are used.
  type(t_scalarTrafoInfo), intent(in), target, optional :: rtrafoInfo

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<output>
  ! Calculated integral
  real(DP), intent(out) :: dintvalue
!</output>

!</subroutine>

    ! local variables
    type(t_scalarCubatureInfo), target :: rtempCubatureInfo
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo
    type(t_bmaIntegralAssembly) :: rintegralAssembly
    integer :: ielementDistr,icubatureBlock,NEL
    integer, dimension(:), pointer :: p_IelementList
    integer(I32) :: ccubType,celement,ctrafoType
    type(t_spatialDiscretisation), pointer :: p_rdiscr
    type(t_scalarTrafoInfo), pointer :: p_rtrafoInfo
    type(t_scalarTrafoInfo), target :: rtempTrafoInfo
    integer :: i
    real(DP) :: dinttemp
    type(t_triangulation), pointer :: p_rtriangulation

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => bma_perfconfig
    end if
    
    ! Try to find a discretisation structure.
    nullify(p_rdiscr)
    
    if (present(revalVectors)) then
      do i=1,revalVectors%ncount
        if (associated(revalVectors%p_RvectorData(1)%p_rvector)) then
          p_rdiscr => revalVectors%p_RvectorData(1)%p_rvector%p_rspatialDiscr
          exit
        end if
      end do
    end if
    
    ! Mesh?
    if (present(rtriangulation)) then
      ! Given by the parameter
      p_rtriangulation => rtriangulation
    else if (associated(p_rdiscr)) then
      ! Taken from the discretisation
      p_rtriangulation => p_rdiscr%p_rtriangulation
    end if
    
    if (.not. associated(p_rtriangulation)) then
      call output_line ("Cannot determine underlying mesh.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_buildIntegral")
      call sys_halt()
    end if
    
    ! If we do not have it, create a cubature info structure that
    ! defines how to do the assembly.
    if (.not. present(rcubatureInfo)) then
      if (associated(p_rdiscr)) then
        call spdiscr_createDefCubStructure(p_rdiscr,&
            rtempCubatureInfo,CUB_GEN_DEPR_BILFORM)
      else
        ! If there is no cubature and no discretisation given,
        ! there is no chance to determine
        call output_line ("No cubature rule and no discretisation given.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_buildIntegral")
        call sys_halt()
      end if
      p_rcubatureInfo => rtempCubatureInfo
    else
      p_rcubatureInfo => rcubatureInfo
    end if

    ! Transformation structure?
    if (present(rtrafoInfo)) then
      p_rtrafoInfo => rtrafoInfo
    else
      ! Set up a default transformation structure.
      ! Do we have a discretisation that allows us to set up the structure?
      if (associated(p_rdiscr)) then
        ! We found a discretisation. Set up the transformation
        call spdiscr_createDefTrafoStructure (p_rdiscr, rtempTrafoInfo)
      else
        ! Fallback procedure: No discretisation found.
        ! Create a standard transformation based on the cubature formula.
        call spdiscr_createDefTrafoByCubInfo (rcubatureInfo, rtempTrafoInfo)
      end if
      
      p_rtrafoInfo => rtempTrafoInfo
    end if

    ! Loop over the cubature blocks to calculate
    dintvalue = 0.0_DP
    do icubatureBlock = 1,p_rcubatureInfo%ninfoBlockCount

      ! Get information about that block as well as an appropriate cubature formula
      call spdiscr_getStdDiscrInfo(icubatureBlock,p_rcubatureInfo,&
          p_rdiscr,ielementDistr,celement,ccubType,NEL,p_IelementList,&
          p_rtrafoInfo,ctrafoType=ctrafoType)

      ! Check if element distribution is empty
      if (NEL .le. 0 ) cycle

      ! Initialise a vector assembly structure for that element distribution
      call bma_initIntAssembly(rintegralAssembly,ccubType,ctrafoType,cflags,&
          p_rtriangulation,rboundary,revalVectors,ielementDistr,rperfconfig)

      ! Assemble the data for all elements in this element distribution
      dinttemp = 0.0_DP
      call bma_assembleSubmeshIntegral (dinttemp,rintegralAssembly, p_IelementList(1:NEL),&
          fcalcLocalIntegral, rcollection)
      
      ! Sum up
      dintvalue = dintvalue + dinttemp

      ! Release the assembly structure.
      call bma_doneIntAssembly(rintegralAssembly)
    end do

    ! Release the assembly structure if necessary.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_releaseCubStructure(rtempCubatureInfo)
    end if
    
    ! Release the transformation structure if necessary.
    if (.not. present(rtrafoInfo)) then
      call spdiscr_releaseTrafoStructure (rtempTrafoInfo)
    end if

  end subroutine

end module
