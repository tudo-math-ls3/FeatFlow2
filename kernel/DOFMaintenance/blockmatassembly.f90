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
!# Finally the following set of callback routines realise standard operators.
!# They can directly be used with bma_buildMatrix:
!#
!# 1.) bma_fcalc_Mass
!#     -> Used with bma_buildMatrix, this calculates mass matrices in all
!#        diagonal blocks of a block matrix.
!#
!# 2.) bma_fcalc_Laplace
!#     -> Used with bma_buildMatrix, this calculates Laplace matrices matrices
!#        in all diagonal blocks of a block matrix.
!# </purpose>
!##############################################################################

module blockmatassembly

  use fsystem
  use storage
  use collection, only: t_collection
  use basicgeometry
  use boundaryaux
  use cubature
  use dofmapping
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

  public :: bma_fcalc_Mass
  public :: bma_fcalc_Laplace
  public :: bma_fcalc_rhsOne

  public :: bma_initMatAssembly
  public :: bma_doneMatAssembly
  public :: bma_assembleSubmeshMatrix
  public :: bma_buildMatrix

  public :: bma_initVecAssembly
  public :: bma_doneVecAssembly
  public :: bma_assembleSubmeshVector
  public :: bma_buildVector

contains

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_Mass(RmatrixData,rassemblyData,rmatrixAssembly,&
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
    type(t_collection), intent(inout), optional :: rcollection
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

  subroutine bma_fcalc_Laplace(RmatrixData,rassemblyData,rmatrixAssembly,&
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
    type(t_collection), intent(inout), optional :: rcollection
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

  subroutine bma_fcalc_rhsOne(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f=1 in all components.
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: rvectorData
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
    type(t_collection), intent(inout), optional :: rcollection
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

          ! Outer loop over the DOF's i=1..ndof on our current element,
          ! which corresponds to the (test) basis functions Phi_i:
          do idofe=1,p_rvectorData%ndofTest
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! into dbasI
            dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
            
            ! The value of the RHS is just the 1.0-function here
            dval = 1.0_DP
            
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

  integer function containsDiscr (Rlist,nentries,rentry)
  ! returns <> 0 (the index) if the list Rlist contains rentry
  type(t_fev2FemData), dimension(:), intent(in) :: Rlist
  type(t_spatialDiscretisation), intent(in), target :: rentry
  integer, intent(in) :: nentries

    integer :: i

    do i=1,nentries
      if (associated(Rlist(i)%p_rdiscr,rentry)) then
        containsDiscr = i
        return
      end if
    end do

    containsDiscr = 0

  end function

  !****************************************************************************

!<subroutine>

  subroutine bma_prepareMatrixData(&
      rmatrix,RmatrixData,RfemData,ielementDistr)

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

  ! List of all FEM spaces used during the assembly
  type(t_fev2FemData), dimension(:), intent(in) :: RfemData
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
          RmatrixData(i,j)%iidxFemDataTrial = containsDiscr(RfemData,size(RfemData),p_rdiscrTrial)
          RmatrixData(i,j)%iidxFemDataTest = containsDiscr (RfemData,size(RfemData),p_rdiscrTest)
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
      rmatrix,ielementDistr,RmaxDerivativeTest,RmaxDerivativeTrial)

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
    call fev2_prepareFemDataBMat(rmatrix,rassemblyData%p_RfemData,&
        ielementDistr,RmaxDerivativeTest,RmaxDerivativeTrial)

    ! Initialise the blocks according to the structure of the matrix
    ! and the FEM spaces used.
    allocate(rassemblyData%p_RmatrixData( &
        ubound(rmatrix%RmatrixBlock,1),ubound(rmatrix%RmatrixBlock,2)))

    call bma_prepareMatrixData(&
        rmatrix,rassemblyData%p_RmatrixData,rassemblyData%p_RfemData,ielementDistr)

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
    deallocate(rassemblyData%p_RfemData)
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

    ! Not initialise the dynamic content.
    !
    ! Copy the FEM data and the matrix data
    allocate (rassemblyData%p_RfemData(size(rassemblyDataTemplate%p_RfemData)))
    allocate (rassemblyData%p_RmatrixData(&
        ubound(rassemblyDataTemplate%p_RmatrixData,1),&
        ubound(rassemblyDataTemplate%p_RmatrixData,2)))

    rassemblyData%p_RfemData(:) = rassemblyDataTemplate%p_RfemData(:)
    rassemblyData%p_RmatrixData(:,:) = rassemblyDataTemplate%p_RmatrixData(:,:)

    ! Initialise the FEM evaluation structures.
    call fev2_createFemData(rassemblyData%p_RfemData,&
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
              rassemblyData%p_RfemData(p_rmatrixData%iidxFemDataTrial)%p_Dbas
          p_rmatrixData%p_DbasTest => &
              rassemblyData%p_RfemData(p_rmatrixData%iidxFemDataTest)%p_Dbas

          p_rmatrixData%p_IdofsTrial => &
              rassemblyData%p_RfemData(p_rmatrixData%iidxFemDataTrial)%p_Idofs
          p_rmatrixData%p_IdofsTest => &
              rassemblyData%p_RfemData(p_rmatrixData%iidxFemDataTest)%p_Idofs

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
      revalVectors%p_RvectorData(:) = revalVectorsTemplate%p_RvectorData(:)

      ! Allocate memory for the values in the cubature points
      call fev2_initVectorEval(revalVectors,rassemblyData%p_RfemData,&
          rassemblyData%revalElementSet)
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
    call fev2_releaseFemData(rassemblyData%p_RfemData)

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
    deallocate (rassemblyData%p_RfemData)
    deallocate (rassemblyData%p_RmatrixData)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_initMatAssembly(rmatrixAssembly,ccubType,cflags,&
      rmatrix,ielementDistr,RmaxDerivativeTest,RmaxDerivativeTrial,&
      revalVectors,rperfconfig)

!<description>
  ! Initialise a matrix assembly structure for assembling a bilinear form
  ! for a certain element distribution.
!</description>

!<input>
  ! Cubature formla to use
  integer(I32), intent(in) :: ccubType

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

    integer :: i
    integer(I32) :: cshape

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
        rmatrix,ielementDistr,RmaxDerivativeTest,RmaxDerivativeTrial)

    ! Get the transformation from the first FEM data structure.
    rmatrixAssembly%ctrafoType = &
        rmatrixAssembly%rassemblyDataTemplate%p_RfemData(1)%ctrafoType

    ! Get the cubature formula on the reference element
    rmatrixAssembly%ccubType = ccubType
    rmatrixAssembly%ncubp = cub_igetNumPts(ccubType)

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

    if (.not. associated(rmatrixAssembly%rassemblyDataTemplate%p_RfemData)) then
      call output_line ("Matrix is empty, nothing to assemble.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_initMatAssembly")
      call sys_halt()
    end if

    ! Add the evaluation tags of all FEM spaces to a common evaluation tag.
    do i=1,size(rmatrixAssembly%rassemblyDataTemplate%p_RfemData)
      rmatrixAssembly%cevaluationTag = ior(rmatrixAssembly%cevaluationTag,&
          elem_getEvaluationTag(rmatrixAssembly%rassemblyDataTemplate%p_RfemData(i)%celement))
    end do

    ! Check that all elements have the same NVE. Otherwise,
    ! something is wrong.
    cshape = &
        elem_igetShape(rmatrixAssembly%rassemblyDataTemplate%p_RfemData(1)%celement)

    do i=2,size(rmatrixAssembly%rassemblyDataTemplate%p_RfemData)
      if (cshape .ne. &
          elem_igetShape(rmatrixAssembly%rassemblyDataTemplate%p_RfemData(i)%celement)) then
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
    type(t_fev2FemData), pointer :: p_rfemData
    type(t_bmaMatrixData), pointer :: p_rmatrixData

    ! Remember the element list
    rassemblyData%p_IelementList => IelementList
    rassemblyData%nelements = size(IelementList)

    ! Now, loop through all FEM spaces to calculate the FEM basis functions
    do i=1,size(rassemblyData%p_RfemData)

      p_rfemData => rassemblyData%p_RfemData(i)

      ! Calculates the degrees of freedom on the elements
      call dof_locGlobMapping_mult(p_rfemData%p_rdiscr, &
          IelementList, p_rfemData%p_Idofs)

    end do

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
    call fev2_evaluateFemData(rassemblyData%p_RfemData,rassemblyData%revalElementSet)

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

            do m=1,ubound(p_Dentry,3)
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

            do m=1,ubound(p_DentryIntl,4)
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
      type(t_collection), intent(inout), optional :: rcollection

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

      ! Evaluate the attached vectors in the cubature points.
      call fev2_evaluateVectors(revalVectors,rassemblyData%p_RfemData)

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
      revalVectors,rcubatureInfo,rperfconfig)

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
      type(t_collection), intent(inout), optional :: rcollection

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
    integer(I32) :: ccubType,celement
    type(t_spatialDiscretisation), pointer :: p_rdiscr

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

    ! Loop over the cubature blocks to discretise
    do icubatureBlock = 1,p_rcubatureInfo%ninfoBlockCount

      ! Get information about that block as well as an appropriate cubature formula
      call spdiscr_getStdDiscrInfo(icubatureBlock,p_rcubatureInfo,&
          p_rdiscr,ielementDistr,celement,ccubType,NEL,p_IelementList)

      ! Check if element distribution is empty
      if (NEL .le. 0 ) cycle

      ! Initialise a matrix assembly structure for that element distribution
      call bma_initMatAssembly(rmatrixAssembly,ccubType,cflags,&
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

  end subroutine

  ! #################################################################

!****************************************************************************

!<subroutine>

  subroutine bma_prepareVectorData(&
      rvector,RvectorData,RfemData,ielementDistr)

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

  ! List of all FEM spaces used during the assembly
  type(t_fev2FemData), dimension(:), intent(in) :: RfemData
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
      RvectorData(i)%iidxFemDataTest = containsDiscr (RfemData,size(RfemData),p_rdiscrTest)

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
      rvector,ielementDistr,RmaxDerivativeTest)

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
    call fev2_prepareFemDataBVec(rvector,rassemblyData%p_RfemData,&
        ielementDistr,RmaxDerivativeTest)

    ! Initialise the blocks according to the structure of the vector
    ! and the FEM spaces used.
    allocate(rassemblyData%p_RvectorData(ubound(rvector%RvectorBlock,1)))

    call bma_prepareVectorData(&
        rvector,rassemblyData%p_RvectorData,rassemblyData%p_RfemData,ielementDistr)

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
    deallocate(rassemblyData%p_RfemData)
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

    ! Not initialise the dynamic content.
    !
    ! Copy the FEM data and the vector data
    allocate (rassemblyData%p_RfemData(size(rassemblyDataTemplate%p_RfemData)))
    allocate (rassemblyData%p_RvectorData(&
        ubound(rassemblyDataTemplate%p_RvectorData,1)))

    rassemblyData%p_RfemData(:) = rassemblyDataTemplate%p_RfemData(:)
    rassemblyData%p_RvectorData(:) = rassemblyDataTemplate%p_RvectorData(:)

    ! Initialise the FEM evaluation structures.
    call fev2_createFemData(rassemblyData%p_RfemData,&
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
          rassemblyData%p_RfemData(p_rvectorData%iidxFemDataTest)%p_Dbas

      p_rvectorData%p_IdofsTest => &
          rassemblyData%p_RfemData(p_rvectorData%iidxFemDataTest)%p_Idofs

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
      revalVectors%p_RvectorData(:) = revalVectorsTemplate%p_RvectorData(:)

      ! Allocate memory for the values in the cubature points
      call fev2_initVectorEval(revalVectors,rassemblyData%p_RfemData,&
          rassemblyData%revalElementSet)
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
    call fev2_releaseFemData(rassemblyData%p_RfemData)

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
    deallocate (rassemblyData%p_RfemData)
    deallocate (rassemblyData%p_rvectorData)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_initVecAssembly(rvectorAssembly,ccubType,cflags,&
      rvector,ielementDistr,RmaxDerivativeTest,&
      revalVectors,rperfconfig)

!<description>
  ! Initialise a vector assembly structure for assembling a linear form
  ! for a certain element distribution.
!</description>

!<input>
  ! Cubature formla to use
  integer(I32), intent(in) :: ccubType

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

    integer :: i
    integer(I32) :: cshape

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
        rvector,ielementDistr,RmaxDerivativeTest)

    ! Get the transformation from the first FEM data structure.
    rvectorAssembly%ctrafoType = &
        rvectorAssembly%rassemblyDataTemplate%p_RfemData(1)%ctrafoType

    ! Get the cubature formula on the reference element
    rvectorAssembly%ccubType = ccubType
    rvectorAssembly%ncubp = cub_igetNumPts(ccubType)

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

    if (.not. associated(rvectorAssembly%rassemblyDataTemplate%p_RfemData)) then
      call output_line ("Vector is empty, nothing to assemble.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_initVecAssembly")
      call sys_halt()
    end if

    ! Add the evaluation tags of all FEM spaces to a common evaluation tag.
    do i=1,size(rvectorAssembly%rassemblyDataTemplate%p_RfemData)
      rvectorAssembly%cevaluationTag = ior(rvectorAssembly%cevaluationTag,&
          elem_getEvaluationTag(rvectorAssembly%rassemblyDataTemplate%p_RfemData(i)%celement))
    end do

    ! Check that all elements have the same NVE. Otherwise,
    ! something is wrong.
    cshape = &
        elem_igetShape(rvectorAssembly%rassemblyDataTemplate%p_RfemData(1)%celement)

    do i=2,size(rvectorAssembly%rassemblyDataTemplate%p_RfemData)
      if (cshape .ne. &
          elem_igetShape(rvectorAssembly%rassemblyDataTemplate%p_RfemData(i)%celement)) then
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
    type(t_fev2FemData), pointer :: p_rfemData
    type(t_bmaVectorData), pointer :: p_rvectorData

    ! Remember the element list
    rassemblyData%p_IelementList => IelementList
    rassemblyData%nelements = size(IelementList)

    ! Now, loop through all FEM spaces to calculate the FEM basis functions
    do i=1,size(rassemblyData%p_RfemData)

      p_rfemData => rassemblyData%p_RfemData(i)

      ! Calculates the degrees of freedom on the elements
      call dof_locGlobMapping_mult(p_rfemData%p_rdiscr, &
          IelementList, p_rfemData%p_Idofs)

    end do

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
    call fev2_evaluateFemData(rassemblyData%p_RfemData,rassemblyData%revalElementSet)

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

        do l=1,ubound(p_Dentry,2)
          do k=1,ubound(p_Dentry,1)
            p_Ddata(p_IdofsTest(k,l)) = p_Ddata(p_IdofsTest(k,l)) + p_Dentry(k,l)
          end do
        end do

      else

        ! Interleaved
        p_DentryIntl => p_rvectorData%p_DentryIntl
        nvar = p_rvectorData%nvar

        do l=1,ubound(p_DentryIntl,3)
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
      type(t_collection), intent(inout), optional :: rcollection

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

      ! Evaluate the attached vectors in the cubature points.
      call fev2_evaluateVectors(revalVectors,rassemblyData%p_RfemData)

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
      RmaxDerivativeTest,revalVectors,rcubatureInfo,rperfconfig)

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
      type(t_collection), intent(inout), optional :: rcollection

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
    integer(I32) :: ccubType,celement
    type(t_spatialDiscretisation), pointer :: p_rdiscr

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

    ! Loop over the cubature blocks to discretise
    do icubatureBlock = 1,p_rcubatureInfo%ninfoBlockCount

      ! Get information about that block as well as an appropriate cubature formula
      call spdiscr_getStdDiscrInfo(icubatureBlock,p_rcubatureInfo,&
          p_rdiscr,ielementDistr,celement,ccubType,NEL,p_IelementList)

      ! Check if element distribution is empty
      if (NEL .le. 0 ) cycle

      ! Initialise a vector assembly structure for that element distribution
      call bma_initVecAssembly(rvectorAssembly,ccubType,cflags,&
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

  end subroutine

end module
