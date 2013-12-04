! This module contains the core of the compare-application
! It starts at the point where you have 2 FEM Solutions,
! created the according discretisations and parametrisations
! and now want to calculate the differences.

module ExtFEcomparer_core

  use fsystem
  use storage
  use genoutput

  use element
  use cubature
  use spatialdiscretisation
  use linearsystemblock
  use linearsystemscalar

  use bilinearformevaluation
  use derivatives
  use collection
  use feevaluation2
  use feevaluation

  use blockmatassembly
  use blockmatassemblybase
  use blockmatassemblystdop

  use ExtFEcomparer_typedefs

  use fparser
  use paramlist

  use vectorio

implicit none

contains

  ! This routine shall prepare everything for the calculation
  ! It is the one that is called from outside. The other routines
  ! (that do the actual calculation) are called by this routine.

  subroutine calculate_difference_2D(rproblem_1,rproblem_2)

  type(t_problem), intent(inout), target :: rproblem_1, rproblem_2

  ! Level where we want to evaluate
  integer :: lvl1,lvl2

  ! We need 2 Block discretisations
  !type(t_blockDiscretisation), pointer :: rblockDiscr1, rblockDiscr2

  ! and the spatial discretisation on the lvls
  !type(t_spatialDiscretisation), pointer :: rspatialDiscr1, rspatialDiscr2

  ! We also need to create 2 block vectors so we can attatch
  ! the discretisation to it
  type(t_vectorBlock), pointer :: rx1, rx2

  ! We want to calculate integrals, so we need the cubature info
  type(t_scalarCubatureInfo), target:: rcubatureInfoFirst, rcubatureInfoSecond

  ! We need to evaluate our FEM-function in the cubature points
  ! for that purpose we need a vector
  type(t_fev2Vectors) :: rcoeffVectors

 ! We need the collection to pass the second vector to the
 ! actual computation
 type(t_collection) :: rcollection

  ! Local variables
  real(DP) :: dintvalue

  ! For the output
  real(DP), dimension(3) :: results
  integer :: writeOut
  character(LEN=SYS_STRLEN):: outFileName

  !======================================!
  ! We might have everything, so we start!
  !======================================!

  lvl1 = rproblem_1%NLMAX
  lvl2 = rproblem_2%NLMAX


  ! This part is actually neccesary because we want to have the
  ! situation that each block vector has its blockdiscretisation
  ! attatched to it.

  rx1 => rproblem_1%coeffVector
  rx2 => rproblem_2%coeffVector

  ! Attatch the Block-Discritisation to rx1 and rx2
  rx1%p_rblockDiscr => rproblem_1%RlevelInfo(lvl1)%rdiscretisation
  rx2%p_rblockDiscr => rproblem_2%RlevelInfo(lvl2)%rdiscretisation

  !Now we attatch the spatial discretisations to rx1 and rx2
  rx1%RvectorBlock(1)%p_rspatialDiscr => rproblem_1%RlevelInfo(lvl1)%rdiscretisation%RspatialDiscr(1)
  rx1%RvectorBlock(2)%p_rspatialDiscr => rproblem_1%RlevelInfo(lvl1)%rdiscretisation%RspatialDiscr(1)

  rx2%RvectorBlock(1)%p_rspatialDiscr => rproblem_2%RlevelInfo(lvl2)%rdiscretisation%RspatialDiscr(1)
  rx2%RvectorBlock(2)%p_rspatialDiscr => rproblem_2%RlevelInfo(lvl2)%rdiscretisation%RspatialDiscr(1)

  ! In the subroutines we need the parameter list
  ! The lists in both problem structures are the same, so
  ! we can pick any
  rcollection%p_rparlistQuickAccess1 => rproblem_1%rparamlist


  ! ###################################
  ! # FIRST Integral: The common part #
  ! ###################################

  ! We need the cubature informations
  ! First step: only for Integral ||u1-u2|| on common part, we can
  ! actually use the cubature info from any discretisation
  ! So we use the discrectisation from vector 1

  call spdiscr_createDefCubStructure(rx1%RvectorBlock(1)%p_rspatialDiscr,rcubatureInfoFirst,CUB_GEN_AUTO_G3)

  ! first the u component
  call fev2_addVectorToEvalList(rcoeffVectors,rx1%RvectorBlock(1),0)
  ! then the v component
  call fev2_addVectorToEvalList(rcoeffVectors,rx1%RvectorBlock(2),0)

  ! we cannot pass the second vector via this way because it lives
  ! on a different mesh. Therefore we pass it via the collection
  ! The evaluation of this one will be expensive
  rcollection%p_rvectorQuickAccess1=> rx2


  ! We need space in the vector for the second function
  call fev2_addDummyVectorToEvalList(rcoeffVectors)
  call fev2_addDummyVectorToEvalList(rcoeffVectors)

  ! Calculate the L2 error
  call bma_buildIntegral(dintvalue,BMA_CALC_STANDARD,calc_L2error_velocity_2D_common, &
        rcollection=rcollection,revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfoFirst)

  dintvalue = sqrt(dintvalue)
  call output_line ("L2-error on the common part= "//trim(sys_sdEL(dintvalue,10)))

  ! Store it in the result array
  results(1) = dintvalue

  ! Release everything
  call fev2_releaseVectorList(rcoeffVectors)
  call spdiscr_releaseCubStructure(rcubatureInfoFirst)
  rcollection%p_rvectorQuickAccess1 => NULL()
  dintvalue = 0.0_DP




    !######################################
    ! Second integral: Only first function!
    !#####################################!

    ! Pass function 1 to bma_integral
    call fev2_addVectorToEvalList(rcoeffVectors,rx1%RvectorBlock(1),0)
    call fev2_addVectorToEvalList(rcoeffVectors,rx1%RvectorBlock(2),0)

    ! Create the according cubature structure
    call spdiscr_createDefCubStructure(rx1%RvectorBlock(1)%p_rspatialDiscr,rcubatureInfoFirst,CUB_GEN_AUTO_G3)

    ! Calculate the L2 error
    call bma_buildIntegral(dintvalue,BMA_CALC_STANDARD,calc_L2norm_velocity_2D_first, &
            rcollection=rcollection,revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfoFirst)

    dintvalue = sqrt(dintvalue)
    call output_line ("L2-norm of the first function= "//trim(sys_sdEL(dintvalue,10)))

    ! Store it in the result array
    results(2) = dintvalue

    ! Release everything
    call fev2_releaseVectorList(rcoeffVectors)
    call spdiscr_releaseCubStructure(rcubatureInfoFirst)
    dintvalue = 0.0_DP

    !###########################################!
    !# Third integral: Only the second function#!
    !###########################################!

    ! Pass function 2 to bma_build_integral
    call fev2_addVectorToEvalList(rcoeffVectors,rx2%RvectorBlock(1),0)
    call fev2_addVectorToEvalList(rcoeffVectors,rx2%RvectorBlock(2),0)

    ! Create the according cubature structure
    call spdiscr_createDefCubStructure(rx2%RvectorBlock(1)%p_rspatialDiscr,rcubatureInfoSecond,CUB_GEN_AUTO_G3)

   ! Calculate the L2 error
    call bma_buildIntegral(dintvalue,BMA_CALC_STANDARD,calc_L2norm_velocity_2D_second, &
            rcollection=rcollection,revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfoSecond)

    dintvalue = sqrt(dintvalue)
    call output_line ("L2-norm of the second function= "//trim(sys_sdEL(dintvalue,10)))

    ! Write it out in the results array
    results(3) = dintvalue

    ! Release everything
    call fev2_releaseVectorList(rcoeffVectors)
    call spdiscr_releaseCubStructure(rcubatureInfoSecond)
    dintvalue = 0.0_DP



    !##################################################
    !# Postprocessing: write the results out in a file#
    !##################################################

    ! Check what the user set
    call parlst_getvalue_int(rproblem_1%rparamlist,"POSTPROCESSING","writeOut",writeOut,1)
    call parlst_getvalue_string (rproblem_1%rparamList,"POSTPROCESSING",&
                                 "filename",outFileName,bdequote=.true.)

    if(writeOut > 0.5_DP) then
        call vecio_writeArray(results,0,outFileName,'(E20.10)')
    end if

  end subroutine



subroutine calc_L2error_velocity_2D_common(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

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
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dy,dval1,dval2, work
    real(DP), dimension(:), pointer :: dval_u21, dval_u22

    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc11, p_Dfunc12 ,p_Dfunc21, p_Dfunc22
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
    type(t_vectorBlock), pointer :: p_rx
    integer, dimension(:), pointer :: p_Ielements

    ! We need a parser to read in the values
    type(t_fparser) :: rparser
    ! storage for the parser
    character(256) :: domainDescription
    ! And we need the parameter list
    type(t_parlist) :: rparamlist


    ! Read out the area description in the string
    call parlst_getvalue_string(rcollection%p_rparlistQuickAccess1, 'DOMAININFO','common',domainDescription)

    call fparser_create(rparser,1)
    call fparser_parseFunction(rparser,1,domainDescription,(/'x','y'/))

    ! Now we can evaluate the area description with the parser-tools

    ! Get cubature weigths
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values we already have
    p_Dfunc11 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)
    p_Dfunc12 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_FUNC)

    ! and a pointer to the second function
    p_rx => rcollection%p_rvectorQuickAccess1

    ! We need to evaluate the second function in the cubature points
    ! IF we want to evaluate the function in a point where the function
    ! is not defined, we return a 0 as function value. Else we will
    ! get many errors and we do not want that!

    ! first: set a pointer to the vector where we want to store
    ! the values. We want it to look as if we added the vector
    ! to the evaluation list and passed it via bma_build_integral
    p_Dfunc21 => revalVectors%p_RvectorData(3)%p_Ddata(:,:,DER_FUNC)
    p_Dfunc22 => revalVectors%p_RvectorData(4)%p_Ddata(:,:,DER_FUNC)

    ! Now do the evaluation
    ! Because we know we evaluate the function at points where it does
    ! not have function values and we don't want to get error-messages,
    ! we assume zeros there.
    do iel=1,nelements
    call fevl_evaluate(DER_FUNC,p_Dfunc21(:,iel), &
                p_rx%RvectorBlock(1),p_Dpoints(:,:,iel), cnonmeshPoints=FEVL_NONMESHPTS_ZERO)
        call fevl_evaluate(DER_FUNC,p_Dfunc22(:,iel), &
                p_rx%RvectorBlock(2),p_Dpoints(:,:,iel), cnonmeshPoints=FEVL_NONMESHPTS_ZERO)
    end do


    ! We should have everything prepeared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coodinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0

           call fparser_evalFunction(rparser,1,(/dx,dy/),work)

           ! work can only be 0 or 1, so if it is < 0.5 it is 0
           if (abs(work) .le. 0.5_DP) then
                ! We dont need to work
                Dintvalue(1) = Dintvalue(1)
           else
                ! We need to work

                ! We did all evaluation earlier
                dval1 = (p_Dfunc11(icubp,iel) - p_Dfunc21(icubp,iel))**2
                dval2 = (p_Dfunc12(icubp,iel) - p_Dfunc22(icubp,iel))**2
                Dintvalue(1) = Dintvalue(1) + (dval1 + dval2)*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel


end subroutine

subroutine calc_L2norm_velocity_2D_first(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

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
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dy,dval1,dval2, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc1, p_Dfunc2
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
    integer, dimension(:), pointer :: p_Ielements

    ! We need a parser to read in the values
    type(t_fparser) :: rparser
    ! storage for the parser
    character(256) :: domainDescription
    ! And we need the parameter list
    type(t_parlist) :: rparamlist


    ! Read out the area description in the string
    call parlst_getvalue_string(rcollection%p_rparlistQuickAccess1, 'DOMAININFO','onlyFirst',domainDescription)

    call fparser_create(rparser,1)
    call fparser_parseFunction(rparser,1,domainDescription,(/'x','y'/))

    ! Now we can evaluate the area description with the parser-tools

    ! Get cubature weigths
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)
    p_Dfunc2 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_FUNC)

    ! We should have everything prepeared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coodinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0

           call fparser_evalFunction(rparser,1,(/dx,dy/),work)

           ! work can only be 0 or 1, so if it is < 0.5 it is 0
           if (abs(work) .le. 0.5_DP) then
                ! We dont need to work
                Dintvalue(1) = Dintvalue(1)
           else
                ! We need to work
                dval1 = p_Dfunc1(icubp,iel)**2
                dval2 = p_Dfunc2(icubp,iel)**2
                Dintvalue(1) = Dintvalue(1) + (dval1 + dval2)*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel


end subroutine


subroutine calc_L2norm_velocity_2D_second(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

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
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dy,dval1,dval2, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc1, p_Dfunc2
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
    integer, dimension(:), pointer :: p_Ielements

    ! We need a parser to read in the values
    type(t_fparser) :: rparser
    ! storage for the parser
    character(256) :: domainDescription
    ! And we need the parameter list
    type(t_parlist) :: rparamlist

    ! Read out the area description in the string
    call parlst_getvalue_string(rcollection%p_rparlistQuickAccess1, 'DOMAININFO','onlySecond',domainDescription)

    call fparser_create(rparser,1)
    call fparser_parseFunction(rparser,1,domainDescription,(/'x','y'/))

    ! Now we can evaluate the area description with the parser-tools

    ! Get cubature weigths
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)
    p_Dfunc2 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_FUNC)

    ! We should have everything prepeared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coodinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0

           call fparser_evalFunction(rparser,1,(/dx,dy/),work)

           ! work can only be 0 or 1, so if it is < 0.5 it is 0
           if (abs(work) .le. 0.5_DP) then
                ! We dont need to work
                Dintvalue(1) = Dintvalue(1)
           else
                ! We need to work
                dval1 = p_Dfunc1(icubp,iel)**2
                dval2 = p_Dfunc2(icubp,iel)**2
                Dintvalue(1) = Dintvalue(1) + (dval1 + dval2)*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel


end subroutine

end module
