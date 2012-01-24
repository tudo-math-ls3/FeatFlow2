!##############################################################################
!# ****************************************************************************
!# <name> ccobject </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

module ccobject

  use ccbasic
  use basicgeometry
  use cubature
  use derivatives
  use dofmapping
  use element
  use fsystem
  use geometry
  use io
  use linearsystemscalar
  use storage
  use ucd

  implicit none
  
  public

contains

!--------------------------------------------------------------------------------------------------
!
    subroutine cc_cutoff(rgeometryObject,rmatrix,h_IelementList,listsize,nodes_in)   
!
!--------------------------------------------------------------------------------------------------
!
    type(t_geometryObject), pointer, intent(in)  ::  rgeometryObject    
    type(t_matrixscalar), intent(IN), target :: rmatrix 
    integer, intent(in), optional :: nodes_in
    integer, intent (out) :: h_IelementList
    integer, intent (out) :: listsize

 
! <local variables>
    integer :: i,j,icount,iin,iel,ivertices
    real(dp), dimension(:,:), pointer :: p_DvertexCoordinates   
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_Ielements
    real(dp), dimension(2) :: dcoords
! </local>    

    h_IelementList = ST_NOHANDLE
    
    call storage_new('cc_cutoff', 'h_IelementList', rmatrix%p_rspatialdiscrTrial%p_rtriangulation%NEL, &
                     ST_INT, h_IelementList,ST_NEWBLOCK_NOINIT)
    call storage_getbase_int (h_IelementList, p_Ielements)
        
    listsize = 0
    iel = rmatrix%p_rspatialDiscrTrial%p_rtriangulation%NEL
    do i=1,iel
      call storage_getbase_double2d (rmatrix%p_rspatialDiscrTrial%p_rtriangulation%h_DvertexCoords, &
                                     p_DvertexCoordinates)
      call storage_getbase_int2d (rmatrix%p_rspatialDiscrTrial%p_rtriangulation%h_IverticesAtElement,&
                                 p_IverticesAtElement)
    ! Check how many vertices of the element are inside the particle 
      icount = 0
      ivertices = rmatrix%p_rspatialDiscrTrial%p_rtriangulation%NNEE
      do j=1,ivertices
        dcoords(1) = p_DvertexCoordinates(1,p_IverticesAtElement(j,i))
        dcoords(2) = p_DvertexCoordinates(2,p_IverticesAtElement(j,i))                           

        call geom_isInGeometry (rgeometryObject,dcoords, iin)
        if (iin .eq. 1) then
          icount=icount+1
        end if    
      end do
      
      select case (nodes_in)
      
      case (1) ! 0 or 4 nodes inside for normal cubature formula
        if ((icount.eq.0).or.(icount .eq. 4))then
          p_Ielements(listsize+1)=i
          listsize = listsize+1
        end if
        
      case default ! 1-3 nodes inside for adaptive cubature formula
        if ((icount.gt.0).and.(icount .lt. 4))then
          p_Ielements(listsize+1)=i
          listsize = listsize+1
        end if  
        
      end select  
    end do

    if (listsize .gt. 0) then 
      call storage_realloc ("cc_cutoff", listsize, h_IelementList,ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase_int(h_IelementList,p_Ielements)
    end if
  end subroutine

!----------------------------------------------------------------------------------------------------------------------  

!--------------------------------------------------------------------------------------------------
!
    subroutine cc_calculateForces(rsolution,rproblem)   
!
!--------------------------------------------------------------------------------------------------

!<input>
  ! Solution vector to compute the norm/error from.
  type(t_vectorBlock), intent(in), target :: rsolution
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!<local variables>
  type(t_collection) :: rcollection
  ! pointer to the triangulation structure
  type(t_triangulation), pointer :: p_rtriangulation
  ! pointer to the particle collection
  type(t_particleCollection), pointer :: p_rparticleCollection
  ! Object/s geometry
  type(t_geometryObject), pointer :: p_rgeometry
  type(t_geometryObject) :: rgeometryObject1
  ! array of the elements which are cutted by the penalty object
  integer, dimension(:), pointer :: p_conectElements
  ! Vertices at Element
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  ! Coordinates of the vertices
  real(DP), dimension(:,:),pointer :: p_DvertexCoordinates

  integer :: h_IList,h_points,h_numpoints,h_pointsAtElement,iunit,cflag
  integer :: i,j,k,l,iel,ivert,iin,icount,ielind,inum,imaxdim,ipolyhandle,idfl
             
  real(dp) :: dxcenter,dycenter,dradius,dm,da,db,dc,ddiscr,dxsol1,dxsol2,dysol1,dysol2, &
              ddist1,ddist2,dmin,dmax,dabsmin,dabsmax,dpcont,dut,dbdForcesCoeff1,dbdForcesCoeff2
  real(dp) :: px1,px2,py1,py2,pxm,pym,dlh,dlen
  logical :: bfileExists
  real(dp), dimension(NDIM2D) :: dtangential, dnormal
  real(dp) :: dtol
  real(dp), dimension(4) :: DX, DY
  real(dp),dimension(2) :: Dcoord
  real(dp),dimension(2,2) :: locarray
  real(DP), dimension(:,:), pointer :: p_Dvertices
  real(dp), dimension(:,:), pointer :: p_Coords
  integer, dimension(:,:), pointer :: p_pointsAtElement
  integer, dimension(:), pointer :: p_points
  real(dp) :: dfw, daw
  real, dimension(2) :: dpf
  ! Pointer to vector data of solution vector
  real(DP), dimension(:), pointer :: p_DdataUX,p_DdataUY,p_DdataP

  ! 1D cubature formula identifier to use for the line integration
  integer(I32) :: ccub
  ! number/index of cubature points on the reference element
  integer :: ncubp,icubp
  ! Spatial discretisation structure of velocity and pressure
  type(t_spatialDiscretisation), pointer :: p_rdiscrU, p_rdiscrP
  ! Element type identifier for U and P
  integer(I32) :: ielemU, ielemP
  ! Number of local DOF's in U and P
  integer :: idoflocU, idoflocP
  ! An accepting the DOF's of an element
  integer, dimension(EL_MAXNBAS) :: IdofsU, IdofsP
  ! Coordinate system for U and P element
  integer(I32) :: ctrafoU, ctrafoP
  ! Array to tell the element which derivatives to calculate
  logical, dimension(EL_MAXNDER) :: BderU, BderP
  ! Cubature point coordinates on the reference element
  real(dp), dimension(CUB_MAXCUBP,NDIM3D) :: DXi1D, DXi2D
  ! For every cubature point on the reference element, 
  ! the corresponding cubature weight
  real(dp), dimension(CUB_MAXCUBP) :: Domega
  ! Coordinates of the cubature points on reference and real element
  real(dp), dimension(NDIM2D,CUB_MAXCUBP_1D) :: DpointsRef, DpointsReal
  ! Coordinates of the vertices of an element
  real(dp), dimension(2,4) :: Dcoords
  ! Arrays for saving Jacobian determinants and matrices
  real(dp), dimension(CUB_MAXCUBP_1D) :: Ddetj
  real(dp), dimension(EL_NJACENTRIES2D,CUB_MAXCUBP_1D) :: Djac
  ! Value of basis functions
  real(dp), dimension(El_MAXNBAS,EL_MAXNDER,CUB_MAXCUBP_1D) :: DbasU, DbasP
  real(dp), dimension(2) :: DintU, DintP

  type(t_ucdExport) :: rexport
  character(SYS_STRLEN) :: sfile,sfilepoly,sstr,stemp
!</local variables>

  h_IList = ST_NOHANDLE
  h_points = ST_NOHANDLE
  h_numpoints = ST_NOHANDLE
  h_pointsAtElement = ST_NOHANDLE
  dtol = 10e-7

  ! Initialisation
  dfw = 0.0_dp
  daw = 0.0_dp
  dlen = 0.0_dp

  call parlst_getvalue_double (rproblem%rparamlist,'CC-POSTPROCESSING','DBDFORCESCOEFF1',dbdForcesCoeff1,1.0_dp)
  call parlst_getvalue_double (rproblem%rparamlist,'CC-POSTPROCESSING','DBDFORCESCOEFF2',dbdForcesCoeff2,2.0_dp)

  dpf(1) = dbdForcesCoeff1
  dpf(2) = dbdForcesCoeff2


! Get U and P solution
call lsyssc_getbase_double (rsolution%RvectorBlock(1),p_DdataUX)
call lsyssc_getbase_double (rsolution%RvectorBlock(2),p_DdataUY)
call lsyssc_getbase_double (rsolution%RvectorBlock(3),p_DdataP)

! Cubature formula and number of cubature points to evaluate the integral
ccub = CUB_G1_1D
ncubp = 1

! Get the triangulation
  p_rtriangulation => rsolution%p_rblockdiscr%p_rtriangulation

  call storage_new('cc_calculateForces', 'h_IList', p_rtriangulation%NEL, ST_INT, &
                   h_IList,ST_NEWBLOCK_NOINIT)
  call storage_getbase_int (h_IList, p_conectElements)

! Get the penalty object/s
  if (rproblem%rparticlecollection%nparticles .gt. 0) then
    p_rparticleCollection => collct_getvalue_particles(rproblem%rcollection,'particles')
  end if

! 1.) Select all elements which interfears with the penalty object/s
  imaxdim = p_rtriangulation%NEL

! Coordinates and vertices at element
  call storage_getbase_double2d (p_rtriangulation%h_dvertexcoords,p_dvertexcoordinates)
  call storage_getbase_int2d (p_rtriangulation%h_iverticesatelement,p_iverticesatelement)  

  p_rgeometry => p_rparticleCollection%p_rparticles(1)%rgeometryobject
  ! Count the no. of elements which interfears with the particle.
  ielind = 1

  do iel = 1,imaxdim
    icount = 0
    Dcoords(1:2,1:4) = p_DvertexCoordinates(1:2,p_IverticesAtElement (1:4,iel))

    do ivert = 1,p_rtriangulation%NNVE
      Dcoord(1) = Dcoords(1,ivert)
      Dcoord(2) = Dcoords(2,ivert)
      call geom_isInGeometry (p_rgeometry, Dcoord, iin)
      icount = icount + iin
    end do

    if ((icount .gt. 0).and.(icount .lt. 4)) then
      p_conectElements(ielind) = iel
      ielind = ielind+1
    end if
    
  end do ! iel
  
  call storage_realloc ("cc_calculateForces", ielind-1, h_IList,ST_NEWBLOCK_NOINIT, .true.)
  call storage_getbase_int(h_IList,p_conectElements) 

  call intersection_points (ielind-1,p_rtriangulation,p_rgeometry,p_conectElements,p_DvertexCoordinates, &
                            p_IverticesAtElement,h_points,h_numpoints,h_pointsAtElement) 
                            
  call storage_getbase_int (h_numpoints,p_points)                             
  call storage_getbase_int2D (h_pointsAtElement,p_pointsAtElement)                             
  call storage_getbase_double2D (h_points,p_coords)                             

  ! Get U and P solution values
  call lsyssc_getbase_double (rsolution%RvectorBlock(1),p_DdataUX)
  call lsyssc_getbase_double (rsolution%RvectorBlock(2),p_DdataUY)
  call lsyssc_getbase_double (rsolution%RvectorBlock(3),p_DdataP)


  ! Get pointers to the spatial discretisation structures of the
  ! velocity and pressure
  p_rdiscrU => rsolution%RvectorBlock(1)%p_rspatialDiscr
  p_rdiscrP => rsolution%RvectorBlock(3)%p_rspatialDiscr
  ! What is the actual element that is used for the discretisation?
  ielemU = p_rdiscrU%RelementDistr(1)%celement
  ielemP = p_rdiscrP%RelementDistr(1)%celement
  ! Number of local DOF's
  idoflocU = elem_igetNDofLoc(ielemU)
  idoflocP = elem_igetNDofLoc(ielemP)
  ! Derivatves to calculate when evaluating the U and P-element, respectivly
  BderU = .false.
  BderU(DER_DERIV_X) = .true.
  BderU(DER_DERIV_Y) = .true.
  BderP = .false.
  BderP(DER_FUNC) = .true.

  ! Calculate drag and lift
  do i = 1,p_rtriangulation%NEL
    if (p_points(i) .ne. 0) then

      px1 = p_coords(1,p_pointsAtElement(1,i))
      px2 = p_coords(1,p_pointsAtElement(2,i))
      py1 = p_coords(2,p_pointsAtElement(1,i))
      py2 = p_coords(2,p_pointsAtElement(2,i))

      pxm = 0.5_dp*(px1+px2)
      pym = 0.5_dp*(py1+py2)
      dlh = sqrt((px2-px1)**2+(py2-py1)**2)

      dtangential(1) = (px2-px1)/dlh
      dtangential(2) = (py2-py1)/dlh
      dnormal(1) = -dtangential(2)
      dnormal(2) = dtangential(1)
      dlen = dlen + dlh

     ! For the integration, we need the global DOF`s on our element
     ! for U and P:
     call dof_locGlobMapping(p_rdiscrU, i, IdofsU)
     call dof_locGlobMapping(p_rdiscrP, i, IdofsP)

      ! The coordinates of the vertices on current element
      Dcoords (1:2,1:4) = p_DvertexCoordinates(1:2, p_IverticesAtElement (1:4,i))
      ! Real coordinates of the quadrature points
      DpointsReal(1,1) = pxm
      DpointsReal(2,1) = pym

      ! Calculate U and P in all cubature points
      call elem_generic_mult (ielemU, Dcoords, Djac, Ddetj, BderU, DbasU, ncubp, DpointsReal)
      call elem_generic_mult (ielemP, Dcoords, Djac, Ddetj, BderP, DbasP, ncubp, DpointsReal)

      do icubp = 1,ncubp
        dut = 0.0_DP
        do idfl=1,idoflocU
          dut = dut &
              + p_DdataUX(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_X,icubp)*Dtangential(1) &
                                        * Dnormal(1) &
              + p_DdataUY(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_X,icubp)*Dtangential(2) &
                                        * Dnormal(1) &
              + p_DdataUX(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_Y,icubp)*Dtangential(1) &
                                        * Dnormal(2) &
              + p_DdataUY(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_Y,icubp)*Dtangential(2) &
                                        * Dnormal(2)
        end do

        dpcont = 0.0_DP
        do idfl=1,idoflocP
           dpcont = dpcont +p_DdataP(IdofsP(idfl))*DbasP(idfl,DER_FUNC,icubp)
        end do
      end do

      DintU(1) = DintU(1) + dlh * dut * dnormal(2)
      DintU(2) = DintU(2) - dlh * dut * dnormal(1)
      DintP(1) = DintP(1) - dlh * dpcont * dnormal(1)
      DintP(2) = DintP(2) - dlh * dpcont * dnormal(2)
      
      dfw = dfw + dlh*(dpf(1)*dut*dnormal(2)-dpcont*dnormal(1))
      daw = daw - dlh*(dpf(1)*dut*dnormal(1)+dpcont*dnormal(2))

    end if
  end do

  dfw = 2D0*DFW/DPF(2)
  daw = 2D0*DAW/DPF(2)
  call parlst_getvalue_string (rproblem%rparamlist,'CC-POSTPROCESSING','sfilenameBodyForces',sfile,'')
  read(sfile,*) sfile
  call io_openFileForWriting(sfile,iunit,SYS_REPLACE,bfileExists,.true.)
      write (iunit,ADVANCE='YES',FMT='(A)') trim(sys_sdEL(dfw,5)) // ' ' // trim(sys_sdEL(daw,5))
  close (iunit)

  sfile="gmv/u.gmv.poly"
  call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)
    call geom_polygonise(p_rgeometry,ipolyHandle)
    ! Get the vertices
    call storage_getbase_double2D(ipolyHandle, p_Dvertices)
    call ucd_addPolygon(rexport,p_Dvertices,9)
    call storage_free(ipolyHandle)

    do i = 1, imaxdim
      if (p_points(i) .eq. 2) then
        locarray(1,1)=p_coords(1,p_pointsAtElement(1,i))
        locarray(2,1)=p_coords(2,p_pointsAtElement(1,i))
        locarray(1,2)=p_coords(1,p_pointsAtElement(2,i))
        locarray(2,2)=p_coords(2,p_pointsAtElement(2,i))
        call geom_init_polygon(rgeometryObject1,locarray)
        call geom_polygonise(rgeometryObject1,ipolyHandle)
        call storage_getbase_double2D(ipolyHandle, p_Dvertices)
        call ucd_addPolygon(rexport,p_Dvertices,10)
        call storage_free(ipolyHandle)
      else if (p_points(i) .gt. 2) then
        icount = p_points(i)
        do j=1,icount-1
          locarray(1,1)=p_coords(1,p_pointsAtElement(j,i))
          locarray(2,1)=p_coords(2,p_pointsAtElement(j,i))
          locarray(1,2)=p_coords(1,p_pointsAtElement(j+1,i))
          locarray(2,2)=p_coords(2,p_pointsAtElement(j+1,i))
          call geom_init_polygon(rgeometryObject1,locarray)
          call geom_polygonise(rgeometryObject1,ipolyHandle)
          call storage_getbase_double2D(ipolyHandle, p_Dvertices)
          call ucd_addPolygon(rexport,p_Dvertices,10)
          call storage_free(ipolyHandle)
        end do

      end if
    end do

    call ucd_write (rexport)
    call ucd_release (rexport)
    
    call storage_free(h_IList)
    call storage_free(h_points)
    call storage_free(h_numpoints)
    call storage_free(h_pointsAtElement)

  end subroutine

!----------------------------------------------------------------------------------------------------------------------  
!
  subroutine intersection_points (ielements,p_rtriangulation,p_rgeometry,p_Elements,p_DvertexCoordinates, &
                                  p_IverticesAtElement,h_points,h_numpoints,h_pointsAtElement)
!
!----------------------------------------------------------------------------------------------------------------------  
 !<input>
  integer, intent(in) :: ielements
  type(t_triangulation), pointer :: p_rtriangulation
  type(t_geometryObject), pointer, intent(in) :: p_rgeometry
  integer, dimension(:), pointer, intent(in) :: p_Elements
! Vertices at Element
   integer, dimension(:,:), pointer, intent(in) :: p_IverticesAtElement
! Coordinates of the vertices
   real(DP), dimension(:,:),pointer, intent(in) :: p_DvertexCoordinates
 !</input>

 !<output>
  integer, intent(out) :: h_points,h_numpoints,h_pointsAtElement
 !</output>

 !<local>
  integer :: inum,i,j,k,l,z,iel,itest,isize,index,index_points,total_points,iin1,coeff, &
             isemiplan1,isemiplan2,icase
  integer, dimension(2) :: Isiz
  logical :: bcalc
  real(dp) :: ddist1,ddist2,dtol,dxsol1,dxsol2,dysol1,dysol2,dm,da,db,dc,ddiscr, &
              dmax,dmin,dradius,dxcenter,dycenter,test_value,iin,teta0,dchange1,dchange2
  real(dp),dimension(4) :: Dcoord
  real(dp),dimension(2,20), target :: polyPoints, locarray
  real(dp), dimension(:,:), pointer :: p_Coords
  integer, dimension(:,:), pointer :: p_pointsAtElement
  integer, dimension(:), pointer :: p_points
  real(dp), dimension(:),allocatable :: teta
  integer, dimension(:), allocatable :: index_value

  !</local>

  h_points =  ST_NOHANDLE
  h_numpoints =  ST_NOHANDLE
  h_pointsAtElement = ST_NOHANDLE

  dxcenter = p_rgeometry%rcoord2d%dorigin(1)
  dycenter = p_rgeometry%rcoord2d%dorigin(1)
  dradius = p_rgeometry%rcircle%dradius

  Isiz(1) = 2 ! (4 points maximum)
  Isiz(2) = 4*p_rtriangulation%NEL

  call storage_new('intersection_points', 'h_points', Isiz, ST_DOUBLE, h_points,ST_NEWBLOCK_ZERO)
  call storage_getbase_double2D (h_points, p_Coords)

  ! Pointer to the indexes of the intersection points at each element (maximum 4 points)
  Isiz(1) = 4 
  Isiz(2) = p_rtriangulation%NEL 
  call storage_new('intersection_points', 'h_pointsAtElement', Isiz, ST_INT, h_pointsAtElement,ST_NEWBLOCK_ZERO)
  call storage_getbase_int2D (h_pointsAtElement, p_pointsAtelement)
        
  Isiz(2) = p_rtriangulation%NEL
  call storage_new('intersection_points', 'h_numpoints', Isiz(2), ST_INT, h_numpoints,ST_NEWBLOCK_ZERO)
  call storage_getbase_int (h_numpoints, p_points)

! Tolernace for the very close solution to vertices.
  dtol = 10e-7
! Indexing the intersection points
  index_points = 1
! Calculate for each element the intersection points
  do i = 1,p_rtriangulation%NEL
    iel = i
    bcalc = .FALSE.
    do j=1,ielements
      if (iel .eq. p_Elements(j)) bcalc = .TRUE.
!      if ((iel .eq. 130)) bcalc = .TRUE.
    end do
  if(bcalc) then
    inum = 0
    !Init intersection points array
    polypoints(:,:) = 0.0_dp

    do j = 1,4
      if (j .eq. 4) then
        k = 1
      else
        k = j+1
      end if
      ! Coordonatele dreptei:
      Dcoord(1) = p_DvertexCoordinates(1,p_IverticesAtElement(j,iel))
      Dcoord(2) = p_DvertexCoordinates(2,p_IverticesAtElement(j,iel))
      Dcoord(3) = p_DvertexCoordinates(1,p_IverticesAtElement(k,iel))
      Dcoord(4) = p_DvertexCoordinates(2,p_IverticesAtElement(k,iel))
      ddist1 = sqrt((dxcenter-Dcoord(1))**2+(dycenter-Dcoord(2))**2)
      ddist2 = sqrt((dxcenter-Dcoord(3))**2+(dycenter-Dcoord(4))**2)

      if ((ddist1 .le. dradius).or.(ddist2 .le. dradius)) then
        if (abs(Dcoord(3)-Dcoord(1)) .lt. dtol) then
          ! vertical edge
          dxsol1 = Dcoord(1)
          dysol1 = dycenter + sqrt(dradius**2 -(dxsol1-dxcenter)**2)
          dysol2 = dycenter - sqrt(dradius**2 -(dxsol1-dxcenter)**2)
          ! check solution
          dmin = min(Dcoord(2),Dcoord(4)) - dtol
          dmax = max(Dcoord(2),Dcoord(4)) + dtol
          call geom_isInGeometry (p_rgeometry,(/dxsol1,dysol1/), iin1)
          if ((dysol1 .gt. dmin).and.(dysol1 .lt. dmax).and.(iin1 .eq. 1)) then
            if (inum .eq. 0) then
              inum = inum + 1
              polyPoints(1,inum)=dxsol1
              polyPoints(2,inum)=dysol1
            else 
              itest = 0
              do l = 1,inum
                if (abs(dysol1-polyPoints(2,l)) .gt. dtol*(abs(dysol1)+abs(polyPoints(2,l)))) itest = itest + 1 
              end do
              if (itest .eq. inum) then
                inum = inum + 1
                polyPoints(1,inum)=dxsol1
                polyPoints(2,inum)=dysol1
              end if
            end if
          end if

          call geom_isInGeometry (p_rgeometry,(/dxsol1,dysol2/), iin1)
          if ((dysol2 .gt. dmin).and.(dysol2 .lt. dmax).and.(iin1 .eq. 1)) then
            if (inum .eq. 0) then
              inum = inum + 1
              polyPoints(1,inum)=dxsol1
              polyPoints(2,inum)=dysol2
            else 
              itest = 0
              do l = 1,inum
                if (abs(dysol2-polyPoints(2,l)) .gt. dtol*(abs(dysol2)+abs(polyPoints(2,l)))) itest = itest + 1 
              end do
              if (itest .eq. inum) then
                inum = inum + 1
                polyPoints(1,inum)=dxsol1
                polyPoints(2,inum)=dysol2
              end if
            end if
          end if
        else  
          dm = (Dcoord(4)-Dcoord(2))/(Dcoord(3)-Dcoord(1))
          da = 1+dm**2
          db = -2*dxcenter-2*dm**2*Dcoord(1)+2*dm*(Dcoord(2)-dycenter)
          dc = dxcenter**2+dm**2*Dcoord(1)**2-2*dm*Dcoord(1)*(Dcoord(2)-dycenter)+(Dcoord(2)-dycenter)**2-dradius**2
          ddiscr = db**2-4*da*dc
          if (ddiscr .ge. 0) then
            dxsol1 = (-db+sqrt(ddiscr))/(2*da)
            dysol1 = dm*(dxsol1-Dcoord(1))+Dcoord(2)
            dxsol2 = (-db-sqrt(ddiscr))/(2*da)
            dysol2 = dm*(dxsol2-Dcoord(1))+Dcoord(2)

            if(sqrt((dxsol1-dxsol2)**2+(dysol1-dysol2)**2) .ge. &
               sqrt((Dcoord(1)-Dcoord(3))**2+(Dcoord(2)-Dcoord(4))**2)) then

            ! Check solutions
            dmin = min(Dcoord(1),Dcoord(3)) - dtol
            dmax = max(Dcoord(1),Dcoord(3)) + dtol
              
            if ((dxsol1 .gt. dmin) .and. (dxsol1 .lt. dmax)) then
              if (inum .eq. 0)then
                inum = inum + 1
                polyPoints(1,inum)=dxsol1
                polyPoints(2,inum)=dysol1
              else 
                itest = 0
                do l = 1,inum
                  if (abs(dxsol1-polyPoints(1,l)) .gt. dtol*(abs(dxsol1)+abs(polyPoints(1,l)))) itest = itest + 1 
                end do
                if (itest .eq. inum) then
                  inum = inum + 1
                  polyPoints(1,inum)=dxsol1
                  polyPoints(2,inum)=dysol1
                end if
              end if
            end if

            if ((dxsol2 .gt. dmin) .and. (dxsol2 .lt. dmax)) then
              if (inum .eq. 0) then
                inum = inum + 1
                polyPoints(1,inum)=dxsol2
                polyPoints(2,inum)=dysol2
              else 
                itest = 0
                do l = 1,inum
                  if (abs(dxsol2-polyPoints(1,l)) .gt. dtol*(abs(dxsol2)+abs(polyPoints(1,j)))) itest = itest + 1
                end do
                 if (itest .eq. inum) then
                  inum = inum + 1
                  polyPoints(1,inum)=dxsol2
                  polyPoints(2,inum)=dysol2
                end if
              end if
            end if
            end if
          end if
        end if
      end if
    end do ! j

    if (inum .gt. 1) then
      index = i
      p_points(index) = inum
      total_points = sum(p_points(:))
      if (total_points .gt. 0) then
        call storage_realloc ("intersection_points", total_points, h_points,ST_NEWBLOCK_ZERO, .true.)
        call storage_getbase_double2D (h_points,p_coords)
      end if

      if (inum .ge. 2) then
!        allocate(teta(inum))
!        allocate(index_value(inum))
        isemiplan1 = 1
        isemiplan2 = 1
        if (polypoints(2,1)-dycenter .lt. 0.0_dp) isemiplan1 = 2
        if (polypoints(2,2)-dycenter .lt. 0.0_dp) isemiplan2 = 2

        icase = isemiplan1 + isemiplan2

        select case (icase)
          case (2) ! both points in the positive halfplane
            if (polypoints(1,1) .gt. polypoints(1,2)) then 
               dchange1 = polypoints(1,1)
               dchange2 = polypoints(2,1)
               polypoints(1,1) = polypoints(1,2)
               polypoints(2,1) = polypoints(2,2)
               polypoints(1,2) = dchange1
               polypoints(2,2) = dchange2
            end if     
          case (4) ! both points in the negative halfplane
            if (polypoints(1,1) .lt. polypoints(1,2)) then 
               dchange1 = polypoints(1,1)
               dchange2 = polypoints(2,1)
               polypoints(1,1) = polypoints(1,2)
               polypoints(2,1) = polypoints(2,2)
               polypoints(1,2) = dchange1
               polypoints(2,2) = dchange2
             end if     
          case (3) ! one point in the positive halfplane, one in the negative
             if (polypoints(1,1)-polypoints(1,2) .lt. 0.0_dp) then
               if (polypoints(2,1) .gt. polypoints(2,2)) then
                 dchange1 = polypoints(1,1)
                 dchange2 = polypoints(2,1)
                 polypoints(1,1) = polypoints(1,2)
                 polypoints(2,1) = polypoints(2,2)
                 polypoints(1,2) = dchange1
                 polypoints(2,2) = dchange2
               end if
             else
               if (polypoints(2,1) .lt. polypoints(2,2)) then
                 dchange1 = polypoints(1,1)
                 dchange2 = polypoints(2,1)
                 polypoints(1,1) = polypoints(1,2)
                 polypoints(2,1) = polypoints(2,2)
                 polypoints(1,2) = dchange1
                 polypoints(2,2) = dchange2
               end if
             end if
        end select

!        do j = 1,inum 
!          if ((polypoints(1,j) .gt. dxcenter) .and. (polypoints(2,j) .gt. dycenter)) then
!            teta0 = 0.0_dp
!          end if
!          if ((polypoints(1,j) .lt. dxcenter) .and. (polypoints(2,j) .gt. dycenter)) then
!             teta0 = 90.0_dp
!          end if
!          if ((polypoints(1,j) .lt. dxcenter) .and. (polypoints(2,j) .lt. dycenter)) then
!            teta0 = 180.0_dp
!          end if
!          if ((polypoints(1,j) .gt. dxcenter) .and. (polypoints(2,j) .lt. dycenter)) then
!            teta0 = 270.0_dp
!          end if
!          
!          write(*,*) 180.0_dp*asin(((polypoints(2,j)-dycenter)/dradius))/SYS_PI
!          teta(j) = teta0 + 180.0_dp*asin(abs(polypoints(2,j)-dycenter)/dradius)/SYS_PI
!        end do
!
!        dmin = minval(teta)
!        dmax = maxval(teta)
!        do j=1,inum
!          if (teta(j) .eq. dmin) then
!            test_value = j
!            locarray(1,2) = polypoints(1,test_value) 
!            locarray(2,2) = polypoints(2,test_value)
!          end if
!          if (teta(j) .eq. dmax) then
!            test_value = j
!            locarray(1,1) = polypoints(1,test_value) 
!            locarray(2,1) = polypoints(2,test_value)
!          end if
!        end do
        
        if (inum .gt. 2) inum =  inum - 1
!        polypoints(:,:) = locarray(:,:)
        p_points(index) = inum
      
!        deallocate(index_value)
!        deallocate(teta)
      end if
      do j=1,inum
        p_pointsAtElement(j,iel) = index_points
        p_coords(1,p_pointsAtElement(j,iel)) = polypoints(1,j)
        p_coords(2,p_pointsAtElement(j,iel)) = polypoints(2,j)
        index_points = index_points+1
      end do
    end if
  end if ! bcalc
  end do ! i

 end subroutine
!----------------------------------------------------------------------------------------------------------------------  

end module
