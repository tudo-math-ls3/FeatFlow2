!##############################################################################
!# ****************************************************************************
!# <name> discretefbc </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements discrete boundary conditions for fictitious boundary
!# components. Discrete boundary conditions for fictitious boundary components
!# are an element-dependent way to represent analytical 'boundary conditions'
!# in a general sense. After 'discretising', the boundary condition can quickly
!# be implemented into a vector. Therefore, one can see discrete boundary
!# conditions also as a kind of `precalculated boundary conditions of
!# fictitious boundary objects`. This is exploited e.g. in the filter approach,
!# where a filter routine 'applies' a discrete FBC to a vector.
!#
!# Usually, one fictitious boundary object may consist of multiple subobjects
!# like balls or general shapes) which all share the same boundary condition.
!# 'Discretising' this means collecting all the DOF`s that are affected by
!# such an object as well as specifying a strategy how to modify matrix,
!# solution and/or RHS vector.
!#
!# The basic implementation realises the implementation of objects that
!# specify a Dirichlet value in a subdomain of the global domain (e.g. a
!# ball with a specified velocity in a fluid). But this is not
!# the only application: One can even think of a subdomain that exhibits
!# stress to the right hand side of the system.
!#
!# </purpose>
!##############################################################################

module discretefbc

  use fsystem
  use storage
  use boundarycondition
  use dofmapping
  
  implicit none
  
  private

!<constants>

!<constantblock description="General constants concerning filters">

  ! A standard length for arrays holding a set of discretised BC`s
  integer, parameter, public :: DISCFBC_MAXDISCBC         =  32

!</constantblock>

!<constantblock description="The type identifier for discrete boundary conditions">

  ! undefined discrete BC`s
  integer, parameter, public :: DISCFBC_TPUNDEFINED    = 0

  ! Discrete Dirichlet boundary conditions
  integer, parameter, public :: DISCFBC_TPDIRICHLET    = 1
  
!</constantblock>

!<constantblock description="Type identifiers for the callback routine during discretisation of FBC`s">
  
  ! Calculate the function value in some point of the domain.
  integer, parameter, public :: DISCFBC_NEEDFUNCGENERAL    = 0
  
  ! Calculate the function value in corner vertices of elements in the object.
  integer, parameter, public :: DISCFBC_NEEDFUNC           = 1

  ! Calculate the function value in a edge midpoints of elements in the object
  integer, parameter, public :: DISCFBC_NEEDFUNCMID        = 2
  
  ! Calculate the integral mean value on edges in the object
  integer, parameter, public :: DISCFBC_NEEDINTMEAN        = 3

  ! Calculate the function value in midpoints of elements in the object
  integer, parameter, public :: DISCFBC_NEEDFUNCELMID      = 4

  ! Calculate the function value in the face midpoints of elements in the object
  integer, parameter, public :: DISCFBC_NEEDFUNCFACEMID    = 5
  
  ! Calculate the integral mean value on the faces in the object
  integer, parameter, public :: DISCFBC_NEEDFACEINTMEAN    = 6
  

!</constantblock>
  
!</constants>

!<types>
  
!<typeblock>
  
  ! This structure describes the typical way, Dirichlet boundary conditions
  ! for fictitious boundary components can be discretised.
  ! This is done by two arrays: one array is a list of all
  ! DOF`s that refer do Dirichlet nodes. The second array refers to the value
  ! that must be imposed in this DOF.
  ! The variable Icomponents receives a list of all components in the solution
  ! vector that are affected by this bonudary condition.
  
  type t_discreteFBCDirichlet
    
    ! Number of boundary components affected by this boundary condition.
    ! E.g. =2 for X- and Y-velocity.
    integer :: ncomponents = 0
    
    ! A list of 1..ncomponents components of the equation, this discrete BC
    ! is specified for (e.g. [1 2] = X-velocity(1) + Y-velocity(2))
    integer, dimension(:), pointer :: Icomponents => null()
    
    ! Number of Dirichlet nodes; may be different from the length of the array!
    integer :: nDOF = 0
    
    ! Handle to array with all DOF`s that refer to Dirichlet nodes
    !   array [1..*] of integer
    ! p_IdirichletDOFs(i) is the number of the i-th DOF that is to be overwritten
    ! by the 'Dirichlet replacement' filter.
    integer :: h_IdirichletDOFs   = ST_NOHANDLE
    
    ! Handle to array with the Dirichlet value that should be imposed in these nodes
    !   array [1..ncomponents,1..nDOF] of double
    integer :: h_DdirichletValues = ST_NOHANDLE
    
  end type
  
  public :: t_discreteFBCDirichlet
  
!</typeblock>
  
!<typeblock>
  
  ! This describes the basic structure for discrete boundary conditions
  ! on fictitious boundary components.
  ! A type identifier decides on which boundary conditions this structure
  ! describes. Depending on the type, one of the information blocks
  ! is filled with data about the discrete BC`s.
  
  type t_discreteFBCEntry
    
    ! The type identifier. Identifies the type of discrete BC`s, this
    ! structure describes.
    integer                             :: itype = DISCFBC_TPUNDEFINED
    
    ! Structure for discrete Dirichlet BC`s.
    ! Only valid if itype=DISCBC_TPDIRICHLET.
    type(t_discreteFBCDirichlet)        :: rdirichletFBCs
    
  end type
  
  public :: t_discreteFBCEntry
  
!</typeblock>

!<typeblock>

  ! The main structure for discrete boundary conditions of fictitious boundary
  ! objects.
  ! This is just an array of t_discreteFBCEntry structures, each describing
  ! a single discrete boundary condition (so to speak, a collection of objects
  ! sharing the same boundary conditions, discretised in the same way).
  type t_discreteFBC
  
    ! Total number of allocated t_discreteFBCEntry structures in p_RdiscFBCList.
    integer :: inumEntriesAlloc = 0
    
    ! Total number of used t_discreteFBCEntry structures in p_RdiscFBCList.
    integer :: inumEntriesUsed = 0

    ! An array of t_discreteFBCEntry structures. Each structure describes
    ! one discrete fictitions boundary boundary condition - so one group of
    ! objects with the same 'boundary' condition discretised
    ! in a special, discretisation-dependent way.
    type(t_discreteFBCEntry), dimension(:), pointer :: p_RdiscFBCList => null()
  
  end type
  
  public :: t_discreteFBC

!</typeblock>

!<typeblock>

  ! A structure of this type is passed to the callback routine for assembling the FB
  ! boundary conditions. This structure specifies the points/locations where the
  ! callback routine for the assembly should evaluate.
  type t_discreteFBCevaluation

    ! An integer tag that defines what to evaluate. One of the DISCFBC_xxxx constants.
    ! The constant in this tag defines what can be found in the other entries of
    ! this structure!
    integer :: cinfoNeeded = 0

    ! Number of values to calculate. What to calculate is depending on cinfoNeeded:
    ! cinfoNeeded = DISCFBC_NEEDFUNC,DISCFBC_NEEDFUNCGENERAL:
    !               nvalues = number of vertices (=length of p_Iwhere)
    !
    ! cinfoNeeded = DISCFBC_NEEDFUNCMID:
    !               nvalues = number of edges (=length of p_Iwhere)
    !
    ! cinfoNeeded = DISCFBC_NEEDINTMEAN:
    !               nvalues = number of edges (=length of p_Iwhere)
    !
    ! cinfoNeeded = DISCFBC_NEEDFUNCMID:
    !               nvalues = number of elements (=length of p_Iwhere)
    integer :: nvalues = 0

    ! An integer array specifying the location of points/edges/elements where
    ! to evaluate. The content is depending on the situation, specifíed by the
    ! cinfoNeeded parameter: \\
    !
    ! cinfoNeeded = DISCFBC_NEEDFUNCGENERAL: \\
    !   p_Iwhere(.) = 1..NEL -> Number of the element in the triangulation that
    !                           contains the point.
    !
    ! cinfoNeeded = DISCFBC_NEEDFUNC: \\
    !   p_Iwhere(.) = 1..NVT -> Number of the corner vertex in the triangulation
    !                           where to evaluate.\\
    !
    ! cinfoNeeded = DISCFBC_NEEDFUNCMID: \\
    !   p_Iwhere(.) = 1..NMT -> Number of the edge in the triangulation in whose
    !                           midpoint to evaluate.\\
    !
    ! cinfoNeeded = DISCFBC_NEEDINTMEAN: \\
    !   p_Iwhere(.) = 1..NMT -> Number of the edge in the triangulation where to
    !                           evaluate the integral mean value.\\
    !
    ! cinfoNeeded = DISCFBC_NEEDFUNCMID: \\
    !   p_Iwhere(.) = 1..NEL -> Number of the element in the triangulation in which
    !                           midpoint to evaluate.\\
    !
    ! Usually, nvalues defines the number of entries in p_Iwhere where to evaluate;
    ! the actual size of p_Iwhere might be larger than nvalues, so the callback
    ! routine should orientate on nvalues.
    integer, dimension(:), pointer :: p_Iwhere => null()
    
    ! Coordinate array that specifies the coordinate where to evaluate
    ! (if available). The content is depending on the situation, specifíed by the
    ! cinfoNeeded parameter: \\
    !
    ! cinfoNeeded = DISCFBC_NEEDFUNCGENERAL,
    !               DISCFBC_NEEDFUNC,
    !               DISCFBC_NEEDFUNCMID,
    !               DISCFBC_NEEDFUNCELMID,
    !               DISCFBC_NEEDFUNCFACEMID :
    !   p_Dwhere(1:ndim,.) = x/y/z-coordinate of the point where to evaluate
    !
    ! cinfoNeeded = DISCFBC_NEEDINTMEAN,
    !               DISCFBC_NEEDFACEINTMEAN :
    !   p_Dwhere(1:ndim,.) = x/y/z-coordinate of the midpoint of the edge/face
    !                        where to evaluate
    real(DP), dimension(:,:), pointer :: p_Dwhere => null()
  
    ! A pointer to an array that accepts calculated values. The callback routine
    ! for evaluating on the fictitious boundary component must fill this
    ! array with values. The content and shape is depending on cinfoNeeded:
    !
    ! cinfoNeeded = DISCFBC_NEEDFUNC: \\
    !   p_Dvalues(1..nvalues,1) -> Point values in the vertices defined in p_Iwhere \\
    !
    ! cinfoNeeded = DISCFBC_NEEDFUNCMID: \\
    !   p_Dvalues(1..nvalues,1) -> Point values in the edge midpoints defined
    !                              in p_Iwhere\\
    !
    ! cinfoNeeded = DISCFBC_NEEDINTMEAN: \\
    !   p_Dvalues(1..nvalues,1) -> Integral mean values in the edges defined
    !                              in p_Iwhere\\
    !
    ! cinfoNeeded = DISCFBC_NEEDFUNCELMID: \\
    !   p_Dvalues(1..nvalues,1) -> Point values in the midpoints of the elements
    !                              defined in p_Iwhere.\\
    !
    ! Remark: The second dimension is not used for now. In a later implementation
    ! when derivatives might be needed, the shape may be used in a different way,
    ! e.g. p_Dvalues(1..nvalues,1..#derivatives) or similar.
    !
    ! Remark: Depending on what is to be evaluated, not all values (point values,...)
    ! must be calculated. If some values are left out, the callback routine can
    ! use the p_Binside array to indicate what is calculated and what not.
    real(DP), dimension(:,:), pointer :: p_Dvalues => null()
    
    ! A pointer to an array of integer values for each value to calculate.
    ! Must be set to 1 by the callback routine for those values that are
    ! calculated. The meaning is again depending on cinfoNeeded:
    !
    ! cinfoNeeded = DISCFBC_NEEDFUNC: \\
    !   Set "p_Iinside(i) = 1" if vertex i is inside of the FBC object and
    !   p_Dvalues(i,1) is calculated.
    !   Set "p_Iinside(i) = 0" if vertex i is outside of the FBC object;
    !   p_Dvalues(i,1) need not to be calculated.\\
    !
    ! cinfoNeeded = DISCFBC_NEEDFUNCMID: \\
    !   Set "p_Iinside(i) = 1" if midpoint of edge i is inside of the FBC object
    !   and p_Dvalues(i,1) is calculated.
    !   Set "p_Iinside(i) = 0" if midpoint of edge i is outside of the FBC object;
    !   p_Dvalues(i,1) need not to be calculated.\\
    !
    ! cinfoNeeded = DISCFBC_NEEDINTMEAN: \\
    !   Set "p_Iinside(i) = 1" if the integral mean value of edge i is inside
    !   of the FBC object and p_Dvalues(i,1) is calculated.
    !   Set "p_Iinside(i) = 0" if the integral mean value of edge i is outside
    !   of the FBC object; p_Dvalues(i,1) need not to be calculated.\\
    !
    ! cinfoNeeded = DISCFBC_NEEDFUNCELMID: \\
    !   Set "p_Iinside(i) = 1" if the midpoint of element i is inside of the FBC
    !   object and p_Dvalues(i,1) is calculated.
    !   Set "p_Iinside(i) = 0" if the midpoint of element i is outside of the FBC
    !   object; p_Dvalues(i,1) need not to be calculated.
    integer, dimension(:), pointer :: p_Iinside
  
  end type
  
  public :: t_discreteFBCevaluation

!</typeblock>

!</types>
  
end module
