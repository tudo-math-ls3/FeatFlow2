************************************************************************
* This file describes the handle parameter list that is used in the
* main program to pass parameters from one subroutine to another.
* This handle parameter list is an array of handles (integers)
* to the different pieces of information as well as information
* about the total size of the system and more:
*
*  HLIST : array [1..*] of integer
*          = THList
*
* The entries are defined as follows:
*
* THList : record
*   
*   // Triangulation and matrix/vector data
*
*   LTRIAS   : handle to      INTEGER TRIAS(SZTRIA,NNLEV)
*   LCURTR   : handle to      INTEGER CURTRI(SZTRIA,NNLEV)
*              Currently used triangulation, copy of TRIAS,
*              probably with grid deformation
*   LMATDAT  : handle to      INTEGER MATDAT(SZN2MI,NNLEV)
*   LVECDAT  : handle to      INTEGER VECDAT(SZN2VI,NNLEV)
*
*   // Parameters of the stationary solver
*
*   LISTPAR  : handle to      INTEGER ISTPAR(SZNSDI)
*   LDSTPAR  : handle to      DOUBLE PRECISION DSTPAR(SZNSDD)
*    
*   // Multigrid structures provide space for one MG-Solver, Smoother
*   // and coarse grid solver:
*    
*   LIMGPAR  : handle to      INTEGER IMGPAR (SZ020I+SZSLVI+SZSLVI+SZSLVI)
*   LDMGPAR  : handle to      DOUBLE PRECISION DMGPAR (SZ020D+SZSLVD+SZSLVD+SZSLVD)
*    
*   // Assembly information
*
*   LIASMBL  : handle to      INTEGER IASMBL(SZASMI)
*   LDASMBL  : handle to      DOUBLE PRECISION DASMBL(SZASMD)
*
*   // Geometry information
*    
*   LIGEOM    : handle to     INTEGER IGEOM(*)
*   LDGEOM    : handle to     DOUBLE PRECISION DGEOM(*)
*
*   // Vector information:
*
*   LUP       : handle to     DOUBLE PRECISION DUP(*)
*   LRHS      : handle to     DOUBLE PRECISION DRHS(*)
*   LAUX      : handle to     DOUBLE PRECISION DAUX(*)
*
*   // Scalar information:
*
*   NEQU      : INTEGER NEQU
*   NEQP      : INTEGER NEQP
*   NUVP      : INTEGER NUVP: total length of solution vector
*   IBDTYP    : INTEGER IBDTYP
*
*   // General information
*
*   IFXOUT    : INTEGER IFXOUT
*               Handle to the file where the output is written to
*
*   LIOPTP    : Handle to the integer parameters of the optimizer
*   LDOPTP    : Handle to the double parameters of the optimizer
* end
*
* The access to the variables is done using the constants defined here,
* e.g.
*       LUP = HLIST(OLUP)
*
************************************************************************

      INTEGER OLTRIAS ,OLMATDAT,OLVECDAT,OLISTPAR,OLDSTPAR,OLIMGPAR
      INTEGER OLDMGPAR,OLIASMBL,OLDASMBL,OLIGEOM ,OLDGEOM ,ODUP    
      INTEGER ODRHS   ,ONUVP   ,OIFXOUT ,OLCURTR 
      INTEGER OLIOPTP ,OLDOPTP

      PARAMETER (OLTRIAS  = 1)
      PARAMETER (OLMATDAT = 2)
      PARAMETER (OLVECDAT = 3)
      PARAMETER (OLISTPAR = 4)
      PARAMETER (OLDSTPAR = 5)
      PARAMETER (OLIMGPAR = 6)
      PARAMETER (OLDMGPAR = 7)
      PARAMETER (OLIASMBL = 8)
      PARAMETER (OLDASMBL = 9)
      PARAMETER (OLIGEOM  = 10)
      PARAMETER (OLDGEOM  = 11)
      PARAMETER (ODUP     = 12)
      PARAMETER (ODRHS    = 13)
      PARAMETER (ONUVP    = 14)
      PARAMETER (OIFXOUT  = 15)
      PARAMETER (OLCURTR  = 16)
      PARAMETER (OLIOPTP  = 17)
      PARAMETER (OLDOPTP  = 18)

      INTEGER SZHLST
      
      PARAMETER (SZHLST  = 32)
      