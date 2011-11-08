************************************************************************
* This files describes the general parameter structure of the grid
* adaption routines.
*
* The grid adaption routines take two parameters: IPARAM and DPARAM.
* These parameter blocks are designed as structures and realized by
* arrays.
* The structures contain substructures for the different
* components of the method. Each structure and substructure contains
* - a set of input variables, which determines the behaviour of
*   the method and the sub-components
* - a set of output variables, which describe the status of
*   the sub-components
* - a set of temporary/unused/reserved variables that can be used
*   by the method
* The method tells its status by changing the output-variables to
* an appropriate state. The input-variables are not changed, so the
* method can be called again with the same parameters if necessary.
* The initial state of the output- and temporary variables can be
* undefined.
*
* As each grid adaption method uses different parameters, a different
* structure is necessary for each method. Nevertheless some parameters
* are global ones and valid for all methods. The global structures
* are defined here in this file. By putting the structures together,
* a "local" structure is created with the parameters necessary
* for the different routines.
*
* The realization of the structures as arrays is made e.g. by:
*
*   IPARAM = array [1..SZGACI] of integer
*          = TIGridAdaptionCommon
*   DPARAM = array [1..SZGACD] of double
*          = TDGridAdaptionCommon
*
* The Oxxxxx-constants below can be used to access the different
* variables inside of the structure-array, e.g.:
*
*   RESULT = IPARAM(ORESULT)
*
* Commonly we have the following two structures for integer- and
* double-variables:
*
* TIGridAdaptionCommon = record
*   // Output parameters:
*
*   RESULT : The result of the grid adaption.
*            =0: grid adaption successful
*            =1: error in the input parameters
*            =2: grid adaption broke down in the solution process
*                of a linear system
*            =3: grid adaption broke down in moving the grid points
*
*   NLSITE : Number of iterations necessary for solving all
*            linear systems on the finest level of the grid adaption
*   NLSYS  : Number of solved linear systems
*   NCODE  : Number of solved ODE subproblems
*   NMFEVL : Number of monitor function evaluations
*
*   // Input parameters:
*
*   NGASTP : Number of grid adaption (makro-)steps. >= 1, otherwise
*            there is no grid adaption performed.
*   IARCTR : Adaptive restart control
*            =0: No adaptive restart control, anways perform NGASTP
*            (macro-)steps
*   IMNREL : Whether to use a relative or absolute monitor function.
*            =0: Monitor function determines absolute cell size (standard)
*            =1: Monitor function determines relative cell size,
*                based on previous cell distribution
*   IGSOUT : Controls the level of output in the grid adaption process
*            to GMV-files. Bitfield.
*   ICOLMC : Whether to colorize fine-grid cells according to
*            coarse-grid cells for GMV-output
*            =0: disable coloring
*            =1: colorize macros
*
*   ISMMTH : Type of smoothing method to use for grid smoothing.
*            =0: Laplace
*            =1: weighted umbrella
*   NPRGSM : Number of grid pre-smoothing steps of laplacian grid
*            smoothing before each macrostep
*   NPOGSM : Number of grid post-smoothing steps of laplacian grid
*            smoothing after each macrostep
*   NFSMTH : Number of steps for smoothing of monitor function
*
*   ICUNGS : Cubature formula to use for normalizing of functions.
*            =-1: don't use cubature, but use quick-and-dirty
*                 weighting with mean area around vertices.
*
* end
*
* TDGridAdaptionCommon = record
*   // Output parameters:
*
*   TGTOT  : Total time necessary for grid adaption
*   TGGAPR : Total time necessary for preparation of grid adaption
*            (setting up monitor function, normalization,...)
*   TGLSPR : Total time necessary for setting up linear systems
*   TGLS   : Total time necessary for solving linear systems
*   TGODE  : Total time necessary for ODE solving + grid deformation
*   TGGSM  : Total time necessary for grid smoothing
*   TGGRC  : Total time necessary for gradient recovery
*   TGMSM  : Total time necessary for monitor function smoothing
*   TGGPR  : Total time for grid prolongation/restriction
*
*   // Input parameters:
*
*   DMNTOL : Tolerance of monitor function for automatic
*            determination of restarts
*            =0: always perform NGASTP grid adaption steps
*   EPS0GS : EPS0 of monitor function
*   DAPOSM : Adaptive determination of postsmoothing-steps
*            =0: no adaptive post-smoothing
*   DMNSFC : Monitor function scaling/regularisation factor
*            =0: no regularization
* end
*
************************************************************************

***********************
* TIGridAdaptionCommon
***********************
  
      INTEGER ORESULT,ONLSITE,ONLSYS,ONCODE,ONMFEVL,OISMMTH,OICOLMC
      INTEGER ONGASTP,OIARCTR,OIGSOUT,ONPRGSM, ONPOGSM, ONFSMTH
      INTEGER OIMNREL,OICUNGS
  
      PARAMETER (ORESULT= 1)
      PARAMETER (ONLSITE= 2)
      PARAMETER (ONLSYS = 3)
      PARAMETER (ONCODE = 4)
      PARAMETER (ONMFEVL= 5)
      
      PARAMETER (ONGASTP= 16)
      PARAMETER (OIARCTR= 17)
      PARAMETER (OIMNREL= 18)
      PARAMETER (OIGSOUT= 19)
      PARAMETER (OICOLMC= 20)
      PARAMETER (OISMMTH= 21)
      PARAMETER (ONPRGSM= 22)
      PARAMETER (ONPOGSM= 23)
      PARAMETER (ONFSMTH= 24)
      PARAMETER (OICUNGS= 25)
  
C  9..15: reserved
C 24..31: user-defineable
  
      INTEGER SZGACI
      PARAMETER (SZGACI=32)
      
***********************
* TDGridAdaptionCommon
***********************

      INTEGER ODMNTOL, OEPS0GS, ODAPOSM, ODMNSFC
      INTEGER OTGTOT,OTGLSPR,OTGLS,OTGODE,OTGGPR, OTGGSM, OTGMSM
      INTEGER OTGGAPR,OTGGRC
      
      PARAMETER (OTGTOT   = 1)
      PARAMETER (OTGGAPR  = 2)
      PARAMETER (OTGLSPR  = 3)
      PARAMETER (OTGLS    = 4)
      PARAMETER (OTGODE   = 5)
      PARAMETER (OTGGSM   = 6)
      PARAMETER (OTGGRC   = 7)
      PARAMETER (OTGMSM   = 8)
      PARAMETER (OTGGPR   = 9)

      PARAMETER (ODMNTOL  = 16)
      PARAMETER (OEPS0GS  = 17)
      PARAMETER (ODAPOSM  = 18)
      PARAMETER (ODMNSFC  = 19)
      
C  5..24: reserved
C 25..31: user-defineable
  
      INTEGER SZGACD
      PARAMETER (SZGACD=32)
