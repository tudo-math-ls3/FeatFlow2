************************************************************************
* This files describes the parameter structure of the static grid
* adaption routine XGASTA.
*
* The routine takes two parameters: IPARAM and DPARAM. These
* parameter blocks are designed as structures and realized by arrays.
*
*   IPARAM = array [1..SZGASI] of integer 
*          = TIGridAdaptionStaticParams
*   DPARAM = array [1..SZGASI] of double
*          = TDGridAdaptionStaticParams
*
* The two structures are organized as follows: At first there
* is a block for the general grid adaption parameters as well as
* a block for the solver parameters. The method specific information
* is appended to the structure after these both blocks.
*
* TIGridAdaptionStaticParams = record
*  
*   // Input parameters of static method:
*
*   ODESLV : ODE solver
*     =1: explicit Euler
*     =3: Runge-Kutta-3
*   NITDT  : Maximum number of ODE pseudo time-steps
*   NGAPCS : Number of post-correction macro-steps
*
* end
*
* TDGridAdaptionStaticParams = record
*  
*   // Input parameters of static method: 
*
* end
*
* The Oxxxxx-constants below can be used to access the different
* variables inside of the structure-array.
* The variables in the GENERAL-blocks can be accessed directly, e.g.
*
*   ODESLV = IPARAM(OODESLV)
* 
* The caller has to fill all "input" variables in all these blocks.
************************************************************************

*****************************
* TIGridAdaptionStaticParams
*****************************

      INTEGER OODESLV,ONITDT,ONGAPCS
      PARAMETER (OODESLV = 1)
      PARAMETER (ONITDT  = 2)
      PARAMETER (ONGAPCS = 3)
      
C 2..24: reserved
C 24..31: user-defineable
      
      INTEGER SZGASI 
      PARAMETER (SZGASI=32)

*****************************
* TDGridAdaptionStaticParams
*****************************

C 1..24: reserved
C 24..31: user-defineable

      INTEGER SZGASD
      PARAMETER (SZGASD=32)
