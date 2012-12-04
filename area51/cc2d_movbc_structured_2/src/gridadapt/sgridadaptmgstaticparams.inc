************************************************************************
* This files describes the parameter structure of the multiple-grid
* static grid adaption routine XSMMGS.
*
* The routine takes a parameter block IPARAM which describes
* the behaviour of the grid adaption over multiple levels. 
* Each of the configurations in this parameter block can be
* transferred to/activated in the static grid adaption block using the
* routine INGASM.
* The structure is realized as array,
*
*   IPARAM = array [1..SZGAMI] of integer 
*          = TIGridAdaptionMGStaticParams
*
* containing the following variables:
*
* TIGridAdaptionMGStaticParams = record
*  
*   // Input parameters of static method:
*
*  IMLGAS : Number of configurations
*  IGALVL : For each of the IMLGAS configurations:
*           number of the level where to perform
*           grid adaption.
*            >0: Grid adaption on level IGALVL(.).
*           <=0: Grid adaption on level NLMAX+IGALVL(.)
*  IGASLV : For each of the IMLGAS configurations:
*           number of the level where to solve the linear system
*            >0: Solve on level IGALVL(.).
*           <=0: Solve on level NLMAX+IGALVL(.)
*  MNGASTP : For each of the IMLGAS configurations: 
*            Number of grid adaption (macro-) steps
*  MNGAPCS : For each of the IMLGAS configurations: 
*            Number of grid adaption post-correction steps
*  MNPRGSM : For each of the IMLGAS configurations: 
*            Number of Laplacian presmoothing steps
*  MNPOGSM : For each of the IMLGAS configurations: 
*            Number of Laplacian postsmoothing steps
*  MNFSMTH : For each of the IMLGAS configurations: 
*            Number of monitor function smoothing steps
*
* end
*
* The Oxxxxx-constants below can be used to access the different
* variables inside of the structure-array.
* The caller has to fill all "input" variables in all these blocks.
************************************************************************

*******************************
* TIGridAdaptionMGStaticParams
*******************************

C     Maximum count of levels in multiple-grid grid-adaption

      INTEGER NGAMXL
      PARAMETER (NGAMXL=16)

      INTEGER OIMLGAS ,OIGALVL ,OMNGASTP,OMNGAPCS,OMNPRGSM
      INTEGER OMNPOGSM,OMNFSMTH,OIGASLV

      PARAMETER (OIMLGAS  = 1)
      PARAMETER (OIGALVL  = 2)
      PARAMETER (OIGASLV  = 2+1*NGAMXL)
      PARAMETER (OMNGASTP = 2+2*NGAMXL)
      PARAMETER (OMNGAPCS = 2+3*NGAMXL)
      PARAMETER (OMNPRGSM = 2+4*NGAMXL)
      PARAMETER (OMNPOGSM = 2+5*NGAMXL)
      PARAMETER (OMNFSMTH = 2+6*NGAMXL)
      
      INTEGER SZGAMI 
      PARAMETER (SZGAMI=OMNFSMTH+NGAMXL)
