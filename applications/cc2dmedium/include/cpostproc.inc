************************************************************************
* This file contains a COMMON block with the parameters
* for the standard postprocessing routines.
* The data here is filled using the DAT file.
************************************************************************

      INTEGER IERANA,IAUSAV,ISAV  ,INSAV ,INSAVN
      INTEGER IFUSAV,IFPSAV,IFXSAV,IGMV  ,IFINIT,IGINIT,IUPSAV
      INTEGER ITRCES
      
      DOUBLE PRECISION DTFILM,DTGMV 

      CHARACTER*60 CFLAUT,CFLDUP,CFLDU ,CFLDV ,CFLDP ,CFLISO
      CHARACTER*60 CFLGMV,CFLFX1,CFLFY1,CFLFX2,CFLFY2,CFLFX3
      CHARACTER*60 CFLFY3,CFLVLX,CFLVLY,CFLPRE,CFLSTR,CFLIPR
      CHARACTER*60 CFLSAV
      
      COMMON /IPOSTP/ DTFILM,DTGMV,
     *                IERANA,IAUSAV,ISAV,INSAV ,INSAVN,IUPSAV,
     *                IFUSAV,IFPSAV,IFXSAV,IGMV  ,IFINIT,IGINIT,
     *                ITRCES,
     *                CFLAUT,CFLDUP,CFLDU ,CFLDV ,CFLDP ,CFLISO,
     *                CFLGMV,CFLFX1,CFLFY1,CFLFX2,CFLFY2,CFLFX3,
     *                CFLFY3,CFLVLX,CFLVLY,CFLPRE,CFLSTR,CFLIPR,
     *                CFLSAV
     
      SAVE /IPOSTP/ 
      

      