************************************************************************
* This file contains a COMMON block with the parameters
* about the discretization. The data here is filled using the DAT file.
************************************************************************

      INTEGER ISTART,ISOL,NLMIN,NLMAX,ISTOK,IRHS,IBDR,IELT,IUPW
      INTEGER IMASS,IMASSL,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN ,ICUBB,ICUBF
      INTEGER IAPRM,IINFCF
      
      DOUBLE PRECISION RE    ,UPSAM,DMTEP ,DPUMAX
      
      CHARACTER*60 CSTART,CSOL

      COMMON /IDISCR/ RE    ,UPSAM,DMTEP ,DPUMAX,
     *                ISTART,ISOL,NLMIN,NLMAX,ISTOK,IRHS,IBDR,IELT,IUPW,
     *                IMASS,IMASSL,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN ,
     *                ICUBB,ICUBF,IAPRM ,IINFCF,
     *                CSTART,CSOL
     
      SAVE /IDISCR/