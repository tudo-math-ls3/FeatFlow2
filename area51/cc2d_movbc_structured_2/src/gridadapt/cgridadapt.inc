************************************************************************
* This file contains COMMON blocks holding the parameters of the grid
* adaption routines.
*
* The include file "sgridadaptmgstaticparams.inc" has to be included
* prior to this one!
************************************************************************

C     parameters in the DAT file:

      INTEGER NITDT,ICUBGS,ICURGS
      INTEGER IGSSLV,IRELGS,NITSGS,IPCSMG,NSMSGS,IGAMET,IGSOUT
      INTEGER ISRTGS,ICGSGS,IRCGGS,IPCGGS,NICGGS,IPCSGS
      INTEGER IGAMAL,NLASMF,NGASTF,IARCTR,ICGSLV,IMNREL
      INTEGER ODESLV,ISMMTH,ICOLMC,IRCGRD,ICUNGS

      DOUBLE PRECISION EPS0GS,EPSGSS,OMGPGS,EPSCGS,OMCGGS,DMNTOL,DAPOSM
      DOUBLE PRECISION DFAPSM,DFPSRF,DMNSFC
      
      COMMON /GSMINT/ IGAMET,NITDT,ICUBGS,
     *                ICURGS,IGSOUT,IGSSLV,IRELGS,NITSGS,IPCSMG,NSMSGS,
     *                ISRTGS,ICGSGS,IRCGGS,IPCGGS,NICGGS,IPCSGS,
     *                IGAMAL,NLASMF,NGASTF,IARCTR,ICGSLV,IMNREL,
     *                ODESLV,ISMMTH,ICOLMC,IRCGRD,ICUNGS
     
      COMMON /GSMDBL/ EPS0GS,EPSGSS,OMGPGS,EPSCGS,OMCGGS,DMNTOL,DAPOSM,
     *                DFAPSM,DFPSRF,DMNSFC
     
      INTEGER IMLGAS,IGALVL,IGASLV,NGASTP,NGAPCS,NPRGSM,NPOGSM,NFSMTH
     
      COMMON /GSMMLV/ IMLGAS,
     *                IGALVL(NGAMXL),IGASLV(NGAMXL),NGASTP(NGAMXL),
     *                NGAPCS(NGAMXL),NPRGSM(NGAMXL),NPOGSM(NGAMXL),
     *                NFSMTH(NGAMXL)
      
      SAVE /GSMINT/,/GSMDBL/
      
