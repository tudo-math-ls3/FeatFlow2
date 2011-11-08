************************************************************************
* This file contains a COMMON block with the parameters
* for the linear sub-solver of the nonlinear solver, CC2D-type.
* The data here is filled using the DAT file.
************************************************************************

      INTEGER ISTAT ,NITNS ,IFRSTP,IEPSAD,IADTIM,IREPIT,IADIN ,IDMESH
      INTEGER IMSCRT,IMSREP,IMSSTP
      
      DOUBLE PRECISION TIMEMX,TIMENS,TSTEP ,THETA ,EPSNS ,DTMIN ,DTMAX
      DOUBLE PRECISION DTFACT,EPSADU,TIMEIN,EPSADI,EPSADL,DMSMTM
      
      COMMON /TIDISC/ TIMEMX,TIMENS,TSTEP ,THETA ,EPSNS ,DTMIN ,DTMAX,
     *                DTFACT,EPSADU,TIMEIN,EPSADI,EPSADL,DMSMTM,
     *                ISTAT ,NITNS ,IFRSTP,IEPSAD,IADTIM,IREPIT,IADIN ,
     *                IDMESH,IMSCRT ,IMSREP,IMSSTP

      SAVE /TIDISC/
      



