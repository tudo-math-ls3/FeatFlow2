************************************************************************
* This file contains a COMMON block with the parameters
* for the linear sub-solver of the nonlinear solver, CC2D-type.
* The data here is filled using the DAT file.
************************************************************************

      INTEGER ILINSL,ICYCLE,ILMIN ,ILMAX ,IINT  ,IAVPR
      INTEGER ISM   ,NSM   ,NSMFAC,ISL   ,NSL   ,IAPR, IPCMG
      
      DOUBLE PRECISION DMPMG ,EPSMG ,DPREP ,AMINMG,AMAXMG
      DOUBLE PRECISION RLXSM ,RLXSL ,DMPSL ,EPSSL ,OMGLIN
      
      COMMON /LINSCC/ DMPMG ,EPSMG ,DPREP ,AMINMG,AMAXMG,
     *                RLXSM ,RLXSL ,DMPSL ,EPSSL, OMGLIN,
     *                ILINSL,ICYCLE,ILMIN ,ILMAX ,IINT  ,IAVPR,
     *                ISM   ,NSM   ,NSMFAC,ISL   ,NSL   ,IAPR, IPCMG

      SAVE /LINSCC/
