************************************************************************
* This file contains a COMMON block with the parameters
* for the linear sub-solver of the nonlinear solver, CC2D-type.
* The data here is filled using the DAT file.
************************************************************************

      INTEGER INLMIN,INLMAX
      
      DOUBLE PRECISION EPSD  ,EPSDIV,EPSUR ,EPSPR 
      DOUBLE PRECISION DMPD  ,OMGMIN,OMGMAX,OMGINI
      
      COMMON /NLISCC/ DMPD  ,OMGMIN,OMGMAX,OMGINI,
     *                EPSD  ,EPSDIV,EPSUR ,EPSPR, 
     *                INLMIN,INLMAX

      SAVE /NLISCC/ 
