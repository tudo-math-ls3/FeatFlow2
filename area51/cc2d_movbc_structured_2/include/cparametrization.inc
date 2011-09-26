!************************************************************************
!* This file contains a COMMON block with the parameters
!* about the discretization. The data here is filled using the DAT file.
!************************************************************************

      INTEGER IMESH,IRMESH
      
      CHARACTER*60 CPARM,CMESH

      COMMON /IPARTR/ IMESH,IRMESH,CPARM,CMESH
     
      SAVE /IPARTR/