************************************************************************
* Warning!
*
* In comparison to previous FEATFLOW versions, this file has slightly
* changed. The routines in this file are now renamed to
*
*   FPARX, FPARY, FTMAX, FTNVBC
*
* instead of PARX,PARY,TMAX,TNBC. This allowes to dynamically switch
* between FEAT and OMEGA parametrization by setting the parameter
* IMESH=0/1 without exchanging the file parq2d.f.
*
* In contrast to the original parq2d.f from earlier FEATFLOW versions,
* an additional routine FTNBC was added that returns the number of
* boundary components NBCT. This makes the parametrization independent
* of any information in the triangulation structures, as NBC was
* preciously stored in the /TRADx/ COMMON blocks!
************************************************************************

************************************************************************
* This file realizes the hardcoded parametrization corresponding to the
* triangulation MESH3. This is just a rectangular mesh in the size of
* bench1 which was used for tests with a horizontally moving ball in
* that domain.
*
* By setting the parametrization to "Feat-Parametrization" in the
* DAT file of the application and by setting the triangulation to
* "mesh3.tri" in the same DAT file, the geometry can be activated.
* Otherwise, the OMEGA-parametrization is used, which uses the .PRM
* and .TRI file specified in the DAT file.
*
* The user can customize this file if it's necessary to analytically
* specify the parametrization.
************************************************************************

************************************************************************
* Get X-coordinate of parameter value.
*
* This routine returns the X-coordinate of a parameter value.
*
* In:
*   T      - parameter value
*   IBCT   - number of boundary component
*
* Return:
*   X-coordinate of parameter value T on boundary component IBCT.
************************************************************************

      DOUBLE PRECISION FUNCTION FPARX(T,IBCT)

      IMPLICIT NONE
       
C     parameters
       
      DOUBLE PRECISION T
      INTEGER IBCT
      
C     Constants and local variables

      INTEGER MAXI,MAXJ,K
      PARAMETER(MAXI=25,MAXJ=73,K=2*MAXI+2*(MAXJ-2))
      
      DOUBLE PRECISION DLL
      PARAMETER(DLL=2.0D0)
       
      DOUBLE PRECISION X(1:K),DX(1:MAXI-1)
      DOUBLE PRECISION TPAR(1:K+1) 

      INTEGER I,J
      DOUBLE PRECISION DT

      X(1)=0.0D0
      TPAR(1)=0.0D0

C      DX=DLL/DBLE(MAXI-1)
      DO I=1,MAXI-1
        DX(I)=1.0/12.0
      ENDDO  
                    
      DT=1.0D0

      J=0
      DO I=2,MAXI-1
        J=J+1
        X(I)=X(I-1)+DX(J)        
      ENDDO
      DO I=MAXI,MAXI+MAXJ-1
        X(I)=DLL        
      ENDDO
      J=0
      DO I=MAXI+MAXJ,2*MAXI+MAXJ-3
        J=J+1
        X(I)=X(I-1)-DX(J)        
      ENDDO       
      DO I=2*MAXI+MAXJ-2,2*MAXI+2*(MAXJ-2)
        X(I)=X(1)        
      ENDDO       

      DO I=2,K+1
        TPAR(I)=TPAR(I-1)+DT        
      ENDDO       

C----------------------------------------------------

       IF(T.GE.TPAR(K))  THEN
         FPARX=X(K)+(X(1)-X(K))/(TPAR(K+1)-TPAR(K))*(T-TPAR(K)) 
         RETURN
       ENDIF   
       
      DO I=K-1,1,-1
        IF(T.GE.TPAR(I))  THEN
          FPARX=X(I)+(X(I+1)-X(I))/(TPAR(I+1)-TPAR(I))*(T-TPAR(I)) 
          RETURN
        ENDIF
      ENDDO
           
99999 END

************************************************************************
* Get Y-coordinate of parameter value.
*
* This routine returns the Y-coordinate of a parameter value.
*
* In:
*   T      - parameter value
*   IBCT   - number of boundary component
*
* Return:
*   Y-coordinate of parameter value T on boundary component IBCT.
************************************************************************

      DOUBLE PRECISION FUNCTION FPARY(T,IBCT)
       
      IMPLICIT NONE
       
C     parameters
       
      DOUBLE PRECISION T
      INTEGER IBCT
      
C     constants and local variables
      
      INTEGER MAXI,MAXJ,K
      PARAMETER(MAXI=25,MAXJ=73,K=2*MAXI+2*(MAXJ-2)) 
      
      DOUBLE PRECISION DHH
      PARAMETER(DHH=6.0D0) 
      
      DOUBLE PRECISION Y(1:K),DY(1:MAXJ-1)
      DOUBLE PRECISION TPAR(1:K+1) 
      
      INTEGER I,J
      DOUBLE PRECISION DT

      Y(1)=0.0D0
      TPAR(1)=0.0D0

C      DY=DHH/DBLE(MAXJ-1)  
      DO I=1,MAXJ-1
        DY(I)=1.0/12.0
      ENDDO 
                    
      DT=1.0D0

      DO I=2,MAXI
        Y(I)=Y(I-1)        
      ENDDO
      J=0
      DO I=MAXI+1,MAXI+MAXJ-2
        J=J+1
        Y(I)=Y(I-1)+DY(J)        
      ENDDO
      DO I=MAXI+MAXJ-1,2*MAXI+MAXJ-2
        Y(I)=DHH        
      ENDDO
      J=0       
      DO I=2*MAXI+MAXJ-1,2*MAXI+2*(MAXJ-2)
        J=J+1
        Y(I)=Y(I-1)-DY(J)        
      ENDDO       

      DO I=2,K+1
        TPAR(I)=TPAR(I-1)+DT        
      ENDDO       

C----------------------------------------------------

      IF (T.GE.TPAR(K)) THEN
        FPARY=Y(K)+(Y(1)-Y(K))/(TPAR(K+1)-TPAR(K))*(T-TPAR(K)) 
        RETURN
      ENDIF   
       
      DO I=K-1,1,-1
        IF (T.GE.TPAR(I)) THEN
          FPARY=Y(I)+(Y(I+1)-Y(I))/(TPAR(I+1)-TPAR(I))*(T-TPAR(I)) 
          RETURN
        ENDIF
      ENDDO
           
99999 END

************************************************************************
* Get maximum parameter value
*
* This routine returns the maximum parameter value on a boundary
* component
*
* In:
*   IBCT   - number of boundary component
*
* Return:
*   Maximum parameter value on that boundary component.
************************************************************************

      DOUBLE PRECISION FUNCTION FTMAX(IBCT)
      
      IMPLICIT NONE
      
      INTEGER IBCT
      
      INTEGER MAXI,MAXJ
      PARAMETER(MAXI=25,MAXJ=73)
      
      FTMAX=DBLE(2*MAXI+2*(MAXJ-2))
       
      END
 
************************************************************************
* Get number of real boundary components
*
* This routine returns the number of real boundary components in the
* parametrization.
*
* In:
*   -
*
* Return:
*   Number of boundary components
************************************************************************

      INTEGER FUNCTION FTNBC ()

      IMPLICIT NONE
      
      FTNBC = 1

      END
