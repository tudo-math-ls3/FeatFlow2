************************************************************************
* This file describes structure information about system matrices
* of the 2D-problem for solving Navier Stokes equations.
************************************************************************

************************************************************************
* With each level of the discretization, there are multple matrices
* associated, most of them appear in the memory normally, some of them
* are frequently updated or rebuild.
*
* - Mass matrix
* - Stokes ("Laplace")-matrix
* - Nonlinear system matrix
*
* The structure itselves is realized as integer array:
*   IPARAM = array [1..SZN2MI] of integer
*          = TNS2DMatrixParams
*
* TNS2DMatrixParams = record
*
*   // The first part of the structure defines generally the dimension
*   // of the matrices:
*
*   NEQA   : number of rows/columns in the velocity submatrix
*   NEQB   : number of columns in the pressure submatrix
*
*   // The second part of the structure defines the mass matrix.
*
*   NM     : If the mass matrix is not lumped:
*              = Number of entries in the mass matrix.
*            If the mass matrix is lumped:
*              = Number of entries in the original (not lumped) mass
*                matrix. The number of entries the the mass matrix
*                vector identified by LM is NEQA!
*   LM     : If the mass matrix is not lumped:
*              = Handle to array [1..NM] of double
*                The entries in the mass matrix.
*            If the mass matrix is lumped:
*              = Handle to array [1..NEQA] of double
*                The entries on the diagonal of the mass matrix.
*            =0, if no mass matrix exists.
*   IMALMP : =0: lumped mass matrix. LCOL=LLD=0.
*                standard.
*                LM is a handle to an array [1..NEQA] of double, which
*                represents the diagonal entries of the lumped mass
*                matrix. LCOLM/LLDM point to the structure of the
*                original (not lumped) mass matrix (which usially
*                coincides with LCOLA1/LLDA1).
*            =1: real mass matrix
*                LCOLM/LLDM might coincide with LCOLA1/LLDA1 of the
*                system matrix, because the structure is identical.
*   LCOLM  : Handle to array [1..NM] of integer
*            Column structure of the real (not lumped) mass matrix
*            (also if the mass matrix is lumped!).
*            Might coincode with LCLA1.
*            =0, if no mass matrix exists.
*   LLDM   : Handle to array [1..NEQ+1] of integer
*            Row pointers of the real (not lumped) mass matrix
*            (also if the mass matrix is lumped!).
*            Might coincode with LLDA1.
*            =0, if no mass matrix exists.
*
*   // The third part of the structure defines the Stokes matrix,
*   // which is the diffusive (Laplacian) part of the nonlinear
*   // system matrix.
*   // LCOLST/LLDST might coincide with LCOLA/LLDA of the System
*   // matrix because the structure is identical.
*
*   NST    : Number of entries in the Stokes matrix (= NA normally)
*   LST    : Handle to array [1..NST] of double
*            The entries in the Stokes matrix.
*            =0, if no Stokes matrix exists.
*   LCOLST : Handle to array [1..NST] of integer
*            Column number identifier;
*            might coincode with LCLA1
*   LLDST  : Handle to array [1..NEQ+1] of integer
*            Row pointers;
*            might coincide with LLDA1
*
*   // Now, the definition of the system matrix follows.
*   // The system matrix has the followinf form:
*   //
*   //    ( A1    C2    B1)
*   //    ( C1    A2    B2)
*   //    ( B1^T  B2^T  0 )
*
*   // The 4th part of this structure defines the velocity part of
*   // the nonlinear system matrix. This is frequently rebuild and
*   // updated during the solution process.
*   // The velocity matrix consists of two matrices, one for each
*   // dimension.
*
*   NA1    : Number of entries in the 1st system matrix
*   LA1    : Handle to array [1..NA1] of double
*            The entries in the 1st system matrix.
*            =0, if no 1st system exists.
*   LCLA1  : Handle to array [1..NA1] of integer
*            Column number identifier
*   LLDA1  : Handle to array [1..NEQ+1] of integer
*            Row pointers
*
*   // Remark that the 2nd system matrix might coincide with the
*   // first one!
*
*   NA2    : Number of entries in the 2nd system matrix (= NA1 normally)
*   LA2    : Handle to array [1..NA2] of double
*            The entries in the 2nd system matrix.
*            =0, if no 2nd system exists.
*            Might coincide with LA1.
*   LCLA2  : Handle to array [1..NA2] of integer
*            Column number identifier
*            Might coincide with LCLA1.
*   LLDA2  : Handle to array [1..NEQ+1] of integer
*            Row pointers
*            Might coincide with LLSA1.
*
*   // The 5th part of the structure defines the coupling matrices
*   // C1 and C2 in the velocity part of the nonlinear system matrix.
*   // Normally this matrix does not exist (LCx=0), but it might
*   // be used in the future.
*
*   NC1    : Number of entries in the 1st coupling matrix
*   LC1    : Handle to array [1..NC1] of double
*            The entries in the 1st system matrix.
*            =0, if no 1st coupling exists.
*   LCLC1  : Handle to array [1..NC1] of integer
*            Column number identifier
*   LLDC1  : Handle to array [1..NEQ+1] of integer
*            Row pointers
*
*   // Remark that the 2nd coupling matrix might coincide with the
*   // first one!
*
*   NC2    : Number of entries in the 2nd coupling matrix (= NA1 normally)
*   LC2    : Handle to array [1..NA2] of double
*            The entries in the 2nd coupling matrix.
*            =0, if no 2nd coupling exists.
*            Might coincide with LC1.
*   LCLC2  : Handle to array [1..NA2] of integer
*            Column number identifier
*            Might coincide with LCLC1.
*   LLDC2  : Handle to array [1..NEQ+1] of integer
*            Row pointers
*            Might coincide with LLSC1.
*
*   // The 6th part of the structure defines the B1/B2 pressure
*   // submatrices. As this matrix is rectangular rather that quadratic,
*   // it's stored in matrix structure 9.
*
*   NB1    : Number of entries in the 1st pressure matrix
*   LB1    : Handle to array [1..NB1] of double
*            The entries in the 1st pressure matrix.
*            =0, if no 1st pressure exists.
*   LCLB1  : Handle to array [1..NA] of integer
*            Column number identifier
*   LLDB1  : Handle to array [1..NEQ+1] of integer
*            Row pointers
*
*   // Remark that the structure of the 2nd pressure matrix might
*   // coincide with the first one!
*
*   NB2    : Number of entries in the 2nd pressure matrix
*            (= NB1 normally)
*   LB2    : Handle to array [1..NA] of double
*            The entries in the 2nd pressure matrix.
*            =0, if no 2nd pressure exists.
*            Might coincide with LB1, but usually does not.
*   LCLB2  : Handle to array [1..NA] of integer
*            Column number identifier
*            Might coincide with LCLB1
*   LLDB2  : Handle to array [1..NEQ+1] of integer
*            Row pointers
*            Might coincide with LLDB1
* end
************************************************************************

      INTEGER ONEQA  ,ONEQB
      INTEGER ONM    ,OLM    ,OIMALMP,OLCOLM ,OLLDM
      INTEGER ONST   ,OLST   ,OLCOLS ,OLLDST
      INTEGER ONA1   ,OLA1   ,OLCLA1 ,OLLDA1
      INTEGER ONA2   ,OLA2   ,OLCLA2 ,OLLDA2
      INTEGER ONB1   ,OLB1   ,OLCLB1 ,OLLDB1
      INTEGER ONB2   ,OLB2   ,OLCLB2 ,OLLDB2
      INTEGER ONC1   ,OLC1   ,OLCLC1 ,OLLDC1
      INTEGER ONC2   ,OLC2   ,OLCLC2 ,OLLDC2
      INTEGER SZN2MI

      PARAMETER (ONEQA   =  1)
      PARAMETER (ONEQB   =  2)
      
      PARAMETER (ONM     =  5)
      PARAMETER (OLM     =  6)
      PARAMETER (OIMALMP =  7)
      PARAMETER (OLCOLM  =  8)
      PARAMETER (OLLDM   =  9)
      
      PARAMETER (ONST    = 10)
      PARAMETER (OLST    = 11)
      PARAMETER (OLCOLS  = 12)
      PARAMETER (OLLDST  = 13)
      
      PARAMETER (ONA1    = 15)
      PARAMETER (OLA1    = 16)
      PARAMETER (OLCLA1  = 17)
      PARAMETER (OLLDA1  = 18)
      
      PARAMETER (ONA2    = 20)
      PARAMETER (OLA2    = 21)
      PARAMETER (OLCLA2  = 22)
      PARAMETER (OLLDA2  = 23)
      
      PARAMETER (ONB1    = 25)
      PARAMETER (OLB1    = 26)
      PARAMETER (OLCLB1  = 27)
      PARAMETER (OLLDB1  = 28)
      
      PARAMETER (ONB2    = 30)
      PARAMETER (OLB2    = 31)
      PARAMETER (OLCLB2  = 32)
      PARAMETER (OLLDB2  = 33)
      
      PARAMETER (ONC1    = 35)
      PARAMETER (OLC1    = 36)
      PARAMETER (OLCLC1  = 37)
      PARAMETER (OLLDC1  = 38)
      
      PARAMETER (ONC2    = 40)
      PARAMETER (OLC2    = 41)
      PARAMETER (OLCLC2  = 42)
      PARAMETER (OLLDC2  = 43)

      PARAMETER (SZN2MI  = 64)
      