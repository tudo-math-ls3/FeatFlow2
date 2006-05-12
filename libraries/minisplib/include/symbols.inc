
c     --------------------
c     Solver method codes:
c     --------------------
      integer sBICG, sCGNE, sCGNR, sCGS, sCGSTAB, sGMRES, sQMR, sCSCGS
      integer sJACOBI, sGAUSS, sSUCR, sORTHM, sNOSLV, sBICGS, sGMRESt
      parameter(sBICG=1, sCGNE=2, sCGNR=3, sCGS=4, sCGSTAB=5, sGMRES=6)
      parameter(sQMR=7, sBICGS=8, sGMRESt=9, sNOSLV=0)
      parameter(sJACOBI=10, sGAUSS=11, sSUCR=12, sORTHM = 13,sCSCGS=14)

c     ----------------------------
c     Preconditioner method codes:
c     ----------------------------
      integer pNOPC, pILUS, pMILU, pILUT, pSSOR, pTRID, pILU0, pEIMGS
      parameter(pNOPC=0, pILUS=1, pMILU=2)
      parameter(pILUT=3, pSSOR=4, pTRID=5, pILU0=6)
      parameter(pEIMGS=7)

c     ---------------------
c     Scaling method codes:
c     ---------------------
      integer NOSC, RWSC, COLSC, RCSC, MXSC
      parameter(NOSC=0, RWSC=2, COLSC=3, RCSC=4, MXSC=1)
