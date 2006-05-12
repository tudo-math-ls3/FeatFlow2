#!/usr/bin/env gmake

FEATFLOW=../..

include $(FEATFLOW)/Globals.mk

LIBNAME=blas

SRCLIST= \
caxpy.f   csscal.f  dger.f    dtpsv.f   sger.f    strmm.f   zherk.f \
ccopy.f   cswap.f   dnrm2.f   dtrmm.f   snrm2.f   strmv.f   zhpmv.f \
cdotc.f   csymm.f   drot.f    dtrmv.f   srot.f    strsm.f   zhpr2.f \
cdotu.f   csyr2k.f  drotg.f   dtrsm.f   srotg.f   strsv.f   zhpr.f \
cgbmv.f   csyrk.f   drotm.f   dtrsv.f   srotm.f   xerbla.f  \
cgemm.f   ctbmv.f   drotmg.f  dzasum.f  srotmg.f  zaxpy.f   zscal.f \
cgemv.f   ctbsv.f   dsbmv.f   dznrm2.f  ssbmv.f   zcopy.f   zswap.f \
cgerc.f   ctpmv.f   dscal.f   icamax.f  sscal.f   zdotc.f   zsymm.f \
cgeru.f   ctpsv.f   dsdot.f   idamax.f  sspmv.f   zdotu.f   zsyr2k.f \
chbmv.f   ctrmm.f   dspmv.f   isamax.f  sspr2.f   zdrot.f   zsyrk.f \
chemm.f   ctrmv.f   dspr2.f   izamax.f  sspr.f    zdscal.f  ztbmv.f \
chemv.f   ctrsm.f   dspr.f    lsame.f   sswap.f   zgbmv.f   ztbsv.f \
cher2.f   ctrsv.f   dswap.f   sasum.f   ssymm.f   zgemm.f   ztpmv.f \
cher2k.f  dsymm.f   saxpy.f   ssymv.f   zgemv.f   ztpsv.f \
cher.f    dasum.f   dsymv.f   scasum.f  ssyr2.f   zgerc.f   ztrmm.f \
cherk.f   daxpy.f   dsyr2.f   scnrm2.f  ssyr2k.f  zgeru.f   ztrmv.f \
chpmv.f   dcabs1.f  dsyr2k.f  scopy.f   ssyr.f    zhbmv.f   ztrsm.f \
chpr2.f   dcopy.f   dsyr.f    sdot.f    ssyrk.f   zhemm.f   ztrsv.f \
chpr.f    ddot.f    dsyrk.f   sdsdot.f  stbmv.f   zhemv.f \
crotg.f   dgbmv.f   dtbmv.f   sgbmv.f   stbsv.f   zher2.f \
cscal.f   dgemm.f   dtbsv.f   sgemm.f   stpmv.f   zher2k.f \
csrot.f   dgemv.f   dtpmv.f   sgemv.f   stpsv.f   zher.f

# zrotg.f  removed from the file list. 
# Not compatible with 64 bit compilation (64 Bit integer, 128 bit double)

include $(FEATFLOW)/Rules_libs.mk

update:
	@rm -f blas.tgz
	@wget http://netlib.org/blas/blas.tgz
	@(cd src ; tar -xvzf ../blas.tgz )
