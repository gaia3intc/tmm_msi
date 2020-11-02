C$Header: /Users/ikriest/CVS/mops/BGC_DIAGNOSTICS.h,v 1.2 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

! arrays for diagnostics
! f8_out added by T.Tanioka (Nov 2020)
      real*8 f1_out(bgc_ktotal),
     &       f2_out(bgc_ktotal),
     &       f3_out(bgc_ktotal),
     &       f4_out(bgc_ktotal),
     &       f5_out(bgc_ktotal),
     &       f6_out(bgc_ktotal),
     &       f7_out(bgc_ktotal),
     &       f8_out(bgc_ktotal)

#ifdef ORGCARBON
! Added By T.Tanioka (Nov 2020)
! fbgc9 = Sediment_C, fbgc10 = Phytoplankton C:P uptake ratio, fbgc11 = Zooplankton C:P uptake ratio
      real*8 f9_out(bgc_ktotal),f10_out(bgc_ktotal),f11_out(bgc_ktotal)
#endif

#ifndef ORGCARBON
      COMMON/DIAGVARS/f1_out,f2_out,f3_out,f4_out,f5_out,f6_out,f7_out,
     &                f8_out
#else
      COMMON/DIAGVARS/f1_out,f2_out,f3_out,f4_out,f5_out,f6_out,f7_out,
     &                f8_out,f9_out,f10_out,f11_out
#endif
