C$Header: /Users/ikriest/CVS/mops/BGC_DIAGNOSTICS.h,v 1.2 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

! KM restructured for backwards campatibility with PFT (8/2021)
      
! MOPS arrays for diagnostics
! f8_out added by T.Tanioka (Nov 2020)
      real*8 f1_out(bgc_ktotal),    ! primary production
     &       f2_out(bgc_ktotal),    ! zooplankton grazing
     &       f3_out(bgc_ktotal),    ! detritus or POP sedimentation
     &       f4_out(bgc_ktotal),    ! POP and DOP remineralization
     &       f5_out(bgc_ktotal),    ! river runoff
     &       f6_out(bgc_ktotal),    ! N fixation
     &       f7_out(bgc_ktotal),    ! denitrification
     &       f8_out(bgc_ktotal)     ! PAR

      
#ifndef PFT
#ifdef ORGCARBON
#ifndef FLEXCP
      ! MOPS + (CARBON assumed) + ORGCARBON
      ! fbgc9 = Sediment_C
      real*8 f9_out(bgc_ktotal)
      COMMON/DIAGVARS/f1_out,f2_out,f3_out,f4_out,f5_out,f6_out,f7_out,
     &                f8_out,f9_out
#else
      ! MOPS + (CARBON assumed) + ORGCARBON + FLEXCP
      ! TT: fbgc9 = Sediment_C, fbgc10 = Phytoplankton C:P uptake ratio, fbgc11 = Zooplankton C:P uptake ratio
      real*8 f9_out(bgc_ktotal)
      real*8 f10_out(bgc_ktotal),f11_out(bgc_ktotal)
      COMMON/DIAGVARS/f1_out,f2_out,f3_out,f4_out,f5_out,f6_out,f7_out,
     &                f8_out,f9_out,f10_out,f11_out
#endif

#else
      ! MOPS
      COMMON/DIAGVARS/f1_out,f2_out,f3_out,f4_out,f5_out,f6_out,f7_out,
     &                f8_out
#endif 
      
#else
! else of ifndef PFT

#ifndef ORGCARBON      
      ! MOPS + PFT
      ! fbgc9 = Photosynthesis for PFT 2
      real*8 f9_out(bgc_ktotal)
      COMMON/DIAGVARS/f1_out,f2_out,f3_out,f4_out,f5_out,f6_out,f7_out,
     &                f8_out,f9_out      
#else

#ifndef FLEXCP      
      ! MOPS + (CARBON assumed) + ORGCARBON + PFT
      ! fbgc9 = Photosynthesis for PFT 2, fbgc10 = Sediment_C
      real*8 f9_out(bgc_ktotal)
      real*8 f10_out(bgc_ktotal)
      COMMON/DIAGVARS/f1_out,f2_out,f3_out,f4_out,f5_out,f6_out,f7_out,
     &                f8_out,f9_out,f10_out
#else
      ! MOPS + (CARBON assumed) + ORGCARBON + PFT + FLEXCP
      ! fbgc9 = Photosynthesis for PFT 2, fbgc10 = Sediment_C
      ! fbgc11 = PFT1 C:P uptake ratio, fbgc12 = PFT2 C:P uptake ratio, ffbgc13 = Zooplankton C:P uptake ratio
      real*8 f9_out(bgc_ktotal)
      real*8 f10_out(bgc_ktotal)
      real*8 f11_out(bgc_ktotal),f12_out(bgc_ktotal),f13_out(bgc_ktotal)
      COMMON/DIAGVARS/f1_out,f2_out,f3_out,f4_out,f5_out,f6_out,f7_out,
     &                f8_out,f9_out,f10_out,f11_out,f12_out,f13_out
#endif

#endif      
! ends ORGCARBON
#endif      
! ends ifndef PFT
      

ckm      
ckm#ifndef ORGCARBON
ckm      COMMON/DIAGVARS/f1_out,f2_out,f3_out,f4_out,f5_out,f6_out,f7_out,
ckm     &                f8_out
ckm#else
ckm#ifndef PFT      
ckm      COMMON/DIAGVARS/f1_out,f2_out,f3_out,f4_out,f5_out,f6_out,f7_out,
ckm     &                f8_out,f9_out,f10_out,f11_out
ckm#else
ckm      COMMON/DIAGVARS/f1_out,f2_out,f3_out,f4_out,f5_out,f6_out,f7_out,
ckm     &     f8_out,f9_out,f10_out,f11_out,f12_out,f13_out
ckm#endif
ckm#endif
