C$Header: /Users/ikriest/CVS/mops/mops_biogeochem_set_params.F,v 1.3 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

      SUBROUTINE MOPS_BIOGEOCHEM_SET_PARAMS(npar,parvec)

      implicit none

#include "BGC_PARAMS.h"

      integer npar
      real*8 parvec(npar)
      
C BGC PARAMETERS FOR SURFACE OCEAN BIOLOGY      
C Edited by T.Tanioka
C Calibrated BGC Parameters (different for MIT2.8 and ECCO)
        detmartin  = parvec(1) !w = a*z+b
        ro2ut      = parvec(2) !redfield -02:P ratio
        nfix       = parvec(3) !Max rate of N fixation [d-1]
        subdin     = parvec(4) !no oxic degradation below this level
        ACkbaco2   = parvec(5) !Half sat.-constant for oxic degradation.
        ACkbacdin  = parvec(6) !Half sat.-constant for suboxic degradation
C Flexible C:P Stoichiometry Parameters
#ifdef FLEXCP
         cp_option = parvec(7)            ! Which C:P model? (1=:Power-law, 2:Galbraith and Martiny, 0:Redfield)
         maxcp = parvec(8)                ! Upper bound C:P
         mincp = parvec(9)                ! Lower bound C:P
c Paramters related to C:P Power-law
         par_bio_pc0 = parvec(10)         ! Referennce P:C (permil) = Redfield P:C 
         par_bio_po4_ref = parvec(11)     ! Reference PO4 (mmol m-3)
         par_bio_no3_ref = parvec(12)     ! Reference NO3 (mmol m-3)
         par_bio_temp_ref = parvec(13)    ! Reference Temperature (deg C)
         par_bio_light_ref = parvec(14)   ! Reference irradiance (W m-2)
         par_bio_spc_p = parvec(15)       ! s_P:C^PO4
         par_bio_spc_n = parvec(16)       ! s_P:C^NO3
         par_bio_spc_i = parvec(17)       ! s_P:C^I
         par_bio_spc_t = parvec(18)       ! s_P:C^T
c Parameter related to C:P of Zooplankton 
         par_zoo_cp_hom = parvec(19)      ! Homeostatic Parameter (0 = Completely homeostatic, 1: C:P Zoo = C:P of prey)       
#endif
        
      RETURN
      END
