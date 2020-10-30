C$Header: /mops2.0/src/CALC_CP.F 2020/10/30 Tatsuro Tanioka $
C$Name: mops-2_0 $
C Script to Calculate C:P of Phytoplankton and Zooplankton

C !INTERFACE: ==========================================================
C ----------------------------------------------------------------------
C CALCULATE C:P of Phytoplankton      
C CALCULATE C:P of Phytoplankton      
C CALCULATE C:P of Phytoplankton      
C CALCULATE C:P of Phytoplankton      
C CALCULATE C:P of Phytoplankton      
C by T. Tanioka (Nov 2020)
C Reference:
C --- For power-law model: Tanioka and Matsumoto (2020), Biogeosciences, 17(11), 2939-2954      
C --- For GM15 model: Galbraith and Martiny (2015), PNAS, 112(27), 8199-8204

C PREPROCESSOR OPTIONS:
 
C CALLED BY:	bgc_model

C CALLS:	

C INPUT/ARGUMENT LIST:	
C Option        Option for C:P model (1:Power-law, 2:GM15, Other:Redfield)
C PO4	        Phosphate Conc. [mmol P/m3]
C DIN	        Nitrate Conc. [mmol P/m3]
C Light         Irradiance [W/m2]      
C Temp          Temperature [deg C]

C INPUT/COMMON BLOCK: 
C Contants for power-law C:P calculation. See BGC_PARAMS.h and BGC_INI.F.
C Maximum and Minimum C:P. See BGC_INI.F

C OUTPUT/ARGUMENT LIST: 
C cp_phyto     C:P of Phytoplankton [molar]

      SUBROUTINE CALC_CP_PHYTO(Option,PO4,DIN,Light,Temp,cp_phyto)

      implicit none

#include "BGC_PARAMS.h"
      integer Option
      real*8 PO4,DIN,Light,Temp,cp_phyto
      real*8 rpc_det,rcp_det

! Calculate C:P          
      if(Option.eq.1) then         ! Option 1: Calculate C:P using Power-law
           rpc_det = par_bio_pc0/1000.0*(PO4/par_bio_po4_ref)**
     &     par_bio_spc_p*(DIN/par_bio_no3_ref)**par_bio_spc_n 
     &       *((Temp+273.15)/(par_bio_temp_ref+273.15))**
     &     par_bio_spc_t*(Light/par_bio_light_ref)**par_bio_spc_i
           rcp_det = 1.0/rpc_det
       elseif(Option.eq.2) then    ! Option 2: Calculate C:P with GM15
           rcp_det=1000.0d0/(6.9*PO4+4)
       else                        ! Other = Refield
           rcp_det = rcp
       endif
! Imposing Maximum and Minimum C:P
       if(rcp_det.gt.maxcp) then
           rcp_det = maxcp
       elseif(rcp_det.lt.mincp) then
           rcp_det = mincp
       endif
! Output         
       cp_phyto = rcp_det

       RETURN

       END
C ----------------------------------------------------------------------
C CALCULATE C:P of Zooplankton      
C CALCULATE C:P of Zooplankton      
C CALCULATE C:P of Zooplankton      
C CALCULATE C:P of Zooplankton      
C CALCULATE C:P of Zooplankton      
C by T. Tanioka (Nov 2020)
C Reference:
C ---  Based on homeostatic C:P parameter : Persson et al. (2010), Oikos, 119(5), 741-751      

C PREPROCESSOR OPTIONS:
C CALLED BY:	bgc_model

C CALLS:	

C INPUT/ARGUMENT LIST:	
C cp_phyto      Phytoplankton C:P

C INPUT/COMMON BLOCK: 
C Homestatic parameter par_zoo_cp_hom. See BGC_PARAMS.h and BGC_INI.F.

C OUTPUT/ARGUMENT LIST: 
C cp_zoo    C:P of Zooplankton [molar]
      SUBROUTINE CALC_CP_ZOO(cp_phyto,cp_zoo)

      implicit none

#include "BGC_PARAMS.h"
      real*8 rpc_zoo_ref,rpc_zoo
      real*8 cp_phyto,cp_zoo

! Calculate C:P of Zooplankton by Power-law equation
! (P:C)zoo = (P:C)ref^(1-H) * (P:C)prey^H, where H = Homeostatic parameter (0-1)
      rpc_zoo_ref = 1.0d0/rcp                           ! Reference P:C  = Redfield P:C
      rpc_zoo = (rpc_zoo_ref**(1.0d0-par_zoo_cp_hom))*
     &              (1.0d0/cp_phyto)**(par_zoo_cp_hom)
      if(rpc_zoo.gt.0.0d0) then
          cp_zoo = 1.0d0/rpc_zoo
      else
          cp_zoo = rcp
      endif

      RETURN

      END
