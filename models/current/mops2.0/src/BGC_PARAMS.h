C$Header: /Users/ikriest/CVS/mops/BGC_PARAMS.h,v 1.2 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

C EVERYTHING RELATED TO THE MAIN BGC TRACER ARRAY, ITS INDICES AND TYPES

c a dummy for the nominal number of vertical layers
      INTEGER bgc_ktotal
      PARAMETER(bgc_ktotal=100)

c the total number of bgc tracers
      INTEGER bgc_ntracer
      
c number of PFTs (maxpft sets array size; npft read in biogem_parameter.txt)
      INTEGER npft, maxpft
      PARAMETER(maxpft=2)

c the indices of tracers
      INTEGER ipo4,idop,iphy,izoo,ioxy,idet,idin
      PARAMETER(ipo4=1,    !PO4
     &          idop=2,    !DOP
     &          ioxy=3,    !Oxygen
     &          iphy=4,    !Phyto-P
     &          izoo=5,    !Zoo-P
     &          idet=6,    !Detritus-P
     &          idin=7)    !DIN
      
#ifndef CARBON

#ifndef PFT
      ! MOPS (P-based, 7 tracers)
      PARAMETER(bgc_ntracer=7)
#else
      ! MOPS + PFT (still P-based)
      INTEGER iphy2=8
      PARAMETER(bgc_ntracer=8)
#endif
! ends PFT
      
#else
! else CARBON
      
      ! MOPS + CARBON
      INTEGER idic,ialk,iphyc,izooc
      PARAMETER(idic=8,ialk=9,iphyc=10,izooc=11)
      ! connect carbon exchange and P-based BGC
      real*8 ocmip_alkfac,ocmip_silfac
      COMMON/CO2SURFACE/ocmip_alkfac,ocmip_silfac

#ifndef ORGCARBON
      
#ifndef PFT
      ! MOPS + CARBON
      PARAMETER(bgc_ntracer=11)
#else
      ! MOPS + CARBON + PFT
      INTEGER iphy2,iphyc2
      PARAMETER(bgc_ntracer=13)
#endif
! ends PFT      

#else
! else ORGCARBON
      
#ifndef PFT
      ! MOPS + CARBON + ORGCARBON
      PARAMETER(bgc_ntracer=13)
      INTEGER idoc,ipoc
      PARAMETER(idoc=12,ipoc=13)
#else
      ! MOPS + CARBON + ORGCARBON + PFT
      PARAMETER(bgc_ntracer=15)
      INTEGER idoc,ipoc,iphy2,iphyc2
      PARAMETER(idoc=12,ipoc=13,iphy2=14,iphyc2=15)
#endif
      
#endif
! ends ORGCARBON

#endif
! ends CARBON

      
c the tracer field
      REAL*8 bgc_tracer

      COMMON/BGC/bgc_tracer(bgc_ktotal,bgc_ntracer)

C EVERYTHING RELATED TO THE BIOGEOCHEMISTRY EVALUATION

c the flux attenuation curve
      real*8 wdet(bgc_ktotal)

c the biogeochemistry constants; array for PFTs (KM, 8/2021)
      real*8 rcp,rnp,ro2ut,subox,tempB,
ckm     &       acmuphy,acik,ackpo4,ackchl,ackw,aclambda,acomni,
     &       rcp_det_ez(maxpft), rcp_zoo,
     &       acmuphy(maxpft),ackpo4(maxpft),
     &       aclambda(maxpft),acomni(maxpft),
     &       acik,ackchl,ackw,
     &       ACMuzoo,ACkphy,AClambdaz,AComniz,ACeff,
ckm     &       graztodop,dlambda,plambda,zlambda,detlambda,
     &       graztodop,dlambda,plambda(maxpft),zlambda,detlambda,
     &       detmartin,detwa,detwb,detwmin,
     &       vsafe,alimit

c the air-sea gas exchange constants
      real*8 sox1,sox2,sox3,sox4,oA0,oA1,oA2,oA3,oA4,oA5,
     &       oB0,oB1,oB2,oB3,oC0

      COMMON/PFTPARAMS/ npft
      COMMON/BGCZ/wdet
      COMMON/BGCPARAMS/rcp,rnp,ro2ut,subox,tempB,
     &       acmuphy,acik,ackpo4,ackchl,ackw,aclambda,acomni,
     &       ACMuzoo,ACkphy,AClambdaz,AComniz,ACeff,
     &       graztodop,dlambda,plambda,zlambda,detlambda,
     &       detmartin,detwa,detwb,detwmin,
     &       vsafe,alimit
 
      COMMON/O2SURFACE/sox1,sox2,sox3,sox4,oA0,oA1,oA2,oA3,oA4,oA5,
     &       oB0,oB1,oB2,oB3,oC0

c sediment burial and O2 sensitivity of OM degradation
      real*8 burdige_fac,burdige_exp,flux_bury,ACkbaco2
      COMMON/BGCSEDPARAMS/burdige_fac,burdige_exp,flux_bury,ACkbaco2

c parameters related to N-Fixation and denitrification
      real*8 tf2,tf1,tf0,tff,nfix,subdin,rhno3ut,ACkbacdin
      COMMON/BGCNPARAMS/tf2,tf1,tf0,tff,nfix,subdin,rhno3ut,ACkbacdin

c burial of sedimentary organic carbon
c added by T.Tanioka (Nov 2020)
      real*8 flux_bury_c
      COMMON/BGCSEDPARAMS/flux_bury_c

c sinking of organic carbon
c added by T.Tanioka (Nov 2020)
      real*8 detmartin_c,detwa_c,detwb_c,detwmin_c
      COMMON/BGCPARAMS/detmartin_c,detwa_c,detwb_c,detwmin_c
      real*8 wdet_c(bgc_ktotal)
      COMMON/BGCZ/wdet_c

c parameters related C:P power-law
c added by T.Tanioka (Nov 2020)
      integer cp_option
ckm      real*8 par_bio_pc0,par_bio_po4_ref,par_bio_no3_ref,
      real*8 par_bio_pc0(maxpft),par_bio_po4_ref,par_bio_no3_ref,
     &       par_bio_temp_ref,par_bio_light_ref,
c     &       par_bio_spc_p,par_bio_spc_n,par_bio_spc_i,par_bio_spc_t,
     &       par_bio_spc_p(maxpft),par_bio_spc_n(maxpft),
     &       par_bio_spc_i(maxpft),par_bio_spc_t(maxpft),
     &       maxcp,mincp,
     &       par_zoo_cp_hom,par_excrtodoc
      COMMON/BGCPLPARAMS/par_bio_pc0,par_bio_po4_ref,par_bio_no3_ref,
     &       par_bio_temp_ref,par_bio_light_ref,
     &       par_bio_spc_p,par_bio_spc_n,par_bio_spc_i,par_bio_spc_t,
     &       maxcp,mincp,
     &       par_zoo_cp_hom,par_excrtodoc,
     &       cp_option
