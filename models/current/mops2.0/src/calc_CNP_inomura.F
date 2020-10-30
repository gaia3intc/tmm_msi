! calc_CNP_inomura.F
      
! A subroutine to calculate C:N:P using model of Inomura et al.,2020
! Pragnostically determines the growth rate and nutrient limitation 
! using MOPS ouptut
! Use parameterization for Synechococcus linearis  
! Tatsuro Tanioka, 200702
! Modified 200722 to match satellite C:N:P retrieval algorithm

      subroutine calc_cnp_inomura(PO4,DIN,PAR,T,mu,
     &     flightlim,fnutlim,cp_phy,cn_phy,np_phy)
      implicit none
      
#include "BGC_PARAMS.h"      
#include "BGC_CONTROL.h"     
! Input variables      
      real*8 PO4,DIN,PAR,T,mu
      real*8 flightlim,fnutlim
      integer PLim,NLim
      real*8 limnut
! coefficients in Inomura model
      ! Size independent parameters (ref. ST4 of I20)
      real*8 E,QCDNA,QPminRNA,QPmax,Pin(2),Nin(2)
      ! Size dependent Parameters for Synechococcus linearis (ref. ST4 of I20)
      real*8 m,vImax,AI,APho,ABio,QCProOther,ARNAP,APhoPChl
      real*8 QPOther0,QNStomax,QCOther0
      ! Stoichometry of key macromolecular pools (ref. T1 of I20)
      real*8 YChl_CNP(3),YPro_CNP(3),YRNA_CNP(3),YDNA_CNP(3)
      real*8 YPlip_CNP(3),YCSto_CNP(3),YNSto_CNP(3),YPSto_CNP(3)
! Variable in Inomura model
      real*8 I,vI,AChl,BChl,QCChl,chlc
      real*8 aM,bM,QCOther,cM,muImax,determinant
      real*8 aN,bN,cN,QN_nonsto
      real*8 aP,bP,cP,QP_nonsto
      real*8 Ccell_N,Ccell_P,Ccell
      real*8 QPPot,QP_Nlim,QPSto,QNPOt,QNmax,QN_Plim,QNSto
      real*8 muImax_chlc,muImax_QN,muImax_QP_nonsto
      real*8 muImax_Ccell_N,muImax_Ccell_P
      real*8 muImax_QPPot,muImax_QP_Nlim,muImax_QPSto
      real*8 muImax_QNPot,muImax_QN_Plim
      real*8 muImax_nc,muImax_pc,muImax_np
      real*8 nc,pc
! output C:N:P
      real*8 pc_phy, nc_phy, cp_phy, cn_phy, np_phy

! Assigning parameters values (hard-wired)
      ! Size independent parameters (ref. ST4 of I20)
      E = 7.74e-1                  ! respiratory cost of synthesis [dimensionless]
      QCDNA = 9.41e-4              ! carbon quota in DNA [mol C molC-1]
      QPminRNA = 2.23e-4           ! minimum P qupta in RNA [mol PmolC-1]
      QPmax = 9.79e-3              ! Maximum P quota [mol P mol C-1]
      Pin = (/2.00e-3, 2.00e-2 /)   ! P input for P limiting and P-replete culture [molP m-3]
      Nin = (/5.00e-2, 2.00e-1 /)   ! N input for N limiting and N-replete culture [molN m-3]
      ! Size dependent Parameters for Synechococcus linearis (ref. ST4 of I20)
      m = 3.93e-1
      vImax = 2.77e2
      AI = 8.63e-3
      APho = 1.60e1
      ABio = 2.71e-1
      QCProOther = 2.40e-1
      ARNAP = 4.23e-3
      APhoPChl = 2.82e-2
      QPOther0 = 6.53e-4
      QNStomax = 3.50e-2
      QCOther0 = 1.82e-2
      ! Stoichometry of key macromolecular pools (ref. T1 of I20)
      YChl_CNP = (/55.0, 4.0, 0.0/)     ! Chlorophyll C:N:P = 55:4:0
      YPro_CNP = (/4.49, 1.0, 0.0/)   ! Protein C:N:P = 4.49:1:0
      YRNA_CNP = (/10.7, 3.8, 1.0/) ! RNA C:N:P = 10.7:2.8:1
      YDNA_CNP = (/11.1, 3.8, 1.0/) ! DNA C:N:P = 11.1:3.8:1
      YPlip_CNP = (/40.0, 0.0, 1.0/)  ! Phopsholipid C:N:P = 40:0:1
      YCSto_CNP = (/1.0, 0.0, 0.0/)     ! Carbohydrate and non-phospholipid C:N:P = 1:0:0
      YNSto_CNP = (/2.0, 1.0, 0.0/)     ! Cyanophycin C:N:P = 2:1:0
      YPSto_CNP = (/0.0, 0.0, 1.0/)     ! Polyphosphate C:N:P = 0:0:1
! ----------------------------------------------------------------------      
c Determining nutrient limitation       
       limnut = MIN(PO4,DIN/rnp)            
      if(PO4.lt.DIN/rnp) then        ! P-limited Case
          PLim = 1
          NLim = 0
      elseif(PO4.ge.DIN/rnp) then     ! N-limited Case
          PLim = 0
          NLim = 1
      endif

      if(flightlim.gt.fnutlim) then ! Light replete case-> P and N limitation
          PLim = 1
          NLim = 1
      endif
c Setting P Lim purposefully
c      PLim = 1
! ----------------------------------------------------------------------      

c Model Equations
! Converting irradiance in W m-2 to PDF in umol m-2 s-1 (Morel and Smith, 1974)
      I = PAR*0.57
!      print*,'I in wm-2,PFD =', PAR, I
! Photosynthesis as a function of light intensity (eqn 11)
      vI = vImax*(1.0-exp(-1.0*AI*I))
! Chl:C ratio as a function of growth rate and light intensity (eqn 17)
      AChl = (1.0+E)/vI
      BChl = m/vI
      QCChl = AChl*mu + BChl
      chlc = QCChl
!      print*,'chl:C=', chlc
! Maximum Growth rate (eqn 30)       
      aM = YRNA_CNP(1)/YRNA_CNP(3)*ARNAP*(APho*AChl+ABio)
      bM = (1.0 + APho +
     &    YPlip_CNP(1)/YPlip_CNP(3)*APhoPChl)*AChl + ABio +
     &    YRNA_CNP(1)/YRNA_CNP(3)*ARNAP*(APho*BChl +
     &    QCProOther)
      QCOther =  QCProOther + QCDNA + QCOther0
      cM = (1.0 + APho +
     &    YPlip_CNP(1)/YPlip_CNP(3)*APhoPChl)*BChl + QCOther +
     &     YRNA_CNP(1)/YRNA_CNP(3)*QPminRNA - 1.0
      if(cM.lt.0.0d0) then ! Real positive solution
         muImax = (-1.0*bM + sqrt(bM**2 - 4.0*aM*cM))/(2.0*aM)
      else
          if (mu.gt.0.0d0) then ! Assign in-situ growth as its maximum
              muImax = mu
          else
              muImax = 0.0d0
          endif
      endif

      muImax = min(max(mu,muImax),3.0d0)        ! Make sure that max growth rate is greater than the actual growth rate but less than 3.0d-1
      mu = min(mu,muImax)                       ! Make sure that growth rate is less than the maximum

! Modifications to accomodate for very low light regime
! AChl and BChl are adjusted to their maximum theoretical values
      if(cM.ge.0.0d0) then
          cM = 0.0d0
          BChl = (1.0 - YRNA_CNP(1)/YRNA_CNP(3)*QPminRNA - QCOther)/
     &  (1.0 + APho + YPlip_CNP(1)/YPlip_CNP(3)*APhoPChl)
          vI = m/BChl
          AChl = (1.0+E)/vI
      endif

!      print*,'mu, mumax =', mu, muImax
!  Tentative N:C ratio under N limitation (i.e. when QNsto = 0.0)(eqn 33)       
      aN = YRNA_CNP(2)/YRNA_CNP(3)*ARNAP*(ABio + APho*AChl)
      bN = YChl_CNP(2)/YChl_CNP(1)*AChl +
     &  YPro_CNP(2)/YPro_CNP(1)*(ABio+APho*AChl) +
     &  YRNA_CNP(2)/YRNA_CNP(3)*ARNAP*(APho*BChl+QCProOther)
      cN = YChl_CNP(2)/YChl_CNP(1)*BChl +
     &  YPro_CNP(2)/YPro_CNP(1)*(APho*BChl+QCProOther) +
     & YRNA_CNP(2)/YRNA_CNP(3)*QPminRNA+YDNA_CNP(2)/YDNA_CNP(1)*QCDNA 
      QN_nonsto = aN*mu**2 + bN*mu + cN
!      print*,'C:N nonsto=',1.0/QN_nonsto
!  Tentative P:C ratio under P limitation (i.e. when QPsto = 0.0)(eqn 37)
      aP = ARNAP*(ABio+ APho*AChl)
      bP = ARNAP*(APho*BChl + QCProOther) +
     &  APhoPChl*AChl
      cP = QPminRNA + YDNA_CNP(3)/YDNA_CNP(1)*QCDNA + APhoPChl*BChl
     & + QPOther0
      QP_nonsto = aP*mu**2 + bP*mu + cP
!      print*,'C:P nonsto=',1.0/QP_nonsto
! Calculating N and P storage
      Ccell_N = Nin(1)/QN_nonsto ! Celular Carbon Quota under N-limitation (eqn 47)
      Ccell_P = Pin(1)/QP_nonsto ! Cellular Carbon Quota and P-limitation (eqn 48)
      Ccell = MIN(Ccell_N,Ccell_P) ! Actual Cellular C concentration (eqn 50)
      QPPot = Pin(2)/Ccell_N ! Potential level of cellular P under N-limited&P-replete condition (eqn 51)
      QP_Nlim = MIN(QPPot,QPmax)  ! Cellular P:C under N-limited&P-replete condition (eqn 52)
      QPSto = QP_Nlim - QP_nonsto ! Cellular P:C of storage component (eqn53)
      QNPot = Nin(2)/Ccell_P ! Potential level of cellular N under P-limited&N-replete condition (eqn 54)
      QNmax = QNStomax+ QN_nonsto ! Maximum N quota per carbon
      QN_Plim = MIN(QNPot,QNmax)        ! Cellular N:C under P-limited&N-replete condition (eqn 55)
      QNSto = QN_Plim - QN_nonsto       ! Cellular N:C of storage component (eqn 56)
! Calculating the true N:C and P:C
      if(NLim.eq.1) then
          nc = QN_nonsto
      else
          nc = QN_Plim
      endif
!     print*,'C:N true=',1.0/nc
      if(PLim.eq.1) then
          pc = QP_nonsto
      else
          pc = QP_Nlim
      endif
!      print*,'C:P true,PLim=',1.0/pc,PLim
! Calculating C:N:P:Chl at maximum growth rate 
      muImax_chlc = AChl*muImax + BChl               ! Chl:C at maximum growth
      muImax_QN = aN*muImax**2 + bN*muImax + cN      ! No N storage at maximum growth 
      muImax_QP_nonsto = aP*muImax**2 + bP*muImax + cP

      muImax_Ccell_N = Nin(1)/muImax_QN
      muImax_Ccell_P = Pin(1)/muImax_QP_nonsto

      muImax_QPPot = Pin(2)/muImax_Ccell_N
      muImax_QP_Nlim = MIN(muImax_QPPot,QPmax)
      muImax_QPSto = muImax_QP_Nlim - muImax_QP_nonsto

      muImax_QNPot = Nin(2)/muImax_Ccell_P
      muImax_QN_Plim = MIN(muImax_QNPot,muImax_QN)
      if(muImax_QNPot.lt.muImax_QN) then                    ! if QN is less than potential, assign potential as QN
          muImax_QN = muImax_QNPot
      endif
      muImax_nc = muImax_QN;
      if(PLim.eq.1) then
          muImax_pc = MIN(QPmax,muImax_QP_nonsto)
      else
          muImax_pc = muImax_QP_Nlim
      endif
! Final check to make sure N:C and P:C are within the theoretical limits
      nc = MIN(muImax_nc,nc)
      pc = MIN(muImax_pc,pc)
! Final Output      
! If fault, output Redfield C:N:P 
      if (mu.ge.0.0) then
        if(pc.gt.0.0) then
          cp_phy = 1.0/pc
        else
          cp_phy = 117.0d0
        endif
        if(nc.gt.0.0d0) then
          cn_phy = 1.0/nc
        else
          cn_phy = 117.0d0/16.0d0
        endif
      else
          cp_phy = 117.0d0
          cn_phy = 117.0d0/16.0d0
      endif
      np_phy = cp_phy/cn_phy
      
!      print*,'C:N,C:P,N:P =',cp_phy,cn_phy,np_phy
      end

      
