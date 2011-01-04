C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE STRPSRB(LATMOS,TS,ROFF,SFCALB,VIS,THETA,ICLOUD,IWP,
     & ISUB,TAUW,TOP,NCL,INTVAL,tstat,rf,zstat,CLOLWC,CDR,TRANS,TDIR)
C
C     calculate transmittance for clear and cloudy sky
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     LATMOS  : atmosphere type
C               1 - tropical
C               2 - midlatitude summer
C               3 - midlatitude winter
C               4 - subartic summer
C               5 - subartic winter
C     TS      : difference between surface and characteristical temperatures
C     ROFF    : difference in water vapor
C     SFCALB  : surface albedo (0/1)
C     VIS     : visibility at ground station [km]
C     THETA   : zenith angle (o)
C     ICLOUD  : cloud parameter
C              0 - no clouds
C              1 - calculation with clouds
C     IWP     : cloud droplet size distribution
C     ISUB    : use subroutine WOLKE1 (ISUB=1) or WOLKE2 (ISUB=2)
C     TAUW    : cloud optical thickness
C              0.0 - no clouds 
C              100.0 - with clouds 
C     TOP     : cloud top (mbar)
C     NCL     : number of cloud layers
C     INTVAL  : spectral interval
C     INTVA   : value of the begin of spectral interval
C     INTVE   : value of the end of spectral interval
C     WH2O    : precipitable water [cm]
C     CLOLWC  : cloud liquid water content
C     CDR     : change degrees to radians
C     TRANS   : parameterized sky transmittance for broadband solar radiation 
C              in the atmosphere 
C     INDEX   : index of spectral intervals
C     NLEV    : number of atmosphere levels
C     SOL     : solar constant for each spectral interval
C     DZ      : thickness of each atmospheric layer
C     ROL     : air density
C     ROH     : water content (humidity)
C     RO3     : ozone content 
C     NLAY    : number of atmospheric layers
C     IH      : data for index of spectral intervals
C     LCLOUD  : logical variable to calculate with or without clouds
C     PRESS   : atmospheric pressure (mbar)
C     ALT     : altitude (m)
C     TEMP    : temperature (oC)
C     SCA     : 
C     RAY     :
C     IGO     : flag for CLOLWC greater than limits
C     TW      : 
C     K1      :
C     DIR     : direct irradiation for each spectral interval (Wh/m2)
C     DIF     : diffuse irradiation for each spectral interval(Wh/m2)
C     DOWN    : global irradiation (Wh/m2)
C     DIRS    : direct irradiation (Wh/m2)
C     DIFS    : diffuse irradiation (Wh/m2)
C     IBEG    : begin of spectral interval
C     ISTO    : end of spectral interval
C     IALL    : counter of spectral interval
C     IDO     : flag to execute the loop of spectral intervals
C     NLA     : NLEV-1
C     NLY     : NLAY-1
C     ALDIF   : surface albedo for cloudy sky
C     ALDR    : surface albedo for clear sky
C     COSZEN  : cossine of zenith angle
C     LEVMAX  : maximum number of atmopheric levels
C     LAYMAX  : maximum number of atmopheric layers
C     MAXWEL  : number of spectral intervals
C     MAXIND  : number of spectral intervals including their subdivisions
C     IDO     :
C     BET     :
C     BET0    :
C     OM      :
C     TAU     :
C     TOHNW   :
C     A       :
C     AA      :
C     DW      :
C     INAERO  :
C     ASIG    : 
C     U       :
C     V       :
C     WH2O1   :
C     WCO2SW  :
C     WO3SW   :
C     WO2SW   :
C     RAYM    :      
C......................................................................
C
      PARAMETER (LEVMAX = 51, LAYMAX = 50,MAXWEL=37 , MAXIND= 135)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION U(LAYMAX,MAXIND),V(LAYMAX,MAXIND)
      DIMENSION INAERO(LAYMAX),ASIG(LAYMAX,MAXWEL)
      DIMENSION SCA(LEVMAX),RAY(LEVMAX)
      DIMENSION DOWN(LEVMAX),DIR(LEVMAX),DIF(LEVMAX),
     &  TEMP(LEVMAX),PRESS(LEVMAX),ALT(LEVMAX),ROL(LEVMAX),
     &  ROH(LEVMAX),RO3(LEVMAX),DIRS(LEVMAX),DIFS(LEVMAX)
      DIMENSION DZ(LAYMAX),TW(LAYMAX),TAU(LAYMAX),TOHNW(LAYMAX),
     &          BET(LAYMAX),BET0(LAYMAX),OM(LAYMAX)
      DIMENSION SOL(MAXIND),INDEX(MAXIND)
      DIMENSION IH(113),INTVA(8),INTVE(8)
      DIMENSION WH2O1(LAYMAX)
      DIMENSION WCO2SW(LAYMAX),WO3SW(LAYMAX),WO2SW(LAYMAX),
     &          RAYM(LAYMAX)
C
      LOGICAL* 1 LCLOUD
C
      DATA IH / 5*23, 3*24, 25, 5*26, 27, 8*28, 29, 9*30, 31, 10*32
     &        , 2*33, 20*34, 35, 36*36, 10*37 /
C
      DATA INTVA / 20, 16,  1,  0,  0,  0,  0,  1/
C      
      DATA INTVE / 29, 30, 37,  0,  0,  0,  0, 37/
C
      DO 10   I =  1 ,  22
   10 INDEX(I) = I
      DO 20   I = 23 , 135
   20 INDEX(I) = IH(I-22)
C
      LCLOUD = .FALSE.                                                 
      IF (ICLOUD .GT. 0) THEN
      LCLOUD = .TRUE.
      ENDIF
      print *, '1'
C
C......................................................................
C
C     subroutine LESEN - determination of atmospheric profiles
C     input  - LATMOS,VIS,INDEX
C     output - NLEV,SOL,DZ,PRESS,ALT,TEMP,ROL,ROH,RO3,SCA,RAY
C
      CALL LESEN(LATMOS,VIS,INDEX,NLEV,SOL,DZ,PRESS,ALT,TEMP,ROL,ROH,
     & RO3,SCA,RAY)          
      print *, '2'               
C
C     MODIFICACOES SUGERIDAS PELO FERNANDO PARA PREVER A AGUA PRECIPITAVEL 
C     A DIFERENTES ALTITUDES
C
C     calculation of the precitable water as function of relative
C     humidity, partial pressure of water vapor and temperature
C
      XD=287.05/1005.0
      PSTAT = 1013.25*EXP(-0.0001184*ZSTAT)
      PVSAT = EXP(26.23 - 5416.0/TSTAT)/100.0
      WMIXRAT=0.622*RF*PVSAT/1013.25
      XM=(1.0-0.26*WMIXRAT)*XD              !4
      TNEW=TSTAT*((PSTAT/1013.25)**XM)      !3
C     WRITE(*,*) 'tnew=',TNEW,'tstat=',TSTAT,'dalt=',dalt
      PVSAT = EXP(26.23 - 5416.0/TNEW)      !2
      WH2O  = 0.493*RF*PVSAT/TNEW           !1
C
C     test for precitable water less than 0.0000001 cm
C
      IF (WH2O .LE. 0.0000001) THEN
      WRITE (*,*)'NEGATIVE PRECIPITABLE WATER VAPOR',WH2O,
     & ' IS SET TO 0.000001 (cm)',TSTAT,ZSTAT
      WH2O  = 0.0000001
      ENDIF
      print *, '3'
C
C     TERMINO DAS MODIFICACOES
C
C     
C......................................................................
C
      NLAY=NLEV-1
C
      IGO = 1
      IF(ISUB.EQ.1) THEN
C
C......................................................................
C
C     subroutine WOLKE1 - determination of clouds properties, input is
C                         CLOLWC
C     input  - NCL,NLAY,TOP,IWP,LCLOUD,CLOLWC,DZ,TEMP,PRESS
C     output - NCL,IGO,TAUW,TW,K1
C
      CALL WOLKE1(NCL,NLAY,IGO,TOP,TAUW,IWP,LCLOUD,CLOLWC,
     &   DZ,TEMP,PRESS,TW,K1)
      print *, '4'
C
C......................................................................
C
      ELSE
C
C......................................................................
C
C     subroutine WOLKE2 - determination of clouds properties, input is
C                         TAUW
C     input  - NCL,NLAY,TOP,TAUW,IWP,LCLOUD,DZ,TEMP,PRESS
C     output - NCL,CLOLWC,TW,K1
C
         CALL WOLKE2(NCL,NLAY,TOP,TAUW,IWP,LCLOUD,CLOLWC,
     &   DZ,TEMP,PRESS,TW,K1)
C
C......................................................................
C
      END IF
C
      IF(IGO.NE.1) RETURN
C
C......................................................................
C
C     subroutine VARIA - changes atmospheric profiles
C     input  - RAY,NCL,NLAY,TS,ROFF,WH2O,K1,LCLOUD,DZ,TEMP,PRESS,ROL,ROH,RO3
C     output - TEMP,ROH,WH2O1,WCO2SW,WO3SW,WO2SW,RAYM
C     
      print *, '5'
      CALL VARIA(RAY,NCL,NLAY,TS,ROFF,WH2O,K1,LCLOUD,DZ,TEMP,PRESS,ROL,
     &           ROH,RO3,WH2O1,WCO2SW,WO3SW,WO2SW,RAYM)
C
C......................................................................
C
C     subroutine AERO - determination of aerosol profiles
C     input  - SCA,NLAY,DZ,ALT
C     output - INAERO,ASIG
C
      print *, '6'
      CALL AERO(SCA,NLAY,DZ,ALT,INAERO,ASIG)
C
C......................................................................
C
C     subroutine ATMOS - determination of absorption and scattering
C     input  - NLAY,INDEX,DZ,WH2O1,WCO2SW,WO3SW,WO2SW,RAYM
C     output - U,V
C
      print *, '7'
      CALL ATMOS(NLAY,INDEX,DZ,WH2O1,WCO2SW,WO3SW,WO2SW,RAYM,U,V)
C
C......................................................................
C
      THETA = DMAX1(THETA,0.1D0)
      COSZEN= DCOS(CDR*THETA) 
C
C     condition to avoid that some parts of the map be in the night      
C                                                                  
      IF(COSZEN.LT.0.01) COSZEN=0.01
C
      DO 30  I=1,LEVMAX
        DIR(I)    = 0.
        DIF(I)    = 0.
        DOWN(I)   = 0.
        DIRS(I)   = 0.
        DIFS(I)   = 0.
   30 CONTINUE
C
C.....loop spectral intervals..........................................
C
      print *, '8'
      IBEG = INTVA(INTVAL)
      ISTO = INTVE(INTVAL)
      IF(IBEG.EQ.0.AND.ISTO.EQ.0) GOTO 100
C
      DO 40 IALL = 1,135
      IDO = INDEX(IALL)
      IF(IDO.LT.IBEG.OR.IDO.GT.ISTO) GOTO 40
C
      NLA = NLEV
      NLY = NLAY
C
      ALDIF = SFCALB
      ALDR  = SFCALB
C
C......................................................................
C
C     subroutine BETAS - determination of spectral absorption, extinction 
C                        and scattering for each atmospheric layer
C     input  - NCL,NLAY,IALL,COSZEN,INDEX,TW,K1,IWP,LCLOUD,INAERO,ASIG,U,V
C     output - BET,BET0,OM,TAU,TOHNW
C
      CALL BETAS (NCL,NLAY,IALL,COSZEN,INDEX,BET,BET0,OM,TAU,TOHNW,
     & TW,K1,IWP,LCLOUD,INAERO,ASIG,U,V)
C
C......................................................................
C
      DIR(1)  = COSZEN
      A       = 0.
      AA      = 0.
C
      DO 50  I = 2 , NLEV
        A       =  A + TOHNW(I-1)
        AA      =  AA + TAU(I-1)
        DW      =  AA/COSZEN
        DW      =  DMIN1(100.0D0,DW)
        DIR(I)  =  COSZEN*DEXP(-DW)
        IF (A .GT. 25.)  GOTO 200
        IF (I .EQ. NLEV)  GOTO 300
   50 CONTINUE
C
  200 NLA     = I-1
      NLY     = NLA-1
      ALDIF   = 0.
      ALDR    = 0.
  300 CONTINUE
      IF (NLY.GT.0) THEN
C
C......................................................................
C
C     subroutine MATBAU - calculation of diffuse irradiation
C                         for each spectral interval
C     input  - NLY,ALDIF,ALDR,COSZEN,BET,BET0,OM,TAU,DIR
C     output - DIF
C
      CALL MATBAU(NLY,ALDIF,ALDR,COSZEN,BET,BET0,OM,TAU,DIR,DIF)
C
C......................................................................
C
      ENDIF
C
C......................................................................
C
C   upward (UP), downward diffuse (DIF)                    
C   and downward direct (DIR) fluxes                       
C   have been computed for unity input fluxes for each      
C   layer boundary for the wave-length interval IALL.      
C   The true fluxes are obtained by multiplication with the 
C   spectral solar constant (input from <LESEN>)  SOL .     
C.......................................................................
C                                                           
C   for each boundary NLA spectral integrated quantities 
C
C     DIF  : downward diffuse flux                     
C     DIR  : downward direct flux                           
C                                                           
C   are computed                                            
C.......................................................................

      DO 60  I = 1,NLA
        DIF(I) = DIF(I)*SOL(IALL)
        DIR(I) = DIR(I)*SOL(IALL)
        DIRS(I) = DIRS(I)+DIR(I)
        DIFS(I) = DIFS(I)+DIF(I)
   60 CONTINUE

   40 CONTINUE    !end of loop - 40 -  spectral intervals

      DO 70   I =  1 , NLEV
       DOWN(I)    = DIRS(I) + DIFS(I)
   70 CONTINUE

  100 CONTINUE           

C     calculation of atmospheric transmittance

      TRANS=DOWN(NLEV)/DOWN(1)
      TDIR=DIRS(NLEV)/DOWN(1)

      RETURN
      END

C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE LESEN(LATMOS,VIS,INDEX,NLEV,SOL,DZ,PRESS,ALT,TEMP,ROL,
     & ROH,RO3,SCA,RAY)

C            determination of atmospheric profiles

C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C     LATMOS  : atmosphere type
C               1 - tropical
C               2 - midlatitude summer
C               3 - midlatitude winter
C               4 - subartic summer
C               5 - subartic winter
C     VIS     : visibility at ground station [km]
C     INDEX   : index of spectral intervals
C     NLEV    : number of atmosphere levels
C     SOL     : solar constant for each spectral interval
C     DZ      : thickness of each atmospheric layer
C     ROL     : air density
C     ROH     : water content (humidity)
C     RO3     : ozone content 
C     PRESS   : atmospheric pressure (mbar)
C     ALT     : altitude (m)
C     TEMP    : temperature (oC)
C     SCA     : 
C     RAY     :
C     RAYL    :
C     SCAL    : 
C     ISOL    : 
C     WEIGHT  :
C     FITTOP  : top of aerosol profile (km)
C     FITBOT  : bottom of aerosol profile (km)
C     JTOP    :
C     JBOT    :
C     H       : altitude of aerosol profile (km)
C     DLAY    : delta of atmospheric layer (km)
C     SCATOP :
C     SCABOT :
C     CSCAL   :
C     ZJ      : altitude (km)
C     NLAY    : number of atmospheric layers
C     M1      :
C     M2      :
C     CH      :
C     CRAY    :
C     CSC     :
C     X       :
C......................................................................

      PARAMETER (LEVMAX = 51, LAYMAX = 50, MAXWEL= 37, MAXIND= 135)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION RAYL(LEVMAX),SCAL(LEVMAX)
      DIMENSION RAY(LEVMAX),SCA(LEVMAX)
      DIMENSION TEMP(LEVMAX),PRESS(LEVMAX),ALT(LEVMAX),
     &          ROL(LEVMAX),ROH(LEVMAX),RO3(LEVMAX)
      DIMENSION DZ(LAYMAX)
      DIMENSION SOL(MAXIND)
      DIMENSION ISOL(37),WEIGHT(135),INDEX(MAXIND)

      DATA RAYL /
     &                               2.1352E+1,2.4178E+1,2.7378E+1,
     & 3.1117E+1,3.5642E+1,4.0885E+1,4.6968E+1,5.4039E+1,6.2270E+1,
     & 7.1869E+1,8.3082E+1,9.6204E+1,1.1159E+2,1.2965E+2,1.5091E+2,
     & 1.7598E+2,2.0559E+2,2.4064E+2,2.8185E+2,3.2835E+2,3.8280E+2,
     & 4.4660E+2,5.2141E+2,6.0920E+2,7.1230E+2,8.3346E+2,9.7598E+2,
     & 1.1437E+3,1.3414E+3,1.5743E+3,1.8487E+3,2.1624E+3,2.5294E+3,
     & 2.9589E+3,3.4614E+3,4.0495E+3,4.7378E+3,5.5433E+3,6.4861E+3,
     & 7.5853E+3,8.5981E+3,9.7116E+3,1.0933E+4,1.2268E+4,1.3726E+4,
     & 1.5313E+4,1.7037E+4,1.8906E+4,2.0929E+4,2.3115E+4,2.5471E+4
     & /

      DATA SCAL /
     &                                              2.10E-6, 2.39E-6,
     & 2.72E-6, 3.10E-6, 3.53E-6, 4.02E-6, 4.61E-6, 5.29E-6, 6.07E-6,
     & 6.97E-6, 8.00E-6, 9.25E-6, 1.07E-5, 1.24E-5, 1.43E-5, 1.65E-5,
     & 1.90E-5, 2.18E-5, 2.51E-5, 2.89E-5, 3.32E-5, 4.49E-5, 6.07E-5,
     & 8.21E-5, 1.11E-4, 1.50E-4, 1.90E-4, 2.42E-4, 2.96E-4, 3.52E-4,
     & 4.23E-4, 4.92E-4, 5.63E-4, 6.01E-4, 6.41E-4, 6.43E-4, 6.45E-4,
     & 6.22E-4, 6.63E-4, 7.14E-4, 7.87E-4, 9.80E-4, 1.41E-3, 2.30E-3,
     & 3.54E-3, 4.85E-3, 6.43E-3, 8.19E-3, 9.70E-3, 2.58E-2, 6.95E-2
     & /

      DATA WEIGHT /
     & 22*1.00000,
     & 0.012690, 0.114010, 0.004610, 0.017630, 0.851060, 0.364000,
     & 0.284000, 0.352000, 1.000000, 0.002963, 0.104164, 0.005804,
     & 0.018354, 0.868714, 1.000000, 0.001495, 0.074495, 0.265915,
     & 0.024395, 0.008525, 0.091175, 0.040895, 0.493105, 1.000000,
     & 0.018860, 0.001920, 0.084790, 0.168890, 0.017250, 0.010730,
     & 0.108690, 0.035210, 0.553660, 1.000000, 0.055750, 0.043240,
     & 0.165090, 0.017240, 0.113830, 0.070110, 0.033240, 0.081810,
     & 0.124930, 0.294760, 0.977800, 0.022200, 0.028100, 0.023680,
     & 0.082260, 0.157890, 0.013810, 0.058030, 0.025950, 0.031770,
     & 0.096170, 0.283030, 0.006990, 0.005900, 0.020470, 0.039300,
     & 0.003440, 0.014450, 0.006460, 0.007910, 0.023940, 0.070450,
     & 1.000000, 0.011870, 0.026320, 0.031590, 0.030910, 0.096520,
     & 0.005280, 0.077550, 0.122630, 0.065720, 0.051650, 0.111230,
     & 0.183620, 0.002140, 0.004740, 0.005690, 0.005560, 0.017380,
     & 0.000950, 0.013960, 0.022080, 0.011830, 0.009300, 0.020020,
     & 0.033060, 0.000560, 0.001240, 0.001490, 0.001460, 0.004550,
     & 0.000250, 0.003650, 0.005780, 0.003100, 0.002430, 0.005240,
     & 0.008650, 0.014070, 0.080650, 0.183930, 0.007960, 0.208700,
     & 0.014700, 0.024790, 0.228760, 0.103390, 0.133050
     & /
! valores antes da adaptacao para calculo do PAR!     & 22*1.00000,
!     & 0.01269, 0.11401, 0.00461, 0.01763, 0.85106, 0.36320, 0.28320,
!     & 0.35120, 1.00000, 0.01059, 0.10226, 0.00390, 0.01645, 0.86681,
!     & 1.00000, 0.00137, 0.07437, 0.26579, 0.02427, 0.00840, 0.09105,
!     & 0.04177, 0.49298, 1.00000, 0.01886, 0.00192, 0.08479, 0.16889,
!     & 0.01725, 0.01073, 0.10869, 0.03521, 0.55366, 1.00000, 0.05575,
!     & 0.04324, 0.16509, 0.01724, 0.11383, 0.07011, 0.03324, 0.08181,
!     & 0.12493, 0.29476, 0.97780, 0.02220, 0.02810, 0.02368, 0.08226,
!     & 0.15789, 0.01381, 0.05803, 0.02595, 0.03177, 0.09617, 0.28303,
!     & 0.00699, 0.00590, 0.02047, 0.03930, 0.00344, 0.01445, 0.00646,
!     & 0.00791, 0.02394, 0.07045, 1.00000, 0.01187, 0.02632, 0.03159,
!     & 0.03091, 0.09652, 0.00528, 0.07755, 0.12263, 0.06572, 0.05165,
!     & 0.11123, 0.18362, 0.00214, 0.00474, 0.00569, 0.00556, 0.01738,
!     & 0.00095, 0.01396, 0.02208, 0.01183, 0.00930, 0.02002, 0.03306,
!     & 0.00056, 0.00124, 0.00149, 0.00146, 0.00455, 0.00025, 0.00365,
!     & 0.00578, 0.00310, 0.00243, 0.00524, 0.00865, 0.01407, 0.08065,
!     & 0.18393, 0.00796, 0.20870, 0.01470, 0.02479, 0.22876, 0.10339,
!     & 0.13305 /

      DATA ISOL       /
     &     0.4 ,   16.0,   37.0 ,   55.0 ,   54.0,   60.0,    86.0,
     &   160.0 ,  177.0,  273.0 ,  757.5 , 1476.0, 1845.4,  2107.7,
     &  3246.8 , 4280.7, 9476.2 , 9601.2 , 9320.5, 8836.5,  7960.5,
     &  3667.6 , 7949.8, 1220.4 , 1791.7 , 8769.6, 1958.5, 10184.4,
     &  5094.8 , 7307.3, 2996.91,10003.79, 2591.6, 6498.2,   796.5,
     &  3813.85, 1128.3 
     & /
! valores antes da adaptacao para calculo do PAR
!     &    0.4,  16.0,  37.0,  55.0,   54.0,  60.0,  86.0, 160.0, 177.0,
!     &  273.0, 757.5,1476.0,1845.4, 2107.7,2102.5,5425.0,9476.2,9601.2,
!     & 9320.5,8836.5,7960.5,5066.2, 6551.2,1220.4,1791.7,8769.6,1958.5,
!     & 10184.4,5094.8,7307.3,2996.9,10003.8,2591.6,6498.2,796.5,3813.0,
!     & 1128.3 /

      DO 10 I=1,135
   10 SOL(I)=ISOL(INDEX(I))*0.01*WEIGHT(I)

C.....................................................................
C      calculate aerosol-profiles for a given surface visibility
C      Load Aerosol-Profile (Mc Clatchey with surface visibility=50km)
C      fit new aerosol profile between fittop(5km) and fitbot(0km) for
C      a surface visibility VIS = 5 - 50km taken from the input data
C.....................................................................

      FITTOP  = 5.0
      FITBOT  = 0.0
      JTOP    = 51 - INT(FITTOP)
      JBOT    = 51 - INT(FITBOT)
      H       = FITTOP - FITBOT
      DLAY    = 1.0
      SCATOP = SCAL(JTOP)
      SCABOT = 3.912/VIS
      CSCAL   = DLOG(SCATOP/SCABOT)/H

      ZJ = FITTOP + DLAY
      DO 20 J = JTOP,JBOT
        ZJ = ZJ - DLAY
        SCAL(J) = SCABOT*DEXP(CSCAL*(ZJ-FITBOT))
   20 CONTINUE

C ....................................................................
C     *** model atmosphere ***

      OPEN (UNIT=9,FILE=
     & 'Srbpres.dat',
     & FORM='FORMATTED',STATUS='OLD')

      READ (9,*) NLEV
      READ (9,*) (PRESS(I),I=1,NLEV)
C      WRITE(*,*) PRESS
      NLAY = NLEV - 1
      CLOSE (9)
C......................................................................

C     subroutine PTZQD - calculation of atmospheric properties for each
C                        level
C     input  - LATMOS,NLEV,PRESS
C     output - PRESS,ALT,TEMP,ROL,ROH,RO3

      CALL PTZQD(LATMOS,NLEV,PRESS,ALT,TEMP,ROL,ROH,RO3)

C......................................................................

      DO 30  I = 1 , NLEV
         ROL(I) = ROL(I)*1000.0D0
         ROH(I) = ROH(I)*1000.0D0
         RO3(I) = RO3(I)*1000.0D0
         PRESS(I) = PRESS(I)*0.001D0
   30 CONTINUE

      DO 40 I=1,NLAY
      M1 = 50 - INT(ALT(I)-0.001)
      M2 = M1 + 1
      CH = ALT(I) - AINT(ALT(I)-0.0001)
      CRAY= DLOG (RAYL(M1)/RAYL(M2))
      CSC = DLOG (SCAL(M1)/SCAL(M2))
      X = CH - 1.
      RAY(I) = RAYL(M1) * DEXP(X*CRAY)
      SCA(I) = SCAL(M1) * DEXP(X*CSC)
   40 CONTINUE

      SCA(NLEV) = SCAL(51)
      RAY(NLEV) = RAYL(51)

      DO 50 I=1,NLAY
   50 DZ(I)=DABS(ALT(I)-ALT(I+1))

      RETURN
      END

C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PTZQD(LATMOS,NLEV,PRESS,ALT,TEMP,ROL,ROH,RO3)

C     calculation of atmospheric properties for each level

C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C     LATMOS  : atmosphere type
C               1 - tropical
C               2 - midlatitude summer
C               3 - midlatitude winter
C               4 - subartic summer
C               5 - subartic winter
C     NLEV    : number of atmosphere levels
C     ROL     : air density
C     ROH     : water content (humidity)
C     RO3     : ozone content 
C     PRESS   : atmospheric pressure (mbar)
C     ALT     : altitude (m)
C     TEMP    : temperature (oC)
C     NLEVD2  :
C     UMSP    :
C     QL      :
C     RM      :
C     PTZQ1   :
C     PTZQ2   :
C     PTZQ3   :
C     DLOGP   :
C     DELZAP  :
C     NINT    : 
C     ZNINT   : 
C     G0      :
C     AA      :
C     G       :
C     ZMASS   :
C     R       :
C     DZ      :
C     HT      :
C     RKn     :
C......................................................................
C
      PARAMETER(LEVMAX=51)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      
      DIMENSION PRESS(LEVMAX),ALT(LEVMAX),TEMP(LEVMAX),ROL(LEVMAX),
     &          ROH(LEVMAX),RO3(LEVMAX)
C     
      DATA RM,DELZAP,R,G0,ZMASS,AA /348.36764D0,0.5D0,8.314D0,9.81D0,
     & 28.9D0,6400.0D0/
C     
      NLEVD2=NLEV*0.5
      DO 10 I=1,NLEVD2
      J=NLEV-I+1
      UMSP=PRESS(J)
      PRESS(J)=PRESS(I)
      PRESS(I)=UMSP
   10 CONTINUE
C
      ALT(1)=0.0D0
      TEMP(1)=PTZQ1(ALT(1),LATMOS)
      QL=PRESS(1)*RM/TEMP(1)
      ROH(1)=PTZQ2(PRESS(1),LATMOS,QL)
      ROL(1)=QL/(1.0D0+0.608D0*ROH(1)/QL)
      RO3(1)=PTZQ3(PRESS(1),LATMOS)
C
      DO 20 N=1,NLEV-1
      DLOGP=7.0D0*DLOG(PRESS(N)/PRESS(N+1))
      NINT=DLOGP/DELZAP+1
      ZNINT=DBLE(NINT)
      G=G0*(1.0D0-2.0D0*ALT(N)/AA)
      DZ=R*DLOGP/(7.0D0*ZMASS*G*ZNINT)
      HT=ALT(N)
      DO 30 M=1,NINT
      RK1=PTZQ1(HT,LATMOS)*DZ
      RK2=PTZQ1(HT+0.5*RK1,LATMOS)*DZ
      RK3=PTZQ1(HT+0.5*RK2,LATMOS)*DZ
      RK4=PTZQ1(HT+RK3,LATMOS)*DZ
      HT=HT+0.16666667D0*(RK1+RK2+RK2+RK3+RK3+RK4)
   30 CONTINUE
      ALT(N+1)=HT
      TEMP(N+1)=PTZQ1(HT,LATMOS)
      QL=PRESS(N+1)*RM/TEMP(N+1)
      ROH(N+1)=PTZQ2(PRESS(N+1),LATMOS,QL)
      ROL(N+1)=QL/(1.D0+0.608D0*ROH(N+1)/QL)
      RO3(N+1)=PTZQ3(PRESS(N+1),LATMOS)
   20 CONTINUE
      DO 40 N=1,NLEV
      ROH(N)=ROH(N)*0.001D0
      ROL(N)=ROL(N)*0.001D0
   40 CONTINUE
      DO 50 I=1,NLEVD2
      J=NLEV-I+1
      UMSP=PRESS(J)
      PRESS(J)=PRESS(I)
      PRESS(I)=UMSP
      UMSP=ALT(J)
      ALT(J)=ALT(I)
      ALT(I)=UMSP
      UMSP=TEMP(J)
      TEMP(J)=TEMP(I)
      TEMP(I)=UMSP
      UMSP=ROL(J)
      ROL(J)=ROL(I)
      ROL(I)=UMSP
      UMSP=ROH(J)
      ROH(J)=ROH(I)
      ROH(I)=UMSP
      UMSP=RO3(J)
      RO3(J)=RO3(I)
      RO3(I)=UMSP
   50 CONTINUE
      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE WOLKE1
     & (NCL,NLAY,IGO,TOP,TAUW,IWP,LCLOUD,CLOLWC,DZ,TEMP,PRESS,TW,K1)
C
C     determination of clouds properties, input is CLOLWC
C     
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     IWP     : cloud droplet size distribution
C     TAUW    : cloud optical thickness
C              0.0 - no clouds 
C              100.0 - with clouds 
C     TOP     : cloud top (mbar)
C     NCL     : number of cloud layers
C     CLOLWC  : cloud liquid water content
C     DZ      : thickness of each atmospheric layer
C     NLAY    : number of atmospheric layers
C     LCLOUD  : logical variable to calculate with or without clouds
C     PRESS   : atmospheric pressure (mbar)
C     TEMP    : temperature (oC)
C     IGO     : flag for CLOLWC greater than limits
C     K1      :
C     ALW     :
C     FLW     :
C     TW      :
C     PWZW    :
C     ICLT    :
C     ICLB    :
C     TTOP    :
C     TBOT    :
C     EXPON   :
C     ESAT    :
C     ROBOT   :
C     ROTOP   :
C     CLIM    :
C     EXMASS  :
C     EWOL    :
C     FIXLWC  :
C     IWO     :
C     KWO     :
C......................................................................
C
      PARAMETER (LEVMAX = 51, LAYMAX = 50, MAXWEL= 37, MAXIND= 135)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION EWOL(16,8)
      DIMENSION TEMP(LEVMAX),PRESS(LEVMAX)
      DIMENSION DZ(LAYMAX)
      DIMENSION FLW(LAYMAX)
      DIMENSION TW(LAYMAX),FIXLWC(8)
      LOGICAL*1 LCLOUD 
C      
      DATA FIXLWC /
     & 0.05, 0.14, 0.28, 0.47, 1.00, 2.50,
     & 0.13, 0.44
     & /
C
      DATA EWOL /
     &  15.2,  15.4,  15.7,  15.3,  15.6,  15.8,  16.0,  15.8,  16.0,
     &  16.0,  15.9,  16.5,  17.1,  17.2,  21.2,  18.1,
     &  39.4,  40.1,  41.0,  40.9,  42.3,  41.9,  41.6,  41.2,  42.1,
     &  42.9,  43.5,  43.2,  44.1,  45.8,  50.7,  47.1,
     &  72.4,  72.7,  74.1,  75.5,  73.8,  73.9,  74.6,  74.6,  75.0,
     &  75.2,  77.5,  78.4,  78.8,  80.0,  85.6,  82.4,
     &  73.9,  76.6,  74.4,  75.1,  75.5,  75.4,  74.8,  76.0,  76.1,
     &  76.2,  77.2,  77.4,  77.7,  80.3,  80.8,  81.2,
     & 127.8, 131.1, 130.3, 129.1, 130.9, 129.8, 130.9, 130.4, 130.6,
     & 130.7, 134.9, 132.8, 136.4, 137.4, 138.6, 140.2,
     & 121.1, 122.1, 121.2, 122.3, 121.0, 122.3, 121.0, 120.4, 121.4,
     & 122.3, 124.1, 123.3, 125.8, 126.4, 125.2, 127.1,
     &   5.5,   5.5,   5.5,   5.5,   5.5,   5.5,   5.5,   5.5,   5.6,
     &   5.6,   5.6,   5.6,   5.6,   5.6,   5.7,   5.7,
     &  12.3,  12.3,  12.3,  12.3,  12.4,  12.4,  12.4,  12.4,  12.4,
     &  12.4,  12.4,  12.5,  12.5,  12.6,  12.6,  12.7
     & /
C
      TAUW    =  0.0
      ALW     =  0.0
      DO 10  I = 1 , NLAY
        FLW(I) = 0.0
        TW(I)  = 0.0
   10 CONTINUE
      K1 = 0
C
      IF (.NOT. LCLOUD)                      GOTO 100
C
      DO 20  I = 1 , NLAY
        PWZW =  PRESS(I)*1000.0D0 + 1.0D0
        IF (TOP.LE.PWZW)                     GOTO 200
   20 CONTINUE
  200 K1 = I-1
C
      ICLT = K1 + 1
      ICLB = ICLT + NCL
      IF(ICLB.GT.NLAY+1) THEN
        WRITE(*,*) ' cloud bottom has to be adjusted'
        ICLB = NLAY + 1
        NCL = ICLB - ICLT
      ENDIF
C
      TTOP = TEMP(ICLT)
      TBOT = TEMP(ICLB)
      EXPON = (7.5D0*(TTOP-273.15D0))/(TTOP-35.85D0)
      ESAT = 6.11D0*10.0D0**EXPON
      ROTOP = 216.67D0*ESAT/TTOP
      EXPON = (7.5D0*(TBOT-273.15D0))/(TBOT-35.85D0)
      ESAT  = 6.11D0*10.0D0**EXPON
      ROBOT = 216.67D0*ESAT/TBOT
      CLIM  = (ROBOT-ROTOP)*0.2D0
      IF(CLIM.LT.CLOLWC) THEN
        IGO = 0
        RETURN
      ENDIF
C
      EXMASS  = EWOL(1,IWP)/FIXLWC(IWP)/1000.0D0
      DO 30  IWO = 1 , NCL
        KWO      = K1 + IWO
        FLW(KWO) = CLOLWC*DZ(KWO)*1000.0D0
        TW(KWO)  = FLW(KWO)*EXMASS
        TAUW     = TAUW+TW(KWO)
        ALW      = ALW+FLW(KWO)
   30 CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE WOLKE2
     & (NCL,NLAY,TOP,TAUW,IWP,LCLOUD,CLOLWC,DZ,TEMP,PRESS,TW,K1)
C
C     determination of clouds properties, input is CLOLWC
C     
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     IWP     : cloud droplet size distribution
C     TAUW    : cloud optical thickness
C              0.0 - no clouds 
C              100.0 - with clouds 
C     TOP     : cloud top (mbar)
C     NCL     : number of cloud layers
C     CLOLWC  : cloud liquid water content
C     DZ      : thickness of each atmospheric layer
C     NLAY    : number of atmospheric layers
C     LCLOUD  : logical variable to calculate with or without clouds
C     PRESS   : atmospheric pressure (mbar)
C     TEMP    : temperature (oC)
C     TW      :
C     K1      :
C     ALW     :
C     FLW     :
C     PWZW    :
C     ICLT    :
C     ICLB    :
C     TTOP    :
C     TBOT    :
C     EXPON   :
C     ESAT    :
C     ROBOT   :
C     ROTOP   :
C     CLIM    :
C     EXMASS  :
C     EWOL    :
C     FIXLWC  :
C     IWO     :
C     KWO     :
C     ZTOTAL  :
C......................................................................
C
      PARAMETER (LEVMAX = 51, LAYMAX = 50, MAXWEL= 37, MAXIND= 135)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EWOL(16,8)
      DIMENSION TEMP(LEVMAX),PRESS(LEVMAX)
      DIMENSION FLW(LAYMAX)
      DIMENSION DZ(LAYMAX)
      DIMENSION TW(LAYMAX),FIXLWC(8)
      LOGICAL*1 LCLOUD
      DATA FIXLWC /
     & 0.05, 0.14, 0.28, 0.47, 1.00, 2.50,
     & 0.13, 0.44
     & /
C
      DATA EWOL /
     &  15.2,  15.4,  15.7,  15.3,  15.6,  15.8,  16.0,  15.8,  16.0,
     &  16.0,  15.9,  16.5,  17.1,  17.2,  21.2,  18.1,
     &  39.4,  40.1,  41.0,  40.9,  42.3,  41.9,  41.6,  41.2,  42.1,
     &  42.9,  43.5,  43.2,  44.1,  45.8,  50.7,  47.1,
     &  72.4,  72.7,  74.1,  75.5,  73.8,  73.9,  74.6,  74.6,  75.0,
     &  75.2,  77.5,  78.4,  78.8,  80.0,  85.6,  82.4,
     &  73.9,  76.6,  74.4,  75.1,  75.5,  75.4,  74.8,  76.0,  76.1,
     &  76.2,  77.2,  77.4,  77.7,  80.3,  80.8,  81.2,
     & 127.8, 131.1, 130.3, 129.1, 130.9, 129.8, 130.9, 130.4, 130.6,
     & 130.7, 134.9, 132.8, 136.4, 137.4, 138.6, 140.2,
     & 121.1, 122.1, 121.2, 122.3, 121.0, 122.3, 121.0, 120.4, 121.4,
     & 122.3, 124.1, 123.3, 125.8, 126.4, 125.2, 127.1,
     &   5.5,   5.5,   5.5,   5.5,   5.5,   5.5,   5.5,   5.5,   5.6,
     &   5.6,   5.6,   5.6,   5.6,   5.6,   5.7,   5.7,
     &  12.3,  12.3,  12.3,  12.3,  12.4,  12.4,  12.4,  12.4,  12.4,
     &  12.4,  12.4,  12.5,  12.5,  12.6,  12.6,  12.7
     & /
C
      DO 10  I = 1 , NLAY
        FLW(I) = 0.0
        TW(I)  = 0.0
   10 CONTINUE
      K1 = 0
C
      IF (.NOT. LCLOUD)                      GOTO 100
C
      DO 20  I = 1 , NLAY
        PWZW =  PRESS(I)*1000.0 + 1.0
        IF (TOP .LE. PWZW)                   GOTO 200
   20 CONTINUE
  200 K1 = I-1
C
      ICLT = K1 + 1
      ICLB = ICLT + NCL
      IF(ICLB.GT.NLAY+1) THEN
        WRITE(*,*) ' cloud bottom has to be adjusted'
        ICLB = NLAY + 1
        NCL  = ICLB - ICLT
      ENDIF
C
      EXMASS  = EWOL(1,IWP)/FIXLWC(IWP)/1000.
      ZTOTAL = 0.0
C
      DO 30  IWO = 1 , NCL
        KWO      = K1 + IWO
        ZTOTAL = ZTOTAL + DZ(KWO)*1000.
   30 CONTINUE
      ALW = TAUW/EXMASS
      CLOLWC = ALW/ZTOTAL
C
      TTOP = TEMP(ICLT)
      TBOT = TEMP(ICLB)
      EXPON = (7.5D0*(TTOP-273.15D0))/(TTOP-35.85D0)
      ESAT = 6.11D0*10.0D0**EXPON
      ROTOP = 216.67D0*ESAT/TTOP
      EXPON = (7.5D0*(TBOT-273.15D0))/(TBOT-35.85D0)
      ESAT = 6.11D0*10.0D0**EXPON
      ROBOT = 216.67D0*ESAT/TBOT
      CLIM = (ROBOT-ROTOP)*0.2D0
      IF(CLIM.LT.CLOLWC) THEN
      WRITE(*,*) ' WARNING: cloud LWC unrealistically large ',TEMP
      read(*,*) 
      ENDIF
C
      DO 40  IWO = 1 , NCL
        KWO      = K1 + IWO
        FLW(KWO) = CLOLWC*DZ(KWO)*1000.
        TW(KWO)  = FLW(KWO)*EXMASS
   40 CONTINUE
C
C
  100 CONTINUE
C
      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE VARIA(RAY,NCL,NLAY,TS,ROFF,WH2O,K1,LCLOUD,DZ,TEMP,
     &                 PRESS,ROL,ROH,RO3,WH2O1,WCO2SW,WO3SW,WO2SW,RAYM)
C
C                 changes atmospheric profiles
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     TS      : difference between surface and characteristical temperatures
C     ROFF    : difference in water vapor
C     NCL     : number of cloud layers
C     WH2O    : precipitable water [cm]
C     DZ      : thickness of each atmospheric layer
C     ROL     : air density
C     ROH     : water content (humidity)
C     RO3     : ozone content 
C     NLAY    : number of atmospheric layers
C     LCLOUD  : logical variable to calculate with or without clouds
C     PRESS   : atmospheric pressure (mbar)
C     TEMP    : temperature (oC)
C     RAY     :
C     K1      :
C     PPMCO2  :
C     CONZ    :
C     NLEV    :
C     QP0     :
C     T0      :
C     DELT    :
C     DYOU    :
C     XCO2    :
C     SH2O    :
C     SH2ONR  :
C     SO3     :
C     SO3NR   :
C     SCO2    :
C     SCO2NR  :
C     KLEV    :
C     DIFT    :
C     NSTOP   :
C     TM      :
C     WAS     :
C     POO     :
C     NGRZ    :
C     KGRZ    :
C     TGRZ    :
C     EXPON   :
C     ESAT    :
C     XR      :
C     XP      :
C     XH      :
C     X3      :
C     XL      :
C     RAYM    :
C     PP      :
C     RO3M    :
C     DI      :
C     TQ      :
C     TT      :
C     WCO2SW  :
C     WCO2    :
C     WO2SW   :
C     WO3SW   :
C     WO3     :
C     WH2O1   :
C     TB      :
C......................................................................
C
      PARAMETER (LEVMAX = 51,LAYMAX = 50, MAXWEL= 37, MAXIND= 135)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL*1 LCLOUD
      DIMENSION RAY(LEVMAX)
      DIMENSION TEMP(LEVMAX),WAS(LEVMAX),PRESS(LEVMAX),
     &          POO(LEVMAX),ROL(LEVMAX),ROH(LEVMAX),RO3(LEVMAX)
      DIMENSION TM(LAYMAX),DZ(LAYMAX),DI(LAYMAX)
      DIMENSION TB(LAYMAX),RM(LAYMAX),WO3(LAYMAX),PP(LAYMAX),
     &                WCO2(LAYMAX),WH2O1(LAYMAX),TT(LAYMAX)
      DIMENSION WCO2SW(LAYMAX),WO3SW(LAYMAX),WO2SW(LAYMAX),
     &          RAYM(LAYMAX)
C
      DATA O3P,O3T / 0.36479D0, 0.0071908D0 /
C
      PPMCO2 = 300.0D0
      CONZ=PPMCO2/330.D0
      NLEV  = NLAY + 1
      QP0   = 1.D0/1.013D0
      T0    = 273.16D0
      DELT  = 0.0279D0
      DYOU  = (6.D0+3.D0*DELT)/(6.D0-7.D0*DELT)
      XCO2  = CONZ*33.D0
      SH2O   = 0.0
      SH2ONR = 0.0
      SO3    = 0.0
      SO3NR  = 0.0
      SCO2   = 0.0
      SCO2NR = 0.0
C
C    surface temperature TS and water vapor adjustment +- 0-1
C
      KLEV = 5
      DIFT = TS/(KLEV-1)
      NSTOP=NLEV - KLEV - 2
c
      DO 10  IUP = NLEV,NSTOP,-1
         KLEV = KLEV - 1
         TEMP(IUP) = TEMP(IUP) + (DIFT*KLEV)
   10 CONTINUE
c
      DO 20  I=1,NLEV
         ROH(I) = ROH(I)*(1.0+ROFF)
   20 CONTINUE
C
      DO 30  I = 1 , NLAY
         TM(I)  = (TEMP(I)+TEMP(I+1))*0.5D0
         WAS(I) = ROH(I)/ROL(I)
         POO(I) = RO3(I)/ROL(I)
   30 CONTINUE
        WAS(NLEV) = ROH(NLEV)/ROL(NLEV)
        POO(NLEV) = RO3(NLEV)/ROL(NLEV)
C
      IF (.NOT. LCLOUD)       GOTO  100
        NGRZ = NCL +1
         DO 40 I = 1 , NGRZ
           KGRZ = K1 + I
           TGRZ = TEMP(KGRZ)
           EXPON= (7.5D0*(TGRZ-273.15D0))/(TGRZ-35.85D0)
           ESAT = 6.11D0*10.0D0**EXPON
           ROH(KGRZ) = 216.67D0*ESAT/TGRZ
   40    CONTINUE
  100 CONTINUE
C
      DO 50 I=1,NLAY
      XR = DLOG(RAY(I)/RAY(I+1))
      XP = DLOG(PRESS(I)/PRESS(I+1))
      XH = DLOG(ROH(I)/ROH(I+1))
      X3 = DLOG(RO3(I)/RO3(I+1))
      XL = DLOG(ROL(I)/ROL(I+1))
C
      RAYM(I) = RAY(I)*DEXP(-0.5D0*XR) * DYOU
      PP(I)   = PRESS(I)*DEXP(-0.5D0*XP) * QP0
      RM(I)   = ROH(I)*DEXP(-0.5D0*XH)
      RO3M    = RO3(I)*DEXP(-0.5D0*X3)
      DI(I)   = ROL(I)*DEXP(-0.5D0*XL)
      TQ=TM(I)
      TT(I)=T0/TQ
      WCO2SW(I) = XCO2*PP(I)**1.75D0*TT(I)**1.35D0
      WCO2(I)   = CONZ*PP(I)**1.66D0*TT(I)
      WO2SW(I)  = PP(I)**1.75D0*TT(I)
      WO3SW(I)  = 46.667D0*RO3M
      WO3(I)    = WO3SW(I)*PP(I)**O3P*TT(I)**O3T
      WH2O1(I)   = RM(I)*0.1D0*PP(I)**0.9D0*TT(I)**0.45D0
      TB(I)     = 1745.D0 * (1.D0/TQ - 1.D0/296.D0)
C
      SH2O   = SH2O   + RM(I)*0.1D0 * DZ(I)
      SH2ONR = SH2ONR + WH2O1(I)     * DZ(I)
      SO3    = SO3    + WO3SW(I)    * DZ(I)
      SO3NR  = SO3NR  + WO3(I)      * DZ(I)
      SCO2   = SCO2   + WCO2SW(I)   * DZ(I)
      SCO2NR = SCO2NR + WCO2(I)     * DZ(I)
   50 CONTINUE
      DO 60 I=1,NLAY
      WH2O1(I) = WH2O1(I)*WH2O/SH2O
   60 CONTINUE
C
      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE AERO(SCA,NLAY,DZ,ALT,INAERO,ASIG)
C
C          determination of aerosol profiles
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     SCA     : 
C     DZ      : thickness of each atmospheric layer
C     NLAY    : number of atmospheric layers
C     ALT     : altitude (m)
C     INAERO  :
C     ASIG    :
C     K       :
C     XS      :
C     SC      :
C     EXAER   :
C     VIS1    :
C......................................................................
C
      PARAMETER (LEVMAX = 51, LAYMAX = 50, MAXWEL= 37, MAXIND= 135)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION INAERO(LAYMAX),ASIG(LAYMAX,MAXWEL)
      DIMENSION EXAER(MAXWEL,3)
      DIMENSION SCA(LEVMAX)
      DIMENSION ALT(LEVMAX)
      DIMENSION DZ(LAYMAX)
C
      DATA EXAER /
     & 1.493, 1.502, 1.514, 1.526, 1.538, 1.550, 1.550, 1.550, 1.550,
     & 1.550, 1.548, 1.534, 1.513, 1.469, 1.424, 1.360, 1.249, 1.123,
     & 1.000, 0.893, 0.789, 0.705, 0.644, 0.605, 0.587, 0.518, 0.458,
     & 0.306, 0.304, 0.245, 0.193, 0.130, 0.086, 0.051, 0.033, 0.016,
     & 0.086,
     & 2.414, 2.366, 2.302, 2.238, 2.174, 2.110, 2.054, 1.998, 1.942,
     & 1.886, 1.818, 1.733, 1.638, 1.559, 1.479, 1.378, 1.252, 1.113,
     & 1.000, 0.909, 0.824, 0.759, 0.712, 0.682, 0.668, 0.614, 0.569,
     & 0.524, 0.459, 0.411, 0.368, 0.311, 0.260, 0.194, 0.158, 0.144,
     & 0.110,
     & 1.316, 1.304, 1.288, 1.272, 1.256, 1.240, 1.226, 1.212, 1.198,
     & 1.184, 1.167, 1.148, 1.128, 1.112, 1.096, 1.076, 1.052, 1.021,
     & 1.000, 0.980, 0.961, 0.949, 0.939, 0.933, 0.930, 0.920, 0.910,
     & 0.899, 0.883, 0.866, 0.850, 0.823, 0.793, 0.748, 0.704, 0.584,
     & 0.665
     & /
C
      DO 10  I = 1 , NLAY
        INAERO(I)  = 1
        IF (ALT(I) .LT. 10.01)    INAERO(I)=2
        IF (ALT(I) .LT.  2.01)    INAERO(I)=3
   10 CONTINUE
C
         DO 20  I = 1 , 37
         DO 30  J = 1 , NLAY
            K         = INAERO(J)
            XS        = DLOG(SCA(J)/SCA(J+1))
            SC        = SCA(J)*DEXP(-0.5*XS)
            ASIG(J,I) = EXAER(I,K)*SC * DZ(J)
   30    CONTINUE
   20    CONTINUE
      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE ATMOS (NLAY,INDEX,DZ,WH2O1,WCO2SW,WO3SW,WO2SW,RAYM,U,V)
C            
C           determination of absorption and scattering
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     INDEX   : index of spectral intervals
C     DZ      : thickness of each atmospheric layer
C     NLAY    : number of atmospheric layers
C     WH2O1   :
C     WCO2SW  :
C     WO3SW   :
C     WO2SW   :
C     RAYM    :      
C     U       :
C     V       :
C     SIHO    :
C     SICO    :
C     SIO3    :
C     BC      :
C     KBC     :
C     B       :
C     RB      :
C     ANUE    :
C     RANUE   :
C     WE      :
C     XWELSW  :
C     X1      :
C     AS      :
C     X2      :
C     SR      :
C......................................................................
C
      PARAMETER (LEVMAX = 51, LAYMAX = 50, MAXWEL= 37, MAXIND= 135)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION U(LAYMAX,MAXIND),V(LAYMAX,MAXIND)
      DIMENSION RB(MAXIND),RANUE(MAXIND)
      DIMENSION KBC(MAXIND)
      DIMENSION DZ(LAYMAX),WH2O1(LAYMAX)
      DIMENSION XWELSW(MAXWEL)
      DIMENSION SICO(50),SIHO(50),SIO3(50),SR(50),
     &          INDEX(MAXIND)
      DIMENSION WCO2SW(LAYMAX),WO3SW(LAYMAX),WO2SW(LAYMAX),
     &          RAYM(LAYMAX)
C
      DATA RB /  22*0.0000E+0,
     & 0.2080E+1,0.9183E-1,0.8318E+1,0.2612E+0,0.0000E+0,0.0000E+0,
     & 0.0000E+0,0.0000E+0,0.0000E+0,0.2023E+1,0.9120E-1,0.8128E+1,
     & 0.2582E+0,0.0000E+0,0.0000E+0,0.1442E+4,0.1178E+1,0.1000E+0,
     & 0.1726E+2,0.1122E+3,0.3112E+0,0.2958E+1,0.0000E+0,0.0000E+0,
     & 0.9728E+1,0.1552E+4,0.1318E+1,0.1035E+0,0.2388E+2,0.1445E+3,
     & 0.3404E+0,0.3589E+1,0.0000E+0,0.0000E+0,0.1057E+1,0.2183E+2,
     & 0.1208E+0,0.2460E+4,0.2291E+1,0.5445E+2,0.2716E+3,0.4613E+0,
     & 0.8337E+1,0.0000E+0,0.0000E+0,0.0000E+0,0.1028E+1,0.2312E+2,
     & 0.2323E+1,0.1138E+0,0.2825E+4,0.5741E+2,0.2897E+3,0.4581E+0,
     & 0.9441E+1,0.0000E+0,0.1028E+1,0.2312E+2,0.2323E+1,0.1138E+0,
     & 0.2825E+4,0.5741E+2,0.2897E+3,0.4581E+0,0.9441E+1,0.0000E+0,
     & 0.0000E+0,0.4732E+2,0.1122E+5,0.1466E+1,0.1282E+4,0.1374E+0,
     & 0.3491E+4,0.3475E+1,0.1072E+3,0.5012E+3,0.5623E+0,0.1977E+2,
     & 0.0000E+0,0.4732E+2,0.1122E+5,0.1466E+1,0.1282E+4,0.1374E+0,
     & 0.3491E+4,0.3475E+1,0.1072E+3,0.5012E+3,0.5623E+0,0.1977E+2,
     & 0.0000E+0,0.4732E+2,0.1122E+5,0.1466E+1,0.1282E+4,0.1374E+0,
     & 0.3491E+4,0.3475E+1,0.1072E+3,0.5012E+3,0.5623E+0,0.1977E+2,
     & 0.0000E+0,0.8531E+2,0.1483E+2,0.1089E+0,0.1941E+4,0.1656E+1,
     & 0.3899E+2,0.2061E+3,0.4009E+0,0.5070E+1,0.0000E+0 /
C
      DATA RANUE /
     &   8.6500,  14.5000,  45.9000, 116.2000, 211.1000, 283.9000,
     & 282.9000, 201.5000, 106.0000,  34.6000,  10.0000,   0.8980,
     &   0.0640,   0.0018,   0.0000,   0.0000,   0.0035,   0.0345,
     &   0.0920,   0.1320,   0.0620,   0.0230,
     &   113*0.0000
     & /
C
      DATA KBC /
     &      0,     0,     0,     0,     0,     0,     0,     0,     0,
     &      0,     0,     0,     0,     0,     0,     0,     0,     0,
     &      0,     0,     0,     0,     0,     0,     0,     0,     0,
     &  20000,     0,     0,     0,     0,     0,     0,     0,     0,
     &      0,     0,     0,     0,     0,     0,     0,     0,     0,
     &      0,     0,     0,     0,     0,     0,     0,     0,     0,
     &      0,     0,    45,    45,    45,    45,    45,    45,    45,
     &     45,    45,    45,    30, 20000,   120,   120,   120,   120,
     &    120,   120,   120,   120,   120,   120, 13100, 13100, 13100,
     &  13100, 13100, 13100, 13100, 13100, 13100, 13100,     0,   230,
     &    230,   230,   230,   230,   230,   230,   230,   230,   230,
     &    230,   230,120000,120000,120000,120000,120000,120000,120000,
     & 120000,120000,120000,120000,120000,900000,900000,900000,900000,
     & 900000,900000,900000,900000,900000,900000,900000,900000,     0,
     &      0,     0,     0,     0,     0,     0,     0,     0,     0
     & /
C
      DATA XWELSW /
     & 0.200, 0.210, 0.220, 0.230, 0.240, 0.250, 0.260, 0.270, 0.280,
     & 0.290, 0.300, 0.320, 0.340, 0.360, 0.380, 0.400, 0.450, 0.500,
     & 0.550, 0.600, 0.650, 0.700, 0.735, 0.765, 0.778, 0.825, 0.877,
     & 0.950, 1.040, 1.150, 1.245, 1.430, 1.625, 1.890, 2.150, 2.590,
     & 3.280
     & /
C
      DO 10  IALL = 1 , 135
        DO 20   I = 1 , NLAY
           SIHO(I) = 0.0
           SICO(I) = 0.0
           SIO3(I) = 0.0
   20   CONTINUE
        BC   =  KBC(IALL)*0.000001
        B    =  RB(IALL)
        ANUE =  RANUE(IALL)
        WE   =  XWELSW(INDEX(IALL))
C
        IF  (IALL        .EQ. 1)             GOTO 100
        IF  (INDEX(IALL) .EQ. INDEX(IALL-1)) GOTO 200
  100   CONTINUE
        X1   =  1./(WE*WE)
        AS   =  (6432.8+2949810./(146.-X1)+25540./(41.-X1))*1.0E-08+1.
        X2   =  (AS*AS - 1.)**2.
        DO 30  I = 1 , NLAY
          SR(I)=1.352*0.1*X2*RAYM(I)/(WE**4.)
   30   CONTINUE
  200   CONTINUE
C
        IF  (KBC(IALL) .EQ. 0)               GOTO 300
        DO 40  I = 1 , NLAY
          IF (WE .GT. 1.) SICO(I)=WCO2SW(I)*BC
          IF (WE .LT. 1.) SICO(I)=WO2SW(I) *BC
   40   CONTINUE
  300   CONTINUE
C
        IF  (B .LT. 0.000001)                GOTO 400
        DO 50  I = 1 , NLAY
          SIHO(I) = WH2O1(I)*B
   50   CONTINUE
  400   CONTINUE
C
        IF  (ANUE .LT. 0.000001)             GOTO 500
        DO 60  I = 1 , NLAY
          SIO3(I) = WO3SW(I)*ANUE
   60   CONTINUE
  500   CONTINUE
C
        DO 70  I = 1 , NLAY
          V(I,IALL) = (SIO3(I) + SIHO(I) + SICO(I)) * DZ(I)
          U(I,IALL) = V(I,IALL) + SR(I) * DZ(I)
   70   CONTINUE
   10 CONTINUE
C
      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE BETAS(NCL,NLAY,IALL,COSZEN,INDEX,BET,BET0,OM,TAU,
     &                 TOHNW,TW,K1,IWP,LCLOUD,INAERO,ASIG,U,V)
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     IWP     : cloud droplet size distribution
C     NCL     : number of cloud layers
C     INDEX   : index of spectral intervals
C     NLAY    : number of atmospheric layers
C     LCLOUD  : logical variable to calculate with or without clouds
C     TW      : 
C     K1      :
C     IALL    : counter of spectral interval
C     COSZEN  : cossine of zenith angle
C     BET     :
C     BET0    :
C     OM      :
C     TAU     :
C     TOHNW   :
C     INAERO  :
C     ASIG    : 
C     U       :
C     V       :
C     UMUE    :
C     KA      :
C     KB      :
C     KC      :
C     XY      :
C     T       :
C     X1      :
C     X2      :
C     OME     :
C     GV      :
C     OMAER   :
C     GAER    :
C     EF      :
C     GE      :
C     LA      :
C     INWOSW  :
C     OWO     :
C     AWOL    :
C     GWO     :
C     GWOL    :
C     IWO     :
C     TWO     :
C     EWOL    :
C     TGES    :
C     H1      :
C     H2      :
C     H3      :
C......................................................................
C
      PARAMETER (LEVMAX = 51, LAYMAX = 50, MAXWEL= 37, MAXIND= 135)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL*1 LCLOUD
      DIMENSION U(LAYMAX,MAXIND),V(LAYMAX,MAXIND)
      DIMENSION INAERO(LAYMAX),ASIG(LAYMAX,MAXWEL)
      DIMENSION AWOL(16,8),EWOL(16,8),GWOL(16,8)
      DIMENSION OMAER(MAXWEL,3),GAER(MAXWEL,3)
      DIMENSION INWOSW(MAXIND)
      DIMENSION TW(LAYMAX)
      DIMENSION INDEX(MAXIND)
      DIMENSION BET(LAYMAX),BET0(LAYMAX),OM(LAYMAX),EF(LAYMAX)
      DIMENSION TAU(LAYMAX),TOHNW(LAYMAX)
C
      DATA AWOL /
     & .000D+0,.602D-5,.110D-4,.900D-5,.140D-4,.220D-4,.195D-3,.121D-3,
     & .300D-3,.430D-3,.501D-2,.273D-2,.317D-2,.607D-2,.398D-1,.417D+0,
     & .000D+0,.799D-5,.150D-4,.130D-4,.230D-4,.280D-4,.226D-3,.147D-3,
     & .332D-3,.513D-3,.608D-2,.332D-2,.404D-2,.793D-2,.553D-1,.438D+0,
     & .000D+0,.900D-5,.180D-4,.150D-4,.230D-4,.330D-4,.234D-3,.165D-3,
     & .400D-3,.578D-3,.777D-2,.394D-2,.473D-2,.954D-2,.688D-1,.459D+0,
     & .000D+0,.150D-4,.290D-4,.300D-4,.460D-4,.530D-4,.365D-3,.261D-3,
     & .585D-3,.909D-3,.113D-1,.607D-2,.728D-2,.157D-1,.114D+0,.481D+0,
     & .000D+0,.170D-4,.310D-4,.270D-4,.410D-4,.610D-4,.450D-3,.319D-3,
     & .800D-3,.116D-2,.135D-1,.740D-2,.883D-2,.176D-1,.130D+0,.474D+0,
     & .101D-5,.280D-4,.470D-4,.450D-4,.520D-4,.106D-3,.779D-3,.638D-3,
     & .150D-2,.232D-2,.296D-1,.175D-1,.215D-1,.424D-1,.255D+0,.468D+0,
     & .256D-5,.269D-4,.349D-4,.574D-4,.785D-4,.168D-3,.278D-3,.732D-3,
     & .886D-3,.250D-2,.614D-2,.799D-1,.615D-1,.768D-1,.124D+0,.456D+0,
     & .373D-5,.393D-4,.519D-4,.852D-4,.115D-3,.296D-3,.415D-3,.107D-2,
     & .132D-2,.366D-2,.905D-2,.113D+0,.875D-1,.108D+0,.170D+0,.452D+0
     & /
C
      DATA GWOL /
     &  .832,  .830,  .811,  .829,  .819,  .812,  .798,  .803,  .799,
     &  .795,  .776,  .771,  .767,  .796,  .881,  .820,
     &  .836,  .828,  .813,  .830,  .820,  .818,  .797,  .805,  .803,
     &  .801,  .789,  .790,  .773,  .798,  .868,  .829,
     &  .848,  .819,  .840,  .836,  .838,  .830,  .817,  .829,  .820,
     &  .812,  .791,  .799,  .776,  .788,  .860,  .852,
     &  .862,  .857,  .825,  .850,  .849,  .851,  .849,  .848,  .847,
     &  .846,  .836,  .838,  .829,  .822,  .845,  .923,
     &  .855,  .834,  .843,  .839,  .838,  .837,  .833,  .837,  .830,
     &  .826,  .814,  .818,  .808,  .814,  .857,  .887,
     &  .860,  .849,  .850,  .847,  .847,  .857,  .843,  .828,  .835,
     &  .842,  .834,  .834,  .824,  .818,  .872,  .904,
     &  .885,  .885,  .885,  .885,  .885,  .885,  .884,  .884,  .884,
     &  .884,  .884,  .897,  .895,  .900,  .923,  .938,
     &  .887,  .888,  .888,  .888,  .888,  .888,  .888,  .888,  .887,
     &  .888,  .889,  .906,  .904,  .910,  .935,  .939
     & /
C
      DATA EWOL /
     &  15.2,  15.4,  15.7,  15.3,  15.6,  15.8,  16.0,  15.8,  16.0,
     &  16.0,  15.9,  16.5,  17.1,  17.2,  21.2,  18.1,
     &  39.4,  40.1,  41.0,  40.9,  42.3,  41.9,  41.6,  41.2,  42.1,
     &  42.9,  43.5,  43.2,  44.1,  45.8,  50.7,  47.1,
     &  72.4,  72.7,  74.1,  75.5,  73.8,  73.9,  74.6,  74.6,  75.0,
     &  75.2,  77.5,  78.4,  78.8,  80.0,  85.6,  82.4,
     &  73.9,  76.6,  74.4,  75.1,  75.5,  75.4,  74.8,  76.0,  76.1,
     &  76.2,  77.2,  77.4,  77.7,  80.3,  80.8,  81.2,
     & 127.8, 131.1, 130.3, 129.1, 130.9, 129.8, 130.9, 130.4, 130.6,
     & 130.7, 134.9, 132.8, 136.4, 137.4, 138.6, 140.2,
     & 121.1, 122.1, 121.2, 122.3, 121.0, 122.3, 121.0, 120.4, 121.4,
     & 122.3, 124.1, 123.3, 125.8, 126.4, 125.2, 127.1,
     &   5.5,   5.5,   5.5,   5.5,   5.5,   5.5,   5.5,   5.5,   5.6,
     &   5.6,   5.6,   5.6,   5.6,   5.6,   5.7,   5.7,
     &  12.3,  12.3,  12.3,  12.3,  12.4,  12.4,  12.4,  12.4,  12.4,
     &  12.4,  12.4,  12.5,  12.5,  12.6,  12.6,  12.7
     & /
C
      DATA OMAER /
     & 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     & 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     & 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     & 1.000, 1.000, 1.000, 1.000, 0.999, 0.995, 0.980, 0.952, 0.775,
     & 0.072,
     & 0.661, 0.679, 0.703, 0.727, 0.751, 0.775, 0.793, 0.811, 0.828,
     & 0.846, 0.867, 0.885, 0.902, 0.902, 0.901, 0.901, 0.899, 0.898,
     & 0.891, 0.889, 0.885, 0.879, 0.869, 0.863, 0.860, 0.849, 0.839,
     & 0.826, 0.807, 0.793, 0.781, 0.761, 0.751, 0.767, 0.769, 0.652,
     & 0.758,
     & 0.848, 0.859, 0.874, 0.890, 0.905, 0.920, 0.931, 0.942, 0.954,
     & 0.965, 0.977, 0.980, 0.984, 0.985, 0.985, 0.986, 0.987, 0.988,
     & 0.988, 0.989, 0.990, 0.990, 0.989, 0.989, 0.989, 0.988, 0.988,
     & 0.987, 0.986, 0.986, 0.986, 0.985, 0.986, 0.985, 0.983, 0.901,
     & 0.764
     & /
C
      DATA GAER /
     & 0.681, 0.683, 0.685, 0.688, 0.690, 0.693, 0.697, 0.701, 0.705,
     & 0.709, 0.714, 0.719, 0.726, 0.730, 0.735, 0.739, 0.737, 0.734,
     & 0.726, 0.716, 0.705, 0.694, 0.682, 0.675, 0.672, 0.659, 0.647,
     & 0.626, 0.597, 0.566, 0.535, 0.480, 0.422, 0.350, 0.289, 0.208,
     & 0.137,
     & 0.724, 0.718, 0.711, 0.703, 0.696, 0.688, 0.682, 0.676, 0.671,
     & 0.665, 0.658, 0.655, 0.651, 0.649, 0.648, 0.645, 0.643, 0.639,
     & 0.637, 0.635, 0.632, 0.631, 0.632, 0.632, 0.632, 0.633, 0.633,
     & 0.632, 0.631, 0.633, 0.636, 0.641, 0.656, 0.697, 0.732, 0.776,
     & 0.782,
     & 0.774, 0.772, 0.769, 0.765, 0.762, 0.759, 0.757, 0.754, 0.752,
     & 0.749, 0.747, 0.745, 0.744, 0.745, 0.745, 0.746, 0.747, 0.747,
     & 0.744, 0.747, 0.750, 0.751, 0.752, 0.753, 0.753, 0.754, 0.755,
     & 0.757, 0.760, 0.762, 0.764, 0.769, 0.775, 0.783, 0.790, 0.835,
     & 0.760
     & /
C     
      DATA INWOSW /
     & 22* 1,  5* 2,  3* 3,  4,  5* 5,  6,  8* 7,  8, 9* 9,  10,
     & 10*11,  2*12, 20*13, 14, 36*15, 10*16
     & /
C
      KA    =  INDEX(IALL)
      DO 10  I = 1 , NLAY
         KB       = INAERO(I)
         XY       = ASIG(I,KA)
         T        = XY + U(I,IALL)
         X1       = XY*OMAER(KA,KB)
         X2       = U(I,IALL) - V(I,IALL)
         OME      = (X1+X2)/T
         GV       = X1*GAER(KA,KB)/(X1+X2)
         EF(I)    = GV*GV
         OM(I)    = OME*(1.-EF(I))/(1.-OME*EF(I))
         TAU(I)   = T*(1.-OME*EF(I))
         TOHNW(I) = TAU(I)
         GE       = (GV-EF(I))/(1.-EF(I))
         BET(I)   = (3.-3.*GE)/8.
         BET0(I)  = 0.5-COSZEN*0.75*GE
   10 CONTINUE
C
      IF (.NOT. LCLOUD)                      GOTO 100
C
       LA       = INWOSW(IALL)
       OWO      = 1. - AWOL(LA,IWP)
       GWO      = GWOL(LA,IWP)
C
      DO 20  IWO = 1 , NCL
         KC       = K1 + IWO
         TWO      = TW(KC)*EWOL(LA,IWP)/EWOL(1,IWP)
         TGES     = ASIG(KC,KA) + U(KC,IALL) + TWO
         H1       = U(KC,IALL) - V(KC,IALL)
         H2       = OWO*TWO
         KB       = INAERO(KC)
         H3       = OMAER(KA,KB)*ASIG(KC,KA)
         OME      = (H1+H2+H3)/TGES
         GV       = (H3*GAER(KA,KB)+GWO*H2)/(H1+H2+H3)
         EF(KC)   = GV*GV
         OM(KC)   = OME*(1.-EF(KC))/(1.-OME*EF(KC))
         GE       = (GV-EF(KC))/(1.-EF(KC))
         TAU(KC)  = TGES*(1.-OME*EF(KC))
         BET(KC)  = (3.-3.*GE)/8.
         BET0(KC) = 0.5-COSZEN*0.75*GE
   20 CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE MATBAU(NLY,ALDIF,ALDR,COSZEN,BET,BET0,OM,TAU,
     &                  DIR,DIF)
C
C      calculation of diffuse irradiation for each spectral interval
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     DIR     : direct irradiation for each spectral interval (Wh/m2)
C     DIF     : diffuse irradiation for each spectral interval (Wh/m2)
C     NLY     : NLAY-1
C     ALDIF   : surface albedo for cloudy sky
C     ALDR    : surface albedo for clear sky
C     COSZEN  : cossine of zenith angle
C     BET     :
C     BET0    :
C     OM      :
C     TAU     :
C     UMUE    :
C     ISS     :
C     WU      :
C     BT      :
C     BT0     :
C     OMEG    :
C     HH      :
C     HMH     :
C     A1      :
C     A2      :
C     A3      :
C     A4      :
C     E       :
C     D1      :
C     D2      :
C     CC      :
C     X       :
C     Y       :
C     D       :
C     F       :
C     DW      :
C     XH      :
C     G       :
C     N       :
C     AM3     :
C     AM4     :
C     AM1     :
C     AM2     :
C     S       :
C     UP      :
C     HHH     :
C     HHD     :
C     XNET    :
C     DX1     :
C     DX      :
C......................................................................
C
      PARAMETER (LEVMAX = 51, LAYMAX = 50, MAXWEL= 37, MAXIND= 135)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION AM1(100),AM2(100),AM3(100),AM4(100),S(100)
      DIMENSION DIR(LEVMAX),DIF(LEVMAX),UP(LEVMAX)
      DIMENSION D(50),F(50),G(50),XH(50),X(50),Y(50)
      DIMENSION BET(LAYMAX),BET0(LAYMAX),OM(LAYMAX),TAU(LAYMAX)
C
      UMUE  =  2.
      ISS  =   NLY + 1
      WU   =   1. / COSZEN
C
      DO 10  I = 1 , NLY
        BT   = BET (I)
        BT0  = BET0(I)
        OMEG = OM  (I)
C
        HH   = DBLE(0.9999999)
        OMEG = DMIN1(OMEG,HH)
        HMH  = DBLE(0.0000001)
        OMEG = DMAX1(OMEG,HMH)
C
        A1   = (1.-OMEG*(1.-BT)) * UMUE
        A2   = BT   * OMEG * UMUE
        A3   = OMEG * BT0
        A4   = OMEG * (1.-BT0)
        E    = SQRT(A1*A1-A2*A2)
        D1   = -A1*A3 - A2*A4 + A3/COSZEN
        D2   = -A1*A4 - A2*A3 - A4/COSZEN
        CC   = WU*WU - E*E
        X(I) = D1/CC
        Y(I) = D2/CC
        D(I) = A2/(A1+E)
        F(I) = A2/(A1-E)
        DW   = E * TAU(I)
        DW   = DMIN1(100.0D0,DW)
        XH(I)= DEXP(-DW)
        G(I) = 1./XH(I)
   10 CONTINUE
C
      N      =   2*NLY
      AM3(1) =   D(1)
      AM4(1) =   F(1)
      AM1(N) =   (1.-ALDIF*D(NLY)) * G (NLY)
      AM2(N) =   (1.-ALDIF*F(NLY)) * XH(NLY)
      S(N)   =   ((Y(NLY)*ALDIF-X(NLY)) * WU  +  ALDR) * DIR(ISS)
      S(1)   =  -Y(1)
C
      DO 20  L = 2 , NLY
        LA      =   2*L-1
        LB      =   2*L-2
        AM1(LB) =   D(L-1)*G(L-1)
        AM2(LB) =   F(L-1)*XH(L-1)
        AM3(LB) =  -D(L)
        AM4(LB) =  -F(L)
        AM1(LA) =   G(L-1)
        AM2(LA) =   XH(L-1)
        S(LB)   =   DIR(L)*WU*(Y(L)-Y(L-1))
        S(LA)   =   DIR(L)*WU*(X(L)-X(L-1))
   20 CONTINUE                             
C
C......................................................................
C
C     subroutine LINGAU - 
C     input  - N,AM1,AM2,AM3,AM4,S
C     output - AM2,AM3,AM4,S
C
      CALL LINGAU(N,AM1,AM2,AM3,AM4,S)
C
C......................................................................
C
      UP (ISS) =  S(2*NLY-1)*G(NLY)+S(2*NLY)*XH(NLY)+X(NLY)*WU*DIR(ISS)
      DIF(ISS) =  S(2*NLY-1)*G(NLY)*D(NLY)+S(2*NLY)*XH(NLY)*F(NLY)+
     &            Y(NLY)*WU*DIR(ISS)
C
      DO 30  I = 1 , NLY
        DIF(I) = S(2*I-1)*D(I)+S(2*I)*F(I)+Y(I)*WU*DIR(I)
        UP (I) = S(2*I-1)+S(2*I)+X(I)* WU*DIR(I)
   30 CONTINUE
C
      HHH  =   DBLE(0.0)
      HHD  =   DIR(1)
      DO 40  I = 2 , ISS
        DIF(I) = DMAX1(DIF(I),HHH)
        XNET = DIF(I) - UP(I)
        DX1  = DIF(I-1) + DIR(I-1)
        DX   = DIF(I)   + DIR(I)
         IF(DX.GT.DX1) THEN
           DIF(I) = DX1 - DIR(I)
           UP(I)  = DIF(I) - XNET
         ENDIF
        IF (UP(I) .GT. HHD)     UP(I) = HHH
        UP(I) = DMAX1(UP(I),HHH)
   40 CONTINUE
C
      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE LINGAU(N,AM1,AM2,AM3,AM4,S)
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     N       :
C     AM3     :
C     AM4     :
C     AM1     :
C     AM2     :
C     S       :
C     N3      :
C     L       :
C     L1      :
C     L2      :
C     A       :
C     N1      :
C     I       :
C     II      :
C     I1      :
C     I2      :
C......................................................................
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION AM1(100),AM2(100),AM3(100),AM4(100),S(100)
C
      N3   =  N - 3
C
      DO 10  L  =  1 , N3 , 2
         L1      =  L  + 1
         L2      =  L1 + 1
         A       = -AM1(L1) / AM3(L)
         S(L1)   =  A * S(L) + S(L1)
         AM2(L1) =  AM2(L1) + A * AM4(L)
         A       = -AM1(L2) / AM3(L)
         S(L2)   =  A * S(L) + S(L2)
         AM2(L2) =  AM2(L2) + A * AM4(L)
         A       = -AM2(L2) / AM2(L1)
         S(L2)   =  S(L2) + A * S(L1)
         AM3(L2) =  A * AM3(L1) - 1.
         AM4(L2) =  A * AM4(L1) - 1.
   10 CONTINUE
C
      N1   =  N - 1
      A    = -AM1(N) / AM3(N1)
      S(N) =  (S(N) + A * S(N1))/ (AM2(N) + A * AM4(N1))
C
      DO 20  II =  1 , N3 , 2
         I       =  N1 - II
         I1      =  I  + 1
         I2      =  I1 + 1
         S(I1)   =  (S(I1) - AM4(I1) * S(I2)) / AM3(I1)
         S(I)    =  (S(I) - AM3(I) * S(I1) - AM4(I) * S(I2)) / AM2(I)
   20 CONTINUE
C
      S(1) = (S(1) - AM4(1) * S(2)) / AM3(1)
C
      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DOUBLE PRECISION FUNCTION PTZQ1(DZ,LATMOS)
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DZB(10,6),C(11,6),DELTA(10,6),TSTAR(6)
C
C----->>> tropical sounding <<<
C
      DATA (DZB(N,1),N=1,10)/2.0D0,3.0D0,16.5D0,21.5D0,45.0D0,51.0D0,
     & 70.0D0,100.0D0,200.0D0,300.0D0/
      DATA (C(N,1),N=1,11)/-6.0D0,-4.0D0,-6.7D0,4.0D0,2.2D0,1.0D0,
     & -2.8D0,-0.27D0,0.0D0,0.0D0,0.0D0/
      DATA (DELTA(N,1),N=1,10)/0.5D0,0.5D0,0.3D0,0.5D0,6*1.0D0/
C
C----->>> midlatitude summer <<<
C
      DATA (DZB(N,2),N=1,10)/1.5D0,6.5D0,13.0D0,18.0D0,26.0D0,36.0D0,
     & 48.0D0,50.0D0,70.0D0,100.0D0/
      DATA (C(N,2),N=1,11)/-4.0D0,-6.0D0,-6.5D0,0.0D0,1.2D0,2.2D0,
     & 2.5D0,0.0D0,-3.0D0,-0.25D0,0.0D0/
      DATA(DELTA(N,2),N=1,10)/0.5D0,1.0D0,0.5D0,0.5D0,1.0D0,1.0D0,
     & 2.5D0,0.5D0,1.0D0,1.0D0/
C
C----->>> midlatitude winter <<<
C
      DATA (DZB(N,3),N=1,10)/3.0D0,10.0D0,19.0D0,25.0D0,32.0D0,44.5D0,
     & 50.0D0,71.0D0,98.0D0,200.0D0/
      DATA (C(N,3),N=1,11)/-3.5D0,-6.0D0,-0.5D0,0.0D0,0.4D0,3.2D0,
     & 1.6D0,-1.8D0,-0.7D0,0.0D0,0.0D0/
      DATA (DELTA(N,3),N=1,10)/0.5D0,0.5D0,8*1.0D0/
C
C----->>> subarctic summer <<<
C
      DATA (DZB(N,4),N=1,10)/4.7D0,10.0D0,23.0D0,31.8D0,44.0D0,50.2D0,
     & 69.2D0,100.0D0,102.0D0,103.0D0/
      DATA (C(N,4),N=1,11)/-5.3D0,-7.0D0,0.0D0,1.4D0,3.0D0,0.7D0,
     & -3.3D0,-0.2D0,3*0.0D0/
      DATA(DELTA(N,4),N=1,10)/0.5D0,0.3D0,1.0D0,1.0D0,2.0D0,1.0D0,
     & 1.5D0,3*1.0D0/
C
C----->>> subarctic winter <<<
C
      DATA (DZB(N,5),N=1,10)/1.0D0,3.2D0,8.5D0,15.5D0,25.0D0,30.0D0,
     & 35.0D0,50.0D0,70.0D0,100.0D0/
      DATA (C(N,5),N=1,11)/3.0D0,-3.2D0,-6.8D0,0.0D0,-0.6D0,1.0D0,
     & 1.2D0,2.5D0,-0.7D0,-1.2D0,0.0D0/
      DATA (DELTA(N,5),N=1,10)/0.4D0,1.5D0,0.3D0,0.5D0,6*1.0D0/
C
C----->>> us standard 1976 <<<
C
      DATA (DZB(N,6),N=1,10)/11.0D0,20.1D0,32.1D0,47.4D0,51.4D0,
     & 71.7D0,85.7D0,90.0D0,91.0D0,92.0D0/
      DATA (C(N,6),N=1,11)/-6.5D0,0.0D0,1.0D0,2.75D0,0.0D0,-2.75D0,
     & -1.97D0,4*0.0D0/
      DATA (DELTA(N,6),N=1,10)/0.3D0,9*1.0D0/
      DATA (TSTAR(I),I=1,6) /
     &  300.0D0,294.0D0,272.2D0,287.0D0,257.1D0,288.15D0/
C     
      NLAST=10
      DTEMP=TSTAR(LATMOS)+C(1,LATMOS)*DZ
C
      DO 10 N=1,NLAST
      DEXPO=(DZ-DZB(N,LATMOS))/DELTA(N,LATMOS)
      DEXPP=DZB(N,LATMOS)/DELTA(N,LATMOS)
      DFAC=DEXP(DEXPP)+DEXP(-DEXPP)
C
      IF((DABS(DEXPO)-100.0D0).LE.0.0) THEN
      DX=DEXP(DEXPO)
      DY=DX+1.0D0/DX
      DZLOG=DLOG(DY)
      ELSE
      DZLOG=DABS(DEXPO)
      END IF
C
      IF((DEXPP-100.D0).LT.0.0) THEN
      DFACLG=DLOG(DFAC)
      ELSE
      DFACLG=DEXPP
      END IF
C
      DTEMP=DTEMP+(C(N+1,LATMOS)-C(N,LATMOS))*0.5D0*
     &      (DZ+DELTA(N,LATMOS)*(DZLOG-DFACLG))
   10 CONTINUE
C
      PTZQ1=DTEMP
      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DOUBLE PRECISION FUNCTION PTZQ2(P,LATMOS,QL)
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(6,6),B(6,6),QS(6),PRAT(6,6),PS(6)
      DIMENSION PMIN(6),QCONST(6)
      DATA (QS(I),I=1,6) /
     &  19.D0,14.D0,3.5D0,9.1D0,1.2D0,14.D0/
      DATA ((A(I,J),I=1,6),J=1,6) /
     &  3.55D0,1.37D0,3.99D0,2.79D0,8.66D0,4.66D0,
     &  3.31D0,4.77D0,4.19D0,1.30D0,8.21D0,4.52D0,
     &  2.72D0,2.11D0,3.41D0,5.44D0,3.91D0,0.00D0,
     &  3.76D0,2.95D0,0.83D0,4.71D0,7.43D0,4.28D0,
     &  0.00D0,-0.59D0,4.99D0,4.18D0,3.30D0,0.00D0,
     &  3.31D0,4.77D0,4.19D0,1.30D0,8.21D0,4.52D0/
      DATA ((B(I,J),I=1,6),J=1,6) /
     &  1.91D0,-7.58D0,-0.41D0,-1.52D0,1.81D0,0.00D0,
     &  -1.83D0,0.67D0,0.00D0,-2.26D0,1.67D0,0.00D0,
     &  0.39D0,-2.02D0,-0.81D0,0.71D0,0.00D0,0.00D0,
     &  2.37D0,-1.01D0,-3.17D0,0.00D0,1.41D0,0.00D0,
     &  0.00D0,-4.46D0,0.86D0,0.38D0,0.00D0,0.00D0,
     &  -1.83D0,0.67D0,0.00D0,-2.26D0,1.67D0,0.00D0/
      DATA ((PRAT(I,J),I=1,6),J=1,6) /
     &  0.795D0,0.694D0,0.339D0,0.172D0,0.110D0,0.0D0,
     &  0.558D0,0.421D0,0.278D0,0.172D0,0.110D0,0.0D0,
     &  0.776D0,0.342D0,0.263D0,0.100D0,0.0D0,0.0D0,
     &  0.787D0,0.375D0,0.294D0,0.145D0,0.107D0,0.0D0,
     &  0.876D0,0.350D0,0.185D0,0.100D0,0.00D0,0.00D0,
     &  0.558D0,0.421D0,0.278D0,0.172D0,0.110D0,0.0D0/
      DATA (PS(I),I=1,6) /
     &  1013.D0,1013.D0,1018.D0,1010.D0,1013.D0,1013.D0/
      DATA (PMIN(I),I=1,6) /
     &  110.99D0,152.99D0,117.79D0,124.99D0,110.29D0,152.99D0/
      DATA (QCONST(I),I=1,6) /
     &  3.25D-06,5*4.0D-06/
      IF(P.GT.PS(LATMOS)) P=PS(LATMOS)
      PR=P/PS(LATMOS)
      N=1
  300 IF(PR.GE.PRAT(N,LATMOS)) GOTO 100
      IF(P.LE.PMIN(LATMOS)) GOTO 200
      N=N+1
      GOTO 300
  100 ETA=A(N,LATMOS)+B(N,LATMOS)*DLOG(PR)
      PTZQ2=QS(LATMOS)*PR**ETA
      RETURN
  200 PTZQ2=QCONST(LATMOS)*QL
      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION PTZQ3(PPP,LATMOS)
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P(33,6),O(33,6)
      DATA (P(I,1),I=1,33)/
     &       1.013D+05,9.04D+04,8.05D+04,7.15D+04,6.33D+04,5.59D+04,
     &       4.92D+04,4.32D+04,3.78D+04,3.29D+04,2.86D+04,2.47D+04,
     &       2.13D+04,1.82D+04,1.56D+04,1.32D+04,1.11D+04,9.37D+03,
     &       7.89D+03,6.66D+03,5.65D+03,4.80D+03,4.09D+03,3.50D+03,
     &       3.00D+03,2.57D+03,1.22D+03,6.00D+02,3.05D+02,1.59D+02,
     &       8.54D+01,5.79D+00,3.00D-02/
      DATA (P(I,2),I=1,33)/
     &       1.013D+05,9.02D+04,8.02D+04,7.10D+04,6.28D+04,5.54D+04,
     &       4.87D+04,4.26D+04,3.72D+04,3.24D+04,2.81D+04,2.43D+04,
     &       2.09D+04,1.79D+04,1.53D+04,1.30D+04,1.11D+04,9.50D+03,
     &       8.12D+03,6.95D+03,5.95D+03,5.10D+03,4.37D+03,3.76D+03,
     &       3.22D+03,2.77D+03,1.32D+03,6.52D+02,3.33D+02,1.76D+02,
     &       9.51D+01,6.71D+00,3.00D-02/
      DATA (P(I,3),I=1,33)/
     &      1.018D+05,8.973D+04,7.897D+04,6.938D+04,6.081D+04,5.313D+04
     &     ,4.627D+04,4.016D+04,3.473D+04,2.992D+04,2.568D+04,2.199D+04
     &     ,1.882D+04,1.610D+04,1.378D+04,1.178D+04,1.007D+04,8.610D+03
     &     ,7.350D+03,6.280D+03,5.370D+03,4.580D+03,3.910D+03,3.340D+03
     &     ,2.860D+03,2.430D+03,1.110D+03,5.180D+02,2.530D+02,1.290D+02
     &     ,6.820D+01,4.670D+00,3.000D-02/
      DATA (P(I,4),I=1,33)/
     &      1.010D+05,8.960D+04,7.929D+04,7.000D+04,6.160D+04,5.410D+04
     &     ,4.730D+04,4.130D+04,3.590D+04,3.107D+04,2.677D+04,2.300D+04
     &     ,1.977D+04,1.700D+04,1.460D+04,1.250D+04,1.080D+04,9.280D+03
     &     ,7.980D+03,6.860D+03,5.890D+03,5.070D+03,4.360D+03,3.750D+03
     &     ,3.227D+03,2.780D+03,1.340D+03,6.610D+02,3.400D+02,1.810D+02
     &     ,9.870D+01,7.070D+00,3.000D-02/
      DATA (P(I,5),I=1,33)/
     &      1.013D+05,8.878D+04,7.775D+04,6.798D+04,5.932D+04,5.158D+04
     &     ,4.467D+04,3.853D+04,3.308D+04,2.829D+04,2.418D+04,2.067D+04
     &     ,1.766D+04,1.510D+04,1.291D+04,1.103D+04,9.431D+03,8.058D+03
     &     ,6.882D+03,5.875D+03,5.014D+03,4.277D+03,3.647D+03,3.109D+03
     &     ,2.649D+03,2.256D+03,1.020D+03,4.701D+02,2.243D+02,1.113D+02
     &     ,5.719D+01,4.016D+00,3.000D-02/
      DATA (P(I,6),I=1,33)/
     &      1.013D+05,8.988D+04,7.950D+04,7.012D+04,6.166D+04,5.405D+04
     &     ,4.722D+04,4.111D+04,3.565D+04,3.080D+04,2.650D+04,2.270D+04
     &     ,1.940D+04,1.658D+04,1.417D+04,1.211D+04,1.035D+04,8.850D+03
     &     ,7.565D+03,6.467D+03,5.529D+03,4.729D+03,4.048D+03,3.467D+03
     &     ,2.972D+03,2.549D+03,1.197D+03,5.746D+02,2.871D+02,1.491D+02
     &     ,7.978D+01,5.220D+00,3.008D-02/
      DATA (O(I,1),I=1,33)/
     & 5.6D-08,5.6D-08,5.4D-08,5.1D-08,4.7D-08,4.5D-08,4.3D-08,4.1D-08,
     & 3.9D-08,3.9D-08,3.9D-08,4.1D-08,4.3D-08,4.5D-08,4.5D-08,4.7D-08,
     & 4.7D-08,6.9D-08,9.0D-08,1.4D-07,1.9D-07,2.4D-07,2.8D-07,3.2D-07,
     & 3.4D-07,3.4D-07,2.4D-07,9.2D-08,4.1D-08,1.3D-08,4.3D-09,8.6D-11,
     & 4.3D-14/
      DATA (O(I,2),I=1,33)/
     & 6.0D-08,6.0D-08,6.0D-08,6.2D-08,6.4D-08,6.6D-08,6.9D-08,7.5D-08,
     & 7.9D-08,8.6D-08,9.0D-08,1.1D-07,1.2D-07,1.5D-07,1.8D-07,1.9D-07,
     & 2.1D-07,2.4D-07,2.8D-07,3.2D-07,3.4D-07,3.6D-07,3.6D-07,3.4D-07,
     & 3.2D-07,3.0D-07,2.0D-07,9.2D-08,4.1D-08,1.3D-08,4.3D-09,8.6D-11,
     & 4.3D-14/
      DATA (O(I,3),I=1,33)/
     & 6.0D-08,5.4D-08,4.9D-08,4.9D-08,4.9D-08,5.8D-08,6.4D-08,7.7D-08,
     & 9.0D-08,1.2D-07,1.6D-07,2.1D-07,2.6D-07,3.0D-07,3.2D-07,3.4D-07,
     & 3.6D-07,3.9D-07,4.1D-07,4.3D-07,4.5D-07,4.3D-07,4.3D-07,3.9D-07,
     & 3.6D-07,3.4D-07,1.9D-07,9.2D-08,4.1D-08,1.3D-08,4.3D-09,8.6D-11,
     & 4.3D-14/
      DATA (O(I,4),I=1,33)/
     & 4.9D-08,5.4D-08,5.6D-08,5.8D-08,6.0D-08,6.4D-08,7.1D-08,7.5D-08,
     & 7.9D-08,1.1D-07,1.3D-07,1.8D-07,2.1D-07,2.6D-07,2.8D-07,3.2D-07,
     & 3.4D-07,3.9D-07,4.1D-07,4.1D-07,3.9D-07,3.6D-07,3.2D-07,3.0D-07,
     & 2.8D-07,2.6D-07,1.4D-07,9.2D-08,4.1D-08,1.3D-08,4.3D-09,8.6D-11,
     & 4.3D-14/
      DATA (O(I,5),I=1,33)/
     & 4.1D-08,4.1D-08,4.1D-08,4.3D-08,4.5D-08,4.7D-08,4.9D-08,7.1D-08,
     & 9.0D-08,1.6D-07,2.4D-07,3.2D-07,4.3D-07,4.7D-07,4.9D-07,5.6D-07,
     & 6.2D-07,6.2D-07,6.2D-07,6.0D-07,5.6D-07,5.1D-07,4.7D-07,4.3D-07,
     & 3.6D-07,3.2D-07,1.5D-07,9.2D-08,4.1D-08,1.3D-08,4.3D-09,8.6D-11,
     & 4.3D-14/
      DATA (O(I,6),I=1,33)/
     & 5.4D-08,5.4D-08,5.4D-08,5.0D-08,4.6D-08,4.6D-08,4.5D-08,4.9D-08,
     & 5.2D-08,7.1D-08,9.0D-08,1.3D-07,1.6D-07,1.7D-07,1.9D-07,2.1D-07,
     & 2.4D-07,2.8D-07,3.2D-07,3.5D-07,3.8D-07,3.8D-07,3.9D-07,3.8D-07,
     & 3.6D-07,3.4D-07,2.0D-07,1.1D-07,4.9D-08,1.7D-08,5.3D-09,4.3D-11,
     & 4.3D-14/
      PP=PPP
      PP=PP*100.0D0
      N=2
  200 IF(PP.GT.P(N,LATMOS)) GOTO 100
      N=N+1
      GOTO 200
  100 IF(PP.EQ.P(N-1,LATMOS)) GOTO 300
      B=(O(N-1,LATMOS)-O(N,LATMOS))/(P(N-1,LATMOS)-P(N,LATMOS))
      A=O(N,LATMOS)-B*P(N,LATMOS)
      PTZQ3=A+B*PP
      RETURN
  300 PTZQ3=O(N-1,LATMOS)
      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
