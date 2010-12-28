C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE D2STR(XLON,YLAT,TSTAT,rf,zstat,SFCALB,VIS,JDAY,XHOR,
     &                 ICLOUD,TAUW,TRANS,TDIR)
C
C     calculate transmittance for clear and cloudy sky
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     XLON    : longitude of the ground station (o)
C     YLAT    : latitude of the ground station (o)
C     TSTAT   : temperature at ground station (K)
C     WH2O    : precipitable water (cm)
C     SFCALB  : surface albedo (0/1)
C     VIS     : visibility at ground station (km)
C     JDAY    : julian day
C     PI      : pi
C     CDR     : change degree to radians
C     SC0     : solar constant (W/m2)
C     E0      : eccentricity correction factor of the earth's orbit
C     DEC     : solar declination (degrees)                            
C     DECR    : solar declination (radians) 
C     YLATR   : latitude of the ground station (radians)
C     CODEC   : cossine of solar declination  
C     SIDEC   : sine of solar declination 
C     SILAT   : sine of latitude of ground station
C     ET      : equation of time (minutes)
C     ZN      : time zone number (3.0 - Florianopolis) (hours)
C     TIMCOR  : time correction (hours)       
C     T1S     : standard time at begin of integration interval (hours)
C     T2S     : standard time at end of integration interval (hours)
C     T1A     : apparent time at begin of integration interval (hours)
C     T2A     : apparent time at end of integration interval (hours)
C     IT      : hour of sattelite image
C     COWI    : cossine of hour angle at middle of integration interval
C     COSZEN  : cossine of zenith angle
C     THETA   : zenith angle (degrees)
C     ROFF    : difference in water vapor
C     IWP     : cloud droplet size distribution
C     ISUB    : use subroutine WOLKE1 (ISUB=1) or WOLKE2 (ISUB=2)
C     TOP     : cloud top (mbar)
C     NCL     : number of cloud layers
C     INTVAL  : spectral interval
C     CLOLWC  : cloud liquid water content
C     W1      : hour angle at begin of interval (degrees)
C     W2      : hour angle at end of interval (degrees)
C     WI      : characteristic hour angle (degrees)
C     LATMOS  : atmosphere type
C               1 - tropical
C               2 - midlatitude summer 
C               3 - midlatitude winter
C               4 - subartic summer 
C               5 - subartic winter
C     TS      : difference between characteristical temperature of LATMOS
C              and temperature at gound station (oC)
C     ICLOUD  : cloud parameter 
C              0 - no clouds
C              1 - calculation with clouds
C     TAUW    : cloud optical thickness
C              0.0 - no clouds 
C              100.0 - with clouds 
C     TRANS   : parameterized sky transmittance for broadband solar radiation 
C              in the atmosphere 
C......................................................................
C                                                                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NHOUR=24)
C
C     definition of constants bellow
C
      PI    = ACOS(-1.0)
      CDR  = PI/180.0
C
C......................................................................
C     subroutine ASTRO - calculation of eccentricity correction,
C     declination an equation of time
C     input  - JDAY
C     output - E0,DEC,ET
C
      CALL ASTRO(JDAY,E0,DEC,ET)
C
C......................................................................
C
C     change degrees to radians and calculate sines and cossines
C
      
      DECR = DEC*CDR
      YLATR= YLAT*CDR
      CODEC = COS(DECR)
      COLAT = COS(YLATR)
      SIDEC = SIN(DECR)
      SILAT = SIN(YLATR)
       
C
C     calculate time correction
C
      ZN = 0.0
      TIMCOR = (4.0*(15.0*ZN+XLON)+ET)/60.0
      
      !write(*,*) 'd2str'
      

C
C......................................................................
C     calculate transmittance for characteristic day
C
C     fixed input parameters for subroutine strpsrb
C
      ROFF   = 0.0
      IWP    = 3
      ISUB   = 2
      TOP    = 500.0
      NCL    = 2
      INTVAL = 3
      CLOLWC = 0.0
C
C......................................................................
C     subroutine ATMOSPHERE - choose atmosphere by surface temperature
C     input  - TSTAT
C     output - LATMOS,TS
C
      CALL ATMOSPHERE(TSTAT,LATMOS,TS)
C
C......................................................................
C     construction of matrix TRANS
C
      TSOLAR = XHOR+TIMCOR    !para entrada com horario em UTC
C      TSOLAR = XHOR+ET        !para entrada com horario local
c      T1S = XHOR-1.5
C      T1A = T1S+TIMCOR
C      T2S = XHOR+1.5
C      T2A = T2S+TIMCOR
C
C     characteristic hour angle WI corresponding to the W1-W2 interval
C
      WSOLAR    = (12.00 - TSOLAR)*15.
      COWI  = COS(WSOLAR*CDR)
C      W1    = (12.00 - T1A)*15.
C      W2    = (12.00 - T2A)*15.
C      WI    = (W1 + W2)*0.5
C      COWI  = COS(WI*CDR)
C
C     characteristic solar zenith angle THETA (o)
C
      COSZEN = SIDEC*SILAT + CODEC*COLAT*COWI
      THETA  = ACOS(COSZEN)/CDR
      
C      IF((YLAT<=-17.75 .AND. YLAT>=-20)) THEN
C         WRITE(*,*) DEC, CODEC, COLAT, SIDEC, SILAT, TIMCOR, COWI, YLAT,
C     &   XLON, COSZEN, THETA, WSOLAR, TSOLAR, ET, XHOR
C      END IF
C
C......................................................................
C     subroutine STRPSRB - calculate transmittance for clear sky
C     input  - LATMOS,TS,ROFF,SFCALB,VIS,THETA,ICLOUD,IWP,ISUB,
C              TAUW,TOP,NCL,INTVAL,WH2O,CLOLWC,CDR
C     output - TRANS
C
c      write(*,*) xlon,ylat,tstat,rf
      IF(THETA.LE.90) THEN
      CALL STRPSRB(LATMOS,TS,ROFF,SFCALB,VIS,THETA,ICLOUD,IWP,ISUB,
     &    TAUW,TOP,NCL,INTVAL,tstat,rf,zstat,CLOLWC,CDR,TRANS,TDIR)
      ELSE
      TRANS=-2.0
      TDIR=-2.0
      ENDIF
C
C......................................................................
C
c      write(*,*) xlon,ylat,jday,xhor,coszen,theta,trans
c      PAUSE
      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
c
