     PROGRAM BRASILSR
! ++++++++++++++++++++++++++++++====+++++++++++++++++++++++++++++++++++++++
!                         Version 02/01/06                                 
! Variavel INTVAL ajustado para selecao da faixa espectral referente a todo
! espectro de radiacao solar de ondas curtas.
! Alterados os blocos data ISOL e WEIGHT na sub strpsrb para permitir 
! calculo do PAR no intervalo de 0,3 a 0,7microns.
! ++++++++++++++++++++++++++++++====+++++++++++++++++++++++++++++++++++++++

C     IMON   : month
C     IDAY   : day of month
C     IYEAR  : year
C     XLON   : longitude of the ground station (o)
C     XLAT   : latitude of the ground station  (o)
C     TSTAT  : temperature at ground station (K)
C     SFALB  : surface albedo                (%)
C     RFSURF : relative humidity             (%)
C     ZSTAT  : surface elevation             (m)
C     SSIW   : visibility                    (km)
C     NL     : number of cells in latitude direction
C     NC     : number of cells in longitude direction

!INPUTS
!  imagen	: Archivo temporal de lectura de datos
!  XLAT		: Latitudes
!  XLON		: Longitudes 
!  XALT		: Alturas (m)
!  XALB		: Albedo (%)
!  XUMI		: HR (%)
!  XTEM		: Temp (K)

!OUTPUTS
!  TRANS  : parameterized sky transmittance for broadband solar radiation in the atmosphere 
!  XTCR=TCLEAR
!  XTDIR=TDIR
!  XTCD=TCLOUD

!SUBROUTINE 
! TRANSMIT 	: Correccion de visivilidad y pasa a D2TR.
!  D2STR	: calculate transmittance for clear and cloudy sky
!   ASTRO	: calculate eccentricity correction, declination e equation of time
!   ATMOSPHERE	: choose atmosphere index (standard atmosphere as function of temperature)
!   STRPSRB	: calculate transmittance for clear and cloudy sky
! 
!**************************************************************************
SUBROUTINE TRANSMIT(JDAY,XLON,YLAT,TSTAT,SFALB,
     &                    RFSURF,ZSTAT,SSIW,XHOR,ICLOUD,TAUW,TRANS,TDIR)

      IF(XUMI(I,J)*100.0.GT.-98.0) THEN
	SSIW=XVIS(I,J)
       CALL TRANSMIT(NDIACAR,XLON(I,J),XLAT(I,J),XTEM(I,J),XALB(I,J),
     &  XUMI(I,J),XALT(I,J),SSIW,XHOR,0,0.0D0,TCLEAR,TDIR)
       XTCR=TCLEAR
       XTDIR=TDIR
       CALL TRANSMIT(NDIACAR,XLON(I,J),XLAT(I,J),XTEM(I,J),XALB(I,J),
     &  XUMI(I,J),XALT(I,J),SSIW,XHOR,1,100.0D0,TCLOUD,TDIR)
       XTCD=TCLOUD
      ELSE
       XTCR=-1.0
       XTCD=-1.0
       XTDIR=-1.0
      ENDIF

!**************************************************************************

C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE TRANSMIT(JDAY,XLON,YLAT,TSTAT,SFALB,
     &                    RFSURF,ZSTAT,SSIW,XHOR,ICLOUD,TAUW,TRANS,TDIR)
C
C      estimates global radiation in cloudy and clear-sky atmospheres
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     IMON   : month
C     IDAY   : day of month
C     XLON   : longitude of the ground station (o)
C     YLAT   : latitude of
C     TSTAT  : temperature at ground station (K)
C     SFALB  : surface albedo                (%)
C     RFSURF : relative humidity             (%)
C     ZSTAT  : surface elevation             (m)
C     SSIW   : visibility   (km)
C     VZGH   : control of exponential increasing
C     ZGH    : change meters to kilometers (m/km)
C     IYEAR  : year
C     JDAY   : julian day
C     RF     : relative humidity (0/1) 
C     PVSAT  : partial pressure of water vapor in saturated air [mbar]
C     WH2O   : precipitable water [cm]
C     AK     : auxiliary variable
C     SFCALB : surface albedo (0/1)
C     VIS    : visibility at ground station [km]
C     JDAY   : julian day                                             
C     TOTRAD : daily total global radiation [Wh.m-2]
C     TRANS  : parameterized sky transmittance for broadband solar radiation 
C              in the atmosphere 
C......................................................................
      CALL D2STR(XLON,YLAT,TSTAT,rf,zstat,SFCALB,VIS,JDAY,XHOR,ICLOUD,
     &           TAUW,TRANS,TDIR)


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
     

