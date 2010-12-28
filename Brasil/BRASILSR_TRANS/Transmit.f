C
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
C
      PARAMETER (NHOUR=24)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      WRITE(*,*) XHOR
C      write(*,*) JDAY,XLON,YLAT,TSTAT,SFALB,RFSURF,ZSTAT,SSIW,XHOR,
C     & ICLOUD,TAUW,TRANS
C      READ(*,*)
C
C     definition of constants bellow
C
      VZGH =  100.0
      ZGH  = 1000.0
C
C......................................................................
C
C     change relative humidity and albedo from percentual to relative
C     values
C
      SFCALB = SFALB/1000.0
      RF     = RFSURF/100.0
cC
cC     calculation of the precitable water as function of relative
cC     humidity, partial pressure of water vapor and temperature
cC
c
c      PVSAT   = EXP(26.23 - 5416.0/TSTAT)
c      WH2O  = 0.493*RF*PVSAT/TSTAT
cC
cC     test for precitable water less than 0.0000001 cm
cC
c      IF (WH2O .LE.0.0000001) THEN
c      WRITE (*,*)'NEGATIVE PRECIPITABLE WATER VAPOR',WH2O,
c     &' IS SET TO 0.000001 (cm)'
c      WH2O  = 0.0000001
c      ENDIF
C
C     calculation of the visibility at the station as function of
C     visibility and altitude
C
      AK    = LOG(VZGH/SSIW)/ZGH
      VIS   = SSIW * EXP( AK*ZSTAT )
C
C     test for visibility between 2 and 150 km
C                        
      IF (VIS .GT. 150.)  VIS = 150.
      IF (VIS .LT.   2.)  VIS =   2.
C
C......................................................................
C
C     subroutine D2STR - calculation of atmosphere transmittance for clear
C                        and cloudy sky
C     input  - XLON,YLAT,TSTAT,WH2O,SFCALB,VIS,JDAY
C     output - TRANS
C
      CALL D2STR(XLON,YLAT,TSTAT,rf,zstat,SFCALB,VIS,JDAY,XHOR,ICLOUD,
     &           TAUW,TRANS,TDIR)
C
C......................................................................
C
c      write(*,*) JDAY,XLON,YLAT,TSTAT,SFALB,RFSURF,ZSTAT,SSIW,IT,ICLOUD,
c     & TAUW,TRANS

      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
