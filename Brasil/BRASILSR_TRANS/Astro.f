C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE ASTRO (JDAY,E0,DEC,ET)
C
C     calculate eccentricity correction, declination e equation of
C     time
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     PI    : pi
C     CDR   : change degree to radians
C     DANG  : day angle (radians)
C     DEC   : declination (o)
C     E0    : eccentricity correction factor of earth
C     ET    : equation of time (minutes)
C     JDAY  : julian day
C     COnDA : cossine of n*day angle
C     SInDA : sine of n*day angle                                      
C......................................................................
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     definition of constants bellow
C
      PI    = ACOS(-1.0)
      CDR  = PI/180.0
C
C     calculate day angle
C
      DANG  = 2.0*PI*(FLOAT(JDAY) - 1.0)/365.0
C
C     calculate sines and cossines of the day angle
C
      CO1DA = DCOS(DANG)
      CO2DA = DCOS(2*DANG)
      CO3DA = DCOS(3*DANG)
      SI1DA = DSIN(DANG)
      SI2DA = DSIN(2*DANG)
      SI3DA = DSIN(3*DANG)
C
C     calculate declination
C
      D0    = 0.006918
      DC1   = 0.399912
      DC2   = 0.006758
      DC3   = 0.002697
      DS1   = 0.070257
      DS2   = 0.000907
      DS3   = 0.001480
      DEC   = (D0 - DC1*CO1DA + DS1*SI1DA
     &            - DC2*CO2DA + DS2*SI2DA
     &            - DC3*CO3DA + DS3*SI3DA)/CDR
C
C     calculate eccentricity correction
C
      E00   = 1.000110
      EC1   = 0.034221
      EC2   = 0.000719
      ES1   = 0.001280
      ES2   = 0.000077
      E0    = E00 + EC1*CO1DA + ES1*SI1DA + EC2*CO2DA + ES2*SI2DA
C
C     calculate equation of time
C
      ET0    = 0.000075
      TC1    = 0.001868
      TC2    = 0.014615
      TS1    = 0.032077
      TS2    = 0.040890
      ET     = (ET0 + TC1*CO1DA - TS1*SI1DA - TC2*CO2DA - TS2*SI2DA)
     &         *229.18
C
      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
