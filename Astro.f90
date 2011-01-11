      SUBROUTINE ASTRO (JDAY,E0,DEC,ET)
!     calculate eccentricity correction, declination e equation of time
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     PI    : pi
!     CDR   : change degree to radians
!     DANG  : day angle (radians)
!     DEC   : declination (o)
!     E0    : eccentricity correction factor of earth
!     ET    : equation of time (minutes)
!     JDAY  : julian day
!     COnDA : cossine of n*day angle
!     SInDA : sine of n*day angle                                      
!......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
!     definition of constants bellow
      PI    = ACOS(-1.0)
      CDR  = PI/180.0

!     calculate day angle
      DANG  = 2.0*PI*(FLOAT(JDAY) - 1.0)/365.0

!     calculate sines and cossines of the day angle
      CO1DA = DCOS(DANG)
      CO2DA = DCOS(2*DANG)
      CO3DA = DCOS(3*DANG)
      SI1DA = DSIN(DANG)
      SI2DA = DSIN(2*DANG)
      SI3DA = DSIN(3*DANG)

!     calculate declination
      D0    = 0.006918
      DC1   = 0.399912
      DC2   = 0.006758
      DC3   = 0.002697
      DS1   = 0.070257
      DS2   = 0.000907
      DS3   = 0.001480
      DEC   = (D0 - DC1*CO1DA + DS1*SI1DA- DC2*CO2DA + DS2*SI2DA- DC3*CO3DA + DS3*SI3DA)/CDR

!     calculate eccentricity correction
      E00   = 1.000110
      EC1   = 0.034221
      EC2   = 0.000719
      ES1   = 0.001280
      ES2   = 0.000077
      E0    = E00 + EC1*CO1DA + ES1*SI1DA + EC2*CO2DA + ES2*SI2DA

!     calculate equation of time
      ET0    = 0.000075
      TC1    = 0.001868
      TC2    = 0.014615
      TS1    = 0.032077
      TS2    = 0.040890
      ET     = (ET0 + TC1*CO1DA - TS1*SI1DA - TC2*CO2DA - TS2*SI2DA) *229.18

      RETURN
      END
