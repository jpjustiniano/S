C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE ATMOSPHERE(TSTAT,LATMOS,TS)
C
C     choose atmosphere index (standard atmosphere as function of
C     temperature)
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     TSTAT   : temperature at ground station (K)
C     LATMOS  : atmosphere index
C               1 - tropical
C               2 - midlatitude summer
C               3 - midlatitude winter
C               4 - subartic summer
C               5 - subartic winter
C     TS      : difference between surface and characteristical temperatures
C     TnSF    : temperature for atmosphere type 'n' (K)                
C......................................................................
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     determine characteristical temperatures for each type of
C     atmosphere
C
      T1SFC = 300.0
      T2SFC = 294.0
      T4SFC = 287.0
      T3SFC = 272.2
      T5SFC = 257.1
C
C     choose atmosphere index and calculate difference between surface
C     and characteristical temperatures
C
         LATMOS = 1
         TS     = TSTAT - T1SFC
      IF(TSTAT.LT.297.0 .AND. TSTAT.GE.290.5) THEN
         LATMOS = 2
         TS     = TSTAT - T2SFC
      ENDIF
      IF(TSTAT.LT.290.5 .AND. TSTAT.GE.279.5) THEN
         LATMOS = 4
         TS     = TSTAT - T4SFC
      ENDIF
      IF(TSTAT.LT.279.5 .AND. TSTAT.GE.264.5) THEN
         LATMOS = 3
         TS     = TSTAT - T3SFC
      ENDIF
      IF(TSTAT.LT.264.5) THEN
         LATMOS = 5
         TS     = TSTAT - T5SFC
      ENDIF
C
      RETURN
      END
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
