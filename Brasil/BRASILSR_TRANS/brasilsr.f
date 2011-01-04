      PROGRAM BRASILSR
! ++++++++++++++++++++++++++++++====+++++++++++++++++++++++++++++++++++++++
!                         Version 02/01/06                                 
! Variavel INTVAL ajustado para selecao da faixa espectral referente a todo
! espectro de radiacao solar de ondas curtas.
! Alterados os blocos data ISOL e WEIGHT na sub strpsrb para permitir 
! calculo do PAR no intervalo de 0,3 a 0,7microns.
!
! VERSÃO DESENVOLVIDA PARA LINUX utilizando as imagens do GOES12.
! 
! Calcula as transmitancias global, difusa e direta para mês em estudo
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

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 CJDAY
      CHARACTER*2 CCJDAY
      CHARACTER*3 CMON,CCCJDAY,CCHORI
      CHARACTER*4 CYEAR,CHORI
      CHARACTER*80 CDIR,CFILE,CDESC,CANAL,CSAT,CLATI,CLONG
      CHARACTER*80 CDIRE,CDIFU,CKD,CTRANS,CVISIV
      CHARACTER*80 CALBE,CALTI,CUMID,CTEMP,CGLOB,CIMAG,CKT,CIMAGEM
      REAL*8 XTCD,XTCR,XTDIR
      REAL*8 XLAT(1784,1180),XLON(1784,1180),XVIS(1784,1180)
      REAL*8 XALT(1784,1180),XALB(1784,1180),XUMI(1784,1180)
      REAL*8 XTEM(1784,1180)
      INTEGER*2 LEITURA(1784,1180)

      INTEGER*2 IHOR(48)
      INTEGER IERRO, IERR

!**************************************************************************
C.....input from description parameters

      OPEN(UNIT=10,FILE='gkss.dat')
      DO
        READ(10,1001,IOSTAT=IERR) CDIR
         IF (IERR /= 0) EXIT
        READ(10,*) CMON,IYEAR
        WRITE(CYEAR,'(I4)') IYEAR
        CIMAGEM='imagens'//CMON//CYEAR//'.dat'
        PRINT*, CIMAGEM
      !ABRINDO LISTA COM IMAGENS	
      OPEN(UNIT=14,FILE=CIMAGEM)
C     OBTENDO DADOS DO ARQUIVO DESCRITOR
      CDESC=TRIM(CDIR)//CYEAR//'/'//CMON//'/description/'//CMON//
     &  CYEAR//'.des'
      PRINT*, CDESC
      OPEN(11,FILE=CDESC)
	      READ(11,*) CANAL
	      READ(11,*) CSAT
	      READ(11,*) ILINI
	      READ(11,*) ICOLI
	      READ(11,*) NL
	      READ(11,*) NC
	      NLC = 2*NL*NC
      CLOSE(11)

C     LEITURA DOS ARQUIVOS DE LATITUDE E LONGITUDE E CLIMATICOS
C
      CLATI=TRIM(CDIR)//CYEAR//'/'//CMON//'/description/'
     & //CMON//CYEAR//'.lat'
      CLONG=TRIM(CDIR)//CYEAR//'/'//CMON//'/description/'
     & //CMON//CYEAR//'.lon'
      CALTI=TRIM(CDIR)//CYEAR//'/'//CMON//'/prop/'
     & //CMON//CYEAR//'.alt'
      CALBE=TRIM(CDIR)//CYEAR//'/'//CMON//'/prop/'
     & //CMON//CYEAR//'.alb'
      CUMID=TRIM(CDIR)//CYEAR//'/'//CMON//'/prop/'
     & //CMON//CYEAR//'.umi'
      CTEMP=TRIM(CDIR)//CYEAR//'/'//CMON//'/prop/'
     & //CMON//CYEAR//'.tem'
      CVISIV=TRIM(CDIR)//CYEAR//'/'//CMON//'/prop/'
     & //CMON//CYEAR//'.vis'


      WRITE(*,*) CVISIV
      OPEN(99,FILE=CVISIV,FORM='unformatted',RECL=NLC,
     & STATUS='OLD',ACCESS='direct')
      READ(99,REC=1) ((LEITURA(I,J),I=1,NC),J=1,NL)
      CALL MOVEBITES(LEITURA,NL,NC)
      XVIS = LEITURA/100.0
      CLOSE (99)
     
      WRITE(*,*) CLATI
      OPEN(1,FILE=CLATI,FORM='unformatted',RECL=NLC,
     & STATUS='OLD',ACCESS='direct')
      READ(1,REC=1) ((LEITURA(I,J),I=1,NC),J=1,NL)
      CALL MOVEBITES(LEITURA,NL,NC)
      XLAT = LEITURA/100.0
      CLOSE (1)

      WRITE(*,*) CLONG
      OPEN(2,FILE=CLONG,FORM='unformatted',RECL=NLC,STATUS='OLD',
     & ACCESS='direct')
      READ(2,REC=1,ERR=1111) ((LEITURA(I,J),I=1,NC),J=1,NL)
      CALL MOVEBITES(LEITURA,NL,NC)
      XLON = LEITURA/100.0
      CLOSE (2)

      WRITE(*,*) CALTI
      OPEN(3,FILE=CALTI,FORM='unformatted',RECL=NLC,STATUS='OLD',
     & ACCESS='direct')
      READ(3,REC=1,ERR=1112) ((LEITURA(I,J),I=1,NC),J=1,NL)
      CALL MOVEBITES(LEITURA,NL,NC)
      XALT = LEITURA/1.0     ! altitude entra em metros
      CLOSE (3)

      WRITE(*,*) CALBE
      OPEN(4,FILE=CALBE,FORM='unformatted',RECL=NLC,STATUS='OLD',
     & ACCESS='direct')
      READ(4,REC=1,ERR=1113) ((LEITURA(I,J),I=1,NC),J=1,NL)
      CALL MOVEBITES(LEITURA,NL,NC)
      XALB = LEITURA/1.0      ! albedo entra em porcentagem...... ORIGINAL_
      CLOSE (4)
 
      WRITE(*,*) CUMID
      OPEN(11,FILE=CUMID,FORM='unformatted',RECL=NLC,STATUS='OLD',
     & ACCESS='direct')
      READ(11,REC=1,ERR=1114) ((LEITURA(I,J),I=1,NC),J=1,NL)
      CALL MOVEBITES(LEITURA,NL,NC)
      XUMI = LEITURA/100.0      !umidade entra em %*100...
      CLOSE (11)

      WRITE(*,*) CTEMP
      OPEN(12,FILE=CTEMP,FORM='unformatted',RECL=NLC,STATUS='OLD',
     & ACCESS='direct')
      READ(12,REC=1,ERR=1115) ((LEITURA(I,J),I=1,NC),J=1,NL)
      CALL MOVEBITES(LEITURA,NL,NC)
      XTEM = (LEITURA/100.0) + 273.15 !temperatura entra em Celsius*100
      CLOSE (12)

      DO  KK=1,48
        IHOR(KK)=-1
      END DO !169
      JJJ=1

      DO
       READ(14,*,IOSTAT=IERRO) JDAY,IYEAR,IHORI
       IF (IERRO /= 0) EXIT
       JJ=JJJ
       DO KK=1,JJ
        IF(IHORI.NE.IHOR(KK)) THEN
         IF(IHOR(KK).EQ.-1) THEN
          IHOR(KK)=IHORI
          JJJ=JJ+1
          EXIT
         ENDIF
        ELSE
         EXIT
        ENDIF
       END DO !179
      END DO !174
      CLOSE (14)


      DO K=1,JJJ-1

      XHOR=DINT(IHOR(K)/100.0D0)
C
C     ABERTURA DOS ARQUIVOS DE TRANSMITANCIAS
C
      IF(IHOR(K).LT.1000) THEN
       WRITE(CCHORI,'(I3)') IHOR(K)
       CHORI='0'//CCHORI
      ELSE
       WRITE(CHORI,'(I4)') IHOR(K)
      ENDIF
      CTRANS=TRIM(CDIR)//CYEAR//'/'//CMON//'/trans/trans'//CMON//
     & CYEAR//'.'//CHORI
      WRITE(*,*) CTRANS, IHOR(K)
      OPEN(18,FILE=CTRANS,STATUS='UNKNOWN')

      XHOR=XHOR+(IHOR(K)-100.0D0*XHOR)/60.0D0
      CALL MONTH(CMON,IMON)
      CALL DIACAR(IMON,IYEAR,NDIACAR)
      DO J=1,NL
      DO I=1,NC
      if (xtem(i,j) > 320) then
            print*, i,j,XTEM(I,J),XALB(I,j),XUMI(I,J),XALT(I,J)
            pause
      end if

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
      WRITE(18,1005) XLON(I,J),XLAT(I,J),XTCR,XTCD,XTDIR

      END DO !235
      END DO !233
      CLOSE(18)

      END DO !212
      END DO !66
      CLOSE(10)              
     
 1001 FORMAT(A80)     
 1002 FORMAT(A3,1X,I4)     
 1003 FORMAT(A2)     
 1004 FORMAT(A3,1x,f8.2,3(I5),1X,2(F8.2))     
 1005 FORMAT(2(1X,F8.2),3(1X,F8.3))     
  104 FORMAT(3(I4),5(F9.3),i6)     
  111 FORMAT(2(I4),3(F9.3),i6,2(i6),2(F9.3))     

C*************************************************************************
C   Mensagens de erro de leitura
 1110 PRINT*, 'Erro de leitura do arquivo de latitude'
      STOP
 1111 PRINT*, 'Erro de leitura do arquivo de longitude'
      STOP
 1112 PRINT*, 'Erro de leitura do arquivo de altitude'
      STOP
 1113 PRINT*, 'Erro de leitura do arquivo de albedo'
      STOP
 1114 PRINT*, 'Erro de leitura do arquivo de umidade relativa'
      STOP
 1115 PRINT*, 'Erro de leitura do arquivo de temperatura'
      STOP
      END
Correccion de visivilidad y pasa a D2TR.
