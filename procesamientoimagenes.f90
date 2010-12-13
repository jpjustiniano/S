! Copyright (C) 2010, Juan Pablo Justiniano 

! Optimizar subrutina de calculo de elevacion

 Program ProcesamientoImagenes
 use netcdf
 implicit none
 Integer :: errorread
 character (4) :: ano
 character (2) :: mes, dia, hora, minu
 Integer :: iano, imes, idia, ihora, iminu, imesp=1
 Integer :: ifoto					! 1:media_hora; 2:south_full
 Real    :: az,el,ha,dec,soldst
 Logical :: 
 character (len = *) :: argument, filename
 character (len = 12) :: string
 real, parameter :: lat = -33., long = -70.

 !*********************************************************************** Fin declaracion Variables
 
 open (unit=6, file='log.txt')
 
 call get_command_argument(1, argument)
 OPEN (unit=8, file=trim(argument), status='old', IOSTAT=errorread)  
 IF(errorread/=0) print *, " error en apertura de archivo con lista de entrada" ; Exit 2
 
 100 read (8,*, IOSTAT=errorread) filename
 IF(errorread == -1) print *, " Terminado el procesamiento de imagenes"
 IF(errorread > 0) print *, " Error en lectura de nombre de archivo en archivo lista ", filename
 
 string = trim(filename)
 ano=string(1:4)
 mes=string(5:6)
 dia=string(7:8)
 hora=string(9:10)
 minu=string(11:12)
 foto=string(18:27)
 READ (iano,'(I4)') ano
 READ (imes,'(I2)') mes
 READ (idia,'(I2)') dia
 READ (ihora,'(I2)') hora
 READ (iminu,'(I2)') minu
 
 
 If (foto=='media_hora') then
	ifoto = 1
 ElseIf (foto=='south_full') then
	ifoto = 2
 Else 
	ifoto = 0
 End If
 
 call diajuliano (idia, imes, iano, diaj)
 call sunae(iano,diaj,ihora,lat,long,az,el,ha,dec,soldst)
 
 if (el < 7.0) then
  print *, " Imagen ", trim(filename)," eliminada, nocturna."
  write (6.*) trim(filename)," eliminada, nocturna."
  goto 100
 end if
  
 If (imesp /= imes) then
 
 
 end if
  
 OPEN (unit=10, file=trim(filename), status='old', IOSTAT=errorread)
 IF(errorread /= 0) print *, " Error en apertura de archivo de lista." ; Exit 2
 
 
 ! Lectura de datos de encabezado de archivo, dimensiones, lineas columnas
 ! Lectura pixel por pixel
 ! Correccion de datos
 
 
 
 
 
 
 
 
 contains
 subroutine check(status)
  integer, intent ( in) :: status
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop 2
  end if
 end subroutine check
 
 end Program ProcesamientoImagenes

!*********************************************************************** /Main

!*********************************************************************** Subrutines

FUNCTION norm_hora(value)
! Cambio de hora a forma decimal
IMPLICIT NONE

INTEGER, INTENT(IN) :: value 
REAL :: norm_hora 
CHARACTER(len=4) :: string
CHARACTER(len=2) :: hh, mm 
INTEGER :: h, m

string = ' '

WRITE (string,'(I4)') value
hh=string(1:2)
mm=string(3:4)
READ (hh,'(I2)') h
READ (mm,'(I2)') m
norm_hora = REAL(h)*3600.+REAL(m)*60.
END FUNCTION norm_hora


SUBROUTINE diajuliano (day, month, year, dayj)   
!This program calculates the day of year corresponding to a
!specified date. It illustrates the use of counting loops
!and the SELECT CASE construct.

IMPLICIT NONE
! Data dictionary: declare variable types, definitions, & units
INTEGER, INTENT(IN):: day          !Day (dd)
INTEGER, INTENT(IN) :: month        !Month (mm)
INTEGER, INTENT(IN) :: year         !Year (yyyy)
INTEGER, INTENT(out) :: dayj 		!Day of year
INTEGER :: i            			!Index,variable
INTEGER :: leap_day     			!Extra day for leap year

! Check for leap year, and add extra day if necessary
IF ( MOD(year,400) == 0 ) THEN
    leap_day = 1    ! Years divisible by 400 are leap years
ELSE IF ( MOD(year,100) == 0 ) THEN
    leap_day = 0    ! Other centuries are not leap years
ELSE IF ( MOD(year,4) == 0 ) THEN
    leap_day = 1    ! Otherwise every 4th year 1S a leap year
ELSE
    leap_day = 0    ! Other years are not leap years
END IF


! Calculate day of year
dayj= day
DO i = 1, month-1
    ! Add days in months from January to last month
    SELECT CASE (i)
    CASE (1,3,5,7,8,10,12)
    dayj = dayj + 31
    CASE (4,6,9,11)
    dayj = dayj + 30
    CASE (2)
    dayj = dayj + 28 + leap_day
    END SELECT
END DO

END SUBROUTINE diajuliano


 SUBROUTINE sunae(year,day,hour,lat,long,az,el,ha,dec,soldst)
 !Real, intent(in) :: year
 !Real, intent(in) :: day
 !Real, intent(in) :: hour
 !Real, intent(in) :: lat
 !Real, intent(in) :: long
 !Real, intent(out) :: az
 !Real, intent(out) :: el
 !Real, intent(out) :: ha
 !Real, intent(out) :: dec
 !Real, intent(out) :: soldst
 
! 	Subroutine to determine the Sun's position using Michalsky paper in 
!	Solar Energy Journal, Volume 40
!  
!  	corrections provided by Michalsky 4/23/07
!   (ftp://ftp.srrb.noaa.gov/pub/users/joe/asked_for_stuff/asunpos.f)
!
!  	Work with real variables and define some constants, including
!  	one to change between degs and radians.
implicit real (a-z)
data twopi,pi,rad/6.2831853,3.1415927,.017453293/
      
!   get the current julian date (actually add 2,400,000 for jd)
      delta=year-1949.
      leap=aint(delta/4.)
      jd=32916.5+delta*365.+leap+day+hour/24.
!   1st no. is mid. 0 jan 1949 minus 2.4e6; leap=leap days since 1949
!  the last yr of century is not leap yr unless divisible by 400
      if (amod(year,100.).eq.0.0.and.amod(year,400.).ne.0.0) jd=jd-1.

!   calculate ecliptic coordinates
      time=jd-51545.0
!   51545.0 + 2.4e6 = noon 1 jan 2000

!   force mean longitude between 0 and 360 degs
      mnlong=280.460+.9856474*time
      mnlong=mod(mnlong,360.)
      if(mnlong.lt.0.)mnlong=mnlong+360.

!   mean anomaly in radians between 0 and 2*pi
      mnanom=357.528+.9856003*time
      mnanom=mod(mnanom,360.)
      if(mnanom.lt.0.)mnanom=mnanom+360.
      mnanom=mnanom*rad

!   compute the ecliptic longitude and obliquity of ecliptic in radians
      eclong=mnlong+1.915*sin(mnanom)+.020*sin(2.*mnanom)
      eclong=mod(eclong,360.)
      if (eclong.lt.0.) eclong=eclong+360.
      oblqec=23.439-.0000004*time
      eclong=eclong*rad
      oblqec=oblqec*rad

!   calculate right ascension and declination
      num=cos(oblqec)*sin(eclong)
      den=cos(eclong)
      ra=atan(num/den)
!   force ra between 0 and 2*pi
      if (den.lt.0) then
          ra=ra+pi
      elseif (num.lt.0) then
          ra=ra+twopi
      endif

!   dec in radians
      dec=asin(sin(oblqec)*sin(eclong))

!   calculate Greenwich mean sidereal time in hours
      gmst=6.697375+.0657098242*time+hour 
!   hour not changed to sidereal time since 'time' includes
!   the fractional day 
      gmst = mod(gmst,24.)
      if(gmst.lt.0.) gmst=gmst+24.

!   calculate local mean sidereal time in radians 
      lmst=gmst+long/15.
      lmst=mod(lmst,24.)
      if(lmst.lt.0.) lmst=lmst+24.
      lmst=lmst*15.*rad

!   calculate hour angle in radians between -pi and pi
      ha=lmst-ra
      if(ha.lt.-pi) ha=ha+twopi
      if(ha.gt.pi) ha=ha-twopi

!   change latitude to radians
      lat=lat*rad

!   calculate azimuth and elevation
      el=asin(sin(dec)*sin(lat)+cos(dec)*cos(lat)*cos(ha))
      az=asin(-cos(dec)*sin(ha)/cos(el))

!   this puts azimuth between 0 and 2*pi radians
      if(sin(dec)-sin(el)*sin(lat).ge.0.) then
		if(sin(az).lt.0.) az=az+twopi
      else
      az=pi-az
      endif
!   if az=90 degs, elcritical=asin(sin(dec)/sin(lat))
!    elc=asin(sin(dec)/sin(lat))
!    if(el.ge.elc)az=pi-az
!    if(el.le.elc.and.ha.gt.0.)az=twopi+az

!   calculate refraction correction for US stan. atmosphere
!   need to have el in degs before calculating correction
      el=el/rad
!
      if(el.ge.19.225) then 
         refrac=.00452*3.51823/tan(el*rad)
      else if (el.gt.-.766.and.el.lt.19.225) then
         refrac=3.51823*(.1594+.0196*el+.00002*el**2)/ &
     &   (1.+.505*el+.0845*el**2)
      else if (el.le.-.766) then
         refrac=0.0
      end if

!   note that 3.51823=1013.25 mb/288 C
      el=el+refrac
!   elevation in degs
!
!   calculate distance to sun in A.U. & diameter in degs
      soldst=1.00014-.01671*cos(mnanom)-.00014*cos(2.*mnanom)
      soldia=.5332/soldst

!   convert az and lat to degs before returning
      az=az/rad
      lat=lat/rad
	 ha=ha/rad
	 dec=dec/rad

!   mnlong in degs, gmst in hours, jd in days if 2.4e6 added;
!   mnanom,eclong,oblqec,ra,and lmst in radians
 End subroutine
