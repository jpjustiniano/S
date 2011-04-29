! Copyright (C) 2011, Juan Pablo Justiniano  <jpjustiniano@gmail.com>

program extraccion

 use netcdf

 implicit none
 
 ! Variables Programa
 integer, parameter :: Nimagenesdiamin = 5
 integer, parameter :: Nimagenesdiamax = 26 
 integer :: horaf, horaini, horasdia, diasmes
 integer :: largofilein, i , j, k
 integer :: ano, mes, dayj
 character(8)  :: fecha
 character(4)  :: cano
 character(2)  :: cmes
 character (50) :: argument 
 character(60) :: filename_out
 integer :: Nhora, NXf, NYf, Ndias, iLat, iLon
 integer,dimension(8) :: tiempo, tiempof, tiempoa
 integer :: TIDmax
 real :: Globfact, Dirfact, Latfact, Lonfact, horad, lat, lon
 real :: az,el,ha,dec,soldst, factormagico
 
 integer, dimension(:,:), allocatable :: Nimagen, Lat_CH1, Lon_CH1
 integer, dimension(:,:), allocatable :: Global, Directa, Difusa, DNI
 real, dimension (:), allocatable :: hora
 
 ! Variables NETCDF archivo entrada
 integer :: start(4), count(4), dimids(2), dimids4d(4), start_hor(4), count_hor(4) 
 integer :: Global_varid, Directa_varid, hora_varid, Lat_CH1_varid, Lon_CH1_varid
 integer :: xf_dimid, yf_dimid, hora_dimid, ncid, dia_dimid
 
 
!********************************************************** Fin declaracion Variables

factormagico = 0.913
 
 
 print *
 print *, '                 Extraccion de datos de radiacion Global, Difusa y Directa'
 print *, '                 *********************************************************'
 print *

!**************************************************** Lectura de archivo de entrada

 call date_and_time(DATE=fecha, VALUES=tiempo)
 tiempoa = tiempo
 
 call get_command_argument(1, argument)
 largofilein=len_trim(argument)
 
 If (.Not.(argument(largofilein-6: largofilein) /= 'prom.nc' .or. argument(largofilein-6: largofilein) /= '.rad.nc' ) ) then
	Write (*,*) '  Archivo de entrada no es base de datos de radiacion, Continuar ?? ', &
		'\',trim(argument(largofilein-6: largofilein)),'\'
	Read  (*,*)
 end if


call check( nf90_open(argument, nf90_nowrite, ncid) )

call check( nf90_inq_varid(ncid, "hora", hora_varid) )
call check( nf90_inq_varid(ncid, "Lat_CH1", Lat_CH1_varid) )
call check( nf90_inq_varid(ncid, "Lon_CH1", Lon_CH1_varid) )

call check( nf90_inq_dimid(ncid, "xf", xf_dimid) )
call check( nf90_inq_dimid(ncid, "yf", yf_dimid) )
call check( nf90_inq_dimid(ncid, "hora", hora_dimid) )
 
call check( nf90_inquire_dimension(ncid, xf_dimid, len = NXf) )
call check( nf90_inquire_dimension(ncid, yf_dimid, len = NYf) )
call check( nf90_inquire_dimension(ncid, hora_dimid, len = Nhora) )

call check( nf90_get_att(ncid, NF90_GLOBAL, "ano", ano))
call check( nf90_get_att(ncid, NF90_GLOBAL, "mes", mes))
call check( nf90_get_att(ncid, Lat_CH1_varid, "scale_factor", Latfact)) 
call check( nf90_get_att(ncid, Lon_CH1_varid, "scale_factor", Lonfact)) 
 
! allocate variables
allocate (hora(Nhora))
allocate (Lat_CH1 (NXf,NYf))       
allocate (Lon_CH1 (NXf,NYf))

 
 
 call check( nf90_get_var(ncid, hora_varid, hora) )
 call check( nf90_get_var(ncid, Lat_CH1_varid, Lat_CH1) )
 call check( nf90_get_var(ncid, Lon_CH1_varid, Lon_CH1) )
 
 Write (*,*) 'Latitud y longitud de punto a extraer datos: (grados decimales) '
 read (*,*) Lat, Lon
 
 iLat = nint(-100. * abs(Lat)); iLon = nint(-100. * abs(Lon)) ! -100 = -1/Latfact
 i=1; j=1
exter: Do 
	If (Lat_CH1(1,i)<=iLat) Then
	If (abs(Lat_CH1(1,i)-iLat)>abs(Lat_CH1(1,i-1)-iLat)) i= i-1
inter: Do 
		If (Lon_CH1(j,i)>=iLon) Then
			If (abs(Lon_CH1(j,i)-iLon)>abs(Lon_CH1(j-1,i)-iLon)) j= j-1
			print *, Lon_CH1(j,i),iLon, j, 'j', Lat_CH1(1,i),iLat, i, 'i'
			exit exter
		end if
		 If (j>= NXf) then
			write (*,*) ' Problemas al tratar de ajustar la Longitud!!' , i,j
			read *
			exit exter
		end if 
		j = j+1	
	end do inter
	end if
	If (i> NYf) then
		write (*,*) ' Problemas al tratar de ajustar la Latitud!!' , i,j
		exit exter
	end if 
	i=i+1	
 end do exter
 
 
 write (cano,'(I4)') ano
 If (mes < 10 ) then
	write(cmes,'(I1)') mes
	cmes='0'//trim(cmes)
 Else
	write(cmes,'(I2)') mes
 End if
 filename_out = argument(1:largofilein-3)//'.datos.txt'
 open (unit=15, file=filename_out)
 write (15,*) '  Latitud: ',Lat_CH1(1,i)/100.,' ,  Longitud: ',Lon_CH1(j,i)/100.
 
 start = (/ i, j, 1, 1 /)	 
	 
 If ( argument(largofilein-6: largofilein) == '.rad.nc') then
	print *, ' Base de datos no interpolada'
	write (15,*) 'Base de datos no interpolada ', argument
	
	call check( nf90_inq_varid(ncid, "Global", Global_varid) )
	call check( nf90_inq_varid(ncid, "Directa", Directa_varid) )
	 
	call check( nf90_get_att(ncid, Global_varid, "scale_factor", Globfact))
	call check( nf90_get_att(ncid, Directa_varid, "scale_factor", Dirfact))
	 
	call check( nf90_get_var(ncid, hora_varid, hora) )
	
	allocate (Global ( Nhora, 1))
	allocate (Directa ( Nhora, 1))
	allocate (Difusa ( Nhora, 1))
	allocate (DNI  ( Nhora, 1))
	
	count = (/ 1, 1, Nhora, 1 /)
	
	call check( nf90_get_var (ncid, Global_varid, Global, start, count) )
	call check( nf90_get_var (ncid, Directa_varid, Directa, start, count) )	
	
	where (Global<0) Global =0
	where (Directa<0) Directa =0
	Where (Global>15000) Global = 0
	where (Directa>15000) Directa =0
	
	do i = 2,Nhora-1
		If  (Global ( i, 1) ==0 .and.  Global ( i-1, 1) >0 .and.  Global ( i+1, 1) >0 ) then
			Global ( i, 1) = (Global ( i-1, 1) + Global ( i+1, 1))/2.
		end if 
		If  (Directa ( i, 1) ==0 .and.  Directa ( i-1, 1) >0 .and.  Directa ( i+1, 1) >0 ) then
			Directa ( i, 1) = (Directa ( i-1, 1) + Directa ( i+1, 1))/2.
		end if 
	end do	
	
	write (15,*)
	write (*,*)  ' Año Mes Dia Hora(UTC) Global   Directa   Difusa   DNI   el'
	write (15,*) 'Año;Mes;Dia;Hora(UTC);Global;Directa;Difusa;DNI'
exter2: Do i = 1 , 31
			call diajuliano (i, mes, ano, dayj)  
inter2:	Do j = 1, Nhora
			if (hora(j)> (i-1)*24 .and. hora(j) < (i)*24) then 
				horad = mod(hora(j),24.)
				call  sunae(ano,dayj,horad,iLat/100.,iLon/100.,az,el,ha,dec,soldst)
				write (*,200)  ano, mes, i, horad, Global(j,1)*Globfact, Directa(j,1)*Dirfact,&
					(Global(j,1)- Directa(j,1))*Globfact, Directa(j,1)*Dirfact/sin(el*3.141592/180.), el
				write (15,215)  ano, mes, i, horad, Global(j,1)*Globfact, Directa(j,1)*Dirfact,&
					(Global(j,1)- Directa(j,1))*Globfact, Directa(j,1)*Dirfact/sin(el*3.141592/180.)	
			end if
		end do inter2
	end do exter2
	
	
	
	
 else ! Prom.
 
	print *, ' Base de datos ajustada a medias horas'
	write (15,*) ' Base de datos ajustada a medias horas ', argument
	
	call check( nf90_inq_varid(ncid, "Global_hor", Global_varid) )
	call check( nf90_inq_varid(ncid, "Directa_hor", Directa_varid) )

	call check( nf90_inq_dimid(ncid, "dia", dia_dimid) )
	call check( nf90_inquire_dimension(ncid, dia_dimid, len = Ndias) )
	
	call check( nf90_get_att(ncid, NF90_GLOBAL, "horaini", horaini))
	call check( nf90_get_att(ncid, NF90_GLOBAL, "horasdia", horasdia))
	call check( nf90_get_att(ncid, NF90_GLOBAL, "diasmes", diasmes))
	 
	call check( nf90_get_att(ncid, Global_varid, "scale_factor", Globfact))
	call check( nf90_get_att(ncid, Directa_varid, "scale_factor", Dirfact))

	if (diasmes/=Ndias) print*, 'Diasmes: ',diasmes,' Ndias: ', Ndias
	allocate (Global ( Nhora, Ndias))
	allocate (Directa ( Nhora, Ndias))
	allocate (Difusa ( Nhora, Ndias))
	allocate (DNI  ( Nhora, Ndias))
	
	
	count = (/ 1, 1, Nhora, Ndias /)
	! dimids4d =  (/ xf_dimid_prom, yf_dimid_prom, hora_dimid_prom, diasmes_dimid_prom /) 
	call check( nf90_get_var (ncid, Global_varid, Global, start, count) )
	call check( nf90_get_var (ncid, Directa_varid, Directa, start, count) )
	
	where (Global<0) Global =0
	where (Directa<0) Directa =0
	Where (Global>15000) Global = 0
	where (Directa>15000) Directa =0
	
	Do 	j=1,Ndias
		do i = 2,Nhora-1
			If  (Global ( i, j) ==0 .and.  Global ( i-1, j) >0 .and.  Global ( i+1, j) >0 ) then
				Global ( i, j) = (Global ( i-1, j) + Global ( i+1, j))/2.
			end if 
			If  (Directa ( i, j) ==0 .and.  Directa ( i-1, j) >0 .and.  Directa ( i+1, j) >0 ) then
				Directa ( i, j) = (Directa ( i-1, j) + Directa ( i+1, j))/2.
			end if 
		end do	
	end do
	
	
	write (*,*) ' Año Mes Dia Hora(UTC) Global   Directa   Difusa   DNI   el'
	write (15,*)
	write (15,*) 'Año;Mes;Dia;Hora(UTC);Global;Directa;Difusa;DNI'
	Do 	j=1,Ndias
		call diajuliano (j, mes, ano, dayj)  
		Do i=1,Nhora
			call  sunae(ano,dayj, hora(i)/10.,iLat/100.,iLon/100.,az,el,ha,dec,soldst)
			write (*,200) ano, mes, j, hora(i)/10., Global(i,j)*Globfact, Directa(i,j)*Dirfact,&
			(Global(i,j)- Directa(i,j))*Globfact, Directa(i,j)*Dirfact/sin(el*3.141592/180.), el

			write (15,215)  ano, mes, j, hora(i)/10., Global(i,j)*Globfact, Directa(i,j)*Dirfact,&
			(Global(i,j)- Directa(i,j))*Globfact, Directa(i,j)*Dirfact/sin(el*3.141592/180.)
		End do
	end do
 End if

200 Format ('  ',I4,'  ',I2,'  ',I3,'  ',F4.1,'  ',F7.2,'  ',F7.2,'  ',F7.2,'  ',F7.2,'  ',F7.2)
215 Format (I5,';',I3,';',I3,';',F5.1,';',F8.2,';',F8.2,';',F8.2,';',F8.2)

 contains
 subroutine check(status)
  integer, intent ( in) :: status
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop 2
  end if
 end subroutine check
 
 
 subroutine diajuliano (day, month, year, dayj)   
!This program calculates the day of year corresponding to a specified date.
implicit none

integer, intent(in):: day          !Day (dd)
integer, intent(in) :: month        !Month (mm)
integer, intent(in) :: year         !Year (yyyy)
integer, intent(out) :: dayj 		!Day of year
integer :: i            			!Index,variable
integer :: leap_day     			!Extra day for leap year
 
! Check for leap year, and add extra day if necessary
IF ( mod(year,400) == 0 ) THEN
    leap_day = 1    ! Years divisible by 400 are leap years
ELSE IF ( mod(year,100) == 0 ) THEN
    leap_day = 0    ! Other centuries are not leap years
ELSE IF ( mod(year,4) == 0 ) THEN
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

 Subroutine sunae(year,day,hour, lat, long,az,el,ha,dec,soldst)  ! Sun's position, Michalsky
 implicit none
 integer, intent(in) :: year
 integer, intent(in) :: day
 Real, intent(in) :: hour
 Real, intent(in) :: lat
 Real, intent(in) :: long
 Real, intent(out) :: az
 Real, intent(out) :: el
 Real, intent(out) :: ha
 Real, intent(out) :: dec
 Real, intent(out) :: soldst
 
 real :: twopi, pi, rad
 real :: delta, leap, jd, time
 real :: mnlong, mnanom, eclong, oblqec, num, den, ra
 real :: gmst, lmst, latrad, elc, refrac, soldia

 data twopi,pi,rad/6.2831853,3.1415927,.017453293/
  
!   get the current julian date (actually add 2,400,000 for jd)
      delta=year-1949.
      leap=aint(delta/4.)
      jd=32916.5+delta*365.+leap+day+hour/24.
!   1st no. is mid. 0 jan 1949 minus 2.4e6; leap=leap days since 1949
!  the last yr of century is not leap yr unless divisible by 400
      if (mod(year,100).eq.0.0 .and. mod(year,400).ne.0.0) jd=jd-1.

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

!   hour not changed to sidereal time since 'time' includes the fractional day 
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
      latrad=lat*rad

!   calculate azimuth and elevation
      el=asin(sin(dec)*sin(latrad)+cos(dec)*cos(latrad)*cos(ha))
      az=asin(-cos(dec)*sin(ha)/cos(el))

!!   this puts azimuth between 0 and 2*pi radians
      if (sin(dec)-sin(el)*sin(latrad).ge.0.) then
		if(sin(az).lt.0.) az=az+twopi
      else
      az=pi-az
      endif
!   if az=90 degs, elcritical=asin(sin(dec)/sin(latrad))
    elc=asin(sin(dec)/sin(latrad))
    if(el.ge.elc)az=pi-az
    if(el.le.elc.and.ha.gt.0.)az=twopi+az

!   calculate refraction correction for US stan. atmosphere
!   need to have el in degs before calculating correction
      el=el/rad

      if(el.ge.19.225) then 
         refrac=.00452*3.51823/tan(el*rad)
      else if (el.gt.-.766.and.el.lt.19.225) then
         refrac=3.51823*(.1594+.0196*el+.00002*el**2)/(1.+.505*el+.0845*el**2)
      else if (el.le.-.766) then
         refrac=0.0
      end if

!   note that 3.51823=1013.25 mb/288 C
      el=el+refrac
!   elevation in degs
!   calculate distance to sun in A.U. & diameter in degs
      soldst=1.00014-.01671*cos(mnanom)-.00014*cos(2.*mnanom)
      soldia=.5332/soldst

!   convert az and lat to degs before returning
    az=az/rad
	ha=ha/rad
	dec=dec/rad

!   mnlong in degs, gmst in hours, jd in days if 2.4e6 added;
!   mnanom,eclong,oblqec,ra,and lmst in radians
End subroutine

end program extraccion
