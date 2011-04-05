! Copyright (C) 2011, Juan Pablo Justiniano  <jpjustiniano@gmail.com>
! Programa para el calculo de la radiacion Global, Difusa y Directa.

! Revisar:
! dimids
! Hora inicial sin datos. 


Program ProcesamientoImagenes_media_hora
 use netcdf
 !$ use OMP_LIb
 implicit none
 
 ! Variables Programa
 integer, parameter :: Nimagenesdiamin = 5
 integer, parameter :: Nimagenesdiamax = 26 
 real(8), parameter :: pi= acos(-1.0), cdr= pi/180.0
 integer :: horaf, horaini, horasdia, diasmes
 integer :: largofilein, lat, lon, i , j, k
 integer :: ano, mes, diaj
 character(8)  :: fecha
 character(4)  :: cano
 character(2)  :: cmes
 character (40) :: argument 
 character(25) :: filename_out
 integer :: Nhora, NXf, NYf
 integer,dimension(8) :: tiempo, tiempof, tiempoa
 integer :: TIDmax, tz
 
 ! Suane
  REAL ::  az,el,ha,dec,soldst
 
 
 integer, dimension(:), allocatable :: Nimagenesdia, ryo
 real, dimension (:), allocatable :: rxo, xo, yo, hora, hora_hor, Ex_h, Eryo
 integer, dimension (:,:), allocatable :: Nimagen, Glo, Dir, Lat_CH1, Lon_CH1
 real, dimension (:,:), allocatable :: horaimagenes 
 integer, dimension (:,:), allocatable :: imagenesdia
 integer, dimension (:,:), allocatable :: Horario 
 integer, dimension (:,:,:), allocatable ::  Global, Directa
 real, dimension (:,:,:), allocatable :: Global_hor, Directa_hor
 
 ! Variables NETCDF archivo entrada
 integer :: start(3), count(3), dimids(2), dimids4d(4), start_hor(4), count_hor(4) 
 integer :: Global_varid, Directa_varid, hora_varid, Lat_CH1_varid, Lon_CH1_varid
 integer :: xf_dimid, yf_dimid, hora_dimid, ncid

 ! Variables NETCDF archivo salida
 integer :: Global_varid_prom, Directa_varid_prom, hora_varid_prom, dia_varid_prom
 integer :: xf_dimid_prom, yf_dimid_prom, hora_dimid_prom, ncid_prom, diasmes_dimid_prom
 integer :: Lat_CH1_varid_prom, Lon_CH1_varid_prom

!********************************************************** Fin declaracion Variables
 
 print *
 print *, '                  Calculo de la interpolacion de Global y Directa'
 print *, '                  ************************************************'
 print *
 
 call date_and_time(DATE=fecha, VALUES=tiempo)
 tiempoa = tiempo

 TIDmax = 3    ! Numero de procesadores maximo a utilizar
 
 !$    TID = omp_get_num_procs()
 !$		If (TID>TIDmax) TID = TIDmax
 !$    call OMP_SET_NUM_THREADS(TID)
 
 call date_and_time(DATE=fecha, VALUES=tiempo)
 tiempoa = tiempo
 
 open (unit=16, file='log_prom.txt')
 
 call get_command_argument(1, argument)
 largofilein=len_trim(argument)
 If ( argument(largofilein-6: largofilein) /= '.rad.nc') then
	Write (*,*) '  Archivo de entrada no es base de datos de radiacion, Continuar ??', argument(largofilein-6: largofilein)
	Read  (*,*)
 End if
 
 call check( nf90_open(argument, nf90_nowrite, ncid) )

 call check( nf90_inq_varid(ncid, "hora", hora_varid) )

 call check( nf90_inq_varid(ncid, "Lat_CH1", Lat_CH1_varid) )
 call check( nf90_inq_varid(ncid, "Lon_CH1", Lon_CH1_varid) )
 call check( nf90_inq_varid(ncid, "Global", Global_varid) )
 call check( nf90_inq_varid(ncid, "Directa", Directa_varid) )

 call check( nf90_inq_dimid(ncid, "xf", xf_dimid) )
 call check( nf90_inq_dimid(ncid, "yf", yf_dimid) )
 call check( nf90_inq_dimid(ncid, "hora", hora_dimid) )
 
 call check( nf90_inquire_dimension(ncid, xf_dimid, len = NXf) )
 call check( nf90_inquire_dimension(ncid, yf_dimid, len = NYf) )
 call check( nf90_inquire_dimension(ncid, hora_dimid, len = Nhora) )

 call check( nf90_get_att(ncid, NF90_GLOBAL, "ano", ano))
 call check( nf90_get_att(ncid, NF90_GLOBAL, "mes", mes))
 
 ! allocate variables
 allocate (hora(Nhora))
 allocate (Lat_CH1 (NXf,NYf))       
 allocate (Lon_CH1 (NXf,NYf))
 
 call check( nf90_get_var(ncid, hora_varid, hora) )
 call check( nf90_get_var(ncid, Lat_CH1_varid, Lat_CH1) )
 call check( nf90_get_var(ncid, Lon_CH1_varid, Lon_CH1) )
 
 
 select case (mes) !! revisar horas finales mayores a 24
 case (1)
 horaini = 10
 horasdia = 14
 diasmes =  31
 case (2)
 horaini = 10
 horasdia = 14
 diasmes = 29 
 case (3)
 horaini = 11
 horasdia = 12
 diasmes = 31
 case (4)
 horaini = 11
 horasdia = 11
 diasmes = 30
 case (5)
 horaini = 12
 horasdia = 9
 diasmes = 31
 case (6)
 horaini = 12
 horasdia = 9
 diasmes = 30
 case (7)
 horaini = 12
 horasdia = 9
 diasmes = 31
 case (8)
 horaini = 11
 horasdia = 11
 diasmes = 31
 case (9)
 horaini = 11
 horasdia = 12
 diasmes = 30
 case (10)
 horaini = 10
 horasdia = 13
 diasmes = 31
 case (11)
 horaini = 10
 horasdia = 13
 diasmes = 30
 case (12)
 horaini = 10
 horasdia = 14
 diasmes = 31
 case default
	print *, ' Mes fuera de rango!! ', mes
	stop 2
 end Select
 

  !Creacion de archivo de salida
 
  write (cano,'(I4)') ano
 If (mes < 10 ) then
	write(cmes,'(I1)') mes
	cmes='0'//trim(cmes)
 Else
	write(cmes,'(I2)') mes
 End if
 filename_out = cano//cmes//'.media_hora.prom.nc'
  
 call check( nf90_create(filename_out,NF90_64BIT_OFFSET, ncid_prom) ) ! ncid, archivo de salida recortado
 ! Dimensiones
 call check( nf90_def_dim(ncid_prom, "xf", NXf, xf_dimid_prom) )  ! Define the dimensions. 
 call check( nf90_def_dim(ncid_prom, "yf", NYf, yf_dimid_prom) )
 call check( nf90_def_dim(ncid_prom, "hora", horasdia, hora_dimid_prom) )
 call check( nf90_def_dim(ncid_prom, "dia", diasmes, diasmes_dimid_prom) )
 
 dimids = (/ xf_dimid_prom, yf_dimid_prom /)
 dimids4d =  (/ xf_dimid_prom, yf_dimid_prom, hora_dimid_prom, diasmes_dimid_prom /) 
 
 ! Variables
 call check( nf90_def_var(ncid_prom,"hora", NF90_SHORT, hora_dimid_prom, hora_varid_prom) ) 
 call check( nf90_def_var(ncid_prom,"dia", NF90_SHORT, diasmes_dimid_prom, dia_varid_prom) ) 
 call check( nf90_def_var(ncid_prom,"Global_hor", NF90_SHORT, dimids4d, Global_varid_prom) )
 call check( nf90_def_var(ncid_prom,"Directa_hor", NF90_SHORT, dimids4d, Directa_varid_prom) )
 call check( nf90_def_var(ncid_prom, "Lat_CH1", NF90_SHORT, dimids, Lat_CH1_varid_prom) )
 call check( nf90_def_var(ncid_prom, "Lon_CH1", NF90_SHORT, dimids, Lon_CH1_varid_prom) )
 
 ! Atributos globales    
 call check( nf90_put_att(ncid_prom, nf90_global, "imagenes", "media_hora") )
 call check( nf90_put_att(ncid_prom, nf90_global, "Procesadox", "JPJ") )
 call check( nf90_put_att(ncid_prom, nf90_global, "date", fecha) )
 call check( nf90_put_att(ncid_prom, nf90_global, "horaini", horaini) )
 call check( nf90_put_att(ncid_prom, nf90_global, "horasdia", horasdia) )
 call check( nf90_put_att(ncid_prom, nf90_global, "diasmes", diasmes) )
 call check( nf90_put_att(ncid_prom, nf90_global, "ano", ano))
 call check( nf90_put_att(ncid_prom, nf90_global, "mes", mes))
 call check( nf90_put_att(ncid_prom, nf90_global, "NXf", NXf) )
 call check( nf90_put_att(ncid_prom, nf90_global, "NYf", NYf) )
 
 ! Atributos variables
 call check( nf90_put_att(ncid_prom, Lat_CH1_varid_prom, "units", "degrees_north"))
 call check( nf90_put_att(ncid_prom, Lat_CH1_varid_prom, "scale_factor", 0.01) )
 call check( nf90_put_att(ncid_prom, Lat_CH1_varid_prom, "_CoordinateAxisType", "Lat_CH1") )

 call check( nf90_put_att(ncid_prom, Lon_CH1_varid_prom, "units", "degrees_east") )
 call check( nf90_put_att(ncid_prom, Lon_CH1_varid_prom, "scale_factor", 0.01) )
 call check( nf90_put_att(ncid_prom, Lon_CH1_varid_prom, "_CoordinateAxisType", "Lon_CH1") )
 
 call check( nf90_put_att(ncid_prom, hora_varid_prom, "units", "UTC_hours"))
 call check( nf90_put_att(ncid_prom, hora_varid_prom, "_CoordinateAxisType", "time"))
 
 call check( nf90_put_att(ncid_prom, Global_varid_prom, "units", "W/m2") )
 call check( nf90_put_att(ncid_prom, Global_varid_prom, "missing_value", -32767) )
 call check( nf90_put_att(ncid_prom, Global_varid_prom, "valid_min", -32767) )
 call check( nf90_put_att(ncid_prom, Global_varid_prom, "valid_max", 32767) )
 call check( nf90_put_att(ncid_prom, Global_varid_prom, "scale_factor", 0.1) )
 call check( nf90_put_att(ncid_prom, Global_varid_prom, "_CoordinateAxes", "time Lat_CH1 Lon_CH1") )
 call check( nf90_put_att(ncid_prom, Global_varid_prom, "standard_name", "Radiación Global") )
 
 call check( nf90_put_att(ncid_prom, Directa_varid_prom, "units", "W/m2") )
 call check( nf90_put_att(ncid_prom, Directa_varid_prom, "missing_value", -32767) )
 call check( nf90_put_att(ncid_prom, Directa_varid_prom, "valid_min", -32767) )
 call check( nf90_put_att(ncid_prom, Directa_varid_prom, "valid_max", 32767) )
 call check( nf90_put_att(ncid_prom, Directa_varid_prom, "scale_factor", 0.1) )
 call check( nf90_put_att(ncid_prom, Directa_varid_prom, "_CoordinateAxes", "time Lat_CH1 Lon_CH1") )
 call check( nf90_put_att(ncid_prom, Directa_varid_prom, "standard_name", "Radiación Global") )
 
 ! End define mode.
 call check( nf90_enddef(ncid_prom) )
 
 ! / Fin creacion archivo de salida
 

 allocate (Nimagenesdia(diasmes))
 allocate (Nimagen(diasmes,Nimagenesdiamax))
 allocate (xo(horasdia))
 allocate (yo(horasdia))
 allocate (horaimagenes(Nimagenesdiamax, diasmes))
 allocate (Ex_h(horasdia))
 
 ! Nimagenesdia (dia) : var. temporal con numero de imagenes x dia i
 ! Nimagen (dia, imagendia) : Numero de imagen del total
 ! horaimagen (dia, imagendia) : hora decimal de la imagen

 Nimagenesdia = 0
 do i =1, diasmes
 inte:Do j =1, Nhora   !! Optimizable.
		if (hora(j)> (i-1)*24 .and. hora(j) < (i)*24) then
			Nimagenesdia(i)=Nimagenesdia(i)+1
			Nimagen(i,Nimagenesdia(i))=j
			horaimagenes(Nimagenesdia(i), i) =  hora(j)
		else if (hora(j) > (i)*24) then
			exit inte
		end if
	end do	inte
 end do
 

 horaf = horaini + horasdia
 tz = 3
 xo(1)= horaini+0.5 ! Hora en UTC
 Do i=2,horasdia
	xo(i)= xo(i-1)+1
 end do
 call check( nf90_put_var(ncid_prom, hora_varid_prom, xo) ) 
 
 start = (/ 1, 1, 1 /)
 count = (/ NXf, NYf, 1 /)
 
 start_hor = (/ 1, 1, 1 ,1/)
 count_hor = (/ NXf, NYf, horasdia, 1 /) !! revisar uso de dimensiones..
 !				Lon  Lat	hora   dia
 
 Do i = 1, diasmes ! Se guarda la informacion de todas las imagenes del dia en las matrices 3D
	call check( nf90_put_var(ncid_prom, dia_varid_prom, i,start=(/i/)) )
	print *, ' Dia : ',i, ', imagenes:', Nimagenesdia(i)
	If (Nimagenesdia(i)<Nimagenesdiamin) cycle
	
	allocate (Global  (NXf, NYf, Nimagenesdia(i)) )
	allocate (Directa (NXf, NYf,Nimagenesdia(i)) )
	allocate (Global_hor (NXf, NYf, horasdia))
	allocate (Directa_hor (NXf, NYf, horasdia))
	allocate (rxo(Nimagenesdia(i)))
	allocate (ryo(Nimagenesdia(i)))
	allocate (Eryo(Nimagenesdia(i)))
	
	Global_hor = 0
	Directa_hor = 0
	

	start (3) = Nimagen(i,1)	
	count (3) = Nimagenesdia(i)
	call check( nf90_get_var(ncid, Global_varid, Global, start=start, count=count))
	call check( nf90_get_var(ncid, Directa_varid, Directa, start, count ))
	
	
	rxo = horaimagenes(:,i)-(i-1)*24
	
	call diajuliano (i, mes, ano, diaj)   ! entrada de reales en ves de enteros.

	
	 exter: Do lon = 1, Nxf			! Se procesa en NXxNY la matriz 3D diaria
		inter: Do lat = 2, Nyf
				if (Global (lon,lat, Nimagen(i,2))<0) cycle inter
				if (Directa (lon,lat, Nimagen(i,2)) < 0) cycle inter
				
				do k = 1, horasdia
				call sunae (ano,diaj,xo(k),Lat_CH1(lon,lat)/100., Lon_CH1(lon,lat)/100., az,el,ha,dec,soldst)
					Ex_h(K) = 1367.*(1+0.033*cos(360.*diaj/365.))*cos((90.-el)*cdr)
				end do	
				do k = 1,Nimagenesdia(i)
				call sunae (ano,diaj,horaimagenes(k,i),Lat_CH1(lon,lat)/100., Lon_CH1(lon,lat)/100., az,el,ha,dec,soldst)
					Eryo(K) = 1367.*(1+0.033*cos(360.*diaj/365.))*cos((90.-el)*cdr)
				end do	
				
				! Agregar puntos de salida y puesta del sol
				
				ryo = Global(lon,lat,Nimagen(i,1):Nimagen(i,Nimagenesdia(i)))
				where (ryo<0) ryo=0
				!print *, Nimagen(i,1),Nimagen(i,Nimagenesdia(i)), el, ((90.-el)*cdr)
				!print *, Global(lon,lat,Nimagen(i,1):Nimagen(i,Nimagenesdia(i)))
				call InterLinealPond( horaini,horasdia, Nimagenesdia(i),rxo,ryo, Eryo, xo, Ex_h, yo )

				Global_hor(lon,lat,:) = yo
				!print *,'yo:', Global_hor(lon,lat,:)
				
				
				ryo = real(Directa(lon,lat,:))
				where (ryo<0) ryo=0
				
				call InterLinealPond( horaini,horasdia, Nimagenesdia(i),rxo,ryo, Eryo, xo, Ex_h, yo )
				
				Directa_hor(lon,lat,:) = yo
				!print *, 'yo:',Directa_hor(lon,lat,:)
				
		end do inter
	 end do exter
	 
	start_hor (4) = i
	where (Global_hor>32700)  Global_hor=32700  	 ! Filtro de maximos
	where (Directa_hor>32700)  Directa_hor=32700
	call check( nf90_put_var(ncid_prom, Global_varid_prom, Global_hor, start_hor,count_hor))
	call check( nf90_put_var(ncid_prom, Directa_varid_prom,Directa_hor,start_hor,count_hor))
	
 
 ! Contador imager ned med hor
 	
 	
 	deallocate (Global)
	deallocate (Directa )
	deallocate (Global_hor )
	deallocate (Directa_hor )
	deallocate (rxo)
	deallocate (ryo)	
	call date_and_time(DATE=fecha, VALUES=tiempof)
	print *
	write (*,200) tiempof(5) - tiempoa(5) ,tiempof(6)- tiempoa(6),tiempof(7) - tiempoa(7)
	print *
	tiempoa = tiempof
 End do
call date_and_time(DATE=fecha, VALUES=tiempof)
print *, 'Tiempo Total:'
write (*,200) tiempof(5) - tiempo(5) ,tiempof(6)- tiempo(6),tiempof(7) - tiempo(7)
print *

200 Format ('   Tiempo procesamiento: ',I3' hr., ',I3 'min., ',I3, 'sec.')
300 Format ('   Tiempo procesamiento: ',I3 'min., ',I3, 'sec.')

 contains
subroutine InterLinealPond ( horaini,horasdia,Nimagenesdia,rxo,ryo, Eryo, xo, Ex_h, yo )
implicit none

integer, intent(in) ::horaini,horasdia,Nimagenesdia
real, dimension (:), intent(in) :: rxo, xo, Ex_h, Eryo
integer, dimension (:), intent(in) :: ryo
real, dimension (:), intent(out) :: yo
integer :: i = 1, j= 1
real :: k1, k2

! write (*,*) 'rxo:',rxo
! write (*,*) 'ryo:',ryo
! write (*,*) 'xo:',xo
! write (*,*) 'Ex_h:',Ex_h
 yo=0.
! read (*,*)
 
 
exter: do i = 1, horasdia
	inter: Do j = 1, Nimagenesdia
		if (xo(i) > rxo(j) ) then
			if (j>1) then 
			 k1 = ryo(j-1)/Eryo(j-1)
			 k2 = ryo(j)/Eryo(j)
			 yo = (((xo(i)-rxo(j-1))/(rxo(j)-rxo(j-1)))*k1 +((rxo(j)-xo(i))/(rxo(j)-rxo(j-1)))*k2)*Ex_h(i)
			else 
			 yo = (ryo(j)/Eryo(j))*Ex_h(i)
			end if
		exit inter
		end if 
	end do inter
end do exter


end subroutine InterLinealPond
 
 
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

End program ProcesamientoImagenes_media_hora
