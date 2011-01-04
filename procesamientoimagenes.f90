! Copyright (C) 2010, Juan Pablo Justiniano  <jpjustiniano@gmail.com>
! Programa para el filtraje y almacenamiento de imagenes satelitales en un archivo mensual.
! Funciona con horas en UTC.

! gfortran -c -I/usr/include "%f"
! gfortran -I/usr/include -L/usr/lib -lnetcdff -lnetcdf -o "%e" "%f"
! gfortran -o "%e" "%f" -L/usr/lib  -lnetcdf   ?? de manual...

! Variable coordinada de lat y long real.
! Coordinar las variables de entrada en las mismas coordenadas
! Procesamiento de solo Chile. Eliminar fuera de fronteras.

!Canal 1: [1728, 9020]
!Canal 4: [432, 2255]

!Imagen Media Hora:
!Canal 1_SI: [1, 3700]
!Canal 1_ID: [1728, 5800]
 
!Canal 4_SI: [1, 925]
!Canal 4_ID: [432, 1450]

 Program ProcesamientoImagenes
 use netcdf
 implicit none
 
 real :: latitud = -33., longitud = -70.
 ! Variables Programa
 Integer :: errorread
 character (4) :: ano
 character (2) :: mes, cdia, hora, minu
 character (2) :: foto
 character(8)  :: fecha
 real :: iano, imes, idia, ihora, ihorat, iminu, imesp=0., diaj, hora_utc
 Integer :: ifoto					! 1:media_hora; 2:south_full
 Real    :: az,el,ha,dec,soldst
 character (30) :: argument  ! Nombre de archivo de lista
 character (30) :: filename  ! Nombres dentro de la lista
 character (12) :: filename_out
 character (len = 27) :: string
 integer :: i=1,j=1, rec=1, noct= 0
 logical :: flag1= .false.
 integer,dimension(8) :: tiempo, tiempof
 
 ! Variables NETCDF 
 integer :: ncid, ncid_in, dia =31, status
 integer, parameter :: NDIMS = 3, NDIMS_IN = 2	      ! We are writing 2D data.
 integer, parameter :: NX = 526, NY = 432 , NXf = 2101, NYf = 1728 

 real :: x(NX), y(NY), xf(NXf), yf(NYf)
 integer, dimension(:,:), allocatable :: CH1_out, CH4_out   ! Matrices de archivo original
 integer, dimension(:,:), allocatable :: CH4_in , CH1_in        !  Matrices de archivo recortado
 real, dimension(:,:), allocatable :: CH1_max, CH1_min, CH4_max,CH4_min  !  Matrices max min
 integer :: x_dimid, y_dimid, xf_dimid, yf_dimid, dia_dimid, hora_dimid
 integer :: x_varid, y_varid, xf_varid, yf_varid, dia_varid, hora_varid
 integer :: CH1_in_varid, CH4_in_varid
 integer :: CH1_varid, Ch4_varid
 integer :: CH1_max_varid, CH1_min_varid, CH4_max_varid, CH4_min_varid
 
 integer :: start(NDIMS_IN), start_fine(NDIMS_IN), count(NDIMS_IN), count_fine(NDIMS_IN), start_hora(1)
 integer :: dimids(NDIMS), dimids_fine(NDIMS), dimids2d(NDIMS_IN), dimids2df(NDIMS_IN)
   
 allocate(CH1_out (9020, 1728))   ! Allocate memory for data.
 allocate(CH4_out (2255,432))     ! Lectura de variables de archivo original
 allocate(CH1_in (NXf,NYf))       ! imagenes invertidas
 allocate(CH4_in(NX,NY))
 
 allocate(CH1_max(NXf,NYf))
 allocate(CH1_min(NXf,NYf))
 allocate(CH4_max(NX,NY))
 allocate(CH4_min(NX,NY))
 CH1_max = 0
 CH1_min = 10000
 CH4_max = 0
 CH4_min = 10000 
!********************************************************** Fin declaracion Variables
 
 print *
 print *, '                        Recorte de imagenes NetCDF v 0.1'
 print *, '                        *********************************'
 print *

 call date_and_time(DATE=fecha, VALUES=tiempo)
 open (unit=16, file='log.txt')

 call get_command_argument(1, argument) ! Nombre de archivo lista.txt
 open (unit=8, file=trim(argument), status='old', ACTION='READ', IOSTAT=errorread)  
 If(errorread/=0) then
    write (*,*) " Error en apertura de archivo con lista de entrada!" 
    write (16,*) " Error en apertura de archivo con lista de entrada!"
    go to 999 
 End If
 
100 read (8,*, IOSTAT=errorread) filename
 if(errorread == -1) then
    write (16,*)
    write (*,*)
    write (16,*) " Terminado exitoso del procesamiento de imagenes"
    write (*,*) " Terminado exitoso del procesamiento de imagenes"
    write (*,*)
    write (16,*)
    call date_and_time(DATE=fecha, VALUES=tiempof)
    write (*,*) '   Tiempo procesamiento: ', tiempof (5) - tiempo (5) ,"horas", tiempof (6) - tiempo (6),"min.", &
					& tiempof (7) - tiempo (7), "sec."
	write (*,*) '    Archivos procesados: ' , rec	
	write (*,*) ' Archivos no procesados: ' , noct, '  (nocturnos)'
	write (*,*)
	write (16,*) 'Tiempo procesamiento: ',tiempof(5)-tiempo(5),'horas',tiempof(6)-tiempo(6),'min.',&
			& tiempof(7)-tiempo(7), 'sec.'
	write (16,*) ' Archivos procesados: ' , rec		
    close (16)
    
    call check( nf90_put_var(ncid, CH1_max_varid, CH1_max) )
    call check( nf90_put_var(ncid, CH1_min_varid, CH1_min) )
    call check( nf90_put_var(ncid, CH4_max_varid, CH4_max) )
    call check( nf90_put_var(ncid, CH4_min_varid, CH4_min) )
        
	call check( nf90_close(ncid) )
	
    go to 999
 End if

 if(errorread > 0) then
    write (16,*) " Error en lectura de nombre de archivo en archivo lista ", filename
    write (*,*) " Error en lectura de nombre de archivo en archivo lista ", filename
    call check( nf90_close(ncid) )
    call check( nf90_close(ncid_in) )
    go to 999
 end if
 
 string = trim(filename)
 ano=string(1:4)
 mes=string(5:6)
 cdia=string(7:8)
 hora=string(9:10)
 minu=string(11:12)
 foto=string(18:27)
 READ (ano,'(F4.0)') iano
 READ (mes,'(F2.0)') imes
 READ (cdia,'(F2.0)') idia
 READ (hora,'(F2.0)') ihora
 READ (minu,'(F2.0)') iminu
 ihora = ihora + iminu/60.

 If (foto=='media_hora') then ! no funciona
    ifoto = 1
 ElseIf (foto=='south_full') then
    ifoto = 2
 Else 
    ifoto = 0
 End If

 call diajuliano (idia, imes, iano, diaj)   ! entrada de reales en ves de enteros.
 call sunae(iano,diaj,ihora, latitud, longitud,az,el,ha,dec,soldst)  ! Comp.

 if (el < 7.0) then
    write (*,*) '   ',cdia,'  ', hora,':', minu, " eliminada, nocturna. "
    write (16,*) '                                   ', trim(filename),"   Eliminada, nocturna. "
    noct = noct +1
    goto 100
 end if

 If (imesp /= imes) then
    imesp = imes 
    !*********************************************************************** Definiciones NetCDF
    If (flag1) then
		call check( nf90_close(ncid) ) ! Cierra el archivo antes de abrir el nuevo.
		call date_and_time(VALUES=tiempof)
		write (*,*) '    Archivos procesados : ' , rec
		write (*,*) ' Archivos no procesados : ' , noct, '  (nocturnos)'
		write (*,*) '   Tiempo procesamiento : ',tiempo(5)-tiempof(5),'horas',tiempo(6)-tiempof(6),'min.',&
		& tiempo(7)-tiempof(7), 'sec.'
	End if	
    flag1 = .true.
    
    If (imes < 10 ) then
		filename_out = ano//'0'//mes//'.nc'
	Else
		filename_out = ano//mes//'.nc'
	End if
    write (*,*) filename_out
    write (16,*) filename_out
    write (16,*)
    call check( nf90_create(filename_out,NF90_64BIT_OFFSET, ncid) ) ! ncid, archivo de salida recortado
    
    ! Dimensiones
    call check( nf90_def_dim(ncid, "x", NX, x_dimid) )  ! Define the dimensions. 
    call check( nf90_def_dim(ncid, "y", NY, y_dimid) )
    call check( nf90_def_dim(ncid, "xf", NXf, xf_dimid) )  ! Define the dimensions. 
    call check( nf90_def_dim(ncid, "yf", NYf, yf_dimid) )
    call check( nf90_def_dim(ncid, "dia", dia,  dia_dimid) )
    call check( nf90_def_dim(ncid, "hora", NF90_UNLIMITED, hora_dimid) )
    ! Variables
    call check( nf90_def_var(ncid,"x", NF90_SHORT, x_dimid, x_varid) )	! Define the coordinate variables
    call check( nf90_def_var(ncid,"y", NF90_SHORT, y_dimid, y_varid) )	! 32 bit
    call check( nf90_def_var(ncid,"xf", NF90_SHORT, xf_dimid, xf_varid) )	
    call check( nf90_def_var(ncid,"yf", NF90_SHORT, yf_dimid, yf_varid) )
    call check( nf90_def_var(ncid,"hora", NF90_FLOAT, hora_dimid, hora_varid) ) ! 32 bit   

    dimids =  (/ x_dimid, y_dimid, hora_dimid  /) ! The dimids array is used to pass the IDs of the dimensions
    dimids_fine =  (/ xf_dimid, yf_dimid, hora_dimid /)
	dimids2d =  (/ x_dimid, y_dimid/)
	dimids2df =  (/ xf_dimid, yf_dimid/)
	
    call check( nf90_def_var(ncid, "CH1", NF90_SHORT, dimids_fine, CH1_varid) )  ! Define the variable and type. 
    call check( nf90_def_var(ncid, "CH4", NF90_SHORT, dimids, CH4_varid) )		!  NF90_Short (2-byte integer)
    call check( nf90_def_var(ncid, "CH1_max", NF90_SHORT, dimids2df, CH1_max_varid) )
    call check( nf90_def_var(ncid, "CH1_min", NF90_SHORT, dimids2df, CH1_min_varid) )
    call check( nf90_def_var(ncid, "CH4_max", NF90_SHORT, dimids2d, CH4_max_varid) )
    call check( nf90_def_var(ncid, "CH4_min", NF90_SHORT, dimids2d, CH4_min_varid) )

    ! Assign units attributes to coordinate variables.
    call check( nf90_put_att(ncid, x_varid, "units", "degrees_north"))
    call check( nf90_put_att(ncid, y_varid, "units", "degrees_east") )
     
    ! Atributos variables
    call check( nf90_put_att(ncid, CH1_varid, "units", "Albedo*100%") )
    call check( nf90_put_att(ncid, CH1_varid, "missing_value", -32768) )
    call check( nf90_put_att(ncid, CH1_varid, "valid_min", -32768) )
    call check( nf90_put_att(ncid, CH1_varid, "valid_max", 32768) )
    call check( nf90_put_att(ncid, CH1_varid, "scale_factor", 0.01) )

    call check( nf90_put_att(ncid, CH4_varid, "units", "temp_deg_C") )
    call check( nf90_put_att(ncid, CH4_varid, "missing_value", -32768) )
    call check( nf90_put_att(ncid, CH4_varid, "valid_min", -32768) )
    call check( nf90_put_att(ncid, CH4_varid, "valid_max", 32768) )
    call check( nf90_put_att(ncid, CH4_varid, "scale_factor", 0.01) )

    ! Atributos globales    
    call check( nf90_put_att(ncid, NF90_GLOBAL, "imagenes", "media_hora") )
    call check( nf90_put_att(ncid, NF90_GLOBAL, "projection_names", "orthographic") )
    call check( nf90_put_att(ncid, NF90_GLOBAL, "satellite", "goes-13") )
    call check( nf90_put_att(ncid, NF90_GLOBAL, "sensor", 12) )
    call check( nf90_put_att(ncid, NF90_GLOBAL, "mes", imes) )  
    call check( nf90_put_att(ncid, NF90_GLOBAL, "ano", iano) )  
    call check( nf90_put_att(ncid, NF90_GLOBAL, "NX", NX) )  
    call check( nf90_put_att(ncid, NF90_GLOBAL, "NY", NY) )  
    call check( nf90_put_att(ncid, NF90_GLOBAL, "NXf", NXf) )
    call check( nf90_put_att(ncid, NF90_GLOBAL, "NYf", NYf) )
    call check( nf90_put_att(ncid, NF90_GLOBAL, "Procesadox", "JPJ") )
    call check( nf90_put_att(ncid, NF90_GLOBAL, "date", fecha) )

    ! End define mode.
    call check( nf90_enddef(ncid) )
    !*********************************************************************** Fin Definiciones NetCDF	
   
   x(1) = 65.
   y(1) = 20. 
   xf(1) = x(1)
   yf(1) = y(1)

    do i = 2, NX 	! Guardar coordenadas de lat y lon.
        x(i) =  x(i-1) + 0.02
    end do
    do j = 2, NY
        y(j) = y(j-1) + 0.13
    end do	
    do i = 2, NXf
        xf(i) = xf(i-1) + 0.005
    end do
    do j = 2, NYf
        yf(j) = yf(j-1) + 0.0325
    end do	

    ! Write the coordinate variable data. This will put the latitudes
    ! and longitudes of our data grid into the netCDF file.
    call check( nf90_put_var(ncid, x_varid, x) )
    call check( nf90_put_var(ncid, y_varid, y) )	
    call check( nf90_put_var(ncid, xf_varid, xf) )
    call check( nf90_put_var(ncid, yf_varid, yf) )

    rec = 0
    noct= 0
    
 end if
 
 rec = rec + 1 
 start_hora =(/1/)
 start_hora(1) = rec
 ihorat = (ihora + idia*24)  ! Hora (hr_mes*100)
 call check( nf90_put_var(ncid, hora_varid, ihorat, start_hora)  )  ! Graba hora
 call check( nf90_open(trim(filename), nf90_nowrite, ncid_in) )  ! Abre archivo de lectura, ncid_in
 
 call check( nf90_inq_varid(ncid_in, "gvar_ch1_fine", CH1_in_varid) )
 call check( nf90_inq_varid(ncid_in, "gvar_ch4", CH4_in_varid) )

		
 call check( nf90_get_var(ncid_in, CH4_in_varid, CH4_out)) 
 call check( nf90_get_var(ncid_in, CH1_in_varid, CH1_out))  
 
 ! Recorte de imagenes
 CH4_in = CH4_out(925:1450,:)
 CH1_in = CH1_out(3700:5800,:)
 
 ! Correccion de datos 
 do i = 1, NXf
    do j = 1, NYf
		if (CH1_in(i,j) > CH1_max(i,j)) then
			CH1_max(i,j) = CH1_in(i,j)
		end if
		if (CH1_in(i,j) < CH1_min(i,j)) then
			CH1_min(i,j) = CH1_in(i,j)
		end if
		if (CH1_in(i,j) > 15000) then
            write (16,*) "CH1_in(", i, ",",j, ") =", CH1_in(i, j) , filename, iano,diaj,ihora,"Corregido, vecinos"
   			CH1_in(i,j) = (CH1_in(i+1,j)+CH1_in(i,j+1)+CH1_in(i-1,j)+CH1_in(i,j-1))/4.
        else if (CH1_in(i,j) > 10000) then
            write (16,*) "CH1_in(", i, ",",j, ") =", CH1_in(i, j) , filename, iano,diaj,ihora,"Corregido"
   			CH1_in(i,j) = 10000
        end if
        if (CH1_in(i, j) < 0 ) then
			write (16,*) "CH1_in(", i, ",", j, ") =", CH1_in(i, j), filename, iano,diaj,ihora,"Pixel malo"
			
			If (CH1_in(i+1,j) > 0 .and. CH1_in(i,j+1) > 0 .and. CH1_in(i-1,j) > 0 .and. CH1_in(i,j-1) > 0) Then
				CH1_in(i,j) = (CH1_in(i+1,j)+CH1_in(i,j+1)+CH1_in(i-1,j)+CH1_in(i,j-1))/4.        
			Else
				CH1_in(i,j) = (CH1_in(i+1,j)+CH1_in(i-1,j))/2.
			End If
        end if
    end do
end do
 
 do i = 1, NX
    do j = 1, NY
		if (CH4_in(i,j) > CH4_max(i,j)) then
			CH4_max(i,j) = CH4_in(i,j)
		end if
		if (CH4_in(i,j) < CH4_min(i,j)) then
			CH4_min(i,j) = CH4_in(i,j)
		end if    
		if (CH4_in(i,j) > 15000) then
            write (16,*) "CH4_in(", i, ",",j, ") =", CH4_in(i, j) , filename, iano,diaj,ihora,"Corregido, vecinos"
   			CH4_in(i,j) = (CH4_in(i+1,j)+CH4_in(i,j+1)+CH4_in(i-1,j)+CH4_in(i,j-1))/4.
        else if (CH4_in(i, j) > 10000) then
            write (16,*) "CH4_in(", i, ",", j, ") =", CH4_in(i, j), filename, iano,diaj,ihora,"Corregido"
   			CH4_in(i,j) = 10000
        end if
        if (CH4_in(j, i) < -10000 ) then
			write (16,*) "CH4_in(", i, ",", j, ") =", CH4_in(i, j), filename, iano,diaj,ihora,"Pixel malo"
			If (CH4_in(i+1,j) > -10000 .and. CH4_in(i,j+1) > -10000 .and. CH4_in(i-1,j) > -10000 .and. CH4_in(i,j-1)> -10000) Then
				CH4_in(i,j) = (CH4_in(i+1,j)+CH4_in(i,j+1)+CH4_in(i-1,j)+CH4_in(i,j-1))/4.        
			Else
				CH4_in(i,j) = (CH4_in(i+1,j)+CH4_in(i-1,j))/2.
			End If
        end if
    end do
 end do
 
 ! Guardado de matriz recortada y revisada
 call check( nf90_put_var(ncid,CH4_varid, CH4_in))
 call check( nf90_put_var(ncid,CH1_varid, CH1_in))
 call check( nf90_close(ncid_in) )
 
 write (*,*) '   ',cdia,'  ', hora,':', minu
 
 Go to 100
 
 999 Stop
 
 contains
 subroutine check(status)
  integer, intent ( in) :: status
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop 2
  end if
 end subroutine check
 
 
 End Program ProcesamientoImagenes

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
!This program calculates the day of year corresponding to a specified date.

IMPLICIT NONE
! Data dictionary: declare variable types, definitions, & units
real, INTENT(IN):: day          !Day (dd)
real, INTENT(IN) :: month        !Month (mm)
real, INTENT(IN) :: year         !Year (yyyy)
real, INTENT(out) :: dayj 		!Day of year
INTEGER :: i            			!Index,variable
INTEGER :: leap_day     			!Extra day for leap year

! Check for leap year, and add extra day if necessary
IF ( MOD(year,400.) == 0 ) THEN
    leap_day = 1    ! Years divisible by 400 are leap years
ELSE IF ( MOD(year,100.) == 0 ) THEN
    leap_day = 0    ! Other centuries are not leap years
ELSE IF ( MOD(year,4.) == 0 ) THEN
    leap_day = 1    ! Otherwise every 4th year 1S a leap year
ELSE
    leap_day = 0    ! Other years are not leap years
END IF


! Calculate day of year
dayj= day
DO i = 1, nint(month)-1
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


 SUBROUTINE sunae(year,day,hour, lat, long,az,el,ha,dec,soldst)
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
!      az=asin(-cos(dec)*sin(ha)/cos(el))

!!   this puts azimuth between 0 and 2*pi radians
!      if(sin(dec)-sin(el)*sin(lat).ge.0.) then
!		if(sin(az).lt.0.) az=az+twopi
!      else
!      az=pi-az
!      endif
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
!!
!!   calculate distance to sun in A.U. & diameter in degs
!      soldst=1.00014-.01671*cos(mnanom)-.00014*cos(2.*mnanom)
!      soldia=.5332/soldst

!!   convert az and lat to degs before returning
!      az=az/rad
!      lat=lat/rad
!	 ha=ha/rad
!	 dec=dec/rad

!!   mnlong in degs, gmst in hours, jd in days if 2.4e6 added;
!!   mnanom,eclong,oblqec,ra,and lmst in radians
 End subroutine
