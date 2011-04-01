program extraccion

 use netcdf

 implicit none
 
 ! Variables Programa
 integer, parameter :: Nimagenesdiamin = 5
 integer, parameter :: Nimagenesdiamax = 26 
 integer :: horaf, horaini, horasdia, diasmes
 integer :: largofilein, lat, lon, i , j, k
 integer :: ano, mes
 character(8)  :: fecha
 character(4)  :: cano
 character(2)  :: cmes
 character (40) :: argument 
 character(25) :: filename_out
 integer :: Nhora, NXf, NYf
 integer,dimension(8) :: tiempo, tiempof, tiempoa
 integer :: TIDmax, diasmes, horasdia, horaini, Globfact, Dirfact, Latfact, Lonfact
 
 integer, dimension(:,:), allocatable :: Nimagen, Lat_CH1, Lon_CH1
 integer, dimension(:,:), allocatable :: Global, Directa, Difusa, DNI
 real, dimension (:), allocatable :: hora
 
 
!********************************************************** Fin declaracion Variables
 
 
 print *
 print *, '                 Extraccion de datos de radiacion Global, Difusa y Directa'
 print *, '                 *********************************************************'
 print *

!**************************************************** Lectura de archivo de entrada

 call date_and_time(DATE=fecha, VALUES=tiempo)
 tiempoa = tiempo
 
 open (unit=16, file='log_extraccion.txt')
  
 call get_command_argument(1, argument)
 largofilein=len_trim(argument)
 If ( argument(largofilein-6: largofilein) /= '.prom.nc') then
	Write (*,*) '  Archivo de entrada no es base de datos de radiacion, Continuar ??', argument(largofilein-6: largofilein)
	Read  (*,*)
 End if
 
 call check( nf90_open(argument, nf90_nowrite, ncid) )

 call check( nf90_inq_varid(ncid, "hora", hora_varid) )
 call check( nf90_inq_varid(ncid, "Lat_CH1", Lat_CH1_varid) )
 call check( nf90_inq_varid(ncid, "Lon_CH1", Lon_CH1_varid) )
 call check( nf90_inq_varid(ncid, "Global_hor", Global_varid) )
 call check( nf90_inq_varid(ncid, "Directa_hor", Directa_varid) )

 call check( nf90_inq_dimid(ncid, "xf", xf_dimid) )
 call check( nf90_inq_dimid(ncid, "yf", yf_dimid) )
 call check( nf90_inq_dimid(ncid, "hora", hora_dimid) )
 call check( nf90_inq_dimid(ncid, "dia", dia_dimid) )
 
 call check( nf90_inquire_dimension(ncid, xf_dimid, len = NXf) )
 call check( nf90_inquire_dimension(ncid, yf_dimid, len = NYf) )
 call check( nf90_inquire_dimension(ncid, hora_dimid, len = Nhora) )
 call check( nf90_inquire_dimension(ncid, dia_dimid, len = Ndias) )

 call check( nf90_get_att(ncid, NF90_GLOBAL, "ano", ano))
 call check( nf90_get_att(ncid, NF90_GLOBAL, "mes", mes))
 call check( nf90_get_att(ncid, NF90_GLOBAL, "horaini", horaini))
 call check( nf90_get_att(ncid, NF90_GLOBAL, "horasdia", horasdia))
 call check( nf90_get_att(ncid, NF90_GLOBAL, "diasmes", diasmes))
 
 call check( nf90_get_att(ncid, Global_varid, "scale_factor", Globfact))
 call check( nf90_get_att(ncid, Directa_varid, "scale_factor", Dirfact))
 call check( nf90_get_att(ncid, Lat_CH1_varid, "scale_factor", Latfact)) 
 call check( nf90_get_att(ncid, Lon_CH1_varid, "scale_factor", Lonfact)) 
 
 ! allocate variables
 allocate (hora(Nhora))
 allocate (Lat_CH1 (NXf,NYf))       
 allocate (Lon_CH1 (NXf,NYf))
 allocate (Global ( Nhora, Ndias))
 allocate (Directa ( Nhora, Ndias))
 allocate (Difusa ( Nhora, Ndias))
 allocate (DNI  ( Nhora, Ndias))
 
 
 
 call check( nf90_get_var(ncid, hora_varid, hora) )
 call check( nf90_get_var(ncid, Lat_CH1_varid, Lat_CH1) )
 call check( nf90_get_var(ncid, Lon_CH1_varid, Lon_CH1) )
 
 Write (*,*) 'Latitud y longitud de punto a extraer datos: (grados decimales) '
 read (*,*) Lat, Lon
 
 iLat = nint(-100. * abs(Lat)); iLon = nint(-100. * abs(Lon)) ! -100 = -1/Latfact
 i=0; j=0
 Do 
	i=i+1
	If (Lat_CH1(i,1)<iLat) Then
		If (abs(Lat_CH1(i,1)-iLat)>abs(Lat_CH1(i-1,1)-iLat) i= i-1
		Do 
			j = j+1
			If (Lon_CH1(i,j)<iLon) Then
				If (abs(Lon_CH1(i,j)-iLon)>abs(Lon_CH1(i,j-1)-iLon) j= j-1
				exit
			end if
		 If (j> NXf) then
			write (*,*) ' Problemas al tratar de ajustar la Longitud!!' , i,j
			exit
		end if 	
		end do
		
	end if
 If (i> NYf) then
	write (*,*) ' Problemas al tratar de ajustar la Latitud!!' , i,j
	exit
 end if 	
 end do
 
 start = (/ i, j, 1, 1 /)
 count = (/ 1, 1, Nhora, Ndias /)
 ! dimids4d =  (/ xf_dimid_prom, yf_dimid_prom, hora_dimid_prom, diasmes_dimid_prom /) 
 call check( nf90_get_var (ncid, Global_varid, Global, start, count) )
 call check( nf90_get_var (ncid, Directa_varid, Directa, start, count) )
 
 open (unit=15, file='datos.txt')
 
 write (15,*)
 write (15,*) '   Dia	Hora		Global		Difusa		Directa		DNI'
 write (15,200)  (











200 Format ('   ',I2,'	'

 contains
 subroutine check(status)
  integer, intent ( in) :: status
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop 2
  end if
 end subroutine check
 
end program extraccion
