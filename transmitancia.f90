! Copyright (C) 2010, Juan Pablo Justiniano  <jpjustiniano@gmail.com>
! Programa para el calculo de la radiacion Global, Difusa y Directa.

! gfortran -c -I/usr/include "%f"
! gfortran -I/usr/include -L/usr/lib -lnetcdff -lnetcdf -o "%e" "%f"
! gfortran -o "%e" "%f" -L/usr/lib  -lnetcdf   ?? de manual...

!Canal 1: [1728, 9020]
!Canal 4: [432, 2255]

!Imagen Media Hora:
!Canal 1_SI: [1, 3700]
!Canal 1_ID: [1728, 5800]
 
!Canal 4_SI: [1, 925]
!Canal 4_ID: [432, 1450]

 Program trasmitancia
 use netcdf
 implicit none
 
 print *
 print *, 'Calculo de la radiacion Global, Difusa y Directa'
 print *

 call date_and_time(DATE=fecha, VALUES=tiempo)
 open (unit=16, file='log_trans.txt')

 call get_command_argument(1, argument)
 call check( nf90_open(argument, nf90_nowrite, ncid) )

      
 ! Get the varids of the latitude and longitude coordinate variables.
 call check( nf90_inq_varid(ncid, "x", x_varid) )
 call check( nf90_inq_varid(ncid, "y", y_varid) )
 call check( nf90_inq_varid(ncid, "xf", xf_varid) )
 call check( nf90_inq_varid(ncid, "yf", yf_varid) )
 call check( nf90_inq_varid(ncid, "hora", hora_varid) )

 call check( nf90_inq_varid(ncid, "CH1", CH1_varid) )
 call check( nf90_inq_varid(ncid, "CH4", CH4_varid) )
 call check( nf90_inq_varid(ncid, "CH1_max", CH1_max_varid) )
 call check( nf90_inq_varid(ncid, "CH1_min", CH1_min_varid) )
 call check( nf90_inq_varid(ncid, "CH4_max", CH4_max_varid) )
 call check( nf90_inq_varid(ncid, "CH4_min", CH4_min_varid) )

 ! Read the latitude and longitude data.
 call check( nf90_get_var(ncid, lat_varid, lats) )
 call check( nf90_get_var(ncid, lon_varid, lons) )

 ! Read 1 record of NLONS*NLATS*NLVLS values, starting at the beginning
 ! of the record (the (1, 1, 1, rec) element in the netCDF file).
 count = (/ NLONS, NLATS, NLVLS, 1 /)
 start = (/ 1, 1, 1, 1 /)

 ! Read the surface pressure and temperature data from the file, one
 ! record at a time.
 do rec = 1, NRECS
	start(4) = rec
	call check( nf90_get_var(ncid, pres_varid, pres_in, start = start, &
						   count = count) )
	call check( nf90_get_var(ncid, temp_varid, temp_in, start, count) )
