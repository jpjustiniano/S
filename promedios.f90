! Copyright (C) 2011, Juan Pablo Justiniano  <jpjustiniano@gmail.com>
! Programa para el calculo de la radiacion Global, Difusa y Directa.
! Revisar:


Program ProcesamientoImagenes_media_hora
 use netcdf
 !$ use OMP_LIb
 implicit none
 
 ! Variables Programa
 
 
 
 ! Variables NETCDF archivo entrada


 ! Variables NETCDF archivo salida


!********************************************************** Fin declaracion Variables
 
 
 print *
 print *, '                  Calculo de la radiacion Global, Difusa y Directa'
 print *, '                  ************************************************'
 print *

!**************************************************** Lectura de archivo de entrada

 TID = 1
 !$    TID = omp_get_num_procs();
 !$		If (TID>6) TID = 6
 !$    call OMP_SET_NUM_THREADS(TID)
 !$ print *, ' Numero de procesadores en uso: ',TID
 call date_and_time(DATE=fecha, VALUES=tiempo)
 open (unit=16, file='log_prom.txt')
 
 tiempoa = tiempo
 
 call get_command_argument(1, argument)
 call check( nf90_open(argument, nf90_nowrite, ncid) )
      
 ! Get the varids of the latitude and longitude coordinate variables.
 call check( nf90_inq_varid(ncid, "hora", hora_varid) )

 call check( nf90_inq_varid(ncid, "CH1", CH1_varid) )
 call check( nf90_inq_varid(ncid, "CH4", CH4_varid) )
 call check( nf90_inq_varid(ncid, "Lat_CH4", Lat_CH4_varid) )
 call check( nf90_inq_varid(ncid, "Lon_CH4", Lon_CH4_varid) )
 call check( nf90_inq_varid(ncid, "Lat_CH1", Lat_CH1_varid) )
 call check( nf90_inq_varid(ncid, "Lon_CH1", Lon_CH1_varid) )
 call check( nf90_inq_varid(ncid, "Global", Global_varid) )
 call check( nf90_inq_varid(ncid, "Directa", Directa_varid) )

 call check( nf90_inq_dimid(ncid, "x", x_dimid) )
 call check( nf90_inq_dimid(ncid, "y", y_dimid) )
 call check( nf90_inq_dimid(ncid, "xf", xf_dimid) )
 call check( nf90_inq_dimid(ncid, "yf", yf_dimid) )
 call check( nf90_inq_dimid(ncid, "hora", hora_dimid) )
 
 call check( nf90_inquire_dimension(ncid, x_dimid, len = NX) )
 call check( nf90_inquire_dimension(ncid, y_dimid, len = NY) )
 call check( nf90_inquire_dimension(ncid, xf_dimid, len = NXf) )
 call check( nf90_inquire_dimension(ncid, yf_dimid, len = NYf) )
 call check( nf90_inquire_dimension(ncid, hora_dimid, len = Nhora) )

 call check( nf90_get_att(ncid, NF90_GLOBAL, "ano", ano))
 call check( nf90_get_att(ncid, NF90_GLOBAL, "mes", mes))
 
 ! allocate variables
 allocate(x(NX))
 allocate(y(NY))
 allocate(xf(NXf))
 allocate(yf(NYf))
 allocate(hora(Nhora))
 
 allocate(Horario (NXf,NYf,NhoraMax)) 
 allocate(Diario (NXf,NYf)) 
 allocate(Mensual (NXf,NYf)) 
 
 
 count = (/ NX, NY, 1 /)
 countf = (/ NXf, NYf, 1 /)
 start = (/ 1, 1, 1 /)
 startm = (/ 1, 1, 1 /)
 call check( nf90_get_var(ncid, hora_varid, hora))
 
 
 write (cano,'(I4)') NInt(ano)
 If (mes < 10 ) then
	write(cmes,'(I1)') NInt(mes)
	filename_out = cano//'0'//cmes//'.media_hora.prom.nc'
 Else
	write(cmes,'(I2)') NInt(mes)
	filename_out = cano//cmes//'.media_hora.prom.nc'
 End if
 
 call check( nf90_create(filename_out,NF90_64BIT_OFFSET, ncid_prom) ) ! ncid, archivo de salida recortado
 ! Dimensiones
 call check( nf90_def_dim(ncid_rad, "x", NX, x_dimid_prom) )  ! Define the dimensions. 
 call check( nf90_def_dim(ncid_rad, "y", NY, y_dimid_prom) )
 call check( nf90_def_dim(ncid_rad, "xf", NXf, xf_dimid_prom) )  ! Define the dimensions. 
 call check( nf90_def_dim(ncid_rad, "yf", NYf, yf_dimid_prom) )
 call check( nf90_def_dim(ncid_rad, "hora", Nhora, hora_dimid_prom) )
 
 ! Variables
 call check( nf90_def_var(ncid_rad,"hora", NF90_FLOAT, hora_dimid_prom, hora_varid_prom) ) ! 32 bit
     

 dimids =  (/ xf_dimid_prom, yf_dimid_prom, hora_dimid_prom /) ! The dimids array is used to pass the IDs of the dimensions
 dimids2d =  (/ x_dimid, y_dimid/)
 dimids2df =  (/ xf_dimid, yf_dimid/)

 call check( nf90_def_var(ncid_rad, "1", NF90_SHORT, dimids, h1_varid) )
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
