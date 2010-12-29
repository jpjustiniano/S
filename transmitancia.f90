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
 
 integer i, j
 
 ! Variables NETCDF 
 integer :: ncid, ncid_in, dia =31, status
 integer :: NDIMS = 3, NDIMS_IN	= 2     
 integer :: NX, NY, NXf, NYf 
 real, dimension(:), allocatable :: x, y, xf, yf, hora
 integer, dimension(:,:), allocatable :: CH1_out, CH4_out   ! Matrices de archivo original
 integer, dimension(:,:), allocatable :: CH4_in , CH1_in        !  Matrices de archivo recortado
 integer, dimension(:,:), allocatable :: CH1_max, CH1_min, CH4_max,CH4_min  !  Matrices max min
 integer :: x_dimid, y_dimid, xf_dimid, yf_dimid, dia_dimid, hora_dimid
 integer :: x_varid, y_varid, xf_varid, yf_varid, dia_varid, hora_varid
 integer :: CH1_in_varid, CH4_in_varid
 integer :: CH1_varid, Ch4_varid
 integer :: CH1_max_varid, CH1_min_varid, CH4_max_varid, CH4_min_varid
 
 integer :: start(NDIMS_IN), start_fine(NDIMS_IN), count(NDIMS_IN), count_fine(NDIMS_IN), start_hora(1)
 integer :: dimids(NDIMS), dimids_fine(NDIMS), dimids2d(NDIMS_IN)
   
 
!********************************************************** Fin declaracion Variables
 
 
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
 
 ! allocate variables
 allocate(x(NX))
 allocate(y(NY))
 allocate(xf(NXf))
 allocate(yf(NYf))
 allocate(hora(Nhora))
 
 allocate(CH1_in (NXf,NYf))       ! imagenes invertidas
 allocate(CH4_in(NX,NY))
 
 allocate(CH1_max(NXf,NYf))
 allocate(CH1_min(NXf,NYf))
 allocate(CH4_max(NX,NY))
 allocate(CH4_min(NX,NY))
 

 ! Lectura de matrices de maximo y minimo
 call check( nf90_get_var(ncid, x_varid, x) )
 call check( nf90_get_var(ncid, y_varid, y) )
 call check( nf90_get_var(ncid, CH1_max_varid, CH1_max) )
 call check( nf90_get_var(ncid, CH1_min_varid, CH1_min) )
 call check( nf90_get_var(ncid, CH4_max_varid, CH4_max) )
 call check( nf90_get_var(ncid, CH4_min_varid, CH4_min) )

  
 ! Lectura de pixeles de CH1 y CH4
 count = (/ NX, NY, 1 /)
 countf = (/ NXf, NYf, 1 /)
 start = (/ 1, 1, 1 /)
 call check( nf90_get_var(ncid, hora_varid, hora))


 do rec = 1, Nhora
	dia = hora(rec)/24
	horad = hora(rec)- dia*24
	start(3) = rec

	call check( nf90_get_var(ncid, CH1_varid, CH1_in, start = start, &
						   count = countf) )						   						   
	call check( nf90_get_var(ncid, CH4_varid, CH4_in, start, count) )
	
	Do i = 1, NXf
		Do j=1, NYf
			! Lectura de archivo con valores de temperatura, HR, vis., albedo y altura, lat y long.
			! En funcion del archivo de altura se procesa o no el pixel.
			 
			Temp = 300.
			Alt = 500.
			HR = 40.
			albedo = 20.
			vis = 10.
			!???? Esta lloviendo ??? negativo ??
			! IF(XUMI(I,J)*100.0.GT.-98.0) THEN 
			!   ...... proceso de matriz
			!ELSE
			!	XTCR=-1.0
			!	XTCD=-1.0
			!	XTDIR=-1.0
			!ENDIF
			
			! test for precitable water less than 0.0000001 cm
			! IF (WH2O .LE.0.0000001) WH2O  = 0.0000001
			
			! calculation of the visibility at the station as function of visibility and altitude
			 VZGH =  100.0
			 ZGH  = 1000.0
			 VIS   = vis * EXP( (LOG(VZGH/vis)/ZGH)*Alt )
			! test for visibility between 2 and 150 km
			 IF (VIS .GT. 150.)  VIS = 150.
			 IF (VIS .LT.   2.)  VIS =   2.
			 
			! change relative humidity and albedo from percentual to relative values
			albedo = albedo/1000.0
			HR     = HR/100.0
			! calculation of the precitable water as function of relative humidity, 
			! partial pressure of water vapor and temperature
			! PVSAT   = EXP(26.23 - 5416.0/Temp)
			! WH2O  = 0.493*hr*PVSAT/Temp
			
			CALL D2STR(x(i),y(j),Temp,HR,Alt,albedo,vis,JDAY,XHOR,0, 0.0D0,
						&  TCLEAR(i,j),TDIR(i,j))
			
			CALL D2STR(x(i),y(j),Temp,rf,zstat,SFCALB,VIS,JDAY,XHOR,1,100.0D0,
						&  TCLOUD(i,j),TDIR)
	
	
 end do







	
	
