! Copyright (C) 2010, Juan Pablo Justiniano  <jpjustiniano@gmail.com>
! Programa para el calculo de la radiacion Global, Difusa y Directa.

! gfortran -c -I/usr/include "%f"
! gfortran -I/usr/include -L/usr/lib -lnetcdff -lnetcdf -o transmitancia transmitancia.f90
! gfortran -g -I/usr/include -L/usr/lib -lnetcdff -lnetcdf -o transmitancia transmitancia.f90 Atmosphe.o Strpsrb.o Astro.o

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
 
 integer i, j, rec
 real ano, mes, diaj, dia 
 Real Alt, HR, Albedo, temp, vis
 Real horad
 character (30) :: argument  ! Nombre de archivo 
 character(8)  :: fecha
 integer,dimension(8) :: tiempo, tiempof
 real, parameter :: pi= acos(-1.0), cdr= pi/180.0
 
 ! Variables NETCDF archivo entrada
 integer :: ncid, ncid_in, status
 integer, parameter :: NDIMS = 3, NDIMS_IN	= 2     
 integer :: NX, NY, NXf, NYf, Nhora
 real, dimension(:), allocatable :: x, y, xf, yf, hora
 integer, dimension(:,:), allocatable :: CH1_out, CH4_out   ! Matrices de archivo original
 integer, dimension(:,:), allocatable :: CH4_in , CH1_in        !  Matrices de archivo recortado
 integer, dimension(:,:), allocatable :: CH1_max, CH1_min, CH4_max,CH4_min  !  Matrices max min
 real, dimension(:,:), allocatable :: TCLEAR,TDIR, TCLOUD
 integer :: x_dimid, y_dimid, xf_dimid, yf_dimid, dia_dimid, hora_dimid
 integer :: x_varid, y_varid, xf_varid, yf_varid, dia_varid, hora_varid
 integer :: CH1_in_varid, CH4_in_varid
 integer :: CH1_varid, Ch4_varid
 integer :: CH1_max_varid, CH1_min_varid, CH4_max_varid, CH4_min_varid
 integer :: TDIRn
 
 integer :: start(NDIMS), count(NDIMS), countf(NDIMS)
 integer :: dimids(NDIMS), dimids_fine(NDIMS), dimids2d(NDIMS_IN)
 
 ! Variables NETCDF archivo salida
 character(4)  :: cano
 character(2)  :: cmes
 character(14) :: filename_out
 integer :: ncid_rad
 integer :: x_dimid_rad, y_dimid_rad,hora_dimid_rad
 integer :: x_varid_rad, y_varid_rad, hora_varid_rad
 integer :: Global_varid, Difusa_varid, Directa_varid

 
 ! Variables Brasil
 Real :: E0,DEC,ET
 Real :: DECR, YLATR, CODEC, COLAT, SIDEC, SILAT
 Real :: ZN, TIMCOR, ROFF, IWP, ISUB, TOP, NCL, INTVAL, CLOLWC
 Real :: LATMOS,TS, TSOLAR, WSOLAR, COWI, COSZEN, THETA
   
!********************************************************** Fin declaracion Variables
 
 
 print *
 print *, '                  Calculo de la radiacion Global, Difusa y Directa'
 print *, '                  ************************************************'
 print *

!**************************************************** Lectura de archivo de entrada
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

 call check( nf90_get_att(ncid, NF90_GLOBAL, "ano", ano))
 call check( nf90_get_att(ncid, NF90_GLOBAL, "mes", mes))
   
 
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
 
 allocate(TCLEAR(NX,NY))
 allocate(TDIR(NX,NY))
 allocate(TCLOUD(NX,NY))

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
 
 !**************************************************** Creacion de archivo de salida
 
 write (cano,'(F4.0)') ano
 If (mes < 10 ) then
	write(cmes,'(F1.0)') mes
	filename_out = cano//'0'//cmes//'_rad.nc'
 Else
	write(cmes,'(F2.0)') mes
	filename_out = cano//cmes//'_rad.nc'
 End if
 
 call check( nf90_create(filename_out,NF90_64BIT_OFFSET, ncid_rad) ) ! ncid, archivo de salida recortado
 ! Dimensiones
 call check( nf90_def_dim(ncid_rad, "x", NX, x_dimid_rad) )  ! Define the dimensions. 
 call check( nf90_def_dim(ncid_rad, "y", NY, y_dimid_rad) )
 call check( nf90_def_dim(ncid_rad, "hora", Nhora, hora_dimid_rad) )
 ! Variables
 call check( nf90_def_var(ncid_rad,"x", NF90_SHORT, x_dimid_rad, x_varid_rad) )	! Define the coordinate variables
 call check( nf90_def_var(ncid_rad,"y", NF90_SHORT, y_dimid_rad, y_varid_rad) )	! 32 bit
 call check( nf90_def_var(ncid_rad,"hora", NF90_SHORT, hora_dimid_rad, hora_varid_rad) ) ! 32 bit    

 dimids =  (/ x_dimid_rad, y_dimid_rad, hora_dimid_rad /) ! The dimids array is used to pass the IDs of the dimensions
 dimids2d =  (/ x_dimid_rad, y_dimid_rad/)

 call check( nf90_def_var(ncid_rad, "Global", NF90_SHORT, dimids, Global_varid) )
 call check( nf90_def_var(ncid_rad, "Difusa", NF90_SHORT, dimids, Difusa_varid) )
 call check( nf90_def_var(ncid_rad, "Directa", NF90_SHORT, dimids, Directa_varid) )
 
 ! Atributos variables
 call check( nf90_put_att(ncid_rad, Global_varid, "units", "W/m2") )
 call check( nf90_put_att(ncid_rad, Global_varid, "missing_value", -32768) )
 call check( nf90_put_att(ncid_rad, Global_varid, "valid_min", -32768) )
 call check( nf90_put_att(ncid_rad, Global_varid, "valid_max", 32768) )
 call check( nf90_put_att(ncid_rad, Global_varid, "scale_factor", 0.01) )
 
 call check( nf90_put_att(ncid_rad, Difusa_varid, "units", "W/m2") )
 call check( nf90_put_att(ncid_rad, Difusa_varid, "missing_value", -32768) )
 call check( nf90_put_att(ncid_rad, Difusa_varid, "valid_min", -32768) )
 call check( nf90_put_att(ncid_rad, Difusa_varid, "valid_max", 32768) )
 call check( nf90_put_att(ncid_rad, Difusa_varid, "scale_factor", 0.01) )
 
 call check( nf90_put_att(ncid_rad, Directa_varid, "units", "W/m2") )
 call check( nf90_put_att(ncid_rad, Directa_varid, "missing_value", -32768) )
 call check( nf90_put_att(ncid_rad, Directa_varid, "valid_min", -32768) )
 call check( nf90_put_att(ncid_rad, Directa_varid, "valid_max", 32768) )
 call check( nf90_put_att(ncid_rad, Directa_varid, "scale_factor", 0.01) )
 
 ! End define mode.
 call check( nf90_enddef(ncid_rad) )
 !**************************************************** Procesamiento de las imagenes Sat.
 
 do rec = 1, Nhora
	dia = int(hora(rec)/24)
	horad = hora(rec)- dia*24
	call diajuliano (dia, mes, ano, diaj) 
	
	
	

	start(3) = rec
	

	call check( nf90_get_var(ncid, CH1_varid, CH1_in, start = start, &
						   count = countf) )						   						   
	call check( nf90_get_var(ncid, CH4_varid, CH4_in, start, count) )
	
	Do i = 1, NXf
		Do j=1, NYf
			! Lectura de archivo con valores de temperatura, HR, vis., albedo y altura, lat y long.
			! En funcion del archivo de altura se procesa o no el pixel.
			 
			Temp = 300. ! (K)
			Alt = 500. 	! (m)
			HR = .4 	! (0-1. )
			albedo = .20! (0-1.)
			vis = 10.	! (km)
			
			IF(Alt > 0. ) THEN 
		
			! calculation of the visibility at the station as function of visibility and altitude
			 VIS   = vis * EXP( (LOG(100.0/vis)/1000.0)*Alt )
			! test for visibility between 2 and 150 km
			 IF (VIS .GT. 150.)  VIS = 150.
			 IF (VIS .LT.   2.)  VIS =   2.
			 
			! change relative humidity and albedo from percentual to relative values
			!albedo = albedo/1000.0
			!HR     = HR/100.0
			! calculation of the precitable water as function of relative humidity, 
			! partial pressure of water vapor and temperature 
			! PVSAT   = EXP(26.23 - 5416.0/Temp)
			! WH2O  = 0.493*hr*PVSAT/Temp
			
			! subroutine ASTRO - calculation of eccentricity correction, declination an equation of time
			! input: JDAY ; output: E0,DEC,ET
			CALL ASTRO(diaj,E0,DEC,ET) 
			DECR = DEC*cdr
			YLATR= y(j)*cdr
			CODEC = COS(DECR)
			COLAT = COS(YLATR)
			SIDEC = SIN(DECR)
			SILAT = SIN(YLATR)
			
			!calculate time correction
			ZN = 0.0
			TIMCOR = (4.0*(15.0*ZN+x(i))+ET)/60.0
			
			! Fixed input parameters for subroutine strpsrb
			ROFF   = 0.0
			IWP    = 3
			ISUB   = 2
			TOP    = 500.0
			NCL    = 2
			INTVAL = 3
			CLOLWC = 0.0
			
			! subroutine ATMOSPHERE - choose atmosphere by surface temperature
			! input: Temp ; output: LATMOS,TS
			CALL ATMOSPHERE(Temp,LATMOS,TS)
			
			TSOLAR = horad+TIMCOR    !para entrada com horario em UTC
			
			! Characteristic hour angle WI corresponding to the W1-W2 interval
			WSOLAR    = (12.00 - TSOLAR)*15.
			COWI  = COS(WSOLAR*CDR)
			
			! Characteristic solar zenith angle THETA (o)
			COSZEN = SIDEC*SILAT + CODEC*COLAT*COWI
			THETA  = ACOS(COSZEN)/CDR
			
			! Subroutine STRPSRB - calculate transmittance for clear sky
			! input  - LATMOS,TS,ROFF,SFCALB,VIS,THETA,ICLOUD,IWP,ISUB,TAUW,TOP,NCL,INTVAL,WH2O,CLOLWC,CDR
			! output - TRANS
			If(Theta .LE. 90) then
				CALL STRPSRB(LATMOS,TS,ROFF,albedo,VIS,THETA,0,IWP,ISUB,&	
				&     0.0D0,TOP,NCL,INTVAL,Temp,HR,Alt,CLOLWC,CDR,TCLEAR(i,j),TDIR(i,j))
				print *, TCLEAR(i,j),TDIR(i,j)
				
				CALL STRPSRB(LATMOS,TS,ROFF,albedo,VIS,THETA,1,IWP,ISUB,&	
				&    100.0D0,TOP,NCL,INTVAL,Temp,HR,Alt,CLOLWC,CDR,TCLOUD(i,j),TDIRn)
				print *, TCLOUD(i,j),TDIRn
			
			else
				TCLEAR(i,j) = -2.0
				TDIR(i,j) = -2.0
				TCLOUD(i,j) = -2.0
				
			end if
			
			!CALL D2STR(x(i),y(j),Temp,HR,Alt,albedo,vis,diaj,horad,0, 0.0D0,&
			!			&  TCLEAR(i,j),TDIR(i,j))

			
			
			!CALL D2STR(x(i),y(j),Temp,HR,Alt,albedo,vis,diaj,horad,1,100.0D0,&
			!			&  TCLOUD(i,j),TDIRn)
			
			Else  ! Si esta fuera de territorio Chileno.
				TCLEAR(i,j) = -1.0
				TDIR(i,j) = -1.0
				TCLOUD(i,j) = -1.0
			End if
		End do
	End do

	call check( nf90_put_var(ncid_rad, Global_varid, TCLEAR))
	call check( nf90_put_var(ncid_rad, Difusa_varid, TCLOUD))
	call check( nf90_put_var(ncid_rad, Directa_varid, TDIR))
	
 End do
 
 call check( nf90_close(ncid_rad) )
 
 
 999 Stop
 
 contains
 subroutine check(status)
  integer, intent ( in) :: status
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop 2
  end if
 end subroutine check
 
 
End Program trasmitancia



SUBROUTINE diajuliano (day, month, year, dayj)   
! Calculo de dia juliano

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

