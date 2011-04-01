! Copyright (C) 2010 - 2011, Juan Pablo Justiniano  <jpjustiniano@gmail.com>
! Programa para el calculo de la radiacion Global, Difusa y Directa.


! Revisar:
!	Revisar que pasa con pixeles en borde donde no es procesado el pixel lateral
!	Dejar la matriz de visibilidad procesado segun altura
!	No esta calculando los atos de la primera linea de NYf

!Canal 1: [1728, 9020]
!Canal 4: [432, 2255]

!Imagen Media Hora:
!Canal 1_SI: [1, 4380]
!Canal 1_ID: [1720, 5472]
 
!Canal 4_SI: [1, 1095]
!Canal 4_ID: [430, 1368]

 Program trasmitancia
 use netcdf
 !$ use OMP_LIb
 implicit none
 
 ! Variables Programa
 integer :: i, j, rec, ii,jj, errorread
 real ano, mes, diaj, dia 
 Real horad
 character (30) :: argument  ! Nombre de archivo 
 character(8)  :: fecha
 integer,dimension(8) :: tiempo, tiempof, tiempoa
 real(8), parameter :: pi= acos(-1.0), cdr= pi/180.0
 Logical :: Flag1 = .true.
 Real(8) :: Tempt,Altt,HRt,Albedot,vist
 
 ! Variables NETCDF archivo entrada
 integer :: ncid, ncid_in, ncid_var, status 
 integer, parameter :: NDIMS = 3, NDIMS_IN	= 2     
 integer :: NX, NY, NXf, NYf, Nhora, nboe
 real, dimension(:), allocatable :: x, y, xf, yf, hora
 integer, dimension(:,:), allocatable :: CH1_out, CH4_out   	! Matrices de archivo original
 integer, dimension(:,:), allocatable :: CH4_in , CH1_in        !  Matrices de archivo recortado
 Integer, dimension(:,:), allocatable :: CH1_max, CH1_min
 integer, dimension(:,:), allocatable :: Lat_CH1 , Lon_CH1, Lat_CH4 , Lon_CH4
 integer, dimension(:,:), allocatable :: Alt, Temp, HR, Vis, Albedo
 integer, dimension(:,:), allocatable :: fboe
 integer, dimension(:), allocatable :: dboe
 integer :: x_dimid, y_dimid, xf_dimid, yf_dimid, dia_dimid, hora_dimid
 integer :: x_varid, y_varid, xf_varid, yf_varid, dia_varid, hora_varid
 integer :: CH1_in_varid, CH4_in_varid
 integer :: CH1_varid, Ch4_varid, CH1_max_varid, CH1_min_varid
 integer :: Lat_CH4_varid, Lon_CH4_varid , Lat_CH1_varid, Lon_CH1_varid
 integer :: Lat_CH4_rad_varid, Lon_CH4_rad_varid , Lat_CH1_rad_varid, Lon_CH1_rad_varid 
 integer :: Alt_varid, Temp_varid, Albedo_varid, Vis_varid, HR_varid
 integer :: max1040_varid, max1110_varid, max1140_varid, max1240_varid, max1310_varid, max1340_varid
 integer :: max1410_varid, max1440_varid, max1610_varid, max1640_varid, max1710_varid, max1740_varid
 integer :: max1840_varid, max1910_varid, max1940_varid, max2010_varid, max2040_varid, max2140_varid
 integer :: max2210_varid, max2240_varid
 integer :: min1040_varid, min1110_varid, min1140_varid, min1240_varid, min1310_varid, min1340_varid
 integer :: min1410_varid, min1440_varid, min1610_varid, min1640_varid, min1710_varid, min1740_varid
 integer :: min1840_varid, min1910_varid, min1940_varid, min2010_varid, min2040_varid, min2140_varid
 integer :: min2210_varid, min2240_varid
 
 integer :: start(NDIMS), count(NDIMS), countf(NDIMS), startm(NDIMS)
 integer :: dimids(NDIMS), dimids_fine(NDIMS), dimids2d(NDIMS_IN), dimids2df(NDIMS_IN)
 
 ! Variables NETCDF archivo salida
 character(4)  :: cano
 character(2)  :: chora,cmin, cmes, cdia
 character(25) :: filename_out, ccname
 character (30) :: filenamepathin, pathin, pathout
 integer :: ncid_rad
 integer :: x_dimid_rad, y_dimid_rad,hora_dimid_rad, xf_dimid_rad, yf_dimid_rad
 integer :: x_varid_rad, y_varid_rad, xf_varid_rad, yf_varid_rad, hora_varid_rad
 integer :: Global_varid, Difusa_varid, Directa_varid
 integer(2) :: valid_max = 15000, valid_min=-3
 
 ! Variables OpenMPI
 integer:: TID, TIDmax

 ! Variables Brasil
 Real(8) :: E0,DEC,ET
 Real(8) :: DECR, YLATR, CODEC, COLAT, SIDEC, SILAT
 Real(8) :: ZN, TIMCOR, ROFF, TOP, CLOLWC
 Real(8) :: TS, TSOLAR, WSOLAR, COWI, COSZEN, THETA
 Real(8) :: T1SFC = 300.0, T2SFC = 294.0, T4SFC = 287.0, T3SFC = 272.2, T5SFC = 257.1
 real(8) :: TDIRn, TCLOUDn, TCLEARn, TDIRn2
 integer :: LATMOS, IWP, ISUB, NCL, INTVAL
 real(8) :: COWSR, TSRA, WSR, TSSA, TSRS, TSSS
 real(8) :: XI0, XIM, RSFCN, G, GLINHA, BETA, TAUW, DTAUW, DELTAW
 Integer, dimension(:,:), allocatable :: Global, Directa, cc !cc de boetto
 real(8) :: Glob, Dir
 
!********************************************************** Fin declaracion Variables

 ! Parametros de uso
 pathin = '/media/Elements/dm/'  ! Directorio de archivos de entrada.. NO implementado aun.
 pathout = '/media/Elements/dm/' ! Directorio de archivos de salida.. NO implementado aun
 TIDmax = 6    ! Numero de procesadores maximo a utilizar
 
 !$    TID = omp_get_num_procs()
 !$		If (TID>TIDmax) TID = TIDmax
 !$    call OMP_SET_NUM_THREADS(TID)

 
 print *
 print *, '                  Calculo de la Irradiacion Global y Directa'
 print *, '                  ******************************************'
 print *
 !$ print *, ' Numero de procesadores en uso: ',TID


 call date_and_time(DATE=fecha, VALUES=tiempo)
 open (unit=16, file='log_trans.txt')
 
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
 
call check( nf90_inq_varid(ncid, "max1040", max1040_varid) )
call check( nf90_inq_varid(ncid, "max1110", max1110_varid) )
call check( nf90_inq_varid(ncid, "max1140", max1140_varid) )
call check( nf90_inq_varid(ncid, "max1240", max1240_varid) )
call check( nf90_inq_varid(ncid, "max1310", max1310_varid) )
call check( nf90_inq_varid(ncid, "max1340", max1340_varid) )
call check( nf90_inq_varid(ncid, "max1410", max1410_varid) )
call check( nf90_inq_varid(ncid, "max1440", max1440_varid) )
call check( nf90_inq_varid(ncid, "max1610", max1610_varid) )
call check( nf90_inq_varid(ncid, "max1640", max1640_varid) )
call check( nf90_inq_varid(ncid, "max1710", max1710_varid) )
call check( nf90_inq_varid(ncid, "max1740", max1740_varid) )
call check( nf90_inq_varid(ncid, "max1840", max1840_varid) )
call check( nf90_inq_varid(ncid, "max1910", max1910_varid) )
call check( nf90_inq_varid(ncid, "max1940", max1940_varid) )
call check( nf90_inq_varid(ncid, "max2010", max2010_varid) )
call check( nf90_inq_varid(ncid, "max2040", max2040_varid) )
call check( nf90_inq_varid(ncid, "max2140", max2140_varid) )
call check( nf90_inq_varid(ncid, "max2210", max2210_varid) )
call check( nf90_inq_varid(ncid, "max2240", max2240_varid) )

call check( nf90_inq_varid(ncid, "min1040", min1040_varid) )
call check( nf90_inq_varid(ncid, "min1110", min1110_varid) )
call check( nf90_inq_varid(ncid, "min1140", min1140_varid) )
call check( nf90_inq_varid(ncid, "min1240", min1240_varid) )
call check( nf90_inq_varid(ncid, "min1310", min1310_varid) )
call check( nf90_inq_varid(ncid, "min1340", min1340_varid) )
call check( nf90_inq_varid(ncid, "min1410", min1410_varid) )
call check( nf90_inq_varid(ncid, "min1440", min1440_varid) )
call check( nf90_inq_varid(ncid, "min1610", min1610_varid) )
call check( nf90_inq_varid(ncid, "min1640", min1640_varid) )
call check( nf90_inq_varid(ncid, "min1710", min1710_varid) )
call check( nf90_inq_varid(ncid, "min1740", min1740_varid) )
call check( nf90_inq_varid(ncid, "min1840", min1840_varid) )
call check( nf90_inq_varid(ncid, "min1910", min1910_varid) )
call check( nf90_inq_varid(ncid, "min1940", min1940_varid) )
call check( nf90_inq_varid(ncid, "min2010", min2010_varid) )
call check( nf90_inq_varid(ncid, "min2040", min2040_varid) )
call check( nf90_inq_varid(ncid, "min2140", min2140_varid) )
call check( nf90_inq_varid(ncid, "min2210", min2210_varid) )
call check( nf90_inq_varid(ncid, "min2240", min2240_varid) )
 
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
 allocate(hora(Nhora))
 
 allocate (CH1_in (NXf,NYf))       ! imagenes invertidas
 allocate (CH4_in (NX,NY))

 allocate(Global(NXf,NYf))
 allocate(Directa(NXf,NYf))
 allocate(cc(NXf,NYf))
 
 Global = -3
 Directa= -3
 
 allocate (Lat_CH1 (NXf, NYf) )
 allocate (Lon_CH1 (NXf, NYf) )
 allocate (Lat_CH4 (Nx, Ny) )
 allocate (Lon_CH4 (Nx, Ny) )
 allocate (CH1_max(NXf,NYf),CH1_min(NXf,NYf))
 allocate (Alt(NXf,NYf))
 allocate (Albedo(NXf,NYf))
 allocate (Temp(NXf,NYf))
 allocate (HR(NXf,NYf))

 ! Lectura de matrices 
 call check( nf90_get_var(ncid, hora_varid, hora) )
 call check( nf90_get_var(ncid, Lat_CH4_varid, Lat_CH4) )
 call check( nf90_get_var(ncid, Lon_CH4_varid, Lon_CH4) )
 call check( nf90_get_var(ncid, Lat_CH1_varid, Lat_CH1) )
 call check( nf90_get_var(ncid, Lon_CH1_varid, Lon_CH1) )

 
 
 ! Lectura de pixeles de CH1 y CH4
 count = (/ NX, NY, 1 /)
 countf = (/ NXf, NYf, 1 /)
 start = (/ 1, 1, 1 /)
 startm = (/ 1, 1, 1 /)
 call check( nf90_get_var(ncid, hora_varid, hora))
 
 !**************************************************** Creacion de archivo de salida
 
 write (cano,'(I4)') NInt(ano)
 If (mes .le. 9 ) then
	write(cmes,'(I1)') NInt(mes)
	cmes = '0'//trim(cmes)
 Else
	write(cmes,'(I2)') NInt(mes)
 End if
 filename_out = cano//cmes//'.media_hora.rad.nc'
 
 call check( nf90_create(filename_out,NF90_64BIT_OFFSET, ncid_rad) ) ! ncid, archivo de salida recortado
 ! Dimensiones
 call check( nf90_def_dim(ncid_rad, "x", NX, x_dimid_rad) )  ! Define the dimensions. 
 call check( nf90_def_dim(ncid_rad, "y", NY, y_dimid_rad) )
 call check( nf90_def_dim(ncid_rad, "xf", NXf, xf_dimid_rad) )  ! Define the dimensions. 
 call check( nf90_def_dim(ncid_rad, "yf", NYf, yf_dimid_rad) )
 call check( nf90_def_dim(ncid_rad, "hora", Nhora, hora_dimid_rad) )
 
 ! Variables
 call check( nf90_def_var(ncid_rad,"hora", NF90_FLOAT, hora_dimid_rad, hora_varid_rad) ) ! 32 bit
     

 dimids =  (/ xf_dimid_rad, yf_dimid_rad, hora_dimid_rad /) ! The dimids array is used to pass the IDs of the dimensions
 dimids2d =  (/ x_dimid, y_dimid/)
 dimids2df =  (/ xf_dimid, yf_dimid/)

 call check( nf90_def_var(ncid_rad, "Global", NF90_SHORT, dimids, Global_varid) )
 call check( nf90_def_var(ncid_rad, "Directa", NF90_SHORT, dimids, Directa_varid) )
 call check( nf90_def_var(ncid_rad, "Lat_CH4", NF90_SHORT, dimids2d, Lat_CH4_rad_varid) )
 call check( nf90_def_var(ncid_rad, "Lon_CH4", NF90_SHORT, dimids2d, Lon_CH4_rad_varid) )
 call check( nf90_def_var(ncid_rad, "Lat_CH1", NF90_SHORT, dimids2df, Lat_CH1_rad_varid) )
 call check( nf90_def_var(ncid_rad, "Lon_CH1", NF90_SHORT, dimids2df, Lon_CH1_rad_varid) )
 
 ! Atributos globales    
 call check( nf90_put_att(ncid_rad, NF90_GLOBAL, "imagenes", "media_hora") )
 call check( nf90_put_att(ncid_rad, NF90_GLOBAL, "Procesadox", "JPJ") )
 call check( nf90_put_att(ncid_rad, NF90_GLOBAL, "date", fecha) )
 call check( nf90_put_att(ncid_rad, NF90_GLOBAL, "ano", ano))
 call check( nf90_put_att(ncid_rad, NF90_GLOBAL, "mes", mes))
 call check( nf90_put_att(ncid_rad, NF90_GLOBAL, "NX", NX) )  
 call check( nf90_put_att(ncid_rad, NF90_GLOBAL, "NY", NY) )  
 call check( nf90_put_att(ncid_rad, NF90_GLOBAL, "NXf", NXf) )
 call check( nf90_put_att(ncid_rad, NF90_GLOBAL, "NYf", NYf) )
 
 ! Atributos variables
 call check( nf90_put_att(ncid_rad, Global_varid, "units", "W/m2") )
 call check( nf90_put_att(ncid_rad, Global_varid, "missing_value", valid_min) )
 call check( nf90_put_att(ncid_rad, Global_varid, "valid_min", valid_min) )
 call check( nf90_put_att(ncid_rad, Global_varid, "valid_max", valid_max) )
 call check( nf90_put_att(ncid_rad, Global_varid, "scale_factor", 0.1) )
 call check( nf90_put_att(ncid_rad, Global_varid, "_CoordinateAxes", "time Lat_CH1 Lon_CH1") )
 call check( nf90_put_att(ncid_rad, Global_varid, "standard_name", "Radiación Global") )
 
 call check( nf90_put_att(ncid_rad, Directa_varid, "units", "W/m2") )
 call check( nf90_put_att(ncid_rad, Directa_varid, "missing_value", valid_min) )
 call check( nf90_put_att(ncid_rad, Directa_varid, "valid_min", valid_min) )
 call check( nf90_put_att(ncid_rad, Directa_varid, "valid_max", valid_max) )
 call check( nf90_put_att(ncid_rad, Directa_varid, "scale_factor", 0.1) )
 call check( nf90_put_att(ncid_rad, Directa_varid, "_CoordinateAxes", "time Lat_CH1 Lon_CH1") )
 call check( nf90_put_att(ncid_rad, Directa_varid, "standard_name", "Radiación Global") )
 
 ! Atributos de Geolocalizacion
call check( nf90_put_att(ncid_rad, Lat_CH4_rad_varid, "units", "degrees_north"))
call check( nf90_put_att(ncid_rad, Lat_CH4_rad_varid, "missing_value", -32767) )
call check( nf90_put_att(ncid_rad, Lat_CH4_rad_varid, "valid_min", -32767) )
call check( nf90_put_att(ncid_rad, Lat_CH4_rad_varid, "valid_max", 32767) )
call check( nf90_put_att(ncid_rad, Lat_CH4_rad_varid, "scale_factor", 0.01) )
call check( nf90_put_att(ncid_rad, Lat_CH4_rad_varid, "_CoordinateAxisType", "Lat_CH4") )

call check( nf90_put_att(ncid_rad, Lon_CH4_rad_varid, "units", "degrees_east") )
!call check( nf90_put_att(ncid_rad, Lon_CH4_rad_varid, "missing_value", -32768) )
call check( nf90_put_att(ncid_rad, Lon_CH4_rad_varid, "valid_min", -32767) )
call check( nf90_put_att(ncid_rad, Lon_CH4_rad_varid, "valid_max", 32767) )
call check( nf90_put_att(ncid_rad, Lon_CH4_rad_varid, "scale_factor", 0.01) )
call check( nf90_put_att(ncid_rad, Lon_CH4_rad_varid, "_CoordinateAxisType", "Lon_CH4") )

call check( nf90_put_att(ncid_rad, Lat_CH1_rad_varid, "units", "degrees_north"))
!call check( nf90_put_att(ncid_rad, Lat_CH1_rad_varid, "missing_value", -32768) )
call check( nf90_put_att(ncid_rad, Lat_CH1_rad_varid, "valid_min", -32767) )
call check( nf90_put_att(ncid_rad, Lat_CH1_rad_varid, "valid_max", 32767) )
call check( nf90_put_att(ncid_rad, Lat_CH1_rad_varid, "scale_factor", 0.01) )
call check( nf90_put_att(ncid_rad, Lat_CH1_rad_varid, "_CoordinateAxisType", "Lat_CH1") )

call check( nf90_put_att(ncid_rad, Lon_CH1_rad_varid, "units", "degrees_east") )
!call check( nf90_put_att(ncid_rad, Lon_CH1_rad_varid, "missing_value", -32768) )
call check( nf90_put_att(ncid_rad, Lon_CH1_rad_varid, "valid_min", -32767) )
call check( nf90_put_att(ncid_rad, Lon_CH1_rad_varid, "valid_max", 32767) )
call check( nf90_put_att(ncid_rad, Lon_CH1_rad_varid, "scale_factor", 0.01) )
call check( nf90_put_att(ncid_rad, Lon_CH1_rad_varid, "_CoordinateAxisType", "Lon_CH1") )


 
 ! End define mode.
 call check( nf90_enddef(ncid_rad) )
 !**************************************************** Procesamiento de las imagenes Sat.
 
 call check( nf90_put_var(ncid_rad, Lat_CH4_rad_varid, Lat_CH4) )
 call check( nf90_put_var(ncid_rad, Lon_CH4_rad_varid, Lon_CH4) )
 call check( nf90_put_var(ncid_rad, Lat_CH1_rad_varid, Lat_CH1) )
 call check( nf90_put_var(ncid_rad, Lon_CH1_rad_varid, Lon_CH1) )
 call check( nf90_put_var(ncid_rad, hora_varid_rad, hora) )

 ! Fixed input parameters for subroutine strpsrb
 ROFF   = 0.0		!     ROFF    : difference in water vapor
 IWP    = 3			!     IWP     : cloud droplet size distribution
 ISUB   = 2			!     ISUB    : use subroutine WOLKE1 (ISUB=1) or WOLKE2 (ISUB=2)
 TOP    = 500.0		!     TOP     : cloud top (mbar)
 NCL    = 2			!     NCL     : number of cloud layers
 INTVAL = 3			!     INTVAL  : spectral interval
 CLOLWC = 0.0		!     CLOLWC  : cloud liquid water content
 
 !calculate Diffuse radiation 
 G = 0.85 						! cloud asymmetry parameter
 GLINHA = (G - G*G)/(1 - G*G)	! corrected G
 BETA = 0.5 - 3.*GLINHA/8. - 7.*GLINHA**3/128. - 9.*GLINHA**5/128. ! cloud backscatter coefficient

 ! Variables meteorologicas
! call check( nf90_open('variables.nc', nf90_nowrite, ncid_var) )
! call check( nf90_inq_varid(ncid_var, "Temp", Temp_varid) )
! call check( nf90_inq_varid(ncid_var, "HR", HR_varid) )
! call check( nf90_inq_varid(ncid_var, "Vis", Vis_varid) )
! call check( nf90_inq_varid(ncid_var, "Albedo", Albedo_varid) )
! call check( nf90_inq_varid(ncid_var, "Alt", Alt_varid) )
 
 !call check( nf90_get_var (ncid_var, Alt_varid, Alt) )
 
! ! Boetto3
! open (15,file='trainningmatrix.txt', status='old', ACTION='READ', IOSTAT=errorread)
! if (errorread/=0) Then
!	print *, ' No se encuentra archivo con matriz de entrenamiento'
!	exit
! end if	
 
! nboe =0
! Do
!	read (15,*, iostat= errorread) 
!	if (errorread/=0) exit
!	nboe = nboe+1
!end do
 
! allocate (nboe, 4) fboe
! allocate (nboe) dboe
! do i=1, nboe
!  read (15,1500) fboe(i), dboe(i)
! end do
 
 
!             C1  C4  SP  mes clase
!1500 format ( I4, I4, I4, I2, I2) ! Lectura de p

 !/BOetto
 
 
 
 do rec = 1, Nhora
	Directa = -1
	Global = -1

	print *,'rec: ',rec
	
	dia = int(hora(rec)/24)
	horad = hora(rec)- dia*24

	call diajuliano (dia, mes, ano, diaj) 
	
	start(3) = rec !!!  mes !
	startm(3) = nint(mes)
	
	!Lectura de datos climatologicos
	!call check( nf90_get_var (ncid_var, Temp_varid, Temp, startm, countf) )
	!call check( nf90_get_var (ncid_var, HR_varid, HR, startm, countf) )
	!call check( nf90_get_var (ncid_var, Vis_varid, vis, startm, count = countf) )
	!call check( nf90_get_var (ncid_var, Albedo_varid, Albedo, startm, count = countf) )
	
!	  Temp = 293. ! (K)
!	  HR = .4 	! (0-1. ) 
!	  albedo = .20! (0-1.)
!	  vis = 10.	! (km)
		
	
	call check( nf90_get_var(ncid, CH1_varid, CH1_in, start, countf) ) 					   		   
	call check( nf90_get_var(ncid, CH4_varid, CH4_in, start, count) )
	
	Select Case (NInt(horad*10)) ! Selector de maximo y minimos
	Case (:107)
		CH1_min_varid = min1040_varid
		CH1_max_varid = max1040_varid
	Case (111:112)
		CH1_min_varid = min1110_varid
		CH1_max_varid = max1110_varid
	Case (116:118)
		CH1_min_varid = min1140_varid
		CH1_max_varid = max1140_varid
	Case (121:127)
		CH1_min_varid = min1240_varid
		CH1_max_varid = max1240_varid
	Case (131:132)
		CH1_min_varid = min1310_varid
		CH1_max_varid = max1310_varid		
	Case (136:137)
		CH1_min_varid = min1340_varid
		CH1_max_varid = max1340_varid
	Case (141:142)
		CH1_min_varid = min1410_varid
		CH1_max_varid = max1410_varid
	Case (146:148)
		CH1_min_varid = min1440_varid
		CH1_max_varid = max1440_varid	
	Case (161:162)
		CH1_min_varid = min1610_varid
		CH1_max_varid = max1610_varid		
	Case (166:168)
		CH1_min_varid = min1640_varid
		CH1_max_varid = max1640_varid	
	Case (171:172)
		CH1_min_varid = min1710_varid
		CH1_max_varid = max1710_varid
	Case (176:178)
		CH1_min_varid = min1740_varid
		CH1_max_varid = max1740_varid
	Case (186:188)
		CH1_min_varid = min1840_varid
		CH1_max_varid = max1840_varid
	Case (191:192)
		CH1_min_varid = min1910_varid
		CH1_max_varid = max1910_varid
	Case (196:198)
		CH1_min_varid = min1940_varid
		CH1_max_varid = max1940_varid	
	Case (201:202)
		CH1_min_varid = min2010_varid
		CH1_max_varid = max2010_varid
	Case (206:208)
		CH1_min_varid = min2040_varid
		CH1_max_varid = max2040_varid		
	Case (216:218)
		CH1_min_varid = min2140_varid
		CH1_max_varid = max2140_varid	
	Case (221:222)
		CH1_min_varid = min2210_varid
		CH1_max_varid = max2210_varid
	Case (226:)
		CH1_min_varid = min2240_varid
		CH1_max_varid = max2240_varid	
	Case Default
		print *, ' Fuera de rango de max y min: ', rec, horad
		write (16,*) ' Fuera de rango de max y min: ',rec, horad 
		stop 2
	End Select	
		
	call check( nf90_get_var(ncid, CH1_min_varid, CH1_min) )
	call check( nf90_get_var(ncid, CH1_max_varid, CH1_max) )
	
! Prueba datos Boetto	
		write(chora,'(I2)') Int(horad) 
		write(cmin,'(I2)') NInt((horad-Int(horad))*60)
		If (dia < 10 ) then
			write(cdia,'(I1)') NInt(dia)+1
			cdia='0'//trim(cdia)
		Else
			write(cdia,'(I2)') NInt(dia)+1
		End if
		ccname=cano//cmes//cdia//chora//cmin//'.txt'
		write (*,*) ccname, cdia,dia, chora,hora(rec), cmin
		open (10,file=ccname, status='old', ACTION='READ', IOSTAT=errorread)
		if (errorread/=0) cycle
		
		read (10,*) ((cc(ii,jj), jj = 1, NYf), ii=1, NXf)

!		call LDA(NXf, Nyf, Nx, Ny,nboe, CH1_in, CH4_in, SP, nint(mes), fboe, dboe, DI )
			
! /Prueba datos Boetto			
	
	

	Do j = 1, NYf
		
		! subroutine ASTRO - calculation of eccentricity correction, declination an equation of time
		! input: JDAY ; output: E0,DEC,ET
		CALL ASTRO(diaj,E0,DEC,ET) 
		
		DECR = DEC*cdr
		YLATR= Lat_CH1(NXf/2,j)/100.*cdr 		! Optimizacion, Uso de latitud media.
		CODEC = cos(DECR)
		COLAT = cos(YLATR)
		SIDEC = sin(DECR)
		SILAT = sin(YLATR)
       

		!calculate time correction
		ZN = 0.0
		TIMCOR = (4.0*(15.0*ZN+Lon_CH1(i,j)/100.)+ET)/60.0
		TSOLAR = horad+TIMCOR    !para entrada con horario en UTC
		
		! Characteristic hour angle WI corresponding to the W1-W2 interval
		WSOLAR    = (12.00 - TSOLAR)*15.
		COWI  = COS(WSOLAR*CDR)		
		! Characteristic solar zenith angle THETA (o)
		COSZEN = SIDEC*SILAT + CODEC*COLAT*COWI
		THETA  = ACOS(COSZEN)/CDR
		XI0   = 1367.00*E0*COSZEN

		
!$omp parallel private(XIM,Tempt,Altt,HRt,Albedot,Vist,Latmos,TS,TIMCOR,TSOLAR,WSOLAR,COWI,&
!$omp			& Dir,Glob,TCLEARn,TDIRn,TCLOUDn)
!$omp do
		Do i=1, NXf 	
		! En funcion del archivo de altura se procesa o no el pixel.
		!Altt  = Alt(i,j)
		Altt  = 400.
			
		  IF ( Altt > 0. ) THEN 
			
			!Tempt = Temp(i,j)/100.+ 273.15
			Tempt = 300.
			!HRt   =  HR(i,j)/10000.
			HRt   = 0.4
			!Albedot = albedo(i,j)/100.
			Albedot = 0.15
			!Vist = vis(i,j)
			Vist = 100.
			! calculation of the visibility at the station as function of visibility and altitude
			vist   = vist * EXP( (LOG(100.0/vist)/1000.0)*Altt )
			! test for visibility between 2 and 150 km
			IF (Vist > 150.)  Vist = 150.
			IF (Vist <   2.)  Vist =   2.

			If (Albedot <0. .or. Albedot>1.) then
				print *,'Albedo', i,j,rec,Albedot 
				Albedot = .25
			end if
			If (HRt < 0. .or. HRt >1.) then
				print *,'HR', i,j,rec,HRt 				! No corregido en archivos de kiko
				HRt  = .2
			end if
			If (Tempt<260. .or. Tempt>315.) then
				print *,'Temp', i,j,rec,Tempt
				Tempt=293.
			end if
			
!			! Calculo de Cobertura de nubes.
			if ( CH1_max(i,j) > 2500 .and. CH1_max(i,j)-CH1_min(i,j)>0.1 ) Then    		! Comprobar valor optimo de minimo de Maximos
					XIM =  (CH1_in(i,j)-CH1_min(i,j))/((CH1_max(i,j)-CH1_min(i,j) ) *1.)
			Else
				!XIM = (CH1_in(i,j)-CH1_min(i,j))/((7000. - CH1_min(i,j))*1.)
				 XIM = 0.000
			End if

! Prueba datos Boetto	
			!XIM =cc(i,j)/10000.
! /Prueba datos Boetto				
			
			! Choose atmosphere by surface temperature
			IF (Tempt < 297.0 .AND. Tempt > 290.5) THEN
				LATMOS = 2
				TS     = Tempt - T2SFC
			Else if (tempt < 290.5 .AND. tempt > 279.5) THEN
				LATMOS = 4
				TS     = tempt - T4SFC
			Else IF(tempt < 279.5 .AND. tempt > 264.5) THEN
				LATMOS = 3
				TS     = tempt - T3SFC
			Else IF(tempt < 264.5) THEN
				LATMOS = 5
				TS     = tempt - T5SFC
			Else
				LATMOS = 1
				TS     = Tempt - T1SFC
			end if
		

			! Subroutine STRPSRB - calculate transmittance for clear sky
			! input  - LATMOS,TS,ROFF,SFCALB,VIS,THETA,ICLOUD,IWP,ISUB,TAUW,TOP,NCL,INTVAL,WH2O,CLOLWC,CDR
			! output - TRANS			
			If((mod(i,5)==0) .or. flag1 .or. i==1 .or. j==1) then		 !Optimizar con un llamado para los tres parametros radiativos
		!write (*,*) LATMOS,TS,ROFF,Albedot,vist,THETA,0,IWP,ISUB,0.0D0,TOP,NCL,INTVAL,Tempt,HRt,Altt
				CALL STRPSRB(LATMOS,TS,ROFF,Albedot,vist,THETA,0,IWP,ISUB,&	
				&     0.0D0,TOP,NCL,INTVAL,Tempt,HRt,Altt,CLOLWC,CDR,TCLEARn,TDIRn)

				CALL STRPSRB(LATMOS,TS,ROFF,albedot,vist,THETA,1,IWP,ISUB,&	
				&    100.0D0,TOP,NCL,INTVAL,Tempt,HRt,Altt,CLOLWC,CDR,TCLOUDn,TDIRn2)
				
				flag1= .false.
			end if 
			
			! Calculo de irradiancia global, difusa y directa a partir de transmitancias calculadas.

			Glob  = XI0 *((1.0 - XIM) * (TCLEARn - TCLOUDn) + TCLOUDn) !Revisar XIM de difusa
			if (Glob < 0.0)  Glob = 0.
			
			XIM = 1.0 - XIM		!!
			
			If ((XIM >= 0.95).AND.(XIM < 1.01)) then
				TAUW = 1.0
				DTAUW = EXP(-(1 - TAUW)/(BETA * TAUW))
				Dir = DTAUW * TDIRn * XI0
			else if ((XIM < 0.95).AND.(XIM >= 0.)) then
				TAUW = XIM + 0.05
				DTAUW = EXP(-(1 - TAUW)/(BETA * TAUW))
				Dir = DTAUW * TDIRn * XI0
			else
				DTAUW = -2.0
				Dir=-2.0
				print *, 'Error en (i,j):',i,j, CH1_in(i,j), CH1_max(i,j), CH1_min(i,j), XIM
			end if
			
			
			IF(Dir < 0.0) Dir = 0.0
						
			Directa(i,j) = NInt(Dir*10)
			Global(i,j)  = NInt(Glob*10)			
			
			Else  ! Si esta fuera de territorio Chileno.
				Directa(i,j) = -1.0
				Global(i,j)  = -1.0
111			End if  ! Fin de procesamiento de pixel con altura > 0.
		
		End do
		!$omp end do 
!$omp end parallel
		If (mod(j,5)==0) Flag1 = .true.
		
	End do


	call date_and_time(DATE=fecha, VALUES=tiempof)
	print *, Global(NXf-100,j-500),  Directa(NXf-100,j-500)  
	write(*,300) tiempof (6) - tiempoa (6), tiempof(7) - tiempoa(7)
	tiempoa = tiempof
	

	call check( nf90_put_var(ncid_rad, Global_varid, Global, start=start,count = countf))
	call check( nf90_put_var(ncid_rad, Directa_varid,Directa, start=start,count = countf))
	
 End do
 
 call check( nf90_close(ncid_var) )
 call check( nf90_close(ncid_rad) )
 
 call date_and_time(DATE=fecha, VALUES=tiempof)
 write (*,200) tiempof(5)-tiempo(5),tiempof(6)-tiempo(6),tiempof(7)-tiempo(7)
 					
200 Format ('   Tiempo procesamiento: ',I3' hr., ',I3 'min., ',I3, 'sec.')
300 Format ('   Tiempo procesamiento: ',I3 'min., ',I3, 'sec.')
 
 999 Stop
 
 contains
 subroutine check(status)
  integer, intent ( in) :: status
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop 2
  end if
 end subroutine check
 
 
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

End Program trasmitancia
