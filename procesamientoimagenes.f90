! Copyright (C) 2010 - 2011, Juan Pablo Justiniano  <jpjustiniano@gmail.com>
! Programa para el filtraje y almacenamiento de imagenes satelitales en un archivo mensual.
! Funciona con horas en UTC.

! Variable coordinada de lat y long real.
! Coordinar las variables de entrada en las mismas coordenadas
! Procesamiento de solo Chile. Eliminar fuera de fronteras.

! Revisar:
! Problemas con maximos y minimos en meses que hay una sola imagen

!Canal 1: [1728, 9020]
!Canal 4: [432, 2255]

!Imagen Media Hora:
!Canal 1_SI: [1, 4380]
!Canal 1_ID: [1720, 5477]
 
!Canal 4_SI: [1, 1095]
!Canal 4_ID: [430, 1369]

 Program ProcesamientoImagenes_media_hora
 use netcdf
 !$ use OMP_LIb
 implicit none
 
 real :: latit = -33., longit = -70.
 ! Variables Programa
 Integer :: errorread
 character (4) :: ano
 character (2) :: mes, cdia, hora, minu
 character (10) :: foto
 character(8)  :: fecha
 real :: iano, imes, idia, ihora, ihorat, iminu, imesp=0., diaj, hora_utc
 Integer :: ifoto					! 1:media_hora; 2:south_full
 Real    :: az,el,ha,dec,soldst
 character (30) :: argument  ! Nombre de archivo de lista
 character (30) :: filename  ! Nombres dentro de la lista
 character (30) :: filenamepathin, pathin, pathout  ! Nombres dentro de la lista
 character (23) :: filename_out
 character (len = 27) :: string
 integer :: i=1,j=1, rec=1, noct= 0, k
 logical :: flag1= .false.
 integer,dimension(8) :: tiempo, tiempof
 integer, parameter :: r4i = 1095, r4f = 1368 , r1i = r4i*4, r1f = r4f*4 ! CH1: 4380 a 5472
 integer, parameter :: NLat= 432, Nlon= 2255
 integer, parameter :: NX=(r4f- r4i)+1, NY=430, NXf=(r4f-r4i)*4+1, NYf=(Ny)*4-3
 integer, dimension (Nlon,NLat) :: Latitud , Longitud
 Integer :: L1, L2, TID, TIDmax
 Real :: mm
 
 ! Variables NETCDF 
 integer :: ncid, ncid_in, dia =31, status
 integer, parameter :: NDIMS = 3, NDIMS_IN = 2	      ! We are writing 2D data.

 real :: x(NX), y(NY), xf(NXf), yf(NYf)
 integer, dimension (9020, 1728) :: CH1_out
 integer, dimension (Nlon,NLat) :: CH4_out   ! Matrices de archivo original
 integer, dimension(NX,NY)   :: CH4_in
 integer, dimension(NXf,NYf) :: CH1_in        !  Matrices de archivo recortado
 Integer, dimension(NXf,NYf) :: max1040,max1110,max1140,max1240,max1310,max1340,max1410,max1440,max1610
 Integer, dimension(NXf,NYf) :: min1040,min1110,min1140,min1240,min1310,min1340,min1410,min1440,min1610
 Integer, dimension(NXf,NYf) :: max1640, max1710,max1740,max1840,max1910,max1940,max2010,max2040,max2140
 Integer, dimension(NXf,NYf) :: min1640, min1710,min1740,min1840,min1910,min1940,min2010,min2040,min2140
 Integer, dimension(NXf,NYf) :: max2210, min2210, max2240, min2240
 integer, dimension((r4f- r4i)*4+1, (Ny-1)*4+1) :: Lat_CH1 , Lon_CH1
 integer, dimension (Nx,Ny)  :: Lat_CH4 , Lon_CH4
 integer :: x_dimid, y_dimid, xf_dimid, yf_dimid, dia_dimid, hora_dimid
 integer :: x_varid, y_varid, xf_varid, yf_varid, dia_varid, hora_varid
 integer :: CH1_in_varid, CH4_in_varid
 integer :: CH1_varid, Ch4_varid
 integer :: CH1_max_varid, CH1_min_varid
 integer :: Lat_CH4_varid, Lon_CH4_varid , Lat_CH1_varid, Lon_CH1_varid
 integer :: max1040_varid, max1110_varid, max1140_varid, max1240_varid, max1310_varid, max1340_varid
 integer :: max1410_varid, max1440_varid, max1610_varid, max1640_varid, max1710_varid, max1740_varid
 integer :: max1840_varid, max1910_varid, max1940_varid, max2010_varid, max2040_varid, max2140_varid
 integer :: max2210_varid, max2240_varid
 integer :: min1040_varid, min1110_varid, min1140_varid, min1240_varid, min1310_varid, min1340_varid
 integer :: min1410_varid, min1440_varid, min1610_varid, min1640_varid, min1710_varid, min1740_varid
 integer :: min1840_varid, min1910_varid, min1940_varid, min2010_varid, min2040_varid, min2140_varid
 integer :: min2210_varid, min2240_varid

 integer :: start(NDIMS), startf(NDIMS), count(NDIMS), countf(NDIMS), start_hora(1)
 integer :: dimids(NDIMS), dimids_fine(NDIMS), dimids2d(NDIMS_IN), dimids2df(NDIMS_IN)

 max1040= 0;max1110= 0;max1140= 0;max1240= 0;max1310= 0;max1340= 0;max1410= 0;max1440= 0;max1610= 0; max1640= 0
 max1710= 0;max1740= 0;max1840= 0;max1910= 0;max1940= 0;max2010= 0;max2040= 0;max2140= 0;max2210= 0; max2240= 0 

 min1040=10000;min1110=10000;min1140=10000;min1240=10000;min1310=10000;min1340=10000
 min1410=10000;min1440=10000;min1610=10000;min1710=10000;min1740=10000;min1840=10000;min1910=10000
 min1940=10000;min2010=10000;min2040=10000;min2140=10000;min2210=10000;min2240=10000
!********************************************************** Fin declaracion Variables
 
 ! Parametros de uso
 pathin = '/media/Elements/dm/'  ! Directorio de archivos de entrada.. NO incorporado aun.
 pathout = '/media/Elements/dm/' ! Directorio de archivos de salida.. NO incorporado aun
 TIDmax = 6    ! Numero de procesadores maximo a utilizar
 
 
 
!$    TID = omp_get_num_procs()
!$		If (TID>TIDmax) TID = TIDmax
!$    call OMP_SET_NUM_THREADS(TID)

 
 print *
 print *, '                          Recorte de imagenes NetCDF'
 print *, '                          **************************'
 print *
!$ print *, ' Numero de procesadores en uso: ',TID
 
 
 call date_and_time(DATE=fecha, VALUES=tiempo)
 open (unit=16, file='log.txt')

! Matrices de Geolocalizacion 
 open (unit=12, file='latitude.CH4.media_hora.txt', status= 'old', Action='read' )
 read (12,*) Latitud
 Close (12)
 open (unit=12, file='longitude.CH4.media_hora.txt', status= 'old', Action='read' )
 read (12,*) Longitud
 Close (12)
 
 Lat_CH4 = Latitud(r4i:r4f,:430)
 Lon_CH4 = Longitud(r4i:r4f,:430)
 
Do i = 1, Nx
	Do j =1, Ny
		Lat_CH1((i-1)*4+1, (j-1)*4+1) =  Lat_CH4(i,j)
		Lon_CH1((i-1)*4+1, (j-1)*4+1) =  Lon_CH4(i,j)
	End do
End do

Do j =1, Ny
	Do i = 1, Nx-1
		L1 =  Lat_CH4(i,j)
		L2 =  Lat_CH4(i+1,j)
		mm = (L2-L1)/4.	
		If ( abs(L1)<10 .or. abs(L2)<10 ) then 
			print *, i,j
		End if	
		Do k = 1, 3
			Lat_CH1((i-1)*4+1+k, (j-1)*4+1) = Nint(L1 + k*mm)
		End do	
			
		L1 =  Lon_CH4(i,j)
		L2 =  Lon_CH4(i+1,j)
		mm = (L2-L1)/4.	
		If ( abs(L1)<10 .or. abs(L2)<10 ) then
			print *, i,j
		End if
		Do k = 1, 3
			Lon_CH1((i-1)*4+1+k, (j-1)*4+1) = Nint(L1 + k*mm)
		End do	
	End do
End do

Do i =1, (r4f- r4i)*4+1
	Do j = 1, Ny-1
		L1 =  Lat_CH1(i,(j-1)*4+1)	
		L2 =  Lat_CH1(i,(j-1)*4+1+4)
		mm = (L2-L1)/4.	
		If ( abs(L1)<10 .or. abs(L2)<10 ) then 
			print *, i,j
		End if	
		Do k = 1, 3
			Lat_CH1(i, (j-1)*4+1+k) = Nint(L1 + k*mm)
		End do		
		
		L1 =  Lon_CH1(i,(j-1)*4+1)	
		L2 =  Lon_CH1(i,(j-1)*4+1+4)
		mm = (L2-L1)/4.	
		If ( abs(L1)<10 .or. abs(L2)<10 ) then 
			print *, i,j
		End if	
		Do k = 1, 3
			Lon_CH1(i, (j-1)*4+1+k) = Nint(L1 + k*mm)
		End do		
	End do
End do
 
! Fin de Geolocalizacion 

 
 call get_command_argument(1, argument) ! Nombre de archivo lista.txt
 open (unit=8, file=trim(argument), status='old', ACTION='READ', IOSTAT=errorread)  
 If(errorread/=0) then
    write (*,*) " Error en apertura de archivo con lista de entrada!" 
    write (16,*) " Error en apertura de archivo con lista de entrada!"
    go to 999 
 End If
 
100 read (8,*, IOSTAT=errorread) filename
filenamepathin = pathin//filename
 if(errorread == -1) then
    write (16,*)
    write (*,*)
    write (16,*) " Terminado exitoso del procesamiento de imagenes"
    write (*,*) " Terminado exitoso del procesamiento de imagenes"
    write (*,*)
    write (16,*)
    
    print *,'Maximas: '
    print *,maxval(max1040),maxval(max1110), maxval(max1140),maxval(max1240),maxval(max1310),maxval(max1410)
    print *,maxval(max1610),maxval(max1640), maxval(max1710),maxval(max1740),maxval(max1840),maxval(max1910)
    print *,maxval(max1940),maxval(max2010), maxval(max2040),maxval(max2140),maxval(max2210),maxval(max2240)
    
    call check( nf90_put_var(ncid, max1040_varid, max1040 ) )
    call check( nf90_put_var(ncid, max1110_varid, max1110 ) )
    call check( nf90_put_var(ncid, max1140_varid, max1140 ) )
    call check( nf90_put_var(ncid, max1240_varid, max1240 ) )
    call check( nf90_put_var(ncid, max1310_varid, max1310 ) )
    call check( nf90_put_var(ncid, max1340_varid, max1340 ) )
    call check( nf90_put_var(ncid, max1410_varid, max1410 ) )
    call check( nf90_put_var(ncid, max1440_varid, max1440 ) )
    call check( nf90_put_var(ncid, max1610_varid, max1610 ) )
    call check( nf90_put_var(ncid, max1640_varid, max1640 ) )
    call check( nf90_put_var(ncid, max1710_varid, max1710 ) )
    call check( nf90_put_var(ncid, max1740_varid, max1740 ) )
    call check( nf90_put_var(ncid, max1840_varid, max1840 ) )
    call check( nf90_put_var(ncid, max1910_varid, max1910 ) )
    call check( nf90_put_var(ncid, max1940_varid, max1940 ) )
    call check( nf90_put_var(ncid, max2010_varid, max2010 ) )
    call check( nf90_put_var(ncid, max2040_varid, max2040 ) )
    call check( nf90_put_var(ncid, max2140_varid, max2140 ) )
    call check( nf90_put_var(ncid, max2210_varid, max2210 ) )
    call check( nf90_put_var(ncid, max2240_varid, max2240 ) )
	print *,'Minimas: '
    print *,maxval(min1040),maxval(min1110), maxval(min1140),maxval(min1240),maxval(min1310),maxval(min1410)
    print *,maxval(min1610),maxval(min1640), maxval(min1710),maxval(min1740),maxval(min1840),maxval(min1910)
    print *,maxval(min1940),maxval(min2010), maxval(min2040),maxval(min2140),maxval(min2210),maxval(min2240)
    print *,
    call check( nf90_put_var(ncid, min1040_varid, min1040 ) )
    call check( nf90_put_var(ncid, min1110_varid, min1110 ) )
    call check( nf90_put_var(ncid, min1140_varid, min1140 ) )
    call check( nf90_put_var(ncid, min1240_varid, min1240 ) )
    call check( nf90_put_var(ncid, min1310_varid, min1310 ) )
    call check( nf90_put_var(ncid, min1340_varid, min1340 ) )
    call check( nf90_put_var(ncid, min1410_varid, min1410 ) )
    call check( nf90_put_var(ncid, min1440_varid, min1440 ) )
    call check( nf90_put_var(ncid, min1610_varid, min1610 ) )
    call check( nf90_put_var(ncid, min1640_varid, min1640 ) )
    call check( nf90_put_var(ncid, min1710_varid, min1710 ) )
    call check( nf90_put_var(ncid, min1740_varid, min1740 ) )
    call check( nf90_put_var(ncid, min1840_varid, min1840 ) )
    call check( nf90_put_var(ncid, min1910_varid, min1910 ) )
    call check( nf90_put_var(ncid, min1940_varid, min1940 ) )
    call check( nf90_put_var(ncid, min2010_varid, min2010 ) )
    call check( nf90_put_var(ncid, min2040_varid, min2040 ) )
    call check( nf90_put_var(ncid, min2140_varid, min2140 ) )
    call check( nf90_put_var(ncid, min2210_varid, min2210 ) )
    call check( nf90_put_var(ncid, min2240_varid, min2240 ) )
        
	call check( nf90_close(ncid) )
	
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

 If (foto=='media_hora') then
    ifoto = 1
 ElseIf (foto=='south_full') then
    ifoto = 2
 Else 
    ifoto = 0
 End If
 
 call diajuliano (idia, imes, iano, diaj)   ! entrada de reales en ves de enteros.
 call sunae(iano,diaj,ihora, latit, longit,az,el,ha,dec,soldst)  

 if (el < 7.0) then
    write (*,*) '   ',cdia,'  ', hora,':', minu, " eliminada, nocturna. ", el
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
    filename_out = ano//mes//'.media_hora.nc'
!    If (imes < 10 ) then
!		filename_out = ano//'0'//mes//'.nc'
!	Else
!		filename_out = ano//mes//'.nc'
!	End if
    write (*,*) filename_out
    write (*,*)
    write (16,*) filename_out
    write (16,*)
    call check( nf90_create(filename_out,NF90_64BIT_OFFSET, ncid) ) ! ncid, archivo de salida recortado
    
    ! Dimensiones
    call check( nf90_def_dim(ncid, "x", NX, x_dimid) )  ! Define the dimensions. 
    call check( nf90_def_dim(ncid, "y", NY, y_dimid) )
    call check( nf90_def_dim(ncid, "xf", NXf, xf_dimid) )  ! Define the dimensions. 
    call check( nf90_def_dim(ncid, "yf", NYf, yf_dimid) )
    call check( nf90_def_dim(ncid, "hora", NF90_UNLIMITED, hora_dimid) )
    
    ! Variables Coordinadas
    call check( nf90_def_var(ncid,"hora", NF90_FLOAT, hora_dimid, hora_varid) ) ! 32 bit   

    dimids =  (/ x_dimid, y_dimid, hora_dimid  /) ! The dimids array is used to pass the IDs of the dimensions
    dimids_fine =  (/ xf_dimid, yf_dimid, hora_dimid /)
	dimids2d =  (/ x_dimid, y_dimid/)
	dimids2df =  (/ xf_dimid, yf_dimid/)
	
    call check( nf90_def_var(ncid, "CH1", NF90_SHORT, dimids_fine, CH1_varid) )  
    call check( nf90_def_var(ncid, "CH4", NF90_SHORT, dimids, CH4_varid) )		!  NF90_Short (2-byte integer)
    call check( nf90_def_var(ncid, "Lat_CH4", NF90_SHORT, dimids2d, Lat_CH4_varid) )
    call check( nf90_def_var(ncid, "Lon_CH4", NF90_SHORT, dimids2d, Lon_CH4_varid) )
    call check( nf90_def_var(ncid, "Lat_CH1", NF90_SHORT, dimids2df, Lat_CH1_varid) )
    call check( nf90_def_var(ncid, "Lon_CH1", NF90_SHORT, dimids2df, Lon_CH1_varid) )
    
    call check( nf90_def_var(ncid, "max1040", NF90_SHORT, dimids2df, max1040_varid) )
    call check( nf90_def_var(ncid, "max1110", NF90_SHORT, dimids2df, max1110_varid) )
    call check( nf90_def_var(ncid, "max1140", NF90_SHORT, dimids2df, max1140_varid) )
    call check( nf90_def_var(ncid, "max1240", NF90_SHORT, dimids2df, max1240_varid) )
    call check( nf90_def_var(ncid, "max1310", NF90_SHORT, dimids2df, max1310_varid) )
    call check( nf90_def_var(ncid, "max1340", NF90_SHORT, dimids2df, max1340_varid) )
    call check( nf90_def_var(ncid, "max1410", NF90_SHORT, dimids2df, max1410_varid) )
    call check( nf90_def_var(ncid, "max1440", NF90_SHORT, dimids2df, max1440_varid) )
    call check( nf90_def_var(ncid, "max1610", NF90_SHORT, dimids2df, max1610_varid) )
    call check( nf90_def_var(ncid, "max1640", NF90_SHORT, dimids2df, max1640_varid) )
    call check( nf90_def_var(ncid, "max1710", NF90_SHORT, dimids2df, max1710_varid) )
    call check( nf90_def_var(ncid, "max1740", NF90_SHORT, dimids2df, max1740_varid) )
    call check( nf90_def_var(ncid, "max1840", NF90_SHORT, dimids2df, max1840_varid) )
    call check( nf90_def_var(ncid, "max1910", NF90_SHORT, dimids2df, max1910_varid) )
    call check( nf90_def_var(ncid, "max1940", NF90_SHORT, dimids2df, max1940_varid) )
    call check( nf90_def_var(ncid, "max2010", NF90_SHORT, dimids2df, max2010_varid) )
    call check( nf90_def_var(ncid, "max2040", NF90_SHORT, dimids2df, max2040_varid) )
    call check( nf90_def_var(ncid, "max2140", NF90_SHORT, dimids2df, max2140_varid) )
    call check( nf90_def_var(ncid, "max2210", NF90_SHORT, dimids2df, max2210_varid) )
    call check( nf90_def_var(ncid, "max2240", NF90_SHORT, dimids2df, max2240_varid) )
    
    call check( nf90_def_var(ncid, "min1040", NF90_SHORT, dimids2df, min1040_varid) )
    call check( nf90_def_var(ncid, "min1110", NF90_SHORT, dimids2df, min1110_varid) )
    call check( nf90_def_var(ncid, "min1140", NF90_SHORT, dimids2df, min1140_varid) )
    call check( nf90_def_var(ncid, "min1240", NF90_SHORT, dimids2df, min1240_varid) )
    call check( nf90_def_var(ncid, "min1310", NF90_SHORT, dimids2df, min1310_varid) )
    call check( nf90_def_var(ncid, "min1340", NF90_SHORT, dimids2df, min1340_varid) )
    call check( nf90_def_var(ncid, "min1410", NF90_SHORT, dimids2df, min1410_varid) )
    call check( nf90_def_var(ncid, "min1440", NF90_SHORT, dimids2df, min1440_varid) )
    call check( nf90_def_var(ncid, "min1610", NF90_SHORT, dimids2df, min1610_varid) )
    call check( nf90_def_var(ncid, "min1640", NF90_SHORT, dimids2df, min1640_varid) )
    call check( nf90_def_var(ncid, "min1710", NF90_SHORT, dimids2df, min1710_varid) )
    call check( nf90_def_var(ncid, "min1740", NF90_SHORT, dimids2df, min1740_varid) )
    call check( nf90_def_var(ncid, "min1840", NF90_SHORT, dimids2df, min1840_varid) )
    call check( nf90_def_var(ncid, "min1910", NF90_SHORT, dimids2df, min1910_varid) )
    call check( nf90_def_var(ncid, "min1940", NF90_SHORT, dimids2df, min1940_varid) )
    call check( nf90_def_var(ncid, "min2010", NF90_SHORT, dimids2df, min2010_varid) )
    call check( nf90_def_var(ncid, "min2040", NF90_SHORT, dimids2df, min2040_varid) )
    call check( nf90_def_var(ncid, "min2140", NF90_SHORT, dimids2df, min2140_varid) )
    call check( nf90_def_var(ncid, "min2210", NF90_SHORT, dimids2df, min2210_varid) )
    call check( nf90_def_var(ncid, "min2240", NF90_SHORT, dimids2df, min2240_varid) )

    
    ! Atributos de Geolocalizacion
    call check( nf90_put_att(ncid, Lat_CH4_varid, "units", "degrees_north"))
    call check( nf90_put_att(ncid, Lat_CH4_varid, "missing_value", -300) )
    call check( nf90_put_att(ncid, Lat_CH4_varid, "valid_min", -32768) )
    call check( nf90_put_att(ncid, Lat_CH4_varid, "valid_max", 32768) )
    call check( nf90_put_att(ncid, Lat_CH4_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, Lat_CH4_varid, "_CoordinateAxisType", "Lat_CH4") )
    call check( nf90_put_att(ncid, Lat_CH4_varid, "standard_name", "Latitud_CH4") )
    
    call check( nf90_put_att(ncid, Lon_CH4_varid, "units", "degrees_east") )
    call check( nf90_put_att(ncid, Lon_CH4_varid, "missing_value", -300) )
    call check( nf90_put_att(ncid, Lon_CH4_varid, "valid_min", -32768) )
    call check( nf90_put_att(ncid, Lon_CH4_varid, "valid_max", 32768) )
    call check( nf90_put_att(ncid, Lon_CH4_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, Lon_CH4_varid, "_CoordinateAxisType", "Lon_CH4") )
    call check( nf90_put_att(ncid, Lon_CH4_varid, "standard_name", "Longitud_CH4") )
        
    call check( nf90_put_att(ncid, Lat_CH1_varid, "units", "degrees_north"))
    call check( nf90_put_att(ncid, Lat_CH1_varid, "missing_value", -300) )
    call check( nf90_put_att(ncid, Lat_CH1_varid, "valid_min", -32768) )
    call check( nf90_put_att(ncid, Lat_CH1_varid, "valid_max", 32768) )
    call check( nf90_put_att(ncid, Lat_CH1_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, Lat_CH1_varid, "_CoordinateAxisType", "Lat_CH1") )
    
    call check( nf90_put_att(ncid, Lon_CH1_varid, "units", "degrees_east") )
    call check( nf90_put_att(ncid, Lon_CH1_varid, "missing_value", -300) )
    call check( nf90_put_att(ncid, Lon_CH1_varid, "valid_min", -32768) )
    call check( nf90_put_att(ncid, Lon_CH1_varid, "valid_max", 32768) )
    call check( nf90_put_att(ncid, Lon_CH1_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, Lon_CH1_varid, "_CoordinateAxisType", "Lon_CH1") )
     
    ! Atributos variables
    call check( nf90_put_att(ncid, CH1_varid, "units", "Albedo*100%") )
    call check( nf90_put_att(ncid, CH1_varid, "missing_value", -300) )
    call check( nf90_put_att(ncid, CH1_varid, "valid_min", -32768) )
    call check( nf90_put_att(ncid, CH1_varid, "valid_max", 32768) )
    call check( nf90_put_att(ncid, CH1_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, CH1_varid, "_CoordinateAxes", "time Lat_CH1 Lon_CH1") )
    call check( nf90_put_att(ncid, CH1_varid, "standard_name", "Canal Visible") )

    call check( nf90_put_att(ncid, CH4_varid, "units", "temp_deg_C") )
    call check( nf90_put_att(ncid, CH4_varid, "missing_value", -300) )
    call check( nf90_put_att(ncid, CH4_varid, "valid_min", -32768) )
    call check( nf90_put_att(ncid, CH4_varid, "valid_max", 32768) )
    call check( nf90_put_att(ncid, CH4_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, CH4_varid, "_CoordinateAxes", "time Lat_CH4 Lon_CH4") )
    call check( nf90_put_att(ncid, CH4_varid, "standard_name", "Canal Infrarojo") )
    
    call check( nf90_put_att(ncid, hora_varid, "units", "UTC_hours_from_day1"))
    call check( nf90_put_att(ncid, hora_varid, "_CoordinateAxisType", "time"))
      
    call check( nf90_put_att(ncid, max1040_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max1110_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max1140_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max1240_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max1310_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max1340_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max1410_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max1440_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max1610_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max1640_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max1710_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max1740_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max1840_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max1910_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max1940_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max2010_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max2040_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max2140_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max2210_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, max2240_varid, "scale_factor", 0.01) )
    
    call check( nf90_put_att(ncid, min1040_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min1110_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min1140_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min1240_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min1310_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min1340_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min1410_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min1440_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min1610_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min1640_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min1710_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min1740_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min1840_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min1910_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min1940_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min2010_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min2040_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min2140_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min2210_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, min2240_varid, "scale_factor", 0.01) )   

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

    ! Write the coordinate variable data. This will put the latitudes
    ! and longitudes of our data grid into the netCDF file.
    call check( nf90_put_var(ncid, Lat_CH4_varid, Lat_CH4) )
    call check( nf90_put_var(ncid, Lon_CH4_varid, Lon_CH4) )
    call check( nf90_put_var(ncid, Lat_CH1_varid, Lat_CH1) )
    call check( nf90_put_var(ncid, Lon_CH1_varid, Lon_CH1) )

    rec = 0
    noct= 0
    
 end if
 
 rec = rec + 1 
 start_hora =(/1/)
 start_hora(1) = rec
 count = (/ NX, NY, 1 /)
 countf = (/ NXf, NYf, 1 /)
 start = (/ 1, 1, 1 /)
 start(3) = rec
 
 ihorat = (ihora + (idia-1)*24)  ! Hora (hr_mes*100)
 call check( nf90_put_var(ncid, hora_varid, ihorat, start_hora)  )  ! Graba hora
 call check( nf90_open(trim(filename), nf90_nowrite, ncid_in) )  ! Abre archivo de lectura, ncid_in
 !call check( nf90_open(trim(filenamepathin), nf90_nowrite, ncid_in) )  ! Abre archivo de lectura en Disco externo! ncid_in
 
 call check( nf90_inq_varid(ncid_in, "gvar_ch1_fine", CH1_in_varid) )
 call check( nf90_inq_varid(ncid_in, "gvar_ch4", CH4_in_varid) )

		
 call check( nf90_get_var(ncid_in, CH4_in_varid, CH4_out)) 
 call check( nf90_get_var(ncid_in, CH1_in_varid, CH1_out))  
 
 ! Recorte de imagenes
 CH4_in = CH4_out(r4i:r4f,:NY)!925:1450
 CH1_in = CH1_out(r1i:r1f,:NYf)!3700:5800
 
 ! Revision de matriz
!$omp parallel 
!$omp do
 do i = 1, NXf
    do j = 1, NYf
    
		! Calculo de maximo y minimo mensual
	Select Case (NInt(ihora*10))
	Case (:107)
		if (CH1_in(i,j) > max1040(i,j)) then
			max1040(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min1040(i,j)) then
			min1040(i,j) = CH1_in(i,j)
		end if
	Case (111:113)
		if (CH1_in(i,j) > max1110(i,j)) then
			max1110(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min1110(i,j)) then
			min1110(i,j) = CH1_in(i,j)
		end if	
	Case (116:118)
		if (CH1_in(i,j) > max1140(i,j)) then
			max1140(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min1140(i,j)) then
			min1140(i,j) = CH1_in(i,j)
		end if	
	Case (121:128)
		if (CH1_in(i,j) > max1240(i,j)) then
			max1240(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min1240(i,j)) then
			min1240(i,j) = CH1_in(i,j)
		end if	
	Case (131:133)
		if (CH1_in(i,j) > max1310(i,j)) then
			max1310(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min1310(i,j)) then
			min1310(i,j) = CH1_in(i,j)
		end if		
	Case (136:138)
		if (CH1_in(i,j) > max1340(i,j)) then
			max1340(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min1340(i,j)) then
			min1340(i,j) = CH1_in(i,j)
		end if	
	Case (141:143)
		if (CH1_in(i,j) > max1410(i,j)) then
			max1410(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min1410(i,j)) then
			min1410(i,j) = CH1_in(i,j)
		end if		
	Case (146:148)
		if (CH1_in(i,j) > max1440(i,j)) then
			max1440(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min1440(i,j)) then
			min1440(i,j) = CH1_in(i,j)
		end if		
	Case (161:163)
		if (CH1_in(i,j) > max1610(i,j)) then
			max1610(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min1610(i,j)) then
			min1610(i,j) = CH1_in(i,j)
		end if		
	Case (166:168)
		if (CH1_in(i,j) > max1640(i,j)) then
			max1640(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min1640(i,j)) then
			min1640(i,j) = CH1_in(i,j)
		end if			
	Case (171:173)
		if (CH1_in(i,j) > max1710(i,j)) then
			max1710(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min1710(i,j)) then
			min1710(i,j) = CH1_in(i,j)
		end if		
	Case (176:178)
		if (CH1_in(i,j) > max1740(i,j)) then
			max1740(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min1740(i,j)) then
			min1740(i,j) = CH1_in(i,j)
		end if		
	Case (186:188)
		if (CH1_in(i,j) > max1840(i,j)) then
			max1840(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min1840(i,j)) then
			min1840(i,j) = CH1_in(i,j)
		end if		
	Case (191:193)
		if (CH1_in(i,j) > max1910(i,j)) then
			max1910(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min1910(i,j)) then
			min1910(i,j) = CH1_in(i,j)
		end if	
	Case (196:198)
		if (CH1_in(i,j) > max1940(i,j)) then
			max1940(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min1940(i,j)) then
			min1940(i,j) = CH1_in(i,j)
		end if			
	Case (201:203)
		if (CH1_in(i,j) > max2010(i,j)) then
			max2010(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min2010(i,j)) then
			min2010(i,j) = CH1_in(i,j)
		end if			
	Case (206:208)
		if (CH1_in(i,j) > max2040(i,j)) then
			max2040(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min2040(i,j)) then
			min2040(i,j) = CH1_in(i,j)
		end if			
	Case (216:218)
		if (CH1_in(i,j) > max2140(i,j)) then
			max2140(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min2140(i,j)) then
			min2140(i,j) = CH1_in(i,j)
		end if			
	Case (221:223)
		if (CH1_in(i,j) > max2210(i,j)) then
			max2210(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min2210(i,j)) then
			min2210(i,j) = CH1_in(i,j)
		end if			
	Case (226:)
		if (CH1_in(i,j) > max2240(i,j)) then
			max2240(i,j) = CH1_in(i,j)
		end if 								
		if (CH1_in(i,j) < min2240(i,j)) then
			min2240(i,j) = CH1_in(i,j)
		end if			
	Case Default
		print *, 'Imagen:',filename,' fuera de rango de max y min.'
		write (16,*) 'Imagen:',filename,' fuera de rango de max y min.'
	End Select	
		
	
	!end if
		! Correccion de datos 
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
!$omp end do 

!$omp do
 do i = 1, NX
    do j = 1, NY 
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
!$omp end do 
!$omp end parallel


 ! Guardado de matriz recortada y revisada
 call check( nf90_put_var(ncid,CH4_varid, CH4_in, start = start, &
									count = count) )
 call check( nf90_put_var(ncid,CH1_varid, CH1_in, start, countf))
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
 
 
 End Program ProcesamientoImagenes_media_hora

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


 SUBROUTINE sunae(year,day,hour, lat, long,az,el,ha,dec,soldst)  ! Sun's position, Michalsky
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
      latrad=lat*rad

!   calculate azimuth and elevation
      el=asin(sin(dec)*sin(latrad)+cos(dec)*cos(latrad)*cos(ha))
!      az=asin(-cos(dec)*sin(ha)/cos(el))

!!   this puts azimuth between 0 and 2*pi radians
!      if(sin(dec)-sin(el)*sin(latrad).ge.0.) then
!		if(sin(az).lt.0.) az=az+twopi
!      else
!      az=pi-az
!      endif
!   if az=90 degs, elcritical=asin(sin(dec)/sin(latrad))
!    elc=asin(sin(dec)/sin(latrad))
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
