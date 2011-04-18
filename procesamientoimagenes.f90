! Copyright (C) 2010 - 2011, Juan Pablo Justiniano  <jpjustiniano@gmail.com>
! Programa para el filtraje y almacenamiento de imagenes satelitales en un archivo mensual.
! Funciona con horas en UTC.

! Variable coordinada de lat y long real.
! Coordinar las variables de entrada en las mismas coordenadas
! Procesamiento de solo Chile. Eliminar fuera de fronteras.

! Revisar:
! Problemas con maximos y minimos en meses que hay una sola imagen
! Eliminar maximos
! Mejoras a Correxion de pixeles
! Coordenadas de primeras 4 lineas estan malas

!Canal 1: [1728, 9020]
!Canal 4: [432, 2255]

!Imagen Media Hora:
!Canal 1_SI: [1, 4380]
!Canal 1_ID: [1720, 5472]
 
!Canal 4_SI: [1, 1095]
!Canal 4_ID: [430, 1368]

 Program ProcesamientoImagenesMediaHora
 use netcdf
 !use Subs
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
 character (60) :: filenamepathin
 character (30) :: pathin, pathout  ! Nombres dentro de la lista
 character (23) :: filename_out
 character (len = 27) :: string
 integer :: i=1,j=1, rec=1, noct= 0, k
 logical :: flag1= .false.
 integer,dimension(8) :: tiempo, tiempof
 integer, parameter :: r4i = 1095, r4f = 1368 , r1i = (r4i-1)*4+1, r1f = (r4f-1)*4+1 ! CH1: 4377 a 5469
 integer, parameter :: NLat= 432, Nlon= 2255
 integer, parameter :: NLatF= 4*NLat, NlonF= 4*Nlon
 integer, parameter :: NX=(r4f- r4i)+1, NY=430, NXf=(r4f-r4i)*4+1, NYf=(Ny-1)*4+1
 integer, dimension (Nlon,NLat) :: Latitud , Longitud, scatter_phase
 Integer :: L1, L2, TID, TIDmax
 Real :: mm
 
 
 !Variables Boetto
 real, dimension(:,:), allocatable :: fboe
 integer, dimension(:), allocatable :: dboe
 integer :: nboe, DeltaSP, SPhase, clase
 integer, dimension (NXf,NYf) :: Alt, DI
 integer, dimension(NXf,NYf,800) :: DI3d
 real, dimension (100,4) :: CI_0, CI_0i
 real, dimension (100) :: CI_cs, CI_csi
 real, dimension (NXf,NYf,100) :: CI_cs_pixel
 real, dimension (100,4) :: Ponderador
 real :: IC_max, clearsky, clearsky_pixel, cloudclass, pond
 real, dimension (NXf,NYf) :: CCI_eff
 integer :: CCI_varid
 
 ! Variables NETCDF 
 integer :: ncid, ncid_in, dia =31, status
 integer, parameter :: NDIMS = 3, NDIMS_IN = 2	      ! We are writing 2D data.


 real :: x(NX), y(NY), xf(NXf), yf(NYf)
 integer, dimension (NlonF, NLatF) :: CH1_out
 integer, dimension (Nlon,NLat) :: CH4_out   ! Matrices de archivo original
 integer, dimension (NX,NY)   :: CH4_in
 integer, dimension(NXf,NYf) :: CH1_in        !  Matrices de archivo recortado
 Integer, dimension(NXf,NYf) :: max1040,max1110,max1140,max1240,max1310,max1340,max1410,max1440,max1610
 Integer, dimension(NXf,NYf) :: min1040,min1110,min1140,min1240,min1310,min1340,min1410,min1440,min1610
 Integer, dimension(NXf,NYf) :: max1640, max1710,max1740,max1840,max1910,max1940,max2010,max2040,max2140
 Integer, dimension(NXf,NYf) :: min1640, min1710,min1740,min1840,min1910,min1940,min2010,min2040,min2140
 Integer, dimension(NXf,NYf) :: max2210, min2210, max2240, min2240
 integer, dimension(NXf,NYf) :: Lat_CH1 , Lon_CH1
 integer, dimension (Nx,Ny)  :: Lat_CH4 , Lon_CH4, SP, SP_month
 integer :: x_dimid, y_dimid, xf_dimid, yf_dimid, dia_dimid, hora_dimid
 integer :: x_varid, y_varid, xf_varid, yf_varid, dia_varid, hora_varid
 integer :: CH1_in_varid, CH4_in_varid
 integer :: CH1_varid, Ch4_varid, Alt_varid, ncid_var
 integer :: CH1_max_varid, CH1_min_varid
 integer :: Lat_CH4_varid, Lon_CH4_varid , Lat_CH1_varid, Lon_CH1_varid, scatter_phase_varid,sp_varid
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
 pathin = '/media/Elements/dm/'  ! Directorio de archivos de entrada.. 
 pathout = '/media/Elements/dm/' ! Directorio de archivos de salida.. NO incorporado aun
 TIDmax = 7    ! Numero de procesadores maximo a utilizar
 
 
 
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

! Matrices de Geolocalizacion y scatter phase
 open (unit=12, file='CH4.latitude.media_hora.txt', status= 'old', Action='read', IOSTAT=errorread )
 If(errorread/=0) print *,' Error de lectura de archivo: ', ' CH4.latitude.media_hora.txt'
 read (12,*) Latitud
 Close (12)
 open (unit=12, file='CH4.longitude.media_hora.txt', status= 'old', Action='read', IOSTAT=errorread )
 If(errorread/=0) print *,' Error de lectura de archivo: ', ' CH4.longitude.media_hora.txt'
 read (12,*) Longitud
 Close (12)
 open (unit=12, file='scatter_phase.txt', status= 'old', Action='read', IOSTAT=errorread )
 If(errorread/=0) print *,' Error de lectura de archivo: ', ' scatter_phase.txt'
 read (12,*) scatter_phase
 Close (12)
 
 Lat_CH4 = Latitud(r4i:r4f,1:430)
 Lon_CH4 = Longitud(r4i:r4f,1:430)
 SP_month = scatter_phase(r4i:r4f,1:430)

 
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

Do i =1, NXf
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
 
!/Fin de Geolocalizacion 

! ! Boetto3
! Variables meteorologicas
 call check( nf90_open('variables.nc', nf90_nowrite, ncid_var) )

 call check( nf90_inq_varid(ncid_var, "Alt", Alt_varid) )
 
 call check( nf90_get_var (ncid_var, Alt_varid, Alt) )
 
 call check( nf90_close(ncid_var) )
 
 CI_0=-9999
 CI_cs = -9999
 CI_cs_pixel = -9999
 

 open (15,file='trainningmatrix.txt', status='old', ACTION='READ', IOSTAT=errorread)
 if (errorread/=0) Then
	print *, ' No se encuentra archivo con matriz de entrenamiento'
	go to 999 
 end if	
 
 nboe =0
 Do
	read (15,*, iostat= errorread) 
	if (errorread/=0) exit
	nboe = nboe+1
 end do
 
 allocate (fboe (nboe, 4)) 
 allocate (dboe (nboe)) 
 do i=1, nboe
  read (15,1500) fboe(i,:), dboe(i)
 end do
 
 
!             C1  C4  SP  mes clase
1500 format ( F6.4, F6.2, F5.2, F2.0, I2) ! Lectura de p

!/BOetto



 
 call get_command_argument(1, argument) ! Nombre de archivo lista.txt
 open (unit=8, file=trim(argument), status='old', ACTION='READ', IOSTAT=errorread)  
 If(errorread/=0) then
    write (*,*) " Error en apertura de archivo con lista de entrada!" 
    write (16,*) " Error en apertura de archivo con lista de entrada!"
    go to 999 
 End If
 
100 read (8,*, IOSTAT=errorread) filename
filenamepathin = trim(pathin)//trim(filename)
 if(errorread == -1) then			! Fin de lista de archivos.
    
    !Reinicio lectura para procedimiento Boetto
    rec =0
    Close(8)
    open (unit=8, file=trim(argument), status='old', ACTION='READ', IOSTAT=errorread)
110 read (8,*, IOSTAT=errorread) filename
	filenamepathin = trim(pathin)//trim(filename)
	if(errorread == -1) go to 120 ! salida en caso de terimnar archivo de lista

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
	 
	 call diajuliano (idia, imes, iano, diaj)   ! entrada de reales en ves de enteros.
	 call sunae(iano,diaj,ihora, latit, longit,az,el,ha,dec,soldst)  

	 if (el < 3.0) goto 110
	 rec = rec + 1 
	 start_hora =(/1/)
	 start_hora(1) = rec
	 count = (/ NX, NY, 1 /)
	 countf = (/ NXf, NYf, 1 /)
	 start = (/ 1, 1, 1 /)
	 start(3) = rec
 
	 ihorat = (ihora + (idia-1)*24)  ! Hora (hr_mes*100)
	 call check( nf90_open(trim(filenamepathin), nf90_nowrite, ncid_in) )  ! Abre archivo de lectura en Disco externo! ncid_in
	 
	 call check( nf90_inq_varid(ncid_in, "gvar_ch1_fine", CH1_in_varid) )
	 status = nf90_inq_varid (ncid_in, "scatter_phase", scatter_phase_varid)
	 if(status /= nf90_noerr) then
		SP = SP_month
	 else 
		call check( nf90_get_var(ncid_in, scatter_phase_varid, scatter_phase))
		SP = scatter_phase(r4i:r4f,:NY)
	 end if	

	 call check( nf90_get_var(ncid_in, CH1_in_varid, CH1_out))  
	 
	 ! Recorte de imagenes
	 CH1_in = CH1_out(r1i:r1f,1:NYf)
	 
	 CCI_eff=0
	 
	do i = 1, NXf
		do j = 1, NYf
			SPhase= nint(SP(Nint((i-1)/4.+1), Nint((j-1)/4.+1))/100.)

			If (SPhase >100) SPhase = 100
			if (SPhase <1) SPhase = 1
			clase=DI3d(i,j,rec)
			clearsky=CI_cs(SPhase)
			clearsky_pixel=CI_cs_pixel(i,j,SPhase)
			Cloudclass=CI_0(SPhase, clase)
			Pond=Ponderador(SPhase,clase)
			If (clase==1) then
				if (clearsky_pixel/=-9999) then
					if (cloudclass==0 .or. cloudclass==clearsky_pixel) then 
						CCI_eff(i,j)=0
					else
						CCI_eff(i,j)=abs((CH1_in(i,j) -clearsky_pixel)/(cloudclass-clearsky_pixel))*pond 
					end if
				else
					if (cloudclass==0 .or. cloudclass==clearsky) then
						CCI_eff(i,j)=0
					else
						CCI_eff(i,j)=abs((CH1_in(i,j) -clearsky)/(cloudclass-clearsky))*pond
					end if
				end if
			else
				if (cloudclass==0 .or. cloudclass==clearsky) then
					CCI_eff(i,j) =0
				else
					CCI_eff(i,j)=abs((CH1_in(i,j) -clearsky)/(cloudclass-clearsky))*pond
				end if
			end if
		end do
	end do

	 go to 110
	 
    
    
    
    
    
120    write (16,*)
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
    write (*,200)  tiempof(6) - tiempo (6), tiempof(7) - tiempo(7)
	write (*,*) '    Archivos procesados: ' , rec	
	write (*,*) ' Archivos no procesados: ' , noct, '  (nocturnos)'
	write (*,*)
	write (16,200) tiempof(6) - tiempo (6), tiempof(7) - tiempo(7)
	write (16,*) ' Archivos procesados: ' , rec		
    close (16)
    
    go to 999
 Else if(errorread > 0) then
    write (16,*) " Error en lectura de nombre de archivo en archivo lista ", filenamepathin
    write (*,*) " Error en lectura de nombre de archivo en archivo lista ", filenamepathin
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

 if (el < 3.0) then !! Cambiado a 5º
    write (*,*) '   ',cdia,'  ', hora,':', minu, '        Imagen eliminada, nocturna. ', el
    write (16,*) '                                   ', trim(filename),'   Eliminada, nocturna. '
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
		write (*,200) tiempof(6) - tiempo (6), tiempof(7) - tiempo(7)
	End if	
    flag1 = .true.
    filename_out = ano//mes//'.media_hora.nc'

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
    call check( nf90_def_var(ncid, "scatter_phase", NF90_SHORT, dimids, sp_varid) )
    call check( nf90_def_var(ncid, "CCI", NF90_SHORT, dimids, CCI_varid) )
    
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
    call check( nf90_put_att(ncid, Lat_CH4_varid, "valid_max", 32767) )
    call check( nf90_put_att(ncid, Lat_CH4_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, Lat_CH4_varid, "_CoordinateAxisType", "Lat_CH4") )
    call check( nf90_put_att(ncid, Lat_CH4_varid, "standard_name", "Latitud_CH4") )
    
    call check( nf90_put_att(ncid, Lon_CH4_varid, "units", "degrees_east") )
    call check( nf90_put_att(ncid, Lon_CH4_varid, "missing_value", -300) )
    call check( nf90_put_att(ncid, Lon_CH4_varid, "valid_min", -32768) )
    call check( nf90_put_att(ncid, Lon_CH4_varid, "valid_max", 32767) )
    call check( nf90_put_att(ncid, Lon_CH4_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, Lon_CH4_varid, "_CoordinateAxisType", "Lon_CH4") )
    call check( nf90_put_att(ncid, Lon_CH4_varid, "standard_name", "Longitud_CH4") )
        
    call check( nf90_put_att(ncid, Lat_CH1_varid, "units", "degrees_north"))
    call check( nf90_put_att(ncid, Lat_CH1_varid, "missing_value", -300) )
    call check( nf90_put_att(ncid, Lat_CH1_varid, "valid_min", -32768) )
    call check( nf90_put_att(ncid, Lat_CH1_varid, "valid_max", 32767) )
    call check( nf90_put_att(ncid, Lat_CH1_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, Lat_CH1_varid, "_CoordinateAxisType", "Lat_CH1") )
    
    call check( nf90_put_att(ncid, Lon_CH1_varid, "units", "degrees_east") )
    call check( nf90_put_att(ncid, Lon_CH1_varid, "missing_value", -300) )
    call check( nf90_put_att(ncid, Lon_CH1_varid, "valid_min", -32768) )
    call check( nf90_put_att(ncid, Lon_CH1_varid, "valid_max", 32767) )
    call check( nf90_put_att(ncid, Lon_CH1_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, Lon_CH1_varid, "_CoordinateAxisType", "Lon_CH1") )
     
    ! Atributos variables
    call check( nf90_put_att(ncid, CH1_varid, "units", "Albedo*100%") )
    call check( nf90_put_att(ncid, CH1_varid, "missing_value", -300) )
    call check( nf90_put_att(ncid, CH1_varid, "valid_min", -32768) )
    call check( nf90_put_att(ncid, CH1_varid, "valid_max", 32767) )
    call check( nf90_put_att(ncid, CH1_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, CH1_varid, "_CoordinateAxes", "time Lat_CH1 Lon_CH1") )
    call check( nf90_put_att(ncid, CH1_varid, "standard_name", "Canal Visible") )

    call check( nf90_put_att(ncid, CH4_varid, "units", "temp_deg_C") )
    call check( nf90_put_att(ncid, CH4_varid, "missing_value", -300) )
    call check( nf90_put_att(ncid, CH4_varid, "valid_min", -32768) )
    call check( nf90_put_att(ncid, CH4_varid, "valid_max", 32767) )
    call check( nf90_put_att(ncid, CH4_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, CH4_varid, "_CoordinateAxes", "time Lat_CH4 Lon_CH4") )
    call check( nf90_put_att(ncid, CH4_varid, "standard_name", "Canal Infrarojo") )
    
    
    call check( nf90_put_att(ncid, sp_varid, "units", "degrees") )
    call check( nf90_put_att(ncid, sp_varid, "valid_min", -32768) )
    call check( nf90_put_att(ncid, sp_varid, "valid_max", 32767) )
    call check( nf90_put_att(ncid, sp_varid, "scale_factor", 0.01) )
    call check( nf90_put_att(ncid, sp_varid, "_CoordinateAxes", "time Lat_CH4 Lon_CH4") )
    call check( nf90_put_att(ncid, sp_varid, "standard_name", "Scatter Phase") )
    
    call check( nf90_put_att(ncid, hora_varid, "units", "UTC_hours_from_day1"))
    call check( nf90_put_att(ncid, hora_varid, "_CoordinateAxisType", "time"))
      
    call check( nf90_put_att(ncid, CCI_varid, "units", "%") ) 
    call check( nf90_put_att(ncid, CCI_varid, "scale_factor", 0.01) ) 
    call check( nf90_put_att(ncid, CCI_varid, "_CoordinateAxes", "time Lat_CH1 Lon_CH1") ) 
    call check( nf90_put_att(ncid, CCI_varid, "standard_name", "Cloud Cover Index") ) 
      
      
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
 !call check( nf90_open(trim(filename), nf90_nowrite, ncid_in) )  ! Abre archivo de lectura, ncid_in
 call check( nf90_open(trim(filenamepathin), nf90_nowrite, ncid_in) )  ! Abre archivo de lectura en Disco externo! ncid_in
 
 call check( nf90_inq_varid(ncid_in, "gvar_ch1_fine", CH1_in_varid) )
 call check( nf90_inq_varid(ncid_in, "gvar_ch4", CH4_in_varid) )
 status = nf90_inq_varid (ncid_in, "scatter_phase", scatter_phase_varid)
 if(status /= nf90_noerr) then
	SP = SP_month
	print *, '                   Imagen sin variable de Scatter Phase.     ', filename
	write (16,*) 'Imagen sin variable de Scatter Phase.    ', filename
 else 
 	call check( nf90_get_var(ncid_in, scatter_phase_varid, scatter_phase))
 	SP = scatter_phase(r4i:r4f,:NY)
 end if	

		
 call check( nf90_get_var(ncid_in, CH4_in_varid, CH4_out)) 
 call check( nf90_get_var(ncid_in, CH1_in_varid, CH1_out))  
 
 ! Recorte de imagenes
 CH4_in = CH4_out(r4i:r4f,1:NY)
 CH1_in = CH1_out(r1i:r1f,1:NYf)
 
 
 ! Revision de matriz
 do i = 1, NXf
    do j = 1, NYf
    

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
        		! Calculo de maximo y minimo mensual !!!
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

    end do
end do


 do i = 1, NX
    do j = 1, NY 
		if (CH4_in(i,j) > 15000) then
            write (16,*) "CH4_in(", i, ",",j, ") =", CH4_in(i, j) , filename, iano,diaj,ihora,"Corregido, vecinos"
   			CH4_in(i,j) = (CH4_in(i+1,j)+CH4_in(i,j+1)+CH4_in(i-1,j)+CH4_in(i,j-1))/4.
   			if (CH4_in(i, j) > 7000) CH4_in(i,j) = 7000
        else if (CH4_in(i, j) > 7000) then
            write (16,*) "CH4_in(", i, ",", j, ") =", CH4_in(i, j), filename, iano,diaj,ihora,"Corregido"
   			CH4_in(i,j) = 7000
        end if
        if (CH4_in(i, j) < -7000 ) then
			write (16,*) "CH4_in(", i, ",", j, ") =", CH4_in(i, j), filename, iano,diaj,ihora,"Pixel malo"
			If (i==1 .or. i==NX .or. j==1 .or. j==Ny) then
				CH4_in(i,j) = 0. ! Poner tag de error
			else If (CH4_in(i+1,j) > -7000 .and. CH4_in(i,j+1) > -7000 .and. CH4_in(i-1,j) > -7000 .and. CH4_in(i,j-1)> -7000) Then
				CH4_in(i,j) = (CH4_in(i+1,j)+CH4_in(i,j+1)+CH4_in(i-1,j)+CH4_in(i,j-1))/4.   
				if (CH4_in(i, j) < -7000) CH4_in(i,j) = -7000     
			Else
				CH4_in(i,j) = (CH4_in(i+1,j)+CH4_in(i-1,j))/2.
				if (CH4_in(i, j) < -7000) CH4_in(i,j) = -7000
			End If
        end if
    end do
 end do

 
! Prueba datos Boetto	
 call LDA(NXf, Nyf, Nx, Ny,nboe, (CH1_in/10000.),CH4_in/(-100.), (SP/100.), nint(imes), fboe, dboe, DI )	

! Calculo de Max y Min Boetto
  do i = 1, NXf
    do j = 1, NYf    
		SPhase= nint(SP(Nint((i-1)/4.+1), Nint((j-1)/4.+1))/100.)

		If (SPhase >100) SPhase = 100
		if (SPhase <1) SPhase = 1
		if (CI_0(SPhase,DI(i,j)) == -9999) then
			CI_0(SPhase,DI(i,j)) = CH1_in(i,j)
		Else if (DI(i,j) /= 1) then	
			if (CH1_in(i,j)>CI_0(SPhase,DI(i,j))) CI_0(SPhase,DI(i,j)) = CH1_in(i,j)
		end if
		
		If(DI(i,j)==1) then
			If (CI_cs(SPhase)==-9999) then
				CI_cs(SPhase)=CH1_in(i,j)
			else
				If(alt(i,j)>0.) then
					If (CH1_in(i,j)<CI_cs(SPhase)) CI_cs(SPhase)= CH1_in(i,j)
				end if
			end if
			
			! Boetto es gay
			If (CI_cs_pixel(i,j,Sphase)==-9999) then
				CI_cs_pixel(i,j,Sphase)=CH1_in(i,j)
			else
				If(CH1_in(i,j)<CI_cs_pixel(i,j,Sphase)) CI_cs_pixel(i,j,Sphase) = CH1_in(i,j)
			end if
		end if
	end do
 end do
 
 Do j=1,100
	If ( .not.(CI_0(j,2)==-9999 .and. CI_0(j,3)==-9999 .and. CI_0(j,4)==-9999)) CI_0(j,1)=Maxval(CI_0(j,2:4))
 end do	
 CI_0i=-9999
 CI_csi=-9999
 
 Do i= 1,100
	if (CI_cs(i)==-9999)CI_cs(i)=9999
 end do
 
 Do j= DeltaSP+1,100-DeltaSP
	Do k=1,4
		if (k==1) CI_csi(j) = minval(CI_cs((j-DeltaSP):(j+DeltaSP)))
		CI_0i(j,k) = maxval(CI_0((j-DeltaSP):(j+DeltaSP),K))
	end do
 end do
 
 CI_0=CI_0i
 CI_cs= CI_csi
 Ponderador=1
 do j= 1,100
	if (CI_0(j,2)/=-9999 .and. CI_0(j,3)/=-9999 .and. CI_0(j,4)/=-9999) then
		IC_max=maxval(CI_0(j,2:4))
		Ponderador(j,1) =1
		Do i=2,4
			Ponderador (j,i)=CI_0(j,i)/IC_max
		end do
	end if
 end do
 
 
 !Calculo del indice de cobertura de nubes
 
		
			
			
		
			
		
! /Prueba datos Boetto			




 ! Guardado de matriz recortada y revisada
 call check( nf90_put_var(ncid,sp_varid, SP, start = start,count = count) )
 call check( nf90_put_var(ncid,CH4_varid, CH4_in, start = start,count = count) )
 call check( nf90_put_var(ncid,CH1_varid, CH1_in, start, countf))
 call check( nf90_close(ncid_in) )
 
 write (*,*) '   ',cdia,'  ', hora,':', minu
 
 Go to 100
 
200 Format ('   Tiempo procesamiento: ',I3 'min., ',I3, 'sec.')

 
 999 Stop
 
 contains
 subroutine check(status)
  integer, intent ( in) :: status
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop 2
  end if
 end subroutine check

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
implicit none

real, intent(in):: day          !Day (dd)
real, intent(in) :: month        !Month (mm)
real, intent(in) :: year         !Year (yyyy)
real, intent(out) :: dayj 		!Day of year
integer :: i            			!Index,variable
integer :: leap_day     			!Extra day for leap year
 
! Check for leap year, and add extra day if necessary
IF ( mod(year,400.) == 0 ) THEN
    leap_day = 1    ! Years divisible by 400 are leap years
ELSE IF ( mod(year,100.) == 0 ) THEN
    leap_day = 0    ! Other centuries are not leap years
ELSE IF ( mod(year,4.) == 0 ) THEN
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
 implicit none
 Real, intent(in) :: year
 Real, intent(in) :: day
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


subroutine LDA (NXf, Nyf, Nx, Ny,nboe, CH1, CH4, SP, mes, fboe, dboe, DI ) 

implicit none

integer, intent (in) :: Nxf, Nyf, Nx, Ny, mes, nboe
real, intent (in), dimension (NXf,NYf) :: CH1
real, intent (in), dimension (NX,NY) :: CH4, SP
real, intent (in), dimension (nboe, 4) :: fboe
integer, intent (in), dimension (nboe) :: dboe
integer, intent (out), dimension (NXf,NYf) :: DI
real, dimension (:,:), allocatable ::f, fc
integer, dimension (:), allocatable :: d, dc


integer :: i, j, n=0, nc=0,k, maxCH1, minCH1, maxCH4, minCH4
real :: meanSP, deltaSP=5.
real, dimension(NXf*NYf, 2) :: Xt
integer, dimension (NXf*NYf) :: DI2


n= count(fboe(:,4) = mes)

allocate (f(n,4))
allocate ( d(n))

k = 0
do i= 1, nboe
	if (fboe(i,4)==mes) then
		k= k+1
		f(k,:) = fboe(i,:)
		d(k) = dboe(i)
	end if
end do

if (n/=k) print *,' n y k distintos!', n, k

meanSP = sum(SP)/real(Nx * Ny)
do i= 1, n
	if ((f(i,3)<=(meanSP+deltaSP)) .and.( f(i,3)>=(meanSP-deltaSP))) then 
		nc=nc+1
	end if
end do

allocate ( fc(nc,4))
allocate ( dc(nc))

k = 0

do i= 1, n
	if (f(i,3)<=(meanSP+deltaSP) .and. f(i,3)>=(meanSP-deltaSP)) then
		K = k+1
		fc(k,:) = f(i,:)
		dc(k) = d(i)
	end if
end do
if (nc/=k) print *,' nc y k distintos!', nc, k

k=0
do i= 1, Nxf
	do j = 1, Nyf
	k=k+1	
		Xt(k,1)= CH1(i,j)
		Xt(k,2) = Ch4(Nint((i-1)/4.+1), Nint((j-1)/4.+1))
	end do
end do

maxCH1= Max(Maxval(Xt(:,1)),Maxval(fc(:,1)))
minCH1= Min(MinVAL(Xt(:,1)),Minval(fc(:,1)))
maxCH4= Max(MaxVAL(Xt(:,2)),Maxval(fc(:,2)))
minCH4= Min(Minval(Xt(:,2)),Minval(fc(:,2)))

fc(:,1) = (fc(:,1)-minCH1)/real(MaxCh1-MinCH1)
fc(:,2) = (fc(:,2)-minCH4)/real(MaxCh4-MinCH4)

Xt(:,1) = (fc(:,1)-minCH1)/real(MaxCh1-MinCH1)
Xt(:,2) = (fc(:,2)-minCH4)/real(MaxCh4-MinCH4)

call LDA2 ( nc, NXf,NYf,fc(:,1:2), dc(:), Xt, DI2)

end subroutine LDA

subroutine LDA2 ( nc, NXf,NYf,fc, dc, Xt, DI2)
 
implicit none

integer, intent(in) :: nc, NXf,NYf
real,dimension(NXf*NYf,2),intent(in) :: Xt
real,dimension(nc,2),intent(in) :: fc
integer,dimension(nc),intent(in) :: dc
integer,dimension(NXf*NYf),intent(out) :: DI2
real,dimension(NXf*NYf,4) :: D
integer, parameter :: K=4
integer ::  i, nn, jj, j, errorflag
integer, dimension (:), allocatable :: ii
real, dimension (:,:), allocatable :: Xk
real, dimension (k) :: L, P=0
real, dimension (2,2) :: cw , ck, cw1
real, dimension (2,k) :: xkm
real, dimension (2,1) :: C1=0
real, dimension (NXf*NYf,1) :: C2, VU

VU=1.
D=0.
L=0.
cw=0.
xkm =0.
do i = 1, k
	nn=count(dc=i)
	allocate (ii (nn))
	allocate (XK (nn,2))
	
	jj=1
	do j =1, nc
		if (dc(j)==i) then
			ii(jj)=j
			Xk(jj,:) = fc (j,:)
			jj=jj+1
		end if
	end do
	L(i)=nn

	Xkm(1,i) = sum(XK(:,1))/real(size(XK(:,1)))
	Xkm(2,i) = sum(XK(:,2))/real(size(XK(:,2)))
	!Ck=cov(Xk) ! CovarianzaQ"!! (2x2) 
	Ck=0
	Do	j =1, nn
		Ck(1,1) = Ck(1,1) + ((Xk(j,1)-Xkm(1,i))**2)
		Ck(1,2) = Ck(1,2) + ((Xk(j,1)-Xkm(1,i))*(Xk(j,2)-Xkm(1,i)))
		Ck(2,1) = Ck(2,1) + ((Xk(j,2)-Xkm(1,i))*(Xk(j,1)-Xkm(1,i)))
		Ck(2,2) = Ck(2,2) + ((Xk(j,2)-Xkm(1,i))**2)
	end do
	Ck(1,1) = Ck(1,1)/(nn-1)
	Ck(1,2) = Ck(1,2)/(nn-1)
	Ck(2,1) = Ck(2,1)/(nn-1)
	Ck(2,2) = Ck(2,2)/(nn-1)
	
	Cw=cw+ck*(nn-1)
	P(i)= nn/real(nc)
end do
Cw = Cw/real(nc-k)
call FINDInv(Cw, Cw1, 2, errorflag)
if (errorflag/=0) print *, 'Error de inversion de matriz'

do i= 1,k
	C1(:,1) = matmul(Cw1,Xkm(:,i))	
	C2 = (matmul(reshape(-0.5*Xkm(:,i),(/1,2/)),C1)+Log(P(i)))*VU
	D(:,i) = reshape(matmul(Xt,C1)+C2,(/NXf*NYf/)) 
end do


end subroutine LDA2

SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
	IMPLICIT NONE
	!Declarations
	INTEGER, INTENT(IN) :: n
	INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
	REAL, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
	REAL, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
	
	LOGICAL :: FLAG = .TRUE.
	INTEGER :: i, j, k, l
	REAL :: m
	REAL, DIMENSION(n,2*n) :: augmatrix !augmented matrix
	
	!Augment input matrix with an identity matrix
	DO i = 1, n
		DO j = 1, 2*n
			IF (j <= n ) THEN
				augmatrix(i,j) = matrix(i,j)
			ELSE IF ((i+n) == j) THEN
				augmatrix(i,j) = 1
			Else
				augmatrix(i,j) = 0
			ENDIF
		END DO
	END DO
	
	!Reduce augmented matrix to upper traingular form
	DO k =1, n-1
		IF (augmatrix(k,k) == 0) THEN
			FLAG = .FALSE.
			DO i = k+1, n
				IF (augmatrix(i,k) /= 0) THEN
					DO j = 1,2*n
						augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
					END DO
					FLAG = .TRUE.
					EXIT
				ENDIF
				IF (FLAG .EQV. .FALSE.) THEN
					PRINT*, "Matrix is non - invertible"
					inverse = 0
					errorflag = -1
					return
				ENDIF
			END DO
		ENDIF
		DO j = k+1, n			
			m = augmatrix(j,k)/augmatrix(k,k)
			DO i = k, 2*n
				augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
			END DO
		END DO
	END DO
	
	!Test for invertibility
	DO i = 1, n
		IF (augmatrix(i,i) == 0) THEN
			PRINT*, "Matrix is non - invertible"
			inverse = 0
			errorflag = -1
			return
		ENDIF
	END DO
	
	!Make diagonal elements as 1
	DO i = 1 , n
		m = augmatrix(i,i)
		DO j = i , (2 * n)				
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		END DO
	END DO
	
	!Reduced right side half of augmented matrix to identity matrix
	DO k = n-1, 1, -1
		DO i =1, k
		m = augmatrix(i,k+1)
			DO j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
			END DO
		END DO
	END DO				
	
	!store answer
	DO i =1, n
		DO j = 1, n
			inverse(i,j) = augmatrix(i,j+n)
		END DO
	END DO
	errorflag = 0
END SUBROUTINE FINDinv


 End Program ProcesamientoImagenesMediaHora
