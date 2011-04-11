! Copyright (C) 2011, Juan Pablo Justiniano  <jpjustiniano@gmail.com>
! Programa para el almacenaje de variables climatologicas en formato NetCDF
  
 program variables
   use netcdf
   implicit none
 
   character(*), parameter :: FILE_NAME = "variables.nc" ! This is the name of the data file
   integer :: ncid
   
   integer, parameter :: NDIMS = 2	      
   !integer, parameter :: NX = 1093, NY = 1717
   integer, parameter :: r4i = 1095, r4f = 1368 , r1i = (r4i-1)*4+1, r1f = (r4f-1)*4+1 ! CH1: 4377 a 5469
   integer, parameter :: NLat= 432, Nlon= 2255
   integer, parameter :: NX=(r4f- r4i)+1, NY=430, NXf=(r4f-r4i)*4+1, NYf=(Ny-1)*4+1
   integer, dimension (NXf,NYf) :: Altura , Temperatura, HR, Visibilidad, Albedo
   integer, dimension (NXf,NYf) :: Lat_CH1 , Lon_CH1
   integer, dimension (Nx, Ny)  :: Lat_CH4 , Lon_CH4
   integer, dimension (Nlon,NLat) :: Latitud , Longitud
   integer :: Lat_CH4_varid, Lon_CH4_varid , Lat_CH1_varid, Lon_CH1_varid
  
   integer :: x_dimid, y_dimid, mes_dimid
   integer :: x_varid, y_varid, mes_varid       
   integer :: Alt_varid, Temp_varid, HR_varid, Vis_varid, Albedo_varid
   integer :: dimids(NDIMS), dimids3d(3), rec(3)
   
   integer :: i=1,j=1, lat, lon, k, count(3), start(3)
   Integer :: L1, L2
   Integer :: ierror, al, errorread
   real :: mm
   Character(1) :: c1i
   Character(2) :: ci
   Character(60) ::filename
	   
   !*********************************************************************** Fin declaracion Variables

   call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) ) ! f90_clobber overwrite this file if it already exists.
 
   call check( nf90_def_dim(ncid, "x", NXf, x_dimid) )  ! Define the dimensions. 
   call check( nf90_def_dim(ncid, "y", NYf, y_dimid) )
   call check( nf90_def_dim(ncid, "mes", 12, mes_dimid) )
   
   dimids =  (/ x_dimid, y_dimid/)
   dimids3d = (/ x_dimid, y_dimid, mes_dimid  /)
   
   call check( nf90_def_var(ncid, "Alt", NF90_SHORT, dimids, Alt_varid) )  ! Define the variable and type. NF90_INT (4-byte integer)
   call check( nf90_def_var(ncid, "Temp", NF90_SHORT, dimids, Temp_varid) )
   call check( nf90_def_var(ncid, "HR", NF90_SHORT, dimids, HR_varid) )
   call check( nf90_def_var(ncid, "Albedo", NF90_SHORT, dimids, Albedo_varid) )
   call check( nf90_def_var(ncid, "Lat_CH4", NF90_SHORT, dimids, Lat_CH4_varid) )
   call check( nf90_def_var(ncid, "Lon_CH4", NF90_SHORT, dimids, Lon_CH4_varid) )
   call check( nf90_def_var(ncid, "Lat_CH1", NF90_SHORT, dimids, Lat_CH1_varid) )
   call check( nf90_def_var(ncid, "Lon_CH1", NF90_SHORT, dimids, Lon_CH1_varid) )
   
   call check( nf90_put_att(ncid, Alt_varid, "units", "meters") )
   call check( nf90_put_att(ncid, Alt_varid, "long_name", "Altitud") )
   call check( nf90_put_att(ncid, Alt_varid, "scale_factor", 1) ) 
   call check( nf90_put_att(ncid, Alt_varid, "_CoordinateAxes", "Lat_CH1 Lon_CH1") )
   
   call check( nf90_put_att(ncid, Temp_varid, "units", "C Deg") ) 
   call check( nf90_put_att(ncid, Temp_varid, "long_name", "Temperatura") ) 
   call check( nf90_put_att(ncid, Temp_varid, "scale_factor", 0.01) ) 
   call check( nf90_put_att(ncid, Temp_varid, "_CoordinateAxes", "Lat_CH1 Lon_CH1") )
   
   call check( nf90_put_att(ncid, HR_varid, "units", "percent") )  
   call check( nf90_put_att(ncid, HR_varid, "long_name", "Humedad Relativa") )
   call check( nf90_put_att(ncid, HR_varid, "scale_factor", 0.01) ) 
   call check( nf90_put_att(ncid, HR_varid, "_CoordinateAxes", "Lat_CH1 Lon_CH1") )  
      
   call check( nf90_put_att(ncid, Albedo_varid, "units", "%") )  
   call check( nf90_put_att(ncid, Albedo_varid, "long_name", "Albedo") ) 
   call check( nf90_put_att(ncid, Albedo_varid, "scale_factor", 0.01) ) 
   call check( nf90_put_att(ncid, Albedo_varid, "_CoordinateAxes", "Lat_CH1 Lon_CH1") )
   
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
	
   call check( nf90_enddef(ncid) )  ! End define metadata mode. 
   !*********************************************************************** Fin Definiciones
   !************* Latitud y longitud
	print *

	print *,'  Geolocalizacion'
    count = (/ NXf, NYf, 1 /)
	start = (/ 1, 1, 1 /)
    
! Matrices de Geolocalizacion 
 open (unit=12, file='latitude.CH4.media_hora.txt', status= 'old', Action='read', IOSTAT=errorread )
 If(errorread/=0) print *,' Error de lectura de archivo: ', 'latitude.CH4.media_hora.txt'
 read (12,*) Latitud
 Close (12)
 open (unit=12, file='longitude.CH4.media_hora.txt', status= 'old', Action='read', IOSTAT=errorread )
 If(errorread/=0) print *,' Error de lectura de archivo: ', 'longitude.CH4.media_hora.txt'
 read (12,*) Longitud
 Close (12)
 
 Lat_CH4 = Latitud(r4i:r4f,1:430)
 Lon_CH4 = Longitud(r4i:r4f,1:430)
 
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
 

    call check( nf90_put_var(ncid, Lat_CH4_varid, Lat_CH4) )
    call check( nf90_put_var(ncid, Lon_CH4_varid, Lon_CH4) )
    call check( nf90_put_var(ncid, Lat_CH1_varid, Lat_CH1) )
    call check( nf90_put_var(ncid, Lon_CH1_varid, Lon_CH1) )
 
! Fin de Geolocalizacion 
!************* /Latitud y longitud
   
   
!************* Altitud
print *,'  Altitud'
	Open(10,FILE='./Altura/Altura.media_hora.txt', IOSTAT=ierror, status= 'old', Action='read') ! (m)
	If (ierror /= 0) Then 
		print *, '   Erorr en leer archivo : ./Altura/Altura.media_hora.txt'
		goto 999
	End if
	read (10,*, IOSTAT=ierror) ((Altura(i,j), j = 1, NYf), i=1, NXf)
	close (10)

   call check( nf90_put_var(ncid, Alt_varid, Altura) )
!************* /Altitud
   
!************* Temperatura
print *,'  Temperatura'
Do k = 1, 12
   start(3) = k 
	  
	If ( k .le. 9) then
		write (c1i,'(I1)') k 
		ci=  '0'//c1i
	else
		write (ci,'(I2)') k
   end if 
   filename= './Temperatura/Temperatura.'//ci//'.media_hora.txt'
   
   Open(10,FILE=filename, IOSTAT=ierror, status= 'old', Action='read') !(c*100)
	If (ierror/=0) Then
		print *,'No se puede abrir:',filename,k
		goto 999
	End if
	
	read(10,*, IOSTAT=ierror) ((Temperatura(i,j), j = 1, NYf), i=1, NXf)
	close (10)

	call check( nf90_put_var(ncid, Temp_varid, Temperatura, start, count) )	
End do
!************* /Temperatura
   
!************* Visibilidad
!print *,'  Visivilidad'
!Do k = 1, 12
!   start(3) = k 
	      
!	   If ( k<10) then
!	    write (c1i,'(I1)') k 
!	    filename= './Visibilidad/Visibilidad.'//c1i//'.media_hora.txt'
!	   else
!	   	write (ci,'(I2)') k 
!	    filename= './Visibilidad/Visibilidad.'//ci//'.media_hora.txt'
!	   end if 
       
!       Open(10,FILE=filename, IOSTAT=ierror, status= 'old', Action='read')
!		If (ierror/=0) Then
!			print *,'No se puede abrir: ',filename,k
!			goto 999
!		End if
		
!		read(10,*, IOSTAT=ierror) ((Visibilidad(i,j), j = 1, NYf), i=1, NXf)
!		close (10)
!Do i =1, NXf
!	Do j = 1, Nyf
!		Visibilidad(i,j) = 10
!	end do
!end do
		!call check( nf90_put_var(ncid, Vis_varid, Visibilidad, start, count) )	
!	End do
!************************** /Visibilidad       
       
!**************************  HR
	print *,'  HR'
	Do k = 1, 12
	   start(3) = k 
	      
		If ( k .le. 9) then
			write (c1i,'(I1)') k 
			ci=  '0'//c1i
		else
			write (ci,'(I2)') k 
		end if 
		filename= './HR/HR.'//ci//'.media_hora.txt'
		
		Open(10,FILE=filename, IOSTAT=ierror)
		If (ierror/=0) Then
			print *,'No se puede abrir: ',filename,k
			goto 999
		End if
		
		read(10,*, IOSTAT=ierror) ((HR(i,j), j = 1, NYf), i=1, NXf)
		close (10)

		call check( nf90_put_var(ncid, HR_varid, HR, start, count) )	
	End do
!**************************  /HR  

!**************************  Albedo
	print *,'  Albedo'
	Do k = 1, 12
	   start(3) = k 
	      
		If ( k .le. 9) then
			write (c1i,'(I1)') k 
			ci=  '0'//c1i
		else
			write (ci,'(I2)') k 
		end if 
	    filename= './Albedo/Albedo.'//ci//'.media_hora.txt'
       
       Open(10,FILE=filename, IOSTAT=ierror, status= 'old', Action='read')
		If (ierror/=0) Then
			print *,'No se puede abrir:',filename,k
			goto 999
		End if
		
		read(10,*, IOSTAT=ierror) ((Albedo(i,j), j = 1, NYf), i=1, NXf)
		close (10)
	
		call check( nf90_put_var(ncid, Albedo_varid, Albedo, start, count) )	
	End do
!**************************   /Albeo  
		
	call check( nf90_close(ncid) ) ! Close the file.
     
	print *, '             Terminado exitoso de la escritura de archivo ', FILE_NAME, ' !!'
      
	999 Stop  
     
	contains
       subroutine check(status)
         integer, intent ( in) :: status
     
         if(status /= nf90_noerr) then
           print *, trim(nf90_strerror(status))
           stop 2
         end if
       end subroutine check
      
     end program variables
