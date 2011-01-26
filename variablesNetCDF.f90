! Copyright (C) 2010, Juan Pablo Justiniano  <jpjustiniano@gmail.com>
! Programa para el almacenaje de variables climatologicas en formato NetCDF
  
 program variables
   use netcdf
   implicit none
 
   character(*), parameter :: FILE_NAME = "variables.nc" ! This is the name of the data file
   integer :: ncid
   
   integer, parameter :: NDIMS = 2	      ! We are writing 2D data.
   !integer, parameter :: NX = 1093, NY = 1717
   integer, parameter :: r4i = 1095, r4f = 1368 , r1i = r4i*4, r1f = r4f*4 !
   integer, parameter :: NLat= 432, Nlon= 2255
   integer, parameter :: NX=(r4f- r4i)+1, NY=430, NXf=(r4f-r4i)*4+1, NYf=(Ny-1)*4+1
   integer, dimension (NX,NY) :: Altura , Temperatura, HR, Visibilidad
   integer, dimension ((r4f- r4i)*4+1, (Ny-1)*4+1) :: Lat_CH1 , Lon_CH1
   integer, dimension (Nx, Ny)  :: Lat_CH4 , Lon_CH4
   integer, dimension (Nlon,NLat) :: Latitud , Longitud
  
   integer :: x_dimid, y_dimid
   integer :: x_varid, y_varid       
   integer :: Alt_varid, Temp_varid, HR_varid, Vis_varid
   integer :: dimids(NDIMS), dimids3d(3), rec(3)
   
   integer :: i=1,j=1, lat, lon, k
   Integer :: L1, L2
   Integer :: ierror
   Real :: mm, al
   Character(2) :: ci
   Character(60) ::filename
	   
   !*********************************************************************** Fin declaracion Variables
   call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) ) ! f90_clobber overwrite this file if it already exists.
 
   call check( nf90_def_dim(ncid, "x", NX, x_dimid) )  ! Define the dimensions. 
   call check( nf90_def_dim(ncid, "y", NY, y_dimid) )
   
   dimids =  (/ x_dimid, y_dimid/) ! The dimids array is used to pass the IDs of the dimensions
   
   call check( nf90_def_var(ncid, "Alt", NF90_INT, dimids, Alt_varid) )  ! Define the variable and type. NF90_INT (4-byte integer)
   call check( nf90_def_var(ncid, "Temp", NF90_INT, dimids, Temp_varid) )
   call check( nf90_def_var(ncid, "HR", NF90_INT, dimids, HR_varid) )
   call check( nf90_def_var(ncid, "Vis", NF90_INT, dimids, Vis_varid) )
   
   call check( nf90_put_att(ncid, Alt_varid, "units", "meters") )
   call check( nf90_put_att(ncid, Alt_varid, "long_name", "Altitud") ) 
   
   call check( nf90_put_att(ncid, Temp_varid, "units", "C Deg") ) 
   call check( nf90_put_att(ncid, Temp_varid, "long_name", "Temperatura") ) 
   
   call check( nf90_put_att(ncid, HR_varid, "units", "percent") )  
   call check( nf90_put_att(ncid, HR_varid, "long_name", "Humedad Relativa") )  
   
   call check( nf90_put_att(ncid, Vis_varid, "units", "meters") )  
   call check( nf90_put_att(ncid, Vis_varid, "long_name", "Visibilidad") ) 
	
   call check( nf90_enddef(ncid) )  ! End define metadata mode. 
   !*********************************************************************** Fin Definiciones
   !************* Latitud y longitud
    count = (/ NX, NY, 1 /)
	start = (/ 1, 1, 1 /)
    
    
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
       !************* /Latitud y longitud
       
       !************* Altitud
       Open(10,FILE='./Altura/Altura.media_hora.txt', IOSTAT=ierror)
        Do i=1, NXf
			Do j = 1, NYf
				read(10,*, IOSTAT=ierror) al
				Altura(i,j) = Nint(al*100)
			End do
		End do
		close (10)

       call check( nf90_put_var(ncid, Alt_varid, Altura) )
       !************* /Altitud
       
       !************* Temperatura
	Do k = 1, 12
	   start(3) = k 
	   
	   If ( k<10) write (ci,'(I1)') k 
	   If ( k>=10) write (ci,'(I2)') k 
	   
       filename= './Temperetura/Temperatura.'//ci//'.media_hora.txt'
       Open(10,FILE=filename, IOSTAT=ierror)
       Do  i=1, NXf
          Do j = 1, NYf
			read(10,*, IOSTAT=ierror) al
            Temperatura(i,j) = Nint(al*100)
          end do
       end do
       close (10)
       call check( nf90_put_var(ncid, temp_varid, Temperatura, start, count) )
	End do
       !************* /Temperatura
       
       !************* Visibilidad
       Open(10,FILE='Visibilidad.media_hora.txt', IOSTAT=ierror)
       Do  i=1, NXf
          Do j = 1, NYf
			read(10,*, IOSTAT=ierror) al
             Visibilidad (i,j) = al
          end do
       end do
       
       close (10)
       call check( nf90_put_var(ncid, Vis_varid, Visibilidad) )
       !************* /Visibilidad       
       
       !************* HR
       Open(10,FILE='HR.media_hora.txt', IOSTAT=ierror)
       Do  i=1, NXf
          Do j = 1, NYf
			read(10,*, IOSTAT=ierror) al
             HR (i,j) = al
          end do
       end do
       
       close (10)
       call check( nf90_put_var(ncid, Vis_varid, HR) )
       !************* /HR  
		
       call check( nf90_close(ncid) ) ! Close the file.
     
       print *, '  SUCCESS writing ', FILE_NAME, '!'
     
     contains
       subroutine check(status)
         integer, intent ( in) :: status
     
         if(status /= nf90_noerr) then
           print *, trim(nf90_strerror(status))
           stop 2
         end if
       end subroutine check
     end program variables
