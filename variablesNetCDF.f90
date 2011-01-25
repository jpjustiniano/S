!gfortran -Wall -I/usr/include -L/usr/lib -lnetcdff -lnetcdf -o <name> archivo.f90   
! Agregar dimension mes para algunas variables
! count = (/ NLONS, NLATS, NLVLS, 1 /)
! start = (/ 1, 1, 1, 1 /) 
! do rec = 1, NRECS
!   start(4) = rec
!   call check( nf90_put_var(ncid, temp_varid, temp_out, start = start,count = count) )


     program variables
       use netcdf
       implicit none
     
       character(*), parameter :: FILE_NAME = "variables.nc" ! This is the name of the data file
       integer :: ncid
       
       integer, parameter :: NDIMS = 2	      ! We are writing 2D data.
       integer, parameter :: NX = 1093, NY = 1717
      
       integer :: x_dimid, y_dimid
       integer :: x_varid, y_varid
       
       integer, dimension(:,:), allocatable :: Alt, Temp, HR, Vis
       integer :: Alt_varid, Temp_varid, HR_varid, Vis_varid
       integer :: dimids(NDIMS)
       integer, dimension (Nx,Ny) :: Latitud , Longitud
       integer :: i=1,j=1, lat, lon
           
       allocate(Alt (Nx, Ny)) ! Allocate memory for data.
       allocate(Temp (Nx, Ny))
       allocate(HR (Nx, Ny))
       allocate(Vis (Nx, Ny))
       !*********************************************************************** Fin declaracion Variables
       call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) ) ! f90_clobber overwrite this file if it already exists.
     
       call check( nf90_def_dim(ncid, "x", NX, x_dimid) )  ! Define the dimensions. 
       call check( nf90_def_dim(ncid, "y", NY, y_dimid) )
       
       dimids =  (/ x_dimid, y_dimid/) ! The dimids array is used to pass the IDs of the dimensions
       
       call check( nf90_def_var(ncid, "Alt", NF90_INT, dimids, Alt_varid) )  ! Define the variable and type. NF90_INT (4-byte integer)
       call check( nf90_def_var(ncid, "Temp", NF90_INT, dimids, Temp_varid) )
       call check( nf90_def_var(ncid, "HR", NF90_INT, dimids, HR_varid) )
       call check( nf90_def_var(ncid, "Vis", NF90_INT, dimids, Vis_varid) )
       
       call check( nf90_put_att(ncid, x_dimid, "units", "degrees_south") )	! Assign units attributes to coordinate var data.
       call check( nf90_put_att(ncid, x_dimid, "long_name", "longitude") )
       call check( nf90_put_att(ncid, y_dimid, "units", "degrees_west") )
       call check( nf90_put_att(ncid, y_dimid, "long_name", "latitude") )
       call check( nf90_put_att(ncid, Alt_varid, "units", "m*100") )
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
       Open(10,FILE='Chile_recortado.dat', IOSTAT=ierror)
        Do i=1, NXf
			Do j = 1, NYf
				read(10,*, IOSTAT=ierror) al
				Altura(i,j) = al
			End do
		End do
		
						
				
				
				
		!if (ierror/=0) then
		!	write(*,*) xx, yy, al
		!	Exit
!		endif
!		i=nint((74.+ xx)/0.030042918+1)
!		j=nint(650-(43.5 + yy)/0.04006163)
!		Alt(j,i) = nint(al*100)					! Guardado como entero por 100.
!       Enddo    
!       Close(10)
       
!       call check( nf90_put_var(ncid, Alt_varid, Alt) )
       !************* /Altitud
       
       !************* Temperatura
       print *, "2"
       Do lon = 1, NX
          Do lat = 1, NY
             Temp(lon, lat) = nint((30.-lat/100.)*100.)
          end do
       end do
       call check( nf90_put_var(ncid, temp_varid, temp) )
       !************* /Temperatura
       !************* Visibilidad
       print *, "3"
       Do lon = 1, NX
          Do lat = 1, NY
             Temp(lon, lat) = nint(int(10.-lat*40/100.)*100.)
          end do
       end do
       print *, "5"
       call check( nf90_put_var(ncid, Vis_varid, Temp) )
       print *, "6"
       !************* /Visibilidad       
       !************* HR
       
       Do lon = 1, NX
          Do lat = 1, NY
             HR(lon, lat) = nint((0.60-(lat-300.)/1000.-(lon-130.)/1000.)*100)
          end do
       end do
       call check( nf90_put_var(ncid, HR_varid, HR) )
       !************* /HR  
		
       call check( nf90_close(ncid) ) ! Close the file.
     
       print *, "*** SUCCESS writing example file simple_xy.nc! "
     
     contains
       subroutine check(status)
         integer, intent ( in) :: status
     
         if(status /= nf90_noerr) then
           print *, trim(nf90_strerror(status))
           stop 2
         end if
       end subroutine check
     end program variables
