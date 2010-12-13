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
       integer, parameter :: NX = 10, NY = 15
       real :: x(NX), y(NY)
       real :: Start_Lat= -17.5,  Start_Lon= -74., Dlat =-0.0400, Dlon=0.0300
       integer :: x_dimid, y_dimid
       integer :: x_varid, y_varid
       
       integer, dimension(:,:), allocatable :: Alt, Temp, HR, Vis
       integer :: Alt_varid, Temp_varid, HR_varid, Vis_varid
       integer :: dimids(NDIMS)
       
       integer :: i=1,j=1, lat, lon
           
       allocate(Alt (NY, NX)) ! Allocate memory for data.
       allocate(Temp (NY, NX))
       allocate(HR (NY, NX))
       allocate(Vis (NY, NX))
       !*********************************************************************** Fin declaracion Variables
       call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) ) ! f90_clobber overwrite this file if it already exists.
     
       call check( nf90_def_dim(ncid, "x", NX, x_dimid) )  ! Define the dimensions. 
       call check( nf90_def_dim(ncid, "y", NY, y_dimid) )
     
       call check( nf90_def_var(ncid,"x", NF90_REAL, x_dimid, x_varid) )	! Define the coordinate variables
       call check( nf90_def_var(ncid, "y", NF90_REAL, y_dimid, y_varid) )
       
       dimids =  (/ y_dimid, x_dimid/) ! The dimids array is used to pass the IDs of the dimensions
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
       print *, "1"
       do lat = 1, NY
          Y(lat) = START_LAT + (lat - 1) * (Dlat)
       end do
       do lon = 1, NX
          X(lon) = START_LON + (lon - 1) * (Dlon)
       end do
       call check( nf90_put_var(ncid, x_varid, x) )
       call check( nf90_put_var(ncid, y_varid, y) )
       !************* /Latitud y longitud
       
       !************* Altitud
!       Open(10,FILE='Chile_recortado.dat', IOSTAT=ierror)
!       Do
!		read(10,*, IOSTAT=ierror) xx, yy, al
!		if (ierror/=0) then
!			write(*,*) xx, yy, al
!			Exit
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
