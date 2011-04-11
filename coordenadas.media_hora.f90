Program coordenadas

integer, parameter :: NLat= 432, Nlon= 2255
integer, parameter :: r4i = 1095, r4f = 1368 , r1i = (r4i-1)*4+1, r1f = (r4f-1)*4+1 ! CH1: 4377 a 5469
integer, parameter :: NX = (r4f- r4i)+1, NY = 430 , NXf = (r4f- r4i)*4+1, NYf = (Ny-1)*4+1
integer, dimension (Nx, Ny) :: Lat_CH4 , Lon_CH4
integer, dimension (NXf,NYf) :: Lat_CH1 , Lon_CH1
integer, dimension (Nlon,NLat) :: Latitud , Longitud
Integer :: L1, L2, errorread
Real :: mm

 open (unit=12, file='CH4.latitude.media_hora.txt', status= 'old', Action='read', IOSTAT=errorread )
 If(errorread/=0) print *,' Error de lectura de archivo: ', '  CH4.latitude.media_hora.txt'
 read (12,*) Latitud
 Close (12)
 open (unit=12, file='CH4.longitude.media_hora.txt', status= 'old', Action='read', IOSTAT=errorread )
 If(errorread/=0) print *,' Error de lectura de archivo: ', '  CH4.longitude.media_hora.txt'
 read (12,*) Longitud
 Close (12)

 print *, 'Inicio'
 
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
 
 open (unit=12, file='latitude.CH4.media_hora.txt')
 write (12,*) Lat_CH4
 Close (12)
 open (unit=12, file='longitude.CH4.media_hora.txt')
 write (12,*) Lon_CH4
 Close (12)
 
 open (unit=12, file='latitude.CH1.media_hora.txt')
 write (12,*) Lat_CH1
 Close (12)
 open (unit=12, file='longitude.CH1.media_hora.txt')
 write (12,*) Lon_CH1
 close (12)
 
 print *, 'Terminado..'
 
! print *
! print *, 'CH4',Lat_CH4 (1,1) , Lat_CH4 (2,1), Lat_CH4 (Nx-1,1), Lat_CH4 (Nx,1)
! print *, Lat_CH1 (1,1) , Lat_CH1 (2,1) ,  Lat_CH1 ((Nx-1)*4+1-1,1) , Lat_CH1 ((Nx-1)*4+1,1) 
! print *, Lat_CH1 (1,2) , Lat_CH1 (2,2) ,  Lat_CH1 ((Nx-1)*4+1-1,2) , Lat_CH1 ((Nx-1)*4+1,2)  
! print *, Lat_CH1 (1,3) , Lat_CH1 (2,3) ,  Lat_CH1 ((Nx-1)*4+1-1,3) , Lat_CH1 ((Nx-1)*4+1,3)
! print *, Lat_CH1 (1,4) , Lat_CH1 (2,4) ,  Lat_CH1 ((Nx-1)*4+1-1,4) , Lat_CH1 ((Nx-1)*4+1,4) 
! print *, 'CH4',Lat_CH4 (1,2) , Lat_CH4 (2,2), Lat_CH4 (Nx-1,2), Lat_CH4 (Nx,2) 
! print *,'CH4', Lat_CH4 (1,Ny-1) , Lat_CH4 (2,Ny-1), Lat_CH4 (Nx-1,Ny-1), Lat_CH4 (Nx,Ny-1)
! print *, Lat_CH1(1,(Ny-1)*4+1-2),Lat_CH1 (2,(Ny-1)*4+1-2),Lat_CH1 ((Nx-1)*4+1-1,(Ny-1)*4+1-2),Lat_CH1 ((Nx-1)*4+1,(Ny-1)*4+1-2)
! print *, Lat_CH1(1,(Ny-1)*4+1-1),Lat_CH1 (2,(Ny-1)*4+1-1),Lat_CH1 ((Nx-1)*4+1-1,(Ny-1)*4+1-1),Lat_CH1 ((Nx-1)*4+1,(Ny-1)*4+1-1) 
! print *, Lat_CH1 (1,(Ny-1)*4+1) , Lat_CH1 (2,(Ny-1)*4+1) , Lat_CH1 ((Nx-1)*4+1,(Ny-1)*4+1) , Lat_CH1 ((Nx-1)*4+1,(Ny-1)*4+1)
! print *,'CH4', Lat_CH4 (1,Ny) , Lat_CH4 (2,Ny), Lat_CH4 (Nx-1,Ny), Lat_CH4 (Nx,Ny)
! print *
! print *, 'CH4',Lon_CH4 (1,1) , Lon_CH4 (2,1), Lon_CH4 (Nx-1,1), Lon_CH4 (Nx,1)
! print *, Lon_CH1 (1,1) , Lon_CH1 (2,1) ,  Lon_CH1 ((Nx-1)*4+1-1,1) , Lon_CH1 ((Nx-1)*4+1,1) 
! print *, Lon_CH1 (1,2) , Lon_CH1 (2,2) ,  Lon_CH1 ((Nx-1)*4+1-1,2) , Lon_CH1 ((Nx-1)*4+1,2)  
! print *, Lon_CH1 (1,3) , Lon_CH1 (2,3) ,  Lon_CH1 ((Nx-1)*4+1-1,3) , Lon_CH1 ((Nx-1)*4+1,3) 
! print *, Lon_CH1 (1,4) , Lon_CH1 (2,4) ,  Lon_CH1 ((Nx-1)*4+1-1,4) , Lon_CH1 ((Nx-1)*4+1,4) 
! print *, 'CH4',Lon_CH4 (1,2) , Lon_CH4 (2,2), Lon_CH4 (Nx-1,2), Lon_CH4 (Nx,2)
! print *, 'CH4',Lon_CH4 (1,Ny-1) , Lon_CH4 (2,Ny-1), Lon_CH4 (Nx-1,Ny-1), Lon_CH4 (Nx,Ny-1)
! print *,Lon_CH1(1,(Ny-1)*4+1-2),Lon_CH1 (2,(Ny-1)*4+1-2),Lon_CH1((Nx-1)*4+1-1,(Ny-1)*4+1-2),Lon_CH1((Nx-1)*4+1,(Ny-1)*4+1-2)
! print *,Lon_CH1(1,(Ny-1)*4+1-1),Lon_CH1 (2,(Ny-1)*4+1-1),Lon_CH1((Nx-1)*4+1-1,(Ny-1)*4+1-1),Lon_CH1 ((Nx-1)*4+1,(Ny-1)*4+1-1) 
! print *, Lon_CH1 (1, (Ny-1)*4+1) , Lon_CH1 (2, (Ny-1)*4+1),Lon_CH1 ((Nx-1)*4+1,(Ny-1)*4+1),Lon_CH1 ((Nx-1)*4+1, (Ny-1)*4+1)
! print *, 'CH4',Lon_CH4 (1,Ny) , Lon_CH4 (2,Ny), Lon_CH4 (Nx-1,Ny), Lon_CH4 (Nx,Ny)
  
end program 
