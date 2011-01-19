Program coordenadas

integer, parameter :: NLat= 1330, Nlon= 2565
integer, parameter :: r4i = 1137, r4f = 1441 , r1i = r4i*4, r1f = r4f*4 ! CH1: 4380 a 5472
integer, parameter :: NX = (r4f- r4i)+1, NY = 479 
integer, dimension (Nx, Ny) :: Lat_CH4 , Lon_CH4
integer, dimension ((r4f- r4i)*4+1, (Ny-1)*4+1) :: Lat_CH1 , Lon_CH1
integer, dimension (Nlon,NLat) :: Latitud , Longitud
Integer :: L1, L2
Real :: mm

 open (unit=12, file='latitude.CH4.south_full.txt', status= 'old', Action='read' )
 read (12,*) Latitud
 Close (12)
 open (unit=12, file='longitude.CH4.south_full.txt', status= 'old', Action='read' )
 read (12,*) Longitud
 Close (12)

 print *, 'Inicio'

 Lat_CH4 = Latitud(r4i:r4f,690:1168)
 Lon_CH4 = Longitud(r4i:r4f,690:1168)
 
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
 
 open (unit=12, file='latitude.CH1.south_full.txt')
 write (12,*) Lat_CH1
 Close (12)
 open (unit=12, file='longitude.CH1.south_full.txt')
 write (12,*) Lon_CH1
 Close (12)
 

 end program 
