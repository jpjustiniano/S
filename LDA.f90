Module Subs

contains

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


n= count(fboe(:,4) == mes)

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
	nn=count(dc==i)
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

end Module Subs
