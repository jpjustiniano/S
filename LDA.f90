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

 subroutine LDA2
implicit none
integer, intent(in) :: nc, NXf,NYf
real,dimension(NXf*NYf,2),intent(in) :: Xt
real,dimension(nc,2),intent(in) :: fc
integer,dimension(nc),intent(in) :: dc
integer,dimension(NXf*NYf),intent(out) :: DI2

integer :: K=4, i, nn, jj
integer, dimension (:), allocatable :: ii
real, dimension (:,:), allocatable :: Xk
real, dimension (k) :: L, P=0
real, dimension (2,2) :: cw , ck, cw1
real, dimension (2,k) :: xkm

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
Xkm(:,i) = transpose((/sum(XK(:,1))/real(size(XK(:,1)));sum(XK(:,2))/real(size(XK(:,2)))/))
Ck=cov(Xk) ! CovarianzaQ"!! (2x2) 
Cw=cw+ck*(nn-1)
p(i)= nn/real(nc)
end do
Cw = Cw/real(nc-k)
Cw1=inv(Cw) !!inverse





end subroutine LDA2

end Module Subs
