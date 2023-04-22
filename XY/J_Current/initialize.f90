!--- PROJECT-DEPENDENT ---------------------------------------------
!==============Initialization ======================================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE initialize
	implicit none

	!-- order should be changed --------------------------------------
	call tst_and_prt
	call set_RNG
	call def_para
	call alloc_arr
	call def_latt
	call def_prob
	call def_conf

	!-- measurement initialization -----------------------------------
	Obs = 0.d0;   Quan = 0.d0
	Ave = 0.d0;   Dev  = 0.d0;     Cor = 0.d0
	return
END SUBROUTINE initialize

!==============Test and Print ======================================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE tst_and_prt
	implicit none

	!-- Test and Print -----------------------------------------------
	if((NBlck>MxBlck).or.(NBlck<MnBlck)) then
		write(6,*) 'MnBlck <= NBlck <=MxBlck?';             stop
	endif

	!if((NBlck>200).and.(NBlck/=MxBlck)) then
	!	write(6,*) '"NBlck>200" is supposed for extensive &
	!		& simulation. "NBlck=MxBlk" is suggested!';         stop
	!endif

	if( D/=1 .and. D/=2 .and. D/=3 ) stop "Err: D"

	if(Ly>MxLy) then
		write(6,*) 'Ly<=MxLx?';                             stop
	endif

	if( D==1 ) then
		Vol = Lx
		nnb = 2
		write(6,40) Lx
		40 format(' AT on ',i4,2x,'chain')
	else if ( D==2 ) then
		Vol = Lx*Ly
		nnb = 4
		write(6,41) Lx,Ly
		41 format(' AT on ',i4,1x,'x',1x,i4,2x,'square lattice')
	else if ( D==3 ) then
		Vol = Lx*Ly*Lz
		nnb = 6
		write(6,42) Lx,Ly,Lz
		42 format(' AT on ',i4,1x,'x',1x,i4,1x,'x',1x,i4,2x,'cubic lattice')
	end if
	bVol= Vol*nnb/2
	nnb2=nnb/2
	gVol=(Lx/2+1)*(Ly/2+1)*(Lz/2+1)-1

	write(6,51) betaJ
	51 format(' beta*J:',f12.8)

	write(6,52) Nsamp*NBlck
	52 format(' Will simulate      ',i10,2x,'steps ')

	write(6,53) NBlck
	53 format(' #Blocks            ',i10)

	write(6,54) Ntoss*NBlck
	54 format(' Throw away         ',i10,2x,'steps')

	return
END SUBROUTINE tst_and_prt

!============== definite parameter ====================================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE def_para
	IMPLICIT NONE
	integer    :: i

	!------------- 1234567890123456 -------
	! <=16
	ParaName( 1) = 'ExVxbeta'
	ParaName( 2) = 'NxN_N'
	ParaName( 3) = 'windingx'
	ParaName( 4) = 'windingy'
	ParaName( 5) = 'windingz'
	ParaName( 6) = 'winding'
	ParaName( 7) = 'step'
	ParaName( 8) = 'current'

	ParaName( NObs_b+1) = 'spec_heatxVxb'

	open(1,file=trim(paralist),action='write')
	do i=1, NObs
		write(1,'(1X,I2,4X,A)') i, ParaName(i)
	end do
	close(1)
	return
END SUBROUTINE def_para
!===================================================================

!==============definite Lattice ====================================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE alloc_arr
	implicit none

	if(Ly<4) then
		write(6,*) 'L>4?'; stop
	endif

	wv  = 1.d0/Vol

	allocate(Site(Vol))
	allocate(Bond(Vol*nnb/2))
	allocate(sCur(nnb))
	allocate(Ine(0:MaxCur))
	allocate(BesselRat(0:MaxCur,0:MaxCur))
	allocate(gr(0:gVol,NBlck))
	allocate(gri(0:gVol))
	allocate(Ms(nnb,2))

	BesselRat = 0.d0

	return
END SUBROUTINE alloc_arr
!===================================================================

!==============definite Lattice ====================================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE def_latt
	implicit none

	integer           :: Vc
	integer           :: ix, iy, iz

	Lxyz = (/Lx,Ly,Lz/)

	if( D==1 ) then
		do ix=1, Lx
			Vc = ix
			Site(Vc)%ns(1) = ix-1
			Site(Vc)%ns(2) = ix+1
			Site(Vc)%nb(1) = ix-1
			Site(Vc)%nb(2) = ix
		end do
		Site( 1)%ns(1) = Lx
		Site(Lx)%ns(2) = 1
		Site( 1)%nb(1) = Lx
		Ms(1,1) =  1
		Ms(2,1) =  1
		Ms(1,2) = -1
		Ms(2,2) =  1
		sCur(1) =  1
		sCur(2) = -1
	else if( D==2 ) then
		! ------- calculate the neighbors of every site ----------
		! y
		! ^ 13 14 15 16
		! |  9 10 11 12
		! |  5  6  7  8
		! |  1  2  3  4
		! ---------> x
		! --- neighbors ---
		!      4
		!      |
		! 3 -- * -- 1
		!      |
		!      2
		!
		!------- bonds of 4*4 ----------------------
		! *: site,  1,2,...: the numbers of bonds
		!
		!  |     |     |     |
		!  5     6     7     8  + Vol       4
		!  |     |     |     |              |
		!--*--5--*--6--*--7--*--8      3 -- * -- 1
		!  |     |     |     |              |
		!  1     2     3     4  + Vol       2
		!  |     |     |     |
		!--*--1--*--2--*--3--*--4
		do iy=1, Ly
		do ix=1, Lx
			Vc = ix + (iy-1)*Lx
			Site(Vc)%ns(1) = mod(ix,Lx)+1 + (iy-1)*Lx
			Site(Vc)%ns(3) = mod(ix-2+Lx,Lx)+1 + (iy-1)*Lx
			Site(Vc)%ns(2) = ix + (mod(iy-2+Ly,Ly)+1-1)*Lx
			Site(Vc)%ns(4) = ix + (mod(iy,Ly)+1-1)*Lx

			Site(Vc)%nb(1) = Vc
			Site(Vc)%nb(4) = Vc+Vol
		end do
		end do
		open(1, file="bond.list",action='write')
		do iy=1, Ly
		do ix=1, Lx
			Vc = ix + (iy-1)*Lx
			Site(Vc)%nb(3) = Site(Site(Vc)%ns(3))%nb(1)
			Site(Vc)%nb(2) = Site(Site(Vc)%ns(2))%nb(4)
			write(1,*) Site(Vc)%nb(1:4)
		end do
		end do
		close(1)
		Ms(1,1) = 1
		Ms(2,1) = 2
		Ms(3,1) = 1
		Ms(4,1) = 2
		Ms(1,2) = 1
		Ms(2,2) =-1
		Ms(3,2) =-1
		Ms(4,2) = 1
		sCur(1:nnb2)     =  1
		sCur(nnb2+1:nnb) = -1
	else if( D==3 ) then
		!  6
		!  |
		!  5
		do iz=1, Lz
		do iy=1, Ly
		do ix=1, Lx
			Vc = ix + (iy-1)*Lx + (iz-1)*Lx*Ly
			Site(Vc)%ns(1) = mod(ix,Lx)+1 + (iy-1)*Lx + (iz-1)*Lx*Ly
			Site(Vc)%ns(3) = mod(ix-2+Lx,Lx)+1 + (iy-1)*Lx + (iz-1)*Lx*Ly
			Site(Vc)%ns(2) = ix + (mod(iy-2+Ly,Ly)+1-1)*Lx + (iz-1)*Lx*Ly
			Site(Vc)%ns(4) = ix + (mod(iy,Ly)+1-1)*Lx + (iz-1)*Lx*Ly
			Site(Vc)%ns(5) = Vc -(iz-1)*Lx*Ly + (mod(iz-2+Lz,Lz)+1-1)*Lx*Ly
			Site(Vc)%ns(6) = Vc -(iz-1)*Lx*Ly + (mod(iz,Lz)+1-1)*Lx*Ly

			Site(Vc)%nb(1) = Vc
			Site(Vc)%nb(4) = Vc+Vol
			Site(Vc)%nb(6) = Vc+Vol*2
		end do
		end do
		end do
		open(1, file="bond.list",action='write')
		do iz=1, Lz
		do iy=1, Ly
		do ix=1, Lx
			Vc = ix + (iy-1)*Lx + (iz-1)*Lx*Ly
			Site(Vc)%nb(3) = Site(Site(Vc)%ns(3))%nb(1)
			Site(Vc)%nb(2) = Site(Site(Vc)%ns(2))%nb(4)
			Site(Vc)%nb(5) = Site(Site(Vc)%ns(5))%nb(6)
			write(1,*) Site(Vc)%nb(1:6)
		end do
		end do
		end do
		close(1)
		Ms(1,1) = 1
		Ms(2,1) = 2
		Ms(3,1) = 1
		Ms(4,1) = 2
		Ms(5,1) = 3
		Ms(6,1) = 3
		Ms(1,2) = 1
		Ms(2,2) =-1
		Ms(3,2) =-1
		Ms(4,2) = 1
		Ms(5,2) =-1
		Ms(6,2) = 1
		sCur(1)= 1; sCur(4)= 1; sCur(6)= 1
		sCur(3)=-1; sCur(2)=-1; sCur(5)=-1
	end if

	return
END SUBROUTINE def_latt
!===================================================================

!==============define shifting probability =========================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE def_prob
	implicit none

	integer       :: i

	do i=0, MaxCur-1
		BesselRat(i+1,i) = InIm(i,i+1)
		BesselRat(i,i+1) = 1.d0/BesselRat(i+1,i)
	end do

	! for the calculation of energy
	do i=0, MaxCur
		Ine(i) = InpIn(i, betaJ)
	end do

	return
END SUBROUTINE def_prob
!===================================================================

!==============define gamma function=========================
!! THIS IS PROJECT-DEPENDENT 
FUNCTION InpIn(m,bj)
	implicit none

	integer       :: m
	real(8)       :: bj, bj2
	real(8)       :: InpIn
	real(8)       :: Inm, Inm1, Inp1

	bj2 = bj/2.d0

	if( m==0 ) then
		InpIn =(     bj2      +bj2**3.d0/2.d0+bj2**5.d0/12.d0) &
		&     /(1.d0+bj2**2.d0+bj2**4.d0/4.d0+bj2**6.d0/36.d0)
		return
	end if
	Inm =6*gm(m+1,m+4) + 6*bj2**2.d0*gm(m+2,m+4) &
		& + 3*bj2**4.d0*gm(m+3,m+4) + gm(m+4,m+4)*bj2**6.d0
	Inm1=6*bj2**(-1.d0)*gm(m,m+4) + 6*bj2*gm(m+1,m+4) &
		& + 3*bj2**3.d0*gm(m+2,m+4) + bj2**5.d0*gm(m+3,m+4)
	Inp1=6*bj2*gm(m+2,m+4) + 6*bj2**3.d0*gm(m+3,m+4) &
		& + 3*bj2**5.d0*gm(m+4,m+4) + bj2**7.d0

	InpIn = (Inm1+Inp1)/Inm/2.d0

	return
END FUNCTION InpIn
!===================================================================

! In(x)/Im(x)
FUNCTION InIm(n,m)
	implicit none

	integer       :: n, m
	integer       :: nmax
	real(8)       :: bj2
	real(8)       :: InIm
	real(8)       :: In0, Im0
	integer       :: i

	if( n>=m ) stop "Err: InIm"

	bj2 = betaJ/2.d0

	In0 = 0.d0;   Im0 = 0.d0
	nmax = 10
	do i=0, nmax
		In0 = In0 + bj2**(n*1.d0+2.d0*i)*gm(i+1,nmax)*gm(n+1+i,n+nmax)
		Im0 = Im0 + bj2**(m*1.d0+2.d0*i)*gm(i+1,nmax)*gm(m+1+i,m+nmax)
	end do

	InIm = In0/Im0*gm(n+nmax+1, m+nmax)

	return
END FUNCTION InIm
!===================================================================

FUNCTION gm(n,m)
	implicit none

	integer        :: n, m
	real(8)        :: gm
	integer        :: i

	if(n==0) n=1

	gm=1.d0
	do i=n, m
		gm = gm*i
	end do

	return
END FUNCTION gm

!==========define spin and bond configuration ======================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE def_conf
	implicit none
	integer :: Vc, Ltc

	!-- initialize bond state --------------------------
	Bond  = 0         ! <= nnb2: positive current,  > nnb2: negative current
	Ira   = 0
	Rim   = 0
	Masha = 0
	gr    = 0.d0
	gri   = 0.d0
	ZG    = .False.
	step  = 0.d0

	return
END SUBROUTINE def_conf

