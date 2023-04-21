!--- PROJECT-DEPENDENT ---------------------------------------------
!==============Initialization ======================================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE initialize
	implicit none

	!-- order should be changed --------------------------------------
	call tst_and_prt
	call set_RNG
	call def_para
	!call def_latt
	call read_latt
	call alloc_arr
	call def_prob
	call init_cnf
	call init_fftw
	call init_vrbls

END SUBROUTINE initialize

!==============Test and Print ======================================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE tst_and_prt
	implicit none
	integer(4)                  :: ir
	character(50)               :: charJ

	call redefine_para(beta)

	!-- Test and Print -----------------------------------------------
	!if((NBlck>MxBlck).or.(NBlck<MnBlck)) then
	!	write(6,*) 'MnBlck <= NBlck <=MxBlck?'
	!	stop
	!endif

	if(Lx>MxLx) then
		write(6,*) 'Lx<=MxLx?'
		stop
	endif

	write(*,'(A13,I2)') "r_cnf_stat = ", r_cnf_stat
	write(*,'(A12,A32)') "Algorithm = ", "worm"

	write(*,'(A18)')  "Spin       = Ising"
	write(*,'(A6,I6,A3,I6,A3,I6,A16)') "AT on ",Lx," x ",Ly, " x ", Lz, " periodic box"
	write(*,'(A8,F10.6)') "beta  = ", beta
	write(*,'(A8,F10.6)') "Jcp   = ", Jcp(1)
	write(*,'(A8,I10)')   "nw    = ", nw
	write(*,'(A8,I10)')   "NBlck = ", NBlck
	write(*,'(A8,I10)')   "Ntoss = ", Ntoss
	write(*,'(A8,I10)')   "Nsamp = ", Nsamp

END SUBROUTINE tst_and_prt

!============== definite parameter ====================================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE def_para
	IMPLICIT NONE
	integer    :: i

	ParaName     = 'XXX'
	!------------- 1234567890123456 -------
	ParaName( 1) = 'magnatization^2'

	open(1,file=trim(paralist),action='write')
	do i=1, NObs
		write(1,'(1X,I2,4X,A)') i, ParaName(i)
	end do
	close(1)
END SUBROUTINE def_para
!===================================================================

SUBROUTINE alloc_arr
	implicit none

	allocate(spin(nnr,Vol))
	allocate(current(nnr,maxval(Nb)))

	allocate(Sr(Vol,3))
	allocate(Sq(subl,Lx,Ly,Lz,3))
	allocate(Aq(subl,Lx,Ly,Lz,3))
	allocate(SqSq(subl,subl,Lx,Ly,Lz))

	allocate(SrSr(Lx))

	Sr = 0.d0
	Sq = CMPLX(0.d0,0.d0)
	SqSq = CMPLX(0.d0,0.d0)
	SrSr = 0.d0

END SUBROUTINE alloc_arr
!===================================================================


!==============define shifting probability =========================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE def_prob
	implicit none
	integer(4)                  :: i

	Pbadd = tanh(betaJ)

END SUBROUTINE def_prob
!===================================================================

!==========define spin and bond configuration ======================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE init_cnf
	implicit none
	integer(4)                  :: ir


	!--- initialize conf ---!
	spin = 1
	current = 0
	IraMasha  = 0
	ZG = .false.

	do ir=1, nnr
		BondNum(ir) = sum(current(ir,:))
	end do

	step = 0
	nZG = 1

END SUBROUTINE init_cnf

SUBROUTINE prepare_simulation
	implicit none
#ifdef MPI
	integer                     :: stat(MPI_STATUS_SIZE)

	itag = 100
	if( taskid==0 ) then
		do dest=1, numprocs-1
			beta = betalist(dest)
			call MPI_SEND(beta, 1, MPI_REAL8, dest, itag, MPI_COMM_WORLD, ierr)
		end do
		beta = betalist(0)
	else
		source = 0
		call MPI_RECV(beta, 1, MPI_REAL8, source, itag, MPI_COMM_WORLD, stat, ierr)
	end if
#endif

	call redefine_para(beta)
END SUBROUTINE prepare_simulation

SUBROUTINE init_fftw
	implicit none

	select case (D)
	case(2)
		call init_fftw2
	case(3)
		call init_fftw3
	case default
		write(*,*) "Err: D/=2,3"
		stop
	end select
END SUBROUTINE init_fftw

SUBROUTINE init_fftw2
	implicit none
	allocate(fftw2_in(Lx,Ly))
	allocate(fftw2_out(Lx,Ly))

	plan = fftw_plan_dft_2d(Ly,Lx, fftw2_in, fftw2_out, FFTW_FORWARD, FFTW_ESTIMATE)
END SUBROUTINE init_fftw2

SUBROUTINE init_fftw3
	implicit none
	allocate(fftw3_in(Lx,Ly,Lz))
	allocate(fftw3_out(Lx,Ly,Lz))

	plan = fftw_plan_dft_3d(Lz,Ly,Lx, fftw3_in, fftw3_out, FFTW_FORWARD, FFTW_ESTIMATE)
END SUBROUTINE init_fftw3

SUBROUTINE init_vrbls
	implicit none
	!-- measurement initialization -----------------------------------
	Obs = 0.d0;   Quan = 0.d0
	Ave = 0.d0;   Dev  = 0.d0;     Cor = 0.d0
	!-- the accepted ratio ------------------------------------------------
	accept = 0.d0
END SUBROUTINE init_vrbls

SUBROUTINE redefine_para(newbeta)
	implicit none
	real(8), intent(in)         :: newbeta

	betaJ  = newbeta*Jcp

	call def_prob
END SUBROUTINE redefine_para
