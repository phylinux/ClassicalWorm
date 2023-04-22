
INCLUDE "my_vrbls.f90"
INCLUDE "pro_vrbls.f90"
!=====Main routine =================================================
PROGRAM main
	use my_vrbls
	use pro_vrbls
	implicit none
	integer :: itoss,isamp,iblck,pblck
	integer :: coar_num                    ! count the coarsen times
	integer :: flg
	integer :: smkv

	write(*,*) 'flg, pvor, D, Lx, Ly, Lz, MaxCur, Ntoss, Nsamp, betaJ, smkv, nw, NBlck, Seed'
	read(*,*)   flg, pvor, D, Lx, Ly, Lz, MaxCur, Ntoss, Nsamp, betaJ, smkv, nw, NBlck, Seed
	write(*,*)  flg, pvor, D, Lx, Ly, Lz, MaxCur, Ntoss, Nsamp, betaJ, smkv, nw, NBlck, Seed
	Totsamp = Nsamp*NBlck/1000
	call system('rm -f configuration.txt vector.txt locvor.txt')

	if((Ly/2)*2/=Ly) then
		write(6,*) 'L be even?';   stop
	endif

	!--- Initialization ----------------------------------------------
	!call system('rm -f '//trim(datafile))
	call set_time_elapse
	call initialize
	!   call time_elapse
	t_init = t_elap
	write(6,50) t_init
	50 format(/'        set up time:',f16.7,2x,'s')

	!--- Thermialization ---------------------------------------------
	do iblck = 1, NBlck
		do isamp = 1, Ntoss
			call markov(nw)
		enddo
	enddo
	call time_elapse
	t_toss = t_elap
	write(6,51) t_toss
	51 format(/'thermalization time:',f16.7,2x,'s')

	!--- Determine the value of nw -----------------------------------
	pblck = 0
	step  = 0.d0
	do while ( step<smkv )
		call markov(nw)
		pblck = pblck+1
	end do
	nw = max(Vol*nw/int(step/pblck),1)
	print*, "nw=", nw

	!--- Simulation --------------------------------------------------
	call time_elapse
	t_simu   = 0.d0;   t_meas = 0.d0
	pblck = 1 ! if need to continue previous, pblck should be modified
	! if count the coarsen times and control the frequency of output result
	coar_num = 0
	gri = 0.d0
	DO
		DO iblck = pblck, NBlck
			DO isamp = 1, Nsamp
				call markov(nw)    ! call time_elapse;    t_simu = t_simu+t_elap
				call measure
				! call time_elapse;    t_meas = t_meas+t_elap
				call coll_data(iblck);
			ENDDO
			call norm_Nsamp(iblck)
			! THIS IS DEPENDENT ON THE INPUT PARAMETER
			!if( coar_num>3 ) then
			!	if( mod((iblck-NBlck/2)),NBlck/8)==0 ) then
			!		call stat_analy(iblck)
			!		call write2file
			!	end if
			!end if
		ENDDO
		call time_elapse; t_simu = t_simu + t_elap

		!--- Statistics --------------------------------------------------
		call stat_analy(NBlck)
		call write2file
		!call stat2file
		!call conf2file

		!--- Coarsen data ------------------------------------------------
		call coarsen_data
		pblck = NBlck/2+1
		coar_num = coar_num + 1
		call print_conf
		STOP
	END DO

	STOP

CONTAINS
	INCLUDE "initialize.f90"
	INCLUDE "markov.f90"
	INCLUDE "measure.f90"
	INCLUDE "statistics.f90"
	INCLUDE "my_rng.f90"
	INCLUDE "pconf.f90"

END PROGRAM main
!=====================================================================
