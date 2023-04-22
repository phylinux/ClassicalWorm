!  Model : Classical XY model using worm algorithm
!  Author: Chun-Jiong Huang
!  Build : Apr 21th, 2023.
!===================================================================

!=====Main routine =================================================
SUBROUTINE XY
	implicit none
	integer(4)                  :: itoss
	integer(4)                  :: iblck,iblck0,isamp,isamp0
	integer(4)                  :: id
	integer(4)                  :: ir

	id=1

	open(id,file="inp00",action="read")
	read(id,*) r_cnf_stat
	read(id,*) PT
	read(id,*) D
	read(id,*) Lx
	read(id,*) Ly
	read(id,*) Lz
	read(id,*) subl
	read(id,*) nnr
	read(id,*) Jcp(1)
	read(id,*) beta
	read(id,*) Maxcur
	read(id,*) nw
	read(id,*) NBlck
	read(id,*) Ntoss
	read(id,*) Nsamp
	read(id,*) Seed
	close(id)

	if( nnr/=1 ) then
		write(*,*) "This version can only simulate the nearest neighbor, nnr=1"
		stop
	end if

	Totsamp = Nsamp*NBlck/1000

	!--- Initialization ----------------------------------------------
	call set_time_elapse
	call initialize
	call time_elapse
	t_init = t_elap
	write(*,'(A23,F10.2,A4)') "Initialization time = ", t_init/60.d0, "(m)"

#ifdef ANNEALING
	call annealing
	stop
#endif

	if( r_cnf_stat==0 ) then
		!--- Thermialization ----------------------------------------
		!call system('rm -f '//trim(datafile))
		write(*,'(A36)') "=== Start new simulation ==="

#ifdef THERMALIZATION
	call thermalization
#endif

		call thermal_equilibrium

		call write_cnf_stat
		call time_elapse
		t_ther = t_elap
		write(*,'(A23,F10.2,A4)') "Thermalization time = ", t_ther/60.d0, "(m)"
	else if ( r_cnf_stat==1 ) then
		!-- read cnf ------------------------------------------------------
		call read_cnf_stat(1,0)
		call prepare_simulation
		!call system('rm -f '//trim(datafile))
		write(*,'(A36)') "=== Start on old configuration ==="
		do isamp = 1, Ntoss
			call markov(1)
		end do
		call write_cnf_stat
		call time_elapse
		t_ther = t_elap
		write(*,'(A23,F10.2,A4)') "Thermalization time = ", t_ther/60.d0, "(m)"
	else if ( r_cnf_stat==2 ) then
		!-- read cnf, stat, rng -------------------------------------------
		call read_cnf_stat(1,1)
		call prepare_simulation
		!call system('rm -f '//trim(datafile))
		write(*,'(A36)') "=== Continue old simulation ==="
	else
		write(*,'(A23)') "Err: r_cnf_stat"
		stop
	end if

	if( r_cnf_stat==0 ) then
		iblck = 0
		step  = 0.d0
		do while ( step<1.d4*Vol )
			call markov(1)
			iblck = iblck+1
		end do
		nw = max(Vol*nw/(step*1.d0/iblck),1.1d0)
		print*, "nw=", nw
		call time_elapse; t_adjt = t_elap
		write(*,'(A23,F10.2,A4)') "Adjust         time = ", t_adjt/60.d0, "(m)"
		call write_cnf_stat

#ifdef MPI
		!-- adjust the temperature in parallel tempering -----------------
		if( PT==1 ) call adjust_parallel_tempering
		call time_elapse; t_prtp = t_elap
		write(*,'(A23,F10.2,A4)') "Parallel tempering  = ", t_prtp/60.d0, "(m)"
		exprob = 0.d0
		call write_cnf_stat
#endif
	end if

#ifdef ACF
	call sample_sequence
#endif

#ifdef HISTOGRAM
	call adjust_histogram
#endif

	!--- Simulation --------------------------------------------------
	write(*,'(A16)') "Start simulation"
	t_simu = 0.d0; t_meas = 0.d0; t_prtp = 0.d0
	accept = 0.d0
	if( r_cnf_stat/=2 ) then
		iblck0 = 1
		isamp0 = 1
	else
		iblck0 = ptb
		isamp0 = pts
		!write(*,*) iblck0, isamp0
	end if
	! if count the coarsen times and control the frequency of output result
	DO iblck = iblck0, NBlck
		DO isamp = isamp0, Nsamp
#ifdef MPI
			if( PT/=0 .and. mod(isamp,npt)==0 ) then
				call parallel_tempering
				call time_elapse; t_prtp = t_prtp+t_elap
			end if
#endif
			call markov(nw)
			call time_elapse; t_simu = t_simu+t_elap
			call measure
			call time_elapse; t_meas = t_meas+t_elap
			call coll_data(iblck);
			pts = isamp
		END DO
#ifdef MPI
		if( PT/=0  ) call exchange_probability
#endif
		call norm_Nsamp(iblck)
		! THIS IS DEPENDENT ON THE INPUT PARAMETER
		isamp0 = 1
		ptb = iblck
		if( mod(ptb,32)==0 ) then
			write(*,'(A8,I8)') "iblck = ", iblck
			call stat_analy(iblck)
			call write2file(iblck)
			call write_cnf_stat
		end if
#ifdef RESTART
		if( mod(iblck,32)==0 ) then
			call restart
			call time_elapse;  t_ther = t_ther+t_elap
		end if
#endif
	END DO

	!--- Statistics --------------------------------------------------
	call stat_analy(NBlck)
	call write2file(NBlck)
	call write_cnf_stat

	call my_fftw_destroy

END SUBROUTINE XY
!=====================================================================
