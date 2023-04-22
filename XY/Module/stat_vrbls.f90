!-- the statistic variables -------------------------------------------------
MODULE stat_vrbls
	IMPLICIT NONE
	!-- constant ----------------------------------------------------------
	real(8), parameter          :: eps    = 1.d-14    ! very small number
	real(8), parameter          :: tol    = 0.2d0     ! tolerance for Cor

	!-- statistic parameters ----------------------------------------
	! THIS IS ALMOST PROJECT-INDEPENDENT 
	integer(4), parameter       :: MxBlck = 2**10     ! maximum number of blocks for statistics
	integer(4), parameter       :: MnBlck = 2**6      ! minimum number of blocks

	integer(4)                  :: NBlck              ! # blocks
	integer(4)                  :: ptb=1              ! # block point
	integer(4)                  :: Nsamp              ! # samples in unit 'NBlck'
	integer(4)                  :: pts=1              ! # sample point
	integer(4)                  :: Totsamp
	integer(4)                  :: Ntoss              ! # samples to be thrown away
	logical                     :: prt
	integer(4)                  :: r_cnf_stat
	!-----------------------------------------------------------------

	!-- parameters and variables -------------------------------------
	!! THIS IS PROJECT-DEPENDENT 
	character(16), parameter    :: ident     = 'XY'          ! identifier
	character(16)               :: paralist  = 'parameter.list' ! identifier
	character(16)               :: datafile  = trim(ident)//'.txt'      ! datafile
	character(16)               :: datafile1 = 'dat.txt_cor'    ! datafile for big correlation
	!-----------------------------------------------------------------

	!-- Observables --------------------------------------------------
	!! THIS IS PROJECT-DEPENDENT 

	! GO TO MODIFY PARAMETER LIST in initialize.f90
	integer(4), parameter   :: NObs_b =  6            ! #basic     observables
	integer(4), parameter   :: NObs_c =  3            ! #composite observables
	integer(4), parameter   :: NObs   = NObs_b+NObs_c ! #total  observables
	!-----------------------------------------------------------------

	!-- Statistics ---------------------------------------------------
	! THIS IS PROJECT-DEPENDENT 
	real(8)                 :: Quan(NObs_b)         ! Measured quantities
	real(8)                 :: Obs(NObs,MxBlck)     ! 1st--#quan.  2nd--#block
	real(8)                 :: Ave(NObs), Dev(NObs), Cor(NObs)   ! average, error bars, and correlation of observables
	character(16)           :: ParaName(NObs)
	!-----------------------------------------------------------------

	!-- the accept ratio --------------------------------------------------
	real(8)                     :: accept(0:10,2)
	!----------------------------------------------------------------------

	!-- Idles -------------------------------------------------------------
	!integer(4)                  :: id               ! the idle of file
	!----------------------------------------------------------------------

END MODULE stat_vrbls
