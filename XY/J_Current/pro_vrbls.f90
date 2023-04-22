MODULE pro_vrbls
	IMPLICIT NONE

	!-- statistic parameters ----------------------------------------
	! THIS IS ALMOST PROJECT-INDEPENDENT 
	real(8),          parameter :: PI=3.141592653589793d0
	real(8),          parameter :: PI2=2.d0*PI
	integer,          parameter :: MxBlck = 2**10          ! maximum number of blocks for statistics
	integer,          parameter :: MnBlck = 2**6           ! minimum number of blocks

	integer            :: NBlck                            ! # blocks
	integer            :: Nsamp                            ! # samples in unit 'NBlck'
	integer            :: Totsamp
	integer            :: Ntoss                            ! # samples to be thrown away
	!-----------------------------------------------------------------

	!-- parameters and variables -------------------------------------
	!! THIS IS PROJECT-DEPENDENT 
	character(8)  :: ident     = 'XYmodel'              ! identifier
	character(16) :: paralist  = 'parameter.list'       ! identifier
	character(12) :: datafile  = 'XYmodel.txt'          ! datafile
	character(11) :: datafile1 = 'dat.txt_cor'          ! datafile for big correlation

	integer               :: D                          ! dimensionality

	integer               :: Lx, Ly, Lz, Vol, bVol
	integer               :: gVol                       ! the volumn of green function
	integer               :: Lxyz(3)
	real(8)               :: wv,we                      ! 1/Vol, 1/E
	real(8)               :: betaJ                      ! = beta*Jcp
	integer               :: Ira, Masha
	logical               :: ZG                         ! F: partition function, T: Green function
	real(8)               :: step
	!-----------------------------------------------------------------

	!-- Lattice and State --------------------------------------------
	!! THIS IS PROJECT-DEPENDENT 
	integer, parameter :: MxLy =1024                    ! maximum linear size

	type siteinfo
		integer(4)        :: ns(6)                      ! note the neighbor site
		integer(4)        :: nb(6)                      ! the numbers of linked bonds
	end type

	integer(1)                   :: nnb, nnb2           ! # the number of neighbor, nnb2=nnb/2
	type (siteinfo), allocatable :: Site(:)
	integer(4),      allocatable :: Bond(:)             ! the the currents
	integer(1),      allocatable :: sCur(:)             ! sign of current of different direction
	real(8),         allocatable :: gr(:,:)
	real(8),         allocatable :: gri(:)
	integer(4)                   :: Rim(3)
	integer(4),      allocatable :: Ms(:,:)

	integer                      :: nw                     ! nw sweeps
	integer                      :: MaxCur
	real(8),         allocatable :: Ine(:)
	real(8),         allocatable :: BesselRat(:,:)     ! probability for update
	!-----------------------------------------------------------------

	!-- configuration flag ----------
	integer                      :: pvor

	!-- Observables --------------------------------------------------
	!! THIS IS PROJECT-DEPENDENT 

	! GO TO MODIFY PARAMETER LIST in initialize.f90
	integer, parameter  :: NObs_b =  8                        ! #basic     observable
	integer, parameter  :: NObs_c =  1                        ! #composite observables
	integer, parameter  :: NObs   = NObs_b+NObs_c             ! #total     observables
	!-----------------------------------------------------------------

	!-- Statistics ---------------------------------------------------
	! THIS IS PROJECT-DEPENDENT 
	real(8)          :: Quan(NObs_b)            ! Measured quantities
	real(8)          :: Obs(NObs,MxBlck)             ! 1st--#quan.  2nd--#block
	real(8)          :: Ave(NObs), Dev(NObs), Cor(NObs)   ! average, error bars, and correlation of observables
	character(16)    :: ParaName(NObs)
	!-----------------------------------------------------------------

END MODULE pro_vrbls
