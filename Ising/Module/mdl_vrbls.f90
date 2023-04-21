!-- the model's variables -------------------------------------------------
MODULE mdl_vrbls
	implicit none
	!-- constant -----------------------------------------------------
	real(8), parameter          :: PI  = 3.1415926535898d0
	real(8), parameter          :: DPI = PI*2.d0
	real(8), parameter          :: PI2 = PI*0.5d0

	!-- Lattice and State --------------------------------------------
	!! THIS IS PROJECT-DEPENDENT
	integer(4)                  :: D                  ! dimensionality
	integer(4)                  :: Lx, Ly, Lz
	integer(4)                  :: subl
	integer(4)                  :: nnr = 1            ! nnr types of neighbors
	integer(4), parameter       :: maxir = 4
	integer(4)                  :: Vol
	real(8)                     :: wv                 ! 1/Vol
	real(8)                     :: beta, Jcp(maxir)
	real(8)                     :: betaJ(maxir)         ! = beta*Jcp
	integer(4), parameter       :: MxLx = 1024        ! maximum linear size
	integer(4)                  :: nw                 ! nw sweeps

	integer(4), allocatable     :: spin(:,:)

	logical                     :: ZG
	integer(4)                  :: IraMasha(2)
	integer(4), allocatable     :: current(:,:)
	integer(4)                  :: BondNum(maxir)=0
	real(8)                     :: step
	integer(4)                  :: nZG

	!-- Lattice ---------------------------------------
	integer(4)                  :: nnb(maxir)           ! the number of neighbor
	integer(4)                  :: Nb(maxir)            ! the number of neighbor
	real(8), allocatable        :: Scoordinate(:,:)   ! site real coordinates
	integer(4), allocatable     :: Icoordinate(:,:)   ! site integer coordinates
	integer(4), allocatable     :: sublattice(:)      ! site belong to which sublattice
	integer(4), allocatable     :: nnsite(:,:,:)      ! neighbors' information
	integer(4), allocatable     :: bond(:,:,:)        ! bond
	integer(4), allocatable     :: sitebond(:,:,:)        ! bond
	integer(4), allocatable     :: backdir(:,:,:)     ! the backward direction
	real(8), allocatable        :: sublatvec(:,:)     ! site real coordinates
	real(8), allocatable        :: reclatvec(:,:)     ! site real coordinates

	!-- Probabilities -----------------------------------------------------
	!! THIS IS PROJECT-DEPENDENT
	real(8)                     :: Pbadd(maxir)
	!----------------------------------------------------------------------

	!-- General observables -----------------------------------------------
	real(8)                     :: totalE
	!--- SSF ---!
	real(8), allocatable        :: Sr(:,:)
	complex(8), allocatable     :: Aq(:,:,:,:,:)
	complex(8), allocatable     :: Sq(:,:,:,:,:)
	complex(8), allocatable     :: SqSq(:,:,:,:,:)
	!--- correlation function ---!
	real(8), allocatable        :: SrSr(:)
	!--- histogram ---!
	real(8)                     :: energy_max, energy_min, energy_unit, energy_ave
	integer(4), parameter       :: energy_n=1000
	real(8)                     :: histogram_energy(-energy_n:energy_n)=0.d0

	!-- parallel tempering FLAG -------------------------------------------
	! PT=0: no parallel tempering
	! PT=1: parallel tempering with    adjusting parameter
	! PT=2: parallel tempering without adjusting parameter
	integer(4)                  :: PT

END MODULE mdl_vrbls
