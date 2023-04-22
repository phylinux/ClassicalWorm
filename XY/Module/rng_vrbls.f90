!============== variable module ====================================
MODULE rng_vrbls
	IMPLICIT NONE

	!-- common parameters and variables ------------------------------
	! THIS IS PROJECT-INDEPENDENT 
	real(8), parameter          :: tm32   = 1.d0/(2.d0**32.d0)
	integer(4), parameter       :: Mxint  = 2147483647! maximum integer
	integer(4), parameter       :: Mnint  =-2147483647! minimum integer

	!-- Random-number generator---------------------------------------
	! THIS IS PROJECT-INDEPENDENT 
	integer, parameter          :: mult=32781
	integer, parameter          :: mod2=2796203, mul2=125
	integer, parameter          :: len1=9689,    ifd1=471
	integer, parameter          :: len2=127,     ifd2=30
	integer, dimension(1:len1)  :: inxt1
	integer, dimension(1:len2)  :: inxt2
	integer, dimension(1:len1)  :: ir1
	integer, dimension(1:len2)  :: ir2
	integer                     :: ipnt1, ipnf1
	integer                     :: ipnt2, ipnf2
	integer, parameter          :: mxrn = 10000
	integer, dimension(1:mxrn)  :: irn(mxrn)

	integer                     :: Seed               ! random-number seed
	integer                     :: nrannr             ! random-number counter
	!-----------------------------------------------------------------

	!-- time-checking variables --------------------------------------
	! THIS IS PROJECT-INDEPENDENT 
	character( 8)               :: date
	character(10)               :: time
	character(5 )               :: zone
	integer, dimension(8)       :: tval
	real(8)                     :: t_prev, t_curr, t_elap
	integer                     :: h_prev, h_curr
	real(8)                     :: t_init, t_ther, t_adjt, t_prtp, t_simu, t_meas
	!-----------------------------------------------------------------
END MODULE rng_vrbls
!===================================================================
