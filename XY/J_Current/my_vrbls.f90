
!*******************************************************************
! Ising model on the square Lattice

! Error bars are calculated using the blocking technique. 
! Let ' \chi^2 = \sum_{t=1}^{T} <(O_t -<O>)^2> ' be the squared chi for
! 'T' blocks of observable 'O'. Assuming each block of data is independent
! of each other, the error bar of 'O' is given by 'Dev = \sqrt{\chi^2/T/(T-1)}.

! Reliabity of the obtained errors is monitored by t=1 correlation,
! for which tolerance is set by variable 'tol' (default: tol=0.20d0).

! Composite quantities like Binder ratios are calculated in each block, and
! the associated error bars are obtained from their fluctuations.

! Results are written into a special file 'dat.***' if the number of
! blocks is less than 125 or correlation is too big. Data in each 
! block will be also printed out in this case.

! Default number of extensive simulation is 'NBlck=1024'.

! For test purpose, for which huge amount of information will be 
! printed out, 'NBlck' should be set smaller but >2.

! Dynamical behavior is not studied.

!  Author: Chunjiong Huang
!  Date  : Oct 12th, 2016.
!*******************************************************************

! Look for 'PROJECT-DEPENDENT' for different projects
!============== variable module ====================================
MODULE my_vrbls
	IMPLICIT NONE

	!-- common parameters and variables ------------------------------
	! THIS IS PROJECT-INDEPENDENT 
	double precision, parameter :: tm32   = 1.d0/(2.d0**32.d0)
	double precision, parameter :: eps    = 1.d-14         ! very small number
	double precision, parameter :: tol    = 0.20d0         ! tolerance for Cor
	logical                     :: prt                     ! flag for write2file
	integer,          parameter :: Mxint  = 2147483647     ! maximum integer
	integer,          parameter :: Mnint  =-2147483647     ! minimum integer

	!-- Random-number generator---------------------------------------
	! THIS IS PROJECT-INDEPENDENT 
	integer, parameter           :: mult=32781
	integer, parameter           :: mod2=2796203, mul2=125
	integer, parameter           :: len1=9689,    ifd1=471
	integer, parameter           :: len2=127,     ifd2=30
	integer, dimension(1:len1)   :: inxt1
	integer, dimension(1:len2)   :: inxt2
	integer, dimension(1:len1)   :: ir1
	integer, dimension(1:len2)   :: ir2
	integer                      :: ipnt1, ipnf1
	integer                      :: ipnt2, ipnf2
	integer, parameter           :: mxrn = 10000
	integer, dimension(1:mxrn)   :: irn(mxrn)

	integer                      :: Seed                   ! random-number seed
	integer                      :: nrannr                 ! random-number counter
	!-----------------------------------------------------------------

	!-- time-checking variables --------------------------------------
	! THIS IS PROJECT-INDEPENDENT 
	character( 8)         :: date
	character(10)         :: time
	character(5 )         :: zone
	integer, dimension(8) :: tval
	double precision      :: t_prev, t_curr, t_elap
	integer               :: h_prev, h_curr
	double precision      :: t_init, t_simu, t_meas, t_toss
	!-----------------------------------------------------------------
END MODULE my_vrbls
!===================================================================
