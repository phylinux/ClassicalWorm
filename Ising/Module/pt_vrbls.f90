MODULE pt_vrbls
	implicit none
	!-- parallel tempering ------------------------------------------------
	real(8), parameter          :: eps_con=0.02
	integer(4)                  :: npt=20
	real(8)                     :: exprob(0:1000) = 0.d0
	real(8)                     :: betalist(0:1000)
	real(8)                     :: converge
END MODULE pt_vrbls
