

SUBROUTINE thermalization
	implicit none
	integer(4)                  :: isamp
	real(8)                     :: unitbeta
	integer(4)                  :: nbeta


	write(*,*) "Start to calculate thermalization time"

	open(1,file="energy.txt",action="write")
	write(1,'(I9,ES16.8)') 0, totalE*wv

	do isamp = 1, Ntoss
		call move
		if( mod(isamp,100)==0 ) then
			totalE = cal_energy()
			write(1,'(I9,ES16.8)') isamp+Ntoss, totalE*wv
		end if
	end do

	close(1)

	stop
END SUBROUTINE thermalization


SUBROUTINE thermal_equilibrium
	implicit none
	integer(4)                  :: isamp
	real(8)                     :: unitbeta, beta0
	integer(4)                  :: nbeta

	write(*,*) "Start thermal equilibrium by gradually reducing the temperature"

	nbeta = 100
	beta0 = beta
	unitbeta = beta/nbeta
	beta = unitbeta
	call redefine_para(beta)

	do isamp = 1, Ntoss
		if( mod(isamp, Ntoss/nbeta)==0 ) then
			beta = unitbeta * (isamp/(Ntoss/nbeta))
			call redefine_para(beta)
		end if

		call markov(1)
	end do

	beta = beta0
	call redefine_para(beta)

	do isamp = 1, Ntoss
		call markov(1)
	end do

END SUBROUTINE thermal_equilibrium


SUBROUTINE restart
	implicit none
	integer(4)                  :: isamp

	call init_cnf

	call thermal_equilibrium
END SUBROUTINE restart
