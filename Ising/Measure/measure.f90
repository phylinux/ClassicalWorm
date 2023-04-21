!==============Measurement =========================================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE measure
	implicit none

	real(8)                     :: energy
	integer(4)                  :: i, j

	totalE = cal_energy()
	energy = totalE*wv

	!call structure_factor

#ifdef HISTOGRAM
	call histogram
#endif

	Quan = 0.d0
	Quan( 1) = energy
	Quan( 2) = energy**2
	Quan( 3) = (step*1.d0/nw)**1.d0
	Quan( 4) = (step*1.d0/nw)**2.d0
	Quan( 5) = BondNum(1)**1.d0
	Quan( 6) = BondNum(1)**2.d0

END SUBROUTINE measure

!============== calculate correlation of spin ========================
!! THIS IS PROJECT-DEPENDENT 
!=====================================================================

!==============Calculate Binder ratio 1================================
!! Q=Ave(b2)/Ave(b1)^epo
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE cal_Q(ib, jb,b1,b2,epo)
	implicit none
	integer, intent(in) :: ib, jb, b1,b2
	real(8), intent(in) :: epo
	integer             :: k
	real(8)             :: tmp

	!-- Obs(j,k) series --------------------------------------------
	do k = 1, ib
		tmp = Obs(b1,k)**epo
		if(dabs(Obs(b2,k))>eps) tmp = tmp/Obs(b2,k)
		Obs(jb,k) = tmp
	enddo
	Ave(jb) = SUM(Obs(jb,1:ib))/ib

END SUBROUTINE cal_Q

!==============Calculate fluctuation quantity ==================
!! THIS IS PROJECT-INDEPENDENT 
!! C = V(Ave(b2)-Ave(b1)^2)
SUBROUTINE cal_fluct(ib, jb,b1,b2, sl)
	implicit none
	integer(4), intent(in)      :: ib
	integer(4), intent(in)      :: jb, b1, b2
	real(8), intent(in)         :: sl
	integer(4)                  :: k

	!-- Average ----------------------------------------------------
	!Ave(jb) = Vol*(Ave(b2)-Ave(b1)**2.d0)

	!-- Obs(j,k) series --------------------------------------------
	do k = 1, ib
		Obs(jb,k) = Obs(b2,k)-Obs(b1,k)**2.d0
		Obs(jb,k) = Obs(jb,k)*sl
	end do
	Ave(jb) = SUM(Obs(jb,1:ib))/ib

END SUBROUTINE cal_fluct


SUBROUTINE cal_binder(ib, jb,b1,b2)
	implicit none
	integer(4), intent(in)      :: ib
	integer(4), intent(in)      :: jb, b1, b2
	integer(4)                  :: k

	!-- Average ----------------------------------------------------
	!Ave(jb) = Vol*(Ave(b2)-Ave(b1)**2.d0)

	!-- Obs(j,k) series --------------------------------------------
	do k = 1, ib
		Obs(jb,k) = Obs(b1,k)**2.d0
		if( Obs(jb,k)>eps ) then
			Obs(jb,k) = Obs(b2,k)/Obs(jb,k)
		else
			Obs(jb,k) = Obs(b2,k)
		end if
	end do
	Ave(jb) = SUM(Obs(jb,1:ib))/ib

END SUBROUTINE cal_binder


SUBROUTINE cal_ratio(ib, jb,b1,b2)
	implicit none
	integer(4), intent(in)      :: ib
	integer(4), intent(in)      :: jb, b1, b2
	integer(4)                  :: k

	!-- Average ----------------------------------------------------
	!Ave(jb) = Vol*(Ave(b2)-Ave(b1)**2.d0)

	!-- Obs(j,k) series --------------------------------------------
	do k = 1, ib
		Obs(jb,k) = Obs(b2,k)/Obs(b1,k)
	end do
	Ave(jb) = SUM(Obs(jb,1:ib))/ib

END SUBROUTINE cal_ratio

SUBROUTINE cal_Cv(ib, jb,b1,b2, sl)
	implicit none
	integer             :: ib
	integer, intent(in) :: jb, b1, b2
	real(8), intent(in) :: sl
	integer             :: k
	integer             :: ir

	ir=1
	do k=1, ib
		Obs(jb,k) =  nnb(ir)/2.d0*Vol + &
			&       -(1.d0+1.d0/(tanh(betaJ(ir))**2)) * Obs(b1,k) &
			&       +(Obs(b2,k)-Obs(b1,k)**2)/(sinh(betaJ(ir))**2)
		Obs(jb,k) = Obs(jb,k)*sl
	end do
	Obs(jb,1:ib) = Obs(jb,1:ib)*(betaJ(ir)**2)/(cosh(betaJ(ir))**2)

	Ave(jb) = SUM(Obs(jb,1:ib))/ib
END SUBROUTINE cal_Cv


!==============Calculate composite observables =====================
!! THIS IS PROJECT-INDEPENDENT 
!! call in 'stat_alan'
SUBROUTINE cal_Obs_comp(ib)
	implicit none
	integer(4), intent(in)      :: ib
	integer(4)                  :: jb, b2, b3, b1

	!-- calculate the average ----------------------------------------

	jb = NObs_b+1;   call cal_Cv(ib, jb, 5, 6, 1.d0/Vol*beta*beta)  ! Cv
	jb = NObs_b+2;   call cal_fluct(ib, jb, 3, 4, 1.d0/Vol)  ! X0
	jb = NObs_b+3;   call cal_fluct(ib, jb, 5, 6, 1.d0/Vol)  ! X0

END SUBROUTINE cal_Obs_comp
!===================================================================
