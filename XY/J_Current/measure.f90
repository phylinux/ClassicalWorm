
!==============Measurement =========================================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE measure
	implicit none

	real(8)          :: energy
	integer          :: i, j, k
	real(8)          :: wx, wy, wz

	energy = cal_energy()
	!do i=1, Vol
	!	k = 0
	!	do j=1, nnb
	!		k = k + Bond(Site(i)%nb(j))*sCur(j)
	!	end do
	!	if( k/=0 ) stop "Err: current is not 0"
	!end do

	wx = winding(1)
	wy = winding(2)
	wz = winding(3)

	Quan( 1) = -energy*wv
	Quan( 2) = energy*energy
	Quan( 3) = wx
	Quan( 4) = wy
	Quan( 5) = wz
	Quan( 6) = (wx+wy+wz)/D
	Quan( 7) = step/nw
	Quan( 8) = sum(abs(Bond))*1.d0*wv

	step = 0.d0

	return
END SUBROUTINE measure

!============== calculate energy ========================
!! THIS IS PROJECT-DEPENDENT 
Function cal_energy()
	implicit none

	real(8)       :: cal_energy
	integer       :: bn
	integer       :: i

	cal_energy = 0.d0
	do i=1, bVol
		bn = abs(Bond(i))
		cal_energy = cal_energy+ Ine(bn)
	end do

	return
END Function cal_energy
!=====================================================================

!============== calculate winding number ========================
!! THIS IS PROJECT-DEPENDENT 
Function winding(xyz)
	implicit none

	integer       :: xyz
	integer       :: winding
	integer       :: i, j

	winding = 0

	if( D==1 ) then
	else if ( D==2 ) then
		if( xyz==1 ) then
			do i=1, Vol, Lx
				winding = winding+Bond(i)
			end do
		else
			do i=1, Lx
				winding = winding+Bond(i+Vol)
			end do
		end if
	else
		if( xyz==1 ) then
			do j=1, Lz
			do i=1, Lx*Ly, Lx
				winding = winding+Bond(i+(j-1)*Lx*Ly)
			end do
			end do
		elseif( xyz==2 ) then
			do j=1, Lz
			do i=1, Lx
				winding = winding+Bond(i+Vol+(j-1)*Lx*Ly)
			end do
			end do
		else
			do i=1, Lx*Ly
				winding = winding+Bond(i+Vol*2)
			end do
		end if
	end if

	winding = winding*winding

	return
END Function winding
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
		tmp = Obs(b1,k)**epo;   if(dabs(tmp)>eps) tmp = Obs(b2,k)/tmp
		Obs(jb,k) = tmp
	enddo
	Ave(jb) = SUM(Obs(jb,1:ib))/ib

	return
END SUBROUTINE cal_Q

!==============Calculate fluctuation quantity ==================
!! THIS IS PROJECT-INDEPENDENT 
!! C = V(Ave(b2)-Ave(b1)^2)
SUBROUTINE cal_fluct(ib, jb,b1,b2)
	implicit none
	integer             :: ib
	integer, intent(in) :: jb, b1, b2
	integer             :: k

	!-- Average ----------------------------------------------------
	!Ave(jb) = Vol*(Ave(b2)-Ave(b1)**2.d0)

	!-- Obs(j,k) series --------------------------------------------
	do k = 1, ib
		Obs(jb,k) = Obs(b2,k)-Obs(b1,k)**2.d0
	end do
	Ave(jb) = SUM(Obs(jb,1:ib))/ib

	return
END SUBROUTINE cal_fluct


!==============Calculate composite observables =====================
!! THIS IS PROJECT-INDEPENDENT 
!! call in 'stat_alan'
SUBROUTINE cal_Obs_comp(ib)
	implicit none
	integer       :: ib
	integer       :: jb, b2, b3, b1

	!-- calculate the average ----------------------------------------

	!jb = NObs_b+1;   call cal_Q(ib, jb, 1, 2, 2.d0)         ! Q1 = <m^4>/<m^2>^2
	jb = NObs_b+1;   call cal_fluct(ib, jb, 1, 2)           ! Magnetic susceptibility

	return
END SUBROUTINE cal_Obs_comp
!===================================================================
