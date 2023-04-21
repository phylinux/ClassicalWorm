!===============Shift register random number generator =============
!  very long period sequential version
!! THIS IS PROJECT-INDEPENDENT 
SUBROUTINE set_RNG
	implicit none
	integer :: i_r,k_r,k1_r
	integer :: iseed

	nrannr = mxrn
	iseed  = iabs(Seed)+1
	k_r    = 3**18+2*iseed
	k1_r   = 1313131*iseed
	k1_r   = k1_r-(k1_r/mod2)*mod2

	do i_r = 1, len1
		k_r  = k_r *mult
		k1_r = k1_r*mul2
		k1_r = k1_r-(k1_r/mod2)*mod2
		ir1(i_r) = k_r+k1_r*8193
	enddo

	do i_r = 1, len2
		k_r  = k_r *mult
		k1_r = k1_r*mul2
		k1_r = k1_r-(k1_r/mod2)*mod2
		ir2(i_r) = k_r+k1_r*4099
	enddo

	do i_r = 1, len1
		inxt1(i_r) = i_r+1
	enddo
	inxt1(len1) = 1
	ipnt1 = 1
	ipnf1 = ifd1+1

	do i_r = 1, len2
		inxt2(i_r) = i_r+1
	enddo
	inxt2(len2) = 1
	ipnt2 = 1
	ipnf2 = ifd2 + 1
	return
END SUBROUTINE set_RNG 
!===================================================================


!===============Calculate next random number =======================
!! THIS IS ALMOST PROJECT-INDEPENDENT 
double precision function rn()
! integer function rn()
	implicit none
	integer   :: i_r, l_r, k_r
	nrannr = nrannr +1
	if(nrannr>=mxrn) then
		nrannr = 1
		do i_r= 1, mxrn
			l_r = ieor(ir1(ipnt1),ir1(ipnf1))
			k_r = ieor(ir2(ipnt2),ir2(ipnf2))
			irn(i_r) = ieor(k_r,l_r)
			ir1(ipnt1)=l_r
			ipnt1 = inxt1(ipnt1)
			ipnf1 = inxt1(ipnf1)
			ir2(ipnt2) = k_r
			ipnt2 = inxt2(ipnt2)
			ipnf2 = inxt2(ipnf2)
		enddo
	endif 
	! rn = irn(nrannr)
	rn = irn(nrannr)*tm32+0.5d0
end function rn
!===================================================================

!==============Trace elapsed time ==================================
!! THIS IS PROJECT-INDEPENDENT 
SUBROUTINE set_time_elapse
	implicit none
	!-- read and calculate time (in seconds) -------------------------
	call date_and_time(date, time, zone, tval)
	t_curr = tval(5)*3600.d0+tval(6)*60.d0+tval(7)+tval(8)*0.001d0 
	h_curr = tval(5)
	t_prev = t_curr
	h_prev = h_curr
	return
END SUBROUTINE set_time_elapse


!==============Trace elapsed time ==================================
!! THIS IS PROJECT-INDEPENDENT 
SUBROUTINE time_elapse
	implicit none

	!-- read and calculate time (in seconds) -------------------------
	call date_and_time(date, time, zone, tval)
	t_curr = tval(5)*3600.d0+tval(6)*60.d0+tval(7)+tval(8)*0.001d0 
	h_curr = tval(5)

	t_elap = t_curr-t_prev
	if(h_curr<h_prev) t_elap = t_elap+24*3600.d0
	t_prev = t_curr
	h_prev = h_curr 
	return
END SUBROUTINE time_elapse
!===================================================================

