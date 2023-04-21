
!! THIS IS PROJECT-DEPENDENT
SUBROUTINE markov(inw)
	implicit none
	integer(4), intent(in)      :: inw
	integer(4)                  :: i

	do i= 1, inw
		call worm
	end do

END SUBROUTINE markov
