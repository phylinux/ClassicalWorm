
SUBROUTINE print_conf()
	implicit none

	integer      :: ib

	open(1, file='current.txt', action='write')
	do ib=1, bVol
		write(1,'(2I10)') ib, Bond(ib)
	end do
	close(1)

	return
END SUBROUTINE print_conf
