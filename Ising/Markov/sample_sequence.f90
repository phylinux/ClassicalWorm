SUBROUTINE sample_sequence
	implicit none

	integer(4)                  :: tolsample
	integer(4)                  :: inw
	integer(4)                  :: i

	write(*,*) "Start to collect samples for autocorrelation function"

	tolsample = 10**7
	inw = 1
	open(1,file="sample_sequence.txt",action="write")
	write(1,'(I9, I9)') inw*Vol, tolsample
	write(1,'(I9    )') 0
	do i = 1, tolsample
		call markov(inw)
		totalE = cal_energy()
		write(1,'(ES16.8)') totalE*wv
		if( mod(i,10000)==0 )  call normalize_spin
	end do
	close(1)
	write(*,*) "Samples for autocorrelation function is done"

#ifdef MPI
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	call MPI_FINALIZE(ierr)
#endif

	stop
END SUBROUTINE sample_sequence
