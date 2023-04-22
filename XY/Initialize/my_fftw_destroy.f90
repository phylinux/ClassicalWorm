SUBROUTINE my_fftw_destroy
	implicit none

	call fftw_destroy_plan(plan)
	select case (D)
	case(2)
		deallocate(fftw2_in)
		deallocate(fftw2_out)
	case(3)
		deallocate(fftw3_in)
		deallocate(fftw3_out)
	end select
END SUBROUTINE my_fftw_destroy
