MODULE fftw_vrbls
	use, intrinsic :: iso_c_binding
	implicit none
#include <fftw3.f03>

	complex(C_DOUBLE_COMPLEX), allocatable :: fftw2_in(:,:), fftw2_out(:,:)
	complex(C_DOUBLE_COMPLEX), allocatable :: fftw3_in(:,:,:), fftw3_out(:,:,:)
	type(C_PTR)                            :: plan
END MODULE fftw_vrbls
