FUNCTION cal_energy()
	implicit none
	real(8)                     :: cal_energy
	integer(4)                  :: ir

	ir=1
	cal_energy = -Jcp(ir)*Pbadd(ir)*(nnb(ir)*Vol/2+BondNum(ir)/(sinh(betaJ(ir))**2.d0))
END FUNCTION cal_energy

!SUBROUTINE structure_factor
!	implicit none
!	integer(4)                  :: sb1, sb2
!	integer(4)                  :: ix, iy, iz
!	real(8), allocatable        :: iSi(:,:,:)
!	complex                     :: phase
!	real(8)                     :: sbvq
!	integer(4)                  :: Vc
!	integer(4)                  :: i, j
!
!	Sr = Sr + spin
!
!	allocate(iSi(Lx,Ly,Lz))
!	Aq = CMPLX(0.d0,0.d0)
!
!	do i=1, spintype
!		do j=1, subl
!			do Vc=1, Vol/subl
!				ix=mod(Vc-1,Lx)+1
!				iy=mod((Vc-1)/Lx,Ly)+1
!				iz=(Vc-1)/(Lx*Ly) + 1
!				if( D==2 .and. iz/=1 ) then
!					write(*,*) "Err: D=2, but iz/=1", iz
!					stop
!				end if
!				iSi(ix,iy,iz) = spin((Vc-1)*subl+j,i)
!			end do
!			select case(D)
!			case(2)
!				!-- set the value of in --!
!				!do iy=1, Ly; do ix=1, Lx
!				!	fftw2_in(ix,iy) = CMPLX(iSi(ix,iy,1),0.d0)
!				!end do; end do
!				fftw2_in(:,:) = CMPLX(iSi(:,:,1),0.d0)
!				!-- FFTW --!
!				call fftw_execute_dft(plan, fftw2_in, fftw2_out)
!				!do iy=1, Ly; do ix=1, Lx
!				!	sbvq = dot_product(sublatvec(j,1:3),reclatvec(1,1:3)/Lx*(ix-1)+reclatvec(2,1:3)/Ly*(iy-1))
!				!	phase = exp(-CMPLX(0.d0,1.d0)*sbvq)
!				!	Aq(j,ix,iy,1,i) = Aq(j,ix,iy,1,i) + fftw2_out(ix,iy)*phase
!				!end do; end do
!				Aq(j,:,:,1,i) = Aq(j,:,:,1,i) + fftw2_out(:,:)
!			case(3)
!				!-- set the value of in --!
!				!do iz=1, Lz; do iy=1, Ly; do ix=1, Lx
!				!	fftw3_in(ix,iy,iz) = CMPLX(iSi(ix,iy,iz),0.d0)
!				!end do; end do; end do
!				fftw3_in(:,:,:) = CMPLX(iSi(:,:,:),0.d0)
!				!-- FFTW --!
!				call fftw_execute_dft(plan, fftw3_in, fftw3_out)
!				!do iz=1, Lz; do iy=1, Ly; do ix=1, Lx
!				!	sbvq = dot_product(sublatvec(j,1:3),reclatvec(1,1:3)/Lx*(ix-1)+reclatvec(2,1:3)/Ly*(iy-1)+reclatvec(3,1:3)/Lz*(iz-1))
!				!	phase = exp(-CMPLX(0.d0,1.d0)*sbvq)
!				!	Aq(j,ix,iy,iz,i) = Aq(j,ix,iy,iz,i) + fftw3_out(ix,iy,iz)*phase
!				!end do; end do; end do
!				Aq(j,:,:,:,i) = Aq(j,:,:,:,i) + fftw3_out(:,:,:)
!			end select
!		end do
!	end do
!
!	Sq = Sq + Aq
!
!	do iz=1, Lz; do iy=1, Ly; do ix=1, Lx
!		do sb1=1, subl;  do sb2=1, subl
!			SqSq(sb1,sb2,ix,iy,iz) = SqSq(sb1,sb2,ix,iy,iz) + dot_product(Aq(sb2,ix,iy,iz,1:spintype),Aq(sb1,ix,iy,iz,1:spintype))
!		end do;  end do
!	end do; end do; end do
!
!	deallocate(iSi)
!END SUBROUTINE structure_factor


SUBROUTINE adjust_histogram
	implicit none
	real(8)                     :: energy
	integer(4)                  :: isamp, Nisamp

	write(*,*) "--- adjust histogram -----------------------------"

	energy = 0.d0
	energy_min =  huge(1.d0)
	energy_max = -huge(1.d0)

	Nisamp = Ntoss/10

	DO isamp = 1, Nisamp
#ifdef MPI
		if( PT/=0 .and. mod(isamp,npt)==0 ) then
			call parallel_tempering
		end if
#endif
		call markov(4)
		energy = cal_energy()*wv
		energy_ave = energy_ave + energy
		if(  energy<energy_min ) energy_min = energy
		if(  energy>energy_max ) energy_max = energy
	END DO

	energy_ave = energy_ave/Nisamp
	if( abs(energy_max-energy_ave)>abs(energy_min-energy_ave) ) then
		energy_unit = abs(energy_max-energy_ave) / energy_n *2.d0
	else
		energy_unit = abs(energy_min-energy_ave) / energy_n *2.d0
	end if

	write(*,'(A30,4ES12.4)') "Energy: ave, max, min, unit:", energy_ave, energy_max, energy_min, energy_unit
	write(*,*) "---- end adjust histogram ------------------------"
END SUBROUTINE adjust_histogram


SUBROUTINE histogram
	implicit none
	integer(4)                  :: k
	real(8)                     :: dis

	dis = (totalE*wv-energy_ave)/energy_unit
	k = int( dis + sign(0.5d0,dis) )
	if( abs(k)>energy_n )   k = sign(energy_n,k)
	histogram_energy(k) = histogram_energy(k) + 1.d0

END SUBROUTINE histogram
