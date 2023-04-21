!============== Write to files =====================================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE write2file(iblck)
	IMPLICIT NONE
	integer(4), intent(in)      :: iblck
	integer(4)                  :: i, j, k, Nwri
	real(8)                     :: t_tot
	integer(4)                  :: id

	id=1

	!-- open file ----------------------------------------------------
	!open (id,file=trim(datafile),  access="append")
	open (id,file=trim(datafile),  action="write")
	write(id, *) "===================================================="

	!-- write to data file ------------------------------------------------
	!-- quantities --
	write(id,"(A16,I6,2F16.8,I9,I3,I9)") trim(ident), Lx, beta, Jcp(1), Totsamp, nw, Seed

	do i = 1, Nobs
		write(id,'(I5,2ES20.8,F12.5)') i, Ave(i), Dev(i), Cor(i)
		write(* ,'(I5,2ES20.8,F12.5)') i, Ave(i), Dev(i), Cor(i)
	enddo

	t_tot  = t_init+t_adjt+t_prtp+t_ther+t_simu+t_meas
	write(id,'(A5,L5,1x,A7,I5,1x,A7,I8,1x,A7,I5)') 'PRT= ', prt, 'NBlck= ', NBlck, 'Nsamp= ', Nsamp, "iblck=", iblck
	write(* ,'(A5,L5,1x,A7,I5,1x,A7,I8,1x,A7,I5)') 'PRT= ', prt, 'NBlck= ', NBlck, 'Nsamp= ', Nsamp, "iblck=", iblck
	write(id,'(A23,F10.2,A4)') "Initialization time = ", t_init/60.d0, "(m)"
	write(id,'(A23,F10.2,A4)') "Thermalization time = ", t_ther/60.d0, "(m)"
	write(id,'(A23,F10.2,A4)') "Adjust         time = ", t_adjt/60.d0, "(m)"
	write(id,'(A23,F10.2,A4)') "Parallel tempering  = ", t_prtp/60.d0, "(m)"
	write(id,'(A23,F10.2,A4)') "Simulation     time = ", t_simu/60.d0, "(m)"
	write(id,'(A23,F10.2,A4)') "Measure        time = ", t_meas/60.d0, "(m)"
	write(id,'(A23,F10.2,A4)') "Total          time = ", t_tot/60.d0, "(m)"
	close(id)

	!-- histogram, correlation, structure factor --
	!call write_structure_factor(iblck)
	!call write_realspace_correlation(iblck)
#ifdef HISTOGRAM
	call write_histogram(iblck)
#endif
	!----------------------------------------------------------------------

	!-- write to output file if #block is too small-------------------
	!if(NBlck<=64) then
	!	write(6,*)
	!	Nwri = NObs_b;       if(Nwri>5) Nwri = 5
	!	do k = 1, NBlck
	!		write(6,42) k,(Obs(j,k),j=1,Nwri)
	!		42 format(i4,5f16.8) 
	!	end do
	!endif

	!!-- the accepted ratio ------------------------------------------------
	!write(*,'(A3)') "---"
	!write(*,'(A23,F10.5)') "Total accepted prob = ", accept(0,2)*1.d0/accept(0,1)
	!do i=1, 2
	!	write(*,'(A4,I2.2,A17,F10.5)') "Op ",i," accepted prob = ", accept(i,2)*1.d0/accept(i,1)
	!end do

	!t_ther = t_ther/60.d0 
	!t_simu = t_simu/60.d0
	!t_meas = t_meas/60.d0
	!t_tot  = t_ther+t_simu+t_meas
	write(*,'(A3)') "---"
	write(*,'(A23,F10.2,A4)') "Initialization time = ", t_init/60.d0, "(m)"
	write(*,'(A23,F10.2,A4)') "Thermalization time = ", t_ther/60.d0, "(m)"
	write(*,'(A23,F10.2,A4)') "Adjust         time = ", t_adjt/60.d0, "(m)"
	write(*,'(A23,F10.2,A4)') "Parallel tempering  = ", t_prtp/60.d0, "(m)"
	write(*,'(A23,F10.2,A4)') "Simulation     time = ", t_simu/60.d0, "(m)"
	write(*,'(A23,F10.2,A4)') "Measure        time = ", t_meas/60.d0, "(m)"
	write(*,'(A23,F10.2,A4)') "Total          time = ", t_tot/60.d0, "(m)"
END SUBROUTINE write2file
!===================================================================

SUBROUTINE write_structure_factor(iblck)
	implicit none
	integer(4), intent(in)      :: iblck
	complex(8), allocatable     :: Sabq(:,:,:,:,:)    ! Sabq(a,b,kx,ky,kz)
	integer(4)                  :: sb1, sb2
	complex(8)                  :: ccmplx
	integer(4)                  :: ix, iy, iz
	real(8)                     :: nor

	nor = 1.d0/iblck
	nor = nor/Nsamp

	allocate(Sabq(subl,subl,Lx,Ly,Lz))
	Sabq = CMPLX(0.d0,0.d0)

	!--- write spin configuration --!
	open(1,file="Sr.txt",action="write")
	do ix=1, Vol
		write(1,"(3ES16.8)") Sr(ix,1:3)*nor
	end do
	close(1)

	!--- write spin configuration -----!
	open(1,file="spinconf.txt",access="append")
	write(1,'(A4,I8)') "# ", iblck
	do ix=1, Vol
		write(1,"(3ES16.8)") spin(ix,1:3)
	end do
	close(1)


	do iz=1, Lz; do iy=1, Ly; do ix=1, Lx
		do sb1=1,subl; do sb2=1, subl
			!ccmplx = SqSq(sb1,sb2,ix,iy,iz) - dot_product(Sq(sb2,ix,iy,iz,1:spintype),(Sq(sb1,ix,iy,iz,1:spintype)) * nor
			ccmplx = SqSq(sb1,sb2, ix,iy,iz)
			Sabq(sb1,sb2, ix,iy,iz) = ccmplx
		end do; end do
	end do; end do; end do

	Sabq(:,:,:,:,:) = Sabq(:,:,:,:,:) * nor / Vol


	!--- write SqSq --!
	open(1,file="SqSq.txt",action="write")
	write(1,'(I4)') 0
	do iz=1, Lz; do iy=1, Ly; do ix=1, Lx
		do sb2=1,subl; do sb1=1, subl
			write(1,'(2ES16.8)') Real(Sabq(sb1,sb2,ix,iy,iz)), Aimag(Sabq(sb1,sb2,ix,iy,iz))
		end do; end do
	end do; end do; end do
	close(1)

	deallocate(Sabq)
END SUBROUTINE write_structure_factor

SUBROUTINE write_realspace_correlation(iblck)
	implicit none
	integer(4), intent(in)      :: iblck
	integer(4)                  :: dr
	real(8)                     :: nor

	nor = 1.d0/iblck
	nor = nor/Nsamp

	open(1,file="SrSr.txt",action="write")
	write(1,'(I4,E16.8)') 0, 1.d0
	do dr=1, Lx/2
		write(1,'(I4,E16.8)') dr, SrSr(dr)*nor
	end do
	close(1)

END SUBROUTINE write_realspace_correlation


SUBROUTINE write_histogram(iblck)
	implicit none
	integer(4), intent(in)      :: iblck
	integer(4)                  :: i
	real(8)                     :: energy
	real(8)                     :: nor

	nor = 1.d0/sum(histogram_energy(-energy_n:energy_n))

	open(1,file="histogram_energy.txt",action="write")
	write(1,'(A20)') "# energy histogram"
	do i = -energy_n, energy_n
		energy = energy_ave + energy_unit*i
		write(1,'(I6,2E16.8)') i, energy, histogram_energy(i)*nor
	end do
	close(1)

END SUBROUTINE write_histogram


!-- write cnf and stat
subroutine write_cnf_stat
	implicit none

	open(1,file="cnf.dat",form="UNFORMATTED",access="sequential",status="replace")
	write(1) D
	write(1) Lx
	write(1) Ly
	write(1) Lz
	write(1) Jcp
	write(1) beta

	write(1) spin
	write(1) current
#ifdef MPI
	write(1) exprob
	write(1) betalist
#endif

	!-- close file --!
	close(1)
	write(*,'(A23)') "write cnf done"

	open(1,file="stat.dat",form="UNFORMATTED",access="sequential",status="replace")
	!-- independent --!
	write(1) ptb
	write(1) NBlck
	write(1) pts
	write(1) Nsamp
	write(1) Obs

	!-- close file --!
	close(1)
	write(*,'(A23)') "write stat done"

	call write_rng
	write(*,'(A23)') "write rng done"
	write(*,'(A23)') "-------------"

end subroutine write_cnf_stat

!-- read cnf and stat
subroutine read_cnf_stat(rc,rs)
	implicit none
	integer(4), intent(in)      :: rc,rs
	integer(4)                  :: D1
	integer(4)                  :: Lx1, Ly1, Lz1
	real(8)                     :: beta1, Jcp1(maxir)

	if( rc==1 ) then
		open(1,file="cnf.dat",form="UNFORMATTED",access="sequential",status="old")
		read(1) D1;      if(D1/=D)       stop "Err: D"
		read(1) Lx1;     if(Lx1/=Lx)     stop "Err: Lx"
		read(1) Ly1;     if(Ly1/=Ly)     stop "Err: Ly"
		read(1) Lz1;     if(Lz1/=Lz)     stop "Err: Lz"
		read(1) Jcp1;    !if(dot_product(Jcp1(1:nnr)-Jcp(1:nnr), Jcp1(1:nnr)-Jcp(1:nnr))>1.d-5)   stop "Err: Jcp"
		read(1) beta1

		read(1) spin
		read(1) current
#ifdef MPI
		read(1) exprob
		read(1) betalist
#endif

		!-- close file --!
		close(1)
		write(*,'(A23)') "read cnf done"
	end if

	if( rs==1 ) then
		open(1,file="stat.dat",form="UNFORMATTED",access="sequential",status="old")
		!-- independent --!
		read(1) ptb
		read(1) NBlck
		read(1) pts
		read(1) Nsamp
		read(1) Obs

		if( ptb==NBlck .and. pts==Nsamp ) then
			call coarsen_data
		else if ( pts<Nsamp ) then
			pts = pts+1
		else if ( pts==Nsamp ) then
			if ( ptb >= NBlck ) stop "Err: ptb"
			ptb = ptb+1
			pts = 1
		else
			stop "Err: pts"
		end if

		!-- close file --!
		close(1)
		write(*,'(A23)') "read stat done"

		call read_rng
		write(*,'(A23)') "read rng done"
		write(*,'(A23)') "-------------"

	end if

end subroutine read_cnf_stat

!---------------------------- store rng ----------------------------------------!
subroutine write_rng
	implicit none
	open(1,file="rng.dat",form="UNFORMATTED",access="sequential",status="replace")
	write(1) Seed
	write(1) nrannr
	write(1) ipnt1, ipnf1
	write(1) ipnt2, ipnf2
	write(1) inxt1
	write(1) inxt2
	write(1) ir1
	write(1) ir2
	write(1) irn
	close(1)
end subroutine write_rng
!-- restore rng --!
subroutine read_rng
	implicit none
	integer                     :: Seed1
	open(1,file="rng.dat",form="UNFORMATTED",access="sequential",status="old")
	read(1) Seed1; if( Seed1/=Seed ) stop "Err: Seed1/=Seed"
	read(1) nrannr
	read(1) ipnt1, ipnf1
	read(1) ipnt2, ipnf2
	read(1) inxt1
	read(1) inxt2
	read(1) ir1
	read(1) ir2
	read(1) irn
	close(1)
end subroutine read_rng
