
!--- PROJECT-INDEPENDENT -------------------------------------------
!==============Collect data ========================================
!! THIS IS PROJECT-INDEPENDENT 
SUBROUTINE coll_data(iblck)
	implicit none
	integer, intent(in) :: iblck
	integer             :: j
	do j = 1, NObs_b
		Obs(j,iblck) = Obs(j,iblck)+ Quan(j)
	enddo
	gr(:,iblck) = gr(:,iblck)+gri(:)
	gri = 0.d0
END SUBROUTINE coll_data 

!==============Normalize by Nsamp ==================================
!! THIS IS PROJECT-INDEPENDENT 
SUBROUTINE norm_Nsamp(iblck)
	implicit none
	integer, intent(in) :: iblck
	integer             :: j
	double precision    :: nor
	nor = 1.d0/(Nsamp*1.d0)
	do j = 1, NObs_b
		Obs(j,iblck) = Obs(j,iblck)*nor
	enddo
	gr(:,iblck) = gr(:,iblck)*nor
END SUBROUTINE norm_Nsamp
!===================================================================

!==============Statistics ==========================================
!! THIS IS PROJECT-INDEPENDENT 
SUBROUTINE stat_analy(iblck)
	implicit none
	integer          :: iblck
	integer          :: j, k, k0
	double precision :: devn, devp, nor

	! -- calculate average -------------------------------------------
	nor  = 1.d0/(iblck*1.d0)
	do j = 1, NObs_b
		Ave(j) = nor*Sum(Obs(j,1:iblck))
	enddo

	!Coarsen: do
		! -- calculate error and t=1 correlation for basics obs.--------
		prt = .true.
		DO j = 1, NObs_b
			devp = 0.d0;  Cor(j) = 0.d0;  Dev(j) = 0.d0
			do k = 1,  iblck
				devn   = Obs(j,k)-Ave(j)
				Dev(j) = Dev(j)+devn*devn
				Cor(j) = Cor(j)+devn*devp
				devp   = devn
			enddo 
			Dev(j)   = Dev(j)*nor;        Cor(j) = Cor(j)*nor
			if(Dev(j)>eps)                Cor(j) = Cor(j)/Dev(j)
			Dev(j)   = dsqrt(Dev(j)/(iblck-1.d0))
			if(dabs(Cor(j))>tol) prt = .false.
		ENDDO 

		!IF(prt)                         EXIT Coarsen 
		!IF(NBlck<=64)    THEN
		!	prt = .false.;                EXIT Coarsen 
		!ENDIF

	!	! -- coarsen blocking ------------------------------------------
	!	nor   = nor*2.d0
	!	call coarsen_data
	!enddo Coarsen 

	! -- define auxillary variables and average of composite obs.-----
	call cal_Obs_comp(iblck)

	! -- calculate error and t=1 correlation for composite obs.-----
	do j = 1+NObs_b, NObs
		devp = 0.d0;  Cor(j) = 0.d0;  Dev(j) = 0.d0
		DO k = 1,  iblck
			devn   = Obs(j,k)-Ave(j)
			Dev(j) = Dev(j)+devn*devn
			Cor(j) = Cor(j)+devn*devp
			devp   = devn
		ENDDO
		Dev(j)   = Dev(j)*nor;        Cor(j) = Cor(j)*nor
		IF(Dev(j)>eps)                Cor(j) = Cor(j)/Dev(j)
		Dev(j)   = dsqrt(Dev(j)/(iblck-1.d0))
	enddo
	return
END SUBROUTINE stat_analy
!===================================================================

!============== Write to files =====================================
!! THIS IS PROJECT-DEPENDENT 
SUBROUTINE write2file 
	IMPLICIT NONE
	integer       :: j, k, Nwri
	double precision :: t_tot

	!-- open file ----------------------------------------------------
	open (1,file=datafile,  access='append') 
	write(1, *) "===================================================="

	!-- write to data file--------------------------------------------
	write(1,40) ident, D, Lx, Ly, Lz, Totsamp, nw, betaJ, Seed
	40 format(a10,i2,3i5,i12,i8,f14.8,i8)

	do j = 1, Nobs
		write(1,41) j, Ave(j), Dev(j), Cor(j)
		write(6,41) j, Ave(j), Dev(j), Cor(j)
		41 format(i5,2es20.8,f12.5)
	enddo
	write(1,42) 'PRT= ', prt, 'NBlck= ', NBlck, 'Nsamp= ', Nsamp
	write(6,42) 'PRT= ', prt, 'NBlck= ', NBlck, 'Nsamp= ', Nsamp
	42 format(A5,l5,1x,A7,I5,1x,A7,I8,1x)
	close(1)
	!-- green function --
	open(1,file='green.txt',action='write')
	write(1,'(I4,3I6,I8)') D, Lx, Ly, Lz, NBlck
	do j=1, NBlck
		do k=0, gVol
			write(1,'(I6,ES16.8)') k, gr(k,j)/gr(0,j)
		end do
	end do

	!-- write to output file if #block is too small-------------------
	!if(NBlck<=64) then
	!	write(6,*)
	!	Nwri = NObs_b;       if(Nwri>5) Nwri = 5
	!	do k = 1, NBlck
	!		write(6,42) k,(Obs(j,k),j=1,Nwri)
	!		42 format(i4,5f16.8) 
	!	end do
	!endif

	!t_toss = t_toss/60.d0 
	!t_simu = t_simu/60.d0
	!t_meas = t_meas/60.d0
	t_tot  = t_toss+t_simu+t_meas
	write(6,50) (t_toss/60.d0), (t_toss/3600.d0)
	50 format( '  equilibrate time:',f14.7,2x,'m. or ',f12.7,'  h.')
	write(6,51) (t_simu/60.d0), (t_simu/3600.d0)
	51 format( '  simulation  time:',f14.7,2x,'m. or ',f12.7,'  h.')
	write(6,52) (t_meas/60.d0), (t_meas/3600.d0)
	52 format( '  measure     time:',f14.7,2x,'m. or ',f12.7,'  h.')
	write(6,53) (t_tot/60.d0),  (t_tot/3600.d0)
	53 format( '  total CPU   time:',f14.7,2x,'m. or ',f12.7,'  h.')
	return
END SUBROUTINE write2file
!===================================================================

!============== Coarsen data =======================================
SUBROUTINE coarsen_data
	IMPLICIT NONE
	integer     :: j, k

	NBlck = NBlck/2
	DO j = 1, NObs_b
		do k   = 1, NBlck
			Obs(j,k) = (Obs(j,2*k-1)+Obs(j,2*k))*0.5d0
		enddo 
	ENDDO 
	Obs(1:NObs_b,NBlck+1:NBlck*2) = 0.d0
	Nsamp = Nsamp*2
	NBlck = NBlck*2
	Totsamp = Nsamp*NBlck/1000

	return
END SUBROUTINE coarsen_data
!===================================================================
