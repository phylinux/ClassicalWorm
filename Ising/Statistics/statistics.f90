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
	ptb = NBlck/2+1
	pts = 1

	return
END SUBROUTINE coarsen_data
!===================================================================
