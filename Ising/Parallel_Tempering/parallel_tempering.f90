SUBROUTINE parallel_tempering
	implicit none
	integer                     :: stat(MPI_STATUS_SIZE)

	integer(4)                  :: root
	real(8)                     :: ibeta, jbeta
	real(8)                     :: Ei, Ej
	integer(4)                  :: exchange(0:1000)
	integer(4)                  :: loop(0:1000)
	integer(4)                  :: lphd(0:1000,2)
	integer(4)                  :: nl, lstart
	integer(4)                  :: prevhead, nexthead
	integer(4)                  :: sendlist(0:1000), recvlist(0:1000)
	integer(4)                  :: flag, useless
	real(8)                     :: p
	integer(4)                  :: i, j

	if( taskid==0 )  then
		loop = -1
		lphd = -1
		exchange = 0
	end if

	!if( taskid==0 ) then
	!	open(1,file="energy.list",action="write")
	!	write(1,'(I4,ES16.8)') taskid, totalE/beta
	!	do source=1, numprocs-1
	!		itag = 101
	!		call MPI_RECV(Ei, 1, MPI_REAL8, source, itag, MPI_COMM_WORLD, stat, ierr)
	!		write(1,'(I4,ES16.8)') source, Ei
	!	end do
	!	close(1)
	!else
	!	!--- send energy to taskid=0 ---------------------------
	!	dest = 0
	!	itag = 101
	!	call MPI_SEND(totalE/beta, 1, MPI_REAL8, dest, itag, MPI_COMM_WORLD, ierr)
	!	!---------------------------------------------------------------------------
	!end if

	!write(*,*) "taskid ", taskid, " arrive"
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	totalE = cal_energy()
	if( taskid==0 ) then
		exprob(0) = exprob(0) + 1.d0
		dest = taskid+1
		!--- send spin, beta, energy to the next process ---------------------------
		itag = 101
		call MPI_SEND(beta, 1, MPI_REAL8, dest, itag, MPI_COMM_WORLD, ierr)
		itag = 102
		call MPI_SEND(totalE/beta, 1, MPI_REAL8, dest, itag, MPI_COMM_WORLD, ierr)
		!---------------------------------------------------------------------------
		!write(*,*) "taskid 0 send data to taskid 1"
	else
		!--- receive spin, beta, energy from the previous process ---------------------------
		source = taskid-1
		!write(*,*) "taskid ",taskid," start to receive data from taskid ",source
		itag = 101
		call MPI_RECV(ibeta, 1, MPI_REAL8, source, itag, MPI_COMM_WORLD, stat, ierr)
		itag = 102
		call MPI_RECV(Ei, 1, MPI_REAL8, source, itag, MPI_COMM_WORLD, stat, ierr)
		!write(*,*) "taskid ",taskid," receive data from taskid ",source,"done"
		!------------------------------------------------------------------------------------
		jbeta = beta
		Ej    = totalE/jbeta
		p = exp( (ibeta-jbeta)*(Ei-Ej) )
		if( rn()<p ) then
			totalE = Ei*jbeta     !!! must update energy
			p = 1.d0
		else
			p = 0.d0
		end if
		!--- send the accepted information to process 0 --------------------------------
		dest = 0
		itag = 103
		!write(*,*) "taskid ",taskid," send p to taskid ",dest
		call MPI_SEND(p, 1, MPI_REAL8, dest, itag, MPI_COMM_WORLD, ierr)
		!-------------------------------------------------------------------------------

		if( taskid<numprocs-1 ) then
			!--- send spin, beta, energy to the next process ---------------------------
			dest = taskid+1
			!write(*,*) "taskid ",taskid," send data to taskid ",dest
			itag = 101
			call MPI_SEND(beta, 1, MPI_REAL8, dest, itag, MPI_COMM_WORLD, ierr)
			itag = 102
			call MPI_SEND(totalE/beta, 1, MPI_REAL8, dest, itag, MPI_COMM_WORLD, ierr)
			!---------------------------------------------------------------------------
		end if
	end if

	do i=0,numprocs-1
		recvlist(i) = i
		sendlist(i) = i
	end do
	if( taskid==0 ) then
		do source=1, numprocs-1
			!--- receive accepted information from other process --------------------------------
			itag = 103
			!write(*,*) "taskid ",taskid," receive p from taskid ",source
			call MPI_RECV(p, 1, MPI_REAL8, source, itag, MPI_COMM_WORLD, stat, ierr)
			!------------------------------------------------------------------------------------
			exprob(source) = exprob(source) + p
			!--- construct the exchange list of conf ---!
			if( p>0.5d0 ) then
				exchange(source) = 1
				i = recvlist(source-1)
				j = recvlist(source)
				sendlist(i) = source
				sendlist(j) = source-1

				i = recvlist(source-1)
				recvlist(source-1) = recvlist(source)
				recvlist(source) = i
			end if
		end do
		!write(*,*) recvlist(0:numprocs-1)
		!write(*,*) sendlist(0:numprocs-1)
	end if

	!--- construct transfer loop ---!
	! -1: single conf; n: the head of a loop; 0: the other conf in the loop
	if( taskid==0 ) then
		nl = 1
		prevhead = -1
		lstart = 0
		do i=0, numprocs-1
			if( lstart==0 .and. exchange(i)==1 ) then
				loop(i-1)=nl            ! the head of loop
				loop(i) = 0
				lphd(i-1,1)=prevhead    ! the head of the previous loop
				if( prevhead/=-1 ) lphd(prevhead,2)=i-1    ! the head of the next loop
				prevhead = i-1
				lstart = 1
				nl = nl+1
			else if( lstart==1 .and. exchange(i)==1 ) then
				loop(i) = 0
			else
				lstart = 0
			end if
		end do
	end if

	root = 0
	call MPI_BCAST(recvlist(0), numprocs, MPI_INTEGER4, root, MPI_COMM_WORLD, ierr)
	call MPI_BCAST(sendlist(0), numprocs, MPI_INTEGER4, root, MPI_COMM_WORLD, ierr)
	call MPI_BCAST(loop(0), numprocs, MPI_INTEGER4, root, MPI_COMM_WORLD, ierr)
	call MPI_BCAST(lphd(0,1), numprocs, MPI_INTEGER4, root, MPI_COMM_WORLD, ierr)
	call MPI_BCAST(lphd(0,2), numprocs, MPI_INTEGER4, root, MPI_COMM_WORLD, ierr)

	!--- exchange conf ---!
	flag = loop(taskid)
	dest   = sendlist(taskid)
	source = recvlist(taskid)
	prevhead = lphd(taskid,1)
	nexthead = lphd(taskid,2)
	itag = 100
	useless = 0
	if( flag>0 ) then
		if( prevhead/=-1 ) call MPI_RECV(useless, 1, MPI_INTEGER4, prevhead, itag, MPI_COMM_WORLD, stat, ierr)
		call MPI_SEND(spin(1,1), Vol*3, MPI_COMPLEX, dest, itag, MPI_COMM_WORLD, ierr)
		call MPI_RECV(spin(1,1), Vol*3, MPI_COMPLEX, source, itag, MPI_COMM_WORLD, stat, ierr)
		if( nexthead/=-1 ) call MPI_SEND(useless, 1, MPI_INTEGER4, nexthead, itag, MPI_COMM_WORLD, ierr)
	else if( flag==0 ) then
		newconf = spin
		call MPI_RECV(spin(1,1), Vol*3, MPI_COMPLEX, source, itag, MPI_COMM_WORLD, stat, ierr)
		call MPI_SEND(newconf(1,1), Vol*3, MPI_COMPLEX, dest, itag, MPI_COMM_WORLD, ierr)
	end if

	!--- in ith parallel tempering, the energy of rank i-1 is not be updated ---!
	totalE = cal_energy()
	
	!----------------------------------------------------------
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	!----------------------------------------------------------
	!write(*,*) "parallel tempering done"
END SUBROUTINE parallel_tempering

SUBROUTINE iteration
	implicit none
	integer                     :: stat(MPI_STATUS_SIZE)
	real(8)                     :: ibeta, prevbeta
	real(8)                     :: avepm, pm

	converge = 0.d0

	if( taskid==0 ) then
		open(1,file="exchange_prob1.txt",action="write")
		do source=1, numprocs-1
			write(1,'(I4,F12.6)') source, exprob(source)/exprob(0)
		end do
		close(1)
	end if

	if( taskid/=0 ) then
		itag = 101
		dest = 0
		call MPI_SEND(beta, 1, MPI_REAL8, dest, itag, MPI_COMM_WORLD, ierr)
		source = 0
		call MPI_RECV(beta, 1, MPI_REAL8, source, itag, MPI_COMM_WORLD, stat, ierr)
		!--- calculate betaJ and energy for the new beta ---!
		call redefine_para(beta)
		totalE = cal_energy()
	else
		betalist(0) = beta
		itag = 101
		do source=1, numprocs-1
			call MPI_RECV(betalist(source), 1, MPI_REAL8, source, itag, MPI_COMM_WORLD, stat, ierr)
		end do

		avepm = sum(exprob(1:numprocs-1))/(numprocs-1.d0) / exprob(0)
		prevbeta = beta
		do dest=1, numprocs-1
			pm = exprob(dest) / exprob(0)
			if( pm<1.d-7 ) then
				ibeta = (betalist(dest-1) + betalist(dest))*0.5d0
			else
				ibeta = betalist(dest-1) + (betalist(dest)-prevbeta)*pm/avepm
			end if
			call MPI_SEND(ibeta, 1, MPI_REAL8, dest, itag, MPI_COMM_WORLD, ierr)
			prevbeta = betalist(dest)
			betalist(dest) = ibeta
			converge = converge + abs(prevbeta-ibeta)
		end do
		converge = converge/(numprocs-1)
	end if

	if( taskid==0 ) then
		open(1,file="betalist.txt",action="write")
		do source=0, numprocs-1
			write(1,'(I4,F12.6)') source, betalist(source)
		end do
		close(1)
	end if

	!----------------------------------------------------------
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	!----------------------------------------------------------
END SUBROUTINE iteration

SUBROUTINE adjust_parallel_tempering
	implicit none
	integer(4)                  :: root=0
	integer(4)                  :: iblck,isamp
	integer(4)                  :: cn

	converge = 1.d0
	cn = 0
	do while( converge>eps_con )
		exprob = 0.d0
		do iblck=1, 400
			do isamp = 1, 50
				call markov(nw)
			end do
			call parallel_tempering
		end do
		call iteration
		cn = cn+1
		call MPI_BCAST(converge, 1, MPI_REAL8, root, MPI_COMM_WORLD, ierr)
		if(taskid==0) write(*,'(A10,I4,F10.4)') "iteration ", cn, converge
	end do
	write(*,'(A7,I3,A20)')  "taskid ", taskid, " --- iteration done."
	exprob = 0.d0
END SUBROUTINE adjust_parallel_tempering

SUBROUTINE exchange_probability
	implicit none
	integer                     :: stat(MPI_STATUS_SIZE)
	real(8)                     :: ibeta, prevbeta
	real(8)                     :: avepm, pm

	if( taskid==0 ) then
		open(1,file="exchange_prob2.txt",action="write")
		do source=1, numprocs-1
			write(1,'(I4,F12.6)') source, exprob(source)/exprob(0)
		end do
		close(1)
	end if
END SUBROUTINE exchange_probability
