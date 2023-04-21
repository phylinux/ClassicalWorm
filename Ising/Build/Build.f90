SUBROUTINE Build
	implicit none
	integer                     :: stat(MPI_STATUS_SIZE)
	integer                     :: stat_chdir

	character(10)               :: dname
#ifdef IFORT
	integer, external           :: chdir
#endif
	character(50)               :: chartemp
	real(8)                     :: temp
	logical                     :: existyn
	integer(4)                  :: id
	character(50)               :: parafile
	integer(4)                  :: i

	id = 10
	parafile="parafile.txt"

	!--- send the parameter to each core ---------------------------------------!
	itag=1
	if( taskid==0 ) then
		inquire(file=trim(parafile), exist=existyn)
		if( .not. existyn ) then
			write(*,*) "Err: no '"//trim(parafile)//"'"
			call MPI_ABORT(MPI_COMM_WORLD,999,ierr)
		end if

		open(id,file=trim(parafile),action="read")
		do i=0,numprocs-1
			read(id,*) dest, temp
			if( dest==0 ) then
				write(chartemp,*) temp
			else
				call MPI_SEND(temp, 1, MPI_REAL8, dest, itag, MPI_COMM_WORLD, ierr)
			end if
		end do
		close(id)
	else
		source=0
		itag=1
		call MPI_RECV(temp, 1, MPI_REAL8, source, itag, MPI_COMM_WORLD, stat, ierr)
		write(chartemp,*) temp
	end if
	!---------------------------------------------------------------------------!

	!--- build cxxx for each core ----------------------------------------------!
	if ( taskid < 10 ) then
		write(dname,'("c00",I1)') taskid
	else if ( taskid < 100 ) then
		write(dname,'("c0",I2)') taskid
	else
		write(dname,'("c",I3)') taskid
	end if

	open(id,file="input",action="read")
	read(id,*) r_cnf_stat
	close(id)
	if( r_cnf_stat==0 ) then
		call system('mkdir '//trim(dname))
	end if
	!---------------------------------------------------------------------------!

#ifdef IFORT
	stat_chdir = chdir(trim(dname))
#endif

#ifdef GFORT
	call chdir(trim(dname), stat_chdir)
#endif

	inquire(file="inp00", exist=existyn)
	if( existyn .and. r_cnf_stat==0 ) return

	if ( stat_chdir /= 0 ) then
		print*, "Err: go to directory "//trim(dname)//" failure."
		stop
	end if

	if( r_cnf_stat==0 ) then
		call system('sed "s#@beta#'//trim(chartemp)//'#" ../input    > ./inp00')

		call system('sed "s#@Jcp1#'//trim(chartemp)//'#"     ./inp00      > temp'); call system('mv temp inp00')
		!------------------------------------!
		open(id,file='inp00',access="append")
		write(id,*) iseed0
		close(id)
	end if

	return
END SUBROUTINE Build
