!========== INCLUDE MODULES ==========! 
#include "Module/mdl_vrbls.f90"
#include "Module/rng_vrbls.f90"
#include "Module/stat_vrbls.f90"
#include "Module/fftw_vrbls.f90"

#ifdef MPI
#include "Module/mpi_vrbls.f90"
! parallel tempering only used in MPI
#include "Module/pt_vrbls.f90"
#endif
!=====================================! 

PROGRAM MAIN
	use mdl_vrbls
	use rng_vrbls
	use stat_vrbls
	use fftw_vrbls
#ifdef MPI
	use mpi
	use mpi_vrbls
	use pt_vrbls
	implicit none

	integer(4)                  :: id=10             ! the idle of file
	integer                     :: stat(MPI_STATUS_SIZE)
	character(LEN=MPI_MAX_PROCESSOR_NAME) :: processor_name

	call MPI_INIT(ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
	call MPI_Get_processor_name(processor_name, namelen, ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierr)
	!call MPI_COMM_PROCESSOR_NAME(processor_name, namelen, ierr)
	call date_and_time(date, time, zone, tval)
	iseed0 = tval(2)*31*24 + tval(3)*24
	iseed0 = (iseed0+tval(5))*3600+tval(6)*60+tval(7)
	iseed0 = iseed0+taskid*3600
	write(*,*) processor_name, namelen
	do itag=1, len(trim(processor_name))
		iseed0 = iseed0 + IACHAR(processor_name(itag:itag))
	end do
	!iseed0 = taskid*100 + 1
	write(*,*) numprocs, taskid, iseed0

	!--- build cxxx for mpi run ---!
	call Build
#endif

!================================= SIMULATION =================================!
	call Ising
!==============================================================================!

#ifdef MPI
	if( taskid /= 0 ) then
		itag = 1
		dest = 0
		call MPI_SEND(taskid, 1, MPI_INTEGER, dest, itag, MPI_COMM_WORLD, ierr)
		!TYPE(*), DIMENSION(..), INTENT(IN) :: buf
		!INTEGER, INTENT(IN) :: count, dest, tag
		!TYPE(MPI_Datatype), INTENT(IN) :: datatype
		!TYPE(MPI_Comm), INTENT(IN) :: comm
		!INTEGER, OPTIONAL, INTENT(OUT) :: ierror
	else
		open(id,file="../log",access="append")
		write(id,'(A5,I3,A4,I3,A17)') "task ", taskid, " of ", numprocs, " process is done."
		do source=1, numprocs-1
			itag = 1
			call MPI_RECV(impi, 1, MPI_INTEGER, MPI_ANY_SOURCE, itag, MPI_COMM_WORLD, stat, ierr)
			!TYPE(*), DIMENSION(..) :: buf
			!INTEGER, INTENT(IN) :: count, source, tag
			!TYPE(MPI_Datatype), INTENT(IN) :: datatype
			!TYPE(MPI_Comm), INTENT(IN) :: comm
			!TYPE(MPI_Status) :: status
			!INTEGER, OPTIONAL, INTENT(OUT) :: ierror
			write(id,'(A5,I3,A4,I3,A17)') "task ", stat(MPI_SOURCE), " of ", numprocs, " process is done."
		end do
		close(id)
	end if
	!open(id,file="../log",access="append")
	!write(id,'(A5,I3,A4,I3,A17)') "task ", taskid, " of ", numprocs, " process is done."
	!close(id)
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	call MPI_FINALIZE(ierr)
#endif

	stop

CONTAINS

#include "Ising.f90"

#ifdef MPI
#include "Build/Build.f90"
#include "Parallel_Tempering/parallel_tempering.f90"
#endif

#include "Initialize/initialize.f90"
#include "Initialize/lattice.f90"
#include "Initialize/my_fftw_destroy.f90"

#include "Markov/markov.f90"
#include "Markov/common.f90"
#include "Markov/worm.f90"
#include "Markov/sample_sequence.f90"

#include "Measure/measure.f90"
#include "Measure/observables.f90"
#include "Measure/specific.f90"

#include "Statistics/statistics.f90"

#include "WriteRead/wr_cnf_stat.f90"
#include "RNG/rng.f90"

END PROGRAM MAIN

