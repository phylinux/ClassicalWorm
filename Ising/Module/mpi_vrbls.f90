MODULE mpi_vrbls
	implicit none
	integer,      parameter           :: master=0
	integer,      parameter           :: one=1
	integer                           :: psum
	integer                           :: ptime
	integer                           :: taskid, numprocs, ierr
	integer                           :: namelen
	integer                           :: iseed0
	integer                           :: dest, source
	integer                           :: itag
	integer                           :: impi

	real(8)                           :: taskparalist(0:1000)
END MODULE mpi_vrbls
