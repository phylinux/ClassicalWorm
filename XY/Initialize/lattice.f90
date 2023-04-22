SUBROUTINE read_latt
	implicit none

	integer                     :: stat
#ifdef IFORT
	integer, external           :: chdir
	integer, external           :: getcwd
#endif
	character(200)              :: current_dir
	character(200)              :: lattice_dir
	integer(4)                  :: D1, Lx1, Ly1, Lz1
	integer(4)                  :: subl1
	character(50)               :: charnnsite, charbackdir, charbond
	integer(4)                  :: Vc
	integer(4)                  :: nr
	integer(4)                  :: i, j
	integer(4)                  :: si, sj
	integer(4), allocatable     :: nnm(:), bd(:)
	real(8)                     :: scoord(3)
	integer(4)                  :: ii(6)
	integer(4), allocatable     :: bondcount(:)

#ifdef IFORT
	stat = getcwd(current_dir)
#endif

#ifdef GFORT
	call getcwd(current_dir,stat)
#endif

	write(*,*) "current dir = ", trim(current_dir)
	if ( stat /= 0 ) then
		print*, "Err: get current directory path failure."
		stop
	end if
	open(1,file="../directory.dat",action="read")
	read(1,'(A)') lattice_dir
	close(1)
	write(*,*) "lattice dir = ", trim(lattice_dir)

#ifdef IFORT
	stat = chdir(trim(lattice_dir))
#endif

#ifdef GFORT
	call chdir(trim(lattice_dir), stat)
#endif

	if ( stat /= 0 ) then
		print*, "Err: go to directory "//trim(lattice_dir)//" failure."
		stop
	end if


	!-- lattice info --------------------------
	open(1,file="lattice",action="read")
	read(1,*) D1, Lx1, Ly1, Lz1, subl1
	do i=1, nnr
		read(1,*) nnb(i)
	end do
	close(1)
	if( D1/=D .or. Lx1/=Lx .or. Ly1/=Ly .or. Lz1/=Lz .or. subl1/=subl ) stop "Err: lattice"

	Vol = Lx*Ly*Lz*subl
	wv  = 1.d0/Vol

	allocate(Scoordinate(Vol,3))
	allocate(Icoordinate(Vol,3))
	allocate(sublattice(Vol))
	allocate(nnsite(nnr,Vol,maxval(nnb)))
	allocate(backdir(nnr,Vol,maxval(nnb)))
	allocate(sublatvec(subl,3))
	allocate(reclatvec(3,3))
	nnsite = -1

	allocate(positivec(nnr,Vol,maxval(nnb)))

	!-- site info -----------------------------
	open(1,file="site.dat",action='read')
	read(1,*) Vc
	if( Vc/=Vol ) stop "Err: Vc/=Vol"
	do i=1, Vol
		read(1,*) Vc, ii(1:6), scoord(1:3)
		Scoordinate(Vc,1:3) = scoord(1:3)
		Icoordinate(Vc,1:3) = ii(4:6)
		sublattice(Vc) = ii(2)
	end do
	close(1)

	allocate(nnm(maxval(nnb)))
	allocate(bd(maxval(nnb)))

	do nr = 1, nnr
		write(charbond, '(A4,I1,A4)') "bond", nr,".dat"
		open(1,file=trim(charbond),action='read')
		read(1,*) Nb(nr)
		close(1)
	end do
	allocate(bond(nnr,maxval(Nb),2))
	allocate(sitebond(nnr,Vol,maxval(nnb)))
	allocate(bondcount(Vol))

	do nr = 1, nnr
		write(charnnsite, '(A6,I1,A4)') "nnsite", nr,".dat"
		write(charbackdir,'(A7,I1,A4)') "backdir",nr,".dat"
		open(1,file=trim(charnnsite), action='read')
		open(2,file=trim(charbackdir),action='read')

		do i=1, Vol
			read(1,*) Vc, nnm(1:nnb(nr))
			read(2,*) Vc,  bd(1:nnb(nr))

			nnsite(nr,Vc,1:nnb(nr)) = nnm(1:nnb(nr))
			backdir(nr,Vc,1:nnb(nr)) = bd(1:nnb(nr))
		end do
		close(1)
		close(2)

		write(charbond, '(A4,I1,A4)') "bond", nr,".dat"
		open(1,file=trim(charbond),action='read')
		read(1,*) i
		if( i/=Nb(nr)  ) then
			write(*,*) "Err: bond number"
			stop
		end if
		do i=1, Nb(nr)
			read(1,*) bond(nr,i,1:2)
		end do
		close(1)

		!--- arrange nnsite, backdir, sitebond ---!
		nnsite(nr,:,:) = 0
		backdir(nr,:,:) = 0
		bondcount = 0
		do i=1, Nb(nr)
			si = bond(nr,i,1)
			sj = bond(nr,i,2)
			bondcount(si) = bondcount(si)+1
			bondcount(sj) = bondcount(sj)+1

			nnsite(nr,si,bondcount(si)) = sj
			nnsite(nr,sj,bondcount(sj)) = si

			sitebond(nr,si,bondcount(si)) = i
			sitebond(nr,sj,bondcount(sj)) = i

			backdir(nr,si,bondcount(si)) = bondcount(sj)
			backdir(nr,sj,bondcount(sj)) = bondcount(si)

			positivec(nr,si,bondcount(si)) =  1
			positivec(nr,sj,bondcount(sj)) = -1
		end do
	end do

	!--- check positivec ---!
	do i=1, Vol
		do j=1, nnb(1)
			if( positivec(1,i,j)*positivec(1,nnsite(1,i,j),backdir(1,i,j))>0 ) then
				write(*,'(I4, 4I4, 4I4)') i, nnsite(1,i,1:4), backdir(1,i,1:4)
				write(*,'(I4, 4I4, 4I4)') nnsite(1,i,j), nnsite(1,nnsite(1,i,j),1:4), backdir(1,nnsite(1,i,j),1:4)
				write(*,*) "Err: positivec"
				stop
			end if
		end do
	end do

	!open(1,file="positivec.dat",action="write")
	!do i=1, Vol
	!	do j=1, nnb(1)
	!		write(1,'(I4,4I4)') i, positivec(1,i,1:4)
	!	end do
	!end do
	!close(1)

	open(1,file="sublatvec.dat",action="read")
	do i=1, subl
		read(1,*) sublatvec(i,1:3)
	end do
	close(1)

	open(1,file="reclatvec.dat",action="read")
	do i=1, D
		read(1,*) reclatvec(i,1:3)
	end do
	close(1)


#ifdef IFORT
	stat = chdir(trim(current_dir))
#endif

#ifdef GFORT
	call chdir(trim(current_dir), stat)
#endif

	if ( stat /= 0 ) then
		print*, "Err: go to directory "//trim(current_dir)//" failure."
		stop
	end if

	deallocate(nnm)
	deallocate(bd)
END SUBROUTINE read_latt
