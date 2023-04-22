PROGRAM analyse_data
	implicit none

	character(20)         :: nfile                 ! the data file name
	character(10)         :: charflg
	integer               :: flg
	integer               :: D
	integer               :: Lx, Ly, Lz
	integer               :: Vol, gVol
	integer               :: nblock                ! the number of blocks
	integer               :: rb
	real(8),  allocatable :: gr(:,:)
	real(8)               :: ave, cor, erb         ! the average of data
	real(8)               :: devn, devp, dev
	real(8)               :: nda                   ! current data
	integer               :: stat
	integer               :: i, j, k

	call getarg(1,nfile)
	call getarg(2,charflg)
	read(charflg,*) flg
	if( len_trim(nfile)==0 ) then
		print *, "error: no data file name"
		stop
	end if

	nblock   = 0
	open(1, file=trim(nfile), action='read')
	read(1,*) D, Lx, Ly, Lz, nblock
	Vol = Lx*Ly*Lz
	if( flg==1 ) then
		gVol=(Lx/2+1)*(Ly/2+1)*(Lz/2+1)-1
	else
		gVol=Lx/2
	end if

	if( nblock>=1024 ) then
		if( nblock/1024<4 ) then
			rb = nblock/512
			nblock = 512
		else
			rb = nblock/1024
			nblock = 1024
		end if
	else
		rb = 1
	end if
	allocate(gr(0:gVol, nblock))
	gr = 0.d0

	open(1, file=trim(nfile), action='read')
	do i=0, nblock*rb-1
		do j=0, gVol
			read(1,*) k, nda
			gr(k,i/rb+1) = gr(k,i/rb+1)+nda
		end do
	end do
	close(1)
	gr = gr/rb


	do k=0, gVol
		ave  = sum(gr(k,:))/nblock
		devn = 0.d0; devp = 0.d0; dev = 0.d0;  cor = 0.d0
		do i=1, nblock
			devn = gr(k,i) - ave
			dev  = dev+ devn*devn
			cor  = cor+ devn*devp
			devp = devn
		end do
		cor = cor/dev
		erb = dsqrt(dev/nblock/(nblock-1))

		write(*,'(I8, 3ES20.8)') k, ave, erb, cor
	end do

	stop
END PROGRAM analyse_data
