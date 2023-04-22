
!============ markov chain =========================!
!! THIS IS PROJECT-DEPENDENT
SUBROUTINE markov(nw)
	implicit none

	integer :: nw
	integer :: istep

	do istep= 1, nw
		call worm()
	enddo

	return
END SUBROUTINE markov

!============== worm  ===========================!
!! THIS IS PROJECT-DEPENDENT
SUBROUTINE worm
	implicit none

	if( .not. ZG ) then
		Ira   = int(rn()*Vol)+1
		Masha = Ira
		ZG    = .True.
		sCur  = (int(rn()*2.d0)*2-1)*sCur
		call move()
		do while ( ZG )
			call move()
		end do
	end if

	return
END SUBROUTINE worm

SUBROUTINE move
	implicit none

	integer       :: Vb
	integer       :: dir
	integer       :: i, j
	integer       :: icur
	real(8)       :: p

	step = step+1.d0
	dir  = int(rn()*nnb)+1
	Vb   = Site(Ira)%nb(dir)
	icur = sCur(dir)

	p = BesselRat(abs(Bond(Vb)),abs(Bond(Vb)+icur))
	if( p>=1.d0 ) then
		call update(dir, Vb, icur)
	else if ( rn()<p ) then
		call update(dir, Vb, icur)
	end if
	call up_gr(dir)

	if( Ira==Masha ) then
		ZG = .False.
		Ira=0
		Masha=0
		if( sum(abs(Rim))/=0 ) stop "Err: Rim"
	end if

	return
END SUBROUTINE move

SUBROUTINE update(udir, uVb, uicur)
	implicit none
	integer      :: udir, uVb, uicur
	integer      :: r

	Bond(uVb) = Bond(uVb)+uicur
	Ira = Site(Ira)%ns(udir)
	Rim(Ms(udir,1)) = Rim(Ms(udir,1))+Ms(udir,2)
	if( abs(Bond(uVb))>MaxCur ) stop ">MaxCur"

	return
END SUBROUTINE update

SUBROUTINE up_gr(dir)
	implicit none
	integer      :: dir
	integer      :: r
	real(8)      :: inc

	inc = 1.d0
	r   = 0.d0

	if( abs(Rim(Ms(dir,1)))>=Lxyz(Ms(dir,1)) ) then
		Rim(Ms(dir,1)) = mod(Rim(Ms(dir,1)), Lxyz(Ms(dir,1)))
	end if

	!-- Calculate the distance between Ira and Masha
	if ( D==1 ) then
		!if ( mod(abs(Rim(1)),Lxyz(1)/2)==0 ) inc = inc*2.d0

		if ( abs(Rim(1))>Lxyz(1)/2 ) then
			r = Lxyz(1)-abs(Rim(1))
		else
			r = abs(Rim(1))
		end if
	else if ( D==2 ) then
		if ( mod(abs(Rim(1)),Lxyz(1)/2)==0 ) inc = inc*2.d0
		!if ( abs(Rim(1))==0 ) inc = inc*2.d0

		if ( abs(Rim(1))>Lxyz(1)/2 ) then
			r = Lxyz(1)-abs(Rim(1))
		else
			r = abs(Rim(1))
		end if

		if ( mod(abs(Rim(2)),Lxyz(2)/2)==0 ) inc = inc*2.d0
		!if ( abs(Rim(2))==0 ) inc = inc*2.d0

		if ( abs(Rim(2))>Lxyz(2)/2 ) then
			r = r+(Lxyz(2)-abs(Rim(2)))*(Lxyz(1)/2+1)
		else if ( abs(Rim(2))>0 ) then
			r = r+abs(Rim(2))*(Lxyz(1)/2+1)
		end if
	else
		if ( mod(abs(Rim(1)),Lxyz(1)/2)==0 ) inc = inc*2.d0
		!if ( abs(Rim(1))==0 ) inc = inc*2.d0

		if ( abs(Rim(1))>Lxyz(1)/2 ) then
			r = Lxyz(1)-abs(Rim(1))
		else
			r = abs(Rim(1))
		end if

		if ( mod(abs(Rim(2)),Lxyz(2)/2)==0 ) inc = inc*2.d0
		!if ( abs(Rim(2))==0 ) inc = inc*2.d0

		if ( abs(Rim(2))>Lxyz(2)/2 ) then
			r = r+(Lxyz(2)-abs(Rim(2)))*(Lxyz(1)/2+1)
		else if ( abs(Rim(2))>0 ) then
			r = r+abs(Rim(2))*(Lxyz(1)/2+1)
		end if

		if ( mod(abs(Rim(3)),Lxyz(3)/2)==0 ) inc = inc*2.d0
		!if ( abs(Rim(3))==0 ) inc = inc*2.d0

		if ( abs(Rim(3))>Lxyz(3)/2 ) then
			r = r+(Lxyz(3)-abs(Rim(3)))*(Lxyz(1)/2+1)*(Lxyz(2)/2+1)
		else if ( abs(Rim(3))>0 ) then
			r = r+abs(Rim(3))*(Lxyz(1)/2+1)*(Lxyz(2)/2+1)
		end if
	end if

	gri(r) = gri(r)+inc

	return
END SUBROUTINE up_gr
