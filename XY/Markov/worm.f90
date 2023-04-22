
SUBROUTINE worm
	implicit none
	integer(4)                  :: ir

	!--- decide the flow of this update ---!
	ir=1
	if( rn()<0.5d0 ) then
		positivec(ir,:,:) = -positivec(ir,:,:)
	end if

	do
		call move
		if( .not. ZG ) exit
	end do
	call check
END SUBROUTINE worm

SUBROUTINE move
	implicit none
	integer(4)                  :: dir
	integer(4)                  :: ib     ! bond label
	integer(4)                  :: curr   ! bond state
	integer(4)                  :: flow
	integer(4)                  :: im, moveIM
	real(8)                     :: p
	integer(4)                  :: ir

	step = step+1.d0

	if ( .not. ZG ) then
		IraMasha   = INT(rn()*Vol)+1; if (IraMasha(1)>Vol) stop "Err: in move, rn error"
		ZG    = .true.
	end if

	ir=1

	!--- Move Ira or Masha ---!
	if( rn()<0.5d0 ) then
		im = 1
	else
		im = 2
	end if
	moveIM = IraMasha(im)


	dir = floor(rn()*nnb(ir))+1
	ib  = sitebond(ir,moveIM,dir)
	curr = current(ir,ib)
	flow = positivec(ir,moveIM,dir)
	!--- Ira and Masha carry different sign ---!
	if( im==2 )  flow=-flow

	if( abs(curr+flow)>MaxCur ) stop "Err: cur>MaxCur"

	p = Pbadd(curr,flow)
	if( p>1.d0 .or. rn()<p ) then
		current(ir,ib) = curr+flow
		moveIM         = nnsite(ir,moveIM,dir)
	end if
	IraMasha(im) = moveIM
	!====================!

	!===== Judge whether Ira and Masha at same site =====!
	if( IraMasha(1)==IraMasha(2) ) then
		ZG    = .false.
	end if
	!====================================================!
END SUBROUTINE move

SUBROUTINE check
	implicit none
	integer(4)                  :: currdiv
	integer(4)                  :: ir
	integer(4)                  :: dir
	integer(4)                  :: i

	ir=1
	do i=1, Vol
		currdiv = 0
		do dir=1, nnb(ir)
			currdiv = currdiv + current(ir,sitebond(ir,i,dir))*positivec(ir,i,dir)
		end do
		if( currdiv/=0 ) then
			write(*,*) "Err: div Jij/=0", i, currdiv
			stop
		end if
	end do
END SUBROUTINE check
