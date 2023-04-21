
SUBROUTINE worm
	implicit none

	do
		call move
		if( .not. ZG ) exit
	end do
END SUBROUTINE worm

SUBROUTINE move
	implicit none
	integer(4)                  :: dir
	integer(4)                  :: ib     ! bond label
	integer(4)                  :: curr   ! bond state
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
	curr  = current(ir,ib)

	if( curr==0 ) then
		p = Pbadd(ir)
	else if( curr==1 ) then
		p = 1.d0/Pbadd(ir)
	else
		write(*,*) "Err: current", curr
	end if

	if( p>1.d0 .or. rn()<p ) then
		current(ir,ib) = mod(curr+1,2)
		BondNum(ir) = BondNum(ir) - curr*2+1
		moveIM      = nnsite(ir,moveIM,dir)
	end if
	IraMasha(im) = moveIM
	!====================!

	!===== Judge whether Ira and Masha at same site =====!
	if( IraMasha(1)==IraMasha(2) ) then
		ZG    = .false.
	end if
	!====================================================!
END SUBROUTINE move

