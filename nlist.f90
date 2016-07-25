! AM: Module to creates list of neighbour cells and calculates local densities
!****** FUNCTIONS ******!
! dist() - calculates list of cells that are within cutoff
! dens(nstart,nend) - calculates number of cells within given distance for each cell
!***********************!


MODULE distance

use params
use merstwist

implicit none

contains

! identify neighbor lists for fast integrations. 
! cell-cell interaction cutoffs are set here in rcut
SUBROUTINE dist()

! AM: use modules instead of common blocks

real :: rcut,rr,minr
integer :: i,j,k,l

co=0

minr=100.
rcut=7.0 ! cell-cell interaction cutoff

do i=nstart,nend ! loops over cells assigned to current processor
	do k=1,nns ! loops over total cells
		if (i.eq.k) goto 70 ! am: exclude self-interaction
		rr=sqrt((x(i)-x(k))**2.0+(y(i)-y(k))**2.0)
		if (rr.lt.minr) minr=rr
		if (rr.lt.rcut) then
			co=co+1 ! am: interaction counter for same box mass points
			rd(1,co)=i
			rd(2,co)=k
		end if
70  end do
end do

return
END SUBROUTINE dist

! NDIST calculates number of nearest neighbors. 
! used for cell neighbour alignment and is NOT used in the current version of the paper
SUBROUTINE ndist()

! AM: use modules instead of common blocks

real :: rcut,rr,minr
integer :: i,j,k,l

nco=0

minr=100.

do i=nstart,nend ! loops over cells assigned to current processor
	if (cell_type(i).eq.0) then
		rcut=2.5
	else
		rcut=3.0
	end if
	do k=1,nns ! loops over total cells
		if (i.eq.k) goto 70 ! AM: exclude self-interaction
		rr=sqrt((x(i)-x(k))**2.0+(y(i)-y(k))**2.0)
		if (rr.lt.minr) minr=rr
		if (rr.lt.rcut) then
			nco=nco+1 ! AM: interaction counter for same box mass points
			nrd(1,nco)=i
			nrd(2,nco)=k
			nrd(3,nco)=nrd(3,nco)+1
		end if
70  end do
end do

return
END SUBROUTINE ndist

! DENS calculated cell number density within a certain distance for a given cell.
! the offset due to cell growth (effective decrease of cell number density on the growth regions) is accounted for
SUBROUTINE dens(nstart,nend)

integer :: nstart,nend, i, j,counter
real :: dist_cell, rcut, cell_incr1, cell_incr2

rcut=2.5 ! take only cells within rcut range from a given cell

do i=nstart,nend
	av_rho(i)=0.0	
end do
! calculate number of cells within interaction range

do i=nstart,nend
	counter=0
	dist_cell=0.0
	! the increase of the distance due to the growth of the cell is
	if (cell_type(i).eq.0) then
		cell_incr1=radius(i)+delta(i)-tr0
	else
		cell_incr1=radius(i)+delta(i)-tr1
	end if
	do j=1,nns
		! exclude cell itself and exclude the growth related increase of distance
		if (cell_type(j).eq.0) then
			cell_incr2=radius(j)+delta(j)-tr0
		else
			cell_incr2=radius(j)+delta(j)-tr1
		end if
		if (i.ne.j) dist_cell=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2)-(cell_incr1+cell_incr2)
		if (dist_cell.le.rcut) then
		   if(j.eq.ni_node) then 
		     counter = counter+4
		   else
		     counter = counter+1
		   end if
		end if
	end do
	av_rho(i)=av_rho(i)+counter
end do

return

END SUBROUTINE dens

END MODULE distance
