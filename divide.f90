! module for cell division
! radii of all cells are evaluated and if they are doubled, 
! original cells is replaced by two new cells of the same type as parent. 
! the centres of the masses of two new cells are alogned along randomly oriented line



subroutine divide(rank)

use params
use functions
use merstwist

implicit none

real :: divr0,divr1,divr 
integer :: ran4(100),nnst,iss 
integer :: i,id, counter,rank
integer :: l,nnn,seed 
real :: angle,alpha 

real :: time

real :: rtot ! the combined radius of the cell (both hard and soft cores)

! initialization

! eligible total radii for division are difference for anterior and posterior cells 
divr0=sqrt(2.0)*tr0
divr1=sqrt(2.0)*tr1

counter=1
nnn=8 ! define 8 directions for a "random" orientation of division line
alpha=2*PI/nnn ! 'unit' angle in radians
call cpu_time(time)
seed=5513974+time*10000
call init_genrand(seed)

do i=1,nnn  
	ran4(i) =  int(grnd()*nnn)+1 ! random 'angle' from 0 to 2*pi with the increment of unit angle
end do

ndiv_cell=0 ! reset the number of dividing cells
nnst=nns

! am before dividing and creating new cells i need to see if i have enough space to store those
do id=1,nnst
	rtot=radius(id)+delta(id) ! the current total radius of the cell
	! it is assumed that both soft and hard cells divide at the same size
	if(cell_type(id).eq.0) then
		if((rtot.ge.divr0).and.(id.ne.ni_node)) counter=counter+1
	else if (cell_type(id).eq.1) then
		if((rtot.ge.divr1).and.(id.ne.ni_node)) counter=counter+1
	end if
end do

! am: not enough space to write new data. extend arrays
if(allo_length.le.(nns+counter)) then
	call reallo(nnst,nns+counter) 
	reallocate=1
end if
!otherwise i am good and can use current arrays
counter=1
do id=1,nnst
	rtot=radius(id)+delta(id)
if (id.ne.ni_node) then
if(((cell_type(id).eq.0).and.(rtot.ge.divr0)).or.&
	&((cell_type(id).eq.1).and.(rtot.ge.divr1))) then ! 
	
		counter=counter+1
		ndiv_cell=ndiv_cell+1
		
		! division

		nns=nns+1
		angle=ran4((mod(id,nnn)+1))*alpha! random angle

		! am: calculate new center's of masses
		
		x(nns)=x(id)+0.5*rtot*cos(angle)
		y(nns)=y(id)+0.5*rtot*sin(angle)
		
		x(id)=x(id)+0.5*rtot*cos(angle+pi)
		y(id)=y(id)+0.5*rtot*sin(angle+pi)

		xm(nns)=xm(id)
		ym(nns)=ym(id)

        cell_type(nns)=cell_type(id)		

		if (cell_type(id).eq.0) then
			radius(id)=r0
			!radius_df(id)=r0
			delta(id)=k0*r0
		else
			radius(id)=r1
			!radius_df(id)=r1
			delta(id)=k1*r1
		end if
		if(cell_type(nns).eq.0) then
			radius(nns)=r0
			!radius_df(nns)=r0
			delta(nns)=k0*r0
		else
			radius(nns)=r1
			!radius_df(nns)=r1
			delta(nns)=k1*r1
		end if
		
	end if
end if
end do

return

end
