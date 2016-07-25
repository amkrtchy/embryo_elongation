! module has various useful function
! some functions are used to set simulation parameters (like cell types etc)
! all functions have short description preceding



module functions

use params
use merstwist
use sorting

implicit none

contains

!***********************************************************************
! distribute cells between processors
! estimates the ids of cells that will be assigned to given processor
!
! nproc - array of processes
! rank  - rank of a process calling this function
!***********************************************************************
subroutine distrib_proc(nproc,rank)

integer :: nproc(0:(max_proc-1))
integer :: i, rank

do i=0,nt-1
	if(i.ne.(nt-1)) then 
		nproc(i) = nns/nt
	else 
		nproc(i)=(nns/nt)+mod(nns,nt)
	end if
end do

nstart=(nns/nt)*rank+1 ! start cell id assigned to given processor
nend=(nns/nt)*rank+nproc(rank) ! end cell id assigned to given processor

end subroutine distrib_proc
!***********************************************************************
! estimate which chunk of arrays assigned to given processor.
! used to communicated cell data via MPI between different processors
!
! nproc - array of processes
! pair-proc - data is send/received between pair of processes to avoid deadlocks
!***********************************************************************
subroutine start_end_index(nproc,pair_proc)

integer :: i,k,j,temp_proc
integer :: nproc(0:(max_proc-1)),pair_proc(iter,0:(max_proc-1))
do i=1,iter
	if(i.eq.1) then ! first iteration
		do k=0,nt-1
			start_index_send(i,k)=(nns/nt)*k+1
			end_index_send(i,k)=(nns/nt)*k+nproc(k)
			
			start_index_receive(i,k)=(nns/nt)*pair_proc(i,k)+1
			end_index_receive(i,k)=(nns/nt)*pair_proc(i,k)+nproc(pair_proc(i,k))
		end do
	else ! consequent iterations are calculated based on first one
		do k=0,nt-1 ! loop to calc. send indices
			if(start_index_send(i-1,k).gt.start_index_receive(i-1,k)) then
				start_index_send(i,k)=start_index_receive(i-1,k)
			else
				start_index_send(i,k)=start_index_send(i-1,k)
			end if
			if(end_index_send(i-1,k).lt.end_index_receive(i-1,k)) then
		    	end_index_send(i,k)=end_index_receive(i-1,k)
		    else
		    	end_index_send(i,k)=end_index_send(i-1,k)
		    end if
		
			! for receive indices i need to know which way data goes
			if(pair_proc(i,k).gt.k) then ! data flow is forward
				start_index_receive(i,k)=end_index_send(i,k)+1
				! am on 21/08/14: find out what was the initial process holding that cell
				do j=0,nt-1
					if (start_index_send(1,j).eq.start_index_receive(i,k)) temp_proc=j
				end do
				
				end_index_receive(i,k)=start_index_receive(i,k)-1
				do j=0,2**(i-1)-1
					end_index_receive(i,k)=end_index_receive(i,k)+nproc(temp_proc+j)
				end do

			else ! data flow is backwards
				end_index_receive(i,k)=start_index_send(i,k)-1
				start_index_receive(i,k)=end_index_receive(i,k)+1				
				do j=0,nt-1
					if (end_index_send(1,j).eq.end_index_receive(i,k)) temp_proc=j
				end do
				do j=0,2**(i-1)-1
					start_index_receive(i,k)=start_index_receive(i,k)-nproc(temp_proc-j)
				end do
			end if
		end do ! end loop to calculate indices
		
	end if ! end condition to calc indices for iter>1
	
end do

end subroutine start_end_index

!***********************************************************************
! Read coordinates and radii of cells from configuration file writecont
! Initialize data arrays and assign cell type
!***********************************************************************
SUBROUTINE init_data()

integer :: io, data_lines
integer :: i

data_lines=0

open(unit=75,file='writecont',status='unknown')
do
	read( 75, *, iostat = io)
	if (io < 0) exit
	data_lines=data_lines+1
end do	
close(75)

nns=data_lines ! AM: number of lines in file is the same as number of cells

! AM: allocate data arrays based on the number of cells and fill in the arrays
call allo(nns) 

open(unit=75,file='writecont',status='unknown')
do i=1, nns
	read(75,*) x(i), y(i), radius(i)
end do 
close(75)

! after reading data from input file, assign cells their radii
call assign_cell_type(1,nns)
radius(ni_node)=node_r

call assign_soft_radii(1,nns)

return

END SUBROUTINE init_data

!***********************************************************************
! Allocate array of growing cells. 
!***********************************************************************

subroutine allo_gr_cells(nt)

integer :: nt,i

allocate(growing_cells(nt))
do i=1,nt
	growing_cells(i)=0
end do

end subroutine allo_gr_cells

!***********************************************************************
! assign cell types based on their relative position from the node
! This is the function that defines probability of a cell type
! 
! (nstart, nend) isolates cells for which processor calling the function 
! will attempt to set types
!***********************************************************************

subroutine assign_cell_type(nstart,nend)

integer :: ni, nstart,nend, ntranz, seed,flag,rank ! counter
real :: dist_node ! distance of given cell from the node
real :: a, b, k, c ! coefficients to characterize the weight of distribution function
real :: time
integer :: weight

! parameters of a "linear" probability of a cell being A/P type based on the distance from the node
a=x(ni_node)-softr
b=x(ni_node)
k=100/(b-a)
c=-k*a
	
	do ni=nstart, nend
	
	cell_type(ni)=0 ! anterior. set cell type to 1 if posterior
	
	! uncomment block of code below for wild type embryo here
!	if(x(ni).le.a) cell_type(ni)=0
!	if(x(ni).ge.b) cell_type(ni)=1
!	if((x(ni).lt.b).and.(x(ni).gt.a)) then 
!		weight=k*x(ni)+c
!		cell_type(ni)=wrand(weight)
!	end if
	cell_type(ni_node)=1 ! node is a cell of a posterior type

	end do
	return
end subroutine assign_cell_type


function wrand(weight) result(fn_val)

integer, intent(in) :: weight
integer ::i,seed ! counter
integer, dimension(100) :: weight_array
real :: time,fn_val 

do i=1,weight
weight_array(i)=1
end do
do i=(weight+1),100
weight_array(i)=0
end do

fn_val=weight_array(int(grnd()*100)+1)
return
end function wrand

!***********************************************************************
!***********************************************************************
! Hard core radii for all cells are the same
! Soft core radii depend on the type of the cell. radii are assigned by code below
! 
! (nstart, nend) isolates cells for which processor calling the function 
! will attempt to set types
!***********************************************************************

subroutine assign_soft_radii(nstart,nend)

integer :: ni, nstart, nend, rank ! counter
real :: dist_node ! distance of given cell from the node
real :: a, b, c, d ! coefficients to characterize the changes in radii

	! assign soft radii to the cells based on their type.
	! type 0 cells have smaller soft radii
	! type 1 cells have larger cell radii
	
	do ni=nstart, nend
		if(cell_type(ni).eq.0) then 
			delta(ni)=k0*radius(ni)
		else 
			delta(ni)=k1*radius(ni)
		end if
		delta(ni_node)=knode*node_r ! node always stays the same
	end do

	return
end subroutine assign_soft_radii

!***********************************************************************
! calculate diffusion coefficient of cells based on their type and the number of neighbor cells
!
! dmax, dmin: diffusions of P/A cells in isolation
! ni: id of the cell
!***********************************************************************

function diffc(dmax, dmin, ni) result(fn_val)
	real, intent(in) :: dmin,dmax
	integer, intent(in) :: ni
	real :: rsr,r_attr
	real :: fn_val
	integer :: rank, is,i
	integer :: bonds
	
	! if cells is within attractive interaction region with 3 other similar cell
	! diffusion coefficient slows.
	! first calculate nn of similar cells 
	bonds=0
	do i=1,co
		if (ni.eq.rd(1,i)) then
			if((cell_type(rd(1,i)).eq.0).and.(cell_type(rd(2,i)).eq.0)) then
				r_attr=radius(rd(1,i))+radius(rd(2,i))+delta(rd(1,i))+delta(rd(2,i))
			end if
			rsr=sqrt((x(rd(1,i))-x(rd(2,i)))**2+(y(rd(1,i))-y(rd(2,i)))**2) 
			if (rsr.eq.0.0) rsr=1.0e-08	
			if(rsr.le.r_attr) bonds=bonds+1 
		end if
	end do
	! calculate distances of node from simulation box edges
	if((cell_type(ni).eq.0).and.(bonds.ge.3)) then
		fn_val = dmax/5.0
    else if((cell_type(ni).eq.0).and.(bonds.ge.6)) then
        fn_val = dmin
	else 
		fn_val = dmax
	end if
	return
end function diffc

! allocate arrays for various cell data
SUBROUTINE allo(length)

    integer :: length,i,j,ni

	allo_length = 2*length
	
	allocate(x(allo_length),y(allo_length))
	allocate(vx(allo_length),vy(allo_length))
    allocate(radius(allo_length))
    allocate(cell_type(allo_length))
    allocate(xm(allo_length),ym(allo_length))
    allocate(xp(allo_length),yp(allo_length))
    allocate(c(allo_length),d(allo_length))
    allocate(ctmp(allo_length),dtmp(allo_length))
    
    allocate(fx(allo_length),fy(allo_length))
    
	allocate(rd(2,allo_length*allo_length))
	allocate(nrd(3,allo_length*allo_length))
	
	allocate(delta(allo_length))
	
	allocate(kin(allo_length),pot(allo_length))

	allocate(av_rho(allo_length))
	
	allocate(fluct(2*allo_length))
	
	allocate(av_vx(2*allo_length), av_vy(2*allo_length))

	do ni=1,allo_length ! loop over all cells
		x(ni)=0.0
		y(ni)=0.0
		vx(ni)=0.0
		vy(ni)=0.0
		xm(ni)=0.0
		ym(ni)=0.0
		xp(ni)=0.0
		yp(ni)=0.0				
		fx(ni)=0.0
		fy(ni)=0.0
		kin(ni)=0.0
		pot(ni)=0.0
		radius(ni)=0.0
		delta(ni)=0.0
		c(ni)=0.0
		d(ni)=0.0
		ctmp(ni)=0.0
		dtmp(ni)=0.0
		cell_type(ni)=0
		av_rho(ni)=0.0
		av_vx(ni)=0.0
		av_vy(ni)=0.0
		
	end do ! end loop over all cells
	do ni=1,2*allo_length
		fluct(ni)=0.0
	end do
	do i=1,2
		do ni=1,allo_length*allo_length
			rd(i,ni)=0
		end do
	end do
	do ni=1,allo_length*allo_length
		do i=1,2
			rd(i,ni)=0
			nrd(i,ni)=0
		end do
		nrd(3,ni)=0
	end do
  return

END SUBROUTINE allo

!***********************************************************************
! reallocate cell arrays to increase array sizes dues to increased number of cell in system
! new arrays are allocates, old array data is copied into new arrays and old arrays are dealocated
! 
! oldnns - old array length
! newnns - new array length
!***********************************************************************

SUBROUTINE reallo(oldnns,newnns)
	
	! am on 15/08/14: first create temporary arrays and copy previous arrays into them
	integer :: newnns,oldnns
	real :: xtemp(oldnns),ytemp(oldnns)
	real :: vxtemp(oldnns),vytemp(oldnns)
	real :: radiustemp(oldnns)
	real :: xmtemp(oldnns),ymtemp(oldnns)
	real :: xptemp(oldnns),yptemp(oldnns)
	real :: fxtemp(oldnns),fytemp(oldnns)
	real :: kintemp(oldnns),pottemp(oldnns)

	integer :: cell_type_temp(oldnns)
	integer :: i, ni,j
	integer :: rdtemp(2,oldnns*oldnns)
	integer :: nrdtemp(3,oldnns*oldnns)
	real :: deltatemp(oldnns)
	!real :: radius_df_temp(oldnns)
	real :: av_rho_temp(oldnns)

	! copy old values to the temporary arrays
	do ni=1,oldnns ! loop over all cells
		xtemp(ni)=x(ni)
		ytemp(ni)=y(ni)
		vxtemp(ni)=vx(ni)
		vytemp(ni)=vy(ni)
		xmtemp(ni)=xm(ni)
		ymtemp(ni)=ym(ni)
		xptemp(ni)=xp(ni)
		yptemp(ni)=yp(ni)
		fxtemp(ni)=fx(ni)
		fytemp(ni)=fy(ni)
		kintemp(ni)=kin(ni)
		pottemp(ni)=pot(ni)
		radiustemp(ni)=radius(ni)
		deltatemp(ni)=delta(ni)
		cell_type_temp(ni)=cell_type(ni)
		
		!radius_df_temp(ni)=radius_df(ni)
		
		av_rho_temp(ni)=av_rho(ni)
	end do ! end loop over all cells
	
	do ni=1,oldnns*oldnns
		do i=1,2
			rdtemp(i,ni)=rd(i,ni)
			nrdtemp(i,ni)=nrd(i,ni)
		end do
		nrdtemp(3,ni)=nrd(3,ni)
	end do
	
	! deallocate old arrays
	deallocate(x, y, xm, ym, xp, yp,vx,vy)  
	deallocate(c,d,ctmp,dtmp)
    deallocate(radius,delta)  
    deallocate(fx, fy)  
    deallocate(kin, pot)  
    deallocate(rd,nrd)
    deallocate(cell_type)
    deallocate(av_rho,av_vx,av_vy)
    deallocate(fluct)
    
    
	! am on 15/08/14: allocate new arrays    
    allo_length = 2*newnns
	
	allocate(x(allo_length),y(allo_length))
	allocate(vx(allo_length),vy(allo_length))
	allocate(cell_type(allo_length))
    allocate(radius(allo_length))
    allocate(xm(allo_length),ym(allo_length))
    allocate(xp(allo_length),yp(allo_length))
    allocate(c(allo_length),d(allo_length))
    allocate(ctmp(allo_length),dtmp(allo_length))

    allocate(fx(allo_length),fy(allo_length))
    
    allocate(kin(allo_length),pot(allo_length))
    
    allocate(rd(2,allo_length*allo_length))
    allocate(nrd(3,allo_length*allo_length))
    
    allocate(delta(allo_length))
    allocate(av_rho(allo_length))
    allocate(av_vx(2*allo_length),av_vy(2*allo_length))
    allocate(fluct(2*allo_length))
    
    ! copy back the previous values into new arrays
    do ni=1,oldnns ! loop over all cells
		x(ni)=xtemp(ni)
		y(ni)=ytemp(ni)
		vx(ni)=xtemp(ni)
		vy(ni)=ytemp(ni)
		xm(ni)=xmtemp(ni)
		ym(ni)=ymtemp(ni)
		xp(ni)=xptemp(ni)
		yp(ni)=yptemp(ni)
		fx(ni)=fxtemp(ni)
		fy(ni)=fytemp(ni)
		kin(ni)=kintemp(ni)
		pot(ni)=pottemp(ni)
		radius(ni)=radiustemp(ni)
		delta(ni)=deltatemp(ni)
		c(ni)=0.0
		d(ni)=0.0
		ctmp(ni)=0.0
		dtmp(ni)=0.0
		cell_type(ni)=cell_type_temp(ni)
		av_rho(ni)=av_rho_temp(ni)
		av_vx(ni)=0.0
		av_vy(ni)=0.0
		
	end do ! end loop over all cells
	do ni=1,2*allo_length
		fluct(ni)=0.0
	end do
	do ni=1,oldnns*oldnns
		do i=1,2
			rd(i,ni)=rdtemp(i,ni)
			nrd(i,ni)=nrdtemp(i,ni)
		end do
		nrd(3,ni)=nrdtemp(3,ni)
	end do
	
	do ni=oldnns+1,allo_length ! loop over all cells
		x(ni)=0.0
		y(ni)=0.0
		vx(ni)=0.0
		vy(ni)=0.0
		xm(ni)=0.0
		ym(ni)=0.0
		xp(ni)=0.0
		yp(ni)=0.0
		fx(ni)=0.0
		fy(ni)=0.0
		kin(ni)=0.0
		pot(ni)=0.0
		radius(ni)=0.0
		delta(ni)=0.0
		c(ni)=0.0
		d(ni)=0.0
		ctmp(ni)=0.0
		dtmp(ni)=0.0
		cell_type(ni)=0
		av_rho(ni)=0.0
		av_vx(ni)=0.0
		av_vy(ni)=0.0
	end do ! end loop over all cells

	do ni=(oldnns*oldnns+1),allo_length*allo_length
		do i=1,2
			rd(i,ni)=0
			nrd(i,ni)=0
		end do
		nrd(3,ni)=0
	end do
    	
  return
end subroutine reallo

!***********************************************************************
! initialize buffer arrays to send/receive cell data in MPI
! multidimensional arrays are converted into 1D arrays before data communication
!***********************************************************************
subroutine ini_buffer(startsi,endsi,startri,endri)
integer :: startsi, endsi,startri,endri,l,k

! AM: allocate buffer arrays for send/receive
		allocate(xts(endsi-startsi+1))
		allocate(xmts(endsi-startsi+1))
		allocate(yts(endsi-startsi+1))
		allocate(ymts(endsi-startsi+1))
		allocate(radts(endsi-startsi+1))
		allocate(delts(endsi-startsi+1))
		allocate(vxts(endsi-startsi+1))
		allocate(vyts(endsi-startsi+1))
		allocate(cellts(endsi-startsi+1))
		
		allocate(xtr(endri-startri+1))
		allocate(xmtr(endri-startri+1))
		allocate(ytr(endri-startri+1))
		allocate(ymtr(endri-startri+1))
		allocate(radtr(endri-startri+1))
		allocate(deltr(endri-startri+1))
		allocate(vxtr(endri-startri+1))
		allocate(vytr(endri-startri+1))
		allocate(celltr(endri-startri+1))
		
				
		! assign send array data
		l=0
		do k=startsi,endsi
			l=l+1
			xts(l)=x(k)
			xmts(l)=xm(k)				
			yts(l)=y(k)
			ymts(l)=ym(k)
			radts(l)=radius(k)
			delts(l)=delta(k)	
			vxts(l)=vx(k)
			vyts(l)=vy(k)
			cellts(l)=cell_type(k)
		end do
end subroutine ini_buffer

! after communicating data between different processors, restore data arrays into their initial form
subroutine convert_buffer(startri,endri)
integer :: l, k,startri,endri

l=0
do k=startri,endri
	l=l+1
	x(k)=xtr(l)
	xm(k)=xmtr(l)
	y(k)=ytr(l)
	ym(k)=ymtr(l)
	radius(k)=radtr(l)
	!radius_df(k)=raddftr(l)
	delta(k)=deltr(l)
	vx(k)=vxtr(l)
	vy(k)=vytr(l)
	cell_type(k)=celltr(l)
end do

end subroutine convert_buffer

subroutine del_buffer()
	
	deallocate(xts)
	deallocate(xmts)
	deallocate(yts)
	deallocate(ymts)
	deallocate(radts)
	deallocate(delts)
	deallocate(vxts)
	deallocate(vyts)
	deallocate(cellts)
	deallocate(xtr)
	deallocate(xmtr)
	deallocate(ytr)
	deallocate(ymtr)
	deallocate(radtr)
	deallocate(deltr)
	deallocate(vxtr)
	deallocate(vytr)
	deallocate(celltr)
	
end subroutine del_buffer

! calculate forces between interacting cells
! 
! xx, yy - arrays of cell positions
! nstart, nend - start/end ids of cell arrays assigned to given processor
! 
subroutine calc_f(xx,yy,nstart,nend,tms)

integer :: ni,i,nstart,nend,curri
real :: ll1,ll2,ll3,rsr
real :: xx(allo_length),yy(allo_length)
real :: a_hc, b_hc, a_sc, b_sc, a_attr, b_attr
real :: c_hc, c_sc, c_attr
integer :: tms,cnt
real :: time_step

time_step=0.001

do ni=nstart,nend ! loop over mass points	    
	fx(ni)=0.0
	fy(ni)=0.0
	
	pot(ni)=0.0
	av_vx(ni)=0.0
	av_vy(ni)=0.0
	
	!AM 05/06/2014 activate rigid walls
    if (yy(ni).lt.delly) then
		fy(ni)=fy(ni)+(delly-yy(ni))*0.5
	end if

	if (yy(ni).gt.delry) then
		fy(ni)=fy(ni)+(delry-yy(ni))*0.5
	end if

	if (xx(ni).gt.delrx) then
		fx(ni)=fx(ni)+(delrx-xx(ni))*0.05
	end if

	if (xx(ni).lt.dellx) then
		fx(ni)=fx(ni)+(dellx-xx(ni))*0.05
	end if

!***********************************************************************
! AM: intercellular force calculation
        
	do i=1,co
		! AM: if interacting particles are within the same box cells
		! AM: calculate interaction cutoffs
		! the equivalent force constants are: k_hc=strhc, k_sc=strsc, k_attr=stratr

		if (ni.eq.rd(1,i)) then 

        if((cell_type(rd(1,i)).eq.0).and.(cell_type(rd(2,i)).eq.0)) then
			ll1=radius(rd(1,i))+radius(rd(2,i))+delta(rd(1,i))+delta(rd(2,i)) 
			ll2=radius(rd(1,i))+radius(rd(2,i))+delta(rd(1,i))/2.0+delta(rd(2,i))/2.0 
		else
			ll1=0.0
			ll2=radius(rd(1,i))+radius(rd(2,i))+delta(rd(1,i))+delta(rd(2,i)) 
		end if
		
		ll3=radius(rd(1,i))+radius(rd(2,i)) 
		
		! define interaction coefficients
		! all forces have linear form ax+b

		a_hc=-(strhc*ll3-strsc*ll2)/ll3
		b_hc=strhc*ll3
		
		a_sc=-(strsc*ll2)/(ll2-ll3)
		b_sc=(strsc*(ll2**2))/(ll2-ll3)
		
		a_attr=stratr*ll1/(ll1-ll2)
		b_attr=-(stratr*(ll1**2))/(ll1-ll2)

		
		! add the constants for potential energy (from theoretical calculations)
		c_attr=0.5*a_attr*(ll1**2)+b_attr*ll1
		c_sc=0.5*a_sc*(ll2**2)+b_sc*ll2-0.5*a_attr*(ll2**2)-b_attr*ll2+c_attr
		c_hc=0.5*a_hc*(ll3**2)+b_hc*ll3-0.5*a_sc*(ll3**2)-b_sc*ll3+c_sc
		
	
		rsr=sqrt((xx(rd(1,i))-xx(rd(2,i)))**2+&
			& (yy(rd(1,i))-yy(rd(2,i)))**2) 
		if (rsr.eq.0.0) rsr=1.0e-08
		
		if(rsr.lt.ll3) then 
			! hard core repulsion
			fx(rd(1,i))=fx(rd(1,i))+(a_hc*rsr+b_hc)*(xx(rd(1,i))-xx(rd(2,i)))/rsr
			fy(rd(1,i))=fy(rd(1,i))+(a_hc*rsr+b_hc)*(yy(rd(1,i))-yy(rd(2,i)))/rsr
			pot(rd(1,i))=pot(rd(1,i))-0.5*a_hc*(rsr**2)-b_hc*rsr+c_hc
			
		else if ((rsr.ge.ll3).and.(rsr.lt.ll2)) then 
			! elastic soft core repulsion
			fx(rd(1,i))=fx(rd(1,i))+(a_sc*rsr+b_sc)*(xx(rd(1,i))-xx(rd(2,i)))/rsr
			fy(rd(1,i))=fy(rd(1,i))+(a_sc*rsr+b_sc)*(yy(rd(1,i))-yy(rd(2,i)))/rsr
			pot(rd(1,i))=pot(rd(1,i))-0.5*a_sc*(rsr**2)-b_sc*rsr+c_sc
			
		else if ((rsr.ge.ll2).and.(rsr.lt.ll1)) then 
			! attractive interaction
			fx(rd(1,i))=fx(rd(1,i))+(a_attr*rsr+b_attr)*(xx(rd(1,i))-xx(rd(2,i)))/rsr
			fy(rd(1,i))=fy(rd(1,i))+(a_attr*rsr+b_attr)*(yy(rd(1,i))-yy(rd(2,i)))/rsr
			pot(rd(1,i))=pot(rd(1,i))-0.5*a_attr*(rsr**2)-b_attr*rsr+c_sc
		end if

		end if 
	end do ! finish looping over co
	
 !AM: no neighbour alignment after first 100 steps of division time
 !Since velocities might have big jumps because of the overlap after division
 if (n_align.eq.1) then
	if ((tms.gt.1000).and.(mod(tms,1000).gt.100)) then
	av_vx=0.0
	av_vy=0.0
	cnt=0
	do i=1,nco
	 if ((ni.eq.nrd(1,i)).and.(ni.ne.ni_node)) then 
    	    av_vx(ni)=av_vx(ni)+(xx(nrd(2,i))-xm(nrd(2,i)))/time_step
    	    av_vy(ni)=av_vy(ni)+(yy(nrd(2,i))-ym(nrd(2,i)))/time_step
			cnt=cnt+1
	 end if
	end do
	end if
	if (ni.ne.ni_node) then
    	fx(ni)=fx(ni)+beta*(av_vx(ni)+(xx(ni)-xm(ni))/time_step)/(cnt+1)
    	fy(ni)=fy(ni)+beta*(av_vy(ni)+(yy(ni)-ym(ni))/time_step)/(cnt+1)
    end if
 end if
!***********************************************************************
! AM 08/11/12: attepmt to stabilize system
!	if (abs(fx(ni)).ge.maxfx) then
!		if (fx(ni).gt.0.0) then 
!			fx(ni)=maxfx
!		else 
!			fx(ni)=-maxfx
!		end if
!	end if
!	if (abs(fy(ni)).ge.maxfy) then
!		if (fy(ni).gt.0.0) then 
!			fy(ni)=maxfy
!		else 
!			fy(ni)=-maxfy
!		end if
!	end if
end do ! end looping over ni

end subroutine calc_f

!***********************************************************************

end module functions
