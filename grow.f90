! module for cells growth
! growth is limited to growth regions and also is suppressed by cell number density
! number density at which cells can growth is different for anterior and posterior cells
! since we want for both cases to have number density of ~0.6
! node contributed to ~4 cells to its neighbours cells (rough estimate based on the are of the node) 

module grow

use params
use functions
use merstwist
 
contains 

!***********************************************************************

subroutine grow_cells(nstart, nend, step,rank)

integer :: ni, nstart, nend, seed, nmindist, rank,step
real :: dist_node, gr
real :: mindist ! smallest distance and number of nearest neighbors
real :: time
real :: rtot
real :: gr_hc_thresh, gr_sc_tresh,dens_hc,dens_sc,gr0,agr0,bgr0,agr1,bgr1
real :: dens_upper_0,dens_lower_0,dens_upper_1,dens_lower_1 ! densities for hard and soft cells, thresholds for growth rate
real :: loc_growr

! number of cells in the vicinity of a cell of given type. number is estimated for a rcut defined in dist module (rcut=2.5)
! number of cells is different for anterior and posterior cells since for the same rcut and threshold number density
! we will have fewer posterior cells. 
 
dens_upper_0=8 ! suppress growth is local density > dens_upper_0 for anterior cell
dens_lower_0=0 

dens_upper_1=5 ! suppress growth is local density > dens_upper_1 for posterior cell
dens_lower_1=0

gr0=0.00005 ! growth rate per simulation step

call cpu_time(time)
! generate random fluctuations to mimic asynchronized growth
seed=5513974+time*(rank+1)*10000
call init_genrand(seed)

! AM: Reset number of growing cells
if(mod(step,1000).eq.1) growing_cells(rank+1)=0

!am: increase pressure  
do ni=nstart,nend ! loop over cells

	dist_node=sqrt((x(ni)-x(ni_node))**2+(y(ni)-y(ni_node))**2)
	rtot=radius(ni)+delta(ni)
	if(x(ni).le.x(ni_node)) then
    loc_growr=growra
    else
    loc_growr=growrp
    end if
    if(dist_node.le.loc_growr) then
	  if (cell_type(ni).eq.0) then ! hard cells
		if(av_rho(ni).lt.dens_lower_0) then ! border cells and cells close to node
			gr=gr0*2*grnd()
			radius(ni)=radius(ni)+gr
		else if (av_rho(ni).gt.dens_upper_0) then ! really pressed cells
			gr=0.0
			radius(ni)=radius(ni)+gr
		else ! linear dependency of gr on concentration
        	gr=gr0*2*grnd()
			radius(ni)=radius(ni)+gr
		end if
	  else ! soft cells
	  	if(av_rho(ni).lt.dens_lower_1) then ! border cells and cells close to node
			gr=gr0*2*grnd()
			radius(ni)=radius(ni)+gr
		else if (av_rho(ni).gt.dens_upper_1) then ! really pressed cells
			gr=0.0
			radius(ni)=radius(ni)+gr
		else ! linear dependency of gr on concentration
        	gr=gr0*2*grnd()
			radius(ni)=radius(ni)+gr
		end if
	  end if
		!write(*,*) "growing ni", dist_node, av_rho(ni)
	else ! this part shrinks cells that are out of the growth zone
		if (cell_type(ni).eq.0) then
			if(rtot.gt.tr0) then
				gr=0.0005
				radius(ni)=radius(ni)-gr
			else
				radius(ni)=r0
			end if
		else
			if(rtot.gt.tr1) then
				gr=0.0005
				radius(ni)=radius(ni)-gr
			else
				radius(ni)=r1
			end if
		end if
	end if
	call assign_soft_radii(nstart,nend)
	
	radius(ni_node)=node_r
	delta(ni_node)=knode*node_r
                 
end do ! end loop over cells

return

end subroutine grow_cells


end module

