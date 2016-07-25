! main module: implements MPI communication, integrates equation of motions of cells
! and reads/writes inputs and outputs

program embryo_model

use params ! simulation parameters
use merstwist ! random number generator
use functions	! Module for subroutines 
use grow	! Cell growth subroutine
use distance	! Neighbour list generator + density calculations
use sorting	! Sort arrays

implicit none
include 'mpif.h'

!***********************************************************************
! List of Variables
!***********************************************************************

real :: dt
real :: dst
real :: lx,ly
real :: maxx,minx,maxy,miny

! Random number generator vars
real, dimension (:), allocatable :: ran1
real :: x1, x2, rr
integer :: seed,seed2
integer :: ncc(1000),cc,noa

! Counters
integer :: ni,i,j,k,is,div,imax,jmax
integer :: oldnns
integer :: count_node, count_border, count_dens

! Energy
real :: t1, t2
real :: totpot, totkin, toten

! Node statistics
real :: node_data(1000,2)
real :: av_nodex,std_nodex, av_nodey,std_nodey

! Border statistics
real :: temp_x(100000), temp_y(100000)
real :: xmin_border, xmax_border,ymin_border,ymax_border

! Density statistics
real :: delta_x,delta_y
integer, dimension(:,:), allocatable :: density

! Statistics for number of cells
integer :: na, np, ngr

! MPI vars
integer :: rank,rc,ierr
integer :: source,dest,tag,stat(mpi_status_size)
integer :: nnsproc(0:(max_proc-1)),pair_proc(iter,0:(max_proc-1))

character(len=1024) :: fln

!***********************************************************************
! Initialization of MPI
!***********************************************************************
call mpi_init(ierr)
if (ierr.ne.mpi_success) then
	print *,'error starting mpi program. terminating.'
	call mpi_abort(mpi_comm_world, rc, ierr)
end if

call mpi_comm_rank(mpi_comm_world, rank, ierr) 
call mpi_comm_size(mpi_comm_world, nt, ierr)

write(*,*) "program starts with rank", rank, nt
if(rank.eq.0) call cpu_time(t1)

!***********************************************************************
! Initialization of system
!***********************************************************************


! AM: Set maximum and minimum values for cell coords.
maxx=-1000.0
minx=1000.0
maxy=-1000.0
miny=1000.0 

! Set time steps and time checkpoints 
dt=0.001 
div=1000

! Flags for array reallocation and statistics counters
reallocate = 0
count_node=0
count_border=0
count_dens=0

! AM: Set interval for boundary slicing
delta_x=5.0
delta_y=5.0
do i=1,100000
	temp_x(i)=0.0
	temp_y(i)=0.0
end do

call allo_gr_cells(nt) ! all ranks allocate the counters for growing cells
!***********************************************************************

!***********************************************************************

! AM 02/05/2012: open files to output data
! AM on 19/08/12: master process opens files for output
if (rank.eq.0) then
	if (wr_psfil.eq.1) then ! AM: record movie in postscript format
		open(unit=90,file='prp',status='unknown') 
		open(unit=55,file='psfil',status='unknown') 
	end if
	if (wr_traj.eq.1)&
		&open(unit=12,file='stat',status='unknown') ! AM: record cell coordinates
	if (wr_ncell.eq.1)&
		&open(unit=13,file='ncell',status='unknown') ! AM: record density snapshots
	if (wr_node.eq.1)&
	 	&open(unit=14,file='nodcord',status='unknown') ! AM: record node statistics
	if (wr_border.eq.1)&
		&open(unit=15,file='borders',status='unknown') ! AM : record system borders
	if (wr_gr.eq.1)&
		&open(unit=18,file='growth',status='unknown') ! AM : record system borders
end if 

!***********************************************************************

6    format(4f12.3)

!***********************************************************************

!***********************************************************************

! AM: if there is a continuation file, master processor reads from it
if (rank.eq.0) then

	if(readcont.eq.1) then
		call init_data()
	else
		! AM: generate random numbers to set a system of randomly positioned cells
		seed=5513974+rank*32570
		call init_genrand(seed)
	
		allocate(ran1(100000))
		do i=1,size(ran1)
			ran1(i)=grnd()
		end do

		nns=90 ! AM: start from nns cells
		dst=0.0
		cc=1
		ncc(cc)=1
		lx=25
		ly=50

		! AM: create randomly distributed cell center of masses in the simulation box
		do  i=3,10000,2
			noa=0 
			!this is done to avoid initial overlapping at t=0
			do  j=1,cc
				dst=sqrt(((lx-1.2)*ran1(i)-(lx-1.2)*ran1(ncc(j)))**2 & ! am: cell diameter is 1.0
				& +((ly-1.2)*ran1(i+1)-(ly-1.2)*ran1(ncc(j)+1))**2)
			if (dst.le.1.0) noa=1
			end do

			if (noa.eq.0) then
				cc=cc+1
				ncc(cc)=i
			end if
	
		end do
		
		call allo(nns)   
		do  ni=1,nns ! loop over cells	
			x(ni)=(lx-1.1)*ran1(ncc(ni))
			y(ni)=(ly-1.1)*ran1(ncc(ni)+1)
		end do ! end loop over cells
		
	    ! deallocate random variable array
	    deallocate(ran1)
		! after generating cells, assign specific radii
		call assign_cell_type(1,nns)
	end if ! end readcont=1 condition
	
end if	! end rank=0 condition

if(nt.ne.1) then
if (rank.eq.0) then ! master processor sends total number of cells to all processes
	do i=1, nt-1
		dest = i
		tag=1
		call mpi_send(nns,1,mpi_real,&
		&dest,tag,mpi_comm_world,ierr)
	end do
else ! all processes receive total number of cells
	tag=1
	call mpi_recv(nns,1,mpi_real,&
	&0,tag,mpi_comm_world,stat,ierr)
end if
end if ! nt.eq.1
! AM: all processors allocate arrays for total cell position data
! after allocation exchange data
if (rank.ne.0) call allo(nns) ! rank 0 allocated data arrays in init_data call

if(nt.ne.1) then ! if i=2 and more processors, share the data
if (rank.eq.0) then ! master processor sends coords, etc data to all processes
	do i=1, nt-1
		dest = i

		tag=1
		call mpi_send(x,size(x),mpi_real,&
		&dest,tag,mpi_comm_world,ierr)

		tag=2
		call mpi_send(y,size(y),mpi_real,&
		&dest,tag,mpi_comm_world,ierr)
		
		tag=3
		call mpi_send(radius,size(radius),mpi_real,&
		&dest,tag,mpi_comm_world,ierr)
		
		tag=4
		call mpi_send(delta,size(delta),mpi_real,&
		&dest,tag,mpi_comm_world,ierr)
		
		tag=5
		call mpi_send(cell_type,size(cell_type),mpi_integer,&
		&dest,tag,mpi_comm_world,ierr)

	end do

else ! all the processes receive cell related data

	tag=1
	call mpi_recv(x,size(x),mpi_real,&
	&0,tag,mpi_comm_world,stat,ierr)
		
	tag=2
	call mpi_recv(y,size(y),mpi_real,&
	&0,tag,mpi_comm_world,stat,ierr)
	
	tag=3
	call mpi_recv(radius,size(radius),mpi_real,&
	&0,tag,mpi_comm_world,stat,ierr)
	
	tag=4
	call mpi_recv(delta,size(delta),mpi_real,&
	&0,tag,mpi_comm_world,stat,ierr)
	
	tag=5
	call mpi_recv(cell_type,size(cell_type),mpi_integer,&
	&0,tag,mpi_comm_world,stat,ierr)

end if ! end rank=0 condition
end if ! nt.ne.1	

! Calculate simulation box size after loading data
do ni=1, nns ! loop over all cells
	if (maxx.lt.x(ni)) then
		maxx=x(ni)
	end if
	if (minx.gt.x(ni)) then
		minx=x(ni)
	end if
	if (maxy.lt.y(ni)) then
		maxy=y(ni)
	end if
	if (miny.gt.y(ni)) then
		miny=y(ni)
	end if
end do	! end of loop over all cells

	
! Box size of system
lx=maxx-minx
ly=maxy-miny

!***********************************************************************
  
!***********************************************************************
! AM: initial coordinates
do ni=1, nns ! loop over cells
	xm(ni)=x(ni)! assign previous positions
	ym(ni)=y(ni)
end do ! end loop over cells


!***********************************************************************

!***********************************************************************

! AM: arrays for expand communication
if(nt.ne.1) then
do i=1,iter ! loop over iteration steps
	do k=0,nt-1 ! reset processor pairing for given iteration
		pair_proc(i,k)=-1
	end do
	do k=0,nt-1
		if(pair_proc(i,k).eq.-1) then ! there was not pairing for this process. assign a pair
			pair_proc(i,k)=k+2**(i-1) 
			pair_proc(i,k+2**(i-1))=k
		end if
	end do 
end do

! AM: assign portion of atoms to processors
call distrib_proc(nnsproc,rank)
call start_end_index(nnsproc,pair_proc)

end if ! nt.ne.1

if(nt.eq.1) then
nstart=1
nend=nns
end if

! AM: generate seed for random fluctuations fluctuations

seed2=783257+rank*4855
call init_genrand(seed2)
!***********************************************************************
! System initialization ends here
!***********************************************************************

call mpi_barrier(mpi_comm_world,stat)

!***********************************************************************
! integration steps begin
!***********************************************************************


do is=1,steps ! loop over time
!   
	if(rank.eq.0.and.mod(is,100).eq.0) write(*,*) "step=",is,"nns=",nns
	
	if (mod(is,10).eq.0.or.is.eq.1.or.mod(is,div).eq.1) then ! am: 16/7/12	
		call dist()
		!call ndist()
	end if

	oldnns=nns
	
	! AM on 27/06/14: generate random gaussians for com's x,y components
	! ran2(i) is noise in x an ran2(i+1) is noise in y components correspondingly
	do ni=nstart,nend ! create nns random numbers for x and y components, i.e. 2*nns
	 	rr = 2.d0
    	do while( rr >= 1.d0)
       		x1 = grnd()
       		x2 = grnd()
       		x1 = 2.d0*x1 - 1.d0
       		x2 = 2.d0*x2 - 1.d0
       		rr = x1*x1 + x2*x2
    	end do

    	rr = sqrt((-2.d0*log(rr))/rr)
    	fluct(2*ni-1) = x1*rr
    	fluct(2*ni-0) = x2*rr
	end do
	
	call calc_f(x,y,nstart,nend,is)
	
	do ni=nstart,nend
	
	    ! keep track of old values of X,Y before integration
     	xm(ni)=x(ni) 
		ym(ni)=y(ni)
		
		d(ni)=diffc(diff_max,diff_min,ni)
		c(ni)=1.0/d(ni)
		
		if(ni.eq.ni_node) then
		if (cell_type(ni_node).eq.1) then
            d(ni_node)=diff_max/4.0
			c(ni_node)=1.0/d(ni_node)
        else
            d(ni_node)=diff_min/4.0
            c(ni_node)=1.0/d(ni_node)
		end if
        end if

		x(ni)=x(ni)+sqrt(dt)*sqrt(d(ni))*fluct(2*ni-1) ! prediction: X, Y are temporary values
		y(ni)=y(ni)+sqrt(dt)*sqrt(d(ni))*fluct(2*ni-0) ! prediction: X, Y are temporary values
			
		dtmp(ni)=diffc(diff_max,diff_min,ni) ! updated values for random fluctuation term
		ctmp(ni)=1.0/dtmp(ni)
		
		if(ni.eq.ni_node) then
        if(cell_type(ni_node).eq.1) then
			dtmp(ni_node)=diff_max/4.0
			ctmp(ni_node)=1.0/dtmp(ni_node)
        else
            dtmp(ni_node)=diff_min/4.0
            ctmp(ni_node)=1.0/dtmp(ni_node)
        end if
        end if

		xp(ni)=xm(ni)+(fx(ni)/c(ni))*dt+&
				&0.5*(sqrt(dtmp(ni))+sqrt(d(ni)))*sqrt(dt)*fluct(2*ni-1) ! correction: update old value for x
		yp(ni)=ym(ni)+(fy(ni)/c(ni))*dt+&
				&0.5*(sqrt(dtmp(ni))+sqrt(d(ni)))*sqrt(dt)*fluct(2*ni-0) ! correction: update old value for y

		x(ni)=xp(ni)
		y(ni)=yp(ni)
		
		vx(ni)=(x(ni)-xm(ni))/dt
		vy(ni)=(y(ni)-ym(ni))/dt
		kin(ni)=0.5*(vx(ni)**2+vy(ni)**2);
		
	end do ! end predictor loop

   call mpi_barrier(mpi_comm_world,stat)
   
   call dens(nstart,nend)
   
   call grow_cells(nstart, nend, is,rank)
   
   if(mod(is,1000).eq.0) then
   call mpi_barrier(mpi_comm_world,stat)

   ! AM: master process collects all the growing cell data
		if(nt.ne.1) then
		if (rank.eq.0) then ! master processor receive data
			do i=1, nt-1
				source=i
				tag=i
				
				call mpi_recv(growing_cells(i+1),1,mpi_integer,&
				&source,tag,mpi_comm_world,stat,ierr)
				
			end do

		else ! all the processes send data to master process
            dest =0
			tag=rank
			
		    call mpi_send(growing_cells(rank+1),1,mpi_integer,&
			&0,tag,mpi_comm_world,ierr)

	    end if ! end rank=0 condition
		end if
   end if
   
   if(mod(is,10000).eq.0) then
   		! allocate and record the current radii of all cells
   		call assign_cell_type(nstart,nend)
   		call assign_soft_radii(nstart,nend)
   	
   		radius(ni_node)=node_r
   		delta(ni_node)=knode*node_r

   end if

	! AM: my implementation of expand communication
	if(nt.ne.1) then
	do i=1,iter

		dest=pair_proc(i,rank)
		source=pair_proc(i,rank)
		
		
		
		call ini_buffer(start_index_send(i,rank),end_index_send(i,rank),&
		&start_index_receive(i,rank),end_index_receive(i,rank))
		
	! send and receive data in paired processors. 

	if(rank.lt.dest) then ! message transports forwards first then recieves
		
	    tag=1
		call mpi_send(xts,size(xts),mpi_real,&
			&dest,tag,mpi_comm_world,ierr)
		
		call mpi_recv(xtr,size(xtr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)
		
		tag=2
		call mpi_send(xmts,size(xmts),mpi_real,&
			&dest,tag,mpi_comm_world,ierr)
		
		call mpi_recv(xmtr,size(xmtr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)
		
		tag=3
		call mpi_send(yts,size(yts),mpi_real,&
			&dest,tag,mpi_comm_world,ierr)
		
		call mpi_recv(ytr,size(ytr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)
		
		tag=4
		call mpi_send(ymts,size(ymts),mpi_real,&
			&dest,tag,mpi_comm_world,ierr)
		
		call mpi_recv(ymtr,size(ymtr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)
		
		tag=5
		call mpi_send(radts,size(radts),mpi_real,&
			&dest,tag,mpi_comm_world,ierr)
		
		call mpi_recv(radtr,size(radtr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)
		
		tag=6
		call mpi_send(delts,size(delts),mpi_real,&
			&dest,tag,mpi_comm_world,ierr)
		
		call mpi_recv(deltr,size(deltr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)
		
		tag=7
		call mpi_send(vxts,size(vxts),mpi_real,&
			&dest,tag,mpi_comm_world,ierr)
		
		call mpi_recv(vxtr,size(vxtr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)
		
		tag=8
		call mpi_send(vyts,size(vyts),mpi_real,&
			&dest,tag,mpi_comm_world,ierr)
		
		call mpi_recv(vytr,size(vytr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)
		
		tag=9
		call mpi_send(cellts,size(cellts),mpi_integer,&
			&dest,tag,mpi_comm_world,ierr)
		
		call mpi_recv(celltr,size(celltr),mpi_integer,&
		&source,tag,mpi_comm_world,stat,ierr)

	else ! message receives first and then passed backwards

		tag=1
		call mpi_recv(xtr,size(xtr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)

		call mpi_send(xts,size(xts),mpi_real,&
		&dest,tag,mpi_comm_world,ierr)
		
		tag=2
		call mpi_recv(xmtr,size(xmtr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)

		call mpi_send(xmts,size(xmts),mpi_real,&
		&dest,tag,mpi_comm_world,ierr)
		
		tag=3
		call mpi_recv(ytr,size(ytr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)

		call mpi_send(yts,size(yts),mpi_real,&
		&dest,tag,mpi_comm_world,ierr)
		
		tag=4
		call mpi_recv(ymtr,size(ymtr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)

		call mpi_send(ymts,size(ymts),mpi_real,&
		&dest,tag,mpi_comm_world,ierr)
		
		tag=5
		call mpi_recv(radtr,size(radtr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)

		call mpi_send(radts,size(radts),mpi_real,&
		&dest,tag,mpi_comm_world,ierr)
		
		tag=6
		call mpi_recv(deltr,size(deltr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)

		call mpi_send(delts,size(delts),mpi_real,&
		&dest,tag,mpi_comm_world,ierr)
		
		tag=7
		call mpi_recv(vxtr,size(vxtr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)

		call mpi_send(vxts,size(vxts),mpi_real,&
		&dest,tag,mpi_comm_world,ierr)
		
		tag=8
		call mpi_recv(vytr,size(vytr),mpi_real,&
		&source,tag,mpi_comm_world,stat,ierr)

		call mpi_send(vyts,size(vyts),mpi_real,&
		&dest,tag,mpi_comm_world,ierr)
		
		tag=9
		call mpi_recv(celltr,size(celltr),mpi_integer,&
		&source,tag,mpi_comm_world,stat,ierr)

		call mpi_send(cellts,size(cellts),mpi_integer,&
		&dest,tag,mpi_comm_world,ierr)
		

	end if

		! dump buffer into actual variables
		call convert_buffer(start_index_receive(i,rank),end_index_receive(i,rank))
		
		call mpi_barrier(mpi_comm_world,stat) ! synchronize before next swap of coordinates
		
		call del_buffer() ! deallocate buffer
	end do
	
	endif ! nt.ne.1

! AM: cell division
!********************************************************************
	if (mod(is,div).eq.0) then 
		if (rank.eq.0) call divide(rank)	
		call mpi_barrier(mpi_comm_world,stat) ! synchronize before next swap of coordinates		
		if(nt.ne.1.and.rank.eq.0) then
			! communicate data between processors
			do i=1, nt-1
				dest = i
			
				tag=10
				call mpi_send(reallocate,1,mpi_real,&
				&dest,tag,mpi_comm_world,ierr)
				
				tag=1
				call mpi_send(nns,1,mpi_real,&
				&dest,tag,mpi_comm_world,ierr)

				tag=2
				call mpi_send(x,size(x),mpi_real,&
				&dest,tag,mpi_comm_world,ierr)
				
				tag=3
				call mpi_send(xm,size(xm),mpi_real,&
				&dest,tag,mpi_comm_world,ierr)

				tag=4
				call mpi_send(y,size(y),mpi_real,&
				&dest,tag,mpi_comm_world,ierr)
				
				tag=5
				call mpi_send(ym,size(ym),mpi_real,&
				&dest,tag,mpi_comm_world,ierr)
				
				tag=6
				call mpi_send(radius,size(radius),mpi_real,&
				&dest,tag,mpi_comm_world,ierr)
				
				tag=7
				call mpi_send(delta,size(delta),mpi_real,&
				&dest,tag,mpi_comm_world,ierr)
				
				tag=8
				call mpi_send(cell_type,size(cell_type),mpi_integer,&
				&dest,tag,mpi_comm_world,ierr)

			end do
		end if ! condition nt!=1 and rank=0
		if(nt.ne.1.and.rank.ne.0) then
		
			tag=10
			call mpi_recv(reallocate,1,mpi_real,&
			&0,tag,mpi_comm_world,stat,ierr)
			
			tag=1
			call mpi_recv(nns,1,mpi_real,&
			&0,tag,mpi_comm_world,stat,ierr)
			
			if(reallocate.eq.1) then
			 call reallo(oldnns,nns+1)
			end if
		
			tag=2
			call mpi_recv(x,size(x),mpi_real,&
			&0,tag,mpi_comm_world,stat,ierr)
			
			tag=3
			call mpi_recv(xm,size(xm),mpi_real,&
			&0,tag,mpi_comm_world,stat,ierr)
		
			tag=4
			call mpi_recv(y,size(y),mpi_real,&
			&0,tag,mpi_comm_world,stat,ierr)
			
			tag=5
			call mpi_recv(ym,size(ym),mpi_real,&
			&0,tag,mpi_comm_world,stat,ierr)
						
			tag=6
			call mpi_recv(radius,size(radius),mpi_real,&
			&0,tag,mpi_comm_world,stat,ierr)
			
			tag=7
			call mpi_recv(delta,size(delta),mpi_real,&
			&0,tag,mpi_comm_world,stat,ierr)
			
			tag=8
			call mpi_recv(cell_type,size(cell_type),mpi_integer,&
			&0,tag,mpi_comm_world,stat,ierr)

		end if ! nt.ne.1.and.rank.ne.0	
		call mpi_barrier(mpi_comm_world,stat) ! synchronize before next swap of coordinates
		
		if(nt.ne.1) then
			call distrib_proc(nnsproc,rank)
			call start_end_index(nnsproc,pair_proc)		
        endif ! nt.eq.1
	end if ! end condition for division

!********************************************************************
reallocate=0
	! AM: master process writes output data 
	if(rank.eq.0) then ! check whether ist's master process
	
		if(wr_traj.eq.1) then
		!****** OUTPUT COORDINATES ******
        	if(mod(is,100).eq.0) then
        		do ni=1,nns
          			write(12,*) is,ni,x(ni),y(ni),x(ni_node),y(ni_node)
        		end do
        	end if
        end if
        
        if(writecont.eq.1) then
        !****** OUTPUT CONTINUATION FILE ******
        	if (mod(is,1000).eq.0) then
				open(unit=65,file='writecont',status='unknown')
				do	ni=1,nns
					write(65,*) x(ni),y(ni),radius(ni)
				end do
				close(65)
			end if
		end if
        
        if(wr_psfil.eq.1) then
        !****** OUTPUT PSFIL ******
			if ((mod(is,1000).eq.0)) then		
    			open(unit=45,file='psp',status='unknown')
				do ni=1,nns
					write(45,*) x(ni),y(ni),radius(ni),delta(ni),x(ni_node),y(ni_node),cell_type(ni)
				end do
				close(45)		
				call psnet(nns)
			end if
		end if
		
		if(wr_node.eq.1) then
        !****** OUTPUT NODE POSITIONS ******
			if(mod(is,1000).gt.0) then	
				count_node=count_node+1
				node_data(count_node,1)=x(ni_node)
				node_data(count_node,2)=y(ni_node)		
			else
				! calculate statistics
				av_nodex=0.0
				std_nodex=0.0
				av_nodey=0.0
				std_nodey=0.0
				count_node=0
				! calculate averages
				do i=1,1000
					av_nodex=av_nodex+node_data(i,1)/(1.0*1000)
					av_nodey=av_nodey+node_data(i,2)/(1.0*1000)
				end do
		
				! calculate standard deviation
				do i=1,1000
					std_nodex=std_nodex+(node_data(i,1)-av_nodex)**2/(1.0*(1000-1))
					std_nodey=std_nodey+(node_data(i,2)-av_nodey)**2/(1.0*(1000-1))
				end do
				std_nodex=sqrt(std_nodex)
				std_nodey=sqrt(std_nodey)
			
				if (is.ge.2000) write(14,*) is, av_nodex, std_nodex, av_nodey, std_nodey
				! empty the stat arrays
				do i=1,1000
					node_data(count_node,1)=0.0
					node_data(count_node,2)=0.0
				end do
			end if
		end if
		
		if(wr_border.eq.1) then
		!****** OUTPUT BORDERS OF SYSTEM ******
			if (mod(is,1000).eq.0) then ! AM: calculate the lateral borders of the system
				! calculate the rough positions of the ap border
				! by identifying 10 max and min cell positions and averaging them

				xmin_border=0.0
				xmax_border=0.0
				do i=1,nns
					temp_x(i)=x(i)
				end do
				call sort(temp_x, nns)
				do ni=1,5
					xmin_border=xmin_border+temp_x(ni)/5.0
				end do
				do ni=nns,nns-4,-1
					xmax_border=xmax_border+temp_x(ni)/5.0
				end do
			
				! slice the boundary along the x
				i=1
		
				do while ((temp_x(1)+i*delta_x).le.temp_x(nns))
					ymin_border=0.0
					ymax_border=0.0
					count_border=0
					do ni=1,nns
						if(x(ni).gt.(temp_x(1)+(i-1)*delta_x).and.&
						x(ni).le.(temp_x(1)+i*delta_x)) then
							count_border=count_border+1
							temp_y(count_border)=y(ni)
						end if
					end do
					call sort(temp_y, count_border)
					do ni=1,5
						ymin_border=ymin_border+temp_y(ni)/5.0
					end do
					do ni=count_border,count_border-4,-1
						ymax_border=ymax_border+temp_y(ni)/5.0
					end do
					write(15,*) is, ((temp_x(1)+(i-1)*delta_x)+temp_x(1)+i*delta_x)/2,&
									 &ymin_border, ymax_border,&
									&av_nodex, av_nodey
					i=i+1
				end do
			end if
		end if
		
		if (wr_dens.eq.1) then
			!****** OUTPUT NUMBER DENSITY ******
			if((is.gt.50000).and.(is.le.70000)) then
			if (mod(is,1000).eq.0) then ! AM: calculate the lateral borders of the system
				! calculate the rough positions of the ap border
				! by identifying 10 max and min cell positions and averaging them

				xmin_border=0.0
				xmax_border=0.0
				ymin_border=0.0
				ymax_border=0.0
				
				do i=1,nns
					temp_x(i)=x(i)
					temp_y(i)=y(i)
				end do
				call sort(temp_x, nns)
				call sort(temp_y, nns)
				do ni=1,5
					xmin_border=xmin_border+temp_x(ni)/5.0
				end do
				do ni=nns,nns-4,-1
					xmax_border=xmax_border+temp_x(ni)/5.0
				end do
				do ni=1,5
					ymin_border=ymin_border+temp_y(ni)/5.0
				end do
				do ni=nns,nns-4,-1
					ymax_border=ymax_border+temp_y(ni)/5.0
				end do
				
				! AM: slice the embryo and calculate number of cells in slices
				imax=ceiling((xmax_border-xmin_border)/delta_x)
				jmax=ceiling((ymax_border-ymin_border)/delta_y)
				
				allocate(density(imax,jmax))
				do i=1,imax
					do j=1,jmax
						density(i,j)=0
					end do
				end do
				
				do i=1,imax
					do j=1,jmax
						do ni=1,nns
							if((x(ni).gt.(xmin_border+(i-1)*delta_x)).and.&
							&(x(ni).le.(xmin_border+i*delta_x))) then
								if((y(ni).gt.(ymin_border+(j-1)*delta_y)).and.&
								&(y(ni).le.(ymin_border+j*delta_y))) then
									if(ni.ne.ni_node) then
										density(i,j)=density(i,j)+1
									else
										density(i,j)=density(i,j)+1
									end if
								end if
							end if
						end do
					end do
				end do

				write(fln, "(A4,I5)") "dens", is

				open(unit=16,file=fln,status='UNKNOWN') ! AM 26/09/12: store node coordinates
				write(16,*) xmin_border, ymin_border, imax, jmax
				do i=1,imax
					do j=1,jmax
						write(16,*) i, j, density(i,j)
					end do
				end do 
				close(16)
                deallocate(density)
                
			end if
		end if
		end if
		
		if (wr_ncell.eq.1) then
			if (mod(is,1000).eq.0) then ! every 1000 steps count number of cells in ap
				na=0
				np=0
				do ni=1, nns
					if (x(ni).le.x(ni_node)) then
						na=na+1
					else 
						np=np+1
					end if
				end do		
				write(13,*) is, na, np 
			end if
		end if
		if (wr_gr.eq.1) then
        	if(mod(is,1000).eq.0) then
        		ngr=0
        		do i=1,nt 
        			ngr=ngr+growing_cells(i)/1000.0
        		end do
        		write(18,*) is, ndiv_cell, ngr
        		
        	end if
        end if
		
	end if ! end check for master process
	
	call mpi_barrier(mpi_comm_world,stat)

end do ! AM: times step loop is done

close(12)
close(13)
close(14)
close(15)

if(rank.eq.0) then 
call cpu_time(t2)
write(*,*), "time of run", t2-t1
end if
call mpi_finalize(ierr)
stop
end 














