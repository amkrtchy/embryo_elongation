! Module of simulation parameters

module params

implicit none

! AM: MPI related flags
! inter - defines number of processors when running code in parallel via 2^inter = number_of_processors
integer, parameter :: max_proc=64, iter=3


integer, parameter :: steps=150000 ! sim. steps
integer, parameter :: ni_node = 262 ! node id. coordinated with the writecont so that node ~ in the center of sim. system

! AM: flags to write statistics
integer :: writecont=0, readcont=1  ! record configuration of the system. Read initial configuration of the system as in writecont
integer :: wr_node=1, wr_dens=1 ! record node position, record cell number density
integer :: wr_border=1, wr_traj=1 ! record lateral borders (deprecated), record positions of all cells in the system
integer :: wr_ncell=1, wr_psfil=1 ! record number of anterior/posterior cells in the system
integer :: wr_gr=1 ! record growth rate of system (deprecated)
integer :: n_align=0 

! AM: Define pi
real, parameter :: PI=ACOS(-1.0)

! AM: cell-cell interactions
real,parameter :: stratr=10.0, strsc=1000.0,strhc=10000.0 ! strengths of cell-cell attraction, soft core and hard core repulsion
real, parameter :: softr=25.0, growra=25.0, growrp=15.0 ! length of density/motility gradient, anterior and posterior growth regions
real, parameter :: node_r = 2.0, tnode_r=2.8, knode=(tnode_r-node_r)/node_r ! hard core/ total radii of the node
real, parameter :: r0=0.5, r1=0.5, tr0=0.7, tr1=0.9, k0=(tr0-r0)/r0, k1=(tr1-r1)/r1 ! hard core/total radii of anterior/posterior cells

! AM: Strength of neighbour coupling - same as lower diffusion coefficient
real, parameter :: beta=1.0

! AM: diffusion coefficients
real, parameter :: diff_min=0.01, diff_max=0.1

! AM: System "walls": AP_right, AP-left, lateral_top, lateral_bottom
! walls, initial config of system (writecont) and node position/id are coordinated so that node appears at center.
real, parameter :: delrx=32.0, dellx=-9.0, delry=24.0,delly=-11.0

! AM: MPI communication between processors
integer :: start_index_send(iter,0:(max_proc)),end_index_send(iter,0:(max_proc))
integer :: start_index_receive(iter,0:(max_proc)),end_index_receive(iter,0:(max_proc))

real :: xcm, ycm

real, dimension(:), allocatable :: x,y,xtot,ytot,xm,ym,xp,yp,vx,vy ! cell coordinates/ velocities
real, dimension(:), allocatable :: radius ! cell hard core radii
real, dimension(:), allocatable :: fx,fy,pot,kin ! forces and energies
real, dimension(:), allocatable :: xts,xmts,yts,ymts,xtr,xmtr,ytr,ymtr ! send/receive buffer for MPI
real, dimension(:), allocatable :: radts,radtr ! buffer for send/receive MPI
real, dimension(:), allocatable :: vxts,vxtr,vyts,vytr! buffer for send/receive MPI
integer, dimension(:), allocatable :: cellts,celltr ! buffer for send/receive MPI
real, dimension(:), allocatable :: delta,deltr,delts
real, dimension(:), allocatable :: d,c,dtmp,ctmp ! cell diffusion coefficient

real, dimension(:), allocatable :: fluct ! fluctuations
real, dimension(:), allocatable :: av_rho, av_vx, av_vy ! average local density and velocity for each cell

integer :: co,allo_length,nns,nt,nstart,nend,reallocate
integer :: nco
integer, dimension(:), allocatable :: growing_cells
integer :: ndiv_cell

integer, dimension(:,:), allocatable :: rd, nrd ! vars used for neighbour list calculations
integer, dimension(:), allocatable :: cell_type ! cell type array

end module params
