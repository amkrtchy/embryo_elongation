! GENERAL INFO
! Module generated postscript file which is later used to make a simulation movie
! HOW MODULE WORKS
! Each time step is recorded as a single frame in postscript file. 
! At each time step cells' center of masses, their radii and types are passed to this module
! Based on COMs, circles with hard and soft radii are drawn to represent cells.
! Node is plotted in red color


! PSSETUP sets the header of postscript frame

subroutine pssetup (psunit)

implicit none
!your own postscript setups to file psunit
integer :: psunit

write(psunit,*)'%!PS-Adobe-2.0 EPSF-1.2'
write(psunit,*)'%%BoundingBox: -1500 -800 3500 1500'
write(psunit,*)'200 200 translate'
write(psunit,*)'20.5 20.5 scale'
write(psunit,*)'0.05 setlinewidth'

!a circle macro [adobe: postscript language,  tutorial and cookbook p.137] 
write(psunit,*)'/circledict 8 dict def'
write(psunit,*)'circledict /mtrx matrix put'
write(psunit,*)'/circle { circledict begin'
write(psunit,*)'  /rad exch def'
write(psunit,*)'  /y exch def'
write(psunit,*)'  /x exch def'
write(psunit,*)'  /savematrix mtrx currentmatrix def'
write(psunit,*)'  x y translate'
write(psunit,*)'  rad rad scale'
write(psunit,*)'  0 0 1 0 360 arc'
write(psunit,*)'  savematrix setmatrix'
write(psunit,*)'  end } def'

!network picture title with 12-point times roman
write(psunit,*) '/Times-Roman findfont 0.3 scalefont setfont'
write(psunit,*)'0.015 setlinewidth'

!set label font to normat 6 point. due to the 40 times scaling 0.15
!must be given as the point size
write(psunit,*)'/Times-Roman findfont 0.15 scalefont setfont'

end

! PSCIRCLE plots circles centered at (xxx, yyy) and having hard/soft radii
! (xxx, yyy) cell's center of masses
! radius - cell's hard core radius
! delta - fluctuating shell
! ctp - cell type
! (xnode, ynode) - node's position
 subroutine pscircle(xxx,yyy,radius,delta,xnode,ynode,ctp,psunit)

!prints ps-code for a circle to file psunit 
!fill    0 - only contour, transparent
!        1 - filled
!        2 - only contour, opaque  
implicit none

real :: xxx,yyy,radius,delta,xnode,ynode,ctp
integer :: psunit, fill

! plot node

if(xxx.eq.xnode.and.yyy.eq.ynode) then
write(psunit,*) 'newpath'
write(psunit,*) xxx, yyy, radius+0.5*delta, 0, 360, 'arc'
write(psunit,*) 'closepath'
write(psunit,*) '1 0.83 0.83 setrgbcolor fill'
write(psunit,*) 'stroke'
write(psunit,*) 'newpath'
write(psunit,*) xxx, yyy, radius, 0, 360, 'arc'
write(psunit,*) 'closepath'
write(psunit,*) '1 0.67 0.67 setrgbcolor fill'
write(psunit,*) 'stroke'
write(psunit,*) 'newpath'
write(psunit,*) xxx, yyy, radius-1.0*(radius/5.0), 0, 360, 'arc'
write(psunit,*) 'closepath'
write(psunit,*) '1 0.5 0.5 setrgbcolor fill'
write(psunit,*) 'stroke'
write(psunit,*) 'newpath'
write(psunit,*) xxx, yyy, radius-2.0*(radius/5.0), 0, 360, 'arc'
write(psunit,*) 'closepath'
write(psunit,*) '1 0.33 0.33 setrgbcolor fill'
write(psunit,*) 'stroke'
write(psunit,*) 'newpath'
write(psunit,*) xxx, yyy, radius-3.0*(radius/5.0), 0, 360, 'arc'
write(psunit,*) 'closepath'
write(psunit,*) '1 0.17 0.17 setrgbcolor fill'
write(psunit,*) 'stroke'
write(psunit,*) 'newpath'
write(psunit,*) xxx, yyy, radius-4.0*(radius/5.0), 0, 360, 'arc'
write(psunit,*) 'closepath'
write(psunit,*) '1 0 0 setrgbcolor fill'
write(psunit,*) 'fill'
write(psunit,*) 'stroke'
else 
! plot cells other than node
write(psunit,*) 'newpath'
write(psunit,*) xxx, yyy, radius+0.5*delta, 0, 360, 'arc'
write(psunit,*) 'closepath'
write(psunit,*) '0.83 0.83 0.83 setrgbcolor fill'
write(psunit,*) 'stroke'
write(psunit,*) 'newpath'
write(psunit,*) xxx, yyy, radius, 0, 360, 'arc'
write(psunit,*) 'closepath'
write(psunit,*) '0.67 0.67 0.67 setrgbcolor fill'
write(psunit,*) 'stroke'
write(psunit,*) 'newpath'
write(psunit,*) xxx, yyy, radius-1.0*(radius/5.0), 0, 360, 'arc'
write(psunit,*) 'closepath'
write(psunit,*) '0.5 0.5 0.5 setrgbcolor fill'
write(psunit,*) 'stroke'
write(psunit,*) 'newpath'
write(psunit,*) xxx, yyy, radius-2.0*(radius/5.0), 0, 360, 'arc'
write(psunit,*) 'closepath'
write(psunit,*) '0.33 0.33 0.33 setrgbcolor fill'
write(psunit,*) 'stroke'
write(psunit,*) 'newpath'
write(psunit,*) xxx, yyy, radius-3.0*(radius/5.0), 0, 360, 'arc'
write(psunit,*) 'closepath'
write(psunit,*) '0.17 0.17 0.17 setrgbcolor fill'
write(psunit,*) 'stroke'
write(psunit,*) 'newpath'
write(psunit,*) xxx, yyy, radius-4.0*(radius/5.0), 0, 360, 'arc'
write(psunit,*) 'closepath'
write(psunit,*) '0 0 0 setrgbcolor fill'
write(psunit,*) 'fill'
write(psunit,*) 'stroke'
end if 
     
write(psunit,*) '0 0 0 setrgbcolor'
write(psunit,*) 'newpath'
write(psunit,*) xxx, yyy, radius+delta, 0, 360, 'arc'
write(psunit,*) 'closepath'
write(psunit,*) 'stroke'
write(psunit,*) 'newpath'
write(psunit,*) 'closepath'
write(psunit,*) 'stroke'

end
!*****************************************************************

subroutine  psnet(nn)

!creates a postscript file for the final network.
!local variables

implicit none

integer :: nn, psunit,i
real :: x, y,radius,delta,xnode,ynode,ctp

psunit=55

open(unit=42,file='psp',status='unknown')   
call pssetup(psunit)

do i=1,nn
	read(42,*) x,y,radius,delta,xnode,ynode,ctp
    call pscircle(x,y,radius,delta,xnode,ynode,ctp,psunit) ! am 18/06/2014 added  com coordinates 
end do

write(psunit, *) 'showpage'
close(42)       
end














