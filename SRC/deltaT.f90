!##########################################################
! SUBROUTINE deltaT
! Local time step scheme for steady flows
!  Goal: fill the array dt with the local time step for each cell
!
!        dt =                   A * CFL 
!             -----------------------------------------------
!            Nsides   
!             SUM ( ABS [ u * n_jx + v * n_jy ] + c * L )
!             j=1
!  where A is the area of the cell,
!        CFL is the Courant-Friedrich-Lewy number,
!        u,v are the x,y velocity components inside the cell, 
!        n_jx,n_jy are the x,y components of the normal to each side of the cell (not normalized),
!        c is the wave speed sqrt(g*h)
!        L is the perimeter of the cell
!        Nsides is the number of sides for the cell
!##########################################################
SUBROUTINE deltaT
    USE module_shallow
    IMPLICIT NONE

    INTEGER(ki) :: i, j, max
    REAL(kr)    :: dx, val,valMax, u, v, h, c, surf
    REAL(kr)    :: nx, ny
    REAL(kr)    :: norm, absv
    
    absv   = zero
    valMax = zero
    DO i=1,nbrElem
       h = U0(i*nbvar-2)
       IF (h<eps) h=eps
       u = U0(i*nbvar-1)/h
       v = U0(i*nbvar)/h
       c = SQRT(ggrav*h)
       
       surf = geom_data(i,1)
       val = zero
       DO j=1,nbr_nodes_per_elem(i)
         nx = cell_data_n(i,2*j-1)
         ny = cell_data_n(i,2*j)
         norm = SQRT(nx**2 + ny**2)
         val = val + ABS(u*nx + v*ny) + c*norm
       ENDDO
       
       dt(i*nbvar-2:i*nbvar) = surf*CFL/val
       
    END DO

END SUBROUTINE deltaT
