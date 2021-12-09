!##########################################################
! SUBROUTINE deltaT
! Local time step scheme for steady flows
!  Goal: fill the array dt with the local time step for each cell
!
!        dt =                   A * CFL 
!             -----------------------------------------------
!              3   
!             SUM ( ABS [ u * n_jx + v * n_jy ] + c * L )
!             j=1
!  where A is the area of the cell,
!        CFL is the Courant-Friedrich-Lewy number,
!        u,v are the x,y velocity components inside the cell, 
!        n_jx,n_jy are the x,y components of the normal to each side of the cell (not normalized),
!        c is the wave speed sqrt(g*h)
!        L is the perimeter of the cell
!##########################################################
SUBROUTINE deltaT
    USE module_shallow
    IMPLICIT NONE

    INTEGER(ki) :: i, max
    REAL(kr)    :: dx, val,valMax, u, v, h, c, surf
    REAL(kr)    :: n1_x, n1_y, n2_x, n2_y, n3_x, n3_y
    REAL(kr)    :: norm1, norm2, norm3

    valMax = 0.
    DO i=1,nbrElem
       h = U0(i*nbvar-2)
       IF (h<eps) h=eps
       u = U0(i*nbvar-1)/h
       v = U0(i*nbvar)/h
       c = SQRT(ggrav*h)
       
       surf = geom_data(i,1)
       n1_x = cell_data_n(i,1)
       n1_y = cell_data_n(i,2)
       n2_x = cell_data_n(i,3)
       n2_y = cell_data_n(i,4)
       n3_x = cell_data_n(i,5)
       n3_y = cell_data_n(i,6)
       
       norm1 = SQRT(n1_x**2 + n1_y**2)
       norm2 = SQRT(n2_x**2 + n2_y**2)
       norm3 = SQRT(n3_x**2 + n3_y**2)
       
       val = ABS(u*n1_x+v*n1_y)+ABS(u*n2_x+v*n2_y)+ABS(u*n3_x+v*n3_y)+c*(norm1+norm2+norm3)

       dt(i*nbvar-2:i*nbvar) = surf*CFL/val
       
    END DO

END SUBROUTINE deltaT
