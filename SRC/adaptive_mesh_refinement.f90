!##########################################################
! SUBROUTINE 
!##########################################################
SUBROUTINE get_relative_indicator_amr(gradX,gradY)
    USE module_shallow
    USE OMP_LIB
    IMPLICIT NONE

    ! Subroutine parameters
    REAL(kr), DIMENSION(1:nbvar*nbrElem) :: Urk, Uc, H, gradX, gradY
    
    ! Local parameters
    INTEGER(ki) :: i, j, idL, idR, IDk
    REAL(kr)    :: globmax, tmpmax, tmp
    
    shock_indicator = zero      
    globmax = zero

!$OMP PARALLEL &
!$OMP& default (shared) &
!$OMP& private (tmp,tmpmax,idL,IDk) 
!$OMP DO
    DO i=1,nbrElem
      tmpmax = zero
      
      DO j=1,nbr_nodes_per_elem(i)
        idL = geom_data_ind(i,j)
        IDk = (idL-1)*nbvar + 1  
        tmp = SQRT(gradX(IDk)*gradX(IDk) + gradY(IDk)*gradY(IDk))
        tmpmax = MAX(tmp,tmpmax)
      ENDDO
            
      idL = i
      IDk = (idL-1)*nbvar + 1
      tmp = SQRT(gradX(IDk)*gradX(IDk) + gradY(IDk)*gradY(IDk))
      tmpmax = MAX(tmp,tmpmax)
      
      IF (ABS(tmpmax).LT.1.E-12) tmpmax = zero
      globmax = MAX(globmax,ABS(tmpmax))
      shock_indicator(i) = tmpmax
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    shock_indicator = shock_indicator / globmax

END SUBROUTINE get_relative_indicator_amr
