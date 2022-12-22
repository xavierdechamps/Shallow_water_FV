!##########################################################
! SUBROUTINE get_relative_indicator_amr
!##########################################################
SUBROUTINE amr_get_relative_indicator(gradX,gradY)
    USE module_shallow
    USE OMP_LIB
    IMPLICIT NONE

    ! Subroutine parameters
    REAL(kr), DIMENSION(1:nbvar*nbrElem), INTENT(IN) :: gradX, gradY
    
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

END SUBROUTINE amr_get_relative_indicator

SUBROUTINE amr_split_cells()
    USE module_shallow
    IMPLICIT NONE
    
    ! Local parameters
    INTEGER(ki) :: i, j, typec,newcells,newnodes, inew, ccur, cnew
    INTEGER(ki), DIMENSION(1:nbrElem) :: flag
    REAL(kr)    :: xnew,ynew
    
    REAL(kr), ALLOCATABLE :: noden(:,:) ! new node coordinates
    REAL(kr), ALLOCATABLE :: U0n(:) ! new data
    INTEGER(ki), ALLOCATABLE :: elemn(:,:)  ! new cell data
    INTEGER(ki), ALLOCATABLE    :: nbr_nodes_per_elemn(:)
    
    newcells = 0
    newnodes = 0
    flag     = 0
    j = 0
    DO i=1,nbrElem
      ! Test if the cell must be split
      IF (shock_indicator(i) .GE. 0.1d00) THEN
        typec = nbr_nodes_per_elem(i) ! 3=triange, 4=quadrangle
        IF (typec.EQ.3) THEN
          newcells = newcells + 2
        ELSEIF (typec.EQ.4) THEN
          newcells = newcells + 3
        ENDIF
        newnodes = newnodes + 1
        
        j = j + 1
        flag(j) = i
      ENDIF
    ENDDO
    
    ALLOCATE(noden(1:nbrNodes+newnodes,1:2))
    ALLOCATE(elemn(1:nbrElem+newcells,1:5))
    ALLOCATE(U0n(1:nbvar*(nbrElem+newcells)))
    ALLOCATE(nbr_nodes_per_elemn(1:nbvar*(nbrElem+newcells)))
    
    noden(1:nbrNodes,1:2) = node(1:nbrNodes,1:2)
    elemn(1:nbrElem,1:5)  = elem(1:nbrElem,1:5)
    U0n(1:nbvar*nbrElem)        = U0(1:nbvar*nbrElem)
    nbr_nodes_per_elemn(1:nbrElem) = nbr_nodes_per_elem(1:nbrElem)
    
    cnew = nbrElem+1
        
    DO i=1,newnodes
      inew = nbrNodes+i ! new node index
      ccur = flag(i) ! current cell index
      
      noden(inew,1:2) = geom_data(ccur,3:4) ! coordinates of new node = cell center
      
      typec = nbr_nodes_per_elem(ccur) ! 3=triange, 4=quadrangle
            
      IF (typec.EQ.3) THEN
        elemn(ccur  ,1:5) = (/ elem(ccur,1) , elem(ccur,2) , inew , 0 , elem(ccur,5) /)
        elemn(cnew  ,1:5) = (/ elem(ccur,2) , elem(ccur,3) , inew , 0 , elem(ccur,5)  /)
        elemn(cnew+1,1:5) = (/ elem(ccur,3) , elem(ccur,1) , inew , 0 , elem(ccur,5)  /)
        U0n(nbvar*(cnew-1)+1:nbvar*(cnew-1)+3) = U0(nbvar*(ccur-1)+1:nbvar*(ccur-1)+3)
        U0n(nbvar*(cnew  )+1:nbvar*(cnew  )+3) = U0(nbvar*(ccur-1)+1:nbvar*(ccur-1)+3)
                
        nbr_nodes_per_elemn(cnew:cnew+1) = 3
      
      ELSEIF (typec.EQ.4) THEN
        elemn(ccur  ,1:5) = (/ elem(ccur,1) , elem(ccur,2) , inew , 0 , elem(ccur,5) /)
        elemn(cnew  ,1:5) = (/ elem(ccur,2) , elem(ccur,3) , inew , 0 , elem(ccur,5)  /)
        elemn(cnew+1,1:5) = (/ elem(ccur,3) , elem(ccur,4) , inew , 0 , elem(ccur,5)  /)        
        elemn(cnew+2,1:5) = (/ elem(ccur,4) , elem(ccur,1) , inew , 0 , elem(ccur,5)  /) 
        U0n(nbvar*(cnew-1)+1:nbvar*(cnew-1)+3) = U0(nbvar*(ccur-1)+1:nbvar*(ccur-1)+3)
        U0n(nbvar*(cnew  )+1:nbvar*(cnew  )+3) = U0(nbvar*(ccur-1)+1:nbvar*(ccur-1)+3)       
        U0n(nbvar*(cnew+1)+1:nbvar*(cnew+1)+3) = U0(nbvar*(ccur-1)+1:nbvar*(ccur-1)+3)       
      
        nbr_nodes_per_elemn(ccur) = 3
        nbr_nodes_per_elemn(cnew:cnew+2) = 3
        
      ENDIF
      
      cnew = cnew+typec-1 ! update first new cell index
    ENDDO
        
    CALL write_gmsh(U0n,nbrElem+newcells,"amr.msh",7,noden,elemn,front,nbrNodes+newnodes,nbrElem+newcells,0,0,&
&                   nbrFront,nbr_nodes_per_elemn,0,0)
        
    DEALLOCATE(noden)
    DEALLOCATE(elemn)
    DEALLOCATE(U0n)
    DEALLOCATE(nbr_nodes_per_elemn)
    
END SUBROUTINE amr_split_cells