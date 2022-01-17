!##########################################################
! SUBROUTINE mem_allocate
! 
!##########################################################
SUBROUTINE mem_allocate(node,front,elem,U0,depth,BoundCond,lengU0,nbrNodes,nbrElem,nbrFront)
    USE module_shallow, only : kr,ki,nbvar
    IMPLICIT NONE

!   WARNING !!!!!!!!!!!!
! If you change the input/output parameters of this routine
! you also have to change them in the interface 
!      module_shallow.f90 / module_mem_allocate

    REAL(kr), ALLOCATABLE    :: U0(:)
    REAL(kr), ALLOCATABLE    :: node(:,:)
    REAL(kr), ALLOCATABLE    :: depth(:)
    REAL(kr), ALLOCATABLE    :: BoundCond(:,:)
    INTEGER(ki), ALLOCATABLE    :: elem(:,:)
    INTEGER(ki), ALLOCATABLE    :: front(:,:)
    
    INTEGER(ki), INTENT(IN) :: lengU0,nbrNodes,nbrElem,nbrFront
    
    ALLOCATE(node(1:nbrNodes,1:2))
    ALLOCATE(front(1:nbrFront,1:4))
    ALLOCATE(elem(1:nbrElem,1:4))
    ALLOCATE(U0(1:nbvar*nbrElem))
    ALLOCATE(depth(1:nbrElem))
    ALLOCATE(BoundCond(1:nbrFront,1:3))
    
END SUBROUTINE mem_allocate
