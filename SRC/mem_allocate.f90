!##########################################################
! SUBROUTINE mem_allocate
! 
!##########################################################
SUBROUTINE mem_allocate(node,front,elem,nbr_nodes_per_elem,U0,depth,BoundCond,dt,Source,&
&                       edges,fnormal,geom_data,geom_data_ind,cell_data_n,edges_ind,fnormal_ind,&
&                       lengU0,nbrNodes,nbrElem,nbrFront,nbrInt,only_mesh)
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
    INTEGER(ki), ALLOCATABLE    :: nbr_nodes_per_elem(:)
    
    REAL(kr), ALLOCATABLE :: Source(:)
    REAL(kr), ALLOCATABLE :: dt(:)
    
    REAL(kr), ALLOCATABLE :: edges(:,:)       
    REAL(kr), ALLOCATABLE :: fnormal(:,:)     
    REAL(kr), ALLOCATABLE :: geom_data(:,:)     
    REAL(kr), ALLOCATABLE :: cell_data_n(:,:) 
    INTEGER(ki), ALLOCATABLE :: edges_ind(:,:)  
    INTEGER(ki), ALLOCATABLE :: fnormal_ind(:,:)
    INTEGER(ki), ALLOCATABLE :: geom_data_ind(:,:) 
    
    INTEGER(ki), INTENT(IN) :: lengU0,nbrNodes,nbrElem,nbrFront,nbrInt,only_mesh
    
    WRITE(*,*) "Allocating the memory..."
    
    ! Requested in read_gmsh
    ALLOCATE(node(1:nbrNodes,1:2))
    ALLOCATE(front(1:nbrFront,1:4))
    ALLOCATE(elem(1:nbrElem,1:5))
    ALLOCATE(nbr_nodes_per_elem(1:nbrElem))
    ALLOCATE(U0(1:nbvar*nbrElem))
    ALLOCATE(depth(1:nbrElem))
    ALLOCATE(BoundCond(1:nbrFront,1:3))
        
    IF (only_mesh.NE.1) THEN
    ! Requested in runge_kutta
      ALLOCATE(dt(nbvar*nbrElem))
      ALLOCATE(Source(nbvar*nbrElem))
    
    ! Requested in get_normal_to_cell
      ALLOCATE(geom_data(1:nbrElem,1:4))
      ALLOCATE(geom_data_ind(1:nbrElem,1:4))
      ALLOCATE(fnormal(1:nbrFront,1:5))
      ALLOCATE(fnormal_ind(1:nbrFront,1:4))
      ALLOCATE(cell_data_n(1:nbrElem,1:8))
      
      ! These 2 are allocated in get_normal_to_cell
      ! ALLOCATE(edges(1:nbrInt,1:6))
      ! ALLOCATE(edges_ind(1:nbrInt,1:2))
    ENDIF
END SUBROUTINE mem_allocate

!##########################################################
! SUBROUTINE mem_deallocate
! 
!##########################################################
SUBROUTINE mem_deallocate(node,front,elem,nbr_nodes_per_elem,U0,depth,BoundCond,dt,Source,&
&                       edges,fnormal,geom_data,geom_data_ind,cell_data_n,edges_ind,fnormal_ind)
    USE module_shallow, only : kr,ki
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
    INTEGER(ki), ALLOCATABLE    :: nbr_nodes_per_elem(:)
    
    REAL(kr), ALLOCATABLE :: Source(:)
    REAL(kr), ALLOCATABLE :: dt(:)
    
    REAL(kr), ALLOCATABLE :: edges(:,:)       
    REAL(kr), ALLOCATABLE :: fnormal(:,:)     
    REAL(kr), ALLOCATABLE :: geom_data(:,:)   
    REAL(kr), ALLOCATABLE :: cell_data_n(:,:) 
    INTEGER(ki), ALLOCATABLE :: edges_ind(:,:)  
    INTEGER(ki), ALLOCATABLE :: fnormal_ind(:,:)
    INTEGER(ki), ALLOCATABLE :: geom_data_ind(:,:) 
        
    WRITE(*,*) "Deallocating the memory..."
    
    IF (ALLOCATED(node))  DEALLOCATE(node)
    IF (ALLOCATED(front)) DEALLOCATE(front)
    IF (ALLOCATED(elem))  DEALLOCATE(elem)
    IF (ALLOCATED(nbr_nodes_per_elem))  DEALLOCATE(nbr_nodes_per_elem)
    IF (ALLOCATED(U0))    DEALLOCATE(U0)
    IF (ALLOCATED(depth)) DEALLOCATE(depth)
    IF (ALLOCATED(BoundCond)) DEALLOCATE(BoundCond)
    
    IF (ALLOCATED(Source)) DEALLOCATE(Source)
    IF (ALLOCATED(dt))     DEALLOCATE(dt)
    
    IF (ALLOCATED(geom_data))   DEALLOCATE(geom_data)
    IF (ALLOCATED(geom_data_ind))   DEALLOCATE(geom_data_ind)
    IF (ALLOCATED(fnormal))     DEALLOCATE(fnormal)
    IF (ALLOCATED(fnormal_ind)) DEALLOCATE(fnormal_ind)
    IF (ALLOCATED(cell_data_n)) DEALLOCATE(cell_data_n)
    IF (ALLOCATED(edges))       DEALLOCATE(edges)
    IF (ALLOCATED(edges_ind))   DEALLOCATE(edges_ind)
    
END SUBROUTINE mem_deallocate