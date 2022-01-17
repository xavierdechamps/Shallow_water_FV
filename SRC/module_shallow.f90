!##########################################################
! MODULE module_shallow
!  Goal: contains all the data used for the simulation
!##########################################################
MODULE module_shallow
    IMPLICIT NONE
    SAVE

    INTEGER, PARAMETER    :: kr = 8
    INTEGER, PARAMETER    :: ki = 4
    
    INTEGER(ki), PARAMETER:: nbvar = 3
    REAL(kr), PARAMETER   :: zero = 0.0d00
    
    REAL(kr), ALLOCATABLE :: node(:,:) ! Node coordinates
    REAL(kr), ALLOCATABLE :: U0(:)     ! Solution ( ... h hu hv ...)
    REAL(kr), ALLOCATABLE :: Source(:) ! source term
    REAL(kr), ALLOCATABLE :: dt(:)     ! time increment in each cell
    REAL(kr), ALLOCATABLE :: depth(:)  ! depth in each cell
    REAL(kr), ALLOCATABLE :: edges(:,:)       ! see get_normal_to_cell
    REAL(kr), ALLOCATABLE :: fnormal(:,:)     ! see get_normal_to_cell
    REAL(kr), ALLOCATABLE :: geom_data(:,:)   ! see get_normal_to_cell
    REAL(kr), ALLOCATABLE :: cell_data_n(:,:) ! see get_normal_to_cell
    REAL(kr), ALLOCATABLE :: BoundCond(:,:)   ! Imposed boundary conditions ( boundary_edge , hBC , uBC , vBC)
    
    INTEGER(ki), ALLOCATABLE :: front(:,:)      ! see read_gmsh
    INTEGER(ki), ALLOCATABLE :: elem(:,:)       ! see read_gmsh
    INTEGER(ki), ALLOCATABLE :: edges_ind(:,:)  ! see get_normal_to_cell
    INTEGER(ki), ALLOCATABLE :: fnormal_ind(:,:)! see get_normal_to_cell
    
    INTEGER(ki), PARAMETER :: length_names = 40
    CHARACTER(LEN=length_names) :: mesh_file, file_restart, file_gmsh, file_dat
    
    INTEGER(ki) :: nbrNodes, nbrElem, nbrFront, nbrInt, nbrBC
    INTEGER(ki) :: nTime, shownTime, savenTime
    INTEGER(ki) :: time_begin, time_end
    INTEGER(ki), DIMENSION(2,5) :: CLTable
    
    LOGICAL   :: restart,steady
    
    REAL(kr)  :: CFL, ggrav, deltaTfixed, eps, manning
    
END MODULE module_shallow

MODULE module_mem_allocate
  interface 
  
  SUBROUTINE mem_allocate(node,front,elem,U0,depth,BoundCond,lengU0,nbrNodes,nbrElem,nbrFront)
    USE module_shallow, only : kr,ki,nbvar
    IMPLICIT NONE

    REAL(kr), ALLOCATABLE    :: U0(:)
    REAL(kr), ALLOCATABLE    :: node(:,:)
    REAL(kr), ALLOCATABLE    :: depth(:)
    REAL(kr), ALLOCATABLE    :: BoundCond(:,:)
    INTEGER(ki), ALLOCATABLE    :: elem(:,:)
    INTEGER(ki), ALLOCATABLE    :: front(:,:)
    
    INTEGER(ki), INTENT(IN) :: lengU0,nbrNodes,nbrElem,nbrFront
        
  END SUBROUTINE mem_allocate
  
  end interface 
  
END MODULE module_mem_allocate
