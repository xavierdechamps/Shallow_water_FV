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
    
    CHARACTER(LEN=40) :: mesh_file, file_restart, file_gmsh, file_dat
    
    INTEGER(ki) :: nbrNodes, nbrElem, nbrFront, nbrInt, nbrBC
    INTEGER(ki) :: nTime, shownTime, savenTime
    INTEGER(ki) :: time_begin, time_end
    INTEGER(ki), DIMENSION(2,5) :: CLTable
    
    LOGICAL   :: restart,steady
    
    REAL(kr)  :: CFL, ggrav, deltaTfixed, eps, manning
    
END MODULE module_shallow
!##########################################################
! MODULE signal_handler
!  Goal: contains all the interface to the routines that
!        manage the interception of stop signals
!##########################################################
MODULE signal_handler
!-----------------------------------------------------------------------
  USE, INTRINSIC :: iso_c_binding, only: C_INT, C_CHAR

  IMPLICIT NONE (type, external)

  INTERFACE

  INTEGER(C_INT) FUNCTION watchsignalname(signame, response) bind(C)
    import C_INT, C_CHAR
    CHARACTER(kind=c_char), INTENT(IN) :: signame(*)
    INTEGER(C_INT), INTENT(IN) :: response
  END FUNCTION watchsignalname

  INTEGER(C_INT) FUNCTION watchsignal(sig) bind(C)
    import C_INT
    INTEGER(C_INT), INTENT(IN) :: sig
  END FUNCTION watchsignal

  INTEGER(C_INT) FUNCTION getlastsignal() bind(C)
    import C_INT
  END FUNCTION getlastsignal

  END INTERFACE
!-----------------------------------------------------------------------
END MODULE signal_handler