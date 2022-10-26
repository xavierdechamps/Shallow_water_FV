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
       
    INTEGER(ki), ALLOCATABLE :: geom_data_ind(:,:)    ! see read_gmsh
    INTEGER(ki), ALLOCATABLE :: front(:,:)            ! see read_gmsh
    INTEGER(ki), ALLOCATABLE :: elem(:,:)             ! see read_gmsh
    INTEGER(ki), ALLOCATABLE :: nbr_nodes_per_elem(:) ! see read_gmsh
    INTEGER(ki), ALLOCATABLE :: edges_ind(:,:)  ! see get_normal_to_cell
    INTEGER(ki), ALLOCATABLE :: fnormal_ind(:,:)! see get_normal_to_cell
    
    INTEGER(ki), PARAMETER :: length_names = 40
    CHARACTER(LEN=length_names) :: mesh_file, file_restart, file_gmsh, file_dat
    
    INTEGER(ki) :: nbrNodes, nbrElem, nbrTris, nbrQuads, nbrFront, nbrInt, nbrBC
    INTEGER(ki) :: nTime, shownTime, savenTime
    INTEGER(ki) :: time_begin, time_end, steady, restart
    INTEGER(ki), DIMENSION(2,5) :: CLTable
    
    
    REAL(kr)  :: CFL, ggrav, deltaTfixed, eps, manning_b, manning_w
    
END MODULE module_shallow

MODULE module_mem_allocate
  INTERFACE 
  
! -----------------------
! Subroutine mem_allocate
! -----------------------
  SUBROUTINE mem_allocate(node,front,elem,nbr_nodes_per_elem,U0,depth,BoundCond,dt,Source,&
&                         edges,fnormal,geom_data,geom_data_ind,cell_data_n,edges_ind,fnormal_ind,&
&                         lengU0,nbrNodes,nbrElem,nbrFront,nbrInt,only_mesh)
    USE module_shallow, only : kr,ki,nbvar
    IMPLICIT NONE

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
        
  END SUBROUTINE mem_allocate
  
! -------------------------
! Subroutine mem_deallocate
! -------------------------
  SUBROUTINE mem_deallocate(node,front,elem,nbr_nodes_per_elem,U0,depth,BoundCond,dt,Source,&
&                       edges,fnormal,geom_data,geom_data_ind,cell_data_n,edges_ind,fnormal_ind)
    USE module_shallow, only : kr,ki
    IMPLICIT NONE

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
        
  END SUBROUTINE mem_deallocate
  
  END INTERFACE 
  
END MODULE module_mem_allocate

!##########################################################
! MODULE signal_handler
!  Goal: contains all the interface to the routines that
!        manage the interception of stop signals
!##########################################################
#ifdef WITHSIGWATCH
MODULE signal_handler
!-----------------------------------------------------------------------
!DEC$ IF .NOT. DEFINED( SIGNAL_NAMES )
!DEC$ DEFINE SIGNAL_NAMES

#ifdef UPPERCASE
#define Watch_signal_name WATCHSIGNALNAME
#define Watch_signal WATCHSIGNAL
#define Get_last_signal GETLASTSIGNAL
#endif

#ifdef UPPERCASE_
#define Watch_signal_name WATCHSIGNALNAME_
#define Watch_signal WATCHSIGNAL_
#define Get_last_signal GETLASTSIGNAL_
#endif

#ifdef lowercase_
#define Watch_signal_name watchsignalname_
#define Watch_signal watchsignal_
#define Get_last_signal getlastsignal_
#endif

#ifdef lowercase
#define Watch_signal_name watchsignalname
#define Watch_signal watchsignal
#define Get_last_signal getlastsignal
#endif

!DEC$ ENDIF

  USE, INTRINSIC :: iso_c_binding, only: C_INT, C_CHAR

    IMPLICIT NONE !(type, external)

  INTERFACE

    INTEGER(C_INT) FUNCTION Watch_signal_name(signame, response) bind(C)
    import C_INT, C_CHAR
    CHARACTER(kind=c_char), INTENT(IN) :: signame(*)
    INTEGER(C_INT), INTENT(IN) :: response
    END FUNCTION Watch_signal_name

    INTEGER(C_INT) FUNCTION Watch_signal(sig) bind(C)
    import C_INT
    INTEGER(C_INT), INTENT(IN) :: sig
    END FUNCTION Watch_signal

    INTEGER(C_INT) FUNCTION Get_last_signal() bind(C)
    import C_INT
    END FUNCTION Get_last_signal

  END INTERFACE
!-----------------------------------------------------------------------
END MODULE signal_handler
#endif
