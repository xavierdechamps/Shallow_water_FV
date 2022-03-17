PROGRAM BUILD_INITIAL_SOLUTION
    USE module_shallow
    USE module_mem_allocate
    IMPLICIT NONE
    
    INTEGER(ki) :: ok
    CHARACTER(LEN=100)::param_file
    
    ! Get the number of arguments
    IF(COMMAND_ARGUMENT_COUNT().NE.2)THEN
      WRITE(*,*) "Incorrect number of arguments. Please launch the program as"
      WRITE(*,*) "build_initial_condition   parameter_file   initial_condition_file.msh"
      GOTO 200
    ENDIF 
    CALL get_command_argument(1,param_file)
    CALL get_command_argument(2,file_gmsh)
    
    ! ------- Parameters of the flow ----------
    CALL read_parameters(param_file,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "Problem reading the parameters"
      GOTO 200
    ENDIF
    CALL get_command_argument(2,file_gmsh)
    
    ! Browse the mesh to get the size of the arrays
    CALL browse_gmsh(mesh_file,length_names,nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the browsing of the mesh"
      GOTO 200
    endif
    
    ! Allocate the memory for the arrays
    CALL mem_allocate(node,front,elem,nbr_nodes_per_elem,U0,depth,BoundCond,dt,Source,&
&                     edges,fnormal,geom_data,cell_data_n,edges_ind,fnormal_ind,&
&                     nbvar*nbrElem,nbrNodes,nbrElem,nbrFront,nbrInt,0)

    ! Read the mesh and the initial solution / boundary conditions
    CALL read_gmsh(U0,nbvar*nbrElem,mesh_file,length_names,node,elem,nbr_nodes_per_elem,front,&
&                  depth,BoundCond,nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,1,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the reading of the mesh"
      GOTO 200
    ELSE
      CALL write_initial_condition_gmsh()
    ENDIF
    
    ! Deallocate the memory for the arrays
    CALL mem_deallocate(node,front,elem,nbr_nodes_per_elem,U0,depth,BoundCond,dt,Source,&
&                       edges,fnormal,geom_data,cell_data_n,edges_ind,fnormal_ind)
    
200 CONTINUE
    WRITE(*,*) "End of the program"
    
END PROGRAM

! ******************************************************************************
SUBROUTINE write_initial_condition_gmsh()
    USE module_shallow
    IMPLICIT NONE

    INTEGER(ki) :: ierr, i, wall_type, edge_id, numdigits
    CHARACTER(LEN=2) :: numdig
    CHARACTER(LEN=9) :: formatreal
    REAL(kr) :: hi, Froude, g, ui
    REAL(kr) :: xMax,xMin,xCenter,circleRadius, circleShiftY, x1, x2, xe
    REAL(kr) :: slope
    REAL(kr) :: height_init(nbrElem),depth_init(nbrElem)
    REAL(kr) :: velocity_init(nbrElem,2)
    REAL(kr) :: height_BC(nbrFront),velocity_BC(nbrFront,2)
            
    CALL get_number_digits_integer(nbrNodes,numdigits)
    WRITE(numdig,'(A,I1)') 'I',numdigits+1

    formatreal = 'ES24.15E3'

    !************************************* INITIAL HEIGHT
    hi     = 1.0d00
    Froude = 2.7d00
    g      = 9.81d00
    ui     = Froude * sqrt(g * hi)
    
    height_init = 0.0d00
    DO i=1,nbrElem
      height_init(i) = hi         
    ENDDO
    
    !************************************* INITIAL VELOCITY
    velocity_init = 0.0d00
    DO i=1,nbrElem
      velocity_init(i,1) = ui
    END DO
    
    !************************************* Bathymetric depth
    depth_init = 0.0d00
    
    !************************************* Boundary Condition - Height 
    height_BC = 0.0d00
    DO i=1,nbrFront
      IF ( front(i,3).eq.1 .or.front(i,3).eq.2) THEN ! Inlet 
        height_BC(i) = hi
      ENDIF
    END DO
    
    !************************************* Boundary Condition - Velocity 
    velocity_BC = 0.0d00
    DO i=1,nbrFront
      IF ( front(i,3).eq.1 ) THEN ! Inlet 
        velocity_BC(i,1) = ui
      ENDIF
    END DO
    
    !*************************************
    
    CALL write_gmsh_initial_solution(height_init, velocity_init, depth_init, &
&                                    height_BC, velocity_BC)

END SUBROUTINE write_initial_condition_gmsh

! *******************************************************************
SUBROUTINE sampletime(counter)
    USE module_shallow, ONLY : kr,ki
    IMPLICIT NONE
          
    ! variables passed through header
    INTEGER(ki) ::counter
      
    ! variables declared locally
    INTEGER(ki) ::rate, contmax
 
    ! Determine CPU time
    CALL SYSTEM_CLOCK(counter, rate, contmax )
!-----------------------------------------------------------------------
END SUBROUTINE sampletime
!-----------------------------------------------------------------------

! *******************************************************************
SUBROUTINE time_display
    USE module_shallow
    IMPLICIT NONE

    REAL(kr) :: time
    INTEGER(ki) :: job, cont3
    INTEGER(ki) :: rate, contmax, itime

    CALL SYSTEM_CLOCK(cont3, rate, contmax )
    IF (time_end .ge. time_begin) THEN
       itime=time_end-time_begin
    ELSE
       itime=(contmax - time_begin) + (time_end + 1)
    ENDIF

    time = DFLOAT(itime) / DFLOAT(rate)

    WRITE(*,'(a,f10.4)') " Time needed (s) :            ",time

END SUBROUTINE time_display