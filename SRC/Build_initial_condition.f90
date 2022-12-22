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
&                     edges,fnormal,geom_data,geom_data_ind,cell_data_n,edges_ind,fnormal_ind,&
&                     shock_indicator,nbvar*nbrElem,nbrNodes,nbrElem,nbrFront,nbrInt,0)

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
&                       edges,fnormal,geom_data,geom_data_ind,cell_data_n,edges_ind,fnormal_ind,&
&                       shock_indicator)
    
200 CONTINUE
    WRITE(*,*) "End of the program"
    
END PROGRAM

! ******************************************************************************
SUBROUTINE write_initial_condition_gmsh()
    USE module_shallow
    IMPLICIT NONE

    INTEGER(ki) :: ierr, i, numdigits
    CHARACTER(LEN=2) :: numdig
    CHARACTER(LEN=9) :: formatreal
    REAL(kr) :: height_init(nbrElem),depth_init(nbrElem)
    REAL(kr) :: velocity_init(nbrElem,2)
    REAL(kr) :: height_BC(nbrFront),velocity_BC(nbrFront,2)
    REAL(kr) :: Q, Fr_i, h_i, u_i, g, yMax, xMin, xMax, hleft, slope, xe
    
    CALL get_number_digits_integer(nbrNodes,numdigits)
    WRITE(numdig,'(A,I1)') 'I',numdigits+1

    formatreal = 'ES24.15E3'

    !************************************* INITIAL HEIGHT    
    ! Q0 = 2.5 or 5. L/s
    Q    = 5.0d0 / 1000.0d0
    Fr_i = 4.d0
    g    = 9.81d0
    yMax = maxval(node(:,2) )
    h_i  = ( Q*Q / (yMax*yMax*Fr_i*Fr_i*g) )**(1.d0/3.d0)
    
    height_init = h_i
        
    !************************************* INITIAL VELOCITY
    u_i = Q / ( yMax * h_i )
    velocity_init(:,1) = u_i
    velocity_init(:,2) = 0.0d00

    !************************************* Bathymetric depth
    ! Linear slope
    xMax = maxval(node(:,1) )
    xMin = minval(node(:,1) )
    hleft = 0.0d0
    slope = 0.015d0 ! = Dy / Dx
            
    DO i=1,nbrElem
      xe =  (node(elem(i,1),1)+node(elem(i,2),1)+node(elem(i,3),1) ) /3.d0
      
      ! Linear slope
      depth_init(i) = hleft - slope*(xe-xMin)
    ENDDO
    
    !************************************* Boundary Condition - Height 
    height_BC = zero
    DO i=1,nbrFront
        IF ( front(i,3).eq.1 .or.front(i,3).eq.2) THEN
        ! Inlet 
            height_BC(i) = h_i
        ENDIF
    END DO
    
    !************************************* Boundary Condition - Velocity 
    velocity_BC = zero
    DO i=1,nbrFront
        IF ( front(i,3).eq.1) THEN
        ! Inlet 
            velocity_BC(i,1) = u_i
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