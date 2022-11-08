!##########################################################
! Main program for the 2D shallow water equations
!##########################################################
PROGRAM MAIN
    USE module_shallow
    USE module_mem_allocate
    IMPLICIT NONE

    ! Local parameters
    INTEGER(ki) :: ok, i, j
    REAL(kr)    :: h, incidence, normeU, u, v, hgauche, hdroite
    REAL(kr)    :: hinit, uinit, vinit
    CHARACTER(LEN=100)::param_file 
    
    ! Get the number of arguments
    IF(COMMAND_ARGUMENT_COUNT().NE.1)THEN
      WRITE(*,*) "Incorrect number of arguments. Please launch the program as"
      WRITE(*,*) "shallow name_of_parameter_file"
      GOTO 200
    ENDIF 
    CALL get_command_argument(1,param_file)
    
    ! Get the flow parameters
    CALL read_parameters(param_file,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "Problem reading the parameters"
      GOTO 200
    ENDIF
    eps = 1.E-12
    
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
&                  depth,BoundCond,nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,0,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the reading of the mesh"
      GOTO 200
    ELSE
      ! Compute the required geometrical data for the finite volume method
      CALL get_normal_to_cell()
      IF (restart.EQ.1) THEN
         WRITE(*,*) "Starting from the previous solution ",trim(file_restart)
         CALL read_solution(U0,nbvar*nbrElem,file_restart,length_names,ok)
         IF (ok == 0) THEN
            WRITE(*,*) "The program hasn't started because of a problem during the reading "
            WRITE(*,*) "of the former solution"
            GOTO 200
         ENDIF
      ENDIF
      
      ! Enter the time loop
      CALL runge_kutta()

    ENDIF
    
    ! Deallocate the memory for the arrays
    CALL mem_deallocate(node,front,elem,nbr_nodes_per_elem,U0,depth,BoundCond,dt,Source,&
&                       edges,fnormal,geom_data,geom_data_ind,cell_data_n,edges_ind,fnormal_ind,&
&                       shock_indicator)
    
200 CONTINUE
    WRITE(*,*) "End of the simulation"

END PROGRAM
