!##########################################################
! Main program for the 2D shallow water equations
!##########################################################
PROGRAM MAIN
    USE module_shallow
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
    
    ! Read the mesh and the initial solution
    CALL read_gmsh(ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the reading of the mesh"
      GOTO 200
    ELSE
      ! Compute the required geometrical data for the finite volume method
      CALL get_normal_to_cell

      IF (restart) THEN
         WRITE(*,*) "Starting from the previous solution ",trim(file_restart)
         CALL read_solution(ok)
         IF (ok == 0) THEN
            WRITE(*,*) "The program hasn't started because of a problem during the reading "
            WRITE(*,*) "of the former solution"
            GOTO 200
         ENDIF
      ENDIF
      
      ! Enter the time loop
      CALL runge_kutta

    ENDIF

    IF (ALLOCATED(node)) DEALLOCATE(node)
    IF (ALLOCATED(U0)) DEALLOCATE(U0)
    IF (ALLOCATED(geom_data)) DEALLOCATE(geom_data)
    IF (ALLOCATED(fnormal)) DEALLOCATE(fnormal)
    IF (ALLOCATED(cell_data_n)) DEALLOCATE(cell_data_n)
    IF (ALLOCATED(front)) DEALLOCATE(front)
    IF (ALLOCATED(elem)) DEALLOCATE(elem)
    IF (ALLOCATED(dt)) DEALLOCATE(dt)
    IF (ALLOCATED(dt)) DEALLOCATE(depth)
    IF (ALLOCATED(fnormal_ind)) DEALLOCATE(fnormal_ind)
    IF (ALLOCATED(BoundCond)) DEALLOCATE(BoundCond)
    IF (ALLOCATED(Source)) DEALLOCATE(Source)

200 CONTINUE
    WRITE(*,*) "End of the simulation"

END PROGRAM
