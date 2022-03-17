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
    REAL(kr) :: hleft,hright, xMax,xMin, xe, ye
    REAL(kr) :: hi, he, B0, q, alpha, a, b, c, d, xm
    REAL(kr) :: slope1,slope2,slope
    REAL(kr) :: height_init(nbrElem),depth_init(nbrElem)
    REAL(kr) :: velocity_init(nbrElem,2)
    REAL(kr) :: height_BC(nbrFront),velocity_BC(nbrFront,2)
        
    CALL get_number_digits_integer(nbrNodes,numdigits)
    WRITE(numdig,'(A,I1)') 'I',numdigits+1

    formatreal = 'ES24.15E3'

    !************************************* INITIAL HEIGHT
    hi = 0.7d0
    B0 = 2.0d00
    q  = 4.0d00
    
    height_init = 0.0d00
    DO i=1,nbrElem
      height_init(i) = hi
    ENDDO
    
    !************************************* INITIAL VELOCITY    
    velocity_init = 0.0d00
    DO i=1,nbrElem
! Uniform velocity profile
        velocity_init(i,1) = q / hi

! Y-coordinate of the cell center
!       IF (nbr_nodes_per_elem(i) .EQ. 3) THEN
!         ye = (node(elem(i,1),2)+node(elem(i,2),2)+node(elem(i,3),2) ) / 3.0d00
!       ELSE IF (nbr_nodes_per_elem(i) .EQ. 4) THEN
!         ye = (node(elem(i,1),2)+node(elem(i,2),2)+node(elem(i,3),2)+node(elem(i,4),2) ) * 0.25d00
!       END IF

! Quadratic velocity profile
!        velocity_init(i,1) = -6.d00 * q * ye * ( ye - B0 ) / ( hi * B0 * B0 )
        
! Fourth order velocity profile
!        alpha = -75.d00
!        a = 5.d00 * (B0*q/hi + alpha*B0*B0*B0/12.d00)/(B0**5)
!        b = -2.d00 * a *B0
!        c = 0.5d00 * alpha
!        d = a*B0**3 - 0.5d00*alpha*B0
!        velocity_init(i,1) = a*ye**4 + b*ye**3 + c*ye**2 + d*ye
    END DO
    
    !************************************* Bathymetric depth    
! Linear slope
    xMax = MAXVAL(node(:,1) )
    xMin = MINVAL(node(:,1) )
    hright = 2.0d0
    slope1 = 0.002d0   ! Slope in the first section of the channel
    slope2 = 0.0005d0  ! Slope in the second section of the channel
    
    depth_init = 0.0d00
    DO i=1,nbrElem
! X-coordinate of the center of the cell
       IF (nbr_nodes_per_elem(i) .EQ. 3) THEN
         xe = (node(elem(i,1),1)+node(elem(i,2),1)+node(elem(i,3),1) ) / 3.0d00
       ELSE IF (nbr_nodes_per_elem(i) .EQ. 4) THEN
         xe = (node(elem(i,1),1)+node(elem(i,2),1)+node(elem(i,3),1)+node(elem(i,4),1) ) * 0.25d00
       END IF
      
! X-coordinate where the slope changes from slope1 to slope2
      xm = xMax * 1.5d00
      
      IF (xe .LE. xm) THEN   ! Steep
! First section of the channel
        slope = slope1
        hleft = hright
        depth_init(i) = hleft - slope*(xe-xMin)
      ELSE                   ! Mild
! Second section of the channel
        slope = slope2
        hleft = hright - slope1*xm
        depth_init(i) = hleft - slope*(xe-xMin-xm)     
      ENDIF      
    ENDDO
    
    !************************************* Boundary Condition - Height
    height_BC = 0.0d00
    DO i=1,nbrFront
        IF ( front(i,3).eq.1 ) THEN ! Inlet 
          height_BC(i) = hi
        ELSE IF ( front(i,3).eq.2 ) THEN ! Outlet
          height_BC(i) = hi
        ENDIF
    END DO
        
    !************************************* Boundary Condition - Velocity
    velocity_BC = 0.0d00
    DO i=1,nbrFront
      IF ( front(i,3).eq.1 ) THEN
! Uniform velocity profile
        velocity_BC(i,1) = q / hi

! Y-coordinate of the edge center
!        ye = 0.5d00 * ( node( front(i,1) ,2) + node( front(i,2) ,2) )

! Quadratic velocity profile
!        velocity_BC(i,1) = -6.d00 * q * ye * ( ye - B0 ) / ( hi * B0 * B0 )
        
! Fourth order velocity profile
!        alpha = -75.d00
!        a = 5.d00 * (B0*q/hi + alpha*B0*B0*B0/12.d00)/(B0**5)
!        b = -2.d00 * a *B0
!        c = 0.5d00 * alpha
!        d = a*B0**3 - 0.5d00*alpha*B0
!        velocity_BC(i,1) = a*ye**4 + b*ye**3 + c*ye**2 + d*ye
      END IF
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