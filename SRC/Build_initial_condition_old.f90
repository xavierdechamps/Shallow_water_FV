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
    CALL browse_gmsh(mesh_file,length_names,nbrNodes,nbrElem,nbrFront,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the browsing of the mesh"
      GOTO 200
    endif
    
    ! Allocate the memory for the arrays
    CALL mem_allocate(node,front,elem,U0,depth,BoundCond,dt,Source,&
&                     edges,fnormal,geom_data,cell_data_n,edges_ind,fnormal_ind,&
&                     nbvar*nbrElem,nbrNodes,nbrElem,nbrFront,nbrInt)

    ! Read the mesh and the initial solution / boundary conditions
    CALL read_gmsh(U0,nbvar*nbrElem,mesh_file,length_names,node,elem,front,depth,BoundCond,nbrNodes,nbrElem,nbrFront,ok)
    
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the reading of the mesh"
      GOTO 200
    ELSE
!      CALL get_normal_to_cell()
      CALL write_initial_condition_gmsh()
    ENDIF
    
    ! Deallocate the memory for the arrays
    CALL mem_deallocate(node,front,elem,U0,depth,BoundCond,dt,Source,&
&                       edges,fnormal,geom_data,cell_data_n,edges_ind,fnormal_ind)
    
200 CONTINUE
    WRITE(*,*) "End of the program"
    
END PROGRAM

! ******************************************************************************
SUBROUTINE write_initial_condition_gmsh()
    use module_shallow
    implicit none

    integer(ki) :: ierr, i, test_case, wall_type, edge_id
    real(kr) :: h, hleft,hright, u, v, h_inlet, u_inlet, v_inlet, h_outlet
    real(kr) :: xMax,xMin,xCenter,circleRadius, circleShiftY, x1, x2, xe
    real(kr) :: slope
    real(kr) :: arrayout(1:nbrElem)
    
!   test_case = 1:   rectilinear channel flow
!             = 2:   flow under bridge with cylinders
!             = 3-4: dam break
!             = 5:   cylindrical dam break
!             = 6:   oblique hydraulic jump
!             = 7:   channel with bump
    test_case = 7
        
    open(unit=10,file=file_gmsh,status="replace",iostat=ierr,form='formatted')
    write(10,'(T1,A11)') "$MeshFormat"
    write(10,'(T1,A7)') "2.2 0 8"
    write(10,'(T1,A14)') "$EndMeshFormat"
    write(10,'(T1,A6)') "$Nodes"
    write(10,'(T1,I9)') nbrNodes
    write(10,'(T1,I9,2ES24.16E2,F4.1)') (i, node(i,:),0.,i=1,nbrNodes)
    write(10,'(T1,A9)') "$EndNodes"
    write(10,'(T1,A9)') "$Elements"
    write(10,'(T1,I9)') nbrElem+nbrFront
    write(10,'(T1,I9,2I2,I9,I2,2I9)') (i,1,2,front(i,3),1,front(i,1:2),i=1,nbrFront)    
    write(10,'(T1,I9,2I2,I9,I2,3I9)') (i+nbrFront,2,2,elem(i,4),1,elem(i,1:3),i=1,nbrElem)
    write(10,'(T1,A12)') "$EndElements"

    !************************************* INITIAL HEIGHT
    write(10,'(T1,A12)') "$ElementData"
    write(10,'(T1,A1)') "1"
    write(10,'(T1,A9)') '"Height"'
    write(10,'(T1,A1)') "1"
    write(10,'(T1,I9)') 0
    write(10,'(T1,A1)') "3"
    write(10,'(T1,I9)') 0
    write(10,'(T1,A1)') "1"
    write(10,'(T1,I9)') nbrElem
    
    arrayout = 0.d00
    
    DO i=1,nbrElem
      select case (test_case)
         case (1 : 2) 
         ! rectilinear flow in channel
         ! flow under bridge with cylinders
            arrayout(i) = 0.03048d00
            
         case (3 : 4) ! dam breaks
            IF (elem(i,4).eq.22) THEN
              arrayout(i) = 5.0d0
            ELSEIF (elem(i,4).eq.23) THEN
              arrayout(i) = 10.0d0
            ENDIF
            
         case (5) ! radial dam break
            IF (elem(i,4).eq.27) THEN
              arrayout(i) = 1.0d0
            ELSEIF (elem(i,4).eq.28) THEN
              arrayout(i) = 1.1d0
            ENDIF
            
         case (6) ! oblique hydraulic jump
            arrayout(i) = 1.0d0
            
         case (7) ! channel with bump
            
            xMax = maxval(node(:,1) )
            xMin = minval(node(:,1) )
            hleft = 4.543260901d0 ! normal depth for b1=40, Q=500, S0=0.002, n=0.0389
            ! hleft = 2.5160369d0 ! critical depth for b1=40, Q=500
            
            arrayout(i) = hleft
            
         case default
            write(*,*) "default test case"
            arrayout(i) = 0.0d0
      end select
         
    ENDDO
    
    write(10,'(T1,I9,ES24.16E2)') (i+nbrFront, arrayout(i),i=1,nbrElem)
    write(10,'(T1,A15)') "$EndElementData"
    
    !************************************* INITIAL VELOCITY
    write(10,'(T1,A12)') "$ElementData"
    write(10,'(T1,A1)') "1"
    write(10,'(T1,A10)') '"Velocity"'
    write(10,'(T1,A1)') "1"
    write(10,'(T1,I9)') 0
    write(10,'(T1,A1)') "3"
    write(10,'(T1,I9)') 0
    write(10,'(T1,A1)') "3"
    write(10,'(T1,I9)') nbrElem
    
    do i=1,nbrElem
    
       select case (test_case)
         case (1 ) 
         ! rectilinear flow in channel
            ! U = 2.195837725d0
            U = 8.4566d0
            V = 0.d0
            
         case (2)
         ! flow under bridge with cylinders
            U = 0.1d0
            V = 0.d0
            
         case (6)
         ! oblique hydraulic jump
            U = 9.d0
            V = 0.d0
            
         case (7)
         ! channel with bump
!            U = 0.18d0 / 0.33d0
!            V = 0.d0

         ! Channel with slope U = q / h
            U = 500.d0 / ( 40d0 * 4.543260901d0 )
            V = 0.d0
            
         case default
            U = 0.d0
            V = 0.d0
            
       end select
    
       write(10,'(T1,I9,2ES24.16E2,F4.1)') i+nbrFront, U, V, 0.
    end do
    
    write(10,'(T1,A15)') "$EndElementData"

    !************************************* Bathymetric depth
    write(10,'(T1,A12)') "$ElementData"
    write(10,'(T1,A1)') "1"
    write(10,'(T1,A9)') '"Depth"'
    write(10,'(T1,A1)') "1"
    write(10,'(T1,I9)') 0
    write(10,'(T1,A1)') "3"
    write(10,'(T1,I9)') 0
    write(10,'(T1,A1)') "1"
    write(10,'(T1,I9)') nbrElem
    
    arrayout = 0.d00
    
! Circular bump
    xMax = maxval(node(:,1) )
    xCenter = xMax *0.5d0
    circleRadius = xMax*0.02d0
    circleShiftY = -circleRadius * 0.5d0
    x1 = xCenter-sqrt( circleRadius*circleRadius - circleShiftY*circleShiftY )
    x2 = xCenter+sqrt( circleRadius*circleRadius - circleShiftY*circleShiftY )
    
! Linear slope
    xMax = maxval(node(:,1) )
    xMin = minval(node(:,1) )
    hleft = 0.0d0
!    slope = 0.05664 ! = Dy / Dx
    slope = 0.002d0 ! = Dy / Dx
            
    DO i=1,nbrElem
      
      xe =  (node(elem(i,1),1)+node(elem(i,2),1)+node(elem(i,3),1) ) /3.d0
      
      select case (test_case)
         
         ! case (1) ! supercritical symm. contraction
            ! arrayout(i) = hleft - slope*(xe-xMin)
         
         case (7) ! channel with bump
            
            ! center of element
!            if (xe.ge.8.d0 .and. xe.le.12.d0) then
!              arrayout(i) = -0.05d0*(xe-10.d0)*(xe-10.d0)
!            else
!              arrayout(i) = -0.2d0
!            endif
            
            ! Circular bump              
            ! if (xe.ge.x1 .and. xe.le.x2) then
              ! arrayout(i) = circleShiftY + sqrt( circleRadius*circleRadius - (xe-xCenter)**2 )
            ! endif
            
            ! Linear slope
            arrayout(i) = hleft - slope*(xe-xMin)
            
         case default
            arrayout(i) = 0.0d0
      end select
         
    ENDDO
    
    write(10,'(T1,I9,ES24.16E2)') (i+nbrFront, arrayout(i),i=1,nbrElem)
    write(10,'(T1,A15)') "$EndElementData"
    
    !************************************* Boundary Condition - Height 
    write(10,'(T1,A12)') "$ElementData"
    write(10,'(T1,A1)') "1"
    write(10,'(T1,A14)') '"Height_inlet"'
    write(10,'(T1,A1)') "1"
    write(10,'(T1,I9)') 0
    write(10,'(T1,A1)') "3"
    write(10,'(T1,I9)') 0
    write(10,'(T1,A1)') "1"
    write(10,'(T1,I9)') nbrFront
    
    do i=1,nbrFront
    
       IF (test_case.eq.1) THEN
          ! h_inlet = 0.03048d0
          h_inlet = 1.0d0
          IF ( front(i,3).eq.1 .or. front(i,3).eq.2 ) THEN
        ! Inlet / outlet
            write(10,'(T1,I9,2X,ES24.16E2)') i, h_inlet
          ELSE
            write(10,'(T1,I9,2X,ES24.16E2)') i, 0.
          ENDIF
          
       ELSE IF (test_case.eq.2) THEN
          h_inlet = 1.d0
          IF ( front(i,3).eq.24 .or. front(i,3).eq.25 ) THEN
        ! Inlet / outlet
            write(10,'(T1,I9,2X,ES24.16E2)') i, h_inlet
          ELSE
            write(10,'(T1,I9,2X,ES24.16E2)') i, 0.
          ENDIF
          
       ELSE IF (test_case.eq.6) THEN
          h_inlet = 1.d0
          IF ( front(i,3).eq.21 .or. front(i,3).eq.22 ) THEN
        ! Inlet / outlet
            write(10,'(T1,I9,2X,ES24.16E2)') i, h_inlet
          ELSE
            write(10,'(T1,I9,2X,ES24.16E2)') i, 0.
          ENDIF
          
       ELSE IF (test_case.eq.7) THEN
          ! Circular bump
          ! h_inlet = 2.d0
          ! IF ( front(i,3).eq.2 .or. front(i,3).eq.3 ) THEN
        !! Inlet / outlet
            ! write(10,'(T1,I9,2X,ES24.16E2)') i, h_inlet
          ! ELSE
            ! write(10,'(T1,I9,2X,ES24.16E2)') i, 0.
          ! ENDIF
          
          ! Linear slope
          h_inlet = 4.543260901d0
          h_outlet = 2.5160369d0
          
!          h_inlet = 0.18d0
!          h_outlet = 0.33d0
          IF ( front(i,3).eq.1 ) THEN
        ! Inlet 
            write(10,'(T1,I9,2X,ES24.16E2)') i, h_inlet
          ELSE IF ( front(i,3).eq.2 ) THEN
        ! Outlet 
            write(10,'(T1,I9,2X,ES24.16E2)') i, h_outlet
          ELSE
            write(10,'(T1,I9,2X,ES24.16E2)') i, 0.
          ENDIF
          
       ELSE
          h_inlet = 0.d0
          write(10,'(T1,I9,2X,ES24.16E2)') i, h_inlet
       ENDIF
    
    end do
    
    write(10,'(T1,A15)') "$EndElementData"
    
    !************************************* Boundary Condition - Velocity 
    write(10,'(T1,A12)') "$ElementData"
    write(10,'(T1,A1)') "1"
    write(10,'(T1,A16)') '"Velocity_inlet"'
    write(10,'(T1,A1)') "1"
    write(10,'(T1,I9)') 0
    write(10,'(T1,A1)') "3"
    write(10,'(T1,I9)') 0
    write(10,'(T1,A1)') "3"
    write(10,'(T1,I9)') nbrFront
    
    do i=1,nbrFront
    
        IF (test_case.eq.1) THEN
          ! u_inlet = 2.195837725d0
          u_inlet = 8.4566d0
          v_inlet = 0.d0
          IF ( front(i,3).eq.1 ) THEN
            ! Inlet 
              write(10,'(T1,I9,1X,2ES24.16E2,F4.1)') i, u_inlet, v_inlet, 0.
           ELSE
              write(10,'(T1,I9,1X,3F4.1)') i, 0. , 0. , 0.
           ENDIF
           
        ELSE IF (test_case.eq.2) THEN
          u_inlet = 0.1d0
          v_inlet = 0.d0
          IF ( front(i,3).eq.24 ) THEN
            ! Inlet 
              write(10,'(T1,I9,1X,2ES24.16E2,F4.1)') i, u_inlet, v_inlet, 0.
           ELSE
              write(10,'(T1,I9,1X,3F4.1)') i, 0. , 0. , 0.
           ENDIF
           
        ELSE IF (test_case.eq.6) THEN
          u_inlet = 9.d0
          v_inlet = 0.d0
          IF ( front(i,3).eq.22 ) THEN
            ! Inlet 
              write(10,'(T1,I9,1X,2ES24.16E2,F4.1)') i, u_inlet, v_inlet, 0.
          ELSE
              write(10,'(T1,I9,1X,3F4.1)') i, 0. , 0. , 0.
          ENDIF
           
        ELSE IF (test_case.eq.7) THEN
          ! Circular bump
          
          ! u_inlet = 1.5d0
          ! v_inlet = 0.d0
          ! IF ( front(i,3).eq.3 ) THEN
         !!   Inlet 
              ! write(10,'(T1,I9,1X,2ES24.16E2,F4.1)') i, u_inlet, v_inlet, 0.
          ! ELSE
              ! write(10,'(T1,I9,1X,3F4.1)') i, 0. , 0. , 0.
          ! ENDIF
          
          
          ! Channel with slope U = q / h
          u_inlet = 500.d0 / ( 40d0 * 4.543260901d0 )
          v_inlet = 0.d0
          
!          u_inlet = 0.18d0 / 0.33d0
!          v_inlet = 0.d0
          IF ( front(i,3).eq.1 ) THEN
        ! Inlet 
            write(10,'(T1,I9,1X,2ES24.16E2,F4.1)') i, u_inlet, v_inlet, 0.
          ELSE
            write(10,'(T1,I9,1X,3F4.1)') i, 0. , 0. , 0.
          ENDIF
          
        ELSE
          u_inlet = 0.d0
          v_inlet = 0.d0
          write(10,'(T1,I9,1X,2ES24.16E2,F4.1)') i, u_inlet, v_inlet, 0.
        ENDIF
    
    end do
    
    write(10,'(T1,A15)') "$EndElementData"
    
    !*************************************
    close(unit=10)

end subroutine write_initial_condition_gmsh

! *******************************************************************
subroutine sampletime(counter)
    use module_shallow, only : kr,ki
    implicit none
          
    ! variables passed through header
    integer(ki) ::counter
      
    ! variables declared locally
    integer(ki) ::rate, contmax
 
    ! Determine CPU time
    call system_clock(counter, rate, contmax )
!-----------------------------------------------------------------------
end subroutine sampletime
!-----------------------------------------------------------------------

! *******************************************************************
SUBROUTINE time_display
    use module_shallow
    implicit none

    real(kr) :: time
    integer(ki) :: job, cont3
    integer(ki) :: rate, contmax, itime

    call system_clock(cont3, rate, contmax )
    if (time_end .ge. time_begin) then
       itime=time_end-time_begin
    else
       itime=(contmax - time_begin) + (time_end + 1)
    endif

    time = dfloat(itime) / dfloat(rate)

    write(*,'(a,f10.4)') " Time needed (s) :            ",time

end subroutine time_display