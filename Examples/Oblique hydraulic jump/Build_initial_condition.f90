PROGRAM BUILD_INITIAL_SOLUTION
    USE module_shallow
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
    
    CALL read_gmsh(ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the reading of the mesh"
      GOTO 200
    ELSE
!      CALL get_normal_to_cell()
      CALL write_initial_condition_gmsh()
    ENDIF
    
200 CONTINUE
    WRITE(*,*) "End of the program"
    
END PROGRAM

! ******************************************************************************
SUBROUTINE write_initial_condition_gmsh()
    use module_shallow
    implicit none

    integer(ki) :: ierr, i, wall_type, edge_id
    real(kr) :: h, hleft,hright, u, v, h_inlet, u_inlet, v_inlet, h_outlet
    real(kr) :: xMax,xMin,xCenter,circleRadius, circleShiftY, x1, x2, xe
    real(kr) :: slope
    real(kr) :: arrayout(1:nbrElem)
            
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
        hleft = 1.0d0
        arrayout(i) = hleft         
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
        U = 9.0d0
        V = 0.d0    
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
        h_inlet = 1.0d0
        IF ( front(i,3).eq.22 .or.front(i,3).eq.21) THEN
        ! Inlet 
            write(10,'(T1,I9,2X,ES24.16E2)') i, h_inlet
        ELSE
            write(10,'(T1,I9,2X,ES24.16E2)') i, 0.
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
        u_inlet = 9.0d0
        v_inlet = 0.d0
        IF ( front(i,3).eq.22 ) THEN
        ! Inlet 
            write(10,'(T1,I9,1X,2ES24.16E2,F4.1)') i, u_inlet, v_inlet, 0.
        ELSE
            write(10,'(T1,I9,1X,3F4.1)') i, 0. , 0. , 0.
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