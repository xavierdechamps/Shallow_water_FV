PROGRAM GET_RELATIVE_ERROR_THEORY
    USE module_shallow
    USE module_mem_allocate
    IMPLICIT NONE
    
    INTEGER(ki) :: ok
    CHARACTER(LEN=100)::param_file
    
    ! Get the number of arguments
    IF(COMMAND_ARGUMENT_COUNT().NE.1)THEN
      WRITE(*,*) "Incorrect number of arguments. Please launch the program as"
      WRITE(*,*) "build_initial_condition   parameter_file"
      GOTO 200
    ENDIF 
    CALL get_command_argument(1,param_file)
    
    ! ------- Parameters of the flow ----------
    CALL read_parameters(param_file,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "Problem reading the parameters"
      GOTO 200
    ENDIF
    
    ! Browse the mesh to get the size of the arrays
    CALL browse_gmsh(mesh_file,length_names,nbrNodes,nbrElem,nbrFront,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the browsing of the mesh"
      GOTO 200
    endif
    
    ! Allocate the memory for the arrays
    CALL mem_allocate(node,front,elem,U0,depth,BoundCond,nbvar*nbrElem,nbrNodes,nbrElem,nbrFront)
    
    ! Read the mesh and the initial solution / boundary conditions
    CALL read_gmsh(U0,nbvar*nbrElem,mesh_file,length_names,node,elem,front,depth,BoundCond,nbrNodes,nbrElem,nbrFront,ok)
    
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the reading of the mesh"
      GOTO 200
    ELSE
      CALL read_solution(U0,nbvar*nbrElem,file_restart,length_names,ok)
      IF (ok == 0) THEN
        WRITE(*,*) "The program hasn't started because of a problem during the reading of the mesh"
        GOTO 200
      ELSE
      
        CALL relative_error(ok)
        IF (ok == 0) THEN
          WRITE(*,*) "Error during the computation of the relative error"
          GOTO 200
        ENDIF
      ENDIF
      
    ENDIF
    
200 CONTINUE
    WRITE(*,*) "End of the program"
    
END PROGRAM

! ******************************************************************************
SUBROUTINE relative_error(ok)
    USE MODULE_SHALLOW
    IMPLICIT NONE
    
    integer(ki) :: ok
    
    INTEGER(ki) :: ierr, i, test_case
    REAL(kr)    :: h1,h2,href,u1,u2,uref,g,theta,beta,U
    REAL(kr)    :: x1,x2,x3,x,y1,y2,y3,y,yjump
    REAL(kr)    :: relerru,relerru1,relerru2
    REAL(kr)    :: relerrh,relerrh1,relerrh2
    REAL(kr)    :: deg_2_rad = 3.1415926535898d00 / 180.0d00
    
    
!   test_case = 1:   rectilinear channel flow
!             = 2:   flow under bridge with cylinders
!             = 3-4: dam break
!             = 5:   cylindrical dam break
!             = 6:   oblique hydraulic jump
!             = 7:   channel with bump
    test_case = 6
    relerru   = zero
    relerru1  = zero
    relerru2  = zero
    relerrh   = zero
    relerrh1  = zero
    relerrh2  = zero
    ok        = 1
    
    select case (test_case)
       case (6) ! oblique hydraulic jump
          h1 = 1.0d00
          u1 = 9.0d00
          g  = 9.81d00
          theta = 10.0d00        * deg_2_rad
          beta  = 29.88171415d00 * deg_2_rad
          h2 = h1 * tan(beta) / tan(beta-theta)
          u2 = u1 * cos(beta) / cos(beta-theta)
          
          DO i=1,nbrElem
            
            x1 = node(elem(i,1),1)
            x2 = node(elem(i,2),1)
            x3 = node(elem(i,3),1)
            y1 = node(elem(i,1),2)
            y2 = node(elem(i,2),2)
            y3 = node(elem(i,3),2)
            x  = (x1+x2+x3)/3.0d00
            y  = (y1+y2+y3)/3.0d00
            yjump = (x-5.0d00)*TAN(beta)
            
            IF (x.GE.5.0d00 .and. y.LE.yjump ) THEN ! downstream of jump
              href = h2
              uref = u2
            ELSE ! upstream of jump
              href = h1
              uref = u1
            ENDIF
            
            relerrh1 = relerrh1 + ( U0(i*nbvar-2) - href )**2
            relerrh2 = relerrh2 + href**2
              
            U = sqrt( (U0(i*nbvar-1)/U0(i*nbvar-2))**2 + (U0(i*nbvar)/U0(i*nbvar-2))**2 )
            relerru1 = relerru1 + ( U - uref )**2
            relerru2 = relerru2 + uref**2
            
          ENDDO
          
          relerrh = sqrt( relerrh1 / relerrh2 )
          relerru = sqrt( relerru1 / relerru2 )
          
          write(*,*) relerrh,relerru
          
       case default
            
    end select
    
300 CONTINUE
    
end subroutine relative_error

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