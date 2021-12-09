!##########################################################
! SUBROUTINE read_parameters
!  Goal: read the parameters from the parameter file given 
!        by the user in the command window
!##########################################################
SUBROUTINE read_parameters(file_name,ok)
    USE module_shallow
    IMPLICIT NONE
    
    ! Subroutine parameters
    CHARACTER(LEN=100),INTENT(IN)  :: file_name
    INTEGER(ki),       INTENT(OUT) :: ok
    
    ! Local parameters
    CHARACTER(LEN=200)::tmp 
    CHARACTER(LEN=256) :: my_iomsg
    INTEGER(ki) :: i, ierr

    ok = 1
    OPEN(unit=20,file=file_name,status="old",iostat=ierr,iomsg=my_iomsg,form='formatted')
    IF (ierr .NE. 0) THEN
        WRITE(*,*) my_iomsg
        ok = 0
        RETURN
    ENDIF
    
    READ(20,iostat=ierr,err=8,fmt='(a40,a100)') mesh_file,tmp
    READ(20,iostat=ierr,err=8,fmt='(I9,a100)') nbrBC,tmp
    IF (nbrBC.le.0 .or. nbrBC.gt.5) THEN
      WRITE(*,*) "read_parameters: the number of boundary types is not in the range [1-5]"
      GOTO 50
    ENDIF
    CLTable = 0
    DO i=1,nbrBC
      READ(20,iostat=ierr,err=8,fmt='(I9,a100)') CLTable(1,i),tmp
      CLTable(2,i) = i-1
    ENDDO
    READ(20,iostat=ierr,err=8,fmt='(a40,a100)') file_gmsh,tmp
    READ(20,iostat=ierr,err=8,fmt='(a40,a100)') file_dat,tmp
    READ(20,iostat=ierr,err=8,fmt='(a40,a100)') file_restart,tmp
    READ(20,iostat=ierr,err=8,fmt='(I9,a100)')  nTime,tmp
    READ(20,iostat=ierr,err=8,fmt='(I9,a100)')  shownTime,tmp
    READ(20,iostat=ierr,err=8,fmt='(I9,a100)')  savenTime,tmp
    READ(20,iostat=ierr,err=8,fmt='(ES15.7E3,a100)')  ggrav,tmp
    READ(20,iostat=ierr,err=8,fmt='(ES15.7E3,a100)')  manning,tmp
    READ(20,iostat=ierr,err=8,fmt='(I9,a100)')  steady,tmp
    READ(20,iostat=ierr,err=8,fmt='(ES15.7E3,a100)')  CFL,tmp
    READ(20,iostat=ierr,err=8,fmt='(ES15.7E3,a100)')  deltaTfixed,tmp
    READ(20,iostat=ierr,err=8,fmt='(I9,a100)')  restart,tmp
    
    CLOSE (20)
! ! normal exit
    CALL print_parameters(file_name)
    RETURN

! ! exit on error during read operation 
 8   WRITE(*,*) "I/O error in parameters : ", ierr
    
 50 CONTINUE   
    ok = 0
    call print_correct_parameters()

END SUBROUTINE read_parameters

!##########################################################
! SUBROUTINE print_correct_parameters
!  Goal: print on the screen the parameters that should 
!        appear in the parameter file
!##########################################################
SUBROUTINE print_correct_parameters()
    USE module_shallow, ONLY : ki
    IMPLICIT NONE
    
    ! Local parameters
    INTEGER(ki) :: i
    
    WRITE(*,*) "The file 'parameters' is located in the current directory and contains the "
    WRITE(*,*) "following lines"
    WRITE(*,*) ('*',i=1,79)
    WRITE(*,'(a12,18x,a)') "CHARACTER*40",": name of mesh file"
    WRITE(*,'(a7,23x,a)')  "INTEGER",": number of boundary types as built in gmsh, see"
    WRITE(*,'(32x,a)')     "the following four physical tags"
    WRITE(*,'(a7,23x,a)')  "INTEGER",": physical tag associated to the Inlet"
    WRITE(*,'(a7,23x,a)')  "INTEGER",": [optional] physical tag associated to the Outlet"
    WRITE(*,'(a7,23x,a)')  "INTEGER",": [optional] physical tag associated to the Wall"
    WRITE(*,'(a7,23x,a)')  "INTEGER",": [optional] physical tag associated to the Symmetry"
    WRITE(*,'(a12,18x,a)') "CHARACTER*40",": name of the output gmsh file"
    WRITE(*,'(a12,18x,a)') "CHARACTER*40",": name of the output dat file"
    WRITE(*,'(a12,18x,a)') "CHARACTER*40",": name of the restart dat file"
    WRITE(*,'(a7,23x,a)')  "INTEGER",": Total number of temporal steps"
    WRITE(*,'(a7,23x,a)')  "INTEGER",": Statistics are printed every X steps"
    WRITE(*,'(a7,23x,a)')  "INTEGER",": Solution is saved every X steps"
    WRITE(*,'(a4,26x,a)')  "REAL",": g : the gravity constant"
    WRITE(*,'(a4,26x,a)')  "REAL",": n : Manning roughness coefficient"
    WRITE(*,'(a7,23x,a)')  "INTEGER",": [1] steady or [0] unsteady flow"
    WRITE(*,'(a4,26x,a)')  "REAL",": CFL number (used only for steady flows)"
    WRITE(*,'(a4,26x,a)')  "REAL",": time step [s] (used only for unsteady flows)"
    WRITE(*,'(a7,23x,a)')  "INTEGER",": restart from previous simulation [0] no or [1] yes"

    write(*,*) ('*',i=1,79)
    
END SUBROUTINE print_correct_parameters

!##########################################################
! SUBROUTINE print_parameters
!  Goal: print on the screen the parameters that were read
!        in the parameter file
!##########################################################
SUBROUTINE print_parameters(file_name)
    USE MODULE_SHALLOW
    IMPLICIT NONE
    
    ! Subroutine parameters
    CHARACTER(LEN=100),INTENT(IN)  :: file_name
    
    ! Local parameters
    INTEGER(ki) :: i
    
    WRITE(*,*) "The following parameters were read from parameter file '",trim(file_name),"'"
    WRITE(*,'(a40,a)')          mesh_file,": name of mesh file"
    WRITE(*,'(i9,31x,a)')       nbrBC,": number of boundary types as built in gmsh, see"
    DO i=1,nbrBC
      WRITE(*,'(i9,31x,a)')     CLTable(1,i),": boundary physical tag"
    ENDDO
    WRITE(*,'(a40,a)')          file_gmsh,": name of the output gmsh file"
    WRITE(*,'(a40,a)')          file_dat,": name of the output dat file"
    WRITE(*,'(a40,a)')          file_restart,": name of the restart dat file"
    WRITE(*,'(i9,31x,a)')       nTime,": Total number of temporal steps"
    WRITE(*,'(i9,31x,a)')       shownTime,": Statistics are printed every X steps"
    WRITE(*,'(i9,31x,a)')       savenTime,": Solution is saved every X steps"
    WRITE(*,'(ES15.7E3,25x,a)') ggrav,": g : the gravity constant"
    WRITE(*,'(ES15.7E3,25x,a)') manning,": n : Manning roughness coefficient"
    WRITE(*,'(i9,31x,a)')       steady,": [1] steady or [0] unsteady flow"
    WRITE(*,'(ES15.7E3,25x,a)') CFL,": CFL number (used only for steady flows)"
    WRITE(*,'(ES15.7E3,25x,a)') deltaTfixed,": time step [s] (used only for unsteady flows)"
    WRITE(*,'(i9,31x,a)')       restart,": restart from previous simulation [0] no or [1] yes"
    
END SUBROUTINE print_parameters
