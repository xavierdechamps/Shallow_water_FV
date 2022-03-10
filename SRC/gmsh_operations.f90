
!##########################################################
! SUBROUTINE write_solution
!  Goal: write the restart solution in the dat file
!##########################################################
SUBROUTINE write_solution(U0,lengU0,file_dat,lengch,ok)
    USE MODULE_SHALLOW, only : kr,ki
    IMPLICIT NONE

    ! Subroutine parameters
    INTEGER(ki), INTENT(OUT) :: ok
    INTEGER(ki), INTENT(IN)  :: lengU0,lengch
    REAL(kr), INTENT(IN)     :: U0(lengU0)
    CHARACTER(LEN=lengch), INTENT(IN) :: file_dat
    
    ! Local parameters
    INTEGER(ki) :: i, ierr

    ok = 1
    OPEN(UNIT=20,FILE=file_dat,STATUS="replace",IOSTAT=ierr,FORM='formatted')
    IF (ierr .NE. 0) THEN
      WRITE(*,*) "Error writing the solution in .dat format"
      ok = 0
    END IF

    WRITE(20,'(ES15.7E3)') (U0(i),i=1,lengU0)

    CLOSE(UNIT=20)

END SUBROUTINE write_solution

!##########################################################
! SUBROUTINE read_solution
!  Goal: read the restart solution in the dat file
!##########################################################
SUBROUTINE read_solution(U0,lengU0,file_restart,lengch,ok)
    USE MODULE_SHALLOW, only : kr,ki,zero
    IMPLICIT NONE

    ! Subroutine parameters
    INTEGER(ki), INTENT(OUT):: ok
    INTEGER(ki), INTENT(IN) :: lengU0,lengch
    REAL(kr), INTENT(OUT)   :: U0(lengU0)
    CHARACTER(LEN=lengch), INTENT(IN) :: file_restart
    
    ! Local parameters
    INTEGER(ki) :: i, ierr, count
    REAL(kr) :: a

    ok = 1
    U0 = zero
    OPEN(UNIT=20,FILE=file_restart,STATUS="old",IOSTAT=ierr,FORM='formatted')
    IF (ierr .NE. 0) THEN
      WRITE(*,*) "Error occuring : file '",trim(file_restart),"' not found "
      ok = 0
      GOTO 200
    END IF
    count = 0

    READ(20,IOSTAT=ierr,FMT='(ES15.7E3)') a
    DO 
       count = count + 1
       READ(UNIT=20,IOSTAT=ierr,FMT='(ES15.7E3)') a
       IF (ierr .LT. 0) EXIT
    END DO

    IF (count .EQ. lengU0) THEN
       REWIND(20)
       READ(20,*) (U0(i), i=1,count)
    ELSE
       ok = 0
       WRITE(*,*) "Error occuring during the reading of the former solution : wrong vector length"
       WRITE(*,*) "Length of the data in the restart file: ",count
       WRITE(*,*) "Expected length of the data in the restart file: ",lengU0
    END IF
    CLOSE(UNIT=20)
    
200 CONTINUE  
END SUBROUTINE read_solution

!##########################################################
!SUBROUTINE write_gmsh
!  Goal: write the solution in the Gmsh .msh format
!##########################################################
SUBROUTINE write_gmsh(U0,lengU0,file_gmsh,lengch,node,elem,front,nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,k,count)
    USE MODULE_SHALLOW, only : kr,ki,ggrav,nbvar,eps
    IMPLICIT NONE

    ! Subroutine parameters
    INTEGER(ki), INTENT(IN) :: k, count, nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront
    INTEGER(ki), INTENT(IN) :: lengU0,lengch
    REAL(kr), INTENT(IN)    :: U0(lengU0)
    REAL(kr), INTENT(IN)    :: node(nbrNodes,2)
    INTEGER(ki), INTENT(IN)    :: elem(nbrElem,5)
    INTEGER(ki), INTENT(IN)    :: front(nbrFront,4)
    CHARACTER(LEN=lengch), INTENT(IN) :: file_gmsh
    
    ! Local parameters
    INTEGER(ki) :: ierr, i, fin
    REAL(kr) :: h, u, v
    REAL(kr), DIMENSION(1) :: Froude(1:nbrElem)
    
    fin = 0
    
    ! Write the mesh data (nodes and elements)
    IF (k==0) THEN
       OPEN(UNIT=10,FILE=file_gmsh,STATUS="replace",IOSTAT=ierr,FORM='formatted')
       WRITE(10,'(T1,A11)') "$MeshFormat"
       WRITE(10,'(T1,A7)') "2.2 0 8"
       WRITE(10,'(T1,A14)') "$EndMeshFormat"
       WRITE(10,'(T1,A6)') "$Nodes"
       WRITE(10,'(T1,I9)') nbrNodes
       WRITE(10,'(T1,I9,2ES25.16E3,F4.1)') (i, node(i,:),0.,i=1,nbrNodes)
       WRITE(10,'(T1,A9)') "$EndNodes"
       WRITE(10,'(T1,A9)') "$Elements"
       WRITE(10,'(T1,I9)') nbrElem+nbrFront
       WRITE(10,'(T1,I9,2I2,I9,I2,2I9)') (i,1,2,front(i,3),1,front(i,1:2),i=1,nbrFront)       
       IF (nbrTris.NE.0) write(10,'(T1,I9,2I2,I9,I2,3I9)') (i+nbrFront,2,2,elem(i,5),1,elem(i,1:3),i=1,nbrTris)
       IF (nbrQuads.NE.0) write(10,'(T1,I9,2I2,I9,I2,4I9)') (i+nbrFront+nbrTris,3,2,elem(i+nbrTris,5),1, &
&                                                            elem(i+nbrTris,1:4),i=1,nbrQuads)
       WRITE(10,'(T1,A12)') "$EndElements"
    ELSE
       OPEN(UNIT=10,FILE=file_gmsh,STATUS="old",ACCESS="append",IOSTAT=ierr,FORM='formatted')
    END IF
    
    !*************************************
    ! Write the height
    WRITE(10,'(T1,A12)') "$ElementData"
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,A9)') '"Height"'
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,I9)') k
    WRITE(10,'(T1,A1)') "3"
    WRITE(10,'(T1,I9)') count
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,I9)') nbrElem
    WRITE(10,'(T1,I9,ES25.16E3)') (i+nbrFront, U0(i*nbvar-2),i=1,nbrElem)
    WRITE(10,'(T1,A15)') "$EndElementData"
    
    !*************************************
    ! Write the velocity vector, compute the Froude number in a loop
    WRITE(10,'(T1,A12)') "$ElementData"
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,A10)') '"Velocity"'
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,I9)') k
    WRITE(10,'(T1,A1)') "3"
    WRITE(10,'(T1,I9)') count
    WRITE(10,'(T1,A1)') "3"
    WRITE(10,'(T1,I9)') nbrElem
    DO i=1,nbrElem
       h = U0(i*nbvar-2)
       IF (h<eps) h = eps
       u = U0(i*nbvar-1)/h
       v = U0(i*nbvar)/h
       Froude(i) = sqrt(u*u+v*v)/sqrt(ggrav*h)
       WRITE(10,'(T1,I9,2ES25.16E3,F4.1)') i+nbrFront, u, v, 0.
    END DO
    WRITE(10,'(T1,A15)') "$EndElementData"
    
    !*************************************
    ! Write the Froude number, calculated in the previous block
    WRITE(10,'(T1,A12)') "$ElementData"
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,A8)') '"Froude"'
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,I9)') k
    WRITE(10,'(T1,A1)') "3"
    WRITE(10,'(T1,I9)') count
    WRITE(10,'(T1,A1)') "1"
    WRITE(10,'(T1,I9)') nbrElem
    WRITE(10,'(T1,I9,ES25.16E3)') (i+nbrFront, Froude(i),i=1,nbrElem)
    WRITE(10,'(T1,A15)') "$EndElementData"
    
    !*************************************
    CLOSE(UNIT=10)

END SUBROUTINE write_gmsh

!##########################################################
! SUBROUTINE read_gmsh
!  Goal: read the mesh and the initial solution in the Gmsh .msh format
!##########################################################
SUBROUTINE read_gmsh(U0,lengU0,mesh_file,lengch,node,elem,nbr_nodes_per_elem,front,&
&                    depth,BoundCond,nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,skip_data,ok)
    USE MODULE_SHALLOW, only : kr,ki,nbvar,time_begin, time_end
    IMPLICIT NONE

    ! Subroutine parameters
    INTEGER(ki), INTENT(OUT) :: ok
    INTEGER(ki), INTENT(IN) :: nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront
    INTEGER(ki), INTENT(IN) :: lengU0,lengch,skip_data
    REAL(kr), INTENT(OUT)    :: U0(lengU0)
    REAL(kr), INTENT(OUT)    :: node(nbrNodes,2)
    REAL(kr), INTENT(OUT)    :: depth(nbrElem)
    REAL(kr), INTENT(OUT)    :: BoundCond(nbrFront,3)
    INTEGER(ki), INTENT(OUT)    :: elem(nbrElem,5)
    INTEGER(ki), INTENT(OUT)    :: front(nbrFront,4)
    INTEGER(ki), INTENT(OUT)    :: nbr_nodes_per_elem(nbrElem)
    CHARACTER(LEN=lengch), INTENT(IN) :: mesh_file

    ! Local parameters
    CHARACTER(len=256)    :: line
    INTEGER(ki) :: i, ierr, a, b, c, nbrElemTot=0,nbrNodeEdge=0
    INTEGER(ki) :: istep
    REAL(kr) :: header,tmp
    LOGICAL :: header21 = .true.
    EXTERNAL :: find_elem_front, time_display

    CALL sampletime(time_begin)
    ok = 1
    WRITE(*,*) "-------------------------------------------------------"
    WRITE(*,*) "Name of the mesh data : ", TRIM(ADJUSTL(mesh_file))
    
    istep = 1
    
    OPEN(UNIT=10,FILE=mesh_file,STATUS="old",IOSTAT=ierr,FORM='formatted')
    IF (ierr .NE. 0) THEN
      WRITE(*,*) "Error occuring : file '",TRIM(mesh_file),"' not found "
      ok = 0
      GOTO 100
    END IF

    READ(10,*) line, header, b, c
    
    IF (header .LT. 2.) THEN
      WRITE(*,*) "Error : wrong format for gmsh data"
      ok = 0
      GOTO 100
    END IF

    IF (header .GT. 2.1) header21 = .false.
    
    istep = 2
    
    DO WHILE (line .NE. "$Nodes")
      READ(10,*) line
    END DO

    READ(10,*) a

    DO i=1,nbrNodes ! Read the node coordinates
      READ(10,*) a, node(i,1), node(i,2), tmp
    END DO
    WRITE(*,101) nbrNodes
    
    istep = 3

    DO WHILE (line .NE. "$Elements")
      READ(10,*) line
    END DO

    READ(10,*) a
    
    ! front (:,1:2) : node IDs composing the boundary edge
    ! front (:,3)   : physical tag of the boundary edge
    ! front (:,4)   : ID of the 2D element linked to the boundary edge
    
    ! elem(:,1:4)   : node IDs composing the 2D element
    ! elem(:,5)     : physical tag of the 2D element
    
    IF (header21) THEN  ! Read the boundary edges
       READ(10,*) (a,a,a,front(i,3),a,a,front(i,1),front(i,2),i=1,nbrFront)
    ELSE
       READ(10,*) (a,a,a,front(i,3),a,front(i,1),front(i,2),i=1,nbrFront)
    END IF
    WRITE(*,102) nbrFront
    
    IF (header21) THEN  ! Read the 2D elements
       IF (nbrTris.GT.0) THEN
         READ(10,*) (a,a,a,elem(i,5),a,a,elem(i,1),elem(i,2),elem(i,3),i=1,nbrTris) ! Triangles
       ENDIF
       IF (nbrQuads.GT.0) THEN
         READ(10,*) (a,a,a,elem(i,5),a,a,elem(i,1),elem(i,2),elem(i,3),elem(i,4),i=nbrTris+1,nbrTris+nbrQuads) ! Quadrangles
       ENDIF
    ELSE
       IF (nbrTris.GT.0) THEN
         READ(10,*) (a,a,a,elem(i,5),a,elem(i,1),elem(i,2),elem(i,3),i=1,nbrTris) ! Triangles
       ENDIF
       IF (nbrQuads.GT.0) THEN
         READ(10,*) (a,a,a,elem(i,5),a,elem(i,1),elem(i,2),elem(i,3),elem(i,4),i=nbrTris+1,nbrTris+nbrQuads) ! Quadrangles
       ENDIF
    END IF
    DO i=1,nbrTris
      nbr_nodes_per_elem(i) = 3
    ENDDO
    DO i=nbrTris+1,nbrTris+nbrQuads
      nbr_nodes_per_elem(i) = 4
    ENDDO
    WRITE(*,103) nbrTris
    WRITE(*,104) nbrQuads

    DO i=1,nbrFront
       CALL find_elem_front(i,front(i,4))
    END DO
    
    IF (skip_data.EQ.1) GOTO 100
    
    ! Initial height
    istep = 4
    
    DO WHILE (line .NE. "$ElementData")
      READ(10,*,END=100) line
    END DO
    DO i=1,7
      READ(10,*,END=100) line
    ENDDO
    READ(10,*,END=100) a
    IF (a.ne.nbrElem) THEN
      WRITE(*,*) "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
      WRITE(*,*) "Error reading the initial height."
      WRITE(*,*) "The number of elements with initial height is",a
      WRITE(*,*) "The number of 2D elements is",nbrElem
      WRITE(*,*) "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
      ok = 0
      GOTO 100
    ENDIF
    READ(10,*,END=100) (a,U0((i-1)*nbvar+1),i=1,nbrElem)
    
    ! Initial velocity
    
    istep = 5
    
    DO WHILE (line .NE. "$ElementData")
      READ(10,*,END=100) line
    END DO
    DO i=1,7
      READ(10,*,END=100) line
    ENDDO
    READ(10,*,END=100) a
    
    IF (a.ne.nbrElem) THEN
      WRITE(*,*) "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
      WRITE(*,*) "Error reading the initial velocity."
      WRITE(*,*) "The number of elements with initial velocity is",a
      WRITE(*,*) "The number of 2D elements is",nbrElem
      WRITE(*,*) "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
      ok = 0
      GOTO 100
    ENDIF
    READ(10,*,END=100) (a,U0((i-1)*nbvar+2),U0((i-1)*nbvar+3),tmp,i=1,nbrElem)
    DO i=1,nbrElem
      U0((i-1)*nbvar+2) = U0((i-1)*nbvar+1)*U0((i-1)*nbvar+2)
      U0((i-1)*nbvar+3) = U0((i-1)*nbvar+1)*U0((i-1)*nbvar+3)
    ENDDO
    
    ! depth
    
    istep = 6
    
    DO WHILE (line .NE. "$ElementData")
      READ(10,*,END=100) line
    END DO
    DO i=1,7
      READ(10,*,END=100) line
    ENDDO
    READ(10,*,END=100) a
    IF (a.ne.nbrElem) THEN
      WRITE(*,*) "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
      WRITE(*,*) "Error reading the depth."
      WRITE(*,*) "The number of elements with depth is",a
      WRITE(*,*) "The number of 2D elements is",nbrElem
      WRITE(*,*) "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
      ok = 0
      GOTO 100
    ENDIF
    READ(10,*,END=100) (a,depth(i),i=1,nbrElem)
    
    ! Boundary condition - height
    
    istep = 7
    
    DO WHILE (line .NE. "$ElementData")
      READ(10,*,END=100) line
    END DO
    DO i=1,7
      READ(10,*,END=100) line
    ENDDO
    READ(10,*,END=100) a
    IF (a.ne.nbrFront) THEN
      WRITE(*,*) "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
      WRITE(*,*) "Error reading the boundary condition on the height."
      WRITE(*,*) "The number of 1D elements with imposed height is",a
      WRITE(*,*) "The number of 1D boundary elements is",nbrFront
      WRITE(*,*) "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
      ok = 0
      GOTO 100
    ENDIF
    READ(10,*,END=100) (a,BoundCond(i,1),i=1,nbrFront)
        
    ! Boundary condition - velocity
    
    istep = 8
    
    DO WHILE (line .NE. "$ElementData")
      READ(10,*,END=100) line
    END DO
    DO i=1,7
      READ(10,*,END=100) line
    ENDDO
    READ(10,*,END=100) a    
    IF (a.ne.nbrFront) THEN
      WRITE(*,*) "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
      WRITE(*,*) "Error reading the boundary condition on the velocity."
      WRITE(*,*) "The number of 1D elements with imposed velocity is",a
      WRITE(*,*) "The number of 1D boundary elements is",nbrFront
      WRITE(*,*) "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
      ok = 0
      GOTO 100
    ENDIF
    READ(10,*,END=100) (a,BoundCond(i,2),BoundCond(i,3),tmp,i=1,nbrFront)
    
100 CALL sampletime(time_end)
    CALL time_display
    WRITE(*,*) "-------------------------------------------------------"
    101 FORMAT(" Number of node : ",T41,I6)
    102 FORMAT(" Number of edge elements :", T41, I6)
    103 FORMAT(" Number of triangular surface elements :",T41,I6)
    104 FORMAT(" Number of rectangular surface elements :",T41,I6)
    CLOSE(UNIT=10)
    
END SUBROUTINE read_gmsh

!##########################################################
! SUBROUTINE read_gmsh
!  Goal: read the mesh and the initial solution in the Gmsh .msh format
!##########################################################
SUBROUTINE browse_gmsh(mesh_file,lengch,nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,ok)
    USE MODULE_SHALLOW, only : ki,kr
    IMPLICIT NONE
    
    ! Subroutine parameters
    INTEGER(ki), INTENT(OUT) :: ok
    INTEGER(ki), INTENT(OUT) :: nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront
    INTEGER(ki), INTENT(IN)  :: lengch
    CHARACTER(LEN=lengch), INTENT(IN) :: mesh_file

    ! Local parameters
    CHARACTER(len=256)    :: line
    INTEGER(ki) :: i, j, ierr, a, b, c, nbrElemTot=0,nbrNodeEdge=0
    INTEGER(ki) :: istep
    REAL(kr) :: header,tmp
    LOGICAL :: header21 = .true.
    EXTERNAL :: find_elem_front, time_display

    ok = 1
    istep = 1
    nbrNodes = 0
    nbrElem = 0
    nbrTris= 0
    nbrQuads = 0
    nbrFront = 0
    
    OPEN(UNIT=10,FILE=mesh_file,STATUS="old",IOSTAT=ierr,FORM='formatted')
    IF (ierr .NE. 0) THEN
      WRITE(*,*) "Error occuring : file '",TRIM(mesh_file),"' not found "
      ok = 0
      GOTO 100
    END IF

    READ(10,*) line, header, b, c
    
    IF (header .LT. 2.) THEN
      WRITE(*,*) "Error : wrong format for gmsh data"
      ok = 0
      GOTO 100
    END IF

    IF (header .GT. 2.1) header21 = .false.
    
    istep = 2
    
    DO WHILE (line .NE. "$Nodes")
      READ(10,*) line
    END DO

    READ(10,*) nbrNodes
    
    DO WHILE (line .NE. "$Elements")
      READ(10,*) line
    END DO

    READ(10,*) nbrElemTot
    
    DO j=1,nbrElemTot
      READ(10,*) a,i
      IF (i==1) THEN
        nbrFront = nbrFront + 1
      ELSE IF (i==2) THEN
        nbrTris = nbrTris + 1
      ELSE IF (i==3) THEN
        nbrQuads = nbrQuads + 1
      ELSE IF (i==15) THEN
        nbrNodeEdge = nbrNodeEdge + 1
      END IF
    END DO
    nbrElem = nbrTris + nbrQuads

100 CLOSE(UNIT=10)
    
END SUBROUTINE browse_gmsh
    