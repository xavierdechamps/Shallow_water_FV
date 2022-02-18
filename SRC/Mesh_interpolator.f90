PROGRAM MESH_INTERPOLATOR
    USE module_shallow
    USE module_mem_allocate
    IMPLICIT NONE
    
    INTEGER(ki) :: ok,i,j,modu,found
    CHARACTER(LEN=length_names)::mesh_old,mesh_new,sol_input,sol_output
    
    REAL(kr), ALLOCATABLE    :: U0_new(:)
    REAL(kr), ALLOCATABLE    :: node_new(:,:)
    REAL(kr), ALLOCATABLE    :: depth_new(:)
    REAL(kr), ALLOCATABLE    :: BoundCond_new(:,:)
    INTEGER(ki), ALLOCATABLE    :: elem_new(:,:)
    INTEGER(ki), ALLOCATABLE    :: front_new(:,:)
    INTEGER(ki), ALLOCATABLE    :: nbr_nodes_per_elem_new(:)
    REAL(kr)    :: x1,x2,x3,y1,y2,y3,xc,yc
    REAL(kr)    :: cross1,cross2,cross3
    INTEGER(ki) :: nbrNodes_new,nbrElem_new,nbrTris_new,nbrQuads_new,nbrFront_new,nbrInt_new
    
    ! Get the number of arguments
    IF(COMMAND_ARGUMENT_COUNT().NE.4)THEN
      WRITE(*,*) "Incorrect number of arguments. Please launch the program as"
      WRITE(*,*) "mesh_interpolator  mesh_old.msh  mesh_new.msh  sol_input.dat  sol_output.dat"
      GOTO 200
    ENDIF 
    CALL get_command_argument(1,mesh_old)
    CALL get_command_argument(2,mesh_new)
    CALL get_command_argument(3,sol_input)
    CALL get_command_argument(4,sol_output)
    
    write(*,*) mesh_old
    write(*,*) mesh_new
    write(*,*) sol_input
    write(*,*) sol_output
        
    ! Browse the old/new mesh to get the size of the arrays
    CALL browse_gmsh(mesh_old,length_names,nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the browsing of the old mesh"
      GOTO 200
    endif
    CALL browse_gmsh(mesh_new,length_names,nbrNodes_new,nbrElem_new,nbrTris_new,nbrQuads_new,nbrFront_new,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the browsing of the new mesh"
      GOTO 200
    endif
    
    ! Allocate the memory for the arrays
    CALL mem_allocate(node,front,elem,nbr_nodes_per_elem,U0,depth,BoundCond,dt,Source,&
&                     edges,fnormal,geom_data,cell_data_n,edges_ind,fnormal_ind,&
&                     nbvar*nbrElem,nbrNodes,nbrElem,nbrFront,nbrInt,1)
    CALL mem_allocate(node_new,front_new,elem_new,nbr_nodes_per_elem_new,U0_new,depth_new,BoundCond_new,dt,Source,&
&                     edges,fnormal,geom_data,cell_data_n,edges_ind,fnormal_ind, &
&                     nbvar*nbrElem_new,nbrNodes_new,nbrElem_new,nbrFront_new,nbrInt_new,1)

    ! Read the mesh and the initial solution / boundary conditions
    CALL read_gmsh(U0,nbvar*nbrElem,mesh_old,length_names,node,elem,nbr_nodes_per_elem,front,&
&                 depth,BoundCond,nbrNodes,nbrElem,nbrTris,nbrQuads,nbrFront,1,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the reading of the mesh"
      GOTO 200
    ENDIF
    CALL read_gmsh(U0_new,nbvar*nbrElem_new,mesh_new,length_names,node_new,elem_new,nbr_nodes_per_elem_new,front_new,&
&                  depth_new,BoundCond_new,nbrNodes_new,nbrElem_new,nbrTris_new,nbrQuads_new,nbrFront_new,1,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the reading of the mesh"
      GOTO 200
    ENDIF
    
    ! Read the solution in the dat file
    CALL read_solution(U0,nbvar*nbrElem,sol_input,length_names,ok)
    IF (ok == 0) THEN
      WRITE(*,*) "The program hasn't started because of a problem during the reading of the mesh"
      GOTO 200
    ENDIF
    
    WRITE(*,*) "Interpolating the mesh data on the new mesh..."
    modu = nbrElem_new / 10
    DO i=1,nbrElem_new
      ! Get area of element
      x1 = node_new(elem_new(i,1),1)
      y1 = node_new(elem_new(i,1),2)
      x2 = node_new(elem_new(i,2),1)
      y2 = node_new(elem_new(i,2),2)
      x3 = node_new(elem_new(i,3),1)
      y3 = node_new(elem_new(i,3),2)
      xc = (x1+x2+x3)/3.0d00
      yc = (y1+y2+y3)/3.0d00
      found = 0
      
      DO j=1,nbrElem
        
        ! Get area of sub element 1
        x1 = node(elem(j,1),1)
        y1 = node(elem(j,1),2)
        x2 = node(elem(j,2),1)
        y2 = node(elem(j,2),2)
        x3 = node(elem(j,3),1)
        y3 = node(elem(j,3),2)
        
        ! AP x AB
        cross1 = ( xc-x1 )*( y2-y1 ) - ( x2-x1 )*( yc-y1 )
        ! BP x BC
        cross2 = ( xc-x2 )*( y3-y2 ) - ( x3-x2 )*( yc-y2 )
        ! CP x CA
        cross3 = ( xc-x3 )*( y1-y3 ) - ( x1-x3 )*( yc-y3 )
        
        IF ( ( cross1.GE.0.0d00 .AND. cross2.GE.0.0d00 .AND. cross3.GE.0.0d00 ) .OR. &
&            ( cross1.LE.0.0d00 .AND. cross2.LE.0.0d00 .AND. cross3.LE.0.0d00 )  ) THEN
          
          U0_new( i*nbvar-2 ) = U0( j*nbvar-2 )
          U0_new( i*nbvar-1 ) = U0( j*nbvar-1 )
          U0_new( i*nbvar   ) = U0( j*nbvar   )
          found = 1
          EXIT
        ENDIF
        
      ENDDO
      
      IF (j.EQ.nbrElem .AND. found.EQ.0) THEN
        WRITE(*,*) "Could not find an adequate match for element ",i
      ENDIF
      IF ( MODULO (i-1,modu).EQ.0 ) WRITE(*,'(i4,a)') i*100/nbrElem_new,"% done"
      
    ENDDO
    
    ! Save the solution in the .mesh and .dat formats
    ! CALL write_gmsh(U0_new,nbvar*nbrElem_new,"test_new.msh",length_names,node_new,elem_new,front_new,nbrNodes_new,nbrElem_new,nbrTris_new,nbrQuads_new,nbrFront_new,0,0)
    CALL write_solution(U0_new,nbvar*nbrElem_new,sol_output,length_names,ok)
    
    ! Deallocate the memory for the arrays
    CALL mem_deallocate(node,front,elem,nbr_nodes_per_elem,U0,depth,BoundCond,dt,Source,&
&                       edges,fnormal,geom_data,cell_data_n,edges_ind,fnormal_ind)
    CALL mem_deallocate(node_new,front_new,elem_new,nbr_nodes_per_elem_new,U0_new,depth_new,BoundCond_new,dt,Source,&
&                       edges,fnormal,geom_data,cell_data_n,edges_ind,fnormal_ind)
    
200 CONTINUE
    WRITE(*,*) "End of the program"
    
END PROGRAM

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