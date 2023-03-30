!##########################################################
! SUBROUTINE get_normal_to_cell
!  Goal: compute the geometrical data required by the finite volume
!        method, after having read the mesh.
!##########################################################
SUBROUTINE get_normal_to_cell()
    USE module_shallow
    USE OMP_LIB
    IMPLICIT NONE

    ! Local parameters
    INTEGER(ki) :: i,j,id,q,p,k,l,front1,front2,idx,goahead
    INTEGER(ki) :: idn,typec,edgeID=0,nID=0
    INTEGER(ki) :: myID,numThreads
    INTEGER(ki), DIMENSION(4) :: pos, pos2
    REAL(kr), DIMENSION(4,2)  :: xy
    REAL(kr) :: xa,xb,xc,xd,ya,yb,yc,yd
    LOGICAL :: test = .false.
    REAL(kr), ALLOCATABLE :: cell_data(:,:)
    REAL(kr), ALLOCATABLE :: internal(:,:)
    INTEGER(ki), ALLOCATABLE :: cell_node_ids(:,:,:)
    INTEGER(ki), ALLOCATABLE :: cell_data_ind(:,:)
    INTEGER(ki), ALLOCATABLE :: internal_ind(:,:)
    
    INTEGER(ki), ALLOCATABLE :: fnormal_ind_omp(:,:,:)
    REAL(kr),    ALLOCATABLE :: fnormal_omp(:,:,:)
    INTEGER(ki), ALLOCATABLE :: nId_omp(:)
    ! INTEGER(ki), ALLOCATABLE :: edgeID_omp(:)
    ! INTEGER(ki), ALLOCATABLE :: internal_ind_omp(:,:,:)
    
    CALL sampletime(time_begin)

    ALLOCATE(cell_data(1:nbrElem,1:10))
    ALLOCATE(cell_data_ind(1:nbrElem,1:4)) ! IDs of the 4 2D elements surrounding each 2D element
    ALLOCATE(internal(1:4*nbrElem,1:6))
    ALLOCATE(internal_ind(1:4*nbrElem,1:2))
    ALLOCATE(cell_node_ids(1:nbrElem,1:4,1:2))
    
    test = .FALSE.
    numThreads = 1
!$OMP PARALLEL
    IF (omp_get_num_threads().GT.1) THEN
      numThreads = omp_get_num_threads()
    ENDIF
!$OMP END PARALLEL

!    IF (numThreads.GT.1) THEN
    ALLOCATE (fnormal_ind_omp(0:numThreads-1,1:nbrFront,1:4))
    ALLOCATE (fnormal_omp(0:numThreads-1,1:nbrFront,1:5))
    ALLOCATE (nId_omp(0:numThreads-1))
    ! ALLOCATE (edgeID_omp(0:numThreads-1))
    ! ALLOCATE (internal_ind_omp(0:numThreads-1,1:4*nbrElem/numThreads,1:2))
    fnormal_ind_omp = 0
    fnormal_omp     = zero
    nId_omp         = 0
    ! edgeID_omp = 0
    ! internal_ind_omp = 0
!    ENDIF
    !
    !     a_______________c
    !      \             / \
    !       \           /   \
    !        \    L    /     \
    !         \       /       \
    !          \     /    M    \
    !           \   /           \
    !            \ /_____________\
    !             b               d
    !
    ! geom_data(:,1)   : area of element
    ! geom_data(:,2)   : perimeter of element
    ! geom_data(:,3:4) : XY coordinates of center of element
    ! geom_data_ind(:,1:4) : IDs of the 4 neighbour cells (col 4 is -1 if cell is a triangle)
    
    ! fnormal(:,1:2)   : XY components of normal to boundary edge
    ! fnormal(:,3:4)   : XY components of center of boundary edge
    ! fnormal(:,5)     : Area of subtriangle made by the border edge and the center of cell
    
    ! fnormal_ind(:,1) : 2D element ID to which the boundary edge is attached
    ! fnormal_ind(:,2) : physical tag of the edge as defined in the msh
    ! fnormal_ind(:,3) : CLTable as defined in flux.f90 (Inlet, Outlet, etc.)
    ! fnormal_ind(:,4) : ID of the edge as defined in the msh
    
    ! internal(:,1:2)  : XY components of normal to edge
    ! internal(:,3:4)  : XY components of center of edge
    ! internal(:,5)    : area of left  subtriangle Lbc if edge bc (requested by upwinding of source term)
    ! internal(:,6)    : area of right subtriangle MBC if edge bc (requested by upwinding of source term)
    
    ! internal_ind(:,1) : 2D element ID left to edge
    ! internal_ind(:,2) : 2D element ID right to edge
    
    cell_data = zero
    internal  = zero
    internal_ind  = 0
    fnormal_ind   = 0
    cell_node_ids = 0
    cell_data_ind = 0
    geom_data     = zero
    geom_data_ind = 0
    
    ! Compute the geometrical data for the cells (geom_data) and prepare the data
    ! for the second loop (edges)
!$OMP PARALLEL &
!$OMP& default (shared) &
!$OMP& private (xy,typec,j,k) 
!$OMP DO
    DO i=1,nbrElem
      xy = zero
      typec = nbr_nodes_per_elem(i)
      ! Coordinates of the corner nodes of the cell
      DO j=1,typec
        xy(j,1) = node(elem(i,j),1)
        xy(j,2) = node(elem(i,j),2)        
      END DO
      
      ! Center of cell
      cell_data(i,1:2) = SUM(xy,DIM=1) / REAL(typec,kind=kr)
      
      ! Get the external normals
      DO j=1,typec
        k = j+1
        IF (j.EQ.typec) k = 1
        cell_data(i,2*j+1 : 2*j+2) = (/xy(k,2)-xy(j,2) , xy(j,1)-xy(k,1)/)
        ! Switch the direction of the normals if pointing inwards
        IF ((cell_data(i,2*j+1)*(xy(j,1)-cell_data(i,1)) + cell_data(i,2*j+2)*(xy(j,2)-cell_data(i,2)))<zero) then
          cell_data(i,2*j+1 : 2*j+2) = - cell_data(i,2*j+1 : 2*j+2)
        END IF
        
        ! Find the three 2D elements surrounding the current 2D element and fill cell_data_ind
        CALL find_elem_surf(i,cell_data_ind(i,j),elem(i,j),elem(i,k))
        cell_node_ids(i,j,1:2) = (/ elem(i,j),elem(i,k) /)        
      END DO
      
      ! Get area + perimeter + center
      IF (typec .EQ. 3) THEN
      geom_data(i,1) = 0.5d00*ABS((xy(2,1)-xy(1,1))*(xy(3,2)-xy(1,2))-(xy(3,1)-xy(1,1))*(xy(2,2)-xy(1,2)))
      geom_data(i,2) = SQRT((xy(1,1)-xy(2,1))**2 + (xy(1,2)-xy(2,2))**2) + SQRT((xy(2,1)-xy(3,1))**2 + (xy(2,2)-xy(3,2))**2) + &
&                      SQRT((xy(3,1)-xy(1,1))**2 + (xy(3,2)-xy(1,2))**2)
      ELSE IF (typec .EQ. 4) THEN
      geom_data(i,1) = 0.5d00*ABS((xy(2,1)-xy(1,1))*(xy(3,2)-xy(1,2))-(xy(3,1)-xy(1,1))*(xy(2,2)-xy(1,2))) &
&                    + 0.5d00*ABS((xy(3,1)-xy(1,1))*(xy(4,2)-xy(1,2))-(xy(4,1)-xy(1,1))*(xy(3,2)-xy(1,2)))
      geom_data(i,2) = SQRT((xy(1,1)-xy(2,1))**2 + (xy(1,2)-xy(2,2))**2) + SQRT((xy(2,1)-xy(3,1))**2 + (xy(2,2)-xy(3,2))**2) + &
&                      SQRT((xy(3,1)-xy(4,1))**2 + (xy(3,2)-xy(4,2))**2) + SQRT((xy(1,1)-xy(4,1))**2 + (xy(1,2)-xy(4,2))**2)
      END IF
      geom_data(i,3:4) = cell_data(i,1:2)
      geom_data_ind(i,1:typec) = cell_data_ind(i,1:typec)
    END DO
!$OMP END DO
!$OMP END PARALLEL
    ! Save the x,y components of the normals in the matrix cell_data_n (module_shallow)
    cell_data_n(1:nbrElem,1:8) = cell_data(1:nbrElem,3:10)
        
    ! Second loop to compute the data for the internal/boundary edges
!$OMP PARALLEL &
!$OMP& default (shared) &
!$OMP& private (myID,i,j,id,idn,xa,ya,xb,yb,xc,yc,front1,front2,idx,p,pos,q,pos2,k,l,test) &
!$OMP& reduction(+: nId)
    myID = omp_get_thread_num()
!$OMP DO
    DO i=1,nbrElem
      DO j=10,10+nbr_nodes_per_elem(i)-1
      ! id is the ID of the cell opposite to the edge, if 0 then it is a boundary edge
        id = cell_data_ind(i,j-9)
        IF (id .eq. 0) THEN ! boundary edge
             nId = nId + 1
             idn = (j-10)*2 + 3
             
             nId_omp(myID) = nId_omp(myID) + 1
              
             ! x,y, components of external normal to the edge
!             fnormal(nId,1) = cell_data(i,idn)
!             fnormal(nId,2) = cell_data(i,idn+1)
             
             fnormal_omp(myID,nId_omp(myID),1) = cell_data(i,idn)
             fnormal_omp(myID,nId_omp(myID),2) = cell_data(i,idn+1)
             
             ! Coordinates of center of current cell
             xa = cell_data(i,1)
             ya = cell_data(i,2)
             
             ! Coordinates of the 2 nodes of the edge
             xb = node(cell_node_ids(i,j-9,1),1)
             yb = node(cell_node_ids(i,j-9,1),2)
             xc = node(cell_node_ids(i,j-9,2),1)
             yc = node(cell_node_ids(i,j-9,2),2)
             
             ! x,y components of the center of the edge
!             fnormal(nId,3:4) = 0.5d00*(/ xb+xc , yb+yc/)
             
             fnormal_omp(myID,nId_omp(myID),3:4) = 0.5d00*(/ xb+xc , yb+yc/)
             
             ! Area of subtriangle made by the border edge and the center of current cell
!             fnormal(nId,5) = 0.5d00*ABS((xb-xa)*(yc-ya)-(xc-xa)*(yb-ya))
             
             fnormal_omp(myID,nId_omp(myID),5) = 0.5d00*ABS((xb-xa)*(yc-ya)-(xc-xa)*(yb-ya))
             
             CALL find_vector(front(:,4),nbrFront,i,front1,front2)
             idx = 0
             IF (front2 .ne. 0) THEN ! 2D element with 2 boundary edges
                ! Find which front1 / front2 is the current edge
                xa = node(front(front1,1),1)
                ya = node(front(front1,1),2)
                xb = node(front(front1,2),1)
                yb = node(front(front1,2),2)
                IF (ABS((xa-xb)*cell_data(i,idn)+(ya-yb)*cell_data(i,idn+1))<eps) THEN
                   idx = front1
                ELSE
                   idx = front2
                END IF
             ELSE
                idx = front1
             END IF
!             fnormal_ind(nId,1) = i
!             fnormal_ind(nId,2) = front(idx,3)
!             CALL find_vector(CLTable(1,:),5,fnormal_ind(nId,2),front1,front2)
!             fnormal_ind(nId,3) = CLTable(2,front1)
!             fnormal_ind(nId,4) = idx
             
!             IF (numThreads.GT.1) THEN
              fnormal_ind_omp(myID,nId_omp(myID),1) = i
              fnormal_ind_omp(myID,nId_omp(myID),2) = front(idx,3)
              CALL find_vector(CLTable(1,:),5,front(idx,3),front1,front2)
              fnormal_ind_omp(myID,nId_omp(myID),3) = CLTable(2,front1)
              fnormal_ind_omp(myID,nId_omp(myID),4) = idx
!             ENDIF
             
        ELSE ! internal edge
        
           ! CALL find_mat2(internal_ind_omp(myID,:,:),4*nbrElem/numThreads,i,id,p)
           
           ! IF (p.eq.0) THEN
              ! edgeID_omp(myID) = edgeID_omp(myID) + 1
              ! idn = (j-10)*2 + 3
              
              ! Fill in internal_ind
              ! internal_ind_omp(myID,edgeID_omp(myID),1) = i
              ! internal_ind_omp(myID,edgeID_omp(myID),2) = id
              
              ! Coordinates of center of current cell
              ! xa = cell_data(i,1)
              ! ya = cell_data(i,2)
              
              ! Coordinates of the two nodes of the current internal edge
              ! xb = node(cell_node_ids(i,j-9,1),1)
              ! yb = node(cell_node_ids(i,j-9,1),2)
              ! xc = node(cell_node_ids(i,j-9,2),1)
              ! yc = node(cell_node_ids(i,j-9,2),2)
              
              ! xy components of external normal to edge
              ! internal(edgeID_omp(myID),1) = cell_data(i,idn)
              ! internal(edgeID_omp(myID),2) = cell_data(i,idn+1)
              
              ! Coordinates of center of internal edge
              ! internal(edgeID_omp(myID),3:4) = 0.5d00*(/ xb+xc , yb+yc /)
              
              ! Area of subtriangle left to internal edge
              ! internal(edgeID_omp(myID),5) = 0.5d00*ABS((xb-xa)*(yc-ya)-(xc-xa)*(yb-ya))
              
              ! Coordinates of center of cell opposite to internal edge
              ! xa = cell_data(id,1)
              ! ya = cell_data(id,2)
              
              ! Area of subtriangle right to internal edge
              ! internal(edgeID_omp(myID),6) = 0.5d00*ABS((xb-xa)*(yc-ya)-(xc-xa)*(yb-ya))
          ! END IF
           
!$omp critical           
           ! Search in internal_ind edges that have neighbour cells i and id
           CALL find_mat2(internal_ind(:,:),4*nbrElem,i,id,p)
           
           ! IF (.not.test) THEN
           IF (p.eq.0) THEN
              edgeId = edgeId + 1
              idn = (j-10)*2 + 3
              
              ! Fill in internal_ind
              internal_ind(edgeId,1) = i
              internal_ind(edgeId,2) = id
              
              ! Coordinates of center of current cell
              xa = cell_data(i,1)
              ya = cell_data(i,2)
              
              ! Coordinates of the two nodes of the current internal edge
              xb = node(cell_node_ids(i,j-9,1),1)
              yb = node(cell_node_ids(i,j-9,1),2)
              xc = node(cell_node_ids(i,j-9,2),1)
              yc = node(cell_node_ids(i,j-9,2),2)
              
              ! xy components of external normal to edge
              internal(edgeId,1) = cell_data(i,idn)
              internal(edgeId,2) = cell_data(i,idn+1)
              
              ! Coordinates of center of internal edge
              internal(edgeId,3:4) = 0.5d00*(/ xb+xc , yb+yc /)
              
              ! Area of subtriangle left to internal edge
              internal(edgeId,5) = 0.5d00*ABS((xb-xa)*(yc-ya)-(xc-xa)*(yb-ya))
              
              ! Coordinates of center of cell opposite to internal edge
              xa = cell_data(id,1)
              ya = cell_data(id,2)
              
              ! Area of subtriangle right to internal edge
              internal(edgeId,6) = 0.5d00*ABS((xb-xa)*(yc-ya)-(xc-xa)*(yb-ya))
           END IF
!$omp end critical

        END IF

      END DO
    END DO
!$OMP END DO

!    nId_omp(myID) = nId_omp(myID) + 1
!    fnormal_ind_omp(myID,nId_omp(myID),4) = 2

!    DO j=0,numThreads-1 ! Loop on other columns
!      IF (j.EQ.myID) CYCLE
!      DO i=1,nId_omp(myID) ! loop on own entries
!        DO k=1,nId_omp(j) ! Loop on entries of other columns
!          IF ( fnormal_ind_omp(myID,i,4).EQ.fnormal_ind_omp(j,k,4) ) fnormal_ind_omp(myID,i,5) = 1
!        ENDDO
!      ENDDO
!    ENDDO
    ! write(*,*) "---",myID,edgeID_omp(myID),sum(edgeID_omp)
    
!$OMP END PARALLEL
    
    idx = 0
    DO j=0,numThreads-1 ! Copy contributions of other threads to a single array
      fnormal(idx+1:idx+nId_omp(j),1:5) = fnormal_omp(j,1:nId_omp(j),1:5)
      fnormal_ind(idx+1:idx+nId_omp(j),1:4) = fnormal_ind_omp(j,1:nId_omp(j),1:4)
      idx = idx + nId_omp(j)
    ENDDO
    
    ! do i=0,numThreads-1
      ! write(*,*) "***",i
       ! do j=1,edgeID_omp(i)
         ! write(*,*) (internal_ind_omp(i,j,k),k=1,2)
       ! enddo
    ! enddo
    
    ! write(*,*) "--------------",nId
    
    ! do i=1,nbrFront
      ! write(*,*) (fnormal(i,k),k=1,5)
    ! enddo
    
    ! edges(:,1:2) XY components of normal to edge
    ! edges(:,3:4) XY components of center of edge
    ! edges(:,5:6) Areas of subtriangles Lbc and Mbc located left and right to edge
    ! edges_ind(:,1:2) IDs of 2D elements located left and right of edge
    nbrInt = edgeId
        
    ALLOCATE(edges(1:nbrInt,1:6))
    ALLOCATE(edges_ind(1:nbrInt,1:2))
    
    edges     = internal(1:nbrInt,1:6)
    edges_ind = internal_ind(1:nbrInt,1:2)
        
    WRITE(*,104) idx
    WRITE(*,105) nbrInt
    
    DEALLOCATE(cell_data, internal)
    DEALLOCATE(cell_data_ind, internal_ind)
    DEALLOCATE(cell_node_ids)
    
    IF (ALLOCATED(fnormal_ind_omp)) DEALLOCATE(fnormal_ind_omp)
    IF (ALLOCATED(fnormal_omp)) DEALLOCATE(fnormal_omp)
    IF (ALLOCATED(nId_omp)) DEALLOCATE(nId_omp)
    ! IF (ALLOCATED(edgeID_omp)) DEALLOCATE(edgeID_omp)
    ! IF (ALLOCATED(internal_ind_omp)) DEALLOCATE(internal_ind_omp)
    
    CALL sampletime(time_end)
    CALL time_display
    
    104 FORMAT(" Number of boundary edges :", T35, I6)
    105 FORMAT(" Number of internal edges :", T35, I6)
    !
    WRITE(*,*) "-------------------------------------------------------"
    
END SUBROUTINE get_normal_to_cell
