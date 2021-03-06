!##########################################################
! SUBROUTINE get_normal_to_cell
!  Goal: compute the geometrical data required by the finite volume
!        method, after having read the mesh.
!##########################################################
SUBROUTINE get_normal_to_cell()
    USE module_shallow
    IMPLICIT NONE

    ! Local parameters
    INTEGER(ki) :: i,j,id,q,p,k,l,front1,front2,idx,goahead
    INTEGER(ki) :: idn,typec,edgeID=0,nID=0
    INTEGER(ki), DIMENSION(4) :: pos, pos2
    REAL(kr), DIMENSION(4,2)  :: xy
    REAL(kr) :: xa,xb,xc,xd,ya,yb,yc,yd
    LOGICAL :: test = .false.
    REAL(kr), ALLOCATABLE :: cell_data(:,:)
    REAL(kr), ALLOCATABLE :: internal(:,:)
    INTEGER(ki), ALLOCATABLE :: cell_node_ids(:,:,:)
    INTEGER(ki), ALLOCATABLE :: cell_data_ind(:,:)
    INTEGER(ki), ALLOCATABLE :: internal_ind(:,:)

    CALL sampletime(time_begin)

    ALLOCATE(cell_data(1:nbrElem,1:10))
    ALLOCATE(cell_data_ind(1:nbrElem,1:4)) ! IDs of the 4 2D elements surrounding each 2D element
    ALLOCATE(internal(1:4*nbrElem,1:6))
    ALLOCATE(internal_ind(1:4*nbrElem,1:2))
    ALLOCATE(cell_node_ids(1:nbrElem,1:4,1:2))
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
    
    ! Compute the geometrical data for the cells (geom_data) and prepare the data
    ! for the second loop (edges)
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
    END DO
    ! Save the x,y components of the normals in the matrix cell_data_n (module_shallow)
    cell_data_n(1:nbrElem,1:8) = cell_data(1:nbrElem,3:10)
        
    ! Second loop to compute the data for the internal/boundary edges
    DO i=1,nbrElem
      DO j=10,10+nbr_nodes_per_elem(i)-1
      ! id is the ID of the cell opposite to the edge, if 0 then it is a boundary edge
        id = cell_data_ind(i,j-9)
        IF (id .eq. 0) THEN ! boundary edge
             nId = nId + 1
             idn = (j-10)*2 + 3
             
             ! x,y, components of external normal to the edge
             fnormal(nId,1) = cell_data(i,idn)
             fnormal(nId,2) = cell_data(i,idn+1)
             
             ! Coordinates of center of current cell
             xa = cell_data(i,1)
             ya = cell_data(i,2)
             
             ! Coordinates of the 2 nodes of the edge
             xb = node(cell_node_ids(i,j-9,1),1)
             yb = node(cell_node_ids(i,j-9,1),2)
             xc = node(cell_node_ids(i,j-9,2),1)
             yc = node(cell_node_ids(i,j-9,2),2)
             
             ! x,y components of the center of the edge
             fnormal(nId,3:4) = 0.5d00*(/ xb+xc , yb+yc/)
             
             ! Area of subtriangle made by the border edge and the center of current cell
             fnormal(nId,5) = 0.5d00*ABS((xb-xa)*(yc-ya)-(xc-xa)*(yb-ya))
                          
             fnormal_ind(nId,1) = i
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
             fnormal_ind(nId,2) = front(idx,3)
             CALL find_vector(CLTable(1,:),5,fnormal_ind(nId,2),front1,front2)
             fnormal_ind(nId,3) = CLTable(2,front1)
             fnormal_ind(nId,4) = idx
             
        ELSE ! internal edge
           ! Search in internal_ind edges that have neighbour cell i
           CALL find_mat(internal_ind(:,1:2),i,p,pos)
           
           ! Search in internal_ind edges that have neighbour cell id
           CALL find_mat(internal_ind(:,1:2),id,q,pos2)
           
           ! Find if current edge already in internal_ind (edge that has both cells i and id as neighbours)
           do1: DO k=1,p
              DO l=1,q
                 IF (pos(k).eq.pos2(l)) THEN
                    test = .true. ! already in internal
                    EXIT do1
                 END IF
              END DO
           END DO do1

           IF (.not.test) THEN
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
           test = .false.
        END IF

      END DO
    END DO
    
    ! edges(:,1:2) XY components of normal to edge
    ! edges(:,3:4) XY components of center of edge
    ! edges(:,5:6) Areas of subtriangles Lbc and Mbc located left and right to edge
    ! edges_ind(:,1:2) IDs of 2D elements located left and right of edge
        
    edges     = internal(1:edgeId,1:6)
    edges_ind = internal_ind(1:edgeId,1:2)
    
    nbrInt = edgeId
    
    WRITE(*,104) nId
    WRITE(*,105) nbrInt

    DEALLOCATE(cell_data, internal)
    DEALLOCATE(cell_data_ind, internal_ind)
    DEALLOCATE(cell_node_ids)

    CALL sampletime(time_end)
    CALL time_display
    
    104 FORMAT(" Number of boundary edges :", T35, I6)
    105 FORMAT(" Number of internal edges :", T35, I6)
    !
    WRITE(*,*) "-------------------------------------------------------"
    
END SUBROUTINE get_normal_to_cell
