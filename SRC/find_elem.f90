!##########################################################
! SUBROUTINE find_elem_front
!  Goal: get the ID of the 2D element linked to a boundary edge elem_front
!##########################################################
SUBROUTINE find_elem_front(elem_front,elem_surf)
    USE module_shallow
    IMPLICIT NONE

    ! Subroutine parameters
    INTEGER(ki), INTENT(IN) :: elem_front
    INTEGER(ki), INTENT(OUT):: elem_surf

    ! Local parameters
    INTEGER(ki) :: i,node1, node2,ind
    INTEGER(ki), DIMENSION(10) :: temp
    
    ind = 0
    ! The 2 nodes that define the boundary edge elem_front
    node1 = front(elem_front,1)
    node2 = front(elem_front,2)
    
    ! Find all the elements that contain the first node node1 and store them in temp
    DO i=1,nbrElem
      IF ((elem(i,1).eq.node1).or.(elem(i,2).eq.node1) .or. (elem(i,3).eq.node1) .or. (elem(i,4).eq.node1)) then
        ind = ind + 1
        temp(ind) = i
      END IF
    END DO
    
    ! Among the previously found elements, find the boundary 2D element that contains the second node
    DO i=1,ind
      IF ((elem(temp(i),1)==node2) .or. (elem(temp(i),2)==node2) .or. (elem(temp(i),3)==node2) .or. (elem(temp(i),4)==node2)) then
        elem_surf = temp(i)
        EXIT
      END IF
    END DO

END SUBROUTINE find_elem_front

!##########################################################
! SUBROUTINE find_elem_surf
!  Goal: get the ID (elem_surf) of the 2D element located on the 
!        other side of the edge of another 2D element (elem_cur)
!##########################################################
SUBROUTINE find_elem_surf(elem_cur,elem_surf,node1,node2)
    USE module_shallow
    IMPLICIT NONE

    ! Subroutine parameters
    INTEGER(ki), INTENT(IN) :: elem_cur, node1, node2
    INTEGER(ki), INTENT(OUT):: elem_surf

    ! Local parameters
    INTEGER(ki) :: i,ind,ind2
    INTEGER(ki), DIMENSION(10) :: temp, temp2

    ind = 0
    ind2 = 0

    ! Find all the elements that contain the first node node1 and store them in temp
    DO i=1,nbrElem
      IF ((elem(i,1).eq.node1).or.(elem(i,2).eq.node1) .or. (elem(i,3).eq.node1) .or. (elem(i,4).eq.node1)) THEN
        ind = ind + 1
        temp(ind) = i
      END IF
    END DO

    ! Among the previously found elements, find the 2D elements that contains the second node
    DO i=1,ind
      IF ((elem(temp(i),1)==node2) .or. (elem(temp(i),2)==node2) .or. (elem(temp(i),3)==node2) .or. (elem(temp(i),4)==node2)) THEN
        ind2 = ind2 + 1
        temp2(ind2) = temp(i)
      END IF
    END DO
    
    ! Set the output
    IF (ind2 == 1) THEN
      elem_surf = 0
    ELSE IF (elem_cur == temp2(1)) THEN
      elem_surf = temp2(2)
    ELSE 
      elem_surf = temp2(1)
    END IF

END SUBROUTINE find_elem_surf

!##########################################################
! SUBROUTINE find_vector
!  Goal: find index 'element' inside array 'vector'
!##########################################################
SUBROUTINE find_vector(vector,element,front1,front2)
    USE module_shallow
    IMPLICIT NONE

    ! Subroutine parameters
    INTEGER(ki), DIMENSION(nbrFront), INTENT(IN) :: vector
    INTEGER(ki), INTENT(IN) :: element
    INTEGER(ki), INTENT(OUT) :: front1,front2
    
    ! Local parameters
    INTEGER(ki) :: i, ind
    
    front1 = 0
    front2 = 0
    ind = 0

    loop : DO i=1,nbrFront
      IF (vector(i) .EQ. element) THEN
        ind = ind + 1
        IF (ind .EQ. 1) THEN
           front1 = i
        ELSE
           front2 = i
           EXIT loop
        END IF
      END IF
    END DO loop

END SUBROUTINE find_vector

!##########################################################
! SUBROUTINE find_mat
!  Goal: find an index inside a 2D array
!##########################################################
SUBROUTINE find_mat(mat,element,ind,pos)
    USE module_shallow
    IMPLICIT NONE

    ! Subroutine parameters
    INTEGER(ki), DIMENSION(4*nbrElem,2), INTENT(IN) :: mat
    INTEGER(ki), INTENT(IN) :: element
    INTEGER(ki), INTENT(OUT) :: ind
    INTEGER(ki), DIMENSION(4), INTENT(OUT) :: pos
    
    ! Local parameters
    INTEGER(ki) :: i
    
    ind = 0 ! number of times that element appears in mat
    pos = (/0,0,0,0/)

    DO i=1,4*nbrElem
       IF ((mat(i,1).eq.element).or.(mat(i,2).eq.element)) THEN
          ind = ind + 1
          pos(ind) = i
       END IF
       IF ((mat(i,1).eq.0).or.(mat(i,2).eq.0)) EXIT
    END DO

END SUBROUTINE find_mat
