!##########################################################
! SUBROUTINE flux
!  Goal: compute the global flux defined by the finite volume method
!##########################################################
SUBROUTINE flux(Hvec,Source_sf, Q)
    USE module_shallow
    IMPLICIT NONE

    ! Subroutine parameters
    REAL(kr), DIMENSION(1:nbvar*nbrElem), INTENT(IN) :: Q
    REAL(kr), DIMENSION(1:nbvar*nbrElem), INTENT(OUT) :: Hvec
    REAL(kr), DIMENSION(1:nbvar*nbrElem), INTENT(OUT) :: Source_sf

    ! Local parameters
    REAL(kr), DIMENSION(1:nbvar) :: temp
    REAL(kr), DIMENSION(1:nbvar) :: qL, qR, Fav, Fup, qAv, sourceloc_f, sourceloc_s
    REAL(kr), DIMENSION(1:nbvar,1:2) :: FL, FR
    REAL(kr), DIMENSION(1:2) :: n
    REAL(kr) :: ds, SL, SR, h, u, v, c, Froude, dij, Hi, Hj
    INTEGER(ki) :: i, idL, idR, error
    LOGICAL  :: iswall

    Hvec = zero
    Source_sf = zero
    iswall = .FALSE.
    ! Loop on the internal edges
    internal : DO i=1,nbrInt
       temp = zero
       Fup  = zero
       idL = edges_ind(i,1) ! ID of 2D element, normal pointing outwards
       idR = edges_ind(i,2) ! ID of 2D element, normal pointing inwards
       SL = geom_data(idL,1)! area of element idL
       SR = geom_data(idR,1)! area of element idR
       qL = Q(idL*nbvar-2:idL*nbvar) ! solution h, hu, hv inside element idL
       qR = Q(idR*nbvar-2:idR*nbvar) ! solution h, hu, hv inside element idR

       IF (ABS(qL(1))<eps) qL(1)=eps
       IF (ABS(qR(1))<eps) qR(1)=eps
       qL(2:3) = qL(2:3)/qL(1)
       qR(2:3) = qR(2:3)/qR(1)

       n(1:2) = edges(i,1:2)            ! x,y components of the external normal to the edge
       ds = SQRT(n(1)*n(1) + n(2)*n(2)) ! length of the edge
       n = n/ds                         ! normalize the normal
       
       CALL getFlux(qL,FL) ! Get the local flux for the cell left of the edge
       CALL getFlux(qR,FR) ! Get the local flux for the cell right of the edge
       
       ! Average the flux
       Fav = ((SL*FL(:,1) + SR*FR(:,1))*n(1) + (SL*FL(:,2) + SR*FR(:,2))*n(2))/(SL+SR)
       qAv = (SL*qL + SR*qR)/(SL+SR) ! Average the solution 
       
       ! Get the distance from center of element to center of edge
       dij = SQRT( ( geom_data(idL,3) - edges(i,3) )**2 + ( geom_data(idL,4) - edges(i,4) )**2 )
       Hi = depth(idL) ! geometric depth, cell left to the edge
       Hj = depth(idR) ! geometric depth, cell right to the edge
       
       ! Get the upwind + source terms
       CALL getUpwind_and_Source(qL,qR,qAv,n,Fup,sourceloc_f,sourceloc_s,Hi,Hj,dij,iswall,ds,SL)
       
       ! Add the contribution to the element idL
       temp = (Fav + Fup)*ds/SL
       Hvec     (idL*nbvar-2:idL*nbvar) = Hvec     (idL*nbvar-2:idL*nbvar) + temp
       Source_sf(idL*nbvar-2:idL*nbvar) = Source_sf(idL*nbvar-2:idL*nbvar) + (sourceloc_f+sourceloc_s) * edges(i,5) / SL
       
       ! Add the contribution to the other element idR, Fav -> -Fav
       n = -n
       
       ! Get the distance from center of element to center of edge
       dij = SQRT( ( geom_data(idR,3) - edges(i,3) )**2 + ( geom_data(idR,4) - edges(i,4) )**2 )
       
       ! Get the upwind + source terms
       CALL getUpwind_and_Source(qR,qL,qAv,n,Fup,sourceloc_f,sourceloc_s,Hj,Hi,dij,iswall,ds,SR)
       
       ! Add the contribution to the element idR
       temp = (- Fav + Fup)*ds/SR
       Hvec     (idR*nbvar-2:idR*nbvar) = Hvec     (idR*nbvar-2:idR*nbvar) + temp
       Source_sf(idR*nbvar-2:idR*nbvar) = Source_sf(idR*nbvar-2:idR*nbvar) + (sourceloc_f+sourceloc_s) * edges(i,6) / SR

    END DO internal

    ! Loop on the border edges to take into account the boundary conditions
    border : DO i=1,nbrFront 
       temp = zero
       Fup  = zero
       idL = fnormal_ind(i,1) ! ID of 2D element
       SL = geom_data(idL,1)  ! area of element idL
       qL = Q(idL*nbvar-2:idL*nbvar) ! solution h, hu, hv inside element idL
       IF (ABS(qL(1))<eps) qL(1)=eps
       qL(2:3) = qL(2:3)/qL(1)

       n(1:2) = fnormal(i,1:2)      ! x,y components of the external normal to the edge
       ds = SQRT(n(1)**2 + n(2)**2) ! length of the edge
       n = n/ds                     ! normalize the normal
       iswall = .FALSE.
       
       ! Select the type of boundary condition to impose
       ! Fill the array qR for the ghost cell
       SELECT CASE (fnormal_ind(i,3))
          CASE (0) ! inlet
             h = BoundCond(fnormal_ind(i,4),1) ! imposed h
             u = BoundCond(fnormal_ind(i,4),2) ! imposed u
             v = BoundCond(fnormal_ind(i,4),3) ! imposed v
             c = SQRT(ggrav*h) ! wave speed
             Froude = SQRT(u*u+v*v)/c
             IF(Froude.GE.1.0D00) THEN ! supercritical flow
                qR(1) = 2*h - qL(1)
                qR(2) = 2*u - qL(2)
                qR(3) = 2*v - qL(3)
             ELSE ! subcritical flow
                qR(1) = qL(1)
                
                ! Force to keep the prescribed unit discharge q = u qL(1)
                u = u * h / qL(1)
                v = v * h / qL(1)
                
                qR(2) = 2*u - qL(2)
                qR(3) = 2*v - qL(3)
             ENDIF
             
          CASE (1) ! outlet
             h = BoundCond(fnormal_ind(i,4),1) ! imposed h
             u = BoundCond(fnormal_ind(i,4),2) ! imposed u
             v = BoundCond(fnormal_ind(i,4),3) ! imposed v
             c = SQRT(ggrav*h) ! wave speed
             Froude = SQRT(u*u+v*v)/c
             IF(Froude>=1) THEN ! supercritical flow
                qR = qL
             ELSE ! subcritical
                qR(1)   = 2*h - qL(1)
                qR(2:3) = qL(2:3)
             ENDIF

          CASE (2) ! wall, cancel the wall-normal velocity component
             qR(1) = qL(1)
             qR(2) = qL(2) - 2*(qL(2)*n(1)+qL(3)*n(2))*n(1)
             qR(3) = qL(3) - 2*(qL(2)*n(1)+qL(3)*n(2))*n(2)
             iswall = .TRUE.
     
          CASE (3) ! symmetry
             qR = qL
             qR(2:3) = -qR(2:3)
             
          CASE (4) ! periodicity
             qR = qL

          CASE DEFAULT
             WRITE(*,*) "Error : unrecognized boundary condition ",fnormal_ind(i,3)
       END SELECT
       
       CALL getFlux(qL,FL) ! Get the local flux for the cell left of the edge
       CALL getFlux(qR,FR) ! Get the local flux for the ghost cell
       
       ! Average the flux and the solution
       Fav = 0.5*(FL(:,1)+FR(:,1))*n(1) + 0.5*(FL(:,2)+FR(:,2))*n(2)
       qAv = 0.5*(qL+qR)
       
       ! Get the distance from center of element to center of edge
       dij = SQRT( ( geom_data(idL,3) - fnormal(i,3) )**2 + ( geom_data(idL,4) - fnormal(i,4) )**2 )
       Hi = depth(idL)
       Hj = depth(idL)
       
       ! Get the upwind term. The source term is not relevant for border edges
       CALL getUpwind_and_Source(qL,qR,qAv,n,Fup,sourceloc_f,sourceloc_s,Hi,Hj,dij,iswall,ds,SL)
       
       ! Add the contribution to the element idL
       temp = (Fav + Fup)*ds/SL
       Hvec(idL*nbvar-2:idL*nbvar) = Hvec(idL*nbvar-2:idL*nbvar) + temp
       Source_sf(idL*nbvar-2:idL*nbvar) = Source_sf(idL*nbvar-2:idL*nbvar) + sourceloc_f * fnormal(i,5) / SL
    END DO border

END SUBROUTINE flux

!##########################################################
! SUBROUTINE getFlux
!  Goal: compute the local flux defined by the finite volume method
!      (h u)                           (h v)
! Fx = (h u^2 + 0.5 g h^2)        Fy = (h u v)
!      (h u v)                         (h v^2 + 0.5 g h^2)
!##########################################################
SUBROUTINE getFlux(q,F)
    USE module_shallow
    IMPLICIT NONE

    ! Subroutine parameters
    REAL(kr), DIMENSION(nbvar), INTENT(IN) :: q
    REAL(kr), DIMENSION(1:nbvar,1:2), INTENT(OUT) :: F
    
    ! Flux for the X-gradient
    F(1,1) = q(1)*q(2)
    F(2,1) = q(1)*q(2)*q(2) + 0.5*ggrav*q(1)*q(1)
    F(3,1) = q(1)*q(2)*q(3)

    ! Flux for the Y-gradient
    F(1,2) = q(1)*q(3)
    F(2,2) = q(1)*q(2)*q(3)
    F(3,2) = q(1)*q(3)*q(3) + 0.5*ggrav*q(1)*q(1)

END SUBROUTINE getFlux

!##########################################################
! SUBROUTINE getUpwind_and_Source
!  Goal: compute the upwind correction flux defined by the finite volume method
!        for both the flux and source terms
!##########################################################
SUBROUTINE getUpwind_and_Source(qL,qR,qAvg,n,Fout,source_f,source_s,Hi,Hj,dij,iswall,ds,omega)
    USE module_shallow, only : kr,ki,ggrav,manning_b,manning_w,nbvar,zero
    USE booklib
    IMPLICIT NONE
    
    ! Subroutine parameters
    REAL(kr), DIMENSION(2), INTENT(IN)        :: n
    REAL(kr), DIMENSION(1:nbvar), INTENT(IN)  :: qL, qR, qAvg
    REAL(kr), DIMENSION(1:nbvar), INTENT(OUT) :: Fout,source_f,source_s
    REAL(kr), INTENT(IN)                      :: Hi, Hj, dij, ds,omega
    LOGICAL, INTENT(IN)                       :: iswall
    
    ! Local parameters
    REAL(kr), DIMENSION(1:nbvar,1:nbvar) :: Qmat,sourceMat
    REAL(kr), DIMENSION(1:nbvar)         :: qLloc, qRloc,tmp
    REAL(kr)                             :: tempval, locdist
    
    source_f = zero
    source_s = zero
    
    qLloc = QL
    qRloc = QR
    ! Set again hu,hv as the variables
    qLloc(2:3) = qLloc(2:3) * qLloc(1)
    qRloc(2:3) = qRloc(2:3) * qRloc(1)
    
    ! Get the matrices    abs(Q) = X abs(lambda) X-1
    !  and             sourceMat = X ( I - abs(lambda) lambda-1 ) X-1
    CALL upwind_term(qAvg,Qmat,sourceMat,n)
    
    ! Get the upwind correction term for the flux
    tmp = MATMUL( Qmat , qRloc - qLloc )
    Fout = - 0.5d00 * tmp
    
    ! Get the upwind corrected source term (bed slope)
    tmp = zero
    tempval = - 0.5d00 * ggrav * ( QL(1)+QR(1) ) * ( Hj - Hi ) / dij
    tmp(2) = tempval*n(1)
    tmp(3) = tempval*n(2)
    source_s = MATMUL( sourceMat , tmp )
    
    ! Friction slope x: -g h n^2 u(u^2 + v^2)^(1/2) / h^(4/3)
    !                y: -g h n^2 v(u^2 + v^2)^(1/2) / h^(4/3)
    tempval = - ggrav * QL(1) * sqrt(ql(2)**2 + qR(2)**2)
    locdist = (manning_b)**(1.5D00) / QL(1) 
    IF (iswall) THEN
      locdist = locdist + (manning_w)**(1.5D0) * ds / omega
    ENDIF
    source_f(2) = tempval * QL(2) * ( locdist**(4.D00 / 3.D00) )
    source_f(3) = tempval * QL(3) * ( locdist**(4.D00 / 3.D00) )
    
    ! source_f(2) = - ggrav * QL(1) * manning_b*manning_b * QL(2) * sqrt(ql(2)**2 + qR(2)**2) / (QL(1)**(4.d00/3.d00))
    ! source_f(3) = - ggrav * QL(1) * manning_b*manning_b * QL(3) * sqrt(ql(2)**2 + qR(2)**2) / (QL(1)**(4.d00/3.d00))
    
    ! IF (iswall) THEN
      ! source_f(2) = source_f(2) - ggrav * QL(1) * manning_w*manning_w * QL(2) * sqrt(ql(2)**2 + qR(2)**2) / ((omega/ds)**(4.d00/3.d00))
      ! source_f(3) = source_f(3) - ggrav * QL(1) * manning_w*manning_w * QL(3) * sqrt(ql(2)**2 + qR(2)**2) / ((omega/ds)**(4.d00/3.d00))
    ! ENDIF
    
END SUBROUTINE getUpwind_and_Source

!##########################################################
! SUBROUTINE upwind_term
!  Goal: get the matrices   abs(Q) = X abs(lambda) X-1
!                        sourceMat = X ( I - abs(lambda) lambda-1 ) X-1
!##########################################################
SUBROUTINE upwind_term(q,Qmat,sourceMat,n)
    USE module_shallow
    USE booklib
    IMPLICIT NONE
    
    ! Subroutine parameters
    REAL(kr), DIMENSION(2), INTENT(IN) :: n
    REAL(kr), DIMENSION(1:nbvar), INTENT(IN) :: q
    REAL(kr), DIMENSION(1:nbvar,1:nbvar), INTENT(OUT) :: Qmat
    REAL(kr), DIMENSION(1:nbvar,1:nbvar), INTENT(OUT) :: sourceMat
    
    ! Local parameters
    REAL(kr) :: one = 1.0d00
    REAL(kr), DIMENSION(1:nbvar,1:nbvar) :: X, Xm1, abslambda,invlambda,identity, temp, temp2
    REAL(kr) :: h, u, v, c,lambda1, lambda2, lambda3
    
    h = q(1)
    u = q(2)
    v = q(3)
    c = SQRT(ggrav*h) ! wave speed

    lambda1 = u*n(1)+v*n(2)      ! first eigenvalue
    lambda2 = lambda1      + c   ! second eigenvalue
    lambda3 = lambda1      - c   ! third eigenvalue
    
    ! Fill in the matrix that contains the absolute eigenvalues on the diagonal
    abslambda = zero
    abslambda(1,1) = ABS(lambda1)
    abslambda(2,2) = ABS(lambda2)
    abslambda(3,3) = ABS(lambda3)
    
    ! Fill in the matrix that contains the inverse of the eigenvalues on the diagonal
    invlambda = zero
    IF (ABS(lambda1).EQ.0.0d0) lambda1=eps
    IF (ABS(lambda2).EQ.0.0d0) lambda2=eps
    IF (ABS(lambda3).EQ.0.0d0) lambda3=eps
    invlambda(1,1) = one/lambda1
    invlambda(2,2) = one/lambda2
    invlambda(3,3) = one/lambda3
    
    ! Identity matrix
    identity = zero
    identity(1,1) = one
    identity(2,2) = one
    identity(3,3) = one
    
    ! Row eigenvector matrix
    X(1,1:3) = (/   zero,      one,     one /)
    X(2,1:3) = (/-c*n(2), u+c*n(1), u-c*n(1)/)
    X(3,1:3) = (/ c*n(1), v+c*n(2), v-c*n(2)/)
    
    ! Column eigenvector matrix
    Xm1(1,1:3) = (/   2*n(2)*u-2*n(1)*v, -2*n(2),  2*n(1)     /)
    Xm1(2,1:3) = (/   c-n(1)*u-n(2)*v,      n(1),    n(2)     /)
    Xm1(3,1:3) = (/   c+n(1)*u+n(2)*v,     -n(1),   -n(2)     /)
    Xm1 = Xm1 * 0.5d0 / c
    
    ! Compute abs(Q) = X abs(lambda) X-1
    temp  = MATMUL(abslambda,Xm1)
    Qmat  = MATMUL(X,temp)
    
    ! Compute sourceMat = X ( I - abs(lambda) lambda-1 ) X-1
    temp       = MATMUL(abslambda,invlambda)
    temp       = identity - temp
    temp2      = MATMUL(temp,Xm1)
    sourceMat  = MATMUL(X,temp2)
    
END SUBROUTINE upwind_term

