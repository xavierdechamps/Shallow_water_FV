MODULE booklib
IMPLICIT NONE
save

! Set access--private by default, and public for specific procedures.
PRIVATE
PUBLIC mat_inv                    ! Matrix inversion
PUBLIC simul                      ! Solve simultaneous eqns

! Declare generic procedures.

INTERFACE mat_inv
   MODULE PROCEDURE mat_inv_sgl
   MODULE PROCEDURE mat_inv_dbl
   MODULE PROCEDURE mat_inv_sgl_cmplx
   MODULE PROCEDURE mat_inv_dbl_cmplx
END INTERFACE

INTERFACE simul
   MODULE PROCEDURE simul_sgl
   MODULE PROCEDURE simul_dbl
   MODULE PROCEDURE simul_sgl_cmplx
   MODULE PROCEDURE simul_dbl_cmplx
END INTERFACE


! Declare module procedures.
CONTAINS


SUBROUTINE mat_inv_sgl ( a, b, ndim, n, error )
!
!  Purpose:
!    Subroutine to n x n matrix using Gaussian elimination
!    and the maximum pivot technique.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/16/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=6) ! Precision
REAL(KIND=kind), PARAMETER :: eps = 10. * EPSILON(0.0_kind)
                                     ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
REAL(KIND=kind),INTENT(IN),DIMENSION(ndim,ndim) :: a
                                     ! Input matrix (N x N).  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
REAL(KIND=kind),INTENT(OUT),DIMENSION(ndim,ndim) :: b
                                     ! Inverse of matrix a.  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local variables:
REAL(KIND=kind),DIMENSION(n,n) :: a1 ! Copy of a to destroy while
                                     ! building the inverse.  Only
                                     ! actual elements in use are
                                     ! duplicated.
REAL(KIND=kind),DIMENSION(n,n) :: b1 ! Array in which to build
                                     ! inverse.  Only acutal elements
                                     ! in use are duplicated.
REAL(KIND=kind) :: factor            ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
REAL(KIND=kind),DIMENSION(n) :: temp ! Scratch array

! Make a copy of the input array.
a1 = a(1:n,1:n)

! Initialize the inverse array.
b1 = 0.
DO irow = 1, n
   b1(irow,irow) = 1.
END DO

! Process n times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to n
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < eps ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp = a1(ipeak,:)
      a1(ipeak,:) = a1(irow,:)     ! Swap rows in a1
      a1(irow,:) = temp
      temp = b1(ipeak,:)
      b1(ipeak,:) = b1(irow,:)     ! Swap rows in b1
      b1(irow,:) = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow), and
   ! add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,:) = a1(irow,:)*factor + a1(jrow,:)
         b1(jrow,:) = b1(irow,:)*factor + b1(jrow,:)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   b1(irow,:) = b1(irow,:) / a1(irow,irow)
   a1(irow,irow) = 1.
END DO divide

! Copy the answer to the output, set error flag
! to 0 and return.
b = 0.
b(1:n,1:n) = b1
error = 0

END SUBROUTINE mat_inv_sgl

SUBROUTINE mat_inv_dbl ( a, b, ndim, n, error )
!
!  Purpose:
!    Subroutine to n x n matrix using Gaussian elimination
!    and the maximum pivot technique.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/16/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=12) ! Precision
REAL(KIND=kind), PARAMETER :: eps = 10. * EPSILON(0.0_kind)
                                     ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
REAL(KIND=kind),INTENT(IN),DIMENSION(ndim,ndim) :: a
                                     ! Input matrix (N x N).  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
REAL(KIND=kind),INTENT(OUT),DIMENSION(ndim,ndim) :: b
                                     ! Inverse of matrix a.  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local variables:
REAL(KIND=kind),DIMENSION(n,n) :: a1 ! Copy of a to destroy while
                                     ! building the inverse.  Only
                                     ! actual elements in use are
                                     ! duplicated.
REAL(KIND=kind),DIMENSION(n,n) :: b1 ! Array in which to build
                                     ! inverse.  Only acutal elements
                                     ! in use are duplicated.
REAL(KIND=kind) :: factor            ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
REAL(KIND=kind),DIMENSION(n) :: temp ! Scratch array

! Make a copy of the input array.
a1 = a(1:n,1:n)

! Initialize the inverse array.
b1 = 0.
DO irow = 1, n
   b1(irow,irow) = 1.
END DO

! Process n times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to n
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < eps ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp = a1(ipeak,:)
      a1(ipeak,:) = a1(irow,:)     ! Swap rows in a1
      a1(irow,:) = temp
      temp = b1(ipeak,:)
      b1(ipeak,:) = b1(irow,:)     ! Swap rows in b1
      b1(irow,:) = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow), and
   ! add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,:) = a1(irow,:)*factor + a1(jrow,:)
         b1(jrow,:) = b1(irow,:)*factor + b1(jrow,:)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   b1(irow,:) = b1(irow,:) / a1(irow,irow)
   a1(irow,irow) = 1.
END DO divide

! Copy the answer to the output, set error flag
! to 0 and return.
b = 0.
b(1:n,1:n) = b1
error = 0

END SUBROUTINE mat_inv_dbl

SUBROUTINE mat_inv_sgl_cmplx ( a, b, ndim, n, error )
!
!  Purpose:
!    Subroutine to n x n matrix using Gaussian elimination
!    and the maximum pivot technique.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/16/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=6) ! Precision
REAL(KIND=kind), PARAMETER :: eps = 10. * EPSILON(0.0_kind)
                                     ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
COMPLEX(KIND=kind),INTENT(IN),DIMENSION(ndim,ndim) :: a
                                     ! Input matrix (N x N).  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
COMPLEX(KIND=kind),INTENT(OUT),DIMENSION(ndim,ndim) :: b
                                     ! Inverse of matrix a.  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local variables:
COMPLEX(KIND=kind),DIMENSION(n,n) :: a1 ! Copy of a to destroy while
                                     ! building the inverse.  Only
                                     ! actual elements in use are
                                     ! duplicated.
COMPLEX(KIND=kind),DIMENSION(n,n) :: b1 ! Array in which to build
                                     ! inverse.  Only acutal elements
                                     ! in use are duplicated.
COMPLEX(KIND=kind) :: factor         ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
COMPLEX(KIND=kind),DIMENSION(n) :: temp ! Scratch array

! Make a copy of the input array.
a1 = a(1:n,1:n)

! Initialize the inverse array.
b1 = 0.
DO irow = 1, n
   b1(irow,irow) = (1.0_kind,0.0_kind)
END DO

! Process n times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to n
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < eps ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp = a1(ipeak,:)
      a1(ipeak,:) = a1(irow,:)     ! Swap rows in a1
      a1(irow,:) = temp
      temp = b1(ipeak,:)
      b1(ipeak,:) = b1(irow,:)     ! Swap rows in b1
      b1(irow,:) = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow), and
   ! add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,:) = a1(irow,:)*factor + a1(jrow,:)
         b1(jrow,:) = b1(irow,:)*factor + b1(jrow,:)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   b1(irow,:) = b1(irow,:) / a1(irow,irow)
   a1(irow,irow) = 1.
END DO divide

! Copy the answer to the output, set error flag
! to 0 and return.
b = 0.
b(1:n,1:n) = b1
error = 0

END SUBROUTINE mat_inv_sgl_cmplx

SUBROUTINE mat_inv_dbl_cmplx ( a, b, ndim, n, error )
!
!  Purpose:
!    Subroutine to n x n matrix using Gaussian elimination
!    and the maximum pivot technique.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    03/16/96    S. J. Chapman        Original code
!
IMPLICIT NONE

! Declare local parameters
INTEGER, PARAMETER :: kind = SELECTED_REAL_KIND(p=12) ! Precision
REAL(KIND=kind), PARAMETER :: eps = 10. * EPSILON(0.0_kind)
                                     ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
COMPLEX(KIND=kind),INTENT(IN),DIMENSION(ndim,ndim) :: a
                                     ! Input matrix (N x N).  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
COMPLEX(KIND=kind),INTENT(OUT),DIMENSION(ndim,ndim) :: b
                                     ! Inverse of matrix a.  This
                                     ! array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
                                     ! The declared dimension ndim
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local variables:
COMPLEX(KIND=kind),DIMENSION(n,n) :: a1 ! Copy of a to destroy while
                                     ! building the inverse.  Only
                                     ! actual elements in use are
                                     ! duplicated.
COMPLEX(KIND=kind),DIMENSION(n,n) :: b1 ! Array in which to build
                                     ! inverse.  Only acutal elements
                                     ! in use are duplicated.
COMPLEX(KIND=kind) :: factor         ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
COMPLEX(KIND=kind),DIMENSION(n) :: temp ! Scratch array

! Make a copy of the input array.
a1 = a(1:n,1:n)

! Initialize the inverse array.
b1 = 0.
DO irow = 1, n
   b1(irow,irow) = (1.0_kind,0.0_kind)
END DO

! Process n times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to n
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < eps ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp = a1(ipeak,:)
      a1(ipeak,:) = a1(irow,:)     ! Swap rows in a1
      a1(irow,:) = temp
      temp = b1(ipeak,:)
      b1(ipeak,:) = b1(irow,:)     ! Swap rows in b1
      b1(irow,:) = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow), and
   ! add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,:) = a1(irow,:)*factor + a1(jrow,:)
         b1(jrow,:) = b1(irow,:)*factor + b1(jrow,:)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   b1(irow,:) = b1(irow,:) / a1(irow,irow)
   a1(irow,irow) = 1.
END DO divide

! Copy the answer to the output, set error flag
! to 0 and return.
b = 0.
b(1:n,1:n) = b1
error = 0

END SUBROUTINE mat_inv_dbl_cmplx


SUBROUTINE simul_sgl ( a, b, soln, ndim, n, error )
!
!  Purpose:
!    Subroutine to solve a set of N linear equations in N
!    unknowns using Gaussian elimination and the maximum
!    pivot technique.  This version of simul has been
!    modified to use array sections and automatic arrays.
!    It DOES NOT DESTROY the original input values.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    10/16/95    S. J. Chapman        Original code
! 1. 10/17/95    S. J. Chapman        Modified to use array ops
! 2. 03/02/96    S. J. Chapman        Add allocatable arrays
!
IMPLICIT NONE

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
REAL, INTENT(IN), DIMENSION(ndim,ndim) :: a
                                     ! Array of coefficients (N x N).
                                     ! This array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
REAL, INTENT(IN), DIMENSION(ndim) :: b
                                     ! Input: Right-hand side of eqns.
REAL, INTENT(OUT), DIMENSION(ndim) :: soln
                                     ! Output: Solution vector.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local parameters
REAL, PARAMETER :: epsilon = 1.0E-6  ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare local variables:
REAL, DIMENSION(n,n) :: a1           ! Copy of "a" which will be
                                     ! destroyed during the solution
REAL :: factor                       ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
REAL :: temp                         ! Scratch value
REAL, DIMENSION(n) :: temp1          ! Scratch array

! Make copies of arrays "a" and "b" for local use
a1 = a(1:n,1:n)
soln = b(1:n)

! Process N times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to N
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < epsilon ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp1 = a1(ipeak,1:n)
      a1(ipeak,1:n) = a1(irow,1:n)   ! Swap rows in a
      a1(irow,1:n) = temp1
      temp = soln(ipeak)
      soln(ipeak) = soln(irow)       ! Swap rows in b
      soln(irow)  = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow),
   ! and add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,1:n) = a1(irow,1:n)*factor + a1(jrow,1:n)
         soln(jrow) = soln(irow)*factor + soln(jrow)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   soln(irow) = soln(irow) / a1(irow,irow)
END DO divide

! Set error flag to 0 and return.
error = 0

END SUBROUTINE simul_sgl

SUBROUTINE simul_dbl ( a, b, soln, ndim, n, error )
!
!  Purpose:
!    Subroutine to solve a set of N linear equations in N
!    unknowns using Gaussian elimination and the maximum
!    pivot technique.  This version of simul has been
!    modified to use array sections and automatic arrays.
!    It uses double precision arithmetic to avoid
!    cumulative roundoff errors.  It DOES NOT DESTROY the
!    original input values.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    10/16/95    S. J. Chapman        Original code
! 1. 10/17/95    S. J. Chapman        Modified to use array ops
! 2. 03/02/96    S. J. Chapman        Add allocatable arrays
! 3. 05/08/96    S. J. Chapman        Double precision
!
IMPLICIT NONE

! Declare parameters
INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(p=12)

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
REAL(KIND=dbl), INTENT(IN), DIMENSION(ndim,ndim) :: a
                                     ! Array of coefficients (N x N).
                                     ! This array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
REAL(KIND=dbl), INTENT(IN), DIMENSION(ndim) :: b
                                     ! Input: Right-hand side of eqns.
REAL(KIND=dbl), INTENT(OUT), DIMENSION(ndim) :: soln
                                     ! Output: Solution vector.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local parameters
REAL(KIND=dbl), PARAMETER :: epsilon = 1.0E-12
                                     ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare local variables:
REAL(KIND=dbl), DIMENSION(n,n) :: a1 ! Copy of "a" which will be
                                     ! destroyed during the solution
REAL(KIND=dbl) :: factor             ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
REAL(KIND=dbl) :: temp               ! Scratch value
REAL(KIND=dbl),DIMENSION(n) :: temp1 ! Scratch array

! Make copies of arrays "a" and "b" for local use
a1 = a(1:n,1:n)
soln = b(1:n)

! Process N times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to N
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < epsilon ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp1 = a1(ipeak,1:n)
      a1(ipeak,1:n) = a1(irow,1:n)   ! Swap rows in a
      a1(irow,1:n) = temp1
      temp = soln(ipeak)
      soln(ipeak) = soln(irow)       ! Swap rows in b
      soln(irow)  = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow),
   ! and add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,1:n) = a1(irow,1:n)*factor + a1(jrow,1:n)
         soln(jrow) = soln(irow)*factor + soln(jrow)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   soln(irow) = soln(irow) / a1(irow,irow)
END DO divide

! Set error flag to 0 and return.
error = 0

END SUBROUTINE simul_dbl

SUBROUTINE simul_sgl_cmplx ( a, b, soln, ndim, n, error )
!
!  Purpose:
!    Subroutine to solve a set of N linear equations in N
!    unknowns using Gaussian elimination and the maximum
!    pivot technique.  This version of simul has been
!    modified to use array sections and automatic arrays.
!    It uses single precision complex arithmetic.  It DOES
!    NOT DESTROY the original input values.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    10/16/95    S. J. Chapman        Original code
! 1. 10/17/95    S. J. Chapman        Modified to use array ops
! 2. 03/02/96    S. J. Chapman        Add allocatable arrays
! 3. 05/08/96    S. J. Chapman        Double precision
! 4. 05/09/96    S. J. Chapman        Complex version
!
IMPLICIT NONE

! Declare parameters
INTEGER, PARAMETER :: sgl = SELECTED_REAL_KIND(p=6)

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
COMPLEX(KIND=sgl), INTENT(IN), DIMENSION(ndim,ndim) :: a
                                     ! Array of coefficients (N x N).
                                     ! This array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
COMPLEX(KIND=sgl), INTENT(IN), DIMENSION(ndim) :: b
                                     ! Input: Right-hand side of eqns.
COMPLEX(KIND=sgl), INTENT(OUT), DIMENSION(ndim) :: soln
                                     ! Output: Solution vector.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local parameters
REAL(KIND=sgl), PARAMETER :: epsilon = 1.0E-6
                                     ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare local variables:
REAL(KIND=sgl), DIMENSION(n,n) :: a1 ! Copy of "a" which will be
                                     ! destroyed during the solution
REAL(KIND=sgl) :: factor             ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
COMPLEX(KIND=sgl) :: temp            ! Scratch value
COMPLEX(KIND=sgl),DIMENSION(n) :: temp1 ! Scratch array

! Make copies of arrays "a" and "b" for local use
a1 = a(1:n,1:n)
soln = b(1:n)

! Process N times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to N
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < epsilon ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp1 = a1(ipeak,1:n)
      a1(ipeak,1:n) = a1(irow,1:n)   ! Swap rows in a
      a1(irow,1:n) = temp1
      temp = soln(ipeak)
      soln(ipeak) = soln(irow)       ! Swap rows in b
      soln(irow)  = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow),
   ! and add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,1:n) = a1(irow,1:n)*factor + a1(jrow,1:n)
         soln(jrow) = soln(irow)*factor + soln(jrow)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   soln(irow) = soln(irow) / a1(irow,irow)
END DO divide

! Set error flag to 0 and return.
error = 0

END SUBROUTINE simul_sgl_cmplx

SUBROUTINE simul_dbl_cmplx ( a, b, soln, ndim, n, error )
!
!  Purpose:
!    Subroutine to solve a set of N linear equations in N
!    unknowns using Gaussian elimination and the maximum
!    pivot technique.  This version of simul has been
!    modified to use array sections and automatic arrays.
!    It uses double precision complex arithmetic.  It DOES
!    NOT DESTROY the original input values.
!
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!    10/16/95    S. J. Chapman        Original code
! 1. 10/17/95    S. J. Chapman        Modified to use array ops
! 2. 03/02/96    S. J. Chapman        Add allocatable arrays
! 3. 05/08/96    S. J. Chapman        Double precision
! 4. 05/09/96    S. J. Chapman        Double complex version
!
IMPLICIT NONE

! Declare parameters
INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(p=12)

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
COMPLEX(KIND=dbl), INTENT(IN), DIMENSION(ndim,ndim) :: a
                                     ! Array of coefficients (N x N).
                                     ! This array is of size ndim x
                                     ! ndim, but only N x N of the
                                     ! coefficients are being used.
COMPLEX(KIND=dbl), INTENT(IN), DIMENSION(ndim) :: b
                                     ! Input: Right-hand side of eqns.
COMPLEX(KIND=dbl), INTENT(OUT), DIMENSION(ndim) :: soln
                                     ! Output: Solution vector.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations

! Declare local parameters
REAL(KIND=dbl), PARAMETER :: epsilon = 1.0E-12
                                     ! A "small" number for comparison
                                     ! when determining singular eqns

! Declare local variables:
REAL(KIND=dbl), DIMENSION(n,n) :: a1 ! Copy of "a" which will be
                                     ! destroyed during the solution
REAL(KIND=dbl) :: factor             ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! currently being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
COMPLEX(KIND=dbl) :: temp            ! Scratch value
COMPLEX(KIND=dbl),DIMENSION(n) :: temp1 ! Scratch array

! Make copies of arrays "a" and "b" for local use
a1 = a(1:n,1:n)
soln = b(1:n)

! Process N times to get all equations...
mainloop: DO irow = 1, n

   ! Find peak pivot for column irow in rows irow to N
   ipeak = irow
   max_pivot: DO jrow = irow+1, n
      IF (ABS(a1(jrow,irow)) > ABS(a1(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot

   ! Check for singular equations.
   singular: IF ( ABS(a1(ipeak,irow)) < epsilon ) THEN
      error = 1
      RETURN
   END IF singular

   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      temp1 = a1(ipeak,1:n)
      a1(ipeak,1:n) = a1(irow,1:n)   ! Swap rows in a
      a1(irow,1:n) = temp1
      temp = soln(ipeak)
      soln(ipeak) = soln(irow)       ! Swap rows in b
      soln(irow)  = temp
   END IF swap_eqn

   ! Multiply equation irow by -a1(jrow,irow)/a1(irow,irow),
   ! and add it to Eqn jrow (for all eqns except irow itself).
   eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a1(jrow,irow)/a1(irow,irow)
         a1(jrow,1:n) = a1(irow,1:n)*factor + a1(jrow,1:n)
         soln(jrow) = soln(irow)*factor + soln(jrow)
      END IF
   END DO eliminate
END DO mainloop

! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   soln(irow) = soln(irow) / a1(irow,irow)
END DO divide

! Set error flag to 0 and return.
error = 0

END SUBROUTINE simul_dbl_cmplx

END MODULE booklib
