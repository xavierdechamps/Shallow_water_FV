!##########################################################
! SUBROUTINE runge_kutta
!  Goal: temporal integration: 4 steps explicit Runge-Kutta
!##########################################################
SUBROUTINE runge_kutta
    USE module_shallow
#ifdef WITHSIGWATCH

!DEC$ IF .NOT. DEFINED( SIGNAL_NAMES )
!DEC$ DEFINE SIGNAL_NAMES

#ifdef UPPERCASE
#define Watch_signal_name WATCHSIGNALNAME
#define Watch_signal WATCHSIGNAL
#define Get_last_signal GETLASTSIGNAL
#endif

#ifdef UPPERCASE_
#define Watch_signal_name WATCHSIGNALNAME_
#define Watch_signal WATCHSIGNAL_
#define Get_last_signal GETLASTSIGNAL_
#endif

#ifdef lowercase_
#define Watch_signal_name watchsignalname_
#define Watch_signal watchsignal_
#define Get_last_signal getlastsignal_
#endif

#ifdef lowercase
#define Watch_signal_name watchsignalname
#define Watch_signal watchsignal
#define Get_last_signal getlastsignal
#endif

!DEC$ ENDIF

    USE signal_handler
#endif
! On Windows, call systemqq instead of system (Linux)
#ifdef WINDOWS
    USE IFPORT, ONLY : systemqq
#endif

    IMPLICIT NONE

    ! Local parameters
    INTEGER(ki) :: i, k, count, ok
    INTEGER(ki) :: statussignal
    REAL(kr), DIMENSION(1:nbvar*nbrElem) :: Urk, Uc, H, gradX, gradY
!    REAL(kr), DIMENSION(4) :: beta
!    REAL(kr), DIMENSION(3) :: alpha
    REAL(kr), ALLOCATABLE  :: error(:,:)
    CHARACTER(20) :: char
    LOGICAL :: kill
    
    LOGICAL(4)           lreturn
    INTEGER ( KIND = 4 ) ierror
    CHARACTER ( LEN = 255 ) :: command 
    CHARACTER ( LEN = 255 ) :: data_filename = 'gnuplot_data.txt'
    CHARACTER ( LEN = 255 ) :: command_filename = 'gnuplot_commands.txt'
    
    ! DAM
    ! REAL(kr), DIMENSION(4) :: dam_points_ref
    ! INTEGER(ki) :: ierr
    ! OPEN(UNIT=105,FILE="dam_points_ref",STATUS="replace",IOSTAT=ierr,FORM='formatted')

#ifdef WITHSIGWATCH
    statussignal = Watch_signal(2)  ! INT
    statussignal = Watch_signal(15) ! KILL
#endif
    kill=.false.
    
    CALL sampletime(time_begin)

!    alpha = (/0.5,0.5,1./)
!    beta = (/6.,3.,3.,6./)
    count = 0
    
    Source = zero
    
    WRITE(*,*) "Begin of time integration"
    
    IF (muscl.NE.0 .OR. amr.NE.0) THEN
      CALL getGradients(U0,gradX,gradY)
      CALL applyTVD_Gradients(gradX,gradY)
    ENDIF
    
    IF (amr.NE.0) THEN
      CALL amr_get_relative_indicator(gradX,gradY)
      CALL amr_split_cells()
    ENDIF
    
    ! Write the initial solution
    CALL write_gmsh(U0,nbvar*nbrElem,file_gmsh,length_names,node,elem,front,nbrNodes,nbrElem,nbrTris,nbrQuads,&
&                   nbrFront,nbr_nodes_per_elem,0,0)
    CALL write_solution(U0,nbvar*nbrElem,file_dat,length_names,ok)
    IF (ok == 0) GOTO 200

    IF (.not.ALLOCATED(error)) THEN
        ALLOCATE( error( nTime ,4 ) ); error = zero
    END IF

    DO i=1,nTime ! Start the temporal loop
       IF (steady.NE.0) THEN
          ! Get the local time step for each cell
          CALL deltaT
       ELSE
          dt = deltaTfixed
       ENDIF
       
!       ! The Runge-Kutta method in itself
!       Urk = U0
!       ! Get the right-hand side term at the first sub-time step
!       CALL flux(H,Source,Urk)
!       ! Uc = U0 - dt*H/beta(1)
!       ! Update the solution at the first sub-time step
!       Uc = U0 + dt*(Source-H)/beta(1)
!       DO k=1,3
!          Urk = U0 + dt*(Source-H)*alpha(k)
!          CALL flux(H,Source,Urk)
!          Uc = Uc + dt*(Source-H)/beta(k+1)
!       END DO

      ! The Runge-Kutta method in itself
       CALL flux(H,Source,U0,gradX,gradY)
       Urk = U0 + dt*(Source-H)
       CALL flux(H,Source,Urk,gradX,gradY)
       Urk = 0.75d00*U0 + 0.25d00*(Urk + dt*(Source-H))
       CALL flux(H,Source,Urk,gradX,gradY)
       Uc  = U0/3.0d00 + 2.0d00/3.0d00 * ( Urk + dt*(Source-H) )
       
       ! DAM
       ! Dam Point 1 [100m] - 1534 1535 1564 1565 
      ! dam_points_ref(1) = 0.25d00*(Uc(1534*nbvar-2) + Uc(1535*nbvar-2) + Uc(1564*nbvar-2) + Uc(1565*nbvar-2))
       ! Dam Point 2 [110m] - 2566 2567 2604 2605
      ! dam_points_ref(2) = 0.25d00*(Uc(2566*nbvar-2) + Uc(2567*nbvar-2) + Uc(2604*nbvar-2) + Uc(2605*nbvar-2))
       ! Dam Point 3 [130m] - 2558 2559 2596 2597
      ! dam_points_ref(3) = 0.25d00*(Uc(2558*nbvar-2) + Uc(2559*nbvar-2) + Uc(2596*nbvar-2) + Uc(2597*nbvar-2))
       ! Dam Point 4 [130m] - 2550 2551 2588 2589
      ! dam_points_ref(4) = 0.25d00*(Uc(2550*nbvar-2) + Uc(2551*nbvar-2) + Uc(2588*nbvar-2) + Uc(2589*nbvar-2))
     
      ! write(105,'(I9,5ES25.16E3)') i,i*deltaTfixed,dam_points_ref(1:4)
       
       ! Print relative errors, useful for steady flows
       error(i,1) = i
       error(i,2) = SQRT(SUM((Uc(1:nbvar*nbrElem:nbvar)-U0(1:nbvar*nbrElem:nbvar))**2))/SQRT(SUM((U0(1:nbvar*nbrElem:nbvar))**2))
       error(i,3) = SQRT(SUM((Uc(2:nbvar*nbrElem:nbvar)-U0(2:nbvar*nbrElem:nbvar))**2))/SQRT(SUM((U0(2:nbvar*nbrElem:nbvar))**2))
       error(i,4) = SQRT(SUM((Uc(3:nbvar*nbrElem:nbvar)-U0(3:nbvar*nbrElem:nbvar))**2))/SQRT(SUM((U0(3:nbvar*nbrElem:nbvar))**2))
       WRITE(char,'(F6.6)') error(i,2)
       char = TRIM(ADJUSTL(char))
       IF (char=="NaN") THEN
          WRITE(*,'(A24)') "DIVERGENCE : error = NaN"
          EXIT
       END IF
       
       ! Every X time steps, print the relative error, plot it with gnuplot
       IF (MOD(i,shownTime)==0) THEN
          WRITE(*,'(A10,3X,I9,A1,I9,A13,3ES15.8)') "Time step ",i,"/",nTime," - Errors : ", error(i,2:4)
    
#ifdef WITHGNUPLOT
          IF (steady.NE.0) THEN

#ifdef WINDOWS
            lreturn = systemqq ( 'start "" /b nircmd win close title "Gnuplot"' ) ! to close the previous gnuplot windows
#else
            ierror = system ( "wmctrl -c Gnuplot" ) ! to close the previous gnuplot windows
#endif
            CALL write_xyy_data ( data_filename, i, i, 4, error(1:i,:), ierror )
            CALL write_xyy_plots ( command_filename, data_filename, 4,ierror )
            CALL run_gnuplot ( command_filename )
          
          ENDIF
#endif
          
       END IF
       
       ! Every Y time steps, save the solution
       count = count + 1
       IF (mod(count,savenTime)==0) THEN
          WRITE(*,*) "Writing the solution in .msh and .dat format"
          CALL write_gmsh(U0,nbvar*nbrElem,file_gmsh,length_names,node,elem,front,nbrNodes,nbrElem,nbrTris,nbrQuads,&
&                         nbrFront,nbr_nodes_per_elem,i,count)
          CALL write_solution(U0,nbvar*nbrElem,file_dat,length_names,ok)
          IF (ok == 0) EXIT
       END IF

       U0 = Uc
       Uc = zero
       
       IF (muscl.NE.0) THEN
         CALL getGradients(U0,gradX,gradY)
         CALL applyTVD_Gradients(gradX,gradY)
       ENDIF
       
       IF (amr.NE.0) THEN
         CALL amr_get_relative_indicator(gradX,gradY)
       ENDIF
       
#ifdef WITHSIGWATCH
       statussignal = Get_last_signal()
       IF (statussignal.eq.2 .or. statussignal.eq.15)   THEN
         kill=.true.
       ENDIF
       
       IF (kill) THEN
         WRITE(*,*) "Caught a kill signal, saving the data and stopping"
         CALL write_gmsh(U0,nbvar*nbrElem,file_gmsh,length_names,node,elem,front,nbrNodes,nbrElem,nbrTris,nbrQuads,&
&                        nbrFront,nbr_nodes_per_elem,i,count)
         CALL write_solution(U0,nbvar*nbrElem,file_dat,length_names,ok)
         EXIT
       ENDIF
#endif

    END DO

    ! DAM
   ! CLOSE(UNIT=105)

    CALL sampletime(time_end)
    CALL time_display

200 CONTINUE

    IF (ALLOCATED(error)) DEALLOCATE( error )

END SUBROUTINE runge_kutta


SUBROUTINE getGradients(q,gradX,gradY)
    USE module_shallow
    USE OMP_LIB
    IMPLICIT NONE
    
    ! Subroutine parameters
    REAL(kr), DIMENSION(1:nbvar*nbrElem), INTENT(IN) :: q
    REAL(kr), DIMENSION(1:nbvar*nbrElem), INTENT(OUT) :: gradX, gradY

    ! Local parameters
    INTEGER(ki) :: i,j,k,IDj
    REAL(kr)    :: Ixx,Iyy,Ixy,D
    REAL(kr), DIMENSION(nbvar) :: Jx,Jy
    
    gradX = zero 
    gradY = zero
    
!$OMP PARALLEL &
!$OMP& default (shared) &
!$OMP& private (Ixx,Iyy,Ixy,Jx,Jy,D,j,k,IDj) 
!$OMP DO
    DO i=1,nbrElem
      Ixx = zero
      Iyy = zero
      Ixy = zero
      Jx  = zero
      Jy  = zero
      DO j=1,nbr_nodes_per_elem(i)
        IDj = geom_data_ind(i,j)
        IF (IDj.EQ.0) CYCLE
        Ixx = Ixx + (geom_data(IDj,3)-geom_data(i,3))*(geom_data(IDj,3)-geom_data(i,3))
        Iyy = Iyy + (geom_data(IDj,4)-geom_data(i,4))*(geom_data(IDj,4)-geom_data(i,4))
        Ixy = Ixy + (geom_data(IDj,3)-geom_data(i,3))*(geom_data(IDj,4)-geom_data(i,4))
        DO k=1,nbvar
          Jx(k) = Jx(k) + (geom_data(IDj,3)-geom_data(i,3))*(q((IDj-1)*nbvar+k)-q((i-1)*nbvar+k))
          Jy(k) = Jy(k) + (geom_data(IDj,4)-geom_data(i,4))*(q((IDj-1)*nbvar+k)-q((i-1)*nbvar+k))
        ENDDO
      ENDDO
      
      D   = Ixx*Iyy - Ixy*Ixy
      DO k=1,nbvar
        gradX((i-1)*nbvar+k) = (Jx(k)*Iyy - Jy(k)*Ixy)/D
        gradY((i-1)*nbvar+k) = (Jy(k)*Ixx - Jx(k)*Ixy)/D
      ENDDO
      
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE getGradients

SUBROUTINE applyTVD_Gradients(gradX,gradY)
    USE module_shallow
    USE OMP_LIB
    IMPLICIT NONE
    
    ! Subroutine parameters
    REAL(kr), DIMENSION(1:nbvar*nbrElem), INTENT(INOUT) :: gradX, gradY

    ! Local parameters
    INTEGER(ki) :: i,j,k,IDj,IDk
    REAL(kr), DIMENSION(1:nbvar*nbrElem) :: gradXlim, gradYlim, signX, signY
    REAL(kr), DIMENSION(1:nbvar)    :: minVx, minVy, minSVx, minSVy, maxSVx, maxSVy
    
    signX = gradX
    signY = gradY
    
    WHERE (signX>zero)
       signX = 1.d00
    ELSEWHERE (signX<zero)
       signX = -1.d00
    END WHERE
    
    WHERE (signY>zero)
       signY = 1.d00
    ELSEWHERE (signY<zero)
       signY = -1.d00
    END WHERE
    
!$OMP PARALLEL &
!$OMP& default (shared) &
!$OMP& private (minVx,minVy,minSVx,minSVy,maxSVx,maxSVy,j,IDj,k,IDk) 
!$OMP DO
    DO i=1,nbrElem
      minVx  = 1.E16
      minVy  = 1.E16
      minSVx = 1.E16
      minSVy = 1.E16
      maxSVx = zero
      maxSVy = zero
      
      DO j=1,nbr_nodes_per_elem(i)
        IDj = geom_data_ind(i,j)
        DO k=1,nbvar
          IDk = (IDj-1)*nbvar+k
          minVx(k) = min( minVx(k) , abs(gradX( IDk ) ) )
          minVy(k) = min( minVy(k) , abs(gradY( IDk ) ) )
          minSVx(k)= min( minSVx(k), signX( IDk ) )
          minSVy(k)= min( minSVy(k), signY( IDk ) )
          maxSVx(k)= max( maxSVx(k), signX( IDk ) )
          maxSVy(k)= max( maxSVy(k), signY( IDk ) )
        ENDDO
      ENDDO
      
      DO k=1,nbvar
        gradXlim((i-1)*nbvar+k) = 0.5d00 * ( minSVx(k) + maxSVx(k) ) * minVx(k)
        gradYlim((i-1)*nbvar+k) = 0.5d00 * ( minSVy(k) + maxSVy(k) ) * minVy(k)
      ENDDO
      
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    
    gradX = gradXlim
    gradY = gradYlim
    
END SUBROUTINE applyTVD_Gradients

!##########################################################
! SUBROUTINE sampletime
!  Goal: get the system time
!##########################################################
SUBROUTINE sampletime(counter)
    USE module_shallow, ONLY : kr,ki
    IMPLICIT NONE
          
    ! variables passed through header
    INTEGER(ki) ::counter
      
    ! variables declared locally
    INTEGER(ki) ::rate, contmax
 
    ! Determine CPU time
    CALL system_clock(counter, rate, contmax )
END SUBROUTINE sampletime

!##########################################################
! SUBROUTINE time_display
!  Goal: print on the screen the time needed by an operation
!##########################################################
SUBROUTINE time_display
    USE module_shallow
    IMPLICIT NONE

    REAL(kr) :: time
    INTEGER(ki) :: job, cont3
    INTEGER(ki) :: rate, contmax, itime

    CALL system_clock(cont3, rate, contmax )
    IF (time_end .ge. time_begin) THEN
       itime=time_end-time_begin
    ELSE
       itime=(contmax - time_begin) + (time_end + 1)
    ENDIF

    time = DFLOAT(itime) / DFLOAT(rate)

    WRITE(*,'(a,f10.4)') " Time needed (s) :            ",time

END SUBROUTINE time_display
