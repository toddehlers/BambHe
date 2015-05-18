CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Finite Element Method for the Helium Production - Diffusion Equation
C______________________________________________________________________
C
C         imax     max number of iterations of linear solver
C         tol     tolerance of linear solver
C
C         ns       size of global stiffness matrix 
C         nz       size of global stiffness matrix after CRS 
C         nt       total number of time steps
C
C         nnode    number of nodes
C         nelem    number of tetrahedarl elements
C         nsfem    number of surface triangle elements
C         nsfnd    number of sufrace nodes
C
C         gnode    global label of trahedral
C
C         tmin     initial time
C         tmax     final time
C         dt       time step length
C
C         xnode    x coordinate of node
C         ynode    y coordinate of node
C         znode    z coordinate of node
C         vloum    volumes of trehedral elements
C
C         qt       diffusivity
C         ft       radiogenic source
C         Tt       temperature history
C
C         D0       diffusivity at infinite
C         Ea       activation energy
C         Gc       gas constant
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE MAINDRIVER_ROUNDDISC (tot_conc,D00,Ea0,nt0,
     *           therm_hist0,ages,conc_hist,nnode0,rad0,
     *           time_out,alpha_ejec0,diffuse,min_diff)
!
! main program 
! 
      USE direct_procedures_rounddisc
      IMPLICIT NONE
!
! A: global stiffness matrix 
! W: global stiffness matrix after compress row storage
! ut: Helium concentration 
! 
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: A,W
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:) :: ut
      REAL(KIND=r8) :: tot_conc(nt0),ages(nt0)

C local variables

      REAL(KIND=r8) :: perim,conc_hist(nnode0)
      REAL(KIND=r8) :: elem_conc,Ea0,D00,time_out(0:nt0)
      REAL(KIND=r8) :: const_rate,var_rate
      REAL(KIND=r8) :: area1,area2,rad0,min_diff
      INTEGER :: err,i,j,k,nt0,concFlag,nnode0,ft_count
      INTEGER :: alpha_ejec0,diffuse,jt,diffuse2
      CHARACTER i_char*4,therm_hist0*100

      write(6,*) 'Processing 2-D Round Disc'

      rad=rad0                  ! Radius of disc
      Ea=Ea0                    ! Activation energy
      D0=D00                    ! Diffusivity at infinity
      nt=nt0                    ! Number of timesteps
      therm_hist=therm_hist0    ! Time-temperature history filename
      alpha_ejec=alpha_ejec0

C parameters initilization

      CALL SETUP_ROUNDDISC
      CALL GKD_ROUNDDISC

      Fmax=1.0E-08
      Emax=1.0E-08

C allocate arrays

      ALLOCATE(A(ns),STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating W!'
         STOP
      END IF

      ALLOCATE(W(nz),STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating W!'
         STOP
      END IF

      ALLOCATE(ut(0:nt,nnode),STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating ut!'
         STOP
      END IF

C forward solver

      ! Set initial Helium concentration to user specified value
      ! or concentrations from previous output file
      ut(0,:)=conc_hist

      jt=0
      DO it=1,nt
         diffuse2=1
         IF(diffuse.eq.0 .or. qt(it).lt.min_diff)  THEN
            diffuse2=0
         END IF
         if (diffuse2 .eq. 1) then
           CALL QEB_ROUNDDISC(A,ut(it,:),ut(it-1,:))
           CALL BDE_ROUNDDISC(A,ut(it,:))
           CALL CRS_ROUNDDISC(A,W)
           CALL QMR_ROUNDDISC(W,ut(it,:))
           jt=it
         else
           do j=1,nnode
           ut(it,1:nnode)=ut(jt,1:nnode)+ft(1:nnode)*(it-jt)*dt
           enddo
         endif
         WRITE(*,*) 'time step:', it 
      END DO

C write helium concentrations and production values for each node at every time step
C added 05/08 by cspath

      ! Scale disc back to original dimensions
      !xnode=xnode*rad
      !ynode=ynode*rad
      area=area*rad**2.
      time=time*L
      ft=ft/L

!       tot_area = pi*rad**2.
      tot_area=0.
      do i=1,nelem
        tot_area=tot_area+area(i)
      enddo

      ft_count = 0                              ! Amount of points less than 20 microns from edge of disc
      const_rate = Prate                        ! Constant production rate (greater than 20 microns from edge)
      var_rate = 0.0                            ! Variable production rate (average of all rates between 20 microns and edge)
      area1 = pi*(rad-20.e-6_r8)**2.            ! Area of inner disc that contains points greater than 20 microns from edge
      area2 = tot_area-area1                    ! Area of outer shell (space between the inner disc and total disc)

      ! Iterate through all nodes to find those that are less than 20 microns from the
      ! edge of the disc
      do i=1,nnode
        if (ft(i).lt.const_rate) then
          var_rate = var_rate+ft(i)
          ft_count = ft_count+1
        endif
      enddo

      ! Divide by the number of nodes in the space to get an average
      ! production rate of the nodes within 20 microns from edge
      var_rate = var_rate/ft_count

      ! Total production rate is the ratio of the inner disc at the constant rate
      ! added with the ratio of the outer surface at the average variable rate
      if (alpha_ejec .eq. 1) then
        tot_prod_rate = (area1/tot_area)*const_rate+
     *                  (area2/tot_area)*var_rate
      else
        tot_prod_rate = const_rate
      endif

      ! Output of Helium concentration and production rate for all nodes
      ! in the disc for each timestep
      DO i=1,nt

        ! Calculate total Helium domain area for timestep by adding
        ! all node concentrations, multiplying by element area (3 nodes
        ! per element), dividing by the number of nodes per element (3), and
        ! adding over all elements
        tot_conc(i) = 0.
        DO j=1,nelem
          elem_conc = 0.
          DO k=1,3
            elem_conc = elem_conc + ut(i,gnode(k,j))
          END DO
          elem_conc = elem_conc/3.*area(j)
          tot_conc(i) = tot_conc(i) + elem_conc
        END DO
		
        ! Calculate the disc age by taking the total Helium area
        ! and dividing by the total disc area to get a concentration
        ! and then dividing by the total production rate in the disc
        ages(i)=(tot_conc(i)/tot_area)/tot_prod_rate

        write(i_char,'(i4)') i
        if (i.lt.10) then
          i_char(1:3)='000'
        else if (i.lt.100) then
          i_char(1:2)='00'
        else if (i.lt.1000) then
          i_char(1:1)='0'
        endif
        OPEN(unit=80,file='rounddisc/He_tec'//i_char//'.dat',
     *       status='unknown')

        write (80,'(A)',ADVANCE="no") 'TITLE = "bambHe Output -'
        write (80,'(A)',ADVANCE="no") ' Round disc - Time: '
        write (80,*) time(i),'s"'
        write (80,'(A)',ADVANCE="no") ' VARIABLES = "X (m)" "Y (m)"'
        write (80,'(A)',ADVANCE="no") ' "He Concentration (nmole/g)"'
        write (80,'(A)') ' "Production Rate (nmole/(g*s))"'
        write (80,'(A)',ADVANCE="no") 'ZONE T = "Round Disc - Time: '
        write (80,*) time(i),'s"'
        write (80,*) 'n=',nnode,', e=',nelem,',et=triangle,f=fepoint'

        ! Write out x, y positions of node, He concentration, and production rate
        DO j=1,nnode
          WRITE(80,'(4es15.7)') xnode(j),ynode(j),ut(i,j),
     *                         ft(j)
        END DO

        ! Write out node-element connectivities
        DO j=1,nelem
          WRITE(80,*) (gnode(k,j),k=1,3)
        END DO

        CLOSE(80)
      END DO

      time_out=time

C clean up variables

      DEALLOCATE(A,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating W!'
         STOP
      END IF

      DEALLOCATE(W,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating W!'
         STOP
      END IF

      DEALLOCATE(ut,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating ut!'
         STOP
      END IF

      CALL CLEANUP_ROUNDDISC

      RETURN
      END SUBROUTINE MAINDRIVER_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SETUP_ROUNDDISC
      USE direct_procedures_rounddisc
      IMPLICIT NONE
      
C local variables

      INTEGER :: err,i,counter
      REAL(KIND=r8) :: q_ROUNDDISC,f_ROUNDDISC,tp

C linear solver parameters

      imax=5000
      tol=1.0E-06_r8

! C time step parameters

!       tmin=0.0_r8
!       tmax=1.0_r8
!       nt=40
!       dt=(tmax-tmin)/nt

! C physical constants

!       D0=1.0_r8
!       Ea=1.0_r8
      Gc=8.314_r8           ! Universal gas constant [J/(K*mol)]

C linear matrices 1D array length

      ns=100E+04
      nz=100E+04

      tot_area=0.0
      pi=3.14159265

C read data

      OPEN(unit=81,file='rounddisc/nnee_rounddisc.dat')
      READ(81,*) nnode,nelem,nedge
      CLOSE(unit=81)

C allocate global matrix

      ALLOCATE(xnode(nnode),ynode(nnode),STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error allocating node corrdinates!'
	  STOP
      END IF

      ALLOCATE(gnode(3,nelem),STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error allocating global node label!'
	  STOP
      END IF

      ALLOCATE(area(nelem),STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error allocating element area!'
	  STOP
      END IF

      ALLOCATE(gedge(2,nedge),STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error allocating edge number!'
	  STOP
      END IF

      ALLOCATE(gboun(nedge),STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error allocating boundary nodal number!'
	  STOP
      END IF

      ALLOCATE(qt(0:nt),STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error allocating qt!'
	  STOP
      END IF

      ALLOCATE(ft(nnode),STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error allocating ft!'
	  STOP
      END IF

      ! Added 05/08 by cspath
      ! Hard coded to be max of 1000 timesteps
      ! This number can be changed if needed
      ALLOCATE(time_init(0:1000),STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error allocating time_init!'
	  STOP
      END IF

      ! Added 05/08 by cspath
      ! Hard coded to be max of 1000 timesteps
      ! This number can be changed if needed
      ALLOCATE(temp_init(0:1000),STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error allocating temp_init!'
	  STOP
      END IF

      ! Added 05/08 by cspath
      ALLOCATE(time(0:nt),STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error allocating time!'
	  STOP
      END IF

      ! Added 05/08 by cspath
      ALLOCATE(temp(0:nt),STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error allocating temp!'
	  STOP
      END IF

C read data from Matlab
	 
      OPEN(unit=82,file='rounddisc/xy_rounddisc.dat')
      DO i=1,nnode
	 READ(82,*) xnode(i),ynode(i)
      END DO
      CLOSE(unit=82)

      OPEN(unit=83,file='rounddisc/node_rounddisc.dat')
      DO i=1,nelem
         READ(83,*) gnode(1,i),gnode(2,i),gnode(3,i)
      END DO
      CLOSE(unit=83)

      OPEN(unit=84,file='rounddisc/edge_rounddisc.dat')
      DO i=1,nedge
	 READ(84,*) gedge(1,i),gedge(2,i)
      END DO
      CLOSE(unit=84)

      OPEN(unit=85,file='rounddisc/area_rounddisc.dat')
      DO i=1,nelem
	 READ(85,*) area(i)
         tot_area = tot_area + area(i)
      END DO
      CLOSE(unit=85)

C allocate global matrix

      ALLOCATE(diag(nnode+1),STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating global label of nodes!'
         STOP
      END IF

      ALLOCATE(colind(nz),STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error allocating colind!'
	  STOP
      END IF

      ALLOCATE(rowptr(2,nnode+1),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating rowptr!'
	  STOP
      END IF

      ALLOCATE(JE(10),STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating JE! '
         STOP
      END IF

      ALLOCATE(JP(20),STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating JP! '
         STOP
      END IF

C read time and temperature histories from data file
C added 05/08 by cspath

      OPEN(unit=86,file=therm_hist)
      READ(86,*)
      counter=0
      DO
        READ(86,*,end=20) time_init(counter),temp_init(counter)
        counter=counter+1
      END DO
20    CLOSE(unit=86)

      counter=counter-1
      L=time_init(counter)              ! Set to last (greatest) time value
      time_init=time_init/L             ! Normalize time history
	  
C calculate tmax, tmin, and dt from read in data
      tmin = time_init(0)
      tmax = time_init(counter)
      dt = (tmax-tmin)/nt

      ! Create a linearly spaced time history based on user
      ! specified number of timesteps
      DO i=0,nt
        time(i) = tmin + i*dt
      END DO

      ! Interpolates the temperature history onto the linearly
      ! spaced time history
      call interp_ROUNDDISC(temp_init,time_init,counter,
     *       temp,time,nt)

C exact qt

      DO i=0,nt
        qt(i)=q_ROUNDDISC(temp(i))
      END DO

!       DO i=0,nt
! 	 tp=tmin+i*dt
! 	 qt(i)=q(tp)
!       END DO

C exact ft

      DO i=1,nnode
         ft(i)=f_ROUNDDISC(xnode(i),ynode(i))
      END DO

C compute the boundary nodal point label

      CALL BOUN_ROUNDDISC

      END SUBROUTINE SETUP_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE BOUN_ROUNDDISC
      USE direct_procedures_rounddisc
      IMPLICIT NONE

C local variable

      INTEGER :: i,i1,i2,j,jb,j1,j2,ip

      gboun=0
      ip=0

      DO i=1,nedge
	 i1=gedge(1,i)
	 i2=gedge(2,i)
	 j1=0
	 j2=0
	 DO j=1,ip
            jb=gboun(j)
	    IF (i1 .EQ. jb) THEN
	        j1=1
	    END IF
	    IF (i2 .EQ. jb) THEN
	        j2=1
	    END IF
	 END DO
         IF (j1 .EQ. 0) THEN
	     ip=ip+1
	     gboun(ip)=i1
	 END IF
	 IF (j2 .EQ. 0) THEN
	     ip=ip+1
	     gboun(ip)=i2
	 END IF
      END DO

      nboun=ip

      END SUBROUTINE BOUN_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GKD_ROUNDDISC
      USE direct_procedures_rounddisc
      IMPLICIT NONE

C     local variables

      INTEGER :: err,L0,LI

      diag(1)=0
      DO L0=1,nnode
         CALL GJP_ROUNDDISC(L0)
         LI=L0-JP(1)
         diag(L0+1)=diag(L0)+LI+1
      END DO
      ns=diag(nnode+1)

      END SUBROUTINE GKD_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GJP_ROUNDDISC(L0)
      USE direct_procedures_rounddisc
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: L0

C local variables    

      INTEGER :: LL,LP,k,ip,LE,IS,kk,i,i0,j0

      LL=0
      LP=0
      CALL GJE_ROUNDDISC(L0,LL)
      DO k=1,3
         JP(k)=gnode(k,JE(1))
      END DO
      ip=3
      DO LE=2,LL
         DO k=1,3
            IS=gnode(k,JE(LE))
            DO kk=1,ip
               IF (IS .NE. JP(kk)) THEN
                  IF (kk .EQ. ip) THEN
                     ip=ip+1
                     JP(ip)=IS
                     GO TO 40
                  END IF
               ELSE
                  GO TO 40
               END IF
            END DO
40       END DO
      END DO
      LP=ip
      DO i=1,LP-1
         i0=JP(i)
         DO k=i+1,LP
            IF (JP(k) .LT. i0) THEN
               j0=k
               i0=JP(k)
            ENDIF
         END DO
         IF (i0 .NE. JP(i)) THEN
            JP(j0)=JP(i)
            JP(i)=i0
         END IF
      END DO
      
      END SUBROUTINE GJP_ROUNDDISC
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GJE_ROUNDDISC(L0,LL)
      USE direct_procedures_rounddisc
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: L0
      INTEGER,INTENT(INOUT) :: LL

C local variables
 
      INTEGER :: i,LE,k

      i=0
      DO LE=1,nelem
         DO k=1,3
            IF (gnode(k,LE) .EQ. L0) THEN
                i=i+1
                JE(i)=LE
                GO TO 30
            END IF
         END DO
30    END DO
      LL=i
      
      END SUBROUTINE GJE_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE QEB_ROUNDDISC(A,un,uo)
      USE direct_procedures_rounddisc
      IMPLICIT NONE

      REAL(KIND=r8),DIMENSION(ns),INTENT(OUT) :: A
      REAL(KIND=r8),DIMENSION(nnode),INTENT(IN) :: uo
      REAL(KIND=r8),DIMENSION(nnode),INTENT(OUT) :: un

C local variables

      INTEGER :: I,J,LE,NI,NJ
      REAL(KIND=r8) :: FE,Aij
      REAL(KIND=r8),DIMENSION(3,3) :: E,F

C initialization

      DO I=1,ns
         A(I)=0.0_r8
      END DO
      DO I=1,nnode
         un(I)=0.0_r8
      END DO

      Aij=0.0_r8
      FE=0.0_r8

      DO LE=1,nelem
         CALL UEF_ROUNDDISC(LE,E,F)
         DO NI=1,3
            CALL GF_ROUNDDISC(NI,LE,FE,uo)
            I=gnode(NI,LE)
            un(I)=un(I)+FE
            DO NJ=1,NI
               J=gnode(NJ,LE)
               CALL GK_ROUNDDISC(NI,NJ,LE,Aij,E,F)
               CALL AK_ROUNDDISC(I,J,Aij,A)
            END DO
         END DO
      END DO
      
      END SUBROUTINE QEB_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE UEF_ROUNDDISC(LE,E,F)
      USE direct_procedures_rounddisc
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: LE
      REAL(KIND=r8),DIMENSION(3,3),INTENT(OUT) :: E,F

C local variables

      INTEGER :: i,j
      REAL(KIND=r8) :: ds
      REAL(KIND=r8),DIMENSION(3) :: x,y,b,c

      ds=area(LE)

      DO i=1,3
         x(i)=xnode(gnode(i,LE))
         y(i)=ynode(gnode(i,LE))
      END DO

      b(1)=y(2)-y(3)
      b(2)=y(3)-y(1)
      b(3)=y(1)-y(2)

      c(1)=x(3)-x(2)
      c(2)=x(1)-x(3)
      c(3)=x(2)-x(1)

      DO i=1,3
         DO j=1,3
            E(i,j)=0.25_r8*(b(i)*b(j)+c(i)*c(j))/ds
            IF (i .EQ. j) THEN
               F(i,j)=ds/6.0_r8
            ELSE
               F(i,j)=ds/12.0_r8
            END IF
         END DO
      END DO   

      END SUBROUTINE UEF_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GF_ROUNDDISC(NI,LE,FE,uo)
      USE direct_procedures_rounddisc
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: NI,LE
      REAL(KIND=r8),INTENT(OUT) :: FE
      REAL(KIND=r8),DIMENSION(nnode),INTENT(IN) :: uo

C local variables

      INTEGER :: i1,i2,i3
      REAL(KIND=r8) :: ds,au,af

      ds=area(LE)

      i1=gnode(1,LE)
      i2=gnode(2,LE)
      i3=gnode(3,LE)

      au=(uo(i1)+uo(i2)+uo(i3))/3.0_r8
      af=(ft(i1)+ft(i2)+ft(i3))/3.0_r8
      FE=(dt*af+au)*ds/3.0_r8


      END SUBROUTINE GF_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GK_ROUNDDISC(NI,NJ,LE,Aij,E,F)
      USE direct_procedures_rounddisc
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: NI,NJ,LE
      REAL(KIND=r8),INTENT(OUT) :: Aij
      REAL(KIND=r8),DIMENSION(3,3),INTENT(IN) :: E,F

C local variables

      IF (ABS(F(NI,NJ)).GT.Fmax) THEN
         Fmax=ABS(F(NI,NJ))
      END IF
      IF (ABS(F(NI,NJ)).LT.Fmin) THEN
         Fmin=ABS(F(NI,NJ))
      END IF
      IF (ABS(E(NI,NJ)).GT.Emax) THEN
         Emax=ABS(E(NI,NJ))
      END IF
      IF (ABS(E(NI,NJ)).LT.Emin) THEN
         Emin=ABS(E(NI,NJ))
      END IF

      Aij=F(NI,NJ)+qt(it)*dt*E(NI,NJ)

      END SUBROUTINE GK_ROUNDDISC
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE AK_ROUNDDISC(i,j,Aij,A)
      USE direct_procedures_rounddisc
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i,j
      REAL(KIND=r8),INTENT(IN) :: Aij
      REAL(KIND=r8),DIMENSION(ns),INTENT(INOUT) :: A

C local variables

      INTEGER :: IGP

      IF (i. GE. j) THEN
         IGP=diag(i+1)-i+j
      ELSE
         IGP=diag(j+1)-j+i
      END IF
      A(IGP)=A(IGP)+Aij
      
      END SUBROUTINE AK_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE CRS_ROUNDDISC(A,W)
      USE direct_procedures_rounddisc
      IMPLICIT NONE
      REAL(KIND=r8),DIMENSION(ns),INTENT(IN) :: A
      REAL(KIND=r8),DIMENSION(nz),INTENT(OUT) :: W

C local variables
 
      INTEGER :: ip,i,Li,mi,j,IGP,LJ,mj

      ip=0
      rowptr(1,1)=1

      DO i=1,nnode
         Li=diag(i+1)-diag(i)
         mi=i-Li+1
         DO j=mi,i
            IGP=diag(i+1)-(i-j)
            IF (ABS(A(IGP)) .GT. tol) THEN
               ip=ip+1
               W(ip)=A(IGP)
               colind(ip)=j
            END IF        
         END DO
         rowptr(2,i)=ip
         DO j=i+1,nnode
            LJ=diag(j+1)-diag(j)
            mj=j-LJ+1
            IF (mj .LE. i) THEN
               IGP=diag(j+1)-(j-i)
               IF (ABS(A(IGP)) .GT. tol) THEN
                  ip=ip+1
                  W(ip)=A(IGP)
                  colind(ip)=j
               END IF
            END IF
         END DO
         rowptr(1,i+1)=ip+1
      END DO

      END SUBROUTINE CRS_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE BDE_ROUNDDISC(A,un)
      USE direct_procedures_rounddisc
      IMPLICIT NONE
      REAL(KIND=r8),DIMENSION(ns),INTENT(INOUT) :: A
      REAL(KIND=r8),DIMENSION(nnode),INTENT(INOUT) :: un

C local variables
     
      INTEGER :: LE,i,IIG,mi,j,IGP,JJG,mj

      DO LE=1,nboun
	 i=gboun(LE)
         un(i)=0.0_r8
         IIG=diag(i+1)-i
         mi=diag(i)-IIG+1
         DO j=mi,i-1
            IGP=IIG+j
            A(IGP)=0.0_r8
         END DO
         A(diag(i+1))=1.0_r8
         DO j=i+1,nnode
            JJG=diag(j+1)-j
            IGP=i+JJG
            mj=diag(j)-JJG+1
            IF (i .GE. mj) THEN
                A(IGP)=0.0_r8
            END IF
         END DO
      END DO

      END SUBROUTINE BDE_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE QMR_ROUNDDISC(W,un)
      USE direct_procedures_rounddisc
      IMPLICIT NONE
      REAL(KIND=r8),DIMENSION(nz),INTENT(IN) :: W
      REAL(KIND=r8),DIMENSION(nnode),INTENT(INOUT) :: un

C local variables

      INTEGER :: i,iter
      REAL(KIND=r8) :: Vm,rho1,rho2,c0,c1,hu0,hu1,Rss
      REAL(KIND=r8) :: Vnorm_ROUNDDISC,Cnorm_ROUNDDISC
      REAL(KIND=r8) :: eps0,eps1,eta,del,eps2
      REAL(KIND=r8),DIMENSION(nnode) :: Av,P,U,R,T     

      iter=0       
      Vm=Vnorm_ROUNDDISC(un)

      DO i=1,nnode
         U(i)=un(i)
         P(i)=0.0_r8
         R(i)=0.0_r8
         un(i)=0.0_r8
      END DO

      rho1=SQRT(Vnorm_ROUNDDISC(U))
      c0=1.0_r8
      eps0=1.0_r8
      hu0=0.0_r8
      eta=-1.0_r8

10    DO i=1,nnode
         U(i)=U(i)/rho1
         T(i)=U(i)/W(rowptr(2,i))
      END DO

      del=Cnorm_ROUNDDISC(U,T)

      DO i=1,nnode
         P(i)=T(i)-(rho1*del/eps0)*P(i)
      END DO
      call Mult_ROUNDDISC(W,P,Av)
      eps1=Cnorm_ROUNDDISC(P,Av)
      eps2=eps1/del
      DO i=1,nnode
         U(i)=Av(i)-eps2*U(i)
      END DO

      rho2=SQRT(Vnorm_ROUNDDISC(U))
      hu1=rho2/(c0*abs(eps2))
      c1=1.0_r8/SQRT(1.0_r8+hu1**2)
      eta=-eta*rho1*c1**2/(eps2*c0**2)

      DO i=1,nnode
         R(i)=eta*P(i)+((hu0*c1)**2)*R(i)
         un(i)=un(i)+R(i)
      END DO

      c0=c1
      hu0=hu1
      rho1=rho2
      eps0=eps1

      iter=iter+1
      Rss=SQRT(Vnorm_ROUNDDISC(R)/Vm)

c      WRITE(6,*) iter,Rss

      IF (Rss .LE. tol) then
c         WRITE(6,*) 'Convergence Achieved! '
      ELSE
         IF (iter .EQ. imax) THEN
            WRITE(6,*) 'imax Exceeded No Convergence!'
         ELSE
            GO TO 10
         END IF
      END IF

      END SUBROUTINE QMR_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE Mult_ROUNDDISC(W,un,Av)
      USE direct_procedures_rounddisc
      IMPLICIT NONE
      REAL(KIND=r8),DIMENSION(nz),INTENT(IN) :: W
      REAL(KIND=r8),DIMENSION(nnode),INTENT(IN):: un
      REAL(KIND=r8),DIMENSION(nnode),INTENT(OUT) :: Av

C local variables

      INTEGER :: i,j

      DO i=1,nnode
         Av(i)=0.0_r8
         DO j=rowptr(1,i),rowptr(1,i+1)-1
            Av(i)=Av(i)+W(j)*un(colind(j))
         END DO
      END DO
      
      END SUBROUTINE Mult_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE CLEANUP_ROUNDDISC
      USE direct_procedures_rounddisc
      IMPLICIT NONE

C local variables
 
      INTEGER :: err


      DEALLOCATE(xnode,ynode,STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error deallocating node corrdinates!'
	  STOP
      END IF

      DEALLOCATE(gnode,STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error deallocating global node label!'
	  STOP
      END IF

      DEALLOCATE(area,STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error deallocating element area!'
	  STOP
      END IF

      DEALLOCATE(gedge,STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error deallocating surface edge number!'
	  STOP
      END IF

      DEALLOCATE(gboun,STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error deallocating boundary nodal number!'
	  STOP
      END IF

      DEALLOCATE(qt,STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error deallocating qt, qn, and dM!'
	  STOP
      END IF

      DEALLOCATE(ft,STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error deallocating ft, fn, and dN!'
	  STOP
      END IF

      DEALLOCATE(diag,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating global label of nodes!'
         STOP
      END IF

      DEALLOCATE(colind,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating colind!'
         STOP
      END IF

      DEALLOCATE(rowptr,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating rowptr!'
         STOP
      END IF

      DEALLOCATE(JE,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating JE! '
         STOP
      END IF

      DEALLOCATE(JP,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating JP! '
         STOP
      END IF

      DEALLOCATE(time,time_init,STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error deallocating time array!'
	  STOP
      END IF

      DEALLOCATE(temp,temp_init,STAT=err)
      IF (err .NE. 0) THEN
	  WRITE(6,*) 'Error deallocating temperature array!'
	  STOP
      END IF

      END SUBROUTINE CLEANUP_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION Vnorm_ROUNDDISC(R)
      USE direct_procedures_rounddisc
      REAL(KIND=r8),DIMENSION(nnode),INTENT(IN) :: R
      REAL(KIND=r8) :: Vnorm_ROUNDDISC

C local variables

      INTEGER :: i
      REAL(KIND=r8) :: sum

      sum=0.0_r8
      DO i=1,nnode
         sum=sum+ABS(R(i))**2
      END DO
      Vnorm_ROUNDDISC=sum
      
      END FUNCTION Vnorm_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION Cnorm_ROUNDDISC(X,Y)
      USE direct_procedures_rounddisc
      REAL(KIND=r8) :: Cnorm_ROUNDDISC
      REAL(KIND=r8),DIMENSION(nnode),INTENT(IN) :: X,Y

C local variables

      INTEGER :: i
      REAL(KIND=r8) :: sum

      sum=0.0_r8
      DO i=1,nnode
         sum=sum+X(i)*Y(i)
      END DO
      Cnorm_ROUNDDISC=sum
      
      END FUNCTION Cnorm_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION q_ROUNDDISC(t)
!
! diffusivity function
! 
      USE direct_procedures_rounddisc
      IMPLICIT NONE
      REAL(KIND=r8) :: t

C local variables
      REAL(KIND=r8) :: q_ROUNDDISC,Te_ROUNDDISC
!       REAL(KIND=r8) :: t0,q

c example 1

c      q=COS(t)

c example 2

!       t0=0.7_r8
! 
!       IF (t .LE. t0) THEN 
!           q=1.0
!       ELSE
!           q=0.5_r8
!       END IF


      q_ROUNDDISC=D0*EXP(-Ea/(Gc*(t+273.15)))
      q_ROUNDDISC=q_ROUNDDISC*(L/rad**2.)

      END FUNCTION q_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION f_ROUNDDISC(x,y)
!
! radiogenic source
!
      USE direct_procedures_rounddisc
      IMPLICIT NONE
      REAL(KIND=r8),INTENT(IN) :: x,y

C local variable
  
      REAL(KIND=r8) :: r,r0,f_ROUNDDISC
      REAL(KIND=r8) :: Um,Thm,Thmc,Umc

      Um=50.                                    ! Uranium concentration [ppm]
      Thm=50.                                   ! Thorium concentration [ppm]
      Umc=Um*(1/238.03)*1e3                     ! Convert U concentration from [ppm]-->[nmol/g]
      Thmc=Thm*(1/232.038)*1e3                  ! Convert Th concentration from [ppm]-->[nmol/g]
      Prate=8.*4.95e-18*Umc+6.*1.56e-18*Thmc    ! He production rate [nmol/(gm*s)]

      r0=20.e-6_r8                                      ! Distance [m] from edge to start decreasing production rate
      r0=r0/rad                                         ! Find same ratio of 20 mircons of actual mesh to unit mesh
      r=SQRT(x*x+y*y)                                   ! Current distance of point from center of disc

      if (alpha_ejec .eq. 1) then
        IF (r .LE. (1.0_r8-r0)) THEN                      ! If the current distance is closer to the center than 20 microns from edge
           f_ROUNDDISC=Prate
           f_ROUNDDISC=f_ROUNDDISC*L                      ! Scale production rate
        ELSE
           f_ROUNDDISC=0.5_r8*Prate*(1.0_r8-r)/r0+        ! Linearly decrease production rate to half as distance approaches edge of disc
     *                 0.5_r8*Prate
           f_ROUNDDISC=f_ROUNDDISC*L                      ! Scale production rate
        END IF
      else
        f_ROUNDDISC=Prate
        f_ROUNDDISC=f_ROUNDDISC*L
      endif

      END FUNCTION f_ROUNDDISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Subroutine that linearly interpolates a temperature history
      ! onto a linearly spaced time history
      SUBROUTINE interp_ROUNDDISC(tp0,ti0,num0,tpf,tif,numf)

      USE direct_procedures_rounddisc
      IMPLICIT NONE

      INTEGER :: num0,numf,i,j
      REAL(KIND=r8) :: tp0(0:num0),ti0(0:num0)
      REAL(KIND=r8) :: tpf(0:numf),tif(0:numf)

      tpf(0)=tp0(0)
      tpf(numf)=tp0(num0)

      DO i=0,numf-1
        DO j=0,num0-1
          IF (ti0(j).eq.tif(i)) THEN
            tpf(i) = tp0(j)
          ELSE IF (ti0(j).lt.tif(i) .and. ti0(j+1).gt.tif(i)) THEN
            tpf(i) = tp0(j) + (tp0(j+1)-tp0(j))/(ti0(j+1)-ti0(j))*
     *               (tif(i)-ti0(j))
          END IF
        END DO
      END DO

      END SUBROUTINE interp_ROUNDDISC





