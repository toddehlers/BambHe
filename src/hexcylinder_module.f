      MODULE direct_procedures_hexcylinder
      IMPLICIT NONE

C r8 is 8-byte (double precision) real

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)

C global variables 

      INTEGER :: imax,ns,nz,nt,it,alpha_ejec
      INTEGER :: nnode,nelem,nsfem,nsfnd
      
      INTEGER,ALLOCATABLE,DIMENSION(:) :: colind,diag,JE,JP,gboun
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: gnode,gsfnd,rowptr

      REAL(KIND=r8) :: tol,tmin,tmax,dt,amp,zamp,side,pi,init_conc
      REAL(KIND=r8) :: D0,Ea,Gc,tot_volum,tot_prod_rate,L,Prate
      
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: volum
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: xnode,ynode,znode
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: qt,ft

      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: time_init,temp_init
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: time,temp

      REAL(KIND=r8) :: Emax,Emin,Fmax,Fmin

      CHARACTER :: therm_hist*100

      END MODULE direct_procedures_hexcylinder