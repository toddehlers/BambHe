
      program bambHe

      implicit none

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)

      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: tot_conc_ball
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: tot_conc_hexcyl,tot_conc_hexcyl2,tot_conc_hexcyl3
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: tot_conc_rdisc
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: tot_conc_hdisc
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: time
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: conc_hist1,ages1
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: conc_hist2,ages2
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: conc_hist3,ages3
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: conc_hist4,ages4
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: conc_hist5,ages5
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: conc_hist6,ages6
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: ages_ball,ages_hexcyl
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: ages_rdisc,ages_hdisc
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: ages_hexcyl2,ages_hexcyl3
      REAL(KIND=r8) :: D0,Ea,init_conc(6),junk,min_diff
      CHARACTER :: therm_hist*100,conc_file*100
      CHARACTER :: line*1024,junk2*2
      INTEGER :: nt,nnode(6),i,diff_flag,j,alpha_ejec,diffuse
      INTEGER :: concFlag(6)
      REAL(KIND=r8) :: t1(61),t2(61),rad(6),height(3)

      ! Read in the number of nodes for each geometry
      ! to allocate the initial concentration history arrays
      open (53,file='ball/nnes_ball.dat',status='old')
      read (53,*) nnode(1),junk,junk
      close (53)

      open (53,file='hexcylinder/nnes_hexcylinder.dat',status='old')
      read (53,*) nnode(2),junk,junk
      close (53)

      open (53,file='hexcylinder_2to1/nnes_hexcylinder.dat',status='old')
      read (53,*) nnode(3),junk,junk
      close (53)

      open (53,file='hexcylinder_3to1/nnes_hexcylinder.dat',status='old')
      read (53,*) nnode(4),junk,junk
      close (53)

      open (53,file='rounddisc/nnee_rounddisc.dat',status='old')
      read (53,*) nnode(5),junk,junk
      close (53)

      open (53,file='hexdisc/nnee_hexdisc.dat',status='old')
      read (53,*) nnode(6),junk,junk
      close (53)

      allocate (conc_hist1(nnode(1)),conc_hist2(nnode(2)))
      allocate (conc_hist3(nnode(3)),conc_hist4(nnode(4)))
      allocate (conc_hist5(nnode(5)),conc_hist6(nnode(6)))

      ! Open the bambHe input file and writes to a scratch file without the comments
      open (54,file='bambHe.in',status='old')
      open (55,status='scratch')
    1 read (54,'(a1024)',end=2) line
      if (line(1:1).ne.'$'.and. line(1:1).ne.' ') write (55,'(a)') line
      goto 1
    2 close (54)
      rewind (55)

      ! Read in time and temperature history filename
      read (55,*) therm_hist

      ! Number of timesteps wanted for run
      read (55,*) nt

      ! Diffusivity and activation energy flag
      ! A value of 1 uses D0 and Ea from Farley, 2000
      ! A value of 2 allows the user to set the D0 and Ea in the input file
      read (55,*) diff_flag

      if (diff_flag.eq.1) then
        D0=50.                                        ! Diffusivity at infinite temp. [cm^2/s]
        D0=D0*0.01*0.01                               ! [m^2/s]
        Ea=33.0e3*4.184                               ! Activation energy in undamaged apatite grain [J/mol]
        read (55,*) junk
        read (55,*) junk
      else if (diff_flag.eq.2) then
        read (55,*) D0
        read (55,*) Ea
      endif

      ! Alpha ejection flag (1=on)
      read (55,*) alpha_ejec

      ! Diffusion flag (1=on)
      read (55,*) diffuse

      ! Read in minimum diffusivity value to allow diffusion to occur
      read (55,*) min_diff

      ! Flag for geometry calculation/initial He concentration
      ! 0 = No calculation on the specified geometry
      ! 1 = Geometry is calculated and user specifies one value
      ! for all nodes to be set to for He initial concentration
      ! 2 = Geometry is calculated and user specifies a file
      ! to read in initial He concentrations for all nodes
      ! Note: This file must be an output file from this program
      ! and should contain the same number of nodes as the current
      ! geometry
      do i=1,6
        conc_file=''
        read (55,*) concFlag(i)

        if (concFlag(i) .eq. 0) then                                            ! No calculation on geometry and basically skips the next two lines of input
          read (55,*) junk
          read (55,*) junk
          if (i.eq.2 .or. i.eq.3 .or. i.eq.4) read (55,*) junk
        else if (concFlag(i) .eq. 1) then
          read (55,*) init_conc(i)
          if (i.eq.1) conc_hist1=init_conc(i)                                   ! Initial concentration of sphere
          if (i.eq.2) conc_hist2=init_conc(i)                                   ! Initial concentration of hexagonal cylinder with 1:1 ratio
          if (i.eq.3) conc_hist3=init_conc(i)                                   ! Initial concentration of hexagonal cylinder with 2:1 ratio
          if (i.eq.4) conc_hist4=init_conc(i)                                   ! Initial concentration of hexagonal cylinder with 3:1 ratio
          if (i.eq.5) conc_hist5=init_conc(i)                                   ! Initial concentration of round disc
          if (i.eq.6) conc_hist6=init_conc(i)                                   ! Initial concentration of hexagonal disc

          read (55,*) rad(i)                                                    ! Radius or amplification factor for the specific geometry
          if (i.eq.2 .or. i.eq.3 .or. i.eq.4) read (55,*) height(i-1)           ! Height of hexagonal cylinder
        else if (concFlag(i) .eq. 2) then
          read (55,'(A100)') conc_file                                              ! Concentration history filename

          open (56,file=conc_file,status="old")
          read (56,*)
          read (56,*)
          read (56,*)
          read (56,*)
            do j=1,nnode(i)
              if (i.eq.1) then
                read (56,*) junk,junk,junk,conc_hist1(j),junk
              else if (i.eq.2) then
                read (56,*) junk,junk,junk,conc_hist2(j),junk
              else if (i.eq.3) then
                read (56,*) junk,junk,junk,conc_hist3(j),junk
              else if (i.eq.4) then
                read (56,*) junk,junk,junk,conc_hist4(j),junk
              else if (i.eq.5) then
                read (56,*) junk,junk,conc_hist5(j),junk
              else if (i.eq.6) then
                read (56,*) junk,junk,conc_hist6(j),junk
              endif
            enddo
          read (55,*) rad(i)
          if (i.eq.2 .or. i.eq.3 .or. i.eq.4) read (55,*) height(i-1)
        endif
      enddo

      close (55)
      close (56)


      allocate (ages_ball(nt),ages_hexcyl(nt))
      allocate (ages_rdisc(nt),ages_hdisc(nt))
      allocate (ages_hexcyl2(nt),ages_hexcyl3(nt))
      allocate (tot_conc_ball(nt),tot_conc_hexcyl(nt))
      allocate (tot_conc_rdisc(nt),tot_conc_hdisc(nt))
      allocate (tot_conc_hexcyl2(nt),tot_conc_hexcyl3(nt))
      allocate (time(0:nt))

      tot_conc_ball=0.0
      tot_conc_hexcyl=0.0
      tot_conc_hexcyl2=0.0
      tot_conc_hexcyl3=0.0
      tot_conc_hdisc=0.0
      tot_conc_rdisc=0.0
      ages_ball=0.0
      ages_hexcyl=0.0
      ages_hexcyl2=0.0
      ages_hexcyl3=0.0
      ages_rdisc=0.0
      ages_hdisc=0.0

      ! Calculate He concentration, ages, production rate, etc
      ! for the sphere if output specified
      if (concFlag(1).ne.0) then
        call MAINDRIVER_BALL(tot_conc_ball,time,D0,Ea,nt,therm_hist,&
             ages_ball,conc_hist1,nnode(1),rad(1),alpha_ejec,diffuse,min_diff)
      endif

      ! Calculate He concentration, ages, production rate, etc
      ! for the hexagonal cylinder (1:1) if output specified
      if (concFlag(2).ne.0) then
        call MAINDRIVER_HEXCYLINDER(tot_conc_hexcyl,D0,Ea,nt,&
             therm_hist,ages_hexcyl,conc_hist2,nnode(2),rad(2),height(1),&
             time,alpha_ejec,diffuse,min_diff)
      endif

      ! Calculate He concentration, ages, production rate, etc
      ! for the hexagonal cylinder (2:1) if output specified
      if (concFlag(3).ne.0) then
        call MAINDRIVER_HEXCYLINDER2(tot_conc_hexcyl2,D0,Ea,nt,&
             therm_hist,ages_hexcyl2,conc_hist3,nnode(3),rad(3),height(2),&
             time,alpha_ejec,diffuse,min_diff)
      endif

      ! Calculate He concentration, ages, production rate, etc
      ! for the hexagonal cylinder (3:1) if output specified
      if (concFlag(4).ne.0) then
        call MAINDRIVER_HEXCYLINDER3(tot_conc_hexcyl3,D0,Ea,nt,&
             therm_hist,ages_hexcyl3,conc_hist4,nnode(4),rad(4),height(3),&
             time,alpha_ejec,diffuse,min_diff)
      endif

      ! Calculate He concentration, ages, production rate, etc
      ! for the round disc if output specified
      if (concFlag(5).ne.0) then
        call MAINDRIVER_ROUNDDISC(tot_conc_rdisc,D0,Ea,nt,&
             therm_hist,ages_rdisc,conc_hist5,nnode(5),rad(5),&
             time,alpha_ejec,diffuse,min_diff)
      endif

      ! Calculate He concentration, ages, production rate, etc
      ! for the hexagonal disc if output specified
      if (concFlag(6).ne.0) then
        call MAINDRIVER_HEXDISC(tot_conc_hdisc,D0,Ea,nt,&
             therm_hist,ages_hdisc,conc_hist6,nnode(6),rad(6),&
             time,alpha_ejec,diffuse,min_diff)
      endif

      ! Output the total Helium domain volume for each geometry
      ! in one unformatted data file and one tecplot formatted file
      call print_He_volume (time,tot_conc_ball,tot_conc_hexcyl,&
           tot_conc_rdisc,tot_conc_hdisc,tot_conc_hexcyl2,tot_conc_hexcyl3,nt,concFlag)

      ! Output the Helium age for each geometry at each timestep in
      ! one unformatted data file and one tecplot formatted file
      call print_He_ages (time,ages_ball,ages_hexcyl,&
           ages_rdisc,ages_hdisc,ages_hexcyl2,ages_hexcyl3,nt,concFlag)

      deallocate (conc_hist1,conc_hist2,conc_hist3,conc_hist4,conc_hist5,conc_hist6)
      deallocate (ages_ball,ages_hexcyl,ages_rdisc,ages_hdisc,ages_hexcyl2,ages_hexcyl3)
      deallocate (tot_conc_ball,tot_conc_hexcyl)
      deallocate (tot_conc_rdisc,tot_conc_hdisc)
      deallocate (tot_conc_hexcyl2,tot_conc_hexcyl3)
      deallocate (time)

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Prints the total amount of Helium in each domain for every timestep
      subroutine print_He_volume (time,tot_conc_ball,tot_conc_hexcyl,&
           tot_conc_rdisc,tot_conc_hdisc,tot_conc_hexcyl2,tot_conc_hexcyl3,nt,concFlag)

      implicit none

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)

      INTEGER :: i,nt,concFlag(6)
      REAL(KIND=r8) :: tot_conc_ball(nt)
      REAL(KIND=r8) :: tot_conc_hexcyl(nt)
      REAL(KIND=r8) :: tot_conc_hexcyl2(nt)
      REAL(KIND=r8) :: tot_conc_hexcyl3(nt)
      REAL(KIND=r8) :: tot_conc_rdisc(nt)
      REAL(KIND=r8) :: tot_conc_hdisc(nt)
      REAL(KIND=r8) :: time(0:nt)

      open (10,file="He_domain_volume_hist.dat",status="unknown")
      write (10,'(A110)') 'Helium Domain Volume'
      write (10,'(A19)',advance="no") 'Time (s)'
      if (concFlag(1).ne.0) write (10,'(A29)',advance="no") ' Sphere (nmol)'
      if (concFlag(2).ne.0) write (10,'(A31)',advance="no") ' Hexagonal Cylinder 1:1 (nmol)'
      if (concFlag(3).ne.0) write (10,'(A)',advance="no") '  Hexagonal Cylinder 2:1 (nmol)'
      if (concFlag(4).ne.0) write (10,'(A)',advance="no") '  Hexagonal Cylinder 3:1 (nmol)'
      if (concFlag(5).ne.0) write (10,'(A)',advance="no") ' Round Disc (nmol)'
      if (concFlag(6).ne.0) write (10,'(A28)',advance="no") ' Hexagonal Disc (nmol)'
      write (10,*)

      open (11,file="He_domain_volume_hist_tec.dat",status="unknown")
      write (11,*) 'TITLE = "He Domain Volume"'
      write (11,'(A)',ADVANCE="no") ' VARIABLES = "Time (s)"'
      if (concFlag(1).ne.0) write (11,'(A)',ADVANCE="no") ' "Sphere (nmol)"'
      if (concFlag(2).ne.0) write (11,'(A)',ADVANCE="no") ' "Hexagonal Cylinder 1:1 (nmol)"'
      if (concFlag(3).ne.0) write (11,'(A)',ADVANCE="no") ' "Hexagonal Cylinder 2:1 (nmol)"'
      if (concFlag(4).ne.0) write (11,'(A)',ADVANCE="no") ' "Hexagonal Cylinder 3:1 (nmol)"'
      if (concFlag(5).ne.0) write (11,'(A)',ADVANCE="no") ' "Round Disc (nmol)"'
      if (concFlag(6).ne.0) write (11,'(A)',ADVANCE="no") ' "Hexagonal Disc (nmol)"'
      write (11,*)
      write (11,*) 'ZONE T="Helium Volume"'
      write (11,*) 'I=',nt,', J=1, K=1, ZONETYPE=Ordered'
      write (11,*) 'DATAPACKING=POINT'
      write (11,'(A)',ADVANCE="no") ' DT=(DOUBLE'

      do i=1,6
        if (concFlag(i).ne.0) write (11,'(A)',ADVANCE="no") ' DOUBLE'
      enddo
      write (11,'(A)') ')'

      do i=1,nt
        write (10,'(es20.5)',ADVANCE="no") time(i)
        if (concFlag(1).ne.0) write (10,'(es27.5)',ADVANCE="no") tot_conc_ball(i)
        if (concFlag(2).ne.0) write (10,'(es27.5)',ADVANCE="no") tot_conc_hexcyl(i)
        if (concFlag(3).ne.0) write (10,'(es27.5)',ADVANCE="no") tot_conc_hexcyl2(i)
        if (concFlag(4).ne.0) write (10,'(es27.5)',ADVANCE="no") tot_conc_hexcyl3(i)
        if (concFlag(5).ne.0) write (10,'(es27.5)',ADVANCE="no") tot_conc_rdisc(i)
        if (concFlag(6).ne.0) write (10,'(es27.5)',ADVANCE="no") tot_conc_hdisc(i)
        write (10,*)

        write (11,'(es15.5)',ADVANCE="no") time(i)
        if (concFlag(1).ne.0) write (11,'(es15.5)',ADVANCE="no") tot_conc_ball(i)
        if (concFlag(2).ne.0) write (11,'(es15.5)',ADVANCE="no") tot_conc_hexcyl(i)
        if (concFlag(3).ne.0) write (11,'(es15.5)',ADVANCE="no") tot_conc_hexcyl2(i)
        if (concFlag(4).ne.0) write (11,'(es15.5)',ADVANCE="no") tot_conc_hexcyl3(i)
        if (concFlag(5).ne.0) write (11,'(es15.5)',ADVANCE="no") tot_conc_rdisc(i)
        if (concFlag(6).ne.0) write (11,'(es15.5)',ADVANCE="no") tot_conc_hdisc(i)
        write (11,*)
      enddo

      close (10)
      close (11)

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Prints the Helium age for each domain at every timestep
      subroutine print_He_ages (time,ages_ball,ages_hexcyl,&
           ages_rdisc,ages_hdisc,ages_hexcyl2,ages_hexcyl3,nt,concFlag)

      implicit none

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)

      INTEGER :: i,nt,concFlag(6)
      REAL(KIND=r8) :: ages_ball(nt)
      REAL(KIND=r8) :: ages_hexcyl(nt)
      REAL(KIND=r8) :: ages_hexcyl2(nt)
      REAL(KIND=r8) :: ages_hexcyl3(nt)
      REAL(KIND=r8) :: ages_rdisc(nt)
      REAL(KIND=r8) :: ages_hdisc(nt)
      REAL(KIND=r8) :: time(0:nt)

      open (10,file="He_domain_ages_hist.dat",status="unknown")
      write (10,'(A110)') 'Helium Domain Ages'
      write (10,'(A19)',advance="no") 'Time (s)'
      if (concFlag(1).ne.0) write (10,'(A26)',advance="no") ' Sphere (Ma)'
      if (concFlag(2).ne.0) write (10,'(A31)',advance="no") ' Hexagonal Cylinder 1:1 (Ma)'
      if (concFlag(3).ne.0) write (10,'(A)',advance="no") ' Hexagonal Cylinder 2:1 (Ma)'
      if (concFlag(4).ne.0) write (10,'(A)',advance="no") ' Hexagonal Cylinder 3:1 (Ma)'
      if (concFlag(5).ne.0) write (10,'(A)',advance="no") ' Round Disc (Ma)'
      if (concFlag(6).ne.0) write (10,'(A25)',advance="no") ' Hexagonal Disc (Ma)'
      write (10,*)

      open (11,file="He_domain_ages_hist_tec.dat",status="unknown")
      write (11,*) 'TITLE = "He Domain Ages"'
      write (11,'(A)',ADVANCE="no") ' VARIABLES = "Time (s)"'
      if (concFlag(1).ne.0) write (11,'(A)',ADVANCE="no") ' "Sphere (Ma)"'
      if (concFlag(2).ne.0) write (11,'(A)',ADVANCE="no") ' "Hexagonal Cylinder 1:1 (Ma)"'
      if (concFlag(3).ne.0) write (11,'(A)',ADVANCE="no") ' "Hexagonal Cylinder 2:1 (Ma)"'
      if (concFlag(4).ne.0) write (11,'(A)',ADVANCE="no") ' "Hexagonal Cylinder 3:1 (Ma)"'
      if (concFlag(5).ne.0) write (11,'(A)',ADVANCE="no") ' "Round Disc (Ma)"'
      if (concFlag(6).ne.0) write (11,'(A)',ADVANCE="no") ' "Hexagonal Disc (Ma)"'
      write (11,*)
      write (11,*) 'ZONE T="Helium Ages"'
      write (11,*) 'I=',nt,', J=1, K=1, ZONETYPE=Ordered'
      write (11,*) 'DATAPACKING=POINT'
      write (11,'(A)',ADVANCE="no") ' DT=(DOUBLE'

      do i=1,6
        if (concFlag(i).ne.0) write (11,'(A)',ADVANCE="no") ' DOUBLE'
      enddo
      write (11,'(A)') ')'

      do i=1,nt
        write (10,'(es20.5)',ADVANCE="no") time(i)
        if (concFlag(1).ne.0) write (10,'(es25.5)',ADVANCE="no") ages_ball(i)/3600./24./365.25e6
        if (concFlag(2).ne.0) write (10,'(es25.5)',ADVANCE="no") ages_hexcyl(i)/3600./24./365.25e6
        if (concFlag(3).ne.0) write (10,'(es25.5)',ADVANCE="no") ages_hexcyl2(i)/3600./24./365.25e6
        if (concFlag(4).ne.0) write (10,'(es25.5)',ADVANCE="no") ages_hexcyl3(i)/3600./24./365.25e6
        if (concFlag(5).ne.0) write (10,'(es25.5)',ADVANCE="no") ages_rdisc(i)/3600./24./365.25e6
        if (concFlag(6).ne.0) write (10,'(es25.5)',ADVANCE="no") ages_hdisc(i)/3600./24./365.25e6
        write (10,*)

        write (11,'(es15.5)',ADVANCE="no") time(i)
        if (concFlag(1).ne.0) write (11,'(es15.5)',ADVANCE="no") ages_ball(i)/3600./24./365.25e6
        if (concFlag(2).ne.0) write (11,'(es15.5)',ADVANCE="no") ages_hexcyl(i)/3600./24./365.25e6
        if (concFlag(3).ne.0) write (11,'(es15.5)',ADVANCE="no") ages_hexcyl2(i)/3600./24./365.25e6
        if (concFlag(4).ne.0) write (11,'(es15.5)',ADVANCE="no") ages_hexcyl3(i)/3600./24./365.25e6
        if (concFlag(5).ne.0) write (11,'(es15.5)',ADVANCE="no") ages_rdisc(i)/3600./24./365.25e6
        if (concFlag(6).ne.0) write (11,'(es15.5)',ADVANCE="no") ages_hdisc(i)/3600./24./365.25e6
        write (11,*)
      enddo
      close (10)
      close (11)

      return
      end subroutine
