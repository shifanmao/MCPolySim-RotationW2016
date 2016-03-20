!---------------------------------------------------------------*
      
PROGRAM wlcsim
      
!     
!     This simulation tracks the dynamics of a single polymer
!     chain modeled as a discrete wormlike chain with bending
!     and stretching energy.
!     
!     Andrew Spakowitz
!     Written 9-2-13
!     
      
!     Kevin J Hou
!
!     Rotation Winter 2016
!     Modifications made to WLC simulation to account for properties of Conjugated Polymers
!         1. Twist Resistance 
!   
!     Search text 'KHOU' to find changes
!


!     Variables within the simulation

  use mt19937, only : grnd, sgrnd, rnorm, mt, mti

  PARAMETER (PI=3.141592654) ! Value of pi
  INTEGER PTON                  ! Parallel Tempering on
  INTEGER NRABOVE                 ! Total number of replicas  
  INTEGER NRBELOW                 ! Total number of replicas  
  INTEGER NRNOW                 ! Total number of replicas  
  INTEGER NRMAX                 ! Total number of replicas  
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: R	 ! Conformation of polymer chains
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: U	 ! Conformation of polymer chains
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: R0	 ! Conformation of polymer chains
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: U0	 ! Conformation of polymer chains
    
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:):: SijU   ! Tensoral Interaction
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:):: SijW   ! Tensoral Interaction
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:):: dSijU   ! Tensoral Interaction
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:):: dSijW   ! Tensoral Interaction
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: SijE        ! Double dot of SijU::SijU

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: PHIA ! Volume fraction of A
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: PHIB ! Volume fraction of B
  INTEGER, ALLOCATABLE, DIMENSION(:):: AB            ! Chemical identity of beads

  INTEGER NT                ! Number of beads in simulation
  INTEGER N                 ! Number of beads in simulation
  INTEGER NB                 ! Number of beads in simulation
  INTEGER NP                ! Number of polymers in simulation
  DOUBLE PRECISION ENERGY   ! Total energy
  DOUBLE PRECISION DT       ! Time step size
  INTEGER I,J,IB            ! Index
  INTEGER INDMAX            ! Maximum index in series
  INTEGER IND               ! Ind in series
  INTEGER TENS              ! Decimal of index
  character*8 fileind       ! Index of output
  character*4 repind        ! Replica index for output
  character*16 snapnm       ! File for output
  INTEGER INDEND            ! Restart index
  logical restart           ! Restart from previous?
  INTEGER WRTON

!     Simulation input variables
  
  INTEGER FRMFILE           ! Initial condition
  INTEGER BROWN             ! Include Brownian forces
  INTEGER INTON             ! Include polymer interactions
  INTEGER NSTEP				! Number of MC steps between save
  INTEGER NINIT				! Number of initialization MC steps
  INTEGER NNOINT			! Number of initialization MC steps without interactions
  INTEGER NCHI			        ! Number of savepoints between chi change
  INTEGER NPT			        ! Number of savepoints between tempering
  INTEGER FRMCHEM           ! Initial chemical sequence

!     Monte Carlo variables (KHOU - > 7)

  DOUBLE PRECISION MCAMP(7) ! Amplitude of random change
  INTEGER MOVEON(7)			! Is the move active 
  INTEGER WINDOW(7)			! Size of window for bead selection 
  INTEGER SUCCESS(7)        ! Number of successes
  DOUBLE PRECISION PHIT(7)     ! % hits per total steps
      
!     Energy variables
      
  DOUBLE PRECISION EELAS(3) ! Elastic force
  DOUBLE PRECISION ECHI   ! CHI energy
  DOUBLE PRECISION EKAP   ! KAP energy
  DOUBLE PRECISION ETWIST ! Twist Energy
  DOUBLE PRECISION EALIGN ! Alignment Energy 
  DOUBLE PRECISION ETOTAL ! Total Energy
      
!     Structure analysis
      
  DOUBLE PRECISION RCOM(3)  ! Center of mass
  DOUBLE PRECISION DELR(3)  ! Mag of gyration tensor
  DOUBLE PRECISION RCOM0(3) ! Init val RCOM
  DOUBLE PRECISION DELR0(3) ! Init val DELR
  DOUBLE PRECISION DRCOM    ! Change in RCOM
  DOUBLE PRECISION SIG(3,3)
  DOUBLE PRECISION COR
  
  INTEGER  SON               !calculate Structure Factors
  INTEGER, PARAMETER:: XNUM = 11
  INTEGER, PARAMETER:: KNUM = (XNUM)**3
  INTEGER NVEC              !number of times calculating SVEC
  DOUBLE PRECISION KVEC(KNUM)
  DOUBLE PRECISION SVEC(KNUM)
      
!     Variables in the simulation
      
  DOUBLE PRECISION PARA(10)
  INTEGER NBIN              ! Number of bins
  INTEGER NBINX              ! Number of bins

!     Variables for the random number generators

  INTEGER IDUM              ! Seed for the generator
  DOUBLE PRECISION MOM(6)

!     Simulation parameters
  
  INTEGER G					! Beads per monomer
  DOUBLE PRECISION LBOX		! Box length (approximate)
  DOUBLE PRECISION DEL      ! Discretization size (approximate)
  DOUBLE PRECISION V		! Bead volume
  DOUBLE PRECISION FA		! Fraction of A beads
  DOUBLE PRECISION LAM		! Chemical correlation parameter
  DOUBLE PRECISION EPS		! Elasticity l0/(2lp)
  DOUBLE PRECISION CHI       ! Chi parameter value
  DOUBLE PRECISION AU,AUV,AV  ! Parameters Describing Alignment
  DOUBLE PRECISION DCHI       ! Chi parameter value
  DOUBLE PRECISION CHI0         ! Initial Chi before ramping
  DOUBLE PRECISION DELCHI       ! Chi ramping rate
  DOUBLE PRECISION KAP		! Incompressibility parameter
  DOUBLE PRECISION TAU    ! Parameter describing twist energy (KHOU)
  DOUBLE PRECISION L0       ! Equilibrium segment length
  INTEGER PTID              ! ID to pair up replicas for PT
  INTEGER ACCBELOW


! KHOU
  ! W tracks twist of chain. W is always orthogonal to corresponding U vector.   
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: W  
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: W0

!     Load in the parameters for the simulation

  open (unit=5, file='input/input')
  read (unit=5, fmt='(4(/))')
  read (unit=5, fmt=*) PTON
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) N
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) G
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) LBOX
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) DEL
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) L0
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) V
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) FA
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) LAM
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) FRMCHEM
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) EPS
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) CHI0
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) DELCHI
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) KAP
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) INDMAX
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) FRMFILE
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) BROWN
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) INTON
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NCHI
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NPT
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NNOINT
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NINIT
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NSTEP
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) WRTON
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) PTID
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NRABOVE
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NRBELOW
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NRNOW
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NRMAX


  !KHOU
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) TAU
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) AU
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) AUV
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) AV
  close(5)

  NP=nint(LBOX**3./(N*G*V))
  LBOX=(V*N*G*NP)**(1./3.)
  call getpara(PARA,EPS,L0,LBOX)
  
  NT=N*NP*G
  ALLOCATE(R(NT,3))
  ALLOCATE(U(NT,3))
  ALLOCATE(R0(NT,3))
  ALLOCATE(U0(NT,3))
  ALLOCATE(AB(NT))

  ! KHOU
  ALLOCATE(W(NT,3))
  ALLOCATE(W0(NT,3))


  NBINX=nint(LBOX/DEL)
  NBIN=NBINX**3.
  DEL=LBOX/NBINX
  
  ALLOCATE(PHIA(NBIN))
  ALLOCATE(PHIB(NBIN))
  ALLOCATE(SijU(NBIN,3,3))
  ALLOCATE(SijW(NBIN,3,3))
  ALLOCATE(dSijU(NBIN,3,3))
  ALLOCATE(dSijW(NBIN,3,3))
  ALLOCATE(SijE(NBIN))

!  Monte-Carlo simulation parameters
  MCAMP(1)=0.5*PI
  MCAMP(2)=0.3*L0
  MCAMP(3)=0.5*PI
  MCAMP(4)=0.5*PI
  MCAMP(5)=0.5*PI
  MCAMP(6)=5.0*L0
  MOVEON(1)=1
  MOVEON(2)=1
  MOVEON(3)=1
  MOVEON(4)=1
  if (INTON.EQ.1) then
     MOVEON(5)=1
     MOVEON(6)=1
  else
     MOVEON(5)=1
     MOVEON(6)=1
  endif

!     Initial segment window for MC moves
  WINDOW(1)=N*G
  WINDOW(2)=N*G
  WINDOW(3)=N*G
  WINDOW(4)=N*G
  WINDOW(5)=N*G
  WINDOW(6)=N*G

! KHOU
  MOVEON(7)=1
  WINDOW(7)=N*G
  MCAMP(7)=0.5*PI

  INQUIRE (FILE = 'data/out1', exist = restart)
  if (.NOT.restart) then

    PRINT*, '-----new simulation-----'

!     Setup the initial condition

    NB=N*G

    ! KHOU
    call initcond(R,U,W,AB,NT,NB,NP,IDUM,FRMFILE,PARA,LBOX)

!     Load in AB sequence
    IF (FRMCHEM.EQ.1) THEN
       OPEN (UNIT = 1, FILE = 'input/ab', STATUS = 'OLD')      
       IB=1
       DO 1 I=1,NP
          DO 2 J=1,NB
             READ(1,"(I2)") AB(IB)
             IB=IB+1
2         CONTINUE
1      CONTINUE 
       CLOSE(1)
    ELSE
       call initchem(AB,NT,N,G,NP,FA,LAM)
    ENDIF 

!     Perform an initialization MC simulation

    !KHOU
    SON=0
    call MCsim(R,U,W,PHIA,PHIB,AB,NT,NB,NP,NBIN,NNOINT,BROWN,0,PARA,V,CHI0,KAP,TAU,AU,AUV,AV,LBOX,L0,DEL,MCAMP,&
         SUCCESS,MOVEON,WINDOW,PHIT,KVEC,SVEC,SON,0,IND,NRABOVE,NRBELOW,NRMAX,NRNOW,ECHI,NPT,PTID, &
         SijU,SijW,dSijU,dSijW,EAlIGN,SijE)
    call MCsim(R,U,W,PHIA,PHIB,AB,NT,NB,NP,NBIN,NINIT,BROWN,INTON,PARA,V,CHI0,KAP,TAU,AU,AUV,AV,LBOX,L0,DEL,MCAMP,&
         SUCCESS,MOVEON,WINDOW,PHIT,KVEC,SVEC,SON,0,IND,NRABOVE,NRBELOW,NRMAX,NRNOW,ECHI,NPT,PTID, &
         SijU,SijW,dSijU,dSijW,EALIGN,SijE)

!     Save initial conformation and PSI angles 
    OPEN (UNIT = 1, FILE = 'data/r0', STATUS = 'NEW')      
    IB=1
    DO 10 I=1,NP
       DO 20 J=1,NB
          R0(IB,1)=R(IB,1)
          R0(IB,2)=R(IB,2)
          R0(IB,3)=R(IB,3)
          U0(IB,1)=U(IB,1)
          U0(IB,2)=U(IB,2)
          U0(IB,3)=U(IB,3)

          ! KHOU
          W0(IB,1)=W(IB,1)
          W0(IB,2)=W(IB,2)
          W0(IB,3)=W(IB,3)
          
          WRITE(1,"(3f8.3,I2)") R(IB,1),R(IB,2),R(IB,3),AB(IB)
          IB=IB+1
20     CONTINUE
10  CONTINUE 
    CLOSE(1)
      
    OPEN (UNIT = 1, FILE = 'data/u0', STATUS = 'NEW')
    IB=1
    DO 30 I=1,NP
       DO 40 J=1,NB
          WRITE(1,"(3f8.3,I2)") U(IB,1),U(IB,2),U(IB,3),AB(IB)
          IB=IB+1
40     CONTINUE
30  CONTINUE 
    CLOSE(1)

!KHOU
    OPEN (UNIT = 1, FILE = 'data/w0', STATUS = 'NEW')
    IB=1
    DO 35 I=1,NP
       DO 45 J=1,NB
          WRITE(1,"(3f8.3,I2)") W(IB,1),W(IB,2),W(IB,3),AB(IB)
          IB=IB+1
45     CONTINUE
35  CONTINUE 
    CLOSE(1)


!     Open the output files
    INDEND = 0
    OPEN (UNIT = 1, FILE = 'data/out1', STATUS = 'NEW')
    OPEN (UNIT = 2, FILE = 'data/out2', STATUS = 'NEW')
    OPEN (UNIT = 3, FILE = 'data/out3', STATUS = 'NEW')
    
    ! KHOU
    CLOSE(1)
    CLOSE(2)
    CLOSE(3)
    

 else

    PRINT*, '-----load simulation-----'
    ! KHOU - TODO: Fix load_from_old
    PRINT*, ' WARNING: Modifications for SCP not implemented for load_from_old -KH 2016 '
    NB=N*G
    call load_from_old(R,U,AB,CHI,DCHI,NT,NB,NP,IDUM,FRMFILE,PARA,LBOX,INDEND)

 endif



!     Begin simulation
  IND=1
  TIME=0.

  !part 7 - PT
  OPEN(unit=1,file='data/ptnow',IOSTAT=IOStatus,status='new')
  CLOSE(unit=1,IOSTAT=IOStatus)
  OPEN(unit=1,file='data/calcpnow',IOSTAT=IOStatus,status='new')
  CLOSE(unit=1,IOSTAT=IOStatus)
  OPEN(unit=1,file='data/swapnow',IOSTAT=IOStatus,status='new')
  CLOSE(unit=1,IOSTAT=IOStatus)
  OPEN(unit=1,file='data/swapend',IOSTAT=IOStatus,status='new')
  CLOSE(unit=1,IOSTAT=IOStatus)

  DO WHILE ((IND+INDEND).LE.INDMAX)
!     Parallel tempering + chi annealing
     call chisched(IND,NR,CHI,DCHI,CHI0,DELCHI,NCHI,NPT,restart)

!     Perform a MC simulation
     SON=0

     call MCsim(R,U,W,PHIA,PHIB,AB,NT,NB,NP,NBIN,NSTEP,BROWN,INTON,PARA,V,CHI,KAP,TAU,AU,AUV,AV,LBOX,L0,DEL,MCAMP,&
          SUCCESS,MOVEON,WINDOW,PHIT,KVEC,SVEC,SON,PTON,IND,NRABOVE,NRBELOW,NRMAX,NRNOW,ECHI,NPT,PTID, &
          SijU,SijW,dSijU,dSijW,EALIGN,SijE)

!     Save the conformation and the metrics
     TENS=nint(log10(1.*INDEND+IND)-0.4999)+1
     write (fileind,'(I4)'), INDEND+IND

      !part 1 - energy
     ECHI=0.
     EKAP=0
     EELAS(1)=0.
     EELAS(2)=0.
     EELAS(3)=0.
     call energy_elas(EELAS(1:3),R,U,NT,NB,NP,PARA)
     call energy_twist(ETWIST,R,U,W,NT,N,NP,TAU)
     call calc_Sij(R,U,W,NT,N,NP,NBIN,&
      V,LBOX,DEL,SijU,SijW,AU,AUV,AV,EALIGN,SijE)

     if (INTON.EQ.1) then
        call energy_int(R,AB,NP,NB,NT,NBIN,V,CHI,KAP,LBOX,DEL,ECHI,EKAP)
     endif

     ETOTAL = EELAS(1)+EELAS(2)+EELAS(3)+ECHI+EKAP+ETWIST+EALIGN

     OPEN (UNIT = 1, FILE = 'data/out1', STATUS = 'OLD')
     DO K = 1,IND+INDEND-1
        READ(1,*)
     ENDDO
     WRITE(1,"(I6,8f20.3)") INDEND+IND,EELAS(1:3),ECHI,EKAP,ETWIST,EALIGN,ETOTAL
     CLOSE(1)
      
     !part 2 - CHI
     OPEN (UNIT = 2, FILE = 'data/out2', STATUS = 'OLD')
     DO K = 1,IND+INDEND-1
        READ(2,*)
     ENDDO
     !  WRITE(2,"(I15,1f10.2)") INDPT,CHI
     WRITE(2,"(I6,1f20.8)") INDEND+IND,CHI*G
     CLOSE(2)

     !part 2.5 - adaptations
     OPEN (UNIT = 3, FILE = 'data/out3', STATUS = 'OLD')
     DO K = 1,IND+INDEND-1
        READ(3,*)
     ENDDO
     WRITE(3,"(I6,18f8.2)") INDEND+IND,REAL(WINDOW(1)),MCAMP(1),PHIT(1), &
          REAL(WINDOW(2)),MCAMP(2),PHIT(2), &
          REAL(WINDOW(3)),MCAMP(3),PHIT(3), &
          REAL(MOVEON(4)),MCAMP(4),PHIT(4), &
          REAL(MOVEON(5)),MCAMP(5),PHIT(5), &
          REAL(MOVEON(6)),MCAMP(6),PHIT(6)
     CLOSE(3)
      
     IF ((WRTON.EQ.1).AND.(MOD(IND+INDEND,250).EQ.1)) THEN

     !part 3 - R
     snapnm= 'data/r'//fileind((4-TENS+1):4)
     OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
     IB=1
     DO 50 I=1,NP
        DO 60 J=1,NB
           WRITE(1,"(3f7.2,I4)") , &
                R(IB,1)-0.*nint(R(IB,1)/LBOX-0.5)*LBOX, &
                R(IB,2)-0.*nint(R(IB,2)/LBOX-0.5)*LBOX, &
                R(IB,3)-0.*nint(R(IB,3)/LBOX-0.5)*LBOX, I
           IB=IB+1
60         CONTINUE
50      CONTINUE
     CLOSE(1)

     !part 4 - U 
    snapnm= 'data/u'//fileind((4-TENS+1):4)
    OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
    IB=1
    DO 70 I=1,NP
       DO 80 J=1,NB
          WRITE(1,"(3f7.3,I2)") U(IB,1),U(IB,2),U(IB,3),AB(IB)
          IB=IB+1
80         CONTINUE
70         CONTINUE 
    CLOSE(1)

    !part 4b - W 
    snapnm= 'data/w'//fileind((4-TENS+1):4)
    OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
    IB=1
    DO 75 I=1,NP
       DO 85 J=1,NB
          WRITE(1,"(3f7.3,I2)") W(IB,1),W(IB,2),W(IB,3),AB(IB)
          IB=IB+1
85         CONTINUE
75         CONTINUE 
    CLOSE(1)


     !part 5 - PHI
!     snapnm= 'data/p'//fileind((4-TENS+1):4)
!     OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
!     DO K=1,NBIN*
!        WRITE(1,"(2f7.2)") PHIA(K),PHIB(K)
!     ENDDO
!     CLOSE(1)

     !TODO: PART 6 - SijU Output
    !snapnm= 'data/SijE'//fileind((4-TENS+1):4)
    !OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
    !DO K=1,NBIN
    !  WRITE(1,"(I6, 1f22.8)") K, SijE(K)
    !ENDDO
    !CLOSE(1)

     ENDIF

     !part 6 - S(k)
     IF (SON.EQ.1) THEN
        snapnm= 'data/s'//fileind((4-TENS+1):4)
        OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
        DO K=1,KNUM
           WRITE(1,"(2f20.3)") KVEC(K),SVEC(K)
        ENDDO
        CLOSE(1)
     ENDIF

     PRINT*, '________________________________________'
     PRINT*, 'Time point ',IND+INDEND, ' out of', INDMAX
     PRINT*, 'Bending energy ', EELAS(1)
     PRINT*, 'Par compression energy ', EELAS(2)
     PRINT*, 'Perp compression energy ', EELAS(3)
     ! PRINT*, 'Chi interaction energy ', ECHI
     PRINT*, 'Kap compression energy ', EKAP
     PRINT*, 'Twist energy', ETWIST
     PRINT*, 'Alignment Energy', EALIGN
     PRINT*, 'Total Energy', ETOTAL
     !   PRINT*, 'Current number of beads ', NT
     !   PRINT*, 'Number of polymers ', NP
         
     IND=IND+1
         
  ENDDO
  
END
      
!---------------------------------------------------------------*
