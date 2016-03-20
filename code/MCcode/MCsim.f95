!---------------------------------------------------------------*
      
!     This subroutine performs a Monte Carlo simulation on the 
!     polymer chain.
      
      SUBROUTINE MCsim(R,U,W,PHIA,PHIB,AB,NT,N,NP,NBIN, &
          NSTEP,BROWN,INTON,PARA,V,CHI,KAP,TAU,AU,AUV,AV,LBOX,L0,& 
          DEL,MCAMP,SUCCESS,MOVEON,WINDOW,PHIT,KVEC,SVEC,SON, &
          PTON,IND,NRABOVE,NRBELOW,NRMAX,NRNOW,ECHI,NPT,PTID, &
          SijU,SijW,dSijU,dSijW,EALIGN,SijE)

      use mt19937, only : grnd, sgrnd, rnorm, mt, mti
      
      PARAMETER (PI=3.141592654) ! Value of pi

      INTEGER SON               !calculate Structure Factors
      INTEGER, PARAMETER:: XNUM = 11
      INTEGER, PARAMETER:: KNUM = (XNUM)**3
      INTEGER, PARAMETER:: NVEC = 50             !number of times calculating SVEC
      DOUBLE PRECISION KVEC(KNUM)
      DOUBLE PRECISION SVEC(KNUM)

      INTEGER PTON                  ! Parallel Tempering on
      INTEGER IND               ! Ind in series
      INTEGER INDPT               ! Ind in series
      INTEGER NRABOVE                 ! Total number of replicas  
      INTEGER NRBELOW                 ! Total number of replicas  
      INTEGER NRNOW                 ! Total number of replicas  
      INTEGER NRMAX                 ! Total number of replicas  
      INTEGER NPT			        ! Number of savepoints between tempering
      INTEGER PTID              ! ID to pair up replicas for PT
      INTEGER ACCBELOW

      DOUBLE PRECISION KMAG(KNUM)
      DOUBLE PRECISION SK(KNUM)
      DOUBLE PRECISION R(NT,3)  ! Bead positions
      DOUBLE PRECISION U(NT,3)  ! Tangent vectors
      DOUBLE PRECISION W(NT,3)  ! Twist vectors (KHOU)
      INTEGER AB(NT)            ! Chemical identity of beads
      DOUBLE PRECISION RP(NT,3)  ! Bead positions
      DOUBLE PRECISION UP(NT,3)  ! Tangent vectors
      DOUBLE PRECISION WP(NT,3)  ! Twist vectors (KHOU)
      INTEGER N,NP,NT           ! Number of beads
      INTEGER NBIN              ! Number of bins
      INTEGER NSTEP             ! Number of MC steps
      INTEGER BROWN            ! Turn on fluctuations
      INTEGER INTON             ! Include polymer interactions
      DOUBLE PRECISION ECHI   ! CHI energy
      DOUBLE PRECISION EKAP   ! KAP energy
      
!     Variables for the simulation
      
      INTEGER ISTEP             ! Current MC step index
      DOUBLE PRECISION PROB     ! Calculated test prob
      DOUBLE PRECISION TEST     ! Random test variable
      INTEGER IB                ! Test bead 
      INTEGER IP                ! Test polymer 
      INTEGER IB1               ! Test bead position 1
      INTEGER IT1               ! Index of test bead 1
      INTEGER IB2               ! Test bead position 2
      INTEGER IT2               ! Index of test bead 2

      INTEGER I,J,K
      DOUBLE PRECISION R0(3)
      
!     Energy variables
      
      DOUBLE PRECISION DEELAS   ! Change in bending energy
      DOUBLE PRECISION DEINT    ! Change in self energy
      DOUBLE PRECISION DEEX     ! Change in external energy
      DOUBLE PRECISION DETWIST  ! Change in twist energy (KHOU)
      DOUBLE PRECISION ENERGY
      DOUBLE PRECISION EALIGN  ! Need this to initialize the Sij calc
      
!     MC adaptation variables (KHOU)
      
      DOUBLE PRECISION MCAMP(7) ! Amplitude of random change      
      INTEGER MCTYPE            ! Type of MC move
      INTEGER NADAPT(7)         ! Num steps btwn adapt
      DOUBLE PRECISION PHIT(7)     ! % hits per total steps
      DOUBLE PRECISION PDESIRE(7) ! Desired hit rate
      INTEGER SUCCESS(7)        ! Number of successes
      DOUBLE PRECISION MINAMP(7) ! Minimum amp to stop
      DOUBLE PRECISION MAXAMP(7) ! Minimum amp to stop
      INTEGER MOVEON(7)			! Is the move active
      INTEGER WINDOW(7)			! Size of window for bead selection
      
!     Variables in the simulation
      
      DOUBLE PRECISION EB,EPAR,EPERP
      DOUBLE PRECISION GAM,ETA,TAU
      DOUBLE PRECISION AU,AUV,AV
      DOUBLE PRECISION XIR,XIU
      DOUBLE PRECISION LHC      ! Length of HC int
      DOUBLE PRECISION VHC      ! HC strength
      DOUBLE PRECISION PARA(10)
      DOUBLE PRECISION VA       ! Volume of monomer A
      DOUBLE PRECISION VB       ! Volume of monomer B
      DOUBLE PRECISION V       ! Volume of monomer
      DOUBLE PRECISION CHI      ! Chi parameter value
      DOUBLE PRECISION KAP      ! Compressibility value
      DOUBLE PRECISION LBOX     ! Simulation box size (approximate)
      DOUBLE PRECISION DEL      ! Discretization size (approximate)
      DOUBLE PRECISION L0       ! Equilibrium segment length

!     Variables for density calculation
      
      DOUBLE PRECISION PHIA(NBIN) ! Volume fraction of A
      DOUBLE PRECISION PHIB(NBIN) ! Volume fraction of B
      DOUBLE PRECISION DPHIA(NBIN) ! Volume fraction of A
      DOUBLE PRECISION DPHIB(NBIN) ! Volume fraction of B
      INTEGER INDPHI(NBIN)      ! Indices of the phi
      INTEGER INDSij(NBIN)      ! Indices of Sij
      INTEGER NPHI		! Number of phi values that change

!     Placeholder Variables for Alignment Energy
      DOUBLE PRECISION DEALIGN   ! Change in Alignment Energy
      DOUBLE PRECISION EALIGNOLD ! Alignment Energy before move
      DOUBLE PRECISION EALIGNNEW ! Alignment Energy after move
      DOUBLE PRECISION SijU(NBIN,3,3)  ! Order Parameter Sij(u) before move
      DOUBLE PRECISION SijW(NBIN,3,3)  ! Order Parameter Sij(v) before move
      DOUBLE PRECISION dSijU(NBIN,3,3) ! Order Parameter Sij(u) after move
      DOUBLE PRECISION dSijW(NBIN,3,3) ! Order Parameter Sij(v) after move
      DOUBLE PRECISION SijE(NBIN)     ! Energy per bin 
      
!     Load the input parameters
      
      EB=PARA(1)
      EPAR=PARA(2)
      EPERP=PARA(3)
      GAM=PARA(4)
      ETA=PARA(5)
      XIR=PARA(6)
      XIU=PARA(7)
      LHC=PARA(9)
      VHC=PARA(10)
	  
      MINAMP(1)=0.2*PI
      MINAMP(2)=0.2*L0
      MINAMP(3)=0.2*PI
      MINAMP(4)=0.2*PI
      MINAMP(5)=0.2*PI
      MINAMP(6)=0.2*L0
	  
      MAXAMP(1)=1.0*PI
      MAXAMP(2)=1.0*L0
      MAXAMP(3)=1.0*PI
      MAXAMP(4)=1.0*PI
      MAXAMP(5)=1.0*PI
      MAXAMP(6)=0.1*LBOX

      NADAPT(1)=1000
      NADAPT(2)=1000
      NADAPT(3)=1000
      NADAPT(4)=1000
      NADAPT(5)=1000
      NADAPT(6)=1000
      if (NSTEP.LE.NADAPT(1)) then
         NADAPT(1)=NSTEP
      endif
      if (NSTEP.LE.NADAPT(2)) then
         NADAPT(2)=NSTEP
      endif
      if (NSTEP.LE.NADAPT(3)) then
         NADAPT(3)=NSTEP
      endif
      if (NSTEP.LE.NADAPT(4)) then
         NADAPT(4)=NSTEP
      endif
      if (NSTEP.LE.NADAPT(5)) then
         NADAPT(5)=NSTEP
      endif
      if (NSTEP.LE.NADAPT(6)) then
         NADAPT(6)=NSTEP
      endif

      !KHOU 
      NADAPT(7)=1000
      if (NSTEP.LE.NADAPT(7)) then
         NADAPT(7)=NSTEP
      endif
      
      PDESIRE(1)=0.5
      PDESIRE(2)=0.5
      PDESIRE(3)=0.5
      PDESIRE(4)=0.5
      PDESIRE(5)=0.5
      PDESIRE(6)=0.5

      SUCCESS(1)=0
      SUCCESS(2)=0
      SUCCESS(3)=0
      SUCCESS(4)=0
      SUCCESS(5)=0
      SUCCESS(6)=0

      !KHOU
      PDESIRE(7)=0.5
      SUCCESS(7)=0
      MAXAMP(7) = 1.0*PI
      MINAMP(7) = 0.2*PI

      DEELAS=0.
      DEINT=0.
      DEALIGN=0.

      if (INTON.EQ.1) then
         call r_to_phi(R,AB,NT,N,NP,NTOT,NBIN, &
              V,CHI,KAP,LBOX,DEL,PHIA,PHIB)

         call calc_Sij(R,U,W,NT,N,NP,NBIN,&
              V,LBOX,DEL,SijU,SijW,AU,AUV,AV,EALIGN,SijE)
      endif

      

!     Begin Monte Carlo simulation
      
      ISTEP=1
      DO K=1,KNUM
         KVEC(K)=0         
         SVEC(K)=0
      ENDDO

      DO WHILE (ISTEP.LE.NSTEP)
         
         !KHOU 
         DO 10 MCTYPE=1,7
            
            if (MOVEON(MCTYPE).EQ.0) then
               goto 60
            endif

            ! KHOU
            call MC_move(R,U,W,RP,UP,WP,NT,N,NP,IP,IB1,IB2,IT1,IT2,MCTYPE,MCAMP,WINDOW)
            
!     Calculate the change in compression and bending energy
            
            call MC_eelas(DEELAS,R,U,RP,UP,NT,N,NP,IP,IB1,IB2,IT1,IT2,EB,EPAR,EPERP,GAM,ETA)
            
!     Calculate the change in the self-interaction energy
            if (INTON.EQ.1) then
               call MC_int(DEINT,R,AB,NT,NBIN, &
                    V,CHI,KAP,LBOX,DEL,PHIA,PHIB,DPHIA,DPHIB, &
                    INDPHI,NPHI,RP,IT1,IT2)

               call MC_align(DEALIGN,R,U,W,RP,UP,WP,NT,N,NP,NTOT,NBIN,V,LBOX,DEL,&
                    AU,AUV,AV,SijU,SijW,dSijU,dSijW,IT1,IT2,INDSij)
            endif

!     KHOU: Calculate the change in Twist Energy

            call MC_twist(DETWIST,R,U,W,RP,UP,WP,NT,N,NP,IT1,IT2,TAU)


!     Change the position if appropriate
            ENERGY=DEELAS+DEINT+DETWIST+DEALIGN
            
            PROB=exp(-ENERGY)

            if (BROWN.EQ.1) then
               TEST=grnd()
            else
               TEST=1.
            endif

            if (TEST.LE.PROB) then ! Move was made
               DO 20 I=IT1,IT2
                  R(I,1)=RP(I,1)
                  R(I,2)=RP(I,2)
                  R(I,3)=RP(I,3)
                  U(I,1)=UP(I,1)
                  U(I,2)=UP(I,2)
                  U(I,3)=UP(I,3)

                  !KHOU
                  W(I,1)=WP(I,1)
                  W(I,2)=WP(I,2)
                  W(I,3)=WP(I,3)
 20            CONTINUE

               if (INTON.EQ.1) then
                  DO 30 I=1,NPHI
                     J=INDPHI(I)
                     PHIA(J)=PHIA(J)+DPHIA(I)
                     PHIB(J)=PHIB(J)+DPHIB(I)	

                     ! Update SijU and SijW if a move was made
                     J = INDSij(I)
                     SijU(J,1,1) = SijU(J,1,1) + dSijU(I,1,1)
                     SijU(J,1,2) = SijU(J,1,2) + dSijU(I,1,2)
                     SijU(J,1,3) = SijU(J,1,3) + dSijU(I,1,3)
                     SijU(J,2,1) = SijU(J,2,1) + dSijU(I,2,1)
                     SijU(J,2,2) = SijU(J,2,2) + dSijU(I,2,2)
                     SijU(J,2,3) = SijU(J,2,3) + dSijU(I,2,3)
                     SijU(J,3,1) = SijU(J,3,1) + dSijU(I,3,1)
                     SijU(J,3,2) = SijU(J,3,2) + dSijU(I,3,2)
                     SijU(J,3,3) = SijU(J,3,3) + dSijU(I,3,3)

                     SijW(J,1,1) = SijW(J,1,1) + dSijW(I,1,1)
                     SijW(J,1,2) = SijW(J,1,2) + dSijW(I,1,2)
                     SijW(J,1,3) = SijW(J,1,3) + dSijW(I,1,3)
                     SijW(J,2,1) = SijW(J,2,1) + dSijW(I,2,1)
                     SijW(J,2,2) = SijW(J,2,2) + dSijW(I,2,2)
                     SijW(J,2,3) = SijW(J,2,3) + dSijW(I,2,3)
                     SijW(J,3,1) = SijW(J,3,1) + dSijW(I,3,1)
                     SijW(J,3,2) = SijW(J,3,2) + dSijW(I,3,2)
                     SijW(J,3,3) = SijW(J,3,3) + dSijW(I,3,3)		  
 30               CONTINUE

                  ECHI=0 
                  EKAP=0
                  DO 70 I=1,NBIN
                     ECHI=ECHI+(DEL**3.)*(CHI/V)*PHIA(I)*PHIB(I)
                     EKAP=EKAP+(DEL**3.)*(KAP/V)*(PHIA(I)+PHIB(I)-1.)**2.
 70               CONTINUE
               else
                  DO 35 I=1,NBIN
                     PHIA(I)=0
                     PHIB(I)=0
 35               CONTINUE
               endif

               SUCCESS(MCTYPE)=SUCCESS(MCTYPE)+1
            endif

!     Adapt the amplitude of step every NADAPT steps
            !restart WINDOW every 20 NADAPT
            if (mod(ISTEP,NADAPT(MCTYPE)*20).EQ.0) then
               WINDOW(MCTYPE) = NINT(N*grnd())
            endif

            !amplitude and window adaptations
            if (mod(ISTEP,NADAPT(MCTYPE)).EQ.0) then
               PHIT(MCTYPE)=real(SUCCESS(MCTYPE))/real(NADAPT(MCTYPE))

               if (PHIT(MCTYPE).GT.PDESIRE(MCTYPE)) then
                  MCAMP(MCTYPE)=MCAMP(MCTYPE)*1.05
               else
                  MCAMP(MCTYPE)=MCAMP(MCTYPE)*0.95
               endif
               
               if (MCAMP(MCTYPE).GT.MAXAMP(MCTYPE)) then
                  MCAMP(MCTYPE)=MAXAMP(MCTYPE)
                  WINDOW(MCTYPE)=WINDOW(MCTYPE)+1
               endif
			   
               if (MCAMP(MCTYPE).LT.MINAMP(MCTYPE)) then
                  MCAMP(MCTYPE)=MINAMP(MCTYPE)
                  WINDOW(MCTYPE)=WINDOW(MCTYPE)-1                
                  if ((MCTYPE.EQ.4).OR.(MCTYPE.EQ.5).OR.(MCTYPE.EQ.6).OR.(MCTYPE.EQ.7)) then
                     MOVEON(MCTYPE)=0
                  endif
               endif
               
               !window limits
               if (WINDOW(MCTYPE).LT.0) then
                  WINDOW(MCTYPE)=1
               endif

               if (WINDOW(MCTYPE).GT.N) then
                  WINDOW(MCTYPE)=N
               endif
		   
               SUCCESS(MCTYPE)=0
               IB=1
               DO 40 I=1,NP
                  R0(1)=nint(R(IB,1)/LBOX-0.5)*LBOX
                  R0(2)=nint(R(IB,2)/LBOX-0.5)*LBOX
                  R0(3)=nint(R(IB,3)/LBOX-0.5)*LBOX
                  DO 50 J=1,N
                     R(IB,1)=R(IB,1)-R0(1)
                     R(IB,2)=R(IB,2)-R0(2)
                     R(IB,3)=R(IB,3)-R0(3)
                     IB=IB+1
 50               CONTINUE
 40             CONTINUE

            endif

 60         CONTINUE

 10      CONTINUE

         IF ((PTON.EQ.1).AND.(INTON.EQ.1)) THEN
            INDPT=ISTEP+NSTEP*(IND-1)
!            call energy_int(R,AB,NP,NB,NT,NBIN,V,CHI,KAP,LBOX,DEL,ECHI,EKAP)
!            PRINT*, ' IN MC Code ', ECHI
            call ptsched(INDPT,IND,NRABOVE,NRBELOW,NRMAX,NRNOW,CHI,ECHI,NPT,PTID)
         ENDIF
         
         if (SON.EQ.1) then
            if (mod(ISTEP,floor(REAL(NSTEP)/NVEC)).EQ.0) then
               call gcalc(KMAG,SK,R,AB,NT)
               KVEC(1:KNUM)=KMAG(1:KNUM)
               SVEC(1:KNUM)=SVEC(1:KNUM)+SK(1:KNUM)
            endif
         endif

         ISTEP=ISTEP+1
      ENDDO


      
      DO K=1,KNUM
         SVEC(K)=REAL(SVEC(K)) / REAL(NVEC)
      ENDDO

      RETURN      
      END
      
!---------------------------------------------------------------*
