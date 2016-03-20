
!---------------------------------------------------------------*
      
!     Calculate the Nematic Order Parameters Sij(u) and Sij(v), and then the change
!     in alignment energy
!     Kevin Hou 2-23-2016
!
!     I basically just made minor modifications to MC_int check that one
!     Discretization code from Shifan Mao
      
      SUBROUTINE MC_align(DEALIGN,R,U,W,RP,UP,WP,NT,N,NP,NTOT,NBIN, &
      V,LBOX,DEL,AU,AUV,AV,SijU,SijW,dSijU,dSijW,IT1,IT2,INDSij)
      
      PARAMETER (PI=3.141592654) ! Value of pi
      
      DOUBLE PRECISION R(NT,3)  ! Conformation of polymer chains
      DOUBLE PRECISION U(NT,3)   ! Tangent Vectors
      DOUBLE PRECISION W(NT,3)   ! Twist (orthogonal) Vectors
      DOUBLE PRECISION RP(NT,3)  ! Conformation of polymer chains
      DOUBLE PRECISION UP(NT,3)   ! Tangent Vectors
      DOUBLE PRECISION WP(NT,3)   ! Twist (orthogonal) Vectors
      DOUBLE PRECISION DR(NT,3)   ! Change in Bead Position
      INTEGER IND(NT,2)         ! Start and end integers      
      INTEGER NT                ! Number of beads
      INTEGER NP                ! Total number of polymers
      INTEGER N                 ! Number of beads per polymer
      INTEGER IT1               ! Index of test bead 1
      INTEGER IT2               ! Index of test bead 2
      INTEGER INDSij(NBIN)      ! Indices
      INTEGER NSij              ! Number of S values that change
      
!     ALIGNMENT VARIABLES
      
      DOUBLE PRECISION DEALIGN
      DOUBLE PRECISION EALIGNOLD
      DOUBLE PRECISION EALIGNNEW
      DOUBLE PRECISION AU,AUV,AV       ! Parameters Describing Strength of interaction
      DOUBLE PRECISION SijU(NBIN,3,3)  ! Order Parameter Sij(u) before move
      DOUBLE PRECISION SijW(NBIN,3,3)  ! Order Parameter Sij(v) before move
      DOUBLE PRECISION dSijU(NBIN,3,3)
      DOUBLE PRECISION dSijW(NBIN,3,3)
      
!     Simulation input variables
      
      DOUBLE PRECISION V       ! Volume of monomer A
      DOUBLE PRECISION LBOX     ! Simulation box size (approximate)
      DOUBLE PRECISION DEL      ! Discretization size (approximate)     
      INTEGER I,J,X,Y

      INTEGER IX(2),IY(2),IZ(2)
      INTEGER IB
      INTEGER NBINX
      
      DOUBLE PRECISION WX(2),WY(2),WZ(2)
      DOUBLE PRECISION WTOT
      DOUBLE PRECISION RBIN(3)
      INTEGER INDBIN
      INTEGER ISX,ISY,ISZ

    !Initialize Change Terms
      do 10 I=1,NBIN
         dSijU(I,1,1)=0.
         dSijU(I,1,2)=0.
         dSijU(I,1,3)=0.
         dSijU(I,2,1)=0.
         dSijU(I,2,2)=0.
         dSijU(I,2,3)=0.
         dSijU(I,3,1)=0.
         dSijU(I,3,2)=0.
         dSijU(I,3,3)=0.
         dSijW(I,1,1)=0.
         dSijW(I,1,2)=0.
         dSijW(I,1,3)=0.
         dSijW(I,2,1)=0.
         dSijW(I,2,2)=0.
         dSijW(I,2,3)=0.
         dSijW(I,3,1)=0.
         dSijW(I,3,2)=0.
         dSijW(I,3,3)=0.
         INDSij(I)=0.
 10   continue
      NBINX=nint(LBOX/DEL)

      DO 15 I=IT1,IT2
         DR(I,1)=RP(I,1)-R(I,1)
         DR(I,2)=RP(I,2)-R(I,2)
         DR(I,3)=RP(I,3)-R(I,3)
 15   CONTINUE
      
      NSij = 0

      ! Set Up Bins
      do 20 IB=IT1,IT2
         RBIN(1)=R(IB,1)-nint(R(IB,1)/LBOX-0.5)*LBOX
         RBIN(2)=R(IB,2)-nint(R(IB,2)/LBOX-0.5)*LBOX
         RBIN(3)=R(IB,3)-nint(R(IB,3)/LBOX-0.5)*LBOX
         
         IX(1)=nint(RBIN(1)/DEL+0.5)
         IY(1)=nint(RBIN(2)/DEL+0.5)
         IZ(1)=nint(RBIN(3)/DEL+0.5)
         
         IX(2)=IX(1)-1
         IY(2)=IY(1)-1
         IZ(2)=IZ(1)-1

!     Calculate the bin weighting
         
         WX(2)=(RBIN(1)-IX(1)*DEL)/(IX(2)*DEL-IX(1)*DEL)
         WX(1)=(IX(2)*DEL-RBIN(1))/(IX(2)*DEL-IX(1)*DEL)
         WY(2)=(RBIN(2)-IY(1)*DEL)/(IY(2)*DEL-IY(1)*DEL)
         WY(1)=(IY(2)*DEL-RBIN(2))/(IY(2)*DEL-IY(1)*DEL)
         WZ(2)=(RBIN(3)-IZ(1)*DEL)/(IZ(2)*DEL-IZ(1)*DEL)
         WZ(1)=(IZ(2)*DEL-RBIN(3))/(IZ(2)*DEL-IZ(1)*DEL)
         
         IX(1)=IX(1)-floor(REAL((IX(1)-1))/REAL(NBINX)) * NBINX
         IX(2)=IX(2)-floor(REAL((IX(2)-1))/REAL(NBINX)) * NBINX
         IY(1)=IY(1)-floor(REAL((IY(1)-1))/REAL(NBINX)) * NBINX
         IY(2)=IY(2)-floor(REAL((IY(2)-1))/REAL(NBINX)) * NBINX
         IZ(1)=IZ(1)-floor(REAL((IZ(1)-1))/REAL(NBINX)) * NBINX
         IZ(2)=IZ(2)-floor(REAL((IZ(2)-1))/REAL(NBINX)) * NBINX

!     Add volume fraction with weighting to each bin
         
         do 30 ISX=1,2
            do 40 ISY=1,2
               do 50 ISZ=1,2
                  WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                  INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX+(IZ(ISZ)-1)*NBINX**2
                  
                  I=1
                  STOP=0
                  do while (STOP.EQ.0)
                     if (INDBIN.EQ.INDSij(I).AND.I.LE.NSij) then
                        STOP=1
                        dSijU(I,1,1) = dSijU(I,1,1) - (U(IB,1)*U(IB,1)-(1/3))*(WTOT*V/DEL**3.)
                        dSijU(I,1,2) = dSijU(I,1,2) - (U(IB,1)*U(IB,2))*(WTOT*V/DEL**3.)
                        dSijU(I,1,3) = dSijU(I,1,3) - (U(IB,1)*U(IB,3))*(WTOT*V/DEL**3.)
                        dSijU(I,2,1) = dSijU(I,2,1) - (U(IB,2)*U(IB,1))*(WTOT*V/DEL**3.)
                        dSijU(I,2,2) = dSijU(I,2,2) - (U(IB,2)*U(IB,2)-(1/3))*(WTOT*V/DEL**3.)
                        dSijU(I,2,3) = dSijU(I,2,3) - (U(IB,2)*U(IB,3))*(WTOT*V/DEL**3.)
                        dSijU(I,3,1) = dSijU(I,3,1) - (U(IB,3)*U(IB,1))*(WTOT*V/DEL**3.)
                        dSijU(I,3,2) = dSijU(I,3,2) - (U(IB,3)*U(IB,2))*(WTOT*V/DEL**3.)
                        dSijU(I,3,3) = dSijU(I,3,3) - (U(IB,3)*U(IB,3)-(1/3))*(WTOT*V/DEL**3.)
 
                        dSijW(I,1,1) = dSijW(I,1,1) - (W(IB,1)*W(IB,1)-(1/3))*(WTOT*V/DEL**3.)
                        dSijW(I,1,2) = dSijW(I,1,2) - (W(IB,1)*W(IB,2))*(WTOT*V/DEL**3.)
                        dSijW(I,1,3) = dSijW(I,1,3) - (W(IB,1)*W(IB,3))*(WTOT*V/DEL**3.)
                        dSijW(I,2,1) = dSijW(I,2,1) - (W(IB,2)*W(IB,1))*(WTOT*V/DEL**3.)
                        dSijW(I,2,2) = dSijW(I,2,2) - (W(IB,2)*W(IB,2)-(1/3))*(WTOT*V/DEL**3.)
                        dSijW(I,2,3) = dSijW(I,2,3) - (W(IB,2)*W(IB,3))*(WTOT*V/DEL**3.)
                        dSijW(I,3,1) = dSijW(I,3,1) - (W(IB,3)*W(IB,1))*(WTOT*V/DEL**3.)
                        dSijW(I,3,2) = dSijW(I,3,2) - (W(IB,3)*W(IB,2))*(WTOT*V/DEL**3.)
                        dSijW(I,3,3) = dSijW(I,3,3) - (W(IB,3)*W(IB,3)-(1/3))*(WTOT*V/DEL**3.)

                     elseif (I.GT.NSij) then
                        STOP=1
                        NSij=NSij+1
                        INDSij(I)=INDBIN

                        dSijU(I,1,1) = dSijU(I,1,1) - (U(IB,1)*U(IB,1)-(1/3))*(WTOT*V/DEL**3.)
                        dSijU(I,1,2) = dSijU(I,1,2) - (U(IB,1)*U(IB,2))*(WTOT*V/DEL**3.)
                        dSijU(I,1,3) = dSijU(I,1,3) - (U(IB,1)*U(IB,3))*(WTOT*V/DEL**3.)
                        dSijU(I,2,1) = dSijU(I,2,1) - (U(IB,2)*U(IB,1))*(WTOT*V/DEL**3.)
                        dSijU(I,2,2) = dSijU(I,2,2) - (U(IB,2)*U(IB,2)-(1/3))*(WTOT*V/DEL**3.)
                        dSijU(I,2,3) = dSijU(I,2,3) - (U(IB,2)*U(IB,3))*(WTOT*V/DEL**3.)
                        dSijU(I,3,1) = dSijU(I,3,1) - (U(IB,3)*U(IB,1))*(WTOT*V/DEL**3.)
                        dSijU(I,3,2) = dSijU(I,3,2) - (U(IB,3)*U(IB,2))*(WTOT*V/DEL**3.)
                        dSijU(I,3,3) = dSijU(I,3,3) - (U(IB,3)*U(IB,3)-(1/3))*(WTOT*V/DEL**3.)
 
                        dSijW(I,1,1) = dSijW(I,1,1) - (W(IB,1)*W(IB,1)-(1/3))*(WTOT*V/DEL**3.)
                        dSijW(I,1,2) = dSijW(I,1,2) - (W(IB,1)*W(IB,2))*(WTOT*V/DEL**3.)
                        dSijW(I,1,3) = dSijW(I,1,3) - (W(IB,1)*W(IB,3))*(WTOT*V/DEL**3.)
                        dSijW(I,2,1) = dSijW(I,2,1) - (W(IB,2)*W(IB,1))*(WTOT*V/DEL**3.)
                        dSijW(I,2,2) = dSijW(I,2,2) - (W(IB,2)*W(IB,2)-(1/3))*(WTOT*V/DEL**3.)
                        dSijW(I,2,3) = dSijW(I,2,3) - (W(IB,2)*W(IB,3))*(WTOT*V/DEL**3.)
                        dSijW(I,3,1) = dSijW(I,3,1) - (W(IB,3)*W(IB,1))*(WTOT*V/DEL**3.)
                        dSijW(I,3,2) = dSijW(I,3,2) - (W(IB,3)*W(IB,2))*(WTOT*V/DEL**3.)
                        dSijW(I,3,3) = dSijW(I,3,3) - (W(IB,3)*W(IB,3)-(1/3))*(WTOT*V/DEL**3.)
                     else
                        I=I+1
                     endif    
                  enddo

 50            continue
 40         continue
 30      continue
 20   continue

      do 60 IB=IT1,IT2
         RBIN(1)=R(IB,1)+DR(IB,1)
         RBIN(2)=R(IB,2)+DR(IB,2)
         RBIN(3)=R(IB,3)+DR(IB,3)

         RBIN(1)=RBIN(1)-nint(RBIN(1)/LBOX-0.5)*LBOX
         RBIN(2)=RBIN(2)-nint(RBIN(2)/LBOX-0.5)*LBOX
         RBIN(3)=RBIN(3)-nint(RBIN(3)/LBOX-0.5)*LBOX
         
         IX(1)=nint(RBIN(1)/DEL+0.5)
         IY(1)=nint(RBIN(2)/DEL+0.5)
         IZ(1)=nint(RBIN(3)/DEL+0.5)
         
         IX(2)=IX(1)-1
         IY(2)=IY(1)-1
         IZ(2)=IZ(1)-1
         
!     Calculate the bin weighting
         
         WX(2)=(RBIN(1)-IX(1)*DEL)/(IX(2)*DEL-IX(1)*DEL)
         WX(1)=(IX(2)*DEL-RBIN(1))/(IX(2)*DEL-IX(1)*DEL)
         WY(2)=(RBIN(2)-IY(1)*DEL)/(IY(2)*DEL-IY(1)*DEL)
         WY(1)=(IY(2)*DEL-RBIN(2))/(IY(2)*DEL-IY(1)*DEL)
         WZ(2)=(RBIN(3)-IZ(1)*DEL)/(IZ(2)*DEL-IZ(1)*DEL)
         WZ(1)=(IZ(2)*DEL-RBIN(3))/(IZ(2)*DEL-IZ(1)*DEL)
         
         IX(1)=IX(1)-floor(REAL((IX(1)-1))/REAL(NBINX)) * NBINX
         IX(2)=IX(2)-floor(REAL((IX(2)-1))/REAL(NBINX)) * NBINX
         IY(1)=IY(1)-floor(REAL((IY(1)-1))/REAL(NBINX)) * NBINX
         IY(2)=IY(2)-floor(REAL((IY(2)-1))/REAL(NBINX)) * NBINX
         IZ(1)=IZ(1)-floor(REAL((IZ(1)-1))/REAL(NBINX)) * NBINX
         IZ(2)=IZ(2)-floor(REAL((IZ(2)-1))/REAL(NBINX)) * NBINX


!     Add volume fraction with weighting to each bin
         
         do 70 ISX=1,2
            do 80 ISY=1,2
               do 90 ISZ=1,2
                  WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                  INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX+(IZ(ISZ)-1)*NBINX**2
                  
                  I=1
                  STOP=0
                  do while (STOP.EQ.0)
                     if (INDBIN.EQ.INDSij(I).AND.I.LE.NSij) then
                        STOP=1
                        dSijU(I,1,1) = dSijU(I,1,1) + (U(IB,1)*U(IB,1)-(1/3))*(WTOT*V/DEL**3.)
                        dSijU(I,1,2) = dSijU(I,1,2) + (U(IB,1)*U(IB,2))*(WTOT*V/DEL**3.)
                        dSijU(I,1,3) = dSijU(I,1,3) + (U(IB,1)*U(IB,3))*(WTOT*V/DEL**3.)
                        dSijU(I,2,1) = dSijU(I,2,1) + (U(IB,2)*U(IB,1))*(WTOT*V/DEL**3.)
                        dSijU(I,2,2) = dSijU(I,2,2) + (U(IB,2)*U(IB,2)-(1/3))*(WTOT*V/DEL**3.)
                        dSijU(I,2,3) = dSijU(I,2,3) + (U(IB,2)*U(IB,3))*(WTOT*V/DEL**3.)
                        dSijU(I,3,1) = dSijU(I,3,1) + (U(IB,3)*U(IB,1))*(WTOT*V/DEL**3.)
                        dSijU(I,3,2) = dSijU(I,3,2) + (U(IB,3)*U(IB,2))*(WTOT*V/DEL**3.)
                        dSijU(I,3,3) = dSijU(I,3,3) + (U(IB,3)*U(IB,3)-(1/3))*(WTOT*V/DEL**3.)
 
                        dSijW(I,1,1) = dSijW(I,1,1) + (W(IB,1)*W(IB,1)-(1/3))*(WTOT*V/DEL**3.)
                        dSijW(I,1,2) = dSijW(I,1,2) + (W(IB,1)*W(IB,2))*(WTOT*V/DEL**3.)
                        dSijW(I,1,3) = dSijW(I,1,3) + (W(IB,1)*W(IB,3))*(WTOT*V/DEL**3.)
                        dSijW(I,2,1) = dSijW(I,2,1) + (W(IB,2)*W(IB,1))*(WTOT*V/DEL**3.)
                        dSijW(I,2,2) = dSijW(I,2,2) + (W(IB,2)*W(IB,2)-(1/3))*(WTOT*V/DEL**3.)
                        dSijW(I,2,3) = dSijW(I,2,3) + (W(IB,2)*W(IB,3))*(WTOT*V/DEL**3.)
                        dSijW(I,3,1) = dSijW(I,3,1) + (W(IB,3)*W(IB,1))*(WTOT*V/DEL**3.)
                        dSijW(I,3,2) = dSijW(I,3,2) + (W(IB,3)*W(IB,2))*(WTOT*V/DEL**3.)
                        dSijW(I,3,3) = dSijW(I,3,3) + (W(IB,3)*W(IB,3)-(1/3))*(WTOT*V/DEL**3.)

                     elseif (I.GT.NSij) then
                        STOP=1
                        NSij=NSij+1
                        INDSij(I)=INDBIN

                        dSijU(I,1,1) = dSijU(I,1,1) + (UP(IB,1)*UP(IB,1)-(1/3))*(WTOT*V/DEL**3.)
                        dSijU(I,1,2) = dSijU(I,1,2) + (UP(IB,1)*UP(IB,2))*(WTOT*V/DEL**3.)
                        dSijU(I,1,3) = dSijU(I,1,3) + (UP(IB,1)*UP(IB,3))*(WTOT*V/DEL**3.)
                        dSijU(I,2,1) = dSijU(I,2,1) + (UP(IB,2)*UP(IB,1))*(WTOT*V/DEL**3.)
                        dSijU(I,2,2) = dSijU(I,2,2) + (UP(IB,2)*UP(IB,2)-(1/3))*(WTOT*V/DEL**3.)
                        dSijU(I,2,3) = dSijU(I,2,3) + (UP(IB,2)*UP(IB,3))*(WTOT*V/DEL**3.)
                        dSijU(I,3,1) = dSijU(I,3,1) + (UP(IB,3)*UP(IB,1))*(WTOT*V/DEL**3.)
                        dSijU(I,3,2) = dSijU(I,3,2) + (UP(IB,3)*UP(IB,2))*(WTOT*V/DEL**3.)
                        dSijU(I,3,3) = dSijU(I,3,3) + (UP(IB,3)*UP(IB,3)-(1/3))*(WTOT*V/DEL**3.)
 
                        dSijW(I,1,1) = dSijW(I,1,1) + (WP(IB,1)*WP(IB,1)-(1/3))*(WTOT*V/DEL**3.)
                        dSijW(I,1,2) = dSijW(I,1,2) + (WP(IB,1)*WP(IB,2))*(WTOT*V/DEL**3.)
                        dSijW(I,1,3) = dSijW(I,1,3) + (WP(IB,1)*WP(IB,3))*(WTOT*V/DEL**3.)
                        dSijW(I,2,1) = dSijW(I,2,1) + (WP(IB,2)*WP(IB,1))*(WTOT*V/DEL**3.)
                        dSijW(I,2,2) = dSijW(I,2,2) + (WP(IB,2)*WP(IB,2)-(1/3))*(WTOT*V/DEL**3.)
                        dSijW(I,2,3) = dSijW(I,2,3) + (WP(IB,2)*WP(IB,3))*(WTOT*V/DEL**3.)
                        dSijW(I,3,1) = dSijW(I,3,1) + (WP(IB,3)*WP(IB,1))*(WTOT*V/DEL**3.)
                        dSijW(I,3,2) = dSijW(I,3,2) + (WP(IB,3)*WP(IB,2))*(WTOT*V/DEL**3.)
                        dSijW(I,3,3) = dSijW(I,3,3) + (WP(IB,3)*WP(IB,3)-(1/3))*(WTOT*V/DEL**3.)
                     else
                        I=I+1
                     endif    
                  enddo

 90            continue
 80         continue
 70      continue
 60   continue
      
      DEALIGN = 0.
      EALIGNOLD = 0.
      EALIGNNEW = 0.

      ! Now, Calculate Energies
      do 100 I=1,NSij
        J = INDSij(I)
        do 110 X=1,3
            do 120 Y =1,3
                EALIGNOLD = EALIGNOLD-(DEL**3./V)*((AU/2)*(SijU(J,X,Y)**2) + &
                AUV*(SijU(J,X,Y)*SijW(J,X,Y))-(AV/2)*(SijW(J,X,Y)**2))

                EALIGNNEW = EALIGNNEW-(DEL**3./V)*((AU/2)*((SijU(J,X,Y)+dSijU(I,X,Y))**2) + &
                AUV*((SijU(J,X,Y)+dSijU(I,X,Y))*(SijW(J,X,Y)+dSijW(I,X,Y))) - &
                AV*0.5*((SijW(J,X,Y)-dSijW(I,X,Y))**2))
 120        continue
 110    continue
 100  continue

      DEALIGN = EALIGNNEW - EALIGNOLD

      RETURN            
      END
      
!---------------------------------------------------------------*
