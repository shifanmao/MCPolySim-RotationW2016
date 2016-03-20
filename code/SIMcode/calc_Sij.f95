
!---------------------------------------------------------------*
      
!     Calculate the Nematic Order Parameters Sij(u) and Sij(v)
!     Kevin Hou 2-23-2016
!
!     I basically just made minor modifications to r_to_phi, check that one
      
      SUBROUTINE calc_Sij(R,U,W,NT,N,NP,NBIN, &
      V,LBOX,DEL,SijU,SijW,AU,AUV,AV,EALIGN,SijE)
      
      PARAMETER (PI=3.141592654) ! Value of pi
      
      DOUBLE PRECISION R(NT,3)	! Conformation of polymer chains
      DOUBLE PRECISION U(NT,3)   ! Tangent Vectors
      DOUBLE PRECISION W(NT,3)   ! Twist (orthogonal) Vectors
      INTEGER IND(NT,2)         ! Start and end integers      
      INTEGER NT                ! Number of beads
      INTEGER NP                ! Total number of polymers
      INTEGER N                 ! Number of beads per polymer
      
!     Variables for density calculation
      
      DOUBLE PRECISION SijU(NBIN,3,3) ! Order Parameter Sij(u)
      DOUBLE PRECISION SijW(NBIN,3,3) ! Order Parameter Sij(v)
      DOUBLE PRECISION SijE(NBIN)     ! Energy per bin  
      DOUBLE PRECISION EALIGN         ! Alignment Energy
      
!     Simulation input variables
      DOUBLE PRECISION AU,AUV,AV       ! Parameters Describing Strength of interaction
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
      
!     Initialize everything to zero
      do 10 I=1,NBIN
         SijE(I) = 0.
         SijU(I,1,1)=0.
         SijU(I,1,2)=0.
         SijU(I,1,3)=0.
         SijU(I,2,1)=0.
         SijU(I,2,2)=0.
         SijU(I,2,3)=0.
         SijU(I,3,1)=0.
         SijU(I,3,2)=0.
         SijU(I,3,3)=0.
         SijW(I,1,1)=0.
         SijW(I,1,2)=0.
         SijW(I,1,3)=0.
         SijW(I,2,1)=0.
         SijW(I,2,2)=0.
         SijW(I,2,3)=0.
         SijW(I,3,1)=0.
         SijW(I,3,2)=0.
         SijW(I,3,3)=0.
 10   continue
      NBINX=nint(LBOX/DEL)
      
!     Cycle through the beads
      
      IB=1
      do 20 I=1,NP
         do 30 J=1,N
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
            
            do 40 ISX=1,2
               do 50 ISY=1,2
                  do 60 ISZ=1,2
                     WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                     INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX+(IZ(ISZ)-1)*NBINX**2

                     SijU(INDBIN,1,1) = SijU(INDBIN,1,1) + (U(IB,1)*U(IB,1)-(1/3))*(WTOT*V/DEL**3.)
                     SijU(INDBIN,1,2) = SijU(INDBIN,1,2) + (U(IB,1)*U(IB,2))*(WTOT*V/DEL**3.)
                     SijU(INDBIN,1,3) = SijU(INDBIN,1,3) + (U(IB,1)*U(IB,3))*(WTOT*V/DEL**3.)
                     SijU(INDBIN,2,1) = SijU(INDBIN,2,1) + (U(IB,2)*U(IB,1))*(WTOT*V/DEL**3.)
                     SijU(INDBIN,2,2) = SijU(INDBIN,2,2) + (U(IB,2)*U(IB,2)-(1/3))*(WTOT*V/DEL**3.)
                     SijU(INDBIN,2,3) = SijU(INDBIN,2,3) + (U(IB,2)*U(IB,3))*(WTOT*V/DEL**3.)
                     SijU(INDBIN,3,1) = SijU(INDBIN,3,1) + (U(IB,3)*U(IB,1))*(WTOT*V/DEL**3.)
                     SijU(INDBIN,3,2) = SijU(INDBIN,3,2) + (U(IB,3)*U(IB,2))*(WTOT*V/DEL**3.)
                     SijU(INDBIN,3,3) = SijU(INDBIN,3,3) + (U(IB,3)*U(IB,3)-(1/3))*(WTOT*V/DEL**3.)

                     SijW(INDBIN,1,1) = SijW(INDBIN,1,1) + (W(IB,1)*W(IB,1)-(1/3))*(WTOT*V/DEL**3.)
                     SijW(INDBIN,1,2) = SijW(INDBIN,1,2) + (W(IB,1)*W(IB,2))*(WTOT*V/DEL**3.)
                     SijW(INDBIN,1,3) = SijW(INDBIN,1,3) + (W(IB,1)*W(IB,3))*(WTOT*V/DEL**3.)
                     SijW(INDBIN,2,1) = SijW(INDBIN,2,1) + (W(IB,2)*W(IB,1))*(WTOT*V/DEL**3.)
                     SijW(INDBIN,2,2) = SijW(INDBIN,2,2) + (W(IB,2)*W(IB,2)-(1/3))*(WTOT*V/DEL**3.)
                     SijW(INDBIN,2,3) = SijW(INDBIN,2,3) + (W(IB,2)*W(IB,3))*(WTOT*V/DEL**3.)
                     SijW(INDBIN,3,1) = SijW(INDBIN,3,1) + (W(IB,3)*W(IB,1))*(WTOT*V/DEL**3.)
                     SijW(INDBIN,3,2) = SijW(INDBIN,3,2) + (W(IB,3)*W(IB,2))*(WTOT*V/DEL**3.)
                     SijW(INDBIN,3,3) = SijW(INDBIN,3,3) + (W(IB,3)*W(IB,3)-(1/3))*(WTOT*V/DEL**3.)
 60               continue
 50            continue
 40         continue
            IB=IB+1
            
 30      continue
 20   continue
      
      EALIGN = 0
      ! Now, Calculate Energies
      do 70 I=1,NBIN
        do 80 X=1,3
            do 90 Y =1,3
                SijE(I) = SijE(I) + SijU(I,X,Y)**2
                EALIGN = EALIGN - (DEL**3./V)*(AU/2.)*(SijU(I,X,Y)**2) 
                !+ AUV*(SijU(I,X,Y)*SijW(I,X,Y)) - (AV/2)*(SijW(I,X,Y)**2)
 90         continue
 80     continue
 70   continue


      RETURN            
      END
      
!---------------------------------------------------------------*
