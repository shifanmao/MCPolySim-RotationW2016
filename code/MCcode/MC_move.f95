!---------------------------------------------------------------*

! Find change in bead position for a crank-shaft type move
     
      SUBROUTINE MC_move(R,U,W,RP,UP,WP,NT,N,NP,IP,IB1,IB2,IT1,IT2,MCTYPE,MCAMP,WINDOW)

      use mt19937, only : grnd, sgrnd, rnorm, mt, mti

      PARAMETER (PI=3.141592654) ! Value of pi
      
      DOUBLE PRECISION R(NT,3)   ! Bead positions
      DOUBLE PRECISION U(NT,3)   ! Tangent vectors
      DOUBLE PRECISION W(NT,3)   ! Twist vectors (KHOU)
      DOUBLE PRECISION RP(NT,3)  ! Bead positions
      DOUBLE PRECISION UP(NT,3)  ! Tangent vectors
      DOUBLE PRECISION WP(NT,3)  ! Twist vectors (KHOU)

      INTEGER N,NP,NT           ! Number of beads

      INTEGER IP                ! Test polymer 
      INTEGER IB1               ! Test bead position 1
      INTEGER IT1               ! Index of test bead 1
      INTEGER IB2               ! Test bead position 2
      INTEGER IT2               ! Index of test bead 2

      INTEGER I,J			! Test indices
      
! Variables for the crank-shaft move

      DOUBLE PRECISION TA(3)    ! Axis of rotation
      DOUBLE PRECISION P1(3)    ! Point on rotation line
      DOUBLE PRECISION MAG      ! Magnitude of vector
      DOUBLE PRECISION ROT(4,4) ! Rotation matrix
            
      DOUBLE PRECISION ALPHA    ! Angle of move
      DOUBLE PRECISION BETA     ! Angle of move
      
!     MC adaptation variables (KHOU)
      
      DOUBLE PRECISION MCAMP(7) ! Amplitude of random change      
      INTEGER MCTYPE            ! Type of MC move
      DOUBLE PRECISION DR(3)    ! Displacement for slide move
      INTEGER TEMP
      INTEGER WINDOW(7)			! Size of window for bead selection

!     Perform crank-shaft move (MCTYPE 1)
          
      if (MCTYPE.EQ.1) then
         IP=nint(0.5+grnd()*NP)
         IB1=nint(0.5+grnd()*N)
         IB2=IB1+nint((grnd()-0.5)*(2*WINDOW(MCTYPE)+1))

         if (IB2.LT.1) then
            IB2=1
         endif
         if (IB2.GT.N) then
            IB2=N
         endif

         if (IB2.LT.IB1) then
            TEMP=IB1
            IB1=IB2
            IB2=TEMP
         endif
         IT1=N*(IP-1)+IB1
         IT2=N*(IP-1)+IB2
         
         if (IB1.EQ.IB2.AND.IB1.EQ.1) then
            TA(1)=R(IT1+1,1)-R(IT1,1)
            TA(2)=R(IT1+1,2)-R(IT1,2)
            TA(3)=R(IT1+1,3)-R(IT1,3)
         elseif (IB1.EQ.IB2.AND.IB1.EQ.N) then
            TA(1)=R(IT1,1)-R(IT1-1,1)
            TA(2)=R(IT1,2)-R(IT1-1,2)
            TA(3)=R(IT1,3)-R(IT1-1,3)
         elseif (IB1.EQ.IB2.AND.IB1.NE.1.AND.IB2.NE.N) then
            TA(1)=R(IT1+1,1)-R(IT1-1,1)
            TA(2)=R(IT1+1,2)-R(IT1-1,2)
            TA(3)=R(IT1+1,3)-R(IT1-1,3)
         else
            TA(1)=R(IT2,1)-R(IT1,1)
            TA(2)=R(IT2,2)-R(IT1,2)
            TA(3)=R(IT2,3)-R(IT1,3)
         endif
         MAG=sqrt(TA(1)**2.+TA(2)**2.+TA(3)**2.)
         TA(1)=TA(1)/MAG
         TA(2)=TA(2)/MAG
         TA(3)=TA(3)/MAG
         P1(1)=R(IT1,1)
         P1(2)=R(IT1,2)
         P1(3)=R(IT1,3)	  
         
         ALPHA=MCAMP(1)*(grnd()-0.5)
         
         ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
         ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
         ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
		 -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)
         
         ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
         ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
         ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
		 -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)
             
         ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
         ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
         ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
         ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
		 -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)
             
         DO 10 I=IT1,IT2
            RP(I,1)=ROT(1,4)+ROT(1,1)*R(I,1)+ROT(1,2)*R(I,2)+ROT(1,3)*R(I,3)
            RP(I,2)=ROT(2,4)+ROT(2,1)*R(I,1)+ROT(2,2)*R(I,2)+ROT(2,3)*R(I,3)
            RP(I,3)=ROT(3,4)+ROT(3,1)*R(I,1)+ROT(3,2)*R(I,2)+ROT(3,3)*R(I,3)
            UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
            UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
            UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)

            !KHOU
            WP(I,1)=ROT(1,1)*W(I,1)+ROT(1,2)*W(I,2)+ROT(1,3)*W(I,3)
            WP(I,2)=ROT(2,1)*W(I,1)+ROT(2,2)*W(I,2)+ROT(2,3)*W(I,3)
            WP(I,3)=ROT(3,1)*W(I,1)+ROT(3,2)*W(I,2)+ROT(3,3)*W(I,3)
 10      CONTINUE

!     Perform slide move (MCTYPE 2)
         
      elseif (MCTYPE.EQ.2) then
         
         IP=nint(0.5+grnd()*NP)
         IB1=nint(0.5+grnd()*N)
         IB2=IB1+nint((grnd()-0.5)*(2*WINDOW(MCTYPE)+1))

         if (IB2.LT.1) then
            IB2=1
         endif
         if (IB2.GT.N) then
            IB2=N
         endif

         if (IB2.LT.IB1) then
            TEMP=IB1
            IB1=IB2
            IB2=TEMP
         endif
         IT1=N*(IP-1)+IB1
         IT2=N*(IP-1)+IB2
         
         DR(1)=MCAMP(2)*(grnd()-0.5)
         DR(2)=MCAMP(2)*(grnd()-0.5)
         DR(3)=MCAMP(2)*(grnd()-0.5)
         
         DO 20 I=IT1,IT2
            RP(I,1)=R(I,1)+DR(1)
            RP(I,2)=R(I,2)+DR(2)
            RP(I,3)=R(I,3)+DR(3)
            UP(I,1)=U(I,1)
            UP(I,2)=U(I,2)
            UP(I,3)=U(I,3)

            !KHOU
            WP(I,1)=W(I,1)
            WP(I,2)=W(I,2)
            WP(I,3)=W(I,3)
 20      CONTINUE    
             
!     Perform pivot move (MCTYPE 3)
             
      elseif (MCTYPE.EQ.3) then

         IP=nint(0.5+grnd()*NP)
         IB1=nint(0.5+grnd()*(2*WINDOW(MCTYPE)))

         if (IB1.LE.WINDOW(MCTYPE)) then
            IB2=IB1
            if (IB2.GT.N) then
               IB2=N
            endif
            IB1=1
            IT1=N*(IP-1)+IB1
            IT2=N*(IP-1)+IB2
            P1(1)=R(IT2,1)
            P1(2)=R(IT2,2)
            P1(3)=R(IT2,3)	  
         else
            IB1=N+WINDOW(MCTYPE)+1-IB1
            if (IB1.LT.1) then
               IB1=1
            endif
            IB2=N
            IT1=N*(IP-1)+IB1
            IT2=N*(IP-1)+IB2
            P1(1)=R(IT1,1)
            P1(2)=R(IT1,2)
            P1(3)=R(IT1,3)	  
         endif             
             
         ALPHA=2.*PI*grnd()
         BETA=acos(2.*grnd()-1.)
         TA(1)=sin(BETA)*cos(ALPHA)
         TA(2)=sin(BETA)*sin(ALPHA)
         TA(3)=cos(BETA)
             
         ALPHA=MCAMP(3)*(grnd()-0.5)
             
         ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
         ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
         ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
		 -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)
             
         ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
         ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
         ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
		 -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)
             
         ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
         ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
         ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
         ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
		 -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)
             
         DO 30 I=IT1,IT2
            RP(I,1)=ROT(1,4)+ROT(1,1)*R(I,1)+ROT(1,2)*R(I,2)+ROT(1,3)*R(I,3)
            RP(I,2)=ROT(2,4)+ROT(2,1)*R(I,1)+ROT(2,2)*R(I,2)+ROT(2,3)*R(I,3)
            RP(I,3)=ROT(3,4)+ROT(3,1)*R(I,1)+ROT(3,2)*R(I,2)+ROT(3,3)*R(I,3)
            UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
            UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
            UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)		

            !KHOU
            WP(I,1)=ROT(1,1)*W(I,1)+ROT(1,2)*W(I,2)+ROT(1,3)*W(I,3)
            WP(I,2)=ROT(2,1)*W(I,1)+ROT(2,2)*W(I,2)+ROT(2,3)*W(I,3)
            WP(I,3)=ROT(3,1)*W(I,1)+ROT(3,2)*W(I,2)+ROT(3,3)*W(I,3)
 30      CONTINUE

!     Perform rotate move (MCTYPE 4)
             
      elseif (MCTYPE.EQ.4) then
         
         IP=nint(0.5+grnd()*NP)
         IB1=nint(0.5+grnd()*N)
         IB2=IB1
         IT1=N*(IP-1)+IB1
         IT2=N*(IP-1)+IB2
         
         ALPHA=2.*PI*grnd()
         BETA=acos(2.*grnd()-1.)
         TA(1)=sin(BETA)*cos(ALPHA)
         TA(2)=sin(BETA)*sin(ALPHA)
         TA(3)=cos(BETA)
             
         ALPHA=MCAMP(4)*(grnd()-0.5)
             
         ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
         ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
         ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
		 -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)
             
         ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
         ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
         ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
		 -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)
             
         ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
         ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
         ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
         ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
		 -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)
             
         I=IT1
         UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
         UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
         UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)

         !KHOU
         WP(I,1)=ROT(1,1)*W(I,1)+ROT(1,2)*W(I,2)+ROT(1,3)*W(I,3)
         WP(I,2)=ROT(2,1)*W(I,1)+ROT(2,2)*W(I,2)+ROT(2,3)*W(I,3)
         WP(I,3)=ROT(3,1)*W(I,1)+ROT(3,2)*W(I,2)+ROT(3,3)*W(I,3)

!     Perform a full chain rotation (MCTYPE 5)

      elseif (MCTYPE.EQ.5) then
         
         IP=nint(0.5+grnd()*NP)
         IB1=1
         IB2=N
         IT1=N*(IP-1)+IB1
         IT2=N*(IP-1)+IB2
         
         ALPHA=2.*PI*grnd()
         BETA=acos(2.*grnd()-1.)
         TA(1)=sin(BETA)*cos(ALPHA)
         TA(2)=sin(BETA)*sin(ALPHA)
         TA(3)=cos(BETA)
         
         ALPHA=MCAMP(5)*(grnd()-0.5)
             
         ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
         ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
         ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
		 -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)
             
         ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
         ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
         ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
		 -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)
             
         ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
         ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
         ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
         ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
		 -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)
             
         DO 40 I=IT1,IT2
            RP(I,1)=ROT(1,4)+ROT(1,1)*R(I,1)+ROT(1,2)*R(I,2)+ROT(1,3)*R(I,3)
            RP(I,2)=ROT(2,4)+ROT(2,1)*R(I,1)+ROT(2,2)*R(I,2)+ROT(2,3)*R(I,3)
            RP(I,3)=ROT(3,4)+ROT(3,1)*R(I,1)+ROT(3,2)*R(I,2)+ROT(3,3)*R(I,3)
            UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
            UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
            UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)		 

           !KHOU
            WP(I,1)=ROT(1,1)*W(I,1)+ROT(1,2)*W(I,2)+ROT(1,3)*W(I,3)
            WP(I,2)=ROT(2,1)*W(I,1)+ROT(2,2)*W(I,2)+ROT(2,3)*W(I,3)
            WP(I,3)=ROT(3,1)*W(I,1)+ROT(3,2)*W(I,2)+ROT(3,3)*W(I,3)
 40      CONTINUE

!     Perform full chain slide move (MCTYPE 6)
         
      elseif (MCTYPE.EQ.6) then
         
         IP=nint(0.5+grnd()*NP)
         IB1=1
         IB2=N
         IT1=N*(IP-1)+IB1
         IT2=N*(IP-1)+IB2
         
         DR(1)=MCAMP(6)*(grnd()-0.5)
         DR(2)=MCAMP(6)*(grnd()-0.5)
         DR(3)=MCAMP(6)*(grnd()-0.5)
         
         DO 50 I=IT1,IT2
            RP(I,1)=R(I,1)+DR(1)
            RP(I,2)=R(I,2)+DR(2)
            RP(I,3)=R(I,3)+DR(3)
            UP(I,1)=U(I,1)
            UP(I,2)=U(I,2)
            UP(I,3)=U(I,3)

            !KHOU
            WP(I,1)=W(I,1)
            WP(I,2)=W(I,2)
            WP(I,3)=W(I,3)
 50      CONTINUE    

!     KHOU: Perform a W-vector rotation (MCTYPE 7)
!     This selects a random bead and rotates it about the position vector U, so as to rotate W
!     without changing the invariant that U and W should be orthogonal.

      elseif (MCTYPE.EQ.7) then
         
         ! Select Random Bead
         IP=nint(0.5+grnd()*NP)
         IB1=nint(0.5+grnd()*N)
         IB2=IB1
         IT1=N*(IP-1)+IB1
         IT2=N*(IP-1)+IB2
         

         ! Set rotation axis to the orientation vector at U
         I = IT1
         TA(1) = U(I,1)
         TA(2) = U(I,2)
         TA(3) = U(I,3)
         
         ! Pick a random rotation angle
         ALPHA=MCAMP(7)*(grnd()-0.5)
             
         ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
         ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
         ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
       -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)
             
         ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
         ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
         ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
         ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
       -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)
             
         ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
         ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
         ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
         ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
       -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)

         UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
         UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
         UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)

         !KHOU
         WP(I,1)=ROT(1,1)*W(I,1)+ROT(1,2)*W(I,2)+ROT(1,3)*W(I,3)
         WP(I,2)=ROT(2,1)*W(I,1)+ROT(2,2)*W(I,2)+ROT(2,3)*W(I,3)
         WP(I,3)=ROT(3,1)*W(I,1)+ROT(3,2)*W(I,2)+ROT(3,3)*W(I,3)

         
      end if
      

      RETURN      
      END
      
!---------------------------------------------------------------!      
